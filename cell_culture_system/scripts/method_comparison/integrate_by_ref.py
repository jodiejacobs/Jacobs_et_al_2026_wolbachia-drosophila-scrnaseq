"""
integrate_by_ref.py
=====================
Integrates query datasets with a reference dataset.

Strategy
--------
1. Load reference (uninfected, has cyclum_stage/pseudotime/leiden).
   Extract log-normalised counts from reference.raw.X.
2. Load each query h5ad, extract raw counts from .raw.X, subset to
   infected cells (treatment != Ctrl).
3. Concatenate reference + query on shared genes.
4. Jointly normalise (already done in ref, redo from raw for consistency)
   -> restrict to reference HVGs -> scale -> PCA -> Harmony.
5. KNN label transfer: for each query cell find k nearest reference
   neighbours in Harmony PCA space, majority-vote their leiden label.
6. Map leiden -> cc_stage / cc_pseudotime via reference majority-vote.
7. Run Q1/Q2/Q3 analyses.

Run with:
    mamba activate scanpy
    python scripts/method_comparison/integrate_by_ref.py \
        --ref   results/integrated/integrated_uninfected_with_cellcycle.h5ad \
        --query results/filtered_h5ad/JW18wMel-SV*_pipseq.h5ad \
                results/filtered_h5ad/JW18wMel-SV*_10x.h5ad \
        --out_path results/integrated/titer_cellcycle.h5ad \
        --fig_dir  figures/titer_vs_cellcycle \
        --harmony_vars method replicate dataset
"""

import os
import glob
import argparse

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.sparse
import scipy.stats
from scipy.stats import kruskal, spearmanr, chi2_contingency, t as t_dist, chi2 as chi2_dist, norm
from itertools import combinations
from statsmodels.stats.multitest import multipletests
from statsmodels.nonparametric.smoothers_lowess import lowess
import anndata as ad
import scanpy as sc
import harmonypy as hm


# -----------------------------------------------------------------------------
# Colours
# -----------------------------------------------------------------------------

CC_ORDER  = ["g0/g1", "s", "g2/m"]
CC_COLORS = {"g0/g1": "#ea546f", "s": "#7bcdca", "g2/m": "#2d9bb4"}   # updated palette

def _cc_palette(stages):
    cmap = matplotlib.colormaps["tab10"]
    return {s: CC_COLORS.get(s, cmap(i % 10)) for i, s in enumerate(stages)}


# -----------------------------------------------------------------------------
# Publication UMAP dimensions
# -----------------------------------------------------------------------------

UMAP_W = 2.75   # inches
UMAP_H = 2.032  # inches


# -----------------------------------------------------------------------------
# P-value helpers -- underflow-safe
# -----------------------------------------------------------------------------

def spearman_p_exact(x, y):
    """Spearman rho + underflow-safe p via t-distribution logsf."""
    rho, _ = spearmanr(x, y)
    n = len(x)
    denom = max(1 - rho**2, 1e-15)
    t_stat = rho * np.sqrt((n - 2) / denom)
    if not np.isfinite(t_stat):
        return rho, 0.0, -np.inf, t_stat
    log_p  = t_dist.logsf(abs(t_stat), df=n - 2) + np.log(2)
    p_val  = np.exp(log_p)
    log10_p = log_p / np.log(10)
    return rho, p_val, log10_p, t_stat

def chi2_p_exact(stat, dof):
    log_p  = chi2_dist.logsf(stat, df=dof)
    p_val  = np.exp(log_p)
    return p_val, log_p / np.log(10)


def kw_p_exact(H, k):
    log_p = chi2_dist.logsf(H, df=k - 1)
    p_val = np.exp(log_p)
    return p_val, log_p / np.log(10)


def z_p_exact(z):
    log_p = norm.logsf(abs(z)) + np.log(2)
    p_val = np.exp(log_p)
    return p_val, log_p / np.log(10)


def format_p(p_val, log10_p):
    if not np.isfinite(log10_p):
        if p_val == 0.0:
            return "p<5e-324"
        return f"p={p_val:.2e}"
    if p_val == 0.0 or p_val < 5e-300:
        exp = int(np.floor(log10_p))
        man = 10 ** (log10_p - exp)
        return f"p={man:.2f}x10^{exp}"
    if p_val < 0.001:
        return f"p={p_val:.2e}"
    return f"p={p_val:.4f}"


def bh_adjust_log10(log10_p_raw_array):
    p_raw = np.array([10**lp for lp in log10_p_raw_array], dtype=float)
    p_raw_safe = np.where(p_raw == 0.0, np.finfo(float).tiny, p_raw)
    _, p_adj, _, _ = multipletests(p_raw_safe, method="fdr_bh")
    log10_p_adj = np.where(
        p_adj == 0.0,
        np.log10(np.finfo(float).tiny),
        np.log10(np.maximum(p_adj, np.finfo(float).tiny)),
    )
    return p_adj, log10_p_adj


def identify_reference_markers(ref_full, fig_dir, sample,
                               groupby="leiden", n_top=5,
                               method="wilcoxon"):
    print(f"\n-- Reference marker genes ({groupby}) --")

    if groupby not in ref_full.obs.columns:
        raise ValueError(f"'{groupby}' not found in reference obs")

    groups = ref_full.obs[groupby].dropna().astype(str).unique().tolist()
    if len(groups) < 2:
        print(f"   SKIP: need at least two {groupby} groups to compute markers")
        return None

    adata = ref_full.copy()
    adata.obs[groupby] = adata.obs[groupby].astype("category")

    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        method=method,
        use_raw=False,
    )

    sc.pl.rank_genes_groups_dotplot(
        adata,
        groupby=groupby,
        n_genes=n_top,
        standard_scale="var",
        show=False,
        save=f"_{sample}_{groupby}_markers_dotplot.pdf",
    )

    marker_rows = []
    for group in adata.uns["rank_genes_groups"]["names"].dtype.names:
        df = sc.get.rank_genes_groups_df(adata, group=group).head(n_top).copy()
        if df.empty:
            continue
        df.insert(0, "cluster", group)
        df.insert(1, "rank", np.arange(1, len(df) + 1))
        marker_rows.append(df)
        print(f"   Cluster {group}: {', '.join(df['names'].astype(str).tolist())}")

    if marker_rows:
        marker_df = pd.concat(marker_rows, ignore_index=True)
        out_csv = os.path.join(fig_dir, f"{sample}_{groupby}_markers.csv")
        marker_df.to_csv(out_csv, index=False)
        print(f"   -> {out_csv}")
    else:
        marker_df = pd.DataFrame()

    return adata, marker_df


# -----------------------------------------------------------------------------
# Step 1 -- Load reference, extract log-normalised counts from .raw
# -----------------------------------------------------------------------------

def load_reference(ref_path, stage_col, pseudotime_col):
    """Load reference and extract raw counts for joint re-processing."""
    print(f"\n-- Loading reference: {ref_path} --")
    ref_full = sc.read_h5ad(ref_path)
    print(f"   {ref_full.n_obs} cells x {ref_full.n_vars} genes")
    print(f"   Leiden clusters : {ref_full.obs['leiden'].nunique()}")

    for col in [stage_col, pseudotime_col, "leiden"]:
        if col not in ref_full.obs.columns:
            raise ValueError(f"'{col}' not in reference obs. "
                             f"Available: {list(ref_full.obs.columns)}")

    print(f"   CC stage distribution (reference, Cyclum):")
    print(ref_full.obs[stage_col].value_counts().to_string())

    if ref_full.raw is None:
        raise ValueError(
            "Reference has no .raw -- re-run filter.py with the updated "
            "analyze_filtered_adata() that saves adata.raw = adata before normalisation."
        )

    raw_X = ref_full.raw.X
    if scipy.sparse.issparse(raw_X):
        raw_X = raw_X.toarray()
    raw_X = raw_X.astype(np.float32)

    ref_raw = ad.AnnData(
        X=scipy.sparse.csr_matrix(raw_X),
        obs=ref_full.obs.copy(),
        var=ref_full.raw.var.copy(),
    )
    ref_raw.obs["dataset"] = "reference"
    ref_raw.obs["cc_stage_source"] = "cyclum"
    print(f"   Reference raw counts: {ref_raw.n_obs} cells x {ref_raw.n_vars} genes")

    return ref_raw, ref_full


# -----------------------------------------------------------------------------
# Step 2 -- Load query files, extract raw counts, subset to infected cells
# -----------------------------------------------------------------------------

def load_query_files(query_paths, ref_condition, titer_col,
                     treatment_col="treatment", infected_only=True):
    adatas = []

    for path in query_paths:
        print(f"\n   Loading query: {path}")
        adata = sc.read_h5ad(path)
        print(f"   {adata.n_obs} cells x {adata.n_vars} genes")

        if adata.raw is None:
            raise ValueError(
                f"{path} has no .raw -- re-run filter.py with the updated "
                "analyze_filtered_adata() that saves adata.raw = adata "
                "before normalisation."
            )

        raw_X = adata.raw.X
        if scipy.sparse.issparse(raw_X):
            raw_X = raw_X.toarray()
        raw_X = raw_X.astype(np.float32)

        a = ad.AnnData(
            X=scipy.sparse.csr_matrix(raw_X),
            obs=adata.obs.copy(),
            var=adata.raw.var.copy(),
        )

        basename = os.path.basename(path).replace(".h5ad", "")
        if "method" not in a.obs.columns:
            method = "pipseq" if "pipseq" in basename.lower() else "10x"
            a.obs["method"] = method
            print(f"   Inferred method='{method}' from filename")
        if "replicate" not in a.obs.columns:
            parts = basename.split("-")
            rep = parts[-1].split("_")[0] if len(parts) >= 3 else "unknown"
            a.obs["replicate"] = rep
            print(f"   Inferred replicate='{rep}' from filename")
        if "cell_line" not in a.obs.columns:
            cell_line = basename.split("-")[0]
            a.obs["cell_line"] = cell_line
            print(f"   Inferred cell_line='{cell_line}' from filename")
        if "treatment" not in a.obs.columns:
            treatment = "Ctrl" if "Ctrl" in basename else "SV"
            a.obs["treatment"] = treatment
            print(f"   Inferred treatment='{treatment}' from filename")

        if infected_only and treatment_col in a.obs.columns:
            mask = (
                ~a.obs["cell_line"].isin(ref_condition) |
                (a.obs[treatment_col] != "Ctrl")
            )
            print(f"   Infected cells: {mask.sum()}/{a.n_obs}")
            a = a[mask].copy()
        else:
            print(f"   Using all {a.n_obs} cells (infected_only=False)")

        if a.n_obs == 0:
            print(f"   WARNING: no cells remain after filtering -- skipping {path}")
            continue

        if titer_col not in a.obs.columns:
            print(f"   WARNING: '{titer_col}' not in obs -- titer analyses will skip this file")

        a.obs["dataset"]          = "query"
        a.obs["source_file"]      = basename
        a.obs["cc_stage_source"]  = "knn_transferred"
        a.obs_names = [f"{basename}_{bc}" for bc in a.obs_names]

        adatas.append(a)
        print(f"   Kept {a.n_obs} cells")

    if not adatas:
        raise ValueError("No query cells loaded -- check paths and --ref_condition")

    print(f"\n-- Concatenating {len(adatas)} query files --")
    query_raw = ad.concat(adatas, join="outer", index_unique=None)
    query_raw.obs_names_make_unique()
    print(f"   Total query cells: {query_raw.n_obs}")

    return query_raw


# -----------------------------------------------------------------------------
# Step 3 -- Joint normalisation -> HVG restriction -> scale -> PCA -> Harmony
# -----------------------------------------------------------------------------

def joint_preprocess_and_harmony(ref_raw, query_raw, ref_full,
                                  harmony_vars, n_pcs=30):
    print(f"\n-- Joint preprocessing --")

    combined = ad.concat(
        [ref_raw, query_raw],
        join="outer",
        index_unique=None,
        label="dataset",
        keys=["reference", "query"],
    )
    combined.obs_names_make_unique()

    if scipy.sparse.issparse(combined.X):
        combined.X = combined.X.toarray()
    combined.X = np.nan_to_num(combined.X.astype(np.float32), nan=0.0)
    combined.X = scipy.sparse.csr_matrix(combined.X)

    ref_mask = combined.obs["dataset"] == "reference"
    print(f"   Combined: {combined.n_obs} cells x {combined.n_vars} genes")
    print(f"   Reference: {ref_mask.sum()}  Query: {(~ref_mask).sum()}")

    print("   Normalising (1e4 per cell) + log1p ...")
    sc.pp.normalize_total(combined, target_sum=1e4)
    sc.pp.log1p(combined)

    ref_hvgs     = ref_full.var_names[ref_full.var["highly_variable"]].tolist()
    hvgs_present = [g for g in ref_hvgs if g in combined.var_names]
    print(f"   Reference HVGs: {len(ref_hvgs)} total, "
          f"{len(hvgs_present)} present in combined dataset")

    if len(hvgs_present) < 100:
        raise ValueError(
            f"Only {len(hvgs_present)} reference HVGs found in combined dataset. "
            "Check that query files share the same gene universe as the reference."
        )

    combined.var["highly_variable"] = combined.var_names.isin(hvgs_present)

    print("   Scaling (max_value=10) ...")
    sc.pp.scale(combined, max_value=10)

    print(f"   PCA ({n_pcs} components) ...")
    sc.tl.pca(combined, n_comps=n_pcs, use_highly_variable=True, svd_solver="arpack")

    missing_vars = [v for v in harmony_vars if v not in combined.obs.columns]
    if missing_vars:
        print(f"   WARNING: harmony_vars {missing_vars} not in obs -- dropping them")
        harmony_vars = [v for v in harmony_vars if v in combined.obs.columns]
    if not harmony_vars:
        raise ValueError("No valid harmony_vars remain -- cannot run Harmony")

    for v in harmony_vars:
        combined.obs[v] = combined.obs[v].astype(str).fillna("unknown")

    print(f"   Harmony correction on: {harmony_vars} ...")
    pca_matrix = combined.obsm["X_pca"]
    meta       = combined.obs[harmony_vars].copy()
    ho = hm.run_harmony(pca_matrix, meta, harmony_vars,
                        max_iter_harmony=30, random_state=42)
    combined.obsm["X_pca_harmony"] = ho.Z_corr.T

    print("   Harmony complete")
    return combined, ref_mask


# -----------------------------------------------------------------------------
# Step 4 -- KNN label transfer
# -----------------------------------------------------------------------------

def knn_label_transfer(combined, ref_mask, ref_full, k=15):
    from sklearn.neighbors import NearestNeighbors

    print(f"\n-- KNN label transfer (k={k}) --")

    ref_pca   = combined.obsm["X_pca_harmony"][ref_mask]
    query_pca = combined.obsm["X_pca_harmony"][~ref_mask]

    ref_leiden = ref_full.obs["leiden"].astype(str).values

    nbrs = NearestNeighbors(n_neighbors=k, metric="euclidean", n_jobs=-1)
    nbrs.fit(ref_pca)
    distances, indices = nbrs.kneighbors(query_pca)

    transferred = []
    confidence  = []
    for row_idx in indices:
        neighbour_labels = ref_leiden[row_idx]
        counts = pd.Series(neighbour_labels).value_counts()
        winner = counts.index[0]
        transferred.append(winner)
        confidence.append(counts.iloc[0] / k)

    query = combined[~ref_mask].copy()
    query.obs["leiden_ref"]     = transferred
    query.obs["knn_confidence"] = confidence
    query.obsm["X_pca_harmony"] = query_pca

    print("   Computing UMAP on combined Harmony embedding ...")
    sc.pp.neighbors(combined, use_rep="X_pca_harmony", n_neighbors=30)
    sc.tl.umap(combined)
    query.obsm["X_umap"] = combined.obsm["X_umap"][~ref_mask]

    print(f"   Transferred leiden_ref distribution:")
    print(pd.Series(transferred).value_counts().sort_index().to_string())
    print(f"   Mean KNN confidence: {np.mean(confidence):.3f}")

    low_conf = np.mean(np.array(confidence) < 0.5)
    if low_conf > 0.2:
        print(f"   WARNING: {low_conf*100:.1f}% of cells have KNN confidence < 0.5 "
              f"-- batch correction may be insufficient")

    return query, combined


# -----------------------------------------------------------------------------
# Step 5 -- Map leiden_ref -> CC stage + pseudotime
# -----------------------------------------------------------------------------

def assign_cc_from_reference(query, ref_full,
                              stage_col="cyclum_stage",
                              pseudotime_col="cyclum_pseudotime",
                              leiden_col="leiden"):
    ref_obs = ref_full.obs[[leiden_col, stage_col]].copy()
    ref_obs[leiden_col] = ref_obs[leiden_col].astype(str)

    def majority(x):
        return x.value_counts().idxmax()
    def purity(x):
        vc = x.value_counts()
        return vc.iloc[0] / vc.sum()

    stage_map = (ref_obs.groupby(leiden_col)[stage_col]
                        .agg(cc_stage=majority, purity=purity)
                        .reset_index()
                        .rename(columns={leiden_col: "cluster"}))

    if pseudotime_col in ref_full.obs.columns:
        pt_map = (ref_full.obs[[leiden_col, pseudotime_col]]
                          .copy()
                          .assign(**{leiden_col: ref_full.obs[leiden_col].astype(str)})
                          .groupby(leiden_col)[pseudotime_col]
                          .median()
                          .reset_index()
                          .rename(columns={leiden_col: "cluster",
                                           pseudotime_col: "median_pseudotime"}))
        stage_map = stage_map.merge(pt_map, on="cluster", how="left")
    else:
        stage_map["median_pseudotime"] = np.nan

    n_ref = ref_obs[leiden_col].value_counts().rename("n_ref_cells").reset_index()
    n_ref.columns = ["cluster", "n_ref_cells"]
    stage_map = stage_map.merge(n_ref, on="cluster", how="left")

    print("\n-- Cluster -> CC stage mapping --")
    print(stage_map.sort_values("cluster").to_string(index=False))

    low_purity = stage_map[stage_map["purity"] < 0.6]
    if len(low_purity):
        print(f"\n   WARNING: {len(low_purity)} clusters have CC stage purity < 60%:")
        print(low_purity[["cluster", "cc_stage", "purity"]].to_string(index=False))

    cluster_to_stage = dict(zip(stage_map["cluster"], stage_map["cc_stage"]))
    cluster_to_pt    = dict(zip(stage_map["cluster"], stage_map["median_pseudotime"]))

    query.obs["leiden_ref"]    = query.obs["leiden_ref"].astype(str)
    query.obs["cc_stage"]      = query.obs["leiden_ref"].map(cluster_to_stage)
    query.obs["cc_pseudotime"] = query.obs["leiden_ref"].map(cluster_to_pt).astype(float)

    n_unmapped = query.obs["cc_stage"].isna().sum()
    if n_unmapped:
        print(f"\n   WARNING: {n_unmapped} query cells did not map to a known cluster")

    print(f"\n   Query CC stage distribution (KNN-transferred):")
    print(query.obs["cc_stage"].value_counts().to_string())

    return query, stage_map


# -----------------------------------------------------------------------------
# Q4 -- Per-sample CC distribution table + titer-shift statistics
# -----------------------------------------------------------------------------

def _infer_sample_label(obs_row):
    cell_line  = str(obs_row.get("cell_line", ""))
    treatment  = str(obs_row.get("treatment", ""))
    timepoint  = str(obs_row.get("timepoint", obs_row.get("day", obs_row.get("dpi", ""))))
    source     = str(obs_row.get("source_file", ""))

    if treatment.lower() in ("ctrl", "control", "uninfected"):
        return f"{cell_line} Ctrl" if cell_line else "Ctrl"
    if timepoint and timepoint not in ("", "nan", "NA"):
        return f"{cell_line} SV {timepoint}" if cell_line else f"SV {timepoint}"
    if source:
        return source
    return cell_line or "unknown"


def q4_cc_distribution_by_sample(
    query_obs,
    ref_full,
    fig_dir,
    sample,
    stage_col="cyclum_stage",
    query_stage_col="cc_stage",
    titer_col="wolbachia_titer",
    sample_label_col=None,
    n_titer_bins=5,
):
    print(f"\n-- Q4: per-sample CC distribution + titer-shift test --")
    os.makedirs(fig_dir, exist_ok=True)

    stages = [s for s in CC_ORDER if True]
    pal    = _cc_palette(CC_ORDER)

    ref_obs = ref_full.obs.copy()
    ref_obs["_sample_label"] = "Uninfected"
    ref_obs["_cc_stage"]     = ref_obs[stage_col].astype(str).str.strip().str.lower()

    qobs = query_obs.copy()
    qobs["_cc_stage"] = qobs[query_stage_col].astype(str).str.strip().str.lower()

    if sample_label_col and sample_label_col in qobs.columns:
        qobs["_sample_label"] = qobs[sample_label_col].astype(str)
    else:
        label_candidates = ["condition", "sample", "timepoint", "day", "dpi",
                            "treatment", "source_file"]
        available = [c for c in label_candidates if c in qobs.columns]
        if available:
            def _make_label(row):
                parts = []
                cl = str(row.get("cell_line", "")).strip()
                tx = str(row.get("treatment", "")).strip()
                tp = str(row.get("timepoint",
                        row.get("day",
                        row.get("dpi", "")))).strip()
                if cl and cl not in ("nan", "NA"):
                    parts.append(cl)
                if tx and tx not in ("nan", "NA", ""):
                    parts.append(tx)
                if tp and tp not in ("nan", "NA", ""):
                    parts.append(tp)
                return " ".join(parts) if parts else str(row.get("source_file", "unknown"))

            qobs["_sample_label"] = qobs.apply(_make_label, axis=1)
        else:
            qobs["_sample_label"] = qobs.get("source_file", "query").astype(str)

    all_rows = pd.concat([
        ref_obs[["_sample_label", "_cc_stage"]],
        qobs[["_sample_label", "_cc_stage"]],
    ], ignore_index=True)

    all_rows = all_rows[all_rows["_cc_stage"].notna() &
                        (all_rows["_cc_stage"] != "nan") &
                        (all_rows["_cc_stage"] != "NA")]

    observed_stages = [s for s in CC_ORDER if s in all_rows["_cc_stage"].unique()]
    observed_stages += [s for s in all_rows["_cc_stage"].unique()
                        if s not in CC_ORDER and s not in observed_stages]

    query_samples = sorted(qobs["_sample_label"].unique().tolist())
    all_samples   = ["Uninfected"] + query_samples

    ct = pd.crosstab(
        all_rows["_cc_stage"], all_rows["_sample_label"]
    ).reindex(index=observed_stages, columns=all_samples, fill_value=0)

    ct["ALL"] = ct.sum(axis=1)
    ct.to_csv(os.path.join(fig_dir, f"q4_cc_distribution_counts_{sample}.csv"))

    col_totals = ct.sum(axis=0)

    def fmt_cell(n, total):
        pct = 100.0 * n / total if total > 0 else 0.0
        return f"{n} ({pct:.1f}%)"

    fmt_rows = {}
    for stage in observed_stages:
        row = {}
        for col in ct.columns:
            row[col] = fmt_cell(ct.loc[stage, col], col_totals[col])
        fmt_rows[stage] = row

    fmt_df = pd.DataFrame(fmt_rows).T
    fmt_df.index.name = "cyclum_stage"
    fmt_df = fmt_df[ct.columns]

    out_dist = os.path.join(fig_dir, f"q4_cc_distribution_{sample}.csv")
    fmt_df.to_csv(out_dist)
    print(f"   -> {out_dist}")
    print(fmt_df.to_string())

    chi_rows = []
    ref_counts = ct.loc[observed_stages, "Uninfected"].values

    for qs in query_samples:
        q_counts = ct.loc[observed_stages, qs].values
        sub_ct   = np.vstack([ref_counts, q_counts])
        nonzero = (sub_ct.sum(axis=0) > 0)
        sub_ct  = sub_ct[:, nonzero]
        if sub_ct.shape[1] < 2:
            continue
        chi2_s, _, dof_s, _ = chi2_contingency(sub_ct)
        p_s, log10_p_s = chi2_p_exact(chi2_s, dof_s)
        chi_rows.append({
            "comparison":   f"Uninfected vs {qs}",
            "chi2":         chi2_s,
            "dof":          dof_s,
            "p_raw":        p_s,
            "log10_p_raw":  log10_p_s,
            "n_ref":        ref_counts.sum(),
            "n_query":      q_counts.sum(),
        })

    ct_core = ct.loc[observed_stages, all_samples]
    nonzero_cols = ct_core.columns[ct_core.sum(axis=0) > 0]
    ct_core = ct_core[nonzero_cols]
    chi2_all, _, dof_all, _ = chi2_contingency(ct_core.values)
    p_all, log10_p_all = chi2_p_exact(chi2_all, dof_all)
    print(f"\n   Overall chi-squared (all samples):")
    print(f"   chi2={chi2_all:.3f}, df={dof_all}, {format_p(p_all, log10_p_all)}")

    if chi_rows:
        chi_df = pd.DataFrame(chi_rows)
        p_adj, log10_p_adj = bh_adjust_log10(chi_df["log10_p_raw"].values)
        chi_df["p_adj_BH"]    = p_adj
        chi_df["log10_p_adj"] = log10_p_adj
        chi_df["significant"] = chi_df["log10_p_adj"] < np.log10(0.05)

        overall_row = pd.DataFrame([{
            "comparison":   "ALL SAMPLES (overall)",
            "chi2":         chi2_all,
            "dof":          dof_all,
            "p_raw":        p_all,
            "log10_p_raw":  log10_p_all,
            "p_adj_BH":     np.nan,
            "log10_p_adj":  np.nan,
            "significant":  p_all < 0.05,
            "n_ref":        ref_counts.sum(),
            "n_query":      ct.loc[observed_stages, query_samples].values.sum(),
        }])
        chi_df = pd.concat([chi_df, overall_row], ignore_index=True)
        out_stat = os.path.join(fig_dir, f"q4_titer_shift_chisq_{sample}.csv")
        chi_df.to_csv(out_stat, index=False)
        print(f"   -> {out_stat}")
        disp = chi_df.copy()
        disp["p_str"] = [format_p(p, lp) for p, lp in
                         zip(disp["p_raw"], disp["log10_p_raw"])]
        print(disp[["comparison", "chi2", "dof", "p_str", "significant"]].to_string(index=False))

    if titer_col in qobs.columns:
        _q4_titer_bin_chisq(qobs, fig_dir, sample,
                            titer_col=titer_col,
                            stage_col="_cc_stage",
                            observed_stages=observed_stages,
                            n_titer_bins=n_titer_bins)

    prop = ct.loc[observed_stages, all_samples].div(
        ct.loc[observed_stages, all_samples].sum(axis=0), axis=1
    ) * 100

    fig, ax = plt.subplots(figsize=(max(8, len(all_samples) * 1.4), 5))
    bottom = np.zeros(len(all_samples))
    for stage in observed_stages:
        vals = prop.loc[stage, all_samples].values.astype(float)
        ax.bar(np.arange(len(all_samples)), vals, bottom=bottom,
               color=pal.get(stage, "#aaaaaa"), label=stage, width=0.7)
        for xi, (v, b) in enumerate(zip(vals, bottom)):
            if v > 5:
                ax.text(xi, b + v / 2, f"{v:.0f}%",
                        ha="center", va="center", fontsize=7,
                        color="white", fontweight="bold")
        bottom += vals

    ax.set_xticks(np.arange(len(all_samples)))
    ax.set_xticklabels(all_samples, rotation=35, ha="right", fontsize=9)
    ax.set_ylabel("% of cells")
    ax.set_ylim(0, 105)
    ax.set_title(
        f"CC phase composition by sample\n"
        f"Overall chi2={chi2_all:.1f}, df={dof_all}, "
        f"{format_p(p_all, log10_p_all)}",
        fontweight="bold",
    )
    ax.legend(title="CC stage", bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=9)

    for xi, s in enumerate(all_samples):
        n = col_totals[s]
        ax.text(xi, 102, f"n={n}", ha="center", va="bottom", fontsize=7, color="grey")

    plt.tight_layout()
    out_fig = os.path.join(fig_dir, f"q4_cc_distribution_{sample}.pdf")
    plt.savefig(out_fig, bbox_inches="tight", dpi=150)
    plt.close()
    print(f"   -> {out_fig}")

    return fmt_df, chi_df if chi_rows else None


def _q4_titer_bin_chisq(qobs, fig_dir, sample,
                         titer_col, stage_col, observed_stages, n_titer_bins):
    print(f"\n   -- Q4 titer-bin analysis --")

    df = qobs[[titer_col, stage_col]].copy()
    df[titer_col] = pd.to_numeric(df[titer_col], errors="coerce")
    df = df.dropna()
    df[stage_col] = df[stage_col].astype(str).str.strip().str.lower()

    if len(df) < 20:
        print("   SKIP: fewer than 20 cells with titer data"); return

    _, bin_edges = pd.qcut(df[titer_col], q=n_titer_bins,
                           retbins=True, duplicates="drop")
    actual_n = len(bin_edges) - 1
    if actual_n < 2:
        print("   SKIP: too few unique titer values"); return

    bin_labels = [f"Q{i+1}" for i in range(actual_n)]
    df["titer_bin"] = pd.cut(df[titer_col], bins=bin_edges,
                             labels=bin_labels, include_lowest=True)
    df = df.dropna(subset=["titer_bin"])

    ct = (pd.crosstab(df["titer_bin"], df[stage_col])
            .reindex(columns=[s for s in observed_stages if s in df[stage_col].unique()],
                     fill_value=0))

    chi2_s, _, dof_s, _ = chi2_contingency(ct.values)
    p_s, log10_p_s = chi2_p_exact(chi2_s, dof_s)
    print(f"   Titer-bin chi-squared: chi2={chi2_s:.3f}, df={dof_s}, "
          f"{format_p(p_s, log10_p_s)}")

    pd.DataFrame({
        "test":    ["chi2_titer_bins"],
        "chi2":    [chi2_s],
        "dof":     [dof_s],
        "p_value": [p_s],
        "log10_p": [log10_p_s],
        "n_bins":  [actual_n],
        "n_cells": [len(df)],
        "interpretation": [
            "Chi-squared tests whether CC phase proportions differ across "
            f"{actual_n} quantile bins of wMel titer"
        ],
    }).to_csv(os.path.join(fig_dir, f"q4_titer_bin_chisq_{sample}.csv"), index=False)

    prop = ct.div(ct.sum(axis=1), axis=0) * 100
    pal  = _cc_palette(observed_stages)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    bottom = np.zeros(len(ct))
    for stage in ct.columns:
        vals = prop[stage].values
        axes[0].bar(np.arange(len(ct)), vals, bottom=bottom,
                    color=pal.get(stage, "#aaaaaa"), label=stage, width=0.7)
        bottom += vals
    bin_ranges = (df.groupby("titer_bin", observed=True)[titer_col]
                    .agg(lo="min", hi="max").reset_index())
    br_map = {row["titer_bin"]: f"{row['lo']:.2f}-{row['hi']:.2f}"
              for _, row in bin_ranges.iterrows()}
    axes[0].set_xticks(np.arange(len(ct)))
    axes[0].set_xticklabels(
        [f"{b}\n({br_map.get(b,'')})" for b in ct.index], fontsize=8)
    axes[0].set_xlabel("wMel titer quantile")
    axes[0].set_ylabel("% of cells")
    axes[0].set_title(
        f"CC phase vs titer (query cells)\n"
        f"chi2={chi2_s:.1f}, df={dof_s}, {format_p(p_s, log10_p_s)}")
    axes[0].legend(title="CC stage", bbox_to_anchor=(1.01, 1), loc="upper left")

    sns.heatmap(prop.T.reindex(observed_stages).dropna(how="all"),
                annot=True, fmt=".1f", cmap="YlOrRd",
                linewidths=0.5, ax=axes[1],
                cbar_kws={"label": "% of cells in titer bin"})
    axes[1].set_xlabel("wMel titer bin")
    axes[1].set_ylabel("CC stage")
    axes[1].set_title("CC proportions per titer bin")

    plt.tight_layout()
    out = os.path.join(fig_dir, f"q4_titer_bin_cc_{sample}.pdf")
    plt.savefig(out, bbox_inches="tight", dpi=150)
    plt.close()
    print(f"   -> {out}")


# -----------------------------------------------------------------------------
# Q1 -- titer vs cell-cycle pseudotime
# -----------------------------------------------------------------------------

def q1_titer_vs_pseudotime(obs, fig_dir, sample,
                            titer_col="wolbachia_titer",
                            pseudotime_col="cc_pseudotime",
                            stage_col="cc_stage",
                            n_bins=8):
    print(f"\n-- Q1: titer vs CC pseudotime --")

    missing = [c for c in [titer_col, pseudotime_col] if c not in obs.columns]
    if missing:
        print(f"   SKIP: {missing} not in obs"); return

    df = obs[[titer_col, pseudotime_col, stage_col]].copy()
    df[titer_col]      = pd.to_numeric(df[titer_col],      errors="coerce")
    df[pseudotime_col] = pd.to_numeric(df[pseudotime_col], errors="coerce")
    df = df.dropna()
    print(f"   {len(df)} cells")

    if len(df) < 20:
        print("   SKIP: fewer than 20 cells"); return

    n_unique_pt = df[pseudotime_col].nunique()
    print(f"   Unique pseudotime values: {n_unique_pt}")

    pt  = df[pseudotime_col].values
    tit = df[titer_col].values

    stages = [s for s in CC_ORDER if s in df[stage_col].unique()]
    stages += [s for s in df[stage_col].unique() if s not in CC_ORDER]
    pal    = _cc_palette(stages)

    rho, p_val, log10_p, t_stat = spearman_p_exact(pt, tit)
    p_str = format_p(p_val, log10_p)
    print(f"   Spearman rho={rho:.3f}, {p_str}  "
          f"(t={t_stat:.2f}, df={len(df)-2}, log10_p={log10_p:.1f})")

    pd.DataFrame({
        "spearman_rho":  [rho],
        "t_stat":        [t_stat],
        "df":            [len(df) - 2],
        "p_value":       [p_val],
        "log10_p":       [log10_p],
        "n_cells":       [len(df)],
        "note": ["cc_pseudotime is cluster-median from reference; "
                 "log10_p exact even when p_value underflows to 0.0"],
    }).to_csv(os.path.join(fig_dir, f"q1_spearman_{sample}.csv"), index=False)

    do_lowess = n_unique_pt >= 5
    if do_lowess:
        order    = np.argsort(pt)
        smoothed = lowess(tit[order], pt[order], frac=0.3, return_sorted=True)

    n_bins_actual = min(n_bins, n_unique_pt)
    do_bins = n_bins_actual >= 2
    if do_bins:
        edges   = np.linspace(pt.min(), pt.max(), n_bins_actual + 1)
        centers = (edges[:-1] + edges[1:]) / 2
        bin_lbl = [f"{e:.2f}" for e in edges[:-1]]
        df["_pt_bin"] = pd.cut(df[pseudotime_col], bins=edges,
                               labels=bin_lbl, include_lowest=True)
        bsum = (df.dropna(subset=["_pt_bin"])
                  .groupby("_pt_bin", observed=True)[titer_col]
                  .agg(n_cells="count", median="median", mean="mean",
                       q25=lambda x: x.quantile(0.25),
                       q75=lambda x: x.quantile(0.75))
                  .reset_index())
        occupied = bsum["_pt_bin"].astype(str).values
        bsum["bin_center"] = [centers[bin_lbl.index(b)]
                              for b in occupied if b in bin_lbl]
        bsum.to_csv(os.path.join(fig_dir, f"q1_bin_summary_{sample}.csv"), index=False)

    pi_ticks  = [0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi]
    pi_labels = ["0", "pi/2", "pi", "3pi/2", "2pi"]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    for stage in stages:
        m = df[stage_col] == stage
        axes[0].scatter(df.loc[m, pseudotime_col], df.loc[m, titer_col],
                        c=[pal[stage]], s=3, alpha=0.4, label=stage, rasterized=True)
    axes[0].set_xlabel("CC pseudotime (cluster-median, radians)")
    axes[0].set_ylabel("wMel titer")
    axes[0].set_title(f"A  Titer vs pseudotime\nSpearman rho={rho:.3f}, {p_str}")
    axes[0].set_xticks(pi_ticks); axes[0].set_xticklabels(pi_labels)
    axes[0].legend(title="CC stage", fontsize=8, markerscale=3)

    if do_lowess:
        axes[1].scatter(pt, tit, c="lightgrey", s=2, alpha=0.25, rasterized=True)
        axes[1].plot(smoothed[:, 0], smoothed[:, 1],
                     color="#d62728", lw=2, label="LOWESS (frac=0.3)")
        axes[1].set_xticks(pi_ticks); axes[1].set_xticklabels(pi_labels)
        axes[1].legend()
    else:
        axes[1].text(0.5, 0.5, f"Too few unique pseudotime\nvalues (n={n_unique_pt})",
                     ha="center", va="center", transform=axes[1].transAxes)
    axes[1].set_xlabel("CC pseudotime (cluster-median, radians)")
    axes[1].set_ylabel("wMel titer")
    axes[1].set_title("B  LOWESS smoothed trend")

    if do_bins and len(bsum) >= 2:
        x_ = bsum["bin_center"].values
        axes[2].plot(x_, bsum["median"].values, "o-", color="#d62728", lw=2, ms=5)
        axes[2].fill_between(x_, bsum["q25"].values, bsum["q75"].values,
                             alpha=0.25, color="#d62728", label="IQR")
        axes[2].set_xticks(x_)
        axes[2].set_xticklabels([f"{v:.2f}" for v in x_], rotation=45)
        axes[2].legend()
    else:
        axes[2].text(0.5, 0.5, "Insufficient pseudotime resolution",
                     ha="center", va="center", transform=axes[2].transAxes)
    axes[2].set_xlabel("Pseudotime bin center (radians)")
    axes[2].set_ylabel("wMel titer")
    axes[2].set_title("C  Median titer +/- IQR")

    plt.suptitle("Q1 -- wMel titer vs cell-cycle pseudotime", fontweight="bold")
    plt.tight_layout()
    out = os.path.join(fig_dir, f"q1_titer_vs_pseudotime_{sample}.pdf")
    plt.savefig(out, bbox_inches="tight", dpi=150)
    plt.close()
    print(f"   -> {out}")


# -----------------------------------------------------------------------------
# Q2 -- CC phase distribution across titer bins
# -----------------------------------------------------------------------------

def q2_phase_distribution_vs_titer(obs, fig_dir, sample,
                                    titer_col="wolbachia_titer",
                                    stage_col="cc_stage",
                                    n_titer_bins=5,
                                    ref_full=None,
                                    stage_col_ref="cyclum_stage"):
    print(f"\n-- Q2: phase distribution vs titer --")

    df = obs[[titer_col, stage_col]].copy()
    df[titer_col] = pd.to_numeric(df[titer_col], errors="coerce")
    df[stage_col] = df[stage_col].astype(str).str.strip().str.lower()
    df = df.dropna()
    print(f"   {len(df)} infected cells")

    stages = [s for s in CC_ORDER if s in df[stage_col].unique()]
    stages += [s for s in df[stage_col].unique() if s not in CC_ORDER]
    pal    = _cc_palette(stages)

    _, bin_edges = pd.qcut(df[titer_col], q=n_titer_bins,
                           retbins=True, duplicates="drop")
    actual_n_bins = len(bin_edges) - 1
    if actual_n_bins < 2:
        print("   SKIP: too few unique titer values to form bins"); return

    bin_labels = [f"Q{i+1}" for i in range(actual_n_bins)]
    df = df.copy()
    df["titer_bin"] = pd.cut(df[titer_col], bins=bin_edges,
                             labels=bin_labels, include_lowest=True)
    df = df.dropna(subset=["titer_bin"])
    actual_bins = df["titer_bin"].cat.categories.tolist()

    if actual_n_bins < n_titer_bins:
        print(f"   NOTE: requested {n_titer_bins} bins but only {actual_n_bins} "
              f"possible given duplicate titer values")

    ct_inf   = pd.crosstab(df["titer_bin"], df[stage_col]).reindex(columns=stages, fill_value=0)
    prop_inf = ct_inf.div(ct_inf.sum(axis=1), axis=0) * 100

    has_ref = ref_full is not None and stage_col_ref in ref_full.obs.columns
    ref_props = None
    ref_n     = 0
    if has_ref:
        ref_stages = (ref_full.obs[stage_col_ref]
                      .astype(str).str.strip().str.lower())
        ref_stages = ref_stages[ref_stages.isin(stages)]
        ref_n      = len(ref_stages)
        ref_vc     = ref_stages.value_counts().reindex(stages, fill_value=0)
        ref_props  = (ref_vc / ref_vc.sum() * 100).to_frame(name="Uninfected").T
        print(f"   {ref_n} uninfected reference cells added for visual comparison")
    else:
        print("   No ref_full supplied -- Uninfected bar omitted from Q2 plot")

    if has_ref:
        prop_full = pd.concat([ref_props, prop_inf], axis=0)
    else:
        prop_full = prop_inf.copy()
    prop_full.to_csv(os.path.join(fig_dir, f"q2_composition_{sample}.csv"))

    chi2_stat, _, dof, _ = chi2_contingency(ct_inf.values)
    p_chi, log10_p_chi   = chi2_p_exact(chi2_stat, dof)
    p_str = format_p(p_chi, log10_p_chi)
    print(f"   Chi-squared (infected bins): chi2={chi2_stat:.3f}, df={dof}, {p_str}  "
          f"(log10_p={log10_p_chi:.1f})")

    pd.DataFrame({
        "chi2":    [chi2_stat],
        "dof":     [dof],
        "p_value": [p_chi],
        "log10_p": [log10_p_chi],
        "n_cells": [len(df)],
    }).to_csv(os.path.join(fig_dir, f"q2_chisq_{sample}.csv"), index=False)

    bin_ranges = (df.groupby("titer_bin", observed=True)[titer_col]
                    .agg(lo="min", hi="max").reset_index())
    bin_annot  = {row["titer_bin"]: f"{row['lo']:.2f}-{row['hi']:.2f}"
                  for _, row in bin_ranges.iterrows()}

    xlabels = []
    if has_ref:
        xlabels.append(f"Uninfected\n(n={ref_n})")
    for b in actual_bins:
        xlabels.append(f"{b}\n({bin_annot.get(b,'')})")

    n_bars = len(prop_full)
    sep_x  = 0.5 if has_ref else None

    fig, axes = plt.subplots(1, 2, figsize=(max(14, n_bars * 1.8 + 4), 5))

    bottom = np.zeros(n_bars)
    for stage in stages:
        vals = prop_full[stage].values if stage in prop_full.columns else np.zeros(n_bars)
        axes[0].bar(np.arange(n_bars), vals, bottom=bottom,
                    color=pal[stage], label=stage, width=0.7)
        for xi, (v, b) in enumerate(zip(vals, bottom)):
            if v >= 7:
                axes[0].text(xi, b + v / 2, f"{v:.0f}%",
                             ha="center", va="center",
                             fontsize=7, color="white", fontweight="bold")
        bottom += vals

    if sep_x is not None:
        axes[0].axvline(sep_x, color="black", lw=1.2, ls="--", alpha=0.6)
        axes[0].text(sep_x + 0.05, 103,
                     "<- Uninfected  |  Infected titer bins ->",
                     fontsize=7, color="dimgrey", va="bottom")

    axes[0].set_xticks(np.arange(n_bars))
    axes[0].set_xticklabels(xlabels, fontsize=8)
    axes[0].set_xlabel(f"wMel titer quantile (Q1=lowest, Q{len(actual_bins)}=highest)"
                       + ("  [Uninfected shown for reference]" if has_ref else ""))
    axes[0].set_ylabel("% of cells")
    axes[0].set_ylim(0, 108)
    axes[0].set_title(f"A  CC phase composition per titer bin\n"
                      f"chi2={chi2_stat:.1f}, df={dof}, {p_str} (infected bins only)")
    axes[0].legend(title="CC stage", bbox_to_anchor=(1.01, 1),
                   loc="upper left", fontsize=9)

    sns.heatmap(prop_full.T.reindex(stages), annot=True, fmt=".1f",
                cmap="YlOrRd", linewidths=0.5, ax=axes[1],
                cbar_kws={"label": "% of cells in bin"})
    axes[1].set_xlabel("Titer bin (Uninfected shown for reference)" if has_ref
                       else "wMel titer bin")
    axes[1].set_ylabel("CC stage")
    axes[1].set_title("B  Phase proportion heatmap")
    if has_ref:
        xlbls = axes[1].get_xticklabels()
        for lbl in xlbls:
            if "Uninfected" in lbl.get_text():
                lbl.set_fontweight("bold")
                lbl.set_color("steelblue")

    plt.suptitle("Q2 -- Cell-cycle phase distribution across wMel titer levels",
                 fontweight="bold")
    plt.tight_layout()
    out = os.path.join(fig_dir, f"q2_phase_dist_vs_titer_{sample}.pdf")
    plt.savefig(out, bbox_inches="tight", dpi=150)
    plt.close()
    print(f"   -> {out}")


# -----------------------------------------------------------------------------
# Standalone publication-sized UMAP helpers (2.75 x 2.032 in)
# -----------------------------------------------------------------------------

def _save_umap_titer(obs, umap_xy, fig_dir, sample, titer_col):
    """Standalone UMAP coloured by wMel titer at publication dimensions."""
    if umap_xy is None:
        return
    tvals = obs[titer_col].values.astype(float)
    valid = ~np.isnan(tvals)

    fig, ax = plt.subplots(figsize=(UMAP_W, UMAP_H))
    sc_ = ax.scatter(umap_xy[valid, 0], umap_xy[valid, 1],
                     c=tvals[valid], cmap="viridis",
                     s=2, alpha=0.7, rasterized=True)
    ax.scatter(umap_xy[~valid, 0], umap_xy[~valid, 1],
               c="lightgrey", s=1, alpha=0.2, rasterized=True)
    plt.colorbar(sc_, ax=ax, label="wMel titer", shrink=0.8)
    ax.set_title("wMel titer", fontsize=7)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("UMAP1", fontsize=6); ax.set_ylabel("UMAP2", fontsize=6)
    plt.tight_layout()
    out = os.path.join(fig_dir, f"umap_titer_{sample}.pdf")
    plt.savefig(out, bbox_inches="tight", dpi=300)
    plt.close()
    print(f"   -> {out}")


def _save_umap_ccstage(obs, umap_xy, fig_dir, sample, stage_col, stages, pal):
    """Standalone UMAP coloured by CC stage at publication dimensions.

    The axes area is exactly UMAP_W x UMAP_H inches; legend and any
    labels extend outside via bbox_inches='tight'.
    """
    if umap_xy is None:
        return
    stage_vals = obs[stage_col].astype(str).str.lower().values
    c_map_vals = [pal.get(s, "grey") for s in stage_vals]

    fig = plt.figure(figsize=(UMAP_W, UMAP_H))
    ax  = fig.add_axes([0, 0, 1, 1])   # axes fills figure exactly
    ax.scatter(umap_xy[:, 0], umap_xy[:, 1],
               c=c_map_vals, s=0.2, alpha=0.5, linewidths=0, rasterized=True)
    for s in stages:
        ax.scatter([], [], color=pal[s], label=s, s=20)
    ax.legend(title="CC stage", fontsize=6, markerscale=2,
              title_fontsize=6, frameon=False,
              bbox_to_anchor=(1.01, 1), loc="upper left")
    ax.set_xticks([]); ax.set_yticks([])
    out = os.path.join(fig_dir, f"umap_ccstage_{sample}.pdf")
    plt.savefig(out, bbox_inches="tight", dpi=300)
    plt.close()
    print(f"   -> {out}")


# -----------------------------------------------------------------------------
# Q3 -- Titer by CC phase + polar cyclicity
# -----------------------------------------------------------------------------

def q3_titer_by_phase(obs, umap_xy, fig_dir, sample,
                      titer_col="wolbachia_titer",
                      stage_col="cc_stage",
                      pseudotime_col="cc_pseudotime",
                      n_sectors=12):
    print(f"\n-- Q3: titer by CC phase --")

    df = obs[[titer_col, stage_col]].copy()
    df[titer_col] = pd.to_numeric(df[titer_col], errors="coerce")
    df[stage_col] = df[stage_col].astype(str).str.strip().str.lower()
    df = df.dropna()
    print(f"   {len(df)} cells")

    stages = [s for s in CC_ORDER if s in df[stage_col].unique()]
    stages += [s for s in df[stage_col].unique() if s not in CC_ORDER]
    pal    = _cc_palette(stages)

    groups = [df.loc[df[stage_col] == s, titer_col].values for s in stages]

    kw_stat, _  = kruskal(*groups)
    kw_p, log10_p_kw = kw_p_exact(kw_stat, len(stages))
    p_str_kw = format_p(kw_p, log10_p_kw)
    print(f"   Kruskal-Wallis: H={kw_stat:.3f}, {p_str_kw}  "
          f"(log10_p={log10_p_kw:.1f})")

    pd.DataFrame({
        "statistic": [kw_stat],
        "p_value":   [kw_p],
        "log10_p":   [log10_p_kw],
        "n_stages":  [len(stages)],
        "n_cells":   [len(df)],
    }).to_csv(os.path.join(fig_dir, f"q3_kruskal_{sample}.csv"), index=False)

    all_ranks  = scipy.stats.rankdata(df[titer_col].values)
    n_total    = len(all_ranks)
    _, counts  = np.unique(df[titer_col].values, return_counts=True)
    tie_factor = np.sum(counts**3 - counts) / (12*(n_total - 1)) if n_total > 1 else 0
    df = df.copy()
    df["_rank"] = all_ranks

    rows = []
    for s_a, s_b in combinations(stages, 2):
        ga = df.loc[df[stage_col] == s_a, "_rank"].values
        gb = df.loc[df[stage_col] == s_b, "_rank"].values
        na, nb = len(ga), len(gb)
        if na < 2 or nb < 2: continue
        se = np.sqrt((n_total*(n_total + 1)/12 - tie_factor) * (1/na + 1/nb))
        if se == 0: continue
        z = (ga.mean() - gb.mean()) / se
        p_raw, log10_p_raw = z_p_exact(z)
        rows.append({
            "stage_A":        s_a,
            "stage_B":        s_b,
            "median_titer_A": np.median(df.loc[df[stage_col] == s_a, titer_col]),
            "median_titer_B": np.median(df.loc[df[stage_col] == s_b, titer_col]),
            "z_stat":         z,
            "p_raw":          p_raw,
            "log10_p_raw":    log10_p_raw,
            "n_A":            na,
            "n_B":            nb,
        })

    dunn_df = pd.DataFrame(rows)
    if len(dunn_df):
        p_adj, log10_p_adj = bh_adjust_log10(dunn_df["log10_p_raw"].values)
        dunn_df["p_adj_BH"]    = p_adj
        dunn_df["log10_p_adj"] = log10_p_adj
        dunn_df["significant"] = dunn_df["log10_p_adj"] < np.log10(0.05)
        dunn_df = dunn_df.sort_values("log10_p_adj")
        dunn_df.to_csv(os.path.join(fig_dir, f"q3_dunn_{sample}.csv"), index=False)

        print(f"   Dunn: {dunn_df['significant'].sum()}/{len(dunn_df)} pairs significant")
        disp = dunn_df[["stage_A", "stage_B", "median_titer_A", "median_titer_B",
                         "log10_p_adj", "significant"]].copy()
        disp["p_adj_str"] = [format_p(pa, lp)
                             for pa, lp in zip(dunn_df["p_adj_BH"], dunn_df["log10_p_adj"])]
        print(disp[["stage_A", "stage_B", "median_titer_A", "median_titer_B",
                    "p_adj_str", "significant"]].to_string(index=False))

    stage_med = df.groupby(stage_col)[titer_col].median().to_dict()
    (df.groupby(stage_col)[titer_col]
       .agg(n_cells="count", median="median", mean="mean", std="std",
            q25=lambda x: x.quantile(0.25),
            q75=lambda x: x.quantile(0.75))
       .reindex(stages).reset_index()
       .rename(columns={stage_col: "stage"})
       .to_csv(os.path.join(fig_dir, f"q3_stage_summary_{sample}.csv"), index=False))

    has_pt = pseudotime_col in obs.columns
    df_pol = None
    if has_pt:
        df_pol = obs[[titer_col, pseudotime_col, stage_col]].copy()
        df_pol[titer_col]      = pd.to_numeric(df_pol[titer_col],      errors="coerce")
        df_pol[pseudotime_col] = pd.to_numeric(df_pol[pseudotime_col], errors="coerce")
        df_pol = df_pol.dropna()

    has_umap = umap_xy is not None

    # -- Combined overview figure (violin + UMAP panels + polar) --------------
    fig = plt.figure(figsize=(22, 5))
    gs  = fig.add_gridspec(1, 4, wspace=0.45)
    ax_vln = fig.add_subplot(gs[0])
    ax_ut  = fig.add_subplot(gs[1])
    ax_uc  = fig.add_subplot(gs[2])
    ax_pol = fig.add_subplot(gs[3], projection="polar")

    sns.violinplot(data=df, x=stage_col, y=titer_col, order=stages,
                   palette=pal, inner=None, linewidth=0.8, ax=ax_vln, cut=0)
    sns.stripplot(data=df, x=stage_col, y=titer_col, order=stages,
                  palette=pal, size=1.5, alpha=0.35, jitter=True, ax=ax_vln)
    for i, s in enumerate(stages):
        ax_vln.scatter(i, stage_med[s], color="white", s=35, zorder=5,
                       edgecolors="black", linewidths=0.8)
    if len(dunn_df):
        sig    = dunn_df[dunn_df["significant"]]
        y_max  = df[titer_col].quantile(0.99)
        y_step = (df[titer_col].quantile(0.99) - df[titer_col].quantile(0.01)) * 0.08
        for k_idx, (_, row) in enumerate(sig.iterrows()):
            if row["stage_A"] not in stages or row["stage_B"] not in stages: continue
            xi = stages.index(row["stage_A"])
            xj = stages.index(row["stage_B"])
            y  = y_max + y_step * (k_idx + 1)
            ax_vln.plot([xi, xj], [y, y], color="black", lw=1)
            ax_vln.text((xi + xj) / 2, y + y_step * 0.15, "*",
                        ha="center", va="bottom", fontsize=10)
    ax_vln.set_xlabel("Cell-cycle stage")
    ax_vln.set_ylabel("wMel titer")
    ax_vln.set_title(f"A  Titer by CC stage\nKW {p_str_kw}")

    if has_umap:
        tvals = obs[titer_col].values.astype(float)
        valid = ~np.isnan(tvals)
        sc_   = ax_ut.scatter(umap_xy[valid, 0], umap_xy[valid, 1],
                              c=tvals[valid], cmap="viridis",
                              s=2, alpha=0.7, rasterized=True)
        ax_ut.scatter(umap_xy[~valid, 0], umap_xy[~valid, 1],
                      c="lightgrey", s=1, alpha=0.2, rasterized=True)
        plt.colorbar(sc_, ax=ax_ut, label="wMel titer", shrink=0.8)
    ax_ut.set_title("B  UMAP -- wMel titer")
    ax_ut.set_xticks([]); ax_ut.set_yticks([])

    if has_umap:
        stage_vals = obs[stage_col].astype(str).str.lower().values
        c_map_vals = [pal.get(s, "grey") for s in stage_vals]
        ax_uc.scatter(umap_xy[:, 0], umap_xy[:, 1],
                      c=c_map_vals, s=2, alpha=0.7, rasterized=True)
        for s in stages:
            ax_uc.scatter([], [], color=pal[s], label=s, s=20)
        ax_uc.legend(title="CC stage", fontsize=8, markerscale=2,
                     bbox_to_anchor=(1.01, 1), loc="upper left")
    ax_uc.set_title("C  UMAP -- CC stage (KNN transferred)")
    ax_uc.set_xticks([]); ax_uc.set_yticks([])

    if has_pt and df_pol is not None and len(df_pol) >= 20:
        edges   = np.linspace(0, 2*np.pi, n_sectors + 1)
        centers = (edges[:-1] + edges[1:]) / 2
        df_pol  = df_pol.copy()
        df_pol["_sector"] = pd.cut(df_pol[pseudotime_col], bins=edges,
                                   labels=range(n_sectors), include_lowest=True)
        sec_med = (df_pol.groupby("_sector", observed=True)[titer_col]
                         .median()
                         .reindex(range(n_sectors))
                         .fillna(0).values)
        vmin, vmax = sec_med.min(), sec_med.max()
        norm_vals  = (sec_med - vmin) / (vmax - vmin + 1e-9)
        width      = 2*np.pi / n_sectors
        for theta, r, nv in zip(centers, sec_med, norm_vals):
            ax_pol.bar(theta, r, width=width * 0.85, bottom=0,
                       color=plt.cm.coolwarm(nv), alpha=0.85)
        r_annot = sec_med.max() * 1.2
        for name, th0, th1 in [("g0/g1", 0, np.pi/2),
                                ("s",     np.pi/2, 3*np.pi/2),
                                ("g2/m",  3*np.pi/2, 2*np.pi)]:
            if name in pal:
                ax_pol.text((th0 + th1) / 2, r_annot, name,
                            ha="center", va="center",
                            fontsize=7, color=pal[name], fontweight="bold")
        ax_pol.set_theta_zero_location("N")
        ax_pol.set_theta_direction(-1)
        ax_pol.set_xticks(np.linspace(0, 2*np.pi, 5)[:-1])
        ax_pol.set_xticklabels(["0", "pi/2", "pi", "3pi/2"], fontsize=8)
        ax_pol.set_title(f"D  Titer cyclicity\n(n={n_sectors} sectors)", pad=18, fontsize=9)
    else:
        ax_pol.set_title("D  Polar: pseudotime not available", pad=15, fontsize=9)

    plt.suptitle("Q3 -- wMel titer across cell-cycle phases", fontweight="bold", y=1.02)
    out = os.path.join(fig_dir, f"q3_titer_by_phase_{sample}.pdf")
    plt.savefig(out, bbox_inches="tight", dpi=150)
    plt.close()
    print(f"   -> {out}")

    # -- Standalone publication-sized UMAP panels (2.75 x 2.032 in) ----------
    _save_umap_titer(obs, umap_xy, fig_dir, sample, titer_col)
    _save_umap_ccstage(obs, umap_xy, fig_dir, sample, stage_col, stages, pal)


# -----------------------------------------------------------------------------
# Sanitise + save
# -----------------------------------------------------------------------------

def _sanitize_obs(adata):
    for col in adata.obs.columns:
        if adata.obs[col].dtype == object or str(adata.obs[col].dtype) == "category":
            adata.obs[col] = adata.obs[col].astype(str).replace("nan", "NA")
    return adata


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def run(ref_path, query_paths, out_path, fig_dir, sample,
        harmony_vars, ref_condition, titer_col, stage_col,
        pseudotime_col, n_pcs, k_neighbors, n_titer_bins, n_sectors,
        infected_only, sample_label_col):

    os.makedirs(fig_dir, exist_ok=True)
    sc.settings.figdir = fig_dir

    expanded = []
    for p in query_paths:
        matched = sorted(glob.glob(p))
        if matched:
            expanded.extend(matched)
        else:
            expanded.append(p)
    query_paths = expanded
    print(f"Query files ({len(query_paths)}):")
    for p in query_paths:
        print(f"  {p}")

    ref_raw, ref_full = load_reference(ref_path, stage_col, pseudotime_col)
    ref_full, marker_df = identify_reference_markers(
        ref_full, fig_dir, sample, groupby="leiden", n_top=25
    )

    print(f"\n-- Loading query files --")
    query_raw = load_query_files(
        query_paths, ref_condition=ref_condition,
        titer_col=titer_col, infected_only=infected_only,
    )

    combined, ref_mask = joint_preprocess_and_harmony(
        ref_raw, query_raw, ref_full,
        harmony_vars=harmony_vars, n_pcs=n_pcs,
    )

    query, combined = knn_label_transfer(
        combined, ref_mask, ref_full, k=k_neighbors,
    )

    query, cluster_map = assign_cc_from_reference(
        query, ref_full,
        stage_col=stage_col, pseudotime_col=pseudotime_col,
    )
    cluster_map.to_csv(
        os.path.join(fig_dir, f"cluster_cc_map_{sample}.csv"), index=False
    )

    sc.pl.umap(combined,
               color=["dataset", "leiden"] if "leiden" in combined.obs.columns
               else ["dataset"],
               save=f"_{sample}_combined_dataset.pdf")

    obs  = query.obs.copy()
    umap = query.obsm.get("X_umap", None)

    q4_cc_distribution_by_sample(
        query_obs=obs,
        ref_full=ref_full,
        fig_dir=fig_dir,
        sample=sample,
        stage_col=stage_col,
        query_stage_col="cc_stage",
        titer_col=titer_col,
        sample_label_col=sample_label_col,
        n_titer_bins=n_titer_bins,
    )

    q1_titer_vs_pseudotime(
        obs, fig_dir, sample,
        titer_col=titer_col, pseudotime_col="cc_pseudotime", stage_col="cc_stage",
    )
    q2_phase_distribution_vs_titer(
        obs, fig_dir, sample,
        titer_col=titer_col, stage_col="cc_stage", n_titer_bins=n_titer_bins,
        ref_full=ref_full, stage_col_ref=stage_col,
    )
    q3_titer_by_phase(
        obs, umap, fig_dir, sample,
        titer_col=titer_col, stage_col="cc_stage",
        pseudotime_col="cc_pseudotime", n_sectors=n_sectors,
    )

    query    = _sanitize_obs(query)
    combined = _sanitize_obs(combined)

    query.write(out_path)
    combined.write(out_path.replace(".h5ad", "_combined.h5ad"))
    print(f"\nSaved query    -> {out_path}")
    print(f"Saved combined -> {out_path.replace('.h5ad', '_combined.h5ad')}")
    print(f"Figures        -> {fig_dir}/")


def main():
    parser = argparse.ArgumentParser(
        description="Project infected cells onto CC reference via Harmony + KNN label transfer"
    )
    parser.add_argument("--ref",            required=True)
    parser.add_argument("--query",          required=True, nargs="+",
                        help="One or more h5ad files (globs ok)")
    parser.add_argument("--out_path",       default="results/integrated/titer_cellcycle.h5ad")
    parser.add_argument("--fig_dir",        default="figures/titer_vs_cellcycle")
    parser.add_argument("--sample",         default="wolbachia_infection")
    parser.add_argument("--harmony_vars",   nargs="+", default=["method", "replicate", "dataset"],
                        help="obs columns to correct in Harmony")
    parser.add_argument("--ref_condition",  nargs="+", default=["JW18DOX"],
                        help="cell_line values that are the uninfected reference")
    parser.add_argument("--titer_col",      default="wolbachia_titer")
    parser.add_argument("--stage_col",      default="cyclum_stage")
    parser.add_argument("--pseudotime_col", default="cyclum_pseudotime")
    parser.add_argument("--n_pcs",          type=int, default=30)
    parser.add_argument("--k_neighbors",    type=int, default=15,
                        help="k for KNN label transfer")
    parser.add_argument("--n_titer_bins",   type=int, default=5)
    parser.add_argument("--n_sectors",      type=int, default=12)
    parser.add_argument("--sample_label_col", default=None,
                        help="obs column to use as sample label in Q4 table "
                             "(default: auto-infer from cell_line + treatment + timepoint)")
    parser.add_argument("--all_cells",      action="store_true",
                        help="Use all cells in query files (default: infected only)")
    args = parser.parse_args()

    run(
        ref_path=args.ref,
        query_paths=args.query,
        out_path=args.out_path,
        fig_dir=args.fig_dir,
        sample=args.sample,
        harmony_vars=args.harmony_vars,
        ref_condition=args.ref_condition,
        titer_col=args.titer_col,
        stage_col=args.stage_col,
        pseudotime_col=args.pseudotime_col,
        n_pcs=args.n_pcs,
        k_neighbors=args.k_neighbors,
        n_titer_bins=args.n_titer_bins,
        n_sectors=args.n_sectors,
        infected_only=not args.all_cells,
        sample_label_col=args.sample_label_col,
    )


if __name__ == "__main__":
    main()