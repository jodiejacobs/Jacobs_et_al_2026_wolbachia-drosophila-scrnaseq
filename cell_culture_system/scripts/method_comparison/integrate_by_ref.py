"""
titer_vs_cellcycle.py
=====================
Tests the relationship between wMel titer and cell cycle in infected Drosophila JW18 cells.

Strategy
--------
1. Load reference (uninfected, has cyclum_stage/pseudotime/leiden).
   Extract log-normalised counts from reference.raw.X.
2. Load each query h5ad, extract raw counts from .raw.X, subset to
   infected cells (treatment != Ctrl).
3. Concatenate reference + query on shared genes.
4. Jointly normalise (already done in ref, redo from raw for consistency)
   → restrict to reference HVGs → scale → PCA → Harmony.
5. KNN label transfer: for each query cell find k nearest reference
   neighbours in Harmony PCA space, majority-vote their leiden label.
6. Map leiden → cc_stage / cc_pseudotime via reference majority-vote.
7. Run Q1/Q2/Q3 analyses.

Run with:
    mamba activate scanpy
    python scripts/method_comparison/titer_vs_cellcycle.py \
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
from scipy.stats import kruskal, spearmanr, chi2_contingency
from itertools import combinations
from statsmodels.stats.multitest import multipletests
from statsmodels.nonparametric.smoothers_lowess import lowess
import anndata as ad
import scanpy as sc
import harmonypy as hm


# ─────────────────────────────────────────────────────────────────────────────
# Colours
# ─────────────────────────────────────────────────────────────────────────────

CC_ORDER  = ["g0/g1", "s", "g2/m"]
CC_COLORS = {"g0/g1": "#4C72B0", "s": "#DD8452", "g2/m": "#55A868"}

def _cc_palette(stages):
    cmap = matplotlib.colormaps["tab10"]
    return {s: CC_COLORS.get(s, cmap(i % 10)) for i, s in enumerate(stages)}


# ─────────────────────────────────────────────────────────────────────────────
# Step 1 — Load reference, extract log-normalised counts from .raw
# ─────────────────────────────────────────────────────────────────────────────

def load_reference(ref_path, stage_col, pseudotime_col):
    """Load reference and extract raw counts for joint re-processing.

    The reference .raw.X contains raw counts (post-QC, pre-normalisation).
    We reconstruct a clean AnnData from .raw so the joint pipeline starts
    from the same point as query samples.

    Returns
    -------
    ref_raw : AnnData
        Raw counts in .X, all reference obs columns preserved.
    ref_full : AnnData
        Original reference (needed for leiden/CC labels and Harmony PCA).
    """
    print(f"\n── Loading reference: {ref_path} ──")
    ref_full = sc.read_h5ad(ref_path)
    print(f"   {ref_full.n_obs} cells × {ref_full.n_vars} genes")
    print(f"   Leiden clusters : {ref_full.obs['leiden'].nunique()}")

    for col in [stage_col, pseudotime_col, "leiden"]:
        if col not in ref_full.obs.columns:
            raise ValueError(f"'{col}' not in reference obs. "
                             f"Available: {list(ref_full.obs.columns)}")

    print(f"   CC stage distribution (reference):")
    print(ref_full.obs[stage_col].value_counts().to_string())

    if ref_full.raw is None:
        raise ValueError(
            "Reference has no .raw — re-run filter.py with the updated "
            "analyze_filtered_adata() that saves adata.raw = adata before normalisation, "
            "then re-run the integration pipeline that built this reference."
        )

    # Reconstruct AnnData from .raw (raw counts, full gene universe)
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
    print(f"   Reference raw counts: {ref_raw.n_obs} cells × {ref_raw.n_vars} genes")

    return ref_raw, ref_full


# ─────────────────────────────────────────────────────────────────────────────
# Step 2 — Load query files, extract raw counts, subset to infected cells
# ─────────────────────────────────────────────────────────────────────────────

def load_query_files(query_paths, ref_condition, titer_col,
                     treatment_col="treatment", infected_only=True):
    """Load and concatenate all query h5ad files.

    Each file must have .raw.X containing raw counts. Infected cells are
    identified as treatment != Ctrl (or all cells if infected_only=False).

    Parameters
    ----------
    query_paths : list[str]
        Paths to filtered h5ad files (may include globs, already expanded).
    ref_condition : list[str]
        cell_line values that define uninfected reference cells — these are
        excluded from the query even if present in the query files.
    titer_col : str
        obs column for wMel titer; must be present in each query file.
    infected_only : bool
        If True, subset to treatment != Ctrl cells only.

    Returns
    -------
    query_raw : AnnData
        Concatenated raw counts from all query files (infected cells only).
    """
    adatas = []

    for path in query_paths:
        print(f"\n   Loading query: {path}")
        adata = sc.read_h5ad(path)
        print(f"   {adata.n_obs} cells × {adata.n_vars} genes")

        if adata.raw is None:
            raise ValueError(
                f"{path} has no .raw — re-run filter.py with the updated "
                "analyze_filtered_adata() that saves adata.raw = adata "
                "before normalisation."
            )

        # Extract raw counts
        raw_X = adata.raw.X
        if scipy.sparse.issparse(raw_X):
            raw_X = raw_X.toarray()
        raw_X = raw_X.astype(np.float32)

        a = ad.AnnData(
            X=scipy.sparse.csr_matrix(raw_X),
            obs=adata.obs.copy(),
            var=adata.raw.var.copy(),
        )

        # Infer metadata from filename if not in obs
        basename = os.path.basename(path).replace(".h5ad", "")
        if "method" not in a.obs.columns:
            method = "pipseq" if "pipseq" in basename.lower() else "10x"
            a.obs["method"] = method
            print(f"   Inferred method='{method}' from filename")
        if "replicate" not in a.obs.columns:
            # e.g. JW18wMel-SV3-1_pipseq → replicate = 1
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

        # Subset to infected cells
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
            print(f"   WARNING: no cells remain after filtering — skipping {path}")
            continue

        if titer_col not in a.obs.columns:
            print(f"   WARNING: '{titer_col}' not in obs — titer analyses will skip this file")

        a.obs["dataset"]   = "query"
        a.obs["source_file"] = basename
        a.obs_names = [f"{basename}_{bc}" for bc in a.obs_names]

        adatas.append(a)
        print(f"   Kept {a.n_obs} cells")

    if not adatas:
        raise ValueError("No query cells loaded — check paths and --ref_condition")

    print(f"\n── Concatenating {len(adatas)} query files ──")
    query_raw = ad.concat(adatas, join="outer", index_unique=None)
    query_raw.obs_names_make_unique()
    print(f"   Total query cells: {query_raw.n_obs}")

    return query_raw


# ─────────────────────────────────────────────────────────────────────────────
# Step 3 — Joint normalisation → HVG restriction → scale → PCA → Harmony
# ─────────────────────────────────────────────────────────────────────────────

def joint_preprocess_and_harmony(ref_raw, query_raw, ref_full,
                                  harmony_vars, n_pcs=30):
    """Jointly normalise and batch-correct reference + query.

    Processing order
    ----------------
    1. Concatenate ref_raw + query_raw on shared genes (outer join,
       missing filled with 0).
    2. Normalise (1e4 per cell) + log1p.
    3. Restrict to reference HVGs (already selected in ref_full.var).
    4. Scale (max_value=10).
    5. PCA (n_pcs components, using HVG subset).
    6. Harmony batch correction on harmony_vars.

    Returns
    -------
    combined : AnnData
        Full combined object with X_pca_harmony in obsm.
    ref_mask : np.ndarray[bool]
        Boolean mask indicating which rows are reference cells.
    """
    print(f"\n── Joint preprocessing ──")

    # ── 1. Concatenate ────────────────────────────────────────────────────────
    combined = ad.concat(
        [ref_raw, query_raw],
        join="outer",      # fills missing genes with NaN → 0 below
        index_unique=None,
        label="dataset",
        keys=["reference", "query"],
    )
    combined.obs_names_make_unique()

    # Fill NaN from outer join
    if scipy.sparse.issparse(combined.X):
        combined.X = combined.X.toarray()
    combined.X = np.nan_to_num(combined.X.astype(np.float32), nan=0.0)
    combined.X = scipy.sparse.csr_matrix(combined.X)

    ref_mask = combined.obs["dataset"] == "reference"
    print(f"   Combined: {combined.n_obs} cells × {combined.n_vars} genes")
    print(f"   Reference: {ref_mask.sum()}  Query: {(~ref_mask).sum()}")

    # ── 2. Normalise + log1p ──────────────────────────────────────────────────
    print("   Normalising (1e4 per cell) + log1p …")
    sc.pp.normalize_total(combined, target_sum=1e4)
    sc.pp.log1p(combined)

    # ── 3. Restrict to reference HVGs ─────────────────────────────────────────
    ref_hvgs = ref_full.var_names[ref_full.var["highly_variable"]].tolist()
    hvgs_present = [g for g in ref_hvgs if g in combined.var_names]
    print(f"   Reference HVGs: {len(ref_hvgs)} total, "
          f"{len(hvgs_present)} present in combined dataset")

    if len(hvgs_present) < 100:
        raise ValueError(
            f"Only {len(hvgs_present)} reference HVGs found in combined dataset. "
            "Check that query files share the same gene universe as the reference."
        )

    combined.var["highly_variable"] = combined.var_names.isin(hvgs_present)

    # ── 4. Scale ──────────────────────────────────────────────────────────────
    print("   Scaling (max_value=10) …")
    sc.pp.scale(combined, max_value=10)

    # ── 5. PCA ────────────────────────────────────────────────────────────────
    print(f"   PCA ({n_pcs} components) …")
    sc.tl.pca(combined, n_comps=n_pcs, use_highly_variable=True, svd_solver="arpack")

    # ── 6. Harmony ────────────────────────────────────────────────────────────
    # Validate harmony_vars exist
    missing_vars = [v for v in harmony_vars if v not in combined.obs.columns]
    if missing_vars:
        print(f"   WARNING: harmony_vars {missing_vars} not in obs — dropping them")
        harmony_vars = [v for v in harmony_vars if v in combined.obs.columns]
    if not harmony_vars:
        raise ValueError("No valid harmony_vars remain — cannot run Harmony")

    # Fill any NaN in harmony vars (e.g. replicate missing for some cells)
    for v in harmony_vars:
        combined.obs[v] = combined.obs[v].astype(str).fillna("unknown")

    print(f"   Harmony correction on: {harmony_vars} …")
    pca_matrix = combined.obsm["X_pca"]
    meta       = combined.obs[harmony_vars].copy()
    ho = hm.run_harmony(pca_matrix, meta, harmony_vars,
                         max_iter_harmony=30, random_state=42)
    combined.obsm["X_pca_harmony"] = ho.Z_corr.T

    print("   Harmony complete")
    return combined, ref_mask


# ─────────────────────────────────────────────────────────────────────────────
# Step 4 — KNN label transfer
# ─────────────────────────────────────────────────────────────────────────────

def knn_label_transfer(combined, ref_mask, ref_full, k=15):
    """Transfer leiden labels from reference to query via KNN in Harmony PCA space.

    For each query cell, find k nearest reference neighbours in
    X_pca_harmony, then assign the majority-vote leiden label.

    Parameters
    ----------
    combined : AnnData
        Output of joint_preprocess_and_harmony.
    ref_mask : np.ndarray[bool]
        Boolean mask for reference cells in combined.
    ref_full : AnnData
        Original reference (provides leiden labels).
    k : int
        Number of nearest neighbours for voting.

    Returns
    -------
    query : AnnData
        Query cells with obs["leiden_ref"] added.
    """
    from sklearn.neighbors import NearestNeighbors

    print(f"\n── KNN label transfer (k={k}) ──")

    ref_pca   = combined.obsm["X_pca_harmony"][ref_mask]
    query_pca = combined.obsm["X_pca_harmony"][~ref_mask]

    # leiden labels aligned to combined reference rows
    ref_leiden = ref_full.obs["leiden"].astype(str).values

    nbrs = NearestNeighbors(n_neighbors=k, metric="euclidean", n_jobs=-1)
    nbrs.fit(ref_pca)
    distances, indices = nbrs.kneighbors(query_pca)

    # Majority vote
    transferred = []
    confidence  = []
    for row_idx in indices:
        neighbour_labels = ref_leiden[row_idx]
        counts = pd.Series(neighbour_labels).value_counts()
        winner = counts.index[0]
        transferred.append(winner)
        confidence.append(counts.iloc[0] / k)

    # Build query AnnData from combined
    query = combined[~ref_mask].copy()
    query.obs["leiden_ref"]        = transferred
    query.obs["knn_confidence"]    = confidence
    query.obsm["X_pca_harmony"]    = query_pca

    # Also store UMAP projection (computed on combined, subset here)
    print("   Computing UMAP on combined Harmony embedding …")
    sc.pp.neighbors(combined, use_rep="X_pca_harmony", n_neighbors=30)
    sc.tl.umap(combined)
    query.obsm["X_umap"] = combined.obsm["X_umap"][~ref_mask]

    print(f"   Transferred leiden_ref distribution:")
    print(pd.Series(transferred).value_counts().sort_index().to_string())
    print(f"   Mean KNN confidence: {np.mean(confidence):.3f}")

    low_conf = np.mean(np.array(confidence) < 0.5)
    if low_conf > 0.2:
        print(f"   WARNING: {low_conf*100:.1f}% of cells have KNN confidence < 0.5 "
              f"— batch correction may be insufficient")

    return query, combined


# ─────────────────────────────────────────────────────────────────────────────
# Step 5 — Map leiden_ref → CC stage + pseudotime
# ─────────────────────────────────────────────────────────────────────────────

def assign_cc_from_reference(query, ref_full,
                              stage_col="cyclum_stage",
                              pseudotime_col="cyclum_pseudotime",
                              leiden_col="leiden"):
    """Assign CC stage and pseudotime to query via reference cluster majority vote."""

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

    print("\n── Cluster → CC stage mapping ──")
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

    print(f"\n   Query CC stage distribution:")
    print(query.obs["cc_stage"].value_counts().to_string())

    return query, stage_map


# ─────────────────────────────────────────────────────────────────────────────
# Q1 — titer vs cell-cycle pseudotime
# ─────────────────────────────────────────────────────────────────────────────

def q1_titer_vs_pseudotime(obs, fig_dir, sample,
                            titer_col="wolbachia_titer",
                            pseudotime_col="cc_pseudotime",
                            stage_col="cc_stage",
                            n_bins=8):
    print(f"\n── Q1: titer vs CC pseudotime ──")

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

    rho, p_val = spearmanr(pt, tit)
    print(f"   Spearman rho={rho:.3f}, p={p_val:.3e}")

    pd.DataFrame({
        "spearman_rho": [rho], "p_value": [p_val], "n_cells": [len(df)],
        "note": ["cc_pseudotime is cluster-median from reference"],
    }).to_csv(os.path.join(fig_dir, f"q1_spearman_{sample}.csv"), index=False)

    # LOWESS
    do_lowess = n_unique_pt >= 5
    if do_lowess:
        order    = np.argsort(pt)
        smoothed = lowess(tit[order], pt[order], frac=0.3, return_sorted=True)

    # Binned summary — cap bins to unique pseudotime values
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
    pi_labels = ["0", "π/2", "π", "3π/2", "2π"]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    for stage in stages:
        m = df[stage_col] == stage
        axes[0].scatter(df.loc[m, pseudotime_col], df.loc[m, titer_col],
                        c=[pal[stage]], s=3, alpha=0.4, label=stage, rasterized=True)
    axes[0].set_xlabel("CC pseudotime (cluster-median, radians)")
    axes[0].set_ylabel("wMel titer")
    axes[0].set_title(f"A  Titer vs pseudotime\nSpearman ρ={rho:.3f}, p={p_val:.2e}")
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
    axes[2].set_title(f"C  Median titer ± IQR")

    plt.suptitle("Q1 — wMel titer vs cell-cycle pseudotime", fontweight="bold")
    plt.tight_layout()
    out = os.path.join(fig_dir, f"q1_titer_vs_pseudotime_{sample}.pdf")
    plt.savefig(out, bbox_inches="tight", dpi=150)
    plt.close()
    print(f"   → {out}")


# ─────────────────────────────────────────────────────────────────────────────
# Q2 — CC phase distribution across titer bins
# ─────────────────────────────────────────────────────────────────────────────

def q2_phase_distribution_vs_titer(obs, fig_dir, sample,
                                    titer_col="wolbachia_titer",
                                    stage_col="cc_stage",
                                    n_titer_bins=5):
    print(f"\n── Q2: phase distribution vs titer ──")

    df = obs[[titer_col, stage_col]].copy()
    df[titer_col] = pd.to_numeric(df[titer_col], errors="coerce")
    df[stage_col] = df[stage_col].astype(str).str.strip().str.lower()
    df = df.dropna()
    print(f"   {len(df)} cells")

    stages = [s for s in CC_ORDER if s in df[stage_col].unique()]
    stages += [s for s in df[stage_col].unique() if s not in CC_ORDER]
    pal    = _cc_palette(stages)

    bin_labels = [f"Q{i+1}" for i in range(n_titer_bins)]
    df = df.copy()
    df["titer_bin"] = pd.qcut(df[titer_col], q=n_titer_bins,
                               labels=bin_labels, duplicates="drop")
    df = df.dropna(subset=["titer_bin"])
    actual_bins = df["titer_bin"].cat.categories.tolist()

    ct   = pd.crosstab(df["titer_bin"], df[stage_col]).reindex(columns=stages, fill_value=0)
    prop = ct.div(ct.sum(axis=1), axis=0) * 100
    prop.to_csv(os.path.join(fig_dir, f"q2_composition_{sample}.csv"))

    chi2, p_chi, dof, _ = chi2_contingency(ct.values)
    print(f"   Chi-squared: χ²={chi2:.3f}, df={dof}, p={p_chi:.3e}")
    pd.DataFrame({
        "chi2": [chi2], "dof": [dof], "p_value": [p_chi], "n_cells": [len(df)],
    }).to_csv(os.path.join(fig_dir, f"q2_chisq_{sample}.csv"), index=False)

    bin_ranges = (df.groupby("titer_bin", observed=True)[titer_col]
                    .agg(lo="min", hi="max").reset_index())
    bin_annot  = {row["titer_bin"]: f"{row['lo']:.2f}–{row['hi']:.2f}"
                  for _, row in bin_ranges.iterrows()}

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    bottom = np.zeros(len(prop))
    for stage in stages:
        vals = prop[stage].values
        axes[0].bar(np.arange(len(prop)), vals, bottom=bottom,
                    color=pal[stage], label=stage, width=0.7)
        bottom += vals
    axes[0].set_xticks(np.arange(len(prop)))
    axes[0].set_xticklabels(
        [f"{b}\n({bin_annot.get(b,'')})" for b in prop.index], fontsize=8)
    axes[0].set_xlabel(f"wMel titer quantile (Q1=lowest, Q{len(actual_bins)}=highest)")
    axes[0].set_ylabel("% of cells")
    axes[0].set_title(f"A  CC phase composition per titer bin\n"
                       f"χ²={chi2:.1f}, df={dof}, p={p_chi:.2e}")
    axes[0].legend(title="CC stage", bbox_to_anchor=(1.01, 1),
                   loc="upper left", fontsize=9)

    sns.heatmap(prop.T.reindex(stages), annot=True, fmt=".1f",
                cmap="YlOrRd", linewidths=0.5, ax=axes[1],
                cbar_kws={"label": "% of cells in titer bin"})
    axes[1].set_xlabel("wMel titer bin")
    axes[1].set_ylabel("CC stage")
    axes[1].set_title("B  Phase proportion heatmap")

    plt.suptitle("Q2 — Cell-cycle phase distribution across wMel titer levels",
                 fontweight="bold")
    plt.tight_layout()
    out = os.path.join(fig_dir, f"q2_phase_dist_vs_titer_{sample}.pdf")
    plt.savefig(out, bbox_inches="tight", dpi=150)
    plt.close()
    print(f"   → {out}")


# ─────────────────────────────────────────────────────────────────────────────
# Q3 — Titer by CC phase + polar cyclicity
# ─────────────────────────────────────────────────────────────────────────────

def q3_titer_by_phase(obs, umap_xy, fig_dir, sample,
                       titer_col="wolbachia_titer",
                       stage_col="cc_stage",
                       pseudotime_col="cc_pseudotime",
                       n_sectors=12):
    print(f"\n── Q3: titer by CC phase ──")

    df = obs[[titer_col, stage_col]].copy()
    df[titer_col] = pd.to_numeric(df[titer_col], errors="coerce")
    df[stage_col] = df[stage_col].astype(str).str.strip().str.lower()
    df = df.dropna()
    print(f"   {len(df)} cells")

    stages = [s for s in CC_ORDER if s in df[stage_col].unique()]
    stages += [s for s in df[stage_col].unique() if s not in CC_ORDER]
    pal    = _cc_palette(stages)

    groups  = [df.loc[df[stage_col] == s, titer_col].values for s in stages]
    kw_stat, kw_p = kruskal(*groups)
    print(f"   Kruskal-Wallis: H={kw_stat:.3f}, p={kw_p:.3e}")

    pd.DataFrame({
        "statistic": [kw_stat], "p_value": [kw_p],
        "n_stages": [len(stages)], "n_cells": [len(df)],
    }).to_csv(os.path.join(fig_dir, f"q3_kruskal_{sample}.csv"), index=False)

    all_ranks  = scipy.stats.rankdata(df[titer_col].values)
    n_total    = len(all_ranks)
    _, counts  = np.unique(df[titer_col].values, return_counts=True)
    tie_factor = np.sum(counts**3 - counts) / (12*(n_total-1)) if n_total > 1 else 0
    df = df.copy()
    df["_rank"] = all_ranks

    rows = []
    for s_a, s_b in combinations(stages, 2):
        ga = df.loc[df[stage_col] == s_a, "_rank"].values
        gb = df.loc[df[stage_col] == s_b, "_rank"].values
        na, nb = len(ga), len(gb)
        if na < 2 or nb < 2: continue
        se = np.sqrt((n_total*(n_total+1)/12 - tie_factor) * (1/na + 1/nb))
        if se == 0: continue
        z = (ga.mean() - gb.mean()) / se
        rows.append({
            "stage_A": s_a, "stage_B": s_b,
            "median_titer_A": np.median(df.loc[df[stage_col]==s_a, titer_col]),
            "median_titer_B": np.median(df.loc[df[stage_col]==s_b, titer_col]),
            "z_stat": z, "p_raw": 2*scipy.stats.norm.sf(abs(z)),
            "n_A": na, "n_B": nb,
        })

    dunn_df = pd.DataFrame(rows)
    if len(dunn_df):
        _, p_adj, _, _ = multipletests(dunn_df["p_raw"], method="fdr_bh")
        dunn_df["p_adj_BH"]    = p_adj
        dunn_df["significant"] = dunn_df["p_adj_BH"] < 0.05
        dunn_df = dunn_df.sort_values("p_adj_BH")
        dunn_df.to_csv(os.path.join(fig_dir, f"q3_dunn_{sample}.csv"), index=False)
        print(f"   Dunn: {dunn_df['significant'].sum()}/{len(dunn_df)} pairs significant")
        print(dunn_df[["stage_A","stage_B","median_titer_A",
                        "median_titer_B","p_adj_BH","significant"]].to_string(index=False))

    stage_med = df.groupby(stage_col)[titer_col].median().to_dict()
    (df.groupby(stage_col)[titer_col]
       .agg(n_cells="count", median="median", mean="mean", std="std",
            q25=lambda x: x.quantile(0.25), q75=lambda x: x.quantile(0.75))
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
    fig = plt.figure(figsize=(22, 5))
    gs  = fig.add_gridspec(1, 4, wspace=0.45)
    ax_vln = fig.add_subplot(gs[0])
    ax_ut  = fig.add_subplot(gs[1])
    ax_uc  = fig.add_subplot(gs[2])
    ax_pol = fig.add_subplot(gs[3], projection="polar")

    # A: violin
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
            ax_vln.text((xi+xj)/2, y+y_step*0.15, "*",
                        ha="center", va="bottom", fontsize=10)
    ax_vln.set_xlabel("Cell-cycle stage")
    ax_vln.set_ylabel("wMel titer")
    ax_vln.set_title(f"A  Titer by CC stage\nKW p={kw_p:.2e}")

    # B: UMAP titer
    if has_umap:
        tvals = obs[titer_col].values.astype(float)
        valid = ~np.isnan(tvals)
        sc_   = ax_ut.scatter(umap_xy[valid,0], umap_xy[valid,1],
                               c=tvals[valid], cmap="viridis",
                               s=2, alpha=0.7, rasterized=True)
        ax_ut.scatter(umap_xy[~valid,0], umap_xy[~valid,1],
                      c="lightgrey", s=1, alpha=0.2, rasterized=True)
        plt.colorbar(sc_, ax=ax_ut, label="wMel titer", shrink=0.8)
    ax_ut.set_title("B  UMAP — wMel titer")
    ax_ut.set_xticks([]); ax_ut.set_yticks([])

    # C: UMAP CC stage
    if has_umap:
        stage_vals = obs[stage_col].astype(str).str.lower().values
        c_map_vals = [pal.get(s, "grey") for s in stage_vals]
        ax_uc.scatter(umap_xy[:,0], umap_xy[:,1],
                      c=c_map_vals, s=2, alpha=0.7, rasterized=True)
        for s in stages:
            ax_uc.scatter([], [], color=pal[s], label=s, s=20)
        ax_uc.legend(title="CC stage", fontsize=8, markerscale=2,
                     bbox_to_anchor=(1.01,1), loc="upper left")
    ax_uc.set_title("C  UMAP — CC stage (KNN transferred)")
    ax_uc.set_xticks([]); ax_uc.set_yticks([])

    # D: polar titer cyclicity
    if has_pt and df_pol is not None and len(df_pol) >= 20:
        edges   = np.linspace(0, 2*np.pi, n_sectors+1)
        centers = (edges[:-1] + edges[1:]) / 2
        df_pol  = df_pol.copy()
        df_pol["_sector"] = pd.cut(df_pol[pseudotime_col], bins=edges,
                                    labels=range(n_sectors), include_lowest=True)
        sec_med = (df_pol.groupby("_sector", observed=True)[titer_col]
                         .median().reindex(range(n_sectors)).fillna(0).values)
        vmin, vmax = sec_med.min(), sec_med.max()
        norm_vals  = (sec_med - vmin) / (vmax - vmin + 1e-9)
        width      = 2*np.pi / n_sectors
        for theta, r, nv in zip(centers, sec_med, norm_vals):
            ax_pol.bar(theta, r, width=width*0.85, bottom=0,
                       color=plt.cm.coolwarm(nv), alpha=0.85)
        r_annot = sec_med.max() * 1.2
        for name, th0, th1 in [("g0/g1", 0, np.pi/2),
                                 ("s",     np.pi/2, 3*np.pi/2),
                                 ("g2/m",  3*np.pi/2, 2*np.pi)]:
            if name in pal:
                ax_pol.text((th0+th1)/2, r_annot, name,
                            ha="center", va="center",
                            fontsize=7, color=pal[name], fontweight="bold")
        ax_pol.set_theta_zero_location("N")
        ax_pol.set_theta_direction(-1)
        ax_pol.set_xticks(np.linspace(0, 2*np.pi, 5)[:-1])
        ax_pol.set_xticklabels(["0", "π/2", "π", "3π/2"], fontsize=8)
        ax_pol.set_title(f"D  Titer cyclicity\n(n={n_sectors} sectors)", pad=18, fontsize=9)
    else:
        ax_pol.set_title("D  Polar: pseudotime not available", pad=15, fontsize=9)

    plt.suptitle("Q3 — wMel titer across cell-cycle phases", fontweight="bold", y=1.02)
    out = os.path.join(fig_dir, f"q3_titer_by_phase_{sample}.pdf")
    plt.savefig(out, bbox_inches="tight", dpi=150)
    plt.close()
    print(f"   → {out}")


# ─────────────────────────────────────────────────────────────────────────────
# Sanitise + save
# ─────────────────────────────────────────────────────────────────────────────

def _sanitize_obs(adata):
    for col in adata.obs.columns:
        if adata.obs[col].dtype == object or str(adata.obs[col].dtype) == "category":
            adata.obs[col] = adata.obs[col].astype(str).replace("nan", "NA")
    return adata


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def run(ref_path, query_paths, out_path, fig_dir, sample,
        harmony_vars, ref_condition, titer_col, stage_col,
        pseudotime_col, n_pcs, k_neighbors, n_titer_bins, n_sectors,
        infected_only):

    os.makedirs(fig_dir, exist_ok=True)
    sc.settings.figdir = fig_dir

    # Expand any globs in query_paths
    expanded = []
    for p in query_paths:
        matched = sorted(glob.glob(p))
        if matched:
            expanded.extend(matched)
        else:
            expanded.append(p)  # keep as-is; will fail gracefully in load_query_files
    query_paths = expanded
    print(f"Query files ({len(query_paths)}):")
    for p in query_paths:
        print(f"  {p}")

    # ── Load ──────────────────────────────────────────────────────────────────
    ref_raw, ref_full = load_reference(ref_path, stage_col, pseudotime_col)

    print(f"\n── Loading query files ──")
    query_raw = load_query_files(
        query_paths, ref_condition=ref_condition,
        titer_col=titer_col, infected_only=infected_only,
    )

    # ── Joint preprocessing + Harmony ─────────────────────────────────────────
    combined, ref_mask = joint_preprocess_and_harmony(
        ref_raw, query_raw, ref_full,
        harmony_vars=harmony_vars, n_pcs=n_pcs,
    )

    # ── KNN label transfer ────────────────────────────────────────────────────
    query, combined = knn_label_transfer(
        combined, ref_mask, ref_full, k=k_neighbors,
    )

    # ── Assign CC stage from cluster mapping ──────────────────────────────────
    query, cluster_map = assign_cc_from_reference(
        query, ref_full,
        stage_col=stage_col, pseudotime_col=pseudotime_col,
    )
    cluster_map.to_csv(
        os.path.join(fig_dir, f"cluster_cc_map_{sample}.csv"), index=False
    )

    # ── Sanity-check UMAP ─────────────────────────────────────────────────────
    sc.pl.umap(combined, color=["dataset", "leiden"] if "leiden" in combined.obs.columns
               else ["dataset"],
               save=f"_{sample}_combined_dataset.pdf")

    # ── Analyses ──────────────────────────────────────────────────────────────
    obs    = query.obs.copy()
    umap   = query.obsm.get("X_umap", None)

    q1_titer_vs_pseudotime(
        obs, fig_dir, sample,
        titer_col=titer_col, pseudotime_col="cc_pseudotime", stage_col="cc_stage",
    )
    q2_phase_distribution_vs_titer(
        obs, fig_dir, sample,
        titer_col=titer_col, stage_col="cc_stage", n_titer_bins=n_titer_bins,
    )
    q3_titer_by_phase(
        obs, umap, fig_dir, sample,
        titer_col=titer_col, stage_col="cc_stage",
        pseudotime_col="cc_pseudotime", n_sectors=n_sectors,
    )

    # ── Save ──────────────────────────────────────────────────────────────────
    query    = _sanitize_obs(query)
    combined = _sanitize_obs(combined)

    query.write(out_path)
    combined.write(out_path.replace(".h5ad", "_combined.h5ad"))
    print(f"\nSaved query    → {out_path}")
    print(f"Saved combined → {out_path.replace('.h5ad', '_combined.h5ad')}")
    print(f"Figures        → {fig_dir}/")


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
    )


if __name__ == "__main__":
    main()