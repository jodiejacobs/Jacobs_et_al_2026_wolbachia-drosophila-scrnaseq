"""
titer_vs_cellcycle.py
=====================
Tests the relationship between wMel titer and cell cycle in infected Drosophila JW18 cells.

Strategy
--------
1. Load pre-built reference (uninfected, Harmony-corrected, Leiden-clustered, Cyclum-annotated).
2. Subset full integrated object to infected cells (treatment != Ctrl).
3. Ingest query onto reference → transfer Leiden cluster labels.
4. Map Leiden → cyclum_stage using the reference's per-cluster majority vote.
5. Optionally assign pseudotime as per-cluster median from reference.
6. Run three biological questions:
   Q1. wMel titer vs cell-cycle pseudotime (scatter + LOWESS + binned ribbon)
   Q2. Cell-cycle phase distribution across titer bins (chi-squared + stacked bar + heatmap)
   Q3. Titer differs by CC phase (KW + Dunn post-hoc + violin + polar cyclicity)

Run with:
    mamba activate scanpy
    python scripts/method_comparison/titer_vs_cellcycle.py \
        --ref     results/integrated/integrated_uninfected_with_cellcycle.h5ad \
        --query   results/integrated/integrated.h5ad \
        --out_path results/integrated/titer_cellcycle.h5ad \
        --fig_dir  figures/titer_vs_cellcycle
"""

import os
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


# ─────────────────────────────────────────────────────────────────────────────
# Colours
# ─────────────────────────────────────────────────────────────────────────────

CC_ORDER  = ["g0/g1", "s", "g2/m"]
CC_COLORS = {"g0/g1": "#4C72B0", "s": "#DD8452", "g2/m": "#55A868"}

def _cc_palette(stages):
    cmap = matplotlib.colormaps["tab10"]
    return {s: CC_COLORS.get(s, cmap(i % 10)) for i, s in enumerate(stages)}


# ─────────────────────────────────────────────────────────────────────────────
# Step 1 — Ingest query onto reference
# ─────────────────────────────────────────────────────────────────────────────

def ingest_query(adata_query, ref, min_genes, fig_dir, sample):
    """Project infected cells onto uninfected CC reference via sc.tl.ingest.

    Transfers leiden cluster label from reference. CC stage and pseudotime are
    then assigned from the reference cluster→stage/pseudotime mapping in the
    calling function.

    Returns
    -------
    query : AnnData
        Query with X_umap projected from reference and obs["leiden_ref"].
    """
    query = adata_query.copy()
    print(f"\n── Ingest: {query.n_obs} query cells ──")
    print(f"   bio_condition : {query.obs['bio_condition'].value_counts().to_dict()}")
    print(f"   methods       : {query.obs['method'].value_counts().to_dict()}")

    # ── Restrict to reference gene universe ───────────────────────────────────
    ref_genes = ref.raw.var_names if ref.raw is not None else ref.var_names
    query = query[:, query.var_names.isin(ref_genes)].copy()

    # ── QC ────────────────────────────────────────────────────────────────────
    if scipy.sparse.issparse(query.X):
        query.X.eliminate_zeros()
    sc.pp.filter_cells(query, min_genes=min_genes)
    sc.pp.filter_cells(query, min_counts=1)
    sc.pp.filter_genes(query, min_counts=1)

    # ── Densify + sanitise ────────────────────────────────────────────────────
    if scipy.sparse.issparse(query.X):
        query.X = query.X.toarray()
    query.X = np.nan_to_num(query.X.astype(np.float64), nan=0.0, posinf=0.0, neginf=0.0)
    query = query[query.X.sum(axis=1) > 0].copy()
    query = query[:, query.X.sum(axis=0) > 0].copy()

    # ── Store raw counts before scaling ───────────────────────────────────────
    query.raw = query

    # ── Pad missing HVGs with zeros + reorder to match reference ──────────────
    missing = ref.var_names[~ref.var_names.isin(query.var_names)]
    if len(missing):
        print(f"   WARNING: {len(missing)} ref HVGs absent from query — zero-padding")
        zero_cols = np.zeros((query.n_obs, len(missing)), dtype=np.float64)
        query = ad.AnnData(
            X=np.hstack([query.X, zero_cols]),
            obs=query.obs.copy(),
            var=pd.DataFrame(index=list(query.var_names) + list(missing)),
        )
    query = query[:, ref.var_names].copy()
    sc.pp.scale(query, max_value=10)
    query.var["highly_variable"] = True

    # ── sc.tl.ingest ──────────────────────────────────────────────────────────
    obs_backup = query.obs.copy()
    print("   Running sc.tl.ingest …")
    sc.tl.ingest(query, ref, obs="leiden", embedding_method="umap")
    leiden_xfer = query.obs["leiden"].copy()
    query.obs = obs_backup
    query.obs["leiden_ref"] = leiden_xfer.values

    print(f"   Transferred leiden_ref distribution:")
    print(query.obs["leiden_ref"].value_counts().sort_index().to_string())

    return query


# ─────────────────────────────────────────────────────────────────────────────
# Step 2 — Map leiden_ref → CC stage + pseudotime from reference
# ─────────────────────────────────────────────────────────────────────────────

def assign_cc_from_reference(query, ref,
                              stage_col="cyclum_stage",
                              pseudotime_col="cyclum_pseudotime",
                              leiden_col="leiden"):
    """Assign CC stage and pseudotime to query cells using reference cluster mappings.

    For each leiden cluster in the reference:
      - CC stage   : majority vote of cyclum_stage within that cluster
      - Pseudotime : median cyclum_pseudotime within that cluster

    These are cluster-level summaries, not single-cell projections — appropriate
    for downstream group comparisons but not for fine-grained pseudotime analysis.
    Both assignments are flagged in the query obs with '_source' columns.

    Parameters
    ----------
    query : AnnData
        Query with obs["leiden_ref"] from ingest.
    ref : AnnData
        Reference with obs[stage_col] and obs[pseudotime_col].

    Returns
    -------
    query : AnnData
        Query with obs["cc_stage"] and obs["cc_pseudotime"] added.
    cluster_map : pd.DataFrame
        Per-cluster mapping table (leiden → cc_stage, median_pseudotime, n_cells, purity).
    """
    if stage_col not in ref.obs.columns:
        raise ValueError(f"'{stage_col}' not in reference obs. "
                         f"Available: {list(ref.obs.columns)}")

    ref_obs = ref.obs[[leiden_col, stage_col]].copy()
    ref_obs[leiden_col] = ref_obs[leiden_col].astype(str)

    # Majority-vote CC stage per cluster
    def majority(x):
        return x.value_counts().idxmax()

    def purity(x):
        vc = x.value_counts()
        return vc.iloc[0] / vc.sum()

    stage_map = (ref_obs.groupby(leiden_col)[stage_col]
                         .agg(cc_stage=majority, purity=purity)
                         .reset_index()
                         .rename(columns={leiden_col: "cluster"}))

    # Median pseudotime per cluster (if available)
    if pseudotime_col in ref.obs.columns:
        pt_map = (ref.obs[[leiden_col, pseudotime_col]]
                     .copy()
                     .assign(**{leiden_col: ref.obs[leiden_col].astype(str)})
                     .groupby(leiden_col)[pseudotime_col]
                     .median()
                     .reset_index()
                     .rename(columns={leiden_col: "cluster",
                                      pseudotime_col: "median_pseudotime"}))
        stage_map = stage_map.merge(pt_map, on="cluster", how="left")
    else:
        stage_map["median_pseudotime"] = np.nan

    # Count reference cells per cluster
    n_ref = ref_obs[leiden_col].value_counts().rename("n_ref_cells").reset_index()
    n_ref.columns = ["cluster", "n_ref_cells"]
    stage_map = stage_map.merge(n_ref, on="cluster", how="left")

    print("\n── Cluster → CC stage mapping (from reference) ──")
    print(stage_map.sort_values("cluster").to_string(index=False))

    # Warn about low-purity clusters
    low_purity = stage_map[stage_map["purity"] < 0.6]
    if len(low_purity):
        print(f"\n   WARNING: {len(low_purity)} clusters have CC stage purity < 60%:")
        print(low_purity[["cluster", "cc_stage", "purity"]].to_string(index=False))

    # Apply to query
    cluster_to_stage = dict(zip(stage_map["cluster"], stage_map["cc_stage"]))
    cluster_to_pt    = dict(zip(stage_map["cluster"], stage_map["median_pseudotime"]))

    query.obs["leiden_ref"]   = query.obs["leiden_ref"].astype(str)
    query.obs["cc_stage"]     = query.obs["leiden_ref"].map(cluster_to_stage)
    query.obs["cc_pseudotime"] = query.obs["leiden_ref"].map(cluster_to_pt).astype(float)

    n_unmapped = query.obs["cc_stage"].isna().sum()
    if n_unmapped:
        print(f"\n   WARNING: {n_unmapped} query cells did not map to a known cluster")

    print(f"\n   Query CC stage distribution (transferred):")
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
        print(f"   SKIP: {missing} not in obs")
        return

    df = obs[[titer_col, pseudotime_col, stage_col]].copy()
    df[titer_col]      = pd.to_numeric(df[titer_col],      errors="coerce")
    df[pseudotime_col] = pd.to_numeric(df[pseudotime_col], errors="coerce")
    df = df.dropna()
    print(f"   {len(df)} cells with titer + pseudotime")

    # Check how many distinct pseudotime values exist (cluster-median = few unique values)
    n_unique_pt = df[pseudotime_col].nunique()
    print(f"   Unique pseudotime values: {n_unique_pt} (cluster-median resolution)")

    if len(df) < 20:
        print("   SKIP: fewer than 20 cells")
        return

    pt  = df[pseudotime_col].values
    tit = df[titer_col].values

    stages = [s for s in CC_ORDER if s in df[stage_col].unique()]
    stages += [s for s in df[stage_col].unique() if s not in CC_ORDER]
    pal    = _cc_palette(stages)

    # Spearman
    rho, p_val = spearmanr(pt, tit)
    print(f"   Spearman rho={rho:.3f}, p={p_val:.3e}")

    pd.DataFrame({
        "spearman_rho": [rho], "p_value": [p_val], "n_cells": [len(df)],
        "note": ["cc_pseudotime is cluster-median from reference; "
                 "Spearman is approximate and at cluster resolution"],
    }).to_csv(os.path.join(fig_dir, f"q1_spearman_{sample}.csv"), index=False)

    # LOWESS — only if enough unique x values
    do_lowess = n_unique_pt >= 5
    if do_lowess:
        order    = np.argsort(pt)
        smoothed = lowess(tit[order], pt[order], frac=0.3, return_sorted=True)

    # Binned summary — cap n_bins to number of unique pseudotime values
    # to avoid empty bins that break the centers alignment
    n_bins_actual = min(n_bins, n_unique_pt)
    if n_bins_actual < 2:
        print("   WARNING: fewer than 2 unique pseudotime values — skipping binned plot")
        do_bins = False
    else:
        do_bins = True
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
        # Align centers to actual non-empty bins
        occupied_labels = bsum["_pt_bin"].astype(str).values
        centers_aligned = np.array([
            centers[bin_lbl.index(b)] for b in occupied_labels if b in bin_lbl
        ])
        bsum["bin_center"] = centers_aligned
        bsum.to_csv(os.path.join(fig_dir, f"q1_bin_summary_{sample}.csv"), index=False)

    pi_ticks  = [0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi]
    pi_labels = ["0", "π/2", "π", "3π/2", "2π"]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # A: scatter coloured by CC stage
    for stage in stages:
        m = df[stage_col] == stage
        axes[0].scatter(df.loc[m, pseudotime_col], df.loc[m, titer_col],
                        c=[pal[stage]], s=3, alpha=0.4, label=stage, rasterized=True)
    axes[0].set_xlabel("CC pseudotime (cluster-median, radians)")
    axes[0].set_ylabel("wMel titer")
    axes[0].set_title(f"A  Titer vs pseudotime\nSpearman ρ={rho:.3f}, p={p_val:.2e}\n"
                      f"({n_unique_pt} unique pseudotime values — cluster resolution)")
    axes[0].set_xticks(pi_ticks); axes[0].set_xticklabels(pi_labels)
    axes[0].legend(title="CC stage", fontsize=8, markerscale=3)

    # B: LOWESS or message
    if do_lowess:
        axes[1].scatter(pt, tit, c="lightgrey", s=2, alpha=0.25, rasterized=True)
        axes[1].plot(smoothed[:, 0], smoothed[:, 1],
                     color="#d62728", lw=2, label="LOWESS (frac=0.3)")
        axes[1].set_xticks(pi_ticks); axes[1].set_xticklabels(pi_labels)
        axes[1].legend()
    else:
        axes[1].text(0.5, 0.5, f"Too few unique pseudotime\nvalues for LOWESS\n(n={n_unique_pt})",
                     ha="center", va="center", transform=axes[1].transAxes, fontsize=10)
    axes[1].set_xlabel("CC pseudotime (cluster-median, radians)")
    axes[1].set_ylabel("wMel titer")
    axes[1].set_title("B  LOWESS smoothed trend")

    # C: binned median ± IQR
    if do_bins and len(bsum) >= 2:
        x_ = bsum["bin_center"].values
        axes[2].plot(x_, bsum["median"].values, "o-", color="#d62728",
                     lw=2, ms=5, label="Median")
        axes[2].fill_between(x_, bsum["q25"].values, bsum["q75"].values,
                             alpha=0.25, color="#d62728", label="IQR")
        axes[2].set_xticks(x_)
        axes[2].set_xticklabels([f"{v:.2f}" for v in x_], rotation=45)
        axes[2].legend()
    else:
        axes[2].text(0.5, 0.5, "Insufficient pseudotime\nresolution for binning",
                     ha="center", va="center", transform=axes[2].transAxes, fontsize=10)
    axes[2].set_xlabel("Pseudotime bin center (radians)")
    axes[2].set_ylabel("wMel titer")
    axes[2].set_title(f"C  Median titer ± IQR  ({n_bins_actual} bins)")

    plt.suptitle("Q1 — wMel titer vs cell-cycle pseudotime (cluster-median)",
                 fontweight="bold")
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
    """Q2: Does cell-cycle phase composition shift with wMel titer?

    Outputs
    -------
    PDF : q2_phase_dist_vs_titer_{sample}.pdf
    CSV : q2_composition_{sample}.csv
          q2_chisq_{sample}.csv
    """
    print(f"\n── Q2: phase distribution vs titer ──")

    df = obs[[titer_col, stage_col]].copy()
    df[titer_col] = pd.to_numeric(df[titer_col], errors="coerce")
    df[stage_col] = df[stage_col].astype(str).str.strip().str.lower()
    df = df.dropna()
    print(f"   {len(df)} cells")

    stages = [s for s in CC_ORDER if s in df[stage_col].unique()]
    stages += [s for s in df[stage_col].unique() if s not in CC_ORDER]
    pal    = _cc_palette(stages)

    # Quantile-bin titer
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
        "chi2": [chi2], "dof": [dof], "p_value": [p_chi],
        "n_cells": [len(df)], "n_titer_bins": [len(actual_bins)],
    }).to_csv(os.path.join(fig_dir, f"q2_chisq_{sample}.csv"), index=False)

    # Titer ranges per bin for x-axis annotation
    bin_ranges = (df.groupby("titer_bin", observed=True)[titer_col]
                    .agg(lo="min", hi="max").reset_index())
    bin_annot  = {row["titer_bin"]: f"{row['lo']:.2f}–{row['hi']:.2f}"
                  for _, row in bin_ranges.iterrows()}

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # A: stacked bar
    bottom = np.zeros(len(prop))
    x_pos  = np.arange(len(prop))
    for stage in stages:
        vals = prop[stage].values
        axes[0].bar(x_pos, vals, bottom=bottom,
                    color=pal[stage], label=stage, width=0.7)
        bottom += vals
    axes[0].set_xticks(x_pos)
    axes[0].set_xticklabels(
        [f"{b}\n({bin_annot.get(b,'')})" for b in prop.index], fontsize=8
    )
    axes[0].set_xlabel(f"wMel titer quantile (Q1=lowest, Q{len(actual_bins)}=highest)")
    axes[0].set_ylabel("% of cells")
    axes[0].set_title(
        f"A  CC phase composition per titer bin\n"
        f"χ²={chi2:.1f}, df={dof}, p={p_chi:.2e}"
    )
    axes[0].legend(title="CC stage", bbox_to_anchor=(1.01, 1),
                   loc="upper left", fontsize=9)

    # B: heatmap
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
    """Q3: Do cells in different CC phases carry different wMel titers,
    and is titer cyclical around the cell cycle?

    Tests
    -----
    Kruskal-Wallis omnibus + Dunn post-hoc pairwise (BH FDR).

    Outputs
    -------
    PDF : q3_titer_by_phase_{sample}.pdf
    CSV : q3_kruskal_{sample}.csv
          q3_dunn_{sample}.csv
          q3_stage_summary_{sample}.csv
    """
    print(f"\n── Q3: titer by CC phase ──")

    df = obs[[titer_col, stage_col]].copy()
    df[titer_col] = pd.to_numeric(df[titer_col], errors="coerce")
    df[stage_col] = df[stage_col].astype(str).str.strip().str.lower()
    df = df.dropna()
    print(f"   {len(df)} cells")

    stages = [s for s in CC_ORDER if s in df[stage_col].unique()]
    stages += [s for s in df[stage_col].unique() if s not in CC_ORDER]
    pal    = _cc_palette(stages)

    # ── Kruskal-Wallis ────────────────────────────────────────────────────────
    groups  = [df.loc[df[stage_col] == s, titer_col].values for s in stages]
    kw_stat, kw_p = kruskal(*groups)
    print(f"   Kruskal-Wallis: H={kw_stat:.3f}, p={kw_p:.3e}")

    pd.DataFrame({
        "statistic": [kw_stat], "p_value": [kw_p],
        "n_stages": [len(stages)], "n_cells": [len(df)],
    }).to_csv(os.path.join(fig_dir, f"q3_kruskal_{sample}.csv"), index=False)

    # ── Dunn post-hoc ─────────────────────────────────────────────────────────
    all_ranks  = scipy.stats.rankdata(df[titer_col].values)
    n_total    = len(all_ranks)
    _, counts  = np.unique(df[titer_col].values, return_counts=True)
    tie_factor = np.sum(counts ** 3 - counts) / (12 * (n_total - 1)) if n_total > 1 else 0
    df = df.copy()
    df["_rank"] = all_ranks

    rows = []
    for s_a, s_b in combinations(stages, 2):
        ga = df.loc[df[stage_col] == s_a, "_rank"].values
        gb = df.loc[df[stage_col] == s_b, "_rank"].values
        na, nb = len(ga), len(gb)
        if na < 2 or nb < 2:
            continue
        se = np.sqrt((n_total * (n_total + 1) / 12 - tie_factor) * (1 / na + 1 / nb))
        if se == 0:
            continue
        z = (ga.mean() - gb.mean()) / se
        rows.append({
            "stage_A": s_a, "stage_B": s_b,
            "median_titer_A": np.median(df.loc[df[stage_col] == s_a, titer_col]),
            "median_titer_B": np.median(df.loc[df[stage_col] == s_b, titer_col]),
            "z_stat": z,
            "p_raw": 2 * scipy.stats.norm.sf(abs(z)),
            "n_A": na, "n_B": nb,
        })

    dunn_df = pd.DataFrame(rows)
    if len(dunn_df):
        _, p_adj, _, _ = multipletests(dunn_df["p_raw"], method="fdr_bh")
        dunn_df["p_adj_BH"]  = p_adj
        dunn_df["significant"] = dunn_df["p_adj_BH"] < 0.05
        dunn_df = dunn_df.sort_values("p_adj_BH")
        dunn_df.to_csv(os.path.join(fig_dir, f"q3_dunn_{sample}.csv"), index=False)
        print(f"   Dunn: {dunn_df['significant'].sum()}/{len(dunn_df)} pairs significant")
        print(dunn_df[["stage_A", "stage_B", "median_titer_A",
                        "median_titer_B", "p_adj_BH", "significant"]].to_string(index=False))

    # Stage summary
    stage_med = df.groupby(stage_col)[titer_col].median().to_dict()
    (df.groupby(stage_col)[titer_col]
       .agg(n_cells="count", median="median", mean="mean", std="std",
            q25=lambda x: x.quantile(0.25), q75=lambda x: x.quantile(0.75))
       .reindex(stages)
       .reset_index()
       .rename(columns={stage_col: "stage"})
       .to_csv(os.path.join(fig_dir, f"q3_stage_summary_{sample}.csv"), index=False))

    # ── Polar data ────────────────────────────────────────────────────────────
    has_pt = pseudotime_col in obs.columns
    df_pol = None
    if has_pt:
        df_pol = obs[[titer_col, pseudotime_col, stage_col]].copy()
        df_pol[titer_col]      = pd.to_numeric(df_pol[titer_col],      errors="coerce")
        df_pol[pseudotime_col] = pd.to_numeric(df_pol[pseudotime_col], errors="coerce")
        df_pol = df_pol.dropna()

    # ── Figure ────────────────────────────────────────────────────────────────
    has_umap = umap_xy is not None
    fig = plt.figure(figsize=(22, 5))
    gs  = fig.add_gridspec(1, 4, wspace=0.45)
    ax_vln = fig.add_subplot(gs[0])
    ax_ut  = fig.add_subplot(gs[1])
    ax_uc  = fig.add_subplot(gs[2])
    ax_pol = fig.add_subplot(gs[3], projection="polar")

    # Panel A — violin + strip + significance bars
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
        for k, (_, row) in enumerate(sig.iterrows()):
            if row["stage_A"] not in stages or row["stage_B"] not in stages:
                continue
            xi = stages.index(row["stage_A"])
            xj = stages.index(row["stage_B"])
            y  = y_max + y_step * (k + 1)
            ax_vln.plot([xi, xj], [y, y], color="black", lw=1)
            ax_vln.text((xi + xj) / 2, y + y_step * 0.15, "*",
                        ha="center", va="bottom", fontsize=10)
    ax_vln.set_xlabel("Cell-cycle stage")
    ax_vln.set_ylabel("wMel titer")
    ax_vln.set_title(f"A  Titer by CC stage\nKW p={kw_p:.2e}")

    # Panel B — UMAP titer
    if has_umap:
        tvals = obs[titer_col].values.astype(float)
        valid = ~np.isnan(tvals)
        sc_   = ax_ut.scatter(umap_xy[valid, 0], umap_xy[valid, 1],
                               c=tvals[valid], cmap="viridis",
                               s=2, alpha=0.7, rasterized=True)
        ax_ut.scatter(umap_xy[~valid, 0], umap_xy[~valid, 1],
                      c="lightgrey", s=1, alpha=0.2, rasterized=True)
        plt.colorbar(sc_, ax=ax_ut, label="wMel titer", shrink=0.8)
    ax_ut.set_title("B  UMAP — wMel titer")
    ax_ut.set_xticks([]); ax_ut.set_yticks([])

    # Panel C — UMAP CC stage
    if has_umap:
        stage_vals  = obs[stage_col].astype(str).str.lower().values
        c_map_vals  = [pal.get(s, "grey") for s in stage_vals]
        ax_uc.scatter(umap_xy[:, 0], umap_xy[:, 1],
                      c=c_map_vals, s=2, alpha=0.7, rasterized=True)
        for s in stages:
            ax_uc.scatter([], [], color=pal[s], label=s, s=20)
        ax_uc.legend(title="CC stage", fontsize=8, markerscale=2,
                     bbox_to_anchor=(1.01, 1), loc="upper left")
    ax_uc.set_title("C  UMAP — CC stage (transferred)")
    ax_uc.set_xticks([]); ax_uc.set_yticks([])

    # Panel D — polar: median titer per pseudotime sector
    if has_pt and df_pol is not None and len(df_pol) >= 20:
        edges   = np.linspace(0, 2 * np.pi, n_sectors + 1)
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
        width      = 2 * np.pi / n_sectors
        for theta, r, nv in zip(centers, sec_med, norm_vals):
            ax_pol.bar(theta, r, width=width * 0.85, bottom=0,
                       color=plt.cm.coolwarm(nv), alpha=0.85)
        # CC stage arc annotations at outer edge
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
        ax_pol.set_xticks(np.linspace(0, 2 * np.pi, 5)[:-1])
        ax_pol.set_xticklabels(["0", "π/2", "π", "3π/2"], fontsize=8)
        ax_pol.set_title(
            f"D  Titer cyclicity\n(median per pseudotime sector, n={n_sectors})",
            pad=18, fontsize=9
        )
    else:
        ax_pol.set_title("D  Polar: pseudotime not available", pad=15, fontsize=9)

    plt.suptitle("Q3 — wMel titer across cell-cycle phases", fontweight="bold", y=1.02)
    out = os.path.join(fig_dir, f"q3_titer_by_phase_{sample}.pdf")
    plt.savefig(out, bbox_inches="tight", dpi=150)
    plt.close()
    print(f"   → {out}")


# ─────────────────────────────────────────────────────────────────────────────
# Sanitise obs
# ─────────────────────────────────────────────────────────────────────────────

def _sanitize_obs(adata):
    for col in adata.obs.columns:
        if adata.obs[col].dtype == object or str(adata.obs[col].dtype) == "category":
            adata.obs[col] = adata.obs[col].astype(str).replace("nan", "NA")
    return adata


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def run(ref_path, query_path, out_path, fig_dir, sample,
        min_genes, titer_col, stage_col, pseudotime_col,
        ref_condition, n_titer_bins, n_sectors):

    os.makedirs(fig_dir, exist_ok=True)
    sc.settings.figdir = fig_dir

    # ── Load ──────────────────────────────────────────────────────────────────
    adata_full = sc.read_h5ad(query_path)
    ref        = sc.read_h5ad(ref_path)

    print(f"Loaded query : {adata_full.n_obs} cells")
    print(f"Loaded ref   : {ref.n_obs} cells, {ref.obs['leiden'].nunique()} clusters")

    # ── Subset to infected (query) cells ─────────────────────────────────────
    ref_conditions = [ref_condition] if isinstance(ref_condition, str) else ref_condition
    query_mask = (
        ~adata_full.obs["cell_line"].isin(ref_conditions) |
        (adata_full.obs["treatment"] != "Ctrl")
    )
    print(f"\nReference conditions : {ref_conditions} + treatment=Ctrl")
    print(f"Query cells          : {query_mask.sum()}/{adata_full.n_obs}")
    print(adata_full.obs.loc[query_mask, "bio_condition"].value_counts().to_string())

    if query_mask.sum() == 0:
        raise ValueError(
            "No query cells found.\n"
            f"  cell_line values : {adata_full.obs['cell_line'].unique().tolist()}\n"
            f"  treatment values : {adata_full.obs['treatment'].unique().tolist()}\n"
            f"  ref_condition    : {ref_conditions}"
        )

    adata_query = adata_full[query_mask].copy()

    # Confirm titer is present
    if titer_col not in adata_query.obs.columns:
        raise ValueError(f"'{titer_col}' not in query obs. "
                         f"Available: {list(adata_query.obs.columns)}")
    n_titer = pd.to_numeric(adata_query.obs[titer_col], errors="coerce").notna().sum()
    print(f"\nQuery cells with non-null {titer_col}: {n_titer}/{adata_query.n_obs}")

    # ── Ingest ────────────────────────────────────────────────────────────────
    query = ingest_query(adata_query, ref=ref,
                         min_genes=min_genes, fig_dir=fig_dir, sample=sample)

    # ── Assign CC stage + pseudotime from reference cluster mapping ───────────
    query, cluster_map = assign_cc_from_reference(
        query, ref,
        stage_col=stage_col,
        pseudotime_col=pseudotime_col,
        leiden_col="leiden",
    )
    # ── Diagnose CC transfer quality ──────────────────────────────────────────
    print("\n── CC transfer diagnostics ──")
    print("Query leiden_ref distribution:")
    print(query.obs["leiden_ref"].value_counts().sort_index().to_string())
    print("\nReference leiden distribution:")
    print(ref.obs["leiden"].value_counts().sort_index().to_string())
    print("\nCluster → stage mapping used:")
    print(cluster_map[["cluster", "cc_stage", "purity", "n_ref_cells"]].to_string(index=False))

    # Check if query cells are piling into a subset of clusters
    query_cluster_counts = query.obs["leiden_ref"].value_counts()
    print(f"\nQuery cells: top 5 clusters by count:")
    print(query_cluster_counts.head(5).to_string())
    print(f"Query cells in clusters mapped to 's': "
        f"{query.obs.loc[query.obs['cc_stage']=='s', 'leiden_ref'].nunique()} distinct clusters")


    cluster_map.to_csv(os.path.join(fig_dir, f"cluster_cc_map_{sample}.csv"), index=False)

    # ── Run analyses ──────────────────────────────────────────────────────────
    obs    = query.obs.copy()
    umap   = query.obsm.get("X_umap", None)

    q1_titer_vs_pseudotime(
        obs, fig_dir, sample,
        titer_col=titer_col,
        pseudotime_col="cc_pseudotime",
        stage_col="cc_stage",
    )
    q2_phase_distribution_vs_titer(
        obs, fig_dir, sample,
        titer_col=titer_col,
        stage_col="cc_stage",
        n_titer_bins=n_titer_bins,
    )
    q3_titer_by_phase(
        obs, umap, fig_dir, sample,
        titer_col=titer_col,
        stage_col="cc_stage",
        pseudotime_col="cc_pseudotime",
        n_sectors=n_sectors,
    )

    # ── Save ──────────────────────────────────────────────────────────────────
    query = _sanitize_obs(query)
    query.write(out_path)
    print(f"\nSaved → {out_path}")
    print(f"Figures → {fig_dir}/")


def main():
    parser = argparse.ArgumentParser(
        description="Ingest infected cells onto CC reference, then test titer vs cell cycle"
    )
    parser.add_argument("--ref",            required=True,
                        help="Reference h5ad: uninfected, Leiden-clustered, cyclum-annotated")
    parser.add_argument("--query",          required=True,
                        help="Full integrated h5ad containing infected cells")
    parser.add_argument("--out_path",       default="results/integrated/titer_cellcycle.h5ad")
    parser.add_argument("--fig_dir",        default="figures/titer_vs_cellcycle")
    parser.add_argument("--sample",         default="wolbachia_infection")
    parser.add_argument("--min_genes",      type=int,   default=200)
    parser.add_argument("--titer_col",      default="wolbachia_titer")
    parser.add_argument("--stage_col",      default="cyclum_stage")
    parser.add_argument("--pseudotime_col", default="cyclum_pseudotime")
    parser.add_argument("--ref_condition",  nargs="+",  default=["JW18DOX"],
                        help="cell_line value(s) that are uninfected reference "
                             "(combined with treatment=Ctrl to define what is NOT query)")
    parser.add_argument("--n_titer_bins",   type=int,   default=5)
    parser.add_argument("--n_sectors",      type=int,   default=12,
                        help="Polar plot sectors for Q3 cyclicity panel")
    args = parser.parse_args()

    run(
        ref_path=args.ref,
        query_path=args.query,
        out_path=args.out_path,
        fig_dir=args.fig_dir,
        sample=args.sample,
        min_genes=args.min_genes,
        titer_col=args.titer_col,
        stage_col=args.stage_col,
        pseudotime_col=args.pseudotime_col,
        ref_condition=args.ref_condition,
        n_titer_bins=args.n_titer_bins,
        n_sectors=args.n_sectors,
    )


if __name__ == "__main__":
    main()