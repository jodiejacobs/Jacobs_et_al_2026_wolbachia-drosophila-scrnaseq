"""
titer_vs_cellcycle.py
=====================
Tests the relationship between wMel titer and cell cycle phase in Drosophila JW18 cells.

Strategy
--------
1. Load pre-built reference (uninfected, Harmony-corrected, Leiden-clustered, UMAP embedded).
   Leiden clusters serve as cell-cycle phase proxies (validated by Cyclum pseudotime).
2. Load query (infected/SV timepoint cells with wolbachia_titer in obs).
3. Ingest query onto reference → transfer Leiden cluster labels as cell-cycle phase proxy.
4. Test three biological questions:
   Q1. Is wMel titer correlated with cell-cycle pseudotime? (linear + LOWESS scatter)
   Q2. Does cell-cycle phase distribution shift with titer? (titer-binned cluster composition)
   Q3. Do cells in different cell-cycle phases carry different titers? (KW + Dunn post-hoc)

Run with:
    mamba activate scanpy
    python scripts/method_comparison/titer_vs_cellcycle.py \
        --ref   results/integrated_by_refintegrated_uninfected_with_cellcycle.h5ad \
        --query results/integrated_by_refintegrated_all_timepoints.h5ad \
        --out_path results/integrated_by_reftiter_cellcycle.h5ad \
        --fig_dir  results/integrated_by_ref/figures/titer_vs_cellcycle \
        --titer_col wolbachia_titer \
        --pseudotime_col cyclum_pseudotime \
        --ref_condition JW18DOX
"""

import os
import argparse

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import scipy.sparse
import scipy.stats
from scipy.stats import kruskal, spearmanr
from itertools import combinations
from statsmodels.stats.multitest import multipletests
from statsmodels.nonparametric.smoothers_lowess import lowess
import anndata as ad
import scanpy as sc


# ─────────────────────────────────────────────────────────────────────────────
# Colour helpers
# ─────────────────────────────────────────────────────────────────────────────

def _cluster_palette(clusters):
    """Return a dict {cluster_label: RGBA} using tab20."""
    cmap = matplotlib.colormaps["tab20"]
    return {cl: cmap(i % 20) for i, cl in enumerate(clusters)}


# ─────────────────────────────────────────────────────────────────────────────
# Step 1 — Ingest query onto reference
# ─────────────────────────────────────────────────────────────────────────────

def ingest_query(adata_query, ref, min_genes, fig_dir, sample):
    """Project infected-timepoint cells onto the uninfected cell-cycle reference.

    Parameters
    ----------
    adata_query : AnnData
        Raw-count query cells (infected / SV timepoints) with wolbachia_titer in obs.
    ref : AnnData
        Pre-built reference: HVGs in var, scaled X, PCA + UMAP in obsm,
        Leiden cluster labels in obs["leiden"].
    min_genes : int
        Minimum genes per cell (QC filter applied to query only).
    fig_dir : str
        Figure output directory.
    sample : str
        Label used in output filenames.

    Returns
    -------
    query : AnnData
        Query with X_umap projected from reference and obs["cc_cluster"] = transferred
        Leiden label (cell-cycle phase proxy).
    """
    query = adata_query.copy()
    print(f"\n── Ingest: {query.n_obs} query cells ──")
    print(f"   Timepoints : {query.obs['timepoint'].value_counts().to_dict()}")
    print(f"   Methods    : {query.obs['method'].value_counts().to_dict()}")

    # ── Gene universe: subset query to HVGs present in reference ──────────────
    ref_genes = ref.raw.var_names if ref.raw is not None else ref.var_names
    query = query[:, query.var_names.isin(ref_genes)].copy()

    # ── QC filtering ──────────────────────────────────────────────────────────
    if scipy.sparse.issparse(query.X):
        query.X.eliminate_zeros()
    sc.pp.filter_cells(query, min_genes=min_genes)
    sc.pp.filter_cells(query, min_counts=1)
    sc.pp.filter_genes(query, min_counts=1)

    # ── Densify and sanitise ──────────────────────────────────────────────────
    if scipy.sparse.issparse(query.X):
        query.X = query.X.toarray()
    query.X = np.nan_to_num(query.X.astype(np.float64), nan=0.0, posinf=0.0, neginf=0.0)
    query = query[query.X.sum(axis=1) > 0].copy()
    query = query[:, query.X.sum(axis=0) > 0].copy()

    # ── Store raw counts before scaling ───────────────────────────────────────
    query.raw = query

    # ── Pad missing reference HVGs with zeros ─────────────────────────────────
    missing = ref.var_names[~ref.var_names.isin(query.var_names)]
    if len(missing):
        print(f"   WARNING: {len(missing)} ref HVGs absent from query — zero-padding")
        zero_cols = np.zeros((query.n_obs, len(missing)), dtype=np.float64)
        query = ad.AnnData(
            X=np.hstack([query.X, zero_cols]),
            obs=query.obs.copy(),
            var=pd.DataFrame(index=list(query.var_names) + list(missing)),
        )

    # Reorder to match reference exactly (required by sc.tl.ingest)
    query = query[:, ref.var_names].copy()
    sc.pp.scale(query, max_value=10)
    query.var["highly_variable"] = True

    # ── sc.tl.ingest: project UMAP + transfer Leiden ──────────────────────────
    obs_backup = query.obs.copy()
    print("   Running sc.tl.ingest …")
    sc.tl.ingest(query, ref, obs="leiden", embedding_method="umap")
    leiden_xfer = query.obs["leiden"].copy()
    query.obs = obs_backup
    query.obs["cc_cluster"] = leiden_xfer.values  # cell-cycle cluster proxy

    print(f"   Cell-cycle cluster distribution (query):")
    print(query.obs["cc_cluster"].value_counts().sort_index().to_string())

    # ── Quick sanity-check UMAP ───────────────────────────────────────────────
    sc.pl.umap(query, color=["cc_cluster", "timepoint", "method"],
               save=f"_{sample}_query_ingested.pdf", ncols=3,
               title=["CC cluster (from ref)", "Timepoint", "Method"])

    # Joint ref + query UMAP
    ref_plot = ref.copy()
    ref_plot.obs["dataset"]    = "uninfected (ref)"
    ref_plot.obs["cc_cluster"] = ref_plot.obs["leiden"]
    ref_plot.obs["timepoint"]  = ref_plot.obs["timepoint"].astype(str).fillna("uninfected")
    query.obs["dataset"]       = "infected (" + query.obs["timepoint"].astype(str).fillna("?") + ")"

    common = ref_plot.var_names.intersection(query.var_names)
    combined_umap = ad.concat(
        [ref_plot[:, common], query[:, common]], join="inner", index_unique="-"
    )
    sc.pl.umap(combined_umap, color=["dataset", "cc_cluster", "timepoint"],
               save=f"_{sample}_combined_ref_query.pdf", ncols=3,
               title=["Dataset", "CC cluster", "Timepoint"])

    return query


# ─────────────────────────────────────────────────────────────────────────────
# Q1 — titer vs cell-cycle pseudotime (linear scatter + LOWESS)
# ─────────────────────────────────────────────────────────────────────────────

def q1_titer_vs_pseudotime(query, fig_dir, sample,
                            titer_col="wolbachia_titer",
                            pseudotime_col="cyclum_pseudotime",
                            cluster_col="cc_cluster",
                            n_bins=8):
    """Q1: Is wMel titer correlated with position in the cell cycle?

    Cyclum pseudotime is circular (0–2π). Linear Spearman is an approximation;
    this is flagged explicitly in the output CSV.

    Outputs
    -------
    PDF : q1_titer_vs_pseudotime_{sample}.pdf   (scatter + LOWESS + binned ribbon)
    CSV : q1_titer_vs_pseudotime_{sample}_spearman.csv
          q1_titer_vs_pseudotime_{sample}_bin_summary.csv
    """
    tag = f"q1_titer_vs_pseudotime_{sample}"
    obs = query.obs.copy()

    missing = [c for c in [titer_col, pseudotime_col] if c not in obs.columns]
    if missing:
        print(f"  Q1 SKIP: {missing} not in query.obs")
        return

    obs[titer_col]      = pd.to_numeric(obs[titer_col],      errors="coerce")
    obs[pseudotime_col] = pd.to_numeric(obs[pseudotime_col], errors="coerce")
    obs = obs.dropna(subset=[titer_col, pseudotime_col])

    if len(obs) < 20:
        print("  Q1 SKIP: fewer than 20 cells with both titer and pseudotime")
        return

    print(f"\n── Q1: titer vs pseudotime  ({len(obs)} cells) ──")

    pt  = obs[pseudotime_col].values
    tit = obs[titer_col].values

    # Spearman
    rho, p_val = spearmanr(pt, tit)
    print(f"   Spearman rho={rho:.3f}, p={p_val:.3e}  (circular pseudotime; interpret with caution)")

    pd.DataFrame({
        "spearman_rho": [rho], "p_value": [p_val], "n_cells": [len(obs)],
        "note": ["cyclum pseudotime is circular (0-2pi); linear Spearman is approximate"],
    }).to_csv(os.path.join(fig_dir, f"{tag}_spearman.csv"), index=False)

    # LOWESS
    order    = np.argsort(pt)
    smoothed = lowess(tit[order], pt[order], frac=0.2, return_sorted=True)

    # Binned summary
    edges   = np.linspace(pt.min(), pt.max(), n_bins + 1)
    centers = (edges[:-1] + edges[1:]) / 2
    labels  = [f"{e:.2f}" for e in edges[:-1]]
    obs["_pt_bin"] = pd.cut(obs[pseudotime_col], bins=edges, labels=labels, include_lowest=True)
    obs_b = obs.dropna(subset=["_pt_bin"])
    bsum = (obs_b.groupby("_pt_bin", observed=True)[titer_col]
            .agg(n_cells="count", median="median", mean="mean",
                 q25=lambda x: x.quantile(0.25), q75=lambda x: x.quantile(0.75))
            .reset_index())
    bsum["bin_center"] = centers
    bsum.to_csv(os.path.join(fig_dir, f"{tag}_bin_summary.csv"), index=False)

    # ── Figure: 3-panel ───────────────────────────────────────────────────────
    # Panel A: scatter coloured by cc_cluster
    # Panel B: LOWESS trend line
    # Panel C: binned median ± IQR ribbon
    has_clusters = cluster_col in obs.columns
    clusters, pal = [], {}
    if has_clusters:
        obs[cluster_col] = obs[cluster_col].astype(str)
        clusters = sorted(obs[cluster_col].unique(), key=lambda x: int(x) if x.isdigit() else x)
        pal = _cluster_palette(clusters)

    pi_ticks = [0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi]
    pi_labels = ["0", "π/2", "π", "3π/2", "2π"]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Panel A — scatter coloured by cell-cycle cluster
    if has_clusters:
        for cl in clusters:
            m = obs[cluster_col] == cl
            axes[0].scatter(pt[m.values], tit[m.values],
                            c=[pal[cl]], s=3, alpha=0.4, label=f"C{cl}", rasterized=True)
        axes[0].legend(fontsize=6, markerscale=2,
                       bbox_to_anchor=(1.01, 1), loc="upper left", title="CC cluster")
    else:
        axes[0].scatter(pt, tit, s=3, alpha=0.4, rasterized=True)
    axes[0].set_xlabel("Cyclum pseudotime (radians)")
    axes[0].set_ylabel(f"wMel titer")
    axes[0].set_title(f"A  Titer vs pseudotime\nSpearman ρ={rho:.3f}, p={p_val:.2e}")
    axes[0].set_xticks(pi_ticks); axes[0].set_xticklabels(pi_labels)

    # Panel B — LOWESS
    axes[1].scatter(pt, tit, c="lightgrey", s=2, alpha=0.25, rasterized=True)
    axes[1].plot(smoothed[:, 0], smoothed[:, 1], color="#d62728", lw=2, label="LOWESS (frac=0.2)")
    axes[1].set_xlabel("Cyclum pseudotime (radians)")
    axes[1].set_ylabel(f"wMel titer")
    axes[1].set_title("B  LOWESS trend")
    axes[1].set_xticks(pi_ticks); axes[1].set_xticklabels(pi_labels)
    axes[1].legend()

    # Panel C — binned median ± IQR
    x_ = bsum["bin_center"].values
    axes[2].plot(x_, bsum["median"].values, "o-", color="#d62728", lw=2, ms=5, label="Median")
    axes[2].fill_between(x_, bsum["q25"].values, bsum["q75"].values,
                         alpha=0.25, color="#d62728", label="IQR")
    axes[2].set_xlabel("Pseudotime bin center (radians)")
    axes[2].set_ylabel(f"wMel titer")
    axes[2].set_title(f"C  Median titer ± IQR across {n_bins} pseudotime bins")
    axes[2].set_xticks(x_)
    axes[2].set_xticklabels([f"{v:.2f}" for v in x_], rotation=45)
    axes[2].legend()

    plt.suptitle("Q1 — wMel titer vs cell-cycle pseudotime", fontweight="bold")
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f"{tag}.pdf"), bbox_inches="tight", dpi=150)
    plt.close()
    print(f"   → {tag}.pdf")


# ─────────────────────────────────────────────────────────────────────────────
# Q2 — Cell-cycle phase distribution vs titer (does titer shift phase?)
# ─────────────────────────────────────────────────────────────────────────────

def q2_phase_distribution_vs_titer(query, fig_dir, sample,
                                    titer_col="wolbachia_titer",
                                    cluster_col="cc_cluster",
                                    n_titer_bins=5):
    """Q2: Does the cell-cycle phase distribution change across titer levels?

    Cells are binned by titer quantile. For each bin, the proportion of cells
    in each cell-cycle cluster is computed (stacked bar). A chi-squared test
    of independence tests whether phase composition is titer-dependent.

    Outputs
    -------
    PDF : q2_phase_dist_vs_titer_{sample}.pdf
    CSV : q2_phase_dist_vs_titer_{sample}_composition.csv
          q2_phase_dist_vs_titer_{sample}_chisq.csv
    """
    tag = f"q2_phase_dist_vs_titer_{sample}"
    obs = query.obs.copy()

    for col in [titer_col, cluster_col]:
        if col not in obs.columns:
            print(f"  Q2 SKIP: '{col}' not in query.obs")
            return

    obs[titer_col]  = pd.to_numeric(obs[titer_col], errors="coerce")
    obs[cluster_col] = obs[cluster_col].astype(str)
    obs = obs.dropna(subset=[titer_col, cluster_col])

    print(f"\n── Q2: phase distribution vs titer  ({len(obs)} cells) ──")

    clusters = sorted(obs[cluster_col].unique(), key=lambda x: int(x) if x.isdigit() else x)
    pal = _cluster_palette(clusters)

    # Quantile-bin titer
    bin_labels = [f"Q{i+1}" for i in range(n_titer_bins)]
    obs["titer_bin"] = pd.qcut(obs[titer_col], q=n_titer_bins,
                                labels=bin_labels, duplicates="drop")
    obs = obs.dropna(subset=["titer_bin"])

    # Contingency table: rows = titer bins, cols = cc clusters
    ct = pd.crosstab(obs["titer_bin"], obs[cluster_col])
    ct = ct.reindex(columns=clusters, fill_value=0)

    # Proportions (row-normalised)
    prop = ct.div(ct.sum(axis=1), axis=0) * 100

    prop.to_csv(os.path.join(fig_dir, f"{tag}_composition.csv"))

    # Chi-squared test of independence
    from scipy.stats import chi2_contingency
    chi2, p_chi, dof, expected = chi2_contingency(ct.values)
    print(f"   Chi-squared: χ²={chi2:.3f}, df={dof}, p={p_chi:.3e}")
    pd.DataFrame({
        "chi2": [chi2], "dof": [dof], "p_value": [p_chi], "n_cells": [len(obs)],
        "note": [f"contingency: {n_titer_bins} titer bins × {len(clusters)} CC clusters"],
    }).to_csv(os.path.join(fig_dir, f"{tag}_chisq.csv"), index=False)

    # ── Figure: 2-panel ───────────────────────────────────────────────────────
    # Panel A: stacked bar — % phase per titer bin
    # Panel B: heatmap of proportions
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Panel A — stacked bar
    bottom = np.zeros(len(prop))
    for cl in clusters:
        vals = prop[cl].values
        axes[0].bar(prop.index.astype(str), vals, bottom=bottom,
                    color=pal[cl], label=f"C{cl}", width=0.7)
        bottom += vals
    axes[0].set_xlabel(f"wMel titer quantile (Q1=lowest, Q{n_titer_bins}=highest)")
    axes[0].set_ylabel("% of cells")
    axes[0].set_title(
        f"A  Cell-cycle phase composition per titer bin\n"
        f"χ²={chi2:.1f}, df={dof}, p={p_chi:.2e}"
    )
    axes[0].legend(title="CC cluster", bbox_to_anchor=(1.01, 1),
                   loc="upper left", fontsize=7)

    # Panel B — heatmap
    sns.heatmap(prop.T, annot=True, fmt=".1f", cmap="YlOrRd",
                linewidths=0.5, ax=axes[1],
                cbar_kws={"label": "% of cells in titer bin"})
    axes[1].set_xlabel("wMel titer bin")
    axes[1].set_ylabel("CC cluster")
    axes[1].set_title("B  Phase proportion heatmap (% per titer bin)")

    plt.suptitle("Q2 — Cell-cycle phase distribution across wMel titer", fontweight="bold")
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f"{tag}.pdf"), bbox_inches="tight", dpi=150)
    plt.close()
    print(f"   → {tag}.pdf")


# ─────────────────────────────────────────────────────────────────────────────
# Q3 — Titer differs by cell-cycle phase (KW + Dunn + violin)
# ─────────────────────────────────────────────────────────────────────────────

def q3_titer_by_phase(query, fig_dir, sample,
                       titer_col="wolbachia_titer",
                       cluster_col="cc_cluster",
                       pseudotime_col="cyclum_pseudotime"):
    """Q3: Do cells in different cell-cycle phases carry different wMel titers?

    Tests
    -----
    - Kruskal-Wallis (omnibus): titer distribution differs across CC clusters?
    - Dunn's post-hoc pairwise (BH FDR): which pairs differ?
    - Spearman: titer vs cluster median titer rank (monotonic summary)
    - Optional: if cyclum_pseudotime present, plots titer on a polar/circular axis
      to reveal whether titer is cyclical with the cell cycle.

    Outputs
    -------
    PDF : q3_titer_by_phase_{sample}.pdf
    CSV : q3_titer_by_phase_{sample}_kruskal.csv
          q3_titer_by_phase_{sample}_dunn.csv
          q3_titer_by_phase_{sample}_cluster_summary.csv
    """
    tag = f"q3_titer_by_phase_{sample}"
    obs = query.obs.copy()

    for col in [titer_col, cluster_col]:
        if col not in obs.columns:
            print(f"  Q3 SKIP: '{col}' not in query.obs")
            return

    obs[titer_col]   = pd.to_numeric(obs[titer_col], errors="coerce")
    obs[cluster_col] = obs[cluster_col].astype(str)
    obs = obs.dropna(subset=[titer_col, cluster_col])

    n_cells  = len(obs)
    clusters = sorted(obs[cluster_col].unique(), key=lambda x: int(x) if x.isdigit() else x)
    pal      = _cluster_palette(clusters)

    print(f"\n── Q3: titer by CC phase  ({n_cells} cells, {len(clusters)} clusters) ──")

    # ── Kruskal-Wallis ────────────────────────────────────────────────────────
    groups = [obs.loc[obs[cluster_col] == cl, titer_col].values for cl in clusters]
    kw_stat, kw_p = kruskal(*groups)
    print(f"   Kruskal-Wallis: H={kw_stat:.3f}, p={kw_p:.3e}")

    pd.DataFrame({
        "statistic": [kw_stat], "p_value": [kw_p],
        "n_clusters": [len(clusters)], "n_cells": [n_cells],
    }).to_csv(os.path.join(fig_dir, f"{tag}_kruskal.csv"), index=False)

    # ── Dunn's post-hoc (manual, BH FDR) ─────────────────────────────────────
    all_ranks  = scipy.stats.rankdata(obs[titer_col].values)
    n_total    = len(all_ranks)
    _, counts  = np.unique(obs[titer_col].values, return_counts=True)
    tie_factor = np.sum(counts ** 3 - counts) / (12 * (n_total - 1)) if n_total > 1 else 0

    obs["_rank"] = all_ranks
    rows = []
    for cl_a, cl_b in combinations(clusters, 2):
        ga = obs.loc[obs[cluster_col] == cl_a, "_rank"].values
        gb = obs.loc[obs[cluster_col] == cl_b, "_rank"].values
        na, nb = len(ga), len(gb)
        if na < 2 or nb < 2:
            continue
        se = np.sqrt((n_total * (n_total + 1) / 12 - tie_factor) * (1 / na + 1 / nb))
        if se == 0:
            continue
        z = (ga.mean() - gb.mean()) / se
        rows.append({"cluster_A": cl_a, "cluster_B": cl_b,
                     "mean_rank_A": ga.mean(), "mean_rank_B": gb.mean(),
                     "z_stat": z, "p_raw": 2 * scipy.stats.norm.sf(abs(z)),
                     "n_A": na, "n_B": nb})

    dunn_df = pd.DataFrame(rows)
    if len(dunn_df):
        _, p_adj, _, _ = multipletests(dunn_df["p_raw"], method="fdr_bh")
        dunn_df["p_adj_BH"] = p_adj
        dunn_df["significant"] = dunn_df["p_adj_BH"] < 0.05
        dunn_df = dunn_df.sort_values("p_adj_BH")
        dunn_df.to_csv(os.path.join(fig_dir, f"{tag}_dunn.csv"), index=False)
        n_sig = dunn_df["significant"].sum()
        print(f"   Dunn post-hoc: {n_sig}/{len(dunn_df)} pairs significant (BH FDR<0.05)")

    # ── Cluster summary ───────────────────────────────────────────────────────
    cluster_median = obs.groupby(cluster_col)[titer_col].median().to_dict()
    summary = (obs.groupby(cluster_col)[titer_col]
               .agg(n_cells="count", median="median", mean="mean", std="std",
                    q25=lambda x: x.quantile(0.25), q75=lambda x: x.quantile(0.75))
               .reset_index().rename(columns={cluster_col: "cluster"}))
    summary.to_csv(os.path.join(fig_dir, f"{tag}_cluster_summary.csv"), index=False)

    # Spearman: titer vs cluster-median-titer (monotonic ordering summary)
    obs["_med"] = obs[cluster_col].map(cluster_median)
    rho, rho_p = spearmanr(obs[titer_col], obs["_med"])
    print(f"   Spearman (titer vs cluster median rank): rho={rho:.3f}, p={rho_p:.3e}")

    # ── Figure: 4-panel ───────────────────────────────────────────────────────
    # Panel A: violin + strip per CC cluster
    # Panel B: UMAP coloured by titer
    # Panel C: UMAP coloured by CC cluster
    # Panel D: polar plot — mean titer by pseudotime phase (if available)
    has_umap = "X_umap" in query.obsm
    has_pt   = pseudotime_col in obs.columns

    ncols = 4 if has_pt else 3
    fig, axes = plt.subplots(1, ncols, figsize=(5 * ncols, 5))

    # Panel A — violin + strip
    ax = axes[0]
    palette = {cl: pal[cl] for cl in clusters}
    sns.violinplot(data=obs, x=cluster_col, y=titer_col, order=clusters,
                   palette=palette, inner=None, linewidth=0.8, ax=ax, cut=0)
    sns.stripplot(data=obs, x=cluster_col, y=titer_col, order=clusters,
                  palette=palette, size=1.5, alpha=0.35, jitter=True, ax=ax)
    for i, cl in enumerate(clusters):
        ax.scatter(i, cluster_median[cl], color="white", s=30, zorder=5,
                   edgecolors="black", linewidths=0.8)
    ax.set_xlabel("Cell-cycle cluster")
    ax.set_ylabel("wMel titer")
    ax.set_title(f"A  Titer by CC phase\nKW p={kw_p:.2e} | ρ={rho:.2f}")
    ax.tick_params(axis="x", rotation=45)

    # Panel B — UMAP titer
    ax = axes[1]
    if has_umap:
        obs_aln  = obs.reindex(query.obs.index)
        tvals    = obs_aln[titer_col].values.astype(float)
        xy       = query.obsm["X_umap"]
        valid    = ~np.isnan(tvals)
        sc_ = ax.scatter(xy[valid, 0], xy[valid, 1], c=tvals[valid],
                         cmap="viridis", s=2, alpha=0.7, rasterized=True)
        ax.scatter(xy[~valid, 0], xy[~valid, 1], c="lightgrey", s=1, alpha=0.2, rasterized=True)
        plt.colorbar(sc_, ax=ax, label="wMel titer", shrink=0.8)
    ax.set_title("B  UMAP — wMel titer")
    ax.set_xticks([]); ax.set_yticks([])

    # Panel C — UMAP CC cluster
    ax = axes[2]
    if has_umap:
        xy       = query.obsm["X_umap"]
        obs_aln  = obs.reindex(query.obs.index)
        cl_vals  = obs_aln[cluster_col].astype(str).values
        int_lbl  = np.array([clusters.index(c) if c in clusters else -1 for c in cl_vals])
        ax.scatter(xy[:, 0], xy[:, 1], c=int_lbl,
                   cmap=matplotlib.colormaps["tab20"].resampled(len(clusters)),
                   s=2, alpha=0.7, rasterized=True)
        for i, cl in enumerate(clusters):
            ax.scatter([], [], color=pal[cl], label=f"C{cl}", s=10)
        ax.legend(fontsize=5, markerscale=2,
                  bbox_to_anchor=(1.01, 1), loc="upper left", title="CC cluster")
    ax.set_title("C  UMAP — CC cluster")
    ax.set_xticks([]); ax.set_yticks([])

    # Panel D — polar plot: is titer cyclical?
    if has_pt:
        ax_pol = fig.add_subplot(1, ncols, ncols, projection="polar")
        obs[pseudotime_col] = pd.to_numeric(obs[pseudotime_col], errors="coerce")
        obs_pt = obs.dropna(subset=[pseudotime_col, titer_col])
        if len(obs_pt) >= 20:
            # Bin pseudotime into equal-width circular sectors
            n_sectors = 12
            edges   = np.linspace(0, 2 * np.pi, n_sectors + 1)
            centers = (edges[:-1] + edges[1:]) / 2
            obs_pt  = obs_pt.copy()
            obs_pt["_sector"] = pd.cut(obs_pt[pseudotime_col], bins=edges,
                                        labels=range(n_sectors), include_lowest=True)
            sector_med = obs_pt.groupby("_sector", observed=True)[titer_col].median()
            sector_med = sector_med.reindex(range(n_sectors)).fillna(0).values

            # Polar bar
            width   = 2 * np.pi / n_sectors
            bars    = ax_pol.bar(centers, sector_med, width=width * 0.85,
                                 bottom=0, alpha=0.8,
                                 color=plt.cm.coolwarm(sector_med / (sector_med.max() + 1e-9)))
            ax_pol.set_theta_zero_location("N")
            ax_pol.set_theta_direction(-1)
            ax_pol.set_xticks(np.linspace(0, 2 * np.pi, 5)[:-1])
            ax_pol.set_xticklabels(["0 (G1?)", "π/2 (S?)", "π (G2?)", "3π/2 (M?)"])
            ax_pol.set_title("D  Titer cyclicity\n(median per pseudotime sector)",
                             pad=15, fontsize=9)
        else:
            ax_pol.set_title("D  Polar: insufficient data", pad=15, fontsize=9)
        # Replace placeholder axis with polar
        axes[3].set_visible(False)

    plt.suptitle("Q3 — wMel titer across cell-cycle phases", fontweight="bold")
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f"{tag}.pdf"), bbox_inches="tight", dpi=150)
    plt.close()
    print(f"   → {tag}.pdf")


# ─────────────────────────────────────────────────────────────────────────────
# Sanitise obs before writing
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
        min_genes, titer_col, pseudotime_col, ref_condition, n_titer_bins):

    os.makedirs(fig_dir, exist_ok=True)
    sc.settings.figdir = fig_dir

    # ── Load ──────────────────────────────────────────────────────────────────
    adata_full = sc.read_h5ad(query_path)
    ref        = sc.read_h5ad(ref_path)

    print(f"Loaded query object : {adata_full.n_obs} cells")
    print(f"Loaded reference    : {ref.n_obs} cells, {ref.obs['leiden'].nunique()} clusters")
    print(f"Query obs columns   : {list(adata_full.obs.columns)}")

    # ── Identify infected / SV-timepoint cells ────────────────────────────────
    ref_conditions = [ref_condition] if isinstance(ref_condition, str) else ref_condition
    query_mask = (
        ~adata_full.obs["cell_line"].isin(ref_conditions) |
        (adata_full.obs["treatment"] != "Ctrl")
    )
    print(f"\nReference conditions : {ref_conditions} (Ctrl only)")
    print(f"Query bio_conditions :\n{adata_full.obs.loc[query_mask, 'bio_condition'].value_counts()}")
    print(f"Query cells          : {query_mask.sum()}")

    if query_mask.sum() == 0:
        raise ValueError(
            "No query cells found. Check --ref_condition and that your query h5ad "
            "contains infected/SV-timepoint cells with treatment != 'Ctrl'.\n"
            f"Available cell_lines : {adata_full.obs['cell_line'].unique().tolist()}\n"
            f"Available treatments : {adata_full.obs['treatment'].unique().tolist()}"
        )

    adata_query = adata_full[query_mask].copy()

    # Verify titer column exists on query
    if titer_col not in adata_query.obs.columns:
        raise ValueError(
            f"'{titer_col}' not found in query obs. "
            f"Available columns: {list(adata_query.obs.columns)}"
        )
    n_with_titer = pd.to_numeric(adata_query.obs[titer_col], errors="coerce").notna().sum()
    print(f"\nQuery cells with non-null {titer_col}: {n_with_titer}/{adata_query.n_obs}")

    # ── Ingest ────────────────────────────────────────────────────────────────
    query = ingest_query(adata_query, ref=ref, min_genes=min_genes,
                         fig_dir=fig_dir, sample=sample)

    # ── Three biological questions ────────────────────────────────────────────
    q1_titer_vs_pseudotime(query, fig_dir, sample,
                           titer_col=titer_col, pseudotime_col=pseudotime_col,
                           cluster_col="cc_cluster")

    q2_phase_distribution_vs_titer(query, fig_dir, sample,
                                   titer_col=titer_col, cluster_col="cc_cluster",
                                   n_titer_bins=n_titer_bins)

    q3_titer_by_phase(query, fig_dir, sample,
                      titer_col=titer_col, cluster_col="cc_cluster",
                      pseudotime_col=pseudotime_col)

    # ── Save outputs ──────────────────────────────────────────────────────────
    query = _sanitize_obs(query)
    ref   = _sanitize_obs(ref)

    query_out = out_path.replace(".h5ad", "_query.h5ad")
    ref_out   = out_path.replace(".h5ad", "_reference.h5ad")
    query.write(query_out)
    ref.write(ref_out)

    print("\n" + "=" * 60)
    print("COMPLETE")
    print("=" * 60)
    print(f"Query   → {query_out}")
    print(f"Reference → {ref_out}")
    print(f"Figures → {fig_dir}/")


def main():
    parser = argparse.ArgumentParser(
        description="Test wMel titer vs cell-cycle phase via reference-anchored ingest"
    )
    parser.add_argument("--ref",            required=True,
                        help="Pre-built reference h5ad (uninfected, Leiden-clustered by CC phase)")
    parser.add_argument("--query",          required=True,
                        help="Full integrated h5ad containing infected/SV-timepoint cells")
    parser.add_argument("--out_path",       default="results/integrated_by_reftiter_cellcycle.h5ad")
    parser.add_argument("--fig_dir",        default="results/integrated_by_ref/figures/titer_vs_cellcycle")
    parser.add_argument("--sample",         default="wolbachia_infection")
    parser.add_argument("--min_genes",      type=int,   default=200)
    parser.add_argument("--titer_col",      default="wolbachia_titer")
    parser.add_argument("--pseudotime_col", default="cyclum_pseudotime")
    parser.add_argument("--ref_condition",  nargs="+",  default=["JW18DOX"],
                        help="cell_line value(s) treated as uninfected reference "
                             "(these + treatment=Ctrl define what is NOT query)")
    parser.add_argument("--n_titer_bins",   type=int,   default=5,
                        help="Number of quantile bins for Q2 titer binning")
    args = parser.parse_args()

    run(
        ref_path=args.ref,
        query_path=args.query,
        out_path=args.out_path,
        fig_dir=args.fig_dir,
        sample=args.sample,
        min_genes=args.min_genes,
        titer_col=args.titer_col,
        pseudotime_col=args.pseudotime_col,
        ref_condition=args.ref_condition,
        n_titer_bins=args.n_titer_bins,
    )


if __name__ == "__main__":
    main()