"""
integrate_by_ref.py
============
Staged scRNA-seq integration pipeline for studying Wolbachia infection dynamics.

Strategy
--------
1. Load a pre-built reference (Harmony-corrected, leiden-clustered, UMAP embedded).
2. Subset the full integrated object to query cells (new-infection timepoints).
3. Map (ingest) query cells onto the reference.
4. Project cluster labels and UMAP coordinates to query cells.

Run with:
mamba activate scanpy 
python scripts/method_comparison/integrate_by_ref.py \
    --ref results/integrated/integrated_uninfected_with_cellcycle.h5ad \
    --query results/integrated/integrated_uninfected.h5ad \
    --out_path results/integrated/integrated_by_cellcycle.h5ad 
    

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
from scipy.stats import kruskal, spearmanr
from itertools import combinations
from statsmodels.stats.multitest import multipletests
import anndata as ad
import scanpy as sc

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _leiden_colors(adata, key="leiden"):
    clusters = sorted(adata.obs[key].unique())
    cmap = matplotlib.colormaps["tab20"]
    return [cmap(i % 20) for i in range(len(clusters))]


# ─────────────────────────────────────────────────────────────────────────────
# Query: map onto reference via ingest
# ─────────────────────────────────────────────────────────────────────────────

def map_query_to_reference(adata_query, ref, min_genes, fig_dir, sample):
    """Project new-infection timepoint cells onto the reference using sc.tl.ingest.

    Parameters
    ----------
    adata_query : AnnData
        Query cells (new-infection timepoints), raw counts in .X.
    ref : AnnData
        Pre-built reference: Harmony-corrected, scaled, PCA/UMAP computed,
        leiden labels in obs["leiden"].
    min_genes : int
        Minimum genes per cell for quality filtering.
    fig_dir : str
        Directory for output figures.
    sample : str
        Sample label used in output file names.

    Returns
    -------
    query : AnnData
        Query object with projected UMAP and transferred leiden_ref labels.
    combined : AnnData
        Concatenation of reference and query for joint plotting.
    """
    query = adata_query.copy()
    print(f"Query cells: {query.n_obs}")
    print(f"Timepoints:  {query.obs['timepoint'].value_counts().to_dict()}")
    print(f"Methods:     {query.obs['method'].value_counts().to_dict()}")

    # Subset to the gene universe present in the reference
    ref_all_genes = ref.raw.var_names if ref.raw is not None else ref.var_names
    query = query[:, query.var_names.isin(ref_all_genes)].copy()

    # Remove explicitly stored zeros before filtering
    if scipy.sparse.issparse(query.X):
        query.X.eliminate_zeros()

    sc.pp.filter_cells(query, min_genes=min_genes)
    sc.pp.filter_cells(query, min_counts=1)
    sc.pp.filter_genes(query, min_counts=1)

    if scipy.sparse.issparse(query.X):
        query.X = query.X.toarray()
    query.X = np.nan_to_num(query.X.astype(np.float64), nan=0.0, posinf=0.0, neginf=0.0)

    # Remove zero-count cells/genes after dense conversion
    cell_sums = query.X.sum(axis=1)
    gene_sums = query.X.sum(axis=0)
    query = query[cell_sums > 0].copy()
    query = query[:, gene_sums > 0].copy()

    query.X = np.nan_to_num(query.X, nan=0.0, posinf=0.0, neginf=0.0)

    # Store pre-scale counts as .raw
    query.raw = query

    # Fill missing reference HVGs with zeros, then reorder vars to match ref exactly.
    # sc.tl.ingest requires identical var_names in the same order.
    missing_genes = ref.var_names[~ref.var_names.isin(query.var_names)]
    if len(missing_genes) > 0:
        print(f"  WARNING: {len(missing_genes)} ref HVGs missing from query — filling with zeros")
        zero_cols = np.zeros((query.n_obs, len(missing_genes)), dtype=np.float64)
        combined_X = np.hstack([query.X, zero_cols])
        combined_var = pd.DataFrame(index=list(query.var_names) + list(missing_genes))
        query = ad.AnnData(X=combined_X, obs=query.obs.copy(), var=combined_var)

    query = query[:, ref.var_names].copy()
    sc.pp.scale(query, max_value=10)
    query.var["highly_variable"] = True

    # ingest projects query into ref PCA + UMAP and transfers leiden labels.
    # Save obs beforehand because sc.tl.ingest overwrites it, then merge back.
    obs_backup = query.obs.copy()
    print("Running sc.tl.ingest …")
    sc.tl.ingest(query, ref, obs="leiden", embedding_method="umap")
    leiden_transferred = query.obs["leiden"].copy()
    query.obs = obs_backup
    query.obs["leiden_ref"] = leiden_transferred.values

    sc.pl.umap(query, color=["leiden_ref", "timepoint", "method"],
               save=f"_{sample}_query_ingested.pdf", ncols=3,
               title=["Cluster (transferred from ref)", "Timepoint", "Method"])

    # Build joint reference + query object for plotting
    ref_plot = ref.copy()
    ref_plot.obs["dataset"]    = "uninfected (reference)"
    ref_plot.obs["leiden_ref"] = ref_plot.obs["leiden"]
    ref_plot.obs["timepoint"]  = ref_plot.obs["timepoint"].astype(str).fillna("uninfected")
    query.obs["timepoint"]     = query.obs["timepoint"].astype(str)
    query.obs["dataset"]       = "new infection (" + query.obs["timepoint"].fillna("?") + ")"

    common_genes = ref_plot.var_names.intersection(query.var_names)
    combined = ad.concat(
        [ref_plot[:, common_genes], query[:, common_genes]],
        join="inner", index_unique="-",
    )

    sc.pl.umap(combined, color=["dataset", "leiden_ref", "timepoint", "method"],
               save=f"_{sample}_combined_ref_query.pdf", ncols=2,
               title=["Dataset (ref vs query)", "Cluster", "Timepoint", "Method"])

    _plot_cluster_by_timepoint(combined, fig_dir, sample)
    _plot_method_validation(ref, query, fig_dir, sample)

    return query, combined


# ─────────────────────────────────────────────────────────────────────────────
# Plotting helpers
# ─────────────────────────────────────────────────────────────────────────────

def _plot_cluster_by_timepoint(combined, fig_dir, sample):
    cluster_col = "leiden_ref" if "leiden_ref" in combined.obs.columns else "leiden"
    colors = _leiden_colors(combined, key=cluster_col)
    timepoints = sorted(combined.obs["timepoint_numeric"].unique())

    comp = pd.crosstab(combined.obs[cluster_col],
                       combined.obs["timepoint_numeric"],
                       normalize="columns") * 100
    fig, ax = plt.subplots(figsize=(10, 6))
    comp.T.plot(kind="bar", stacked=True, ax=ax, color=colors, width=0.8)
    ax.set_xlabel("Timepoint (days post-infection; 0 = uninfected)")
    ax.set_ylabel("% of cells")
    ax.set_title("Cluster composition across infection timepoints")
    ax.legend(title="Cluster", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.xticks(rotation=0)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f"cluster_by_timepoint_{sample}.pdf"), bbox_inches="tight")
    plt.close()

    n_tp = len(timepoints)
    ncols = min(4, n_tp)
    nrows = int(np.ceil(n_tp / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4, nrows * 3.5))
    axes = np.array(axes).flatten()
    umap_xy = combined.obsm["X_umap"]
    all_labels = combined.obs[cluster_col].astype(int).values
    n_cl = len(np.unique(all_labels))

    for i, tp in enumerate(timepoints):
        ax = axes[i]
        mask = combined.obs["timepoint_numeric"].values == tp
        ax.scatter(umap_xy[~mask, 0], umap_xy[~mask, 1],
                   c="lightgrey", s=1, alpha=0.3, rasterized=True)
        ax.scatter(umap_xy[mask, 0], umap_xy[mask, 1],
                   c=all_labels[mask], cmap=matplotlib.colormaps["tab20"].resampled(n_cl),
                   s=2, alpha=0.8, rasterized=True)
        ax.set_title("uninfected (ref)" if tp == 0 else f"D{tp}", fontsize=10)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)
    plt.suptitle("UMAP — cells highlighted per timepoint", fontweight="bold")
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f"umap_by_timepoint_{sample}.pdf"),
                bbox_inches="tight", dpi=150)
    plt.close()

    if "wolbachia_titer" in combined.obs.columns:
        fig, ax = plt.subplots(figsize=(12, 5))
        sns.boxplot(data=combined.obs, x="timepoint_numeric", y="wolbachia_titer",
                    hue=cluster_col, ax=ax, flierprops=dict(markersize=1))
        ax.set_xlabel("Timepoint (days post-infection)")
        ax.set_ylabel("Wolbachia Titer")
        ax.set_title("Wolbachia titer by timepoint and cluster")
        ax.legend(title="Cluster", bbox_to_anchor=(1.05, 1), loc="upper left")
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f"titer_by_timepoint_{sample}.pdf"), bbox_inches="tight")
        plt.close()


def _plot_method_validation(ref, query, fig_dir, sample):
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))

    for row, (adata, title_prefix) in enumerate([(ref, "Reference"), (query, "Query")]):
        umap_xy = adata.obsm["X_umap"]
        cluster_col = "leiden" if "leiden" in adata.obs.columns else "leiden_ref"

        method_cmap = {"10x": "#1f77b4", "pipseq": "#ff7f0e"}
        c = [method_cmap.get(m, "grey") for m in adata.obs["method"].values]
        axes[row, 0].scatter(umap_xy[:, 0], umap_xy[:, 1], c=c, s=1, alpha=0.5, rasterized=True)
        axes[row, 0].set_title(f"{title_prefix}: Method (should overlap post-correction)")
        for m, col in method_cmap.items():
            axes[row, 0].scatter([], [], c=col, label=m)
        axes[row, 0].legend(); axes[row, 0].set_xticks([]); axes[row, 0].set_yticks([])

        conditions = adata.obs["bio_condition"].astype("category")
        cond_cmap = matplotlib.colormaps["Set2"].resampled(len(conditions.cat.categories))
        c = [cond_cmap(conditions.cat.codes.iloc[i]) for i in range(len(conditions))]
        axes[row, 1].scatter(umap_xy[:, 0], umap_xy[:, 1], c=c, s=1, alpha=0.5, rasterized=True)
        axes[row, 1].set_title(f"{title_prefix}: Bio condition (should still separate)")
        for i, cat in enumerate(conditions.cat.categories):
            axes[row, 1].scatter([], [], color=cond_cmap(i), label=cat)
        axes[row, 1].legend(fontsize=7); axes[row, 1].set_xticks([]); axes[row, 1].set_yticks([])

        lbls = adata.obs[cluster_col].astype(int).values
        axes[row, 2].scatter(umap_xy[:, 0], umap_xy[:, 1], c=lbls,
                             cmap=matplotlib.colormaps["tab20"].resampled(len(np.unique(lbls))),
                             s=1, alpha=0.5, rasterized=True)
        axes[row, 2].set_title(f"{title_prefix}: Clusters")
        axes[row, 2].set_xticks([]); axes[row, 2].set_yticks([])

    plt.suptitle("Method correction validation\n"
                 "10X and PIPseq should overlap; biological conditions should separate",
                 fontweight="bold")
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f"method_validation_{sample}.pdf"),
                bbox_inches="tight", dpi=150)
    plt.close()


# ─────────────────────────────────────────────────────────────────────────────
# Titer vs cell-cycle cluster analysis
# ─────────────────────────────────────────────────────────────────────────────

def analyze_titer_vs_clusters(query, fig_dir, sample, titer_col="wolbachia_titer",
                               cluster_col="leiden_ref", n_titer_bins=4):
    """Test whether wMel titer differs across cell-cycle clusters (query cells only).

    Statistical approach
    --------------------
    1. Kruskal-Wallis (omnibus): does titer distribution differ across clusters?
    2. Dunn's post-hoc pairwise test with Benjamini-Hochberg FDR correction.
    3. Spearman correlation between titer and cluster identity encoded as
       median titer rank — summarises the monotonic titer-vs-cluster relationship.

    Outputs
    -------
    CSV  : titer_vs_clusters_{sample}_kruskal.csv
           titer_vs_clusters_{sample}_dunn_posthoc.csv
    PDF  : titer_vs_clusters_{sample}_violin.pdf
           titer_vs_clusters_{sample}_umap_titer.pdf
           titer_vs_clusters_{sample}_titer_bin_composition.pdf

    Parameters
    ----------
    query : AnnData
        Query AnnData with obs[titer_col] and obs[cluster_col].
    fig_dir : str
        Output directory for figures and CSV files.
    sample : str
        Label used in output file names.
    titer_col : str
        obs column containing numeric wMel titer values.
    cluster_col : str
        obs column containing cluster labels (leiden_ref by default).
    n_titer_bins : int
        Number of quantile bins for the titer-bin composition plot.
    """
    obs = query.obs.copy()

    # ── Validate inputs ───────────────────────────────────────────────────────
    if titer_col not in obs.columns:
        print(f"  WARNING: '{titer_col}' not found in query.obs — skipping titer analysis.")
        return
    if cluster_col not in obs.columns:
        print(f"  WARNING: '{cluster_col}' not found in query.obs — skipping titer analysis.")
        return

    obs[titer_col]   = pd.to_numeric(obs[titer_col], errors="coerce")
    obs[cluster_col] = obs[cluster_col].astype(str)
    obs = obs.dropna(subset=[titer_col, cluster_col])

    n_cells = len(obs)
    clusters = sorted(obs[cluster_col].unique(), key=lambda x: int(x) if x.isdigit() else x)
    n_clusters = len(clusters)
    print(f"\nTiter vs clusters: {n_cells} cells, {n_clusters} clusters")

    # ── 1. Kruskal-Wallis omnibus test ────────────────────────────────────────
    groups = [obs.loc[obs[cluster_col] == cl, titer_col].values for cl in clusters]
    kw_stat, kw_p = kruskal(*groups)
    print(f"  Kruskal-Wallis: H={kw_stat:.3f}, p={kw_p:.3e}")

    kw_df = pd.DataFrame({
        "test":      ["Kruskal-Wallis"],
        "statistic": [kw_stat],
        "p_value":   [kw_p],
        "n_clusters":[n_clusters],
        "n_cells":   [n_cells],
    })
    kw_df.to_csv(os.path.join(fig_dir, f"titer_vs_clusters_{sample}_kruskal.csv"), index=False)

    # ── 2. Dunn's post-hoc pairwise tests (BH FDR) ───────────────────────────
    # Dunn's uses the pooled rank sum from the KW test statistic.
    # We compute it manually to avoid adding scikit-posthocs as a dependency.
    all_vals   = obs[titer_col].values
    all_ranks  = scipy.stats.rankdata(all_vals)
    n_total    = len(all_ranks)

    # Tie correction factor
    _, counts  = np.unique(all_vals, return_counts=True)
    tie_factor = np.sum(counts ** 3 - counts) / (12 * (n_total - 1)) if n_total > 1 else 0

    rank_lookup = dict(zip(obs.index, all_ranks))
    obs["_rank"] = obs.index.map(rank_lookup)

    rows = []
    for cl_a, cl_b in combinations(clusters, 2):
        grp_a = obs.loc[obs[cluster_col] == cl_a, "_rank"].values
        grp_b = obs.loc[obs[cluster_col] == cl_b, "_rank"].values
        na, nb = len(grp_a), len(grp_b)
        if na < 2 or nb < 2:
            continue
        mean_rank_a = grp_a.mean()
        mean_rank_b = grp_b.mean()
        se = np.sqrt(
            (n_total * (n_total + 1) / 12 - tie_factor) * (1 / na + 1 / nb)
        )
        if se == 0:
            continue
        z = (mean_rank_a - mean_rank_b) / se
        p = 2 * scipy.stats.norm.sf(abs(z))
        rows.append({
            "cluster_A":     cl_a,
            "cluster_B":     cl_b,
            "mean_rank_A":   mean_rank_a,
            "mean_rank_B":   mean_rank_b,
            "z_stat":        z,
            "p_value_raw":   p,
            "n_A":           na,
            "n_B":           nb,
        })

    dunn_df = pd.DataFrame(rows)
    if len(dunn_df) > 0:
        _, p_adj, _, _ = multipletests(dunn_df["p_value_raw"], method="fdr_bh")
        dunn_df["p_adj_BH"] = p_adj
        dunn_df["significant"] = dunn_df["p_adj_BH"] < 0.05
        dunn_df = dunn_df.sort_values("p_adj_BH")
        dunn_df.to_csv(
            os.path.join(fig_dir, f"titer_vs_clusters_{sample}_dunn_posthoc.csv"), index=False
        )
        n_sig = dunn_df["significant"].sum()
        print(f"  Dunn post-hoc: {n_sig}/{len(dunn_df)} pairs significant (BH FDR < 0.05)")

    # ── 3. Spearman: titer vs median-rank of cluster ──────────────────────────
    # Each cell is assigned the median titer of its cluster; this gives a
    # cluster-level ordering for the correlation.
    cluster_median = obs.groupby(cluster_col)[titer_col].median().to_dict()
    obs["_cluster_median_titer"] = obs[cluster_col].map(cluster_median)
    rho, spear_p = spearmanr(obs[titer_col], obs["_cluster_median_titer"])
    print(f"  Spearman (titer vs cluster median titer rank): rho={rho:.3f}, p={spear_p:.3e}")

    spear_row = pd.DataFrame({
        "rho":     [rho],
        "p_value": [spear_p],
        "n_cells": [n_cells],
    })
    spear_row.to_csv(
        os.path.join(fig_dir, f"titer_vs_clusters_{sample}_spearman.csv"), index=False
    )

    # ── Cluster summary table (appended to KW CSV as a second sheet is not
    #    straightforward in plain CSV; write separately instead) ───────────────
    summary = (
        obs.groupby(cluster_col)[titer_col]
        .agg(n_cells="count", median="median", mean="mean", std="std",
             q25=lambda x: x.quantile(0.25), q75=lambda x: x.quantile(0.75))
        .reset_index()
        .rename(columns={cluster_col: "cluster"})
    )
    summary.to_csv(
        os.path.join(fig_dir, f"titer_vs_clusters_{sample}_cluster_summary.csv"), index=False
    )

    # ── Plots ─────────────────────────────────────────────────────────────────
    leiden_cmap = matplotlib.colormaps["tab20"]
    cluster_colors = {cl: leiden_cmap(i % 20) for i, cl in enumerate(clusters)}

    # --- Plot 1: Violin + strip per cluster -----------------------------------
    fig, ax = plt.subplots(figsize=(max(8, n_clusters * 0.9), 5))
    palette = {cl: cluster_colors[cl] for cl in clusters}
    sns.violinplot(data=obs, x=cluster_col, y=titer_col, order=clusters,
                   palette=palette, inner=None, linewidth=0.8, ax=ax, cut=0)
    sns.stripplot(data=obs, x=cluster_col, y=titer_col, order=clusters,
                  palette=palette, size=1.5, alpha=0.4, jitter=True, ax=ax)
    # Annotate median
    for i, cl in enumerate(clusters):
        med = cluster_median[cl]
        ax.scatter(i, med, color="white", s=30, zorder=5, edgecolors="black", linewidths=0.8)
    ax.set_xlabel("Cluster (cell-cycle phase proxy)")
    ax.set_ylabel("wMel titer")
    ax.set_title(
        f"wMel titer per cell-cycle cluster  |  "
        f"KW p={kw_p:.2e}  |  Spearman rho={rho:.2f}"
    )
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f"titer_vs_clusters_{sample}_violin.pdf"),
                bbox_inches="tight")
    plt.close()

    # --- Plot 2: UMAP colored by titer ----------------------------------------
    if "X_umap" in query.obsm:
        # Re-index obs to align with query cell order
        obs_aligned = obs.reindex(query.obs.index)
        titer_vals  = obs_aligned[titer_col].values.astype(float)
        umap_xy     = query.obsm["X_umap"]

        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        # Left: titer continuous
        valid = ~np.isnan(titer_vals)
        sc_ = axes[0].scatter(
            umap_xy[valid, 0], umap_xy[valid, 1],
            c=titer_vals[valid], cmap="viridis", s=2, alpha=0.7, rasterized=True,
        )
        axes[0].scatter(
            umap_xy[~valid, 0], umap_xy[~valid, 1],
            c="lightgrey", s=1, alpha=0.3, rasterized=True,
        )
        plt.colorbar(sc_, ax=axes[0], label="wMel titer")
        axes[0].set_title("wMel titer (continuous)")
        axes[0].set_xticks([]); axes[0].set_yticks([])

        # Right: cluster labels
        cluster_vals = obs_aligned[cluster_col].values
        int_labels   = np.array([
            clusters.index(str(c)) if str(c) in clusters else -1
            for c in cluster_vals
        ])
        sc2 = axes[1].scatter(
            umap_xy[:, 0], umap_xy[:, 1],
            c=int_labels, cmap=matplotlib.colormaps["tab20"].resampled(n_clusters),
            s=2, alpha=0.7, rasterized=True,
        )
        axes[1].set_title("Cell-cycle clusters (transferred from ref)")
        axes[1].set_xticks([]); axes[1].set_yticks([])
        for i, cl in enumerate(clusters):
            axes[1].scatter([], [], color=leiden_cmap(i % 20), label=f"Cluster {cl}", s=10)
        axes[1].legend(fontsize=6, markerscale=2, bbox_to_anchor=(1.05, 1), loc="upper left")

        plt.suptitle("Query cells: wMel titer on UMAP", fontweight="bold")
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f"titer_vs_clusters_{sample}_umap_titer.pdf"),
                    bbox_inches="tight", dpi=150)
        plt.close()

    # --- Plot 3: Titer-bin stacked bar of cluster composition -----------------
    valid_obs = obs.dropna(subset=[titer_col])
    if len(valid_obs) > n_titer_bins * 10:
        bin_labels = [f"Q{i+1}" for i in range(n_titer_bins)]
        valid_obs = valid_obs.copy()
        valid_obs["titer_bin"] = pd.qcut(
            valid_obs[titer_col], q=n_titer_bins, labels=bin_labels, duplicates="drop"
        )
        comp = pd.crosstab(valid_obs["titer_bin"], valid_obs[cluster_col], normalize="index") * 100
        comp = comp.reindex(columns=clusters, fill_value=0)

        fig, ax = plt.subplots(figsize=(7, 5))
        bottom = np.zeros(len(comp))
        for cl in clusters:
            vals = comp[cl].values
            ax.bar(comp.index.astype(str), vals, bottom=bottom,
                   color=cluster_colors[cl], label=f"Cluster {cl}", width=0.7)
            bottom += vals
        ax.set_xlabel(f"wMel titer quantile bin  (Q1=lowest, Q{n_titer_bins}=highest)")
        ax.set_ylabel("% of cells")
        ax.set_title("Cell-cycle cluster composition across titer bins")
        ax.legend(title="Cluster", bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=7)
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f"titer_vs_clusters_{sample}_titer_bin_composition.pdf"),
                    bbox_inches="tight")
        plt.close()

    print(f"  Titer analysis outputs written to {fig_dir}/")


# ─────────────────────────────────────────────────────────────────────────────
# Titer vs Cyclum pseudotime analysis
# ─────────────────────────────────────────────────────────────────────────────

def analyze_titer_vs_pseudotime(adata, fig_dir, sample, dataset_label,
                                 titer_col="wolbachia_titer",
                                 pseudotime_col="cyclum_pseudotime",
                                 cluster_col="leiden_ref",
                                 n_pseudotime_bins=8):
    """Plot and test the relationship between wMel titer and Cyclum cell-cycle pseudotime.

    Cyclum pseudotime is circular (0–2π), representing position in the cell cycle.
    Titer is treated as a continuous variable. The analysis tests whether titer
    varies systematically across the cell cycle using a circular-aware approach:
    pseudotime is binned into equal-width windows and titer is summarised per bin.

    Statistical approach
    --------------------
    Spearman correlation between titer and pseudotime (linear approximation;
    interpret with caution given circularity — flagged in output).

    Outputs (per dataset_label)
    ---------------------------
    PDF : titer_vs_pseudotime_{sample}_{dataset_label}_scatter.pdf
          titer_vs_pseudotime_{sample}_{dataset_label}_pseudotime_bins.pdf
    CSV : titer_vs_pseudotime_{sample}_{dataset_label}_spearman.csv
          titer_vs_pseudotime_{sample}_{dataset_label}_bin_summary.csv

    Parameters
    ----------
    adata : AnnData
        AnnData object (ref or query) with obs[titer_col] and obs[pseudotime_col].
    fig_dir : str
        Output directory.
    sample : str
        Sample label for file names.
    dataset_label : str
        Short label distinguishing ref vs query (e.g. "reference", "query").
    titer_col : str
        obs column for wMel titer.
    pseudotime_col : str
        obs column for Cyclum pseudotime (expected range 0–2π).
    cluster_col : str
        obs column for cluster labels; used to color scatter points.
    n_pseudotime_bins : int
        Number of equal-width bins across the pseudotime range.
    """
    obs = adata.obs.copy()
    tag = f"{sample}_{dataset_label}"

    # ── Validate ──────────────────────────────────────────────────────────────
    missing = [c for c in [titer_col, pseudotime_col] if c not in obs.columns]
    if missing:
        print(f"  WARNING: {missing} not found in {dataset_label}.obs — skipping pseudotime analysis.")
        return

    obs[titer_col]       = pd.to_numeric(obs[titer_col],       errors="coerce")
    obs[pseudotime_col]  = pd.to_numeric(obs[pseudotime_col],  errors="coerce")
    obs = obs.dropna(subset=[titer_col, pseudotime_col])

    if len(obs) < 20:
        print(f"  WARNING: fewer than 20 cells with both titer and pseudotime in {dataset_label} — skipping.")
        return

    print(f"\nTiter vs pseudotime ({dataset_label}): {len(obs)} cells")

    has_clusters = cluster_col in obs.columns
    if has_clusters:
        obs[cluster_col] = obs[cluster_col].astype(str)
        clusters = sorted(obs[cluster_col].unique(), key=lambda x: int(x) if x.isdigit() else x)
        leiden_cmap = matplotlib.colormaps["tab20"]
        cluster_color_map = {cl: leiden_cmap(i % 20) for i, cl in enumerate(clusters)}
        point_colors = [cluster_color_map.get(c, "grey") for c in obs[cluster_col]]
    else:
        point_colors = obs[titer_col].values  # fallback: color by titer

    pt  = obs[pseudotime_col].values
    tit = obs[titer_col].values

    # ── Spearman (linear approximation — circularity caveat noted) ────────────
    rho, p_val = spearmanr(pt, tit)
    print(f"  Spearman rho={rho:.3f}, p={p_val:.3e}  "
          f"(NOTE: pseudotime is circular; interpret linear correlation with caution)")

    spear_df = pd.DataFrame({
        "dataset":          [dataset_label],
        "spearman_rho":     [rho],
        "p_value":          [p_val],
        "n_cells":          [len(obs)],
        "note":             ["pseudotime is circular (0-2pi); linear Spearman is approximate"],
    })
    spear_df.to_csv(
        os.path.join(fig_dir, f"titer_vs_pseudotime_{tag}_spearman.csv"), index=False
    )

    # ── Plot 1: Scatter — pseudotime vs titer, colored by cluster ─────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Left: raw scatter
    if has_clusters:
        for cl in clusters:
            mask = obs[cluster_col] == cl
            axes[0].scatter(
                pt[mask.values], tit[mask.values],
                c=[cluster_color_map[cl]], s=3, alpha=0.5, label=f"Cluster {cl}",
                rasterized=True,
            )
        axes[0].legend(fontsize=6, markerscale=2, bbox_to_anchor=(1.01, 1),
                       loc="upper left", title="Cluster")
    else:
        sc_ = axes[0].scatter(pt, tit, c=tit, cmap="viridis", s=3, alpha=0.5, rasterized=True)
        plt.colorbar(sc_, ax=axes[0], label="wMel titer")

    axes[0].set_xlabel("Cyclum pseudotime (radians, 0–2π)")
    axes[0].set_ylabel("wMel titer")
    axes[0].set_title(f"{dataset_label}: titer vs pseudotime\nSpearman rho={rho:.3f}, p={p_val:.2e}")

    # Right: LOWESS smoothed trend + scatter (grey)
    from statsmodels.nonparametric.smoothers_lowess import lowess
    order     = np.argsort(pt)
    pt_sort   = pt[order]
    tit_sort  = tit[order]
    smoothed  = lowess(tit_sort, pt_sort, frac=0.2, return_sorted=True)

    axes[1].scatter(pt, tit, c="lightgrey", s=2, alpha=0.3, rasterized=True)
    axes[1].plot(smoothed[:, 0], smoothed[:, 1], color="#d62728", linewidth=2,
                 label="LOWESS (frac=0.2)")
    axes[1].set_xlabel("Cyclum pseudotime (radians, 0–2π)")
    axes[1].set_ylabel("wMel titer")
    axes[1].set_title(f"{dataset_label}: LOWESS trend")
    axes[1].legend()

    # Add π tick labels on x-axis for both panels
    for ax in axes[:2]:
        ax.set_xticks([0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi])
        ax.set_xticklabels(["0", "π/2", "π", "3π/2", "2π"])

    plt.suptitle(f"wMel titer vs Cyclum cell-cycle pseudotime  [{dataset_label}]",
                 fontweight="bold")
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f"titer_vs_pseudotime_{tag}_scatter.pdf"),
                bbox_inches="tight", dpi=150)
    plt.close()

    # ── Plot 2: Titer distribution binned by pseudotime windows ───────────────
    pt_min = obs[pseudotime_col].min()
    pt_max = obs[pseudotime_col].max()
    bin_edges  = np.linspace(pt_min, pt_max, n_pseudotime_bins + 1)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    bin_labels  = [f"{e:.2f}" for e in bin_edges[:-1]]

    obs["_pt_bin"] = pd.cut(
        obs[pseudotime_col], bins=bin_edges, labels=bin_labels, include_lowest=True
    )
    obs_binned = obs.dropna(subset=["_pt_bin"])

    # Per-bin summary for CSV
    bin_summary = (
        obs_binned.groupby("_pt_bin", observed=True)[titer_col]
        .agg(n_cells="count", median="median", mean="mean", std="std",
             q25=lambda x: x.quantile(0.25), q75=lambda x: x.quantile(0.75))
        .reset_index()
        .rename(columns={"_pt_bin": "pseudotime_bin_start"})
    )
    bin_summary["bin_center_radians"] = bin_centers
    bin_summary.to_csv(
        os.path.join(fig_dir, f"titer_vs_pseudotime_{tag}_bin_summary.csv"), index=False
    )

    fig, axes = plt.subplots(1, 2, figsize=(16, 5))

    # Left: violin per bin
    bin_order = [str(b) for b in bin_labels]
    palette   = sns.color_palette("coolwarm", n_colors=n_pseudotime_bins)
    sns.violinplot(data=obs_binned, x="_pt_bin", y=titer_col, order=bin_order,
                   palette=palette, inner=None, linewidth=0.8, ax=axes[0], cut=0)
    sns.stripplot(data=obs_binned, x="_pt_bin", y=titer_col, order=bin_order,
                  color="black", size=1, alpha=0.3, jitter=True, ax=axes[0])
    # Median dots
    for i, bl in enumerate(bin_labels):
        med = obs_binned.loc[obs_binned["_pt_bin"] == bl, titer_col].median()
        axes[0].scatter(i, med, color="white", s=25, zorder=5,
                        edgecolors="black", linewidths=0.8)
    axes[0].set_xlabel("Pseudotime bin (start, radians)")
    axes[0].set_ylabel("wMel titer")
    axes[0].set_title(f"{dataset_label}: titer per pseudotime bin")
    axes[0].tick_params(axis="x", rotation=45)

    # Right: median titer across pseudotime bins (line plot with IQR ribbon)
    medians = bin_summary["median"].values
    q25_    = bin_summary["q25"].values
    q75_    = bin_summary["q75"].values
    x_      = bin_summary["bin_center_radians"].values

    axes[1].plot(x_, medians, "o-", color="#d62728", linewidth=2, markersize=5,
                 label="Median titer")
    axes[1].fill_between(x_, q25_, q75_, alpha=0.25, color="#d62728", label="IQR")
    axes[1].set_xlabel("Pseudotime bin center (radians)")
    axes[1].set_ylabel("wMel titer")
    axes[1].set_title(f"{dataset_label}: median titer ± IQR across pseudotime")
    axes[1].set_xticks(x_)
    axes[1].set_xticklabels([f"{v:.2f}" for v in x_], rotation=45)
    axes[1].legend()

    plt.suptitle(f"wMel titer across cell-cycle pseudotime bins  [{dataset_label}]",
                 fontweight="bold")
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f"titer_vs_pseudotime_{tag}_pseudotime_bins.pdf"),
                bbox_inches="tight", dpi=150)
    plt.close()

    print(f"  Pseudotime analysis outputs written to {fig_dir}/ (tag: {tag})")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def _sanitize_obs(adata):
    """Cast object/category obs columns to str and replace 'nan' with 'NA'."""
    for col in adata.obs.columns:
        if adata.obs[col].dtype == object or str(adata.obs[col].dtype) == "category":
            adata.obs[col] = adata.obs[col].astype(str).replace("nan", "NA")
    return adata


def integrate(
    ref,
    query,
    out_path,
    fig_dir,
    sample,
    batch_key,
    min_cells,
    min_genes,
    n_pcs=30,
    leiden_resolution=0.5,
    ref_condition=None,
):
    """Run staged integration: subset query cells, ingest onto reference, save outputs.

    Parameters
    ----------
    ref : str
        Reference h5ad (pre-built: Harmony-corrected, scaled, PCA/UMAP/leiden).
    query : str
        Query h5ad (to be ingested onto reference).
    out_path : str
        Base path for output h5ad files (_reference, _query, _combined suffixes added).
    fig_dir : str
        Directory for output figures.
    sample : str
        Sample label for figure file names.
    batch_key : str
        Obs column used for batch correction (passed through, not re-run here).
    min_cells : int
        Minimum cells per gene (informational; filtering handled upstream).
    min_genes : int
        Minimum genes per cell for query QC filtering.
    n_pcs : int
        Number of PCs (informational; PCA already computed in reference).
    leiden_resolution : float
        Resolution used when building the reference (used in summary output only).
    ref_condition : str or list[str] or None
        cell_line value(s) that define the reference (default: ["JW18DOX"]).
        Used only to define the query mask on files[0]; files[1] IS the reference.
    """
    os.makedirs(fig_dir, exist_ok=True)
    sc.settings.figdir = fig_dir

    # files[0]: full integrated object — used only to identify query cells
    # files[1]: pre-built reference — passed directly to sc.tl.ingest
    adata_full = sc.read_h5ad(query)
    ref        = sc.read_h5ad(ref)

    print(f"Loaded full object : {adata_full.n_obs} cells")
    print(f"Loaded reference   : {ref.n_obs} cells, "
          f"{ref.obs['leiden'].nunique()} clusters")

    # Resolve reference conditions for query masking
    if ref_condition is None:
        ref_conditions = ["JW18DOX"]
        print(f"No --ref_condition specified. Defaulting to: {ref_conditions}")
    else:
        ref_conditions = ref_condition if isinstance(ref_condition, list) else [ref_condition]

    # Query = cells NOT in ref_conditions OR treatment != Ctrl
    # (i.e. the intermediate SV timepoints)
    query_mask = (
        ~adata_full.obs["cell_line"].isin(ref_conditions) |
        (adata_full.obs["treatment"] != "Ctrl")
    )

    print(f"Reference conditions : {ref_conditions} (treatment=Ctrl only)")
    print(f"Query (SV timepoints):")
    print(adata_full.obs.loc[query_mask, "bio_condition"].value_counts().to_string())
    print(f"Query cells          : {query_mask.sum()}")

    if query_mask.sum() == 0:
        raise ValueError("No query cells found (expected treatment=SV timepoint rows).")

    adata_query = adata_full[query_mask].copy()

    query, combined = map_query_to_reference(
        adata_query, ref=ref,
        min_genes=min_genes,
        fig_dir=fig_dir, sample=sample,
    )

    analyze_titer_vs_clusters(query, fig_dir=fig_dir, sample=sample)

    analyze_titer_vs_pseudotime(ref,   fig_dir=fig_dir, sample=sample,
                                 dataset_label="reference",
                                 cluster_col="leiden")
    analyze_titer_vs_pseudotime(query, fig_dir=fig_dir, sample=sample,
                                 dataset_label="query",
                                 cluster_col="leiden_ref")

    ref      = _sanitize_obs(ref)
    query    = _sanitize_obs(query)
    combined = _sanitize_obs(combined)

    ref_out      = out_path.replace(".h5ad", "_reference.h5ad")
    query_out    = out_path.replace(".h5ad", "_query.h5ad")
    combined_out = out_path.replace(".h5ad", "_combined.h5ad")

    ref.write(ref_out)
    query.write(query_out)
    combined.write(combined_out)

    print("\n" + "=" * 60)
    print("COMPLETE")
    print("=" * 60)
    print(f"Reference  -> {ref_out}")
    print(f"Query      -> {query_out}")
    print(f"Combined   -> {combined_out}")
    print(f"Figures    -> {fig_dir}/")
    print(f"Reference: {ref.n_obs} cells, {ref.obs['leiden'].nunique()} clusters "
          f"at resolution={leiden_resolution}")
    print(f"Query:     {query.n_obs} cells across "
          f"{query.obs['timepoint'].nunique()} timepoints")


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Staged scRNA-seq integration: ingest query cells onto pre-built reference"
    )
    parser.add_argument("--ref",           required=True, help="Pre-built reference h5ad") 
    parser.add_argument("--query",         required=True, help="Query h5ad (to be ingested onto reference)")
    parser.add_argument("--sample",        default="wolbachia_infection")
    parser.add_argument("--batch_key",     default="batch")
    parser.add_argument("--min_cells",     type=int,   default=3)
    parser.add_argument("--min_genes",     type=int,   default=200)
    parser.add_argument("--out_path",      default="all_integrated_by_cellcycle.h5ad")
    parser.add_argument("--fig_dir",       default="figures")
    parser.add_argument("--n_pcs",         type=int,   default=30)
    parser.add_argument("--resolution",    type=float, default=0.5)
    parser.add_argument("--ref_condition", type=str,   nargs="+", default=None,
                        help="cell_line value(s) used to define query mask on files[0] "
                             "(default: JW18DOX). files[1] is always used as the reference.")

    args = parser.parse_args()

    integrate(
        ref=args.ref,
        query=args.query,
        out_path=args.out_path,
        fig_dir=args.fig_dir,
        sample=args.sample,
        batch_key=args.batch_key,
        min_cells=args.min_cells,
        min_genes=args.min_genes,
        n_pcs=args.n_pcs,
        leiden_resolution=args.resolution,
        ref_condition=args.ref_condition,
    )


if __name__ == "__main__":
    main()
