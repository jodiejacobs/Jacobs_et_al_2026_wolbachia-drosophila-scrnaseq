"""
integrate.py
============
Joint clustering pipeline for studying Wolbachia infection dynamics.

Strategy
--------
1. Load all samples, extract metadata from filenames.
2. Preprocess: filter, HVG, PCA.
3. Correct for library-prep method ONLY (BBKNN) — biological signal preserved.
4. Cluster ALL cells together (DOX-Ctrl, wMel-Ctrl, D7, D28, D56).
5. Ask which clusters are enriched at which timepoints/conditions.
6. Ask how Wolbachia titer varies across clusters.
7. Export SCEPTIC-ready files using cluster + condition labels.

This approach is preferred when:
  - You don't expect entirely new cell states to appear upon infection
  - You want to know which existing cell states are most affected by infection
  - You want titer as a continuous variable overlaid on the cluster landscape

Biological interpretation:
  - JW18DOX-Ctrl  = uninfected baseline
  - D7/D28/D56    = infection intermediates
  - JW18wMel-Ctrl = stably infected endpoint
  - Clusters enriched in wMel / high titer = infection-responsive cell states

Usage
-----
python integrate.py \\
    --files results/filtered_h5ad/*.h5ad \\
    --sample wolbachia_infection \\
    --out_path results/integrated/integrated.h5ad \\
    --fig_dir results/integrated/figures \\
    --optimize_resolution \\
    --resolutions 0.1 0.2 0.3 0.4 0.5 0.6 0.8 1.0 1.2 1.5
"""

import os
import re
import argparse

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse
import seaborn as sns
import anndata as ad
import scanpy as sc
import bbknn
from scipy.stats import chi2_contingency, kruskal, mannwhitneyu
from sklearn.metrics import silhouette_score, adjusted_rand_score


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _get_cmap(name, n=None):
    """Matplotlib ≥3.7 compatible colormap helper."""
    cmap = matplotlib.colormaps[name]
    return cmap.resampled(n) if n is not None else cmap


def _leiden_colors(adata, key="leiden"):
    clusters = sorted(adata.obs[key].unique())
    cmap = _get_cmap("tab20")
    return [cmap(i % 20) for i in range(len(clusters))]


def _savefig(fig, path):
    fig.savefig(path, bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  Saved: {path}")


def _sanitize_obs(adata):
    """Convert all object/category obs columns to str for h5py compatibility."""
    for col in adata.obs.columns:
        if adata.obs[col].dtype == object or str(adata.obs[col].dtype) == "category":
            adata.obs[col] = adata.obs[col].astype(str).replace("nan", "NA")
    return adata


# ─────────────────────────────────────────────────────────────────────────────
# Metadata extraction
# ─────────────────────────────────────────────────────────────────────────────

def _extract_timepoint_numeric(row):
    """
    Numeric position on the infection pseudotime axis:
        JW18DOX-Ctrl  -> 0    (uninfected baseline)
        D7/D28/D56    -> 7, 28, 56  (infection intermediates)
        JW18wMel-Ctrl -> 99  (stable infected endpoint)
    """
    tp = row.get("timepoint", None)
    if tp is not None and pd.notna(tp) and str(tp) not in ("nan", "None", ""):
        m = re.search(r"(\d+)", str(tp))
        return int(m.group(1)) if m else 0
    elif str(row.get("cell_line", "")) == "JW18wMel":
        return 99
    else:
        return 0


def add_metadata(adata, batch_key):
    """Extract sample metadata from batch/filename strings."""
    adata.obs["cell_line"]  = adata.obs[batch_key].str.extract(r"(JW18DOX|JW18wMel)")[0]
    adata.obs["treatment"]  = adata.obs[batch_key].str.extract(r"-(Ctrl|SV)")[0]
    adata.obs["timepoint"]  = adata.obs[batch_key].str.extract(r"-(D\d+)-")[0]
    adata.obs["replicate"]  = adata.obs[batch_key].str.extract(r"-(\d+)_")[0]
    adata.obs["method"]     = adata.obs[batch_key].str.extract(r"_(10x|pipseq)$")[0]
    adata.obs["bio_condition"] = adata.obs.apply(
        lambda row: (
            f"{row['cell_line']}-{row['treatment']}-{row['timepoint']}"
            if pd.notna(row["timepoint"])
            else f"{row['cell_line']}-{row['treatment']}"
        ),
        axis=1,
    )
    adata.obs["timepoint_numeric"] = adata.obs.apply(
        _extract_timepoint_numeric, axis=1).astype(int)
    return adata


# ─────────────────────────────────────────────────────────────────────────────
# Preprocessing
# ─────────────────────────────────────────────────────────────────────────────

def preprocess(adata, min_genes, min_cells, n_pcs, n_top_genes=2000):
    """
    Standard scanpy preprocessing pipeline.

    Steps (order matters):
      1.  Remove bacterial genes
      2.  eliminate_zeros() — remove explicitly stored zeros from sparse matrix
      3.  Filter low-quality cells (min_genes, min_counts) and genes (min_cells)
      4.  Convert to dense float64, nan_to_num
      5.  Final zero-count cell/gene removal after dense conversion
      6.  HVG selection on raw counts (flavor='seurat', batch_key='method')
          NOTE: seurat flavor operates on raw counts intentionally here
      7.  normalize_total + log1p
      8.  nan_to_num after log1p
      9.  Store adata.raw (normalized log counts, pre-scale)
      10. Subset to HVGs
      11. scale (zero mean, unit variance)
      12. PCA

    Without normalize + log1p + scale, PCA is driven by library size
    rather than biological variation, producing 1-2 clusters.
    """
    # Remove bacterial / rRNA genes
    bacteria_genes = ["GQX67_00940", "GQX67_05945"] + [
        g for g in adata.var_names if g.startswith("16S_")]
    adata = adata[:, ~adata.var_names.isin(bacteria_genes)].copy()

    # Remove explicitly stored zeros from sparse matrix — these pass
    # min_counts filtering (the entry exists but sums to 0) and cause
    # normalize_total to warn "Some cells have zero counts"
    if scipy.sparse.issparse(adata.X):
        adata.X.eliminate_zeros()

    # Filter cells and genes
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_cells(adata, min_counts=1)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    sc.pp.filter_genes(adata, min_counts=1)

    # Convert to dense — required for nan_to_num and scale
    if scipy.sparse.issparse(adata.X):
        adata.X = adata.X.toarray()
    adata.X = np.nan_to_num(adata.X.astype(np.float64), nan=0.0, posinf=0.0, neginf=0.0)

    # Remove any zero-count cells/genes that appeared after dense conversion
    # (NaN values become 0.0 via nan_to_num, potentially creating empty rows)
    cell_sums = adata.X.sum(axis=1)
    gene_sums = adata.X.sum(axis=0)
    adata = adata[cell_sums > 0].copy()
    adata = adata[:, gene_sums > 0].copy()
    print(f"  After filtering: {adata.n_obs} cells, {adata.n_vars} genes")

    # Remove near-constant genes — these produce NaN dispersion (log(0)) in HVG
    gene_var = adata.X.var(axis=0)
    gene_mean = adata.X.mean(axis=0)
    keep = gene_var > np.maximum(gene_mean * 1e-6, 1e-10)
    adata = adata[:, keep].copy()
    print(f"  After removing near-constant genes: {adata.n_vars} genes")

    # HVG selection on raw counts (before normalization).
    # flavor="seurat" uses mean/dispersion on raw counts — correct here.
    # batch_key="method" selects HVGs that are variable in BOTH 10X and PIPseq,
    # preventing method-specific genes from driving clustering.
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=n_top_genes,
                                 batch_key="method")

    # Normalize + log1p
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    if scipy.sparse.issparse(adata.X):
        adata.X = adata.X.toarray()
    adata.X = np.nan_to_num(adata.X, nan=0.0, posinf=0.0, neginf=0.0)

    # Store normalized log counts before scaling (needed for DE, visualization)
    adata.raw = adata

    # Subset to HVGs, scale, PCA
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)
    sc.pp.pca(adata, n_comps=n_pcs)

    return adata


# ─────────────────────────────────────────────────────────────────────────────
# Resolution optimisation
# ─────────────────────────────────────────────────────────────────────────────

def optimize_leiden_resolution(
    adata,
    resolutions=None,
    n_pcs=30,
    silhouette_subsample=5000,
    ari_n_runs=3,
    fig_dir="figures",
    sample="sample",
    random_state=42,
):
    if resolutions is None:
        resolutions = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0]

    os.makedirs(fig_dir, exist_ok=True)
    X_pca = adata.obsm.get("X_pca_harmony", adata.obsm["X_pca"])[:, :n_pcs]
    rng = np.random.default_rng(random_state)
    sil_idx = rng.choice(adata.n_obs, min(silhouette_subsample, adata.n_obs), replace=False)

    records = []
    cluster_assignments = {}

    print(f"\nSweeping {len(resolutions)} resolutions …")
    for res in resolutions:
        sc.tl.leiden(adata, resolution=res, random_state=random_state,
                     key_added=f"leiden_res{res}")
        labels = adata.obs[f"leiden_res{res}"].astype(int).values
        cluster_assignments[res] = labels
        n_clusters = len(np.unique(labels))
        sil = (silhouette_score(X_pca[sil_idx], labels[sil_idx], metric="euclidean")
               if n_clusters >= 2 else np.nan)

        ari_scores = []
        for seed in range(1, ari_n_runs + 1):
            sc.tl.leiden(adata, resolution=res, random_state=seed, key_added="_tmp")
            ari_scores.append(adjusted_rand_score(
                labels, adata.obs["_tmp"].astype(int).values))
        adata.obs.drop(columns=["_tmp"], inplace=True)

        records.append(dict(
            resolution=res, n_clusters=n_clusters,
            silhouette=sil, ari_stability=float(np.mean(ari_scores))
        ))
        print(f"  res={res:.2f}  n_clusters={n_clusters:3d}  "
              f"silhouette={sil:.4f}  ari={records[-1]['ari_stability']:.4f}")

    scores_df = pd.DataFrame(records)
    scores_df.to_csv(os.path.join(fig_dir, f"resolution_scores_{sample}.csv"), index=False)

    def _norm(v):
        lo, hi = np.nanmin(v), np.nanmax(v)
        return (v - lo) / (hi - lo + 1e-12)

    scores_df["composite"] = (
        _norm(scores_df["silhouette"].values.astype(float)) +
        _norm(scores_df["ari_stability"].values.astype(float))
    ) / 2
    best_res = float(scores_df.loc[scores_df["composite"].idxmax(), "resolution"])
    print(f"\nBest resolution: {best_res}")

    res_vals = scores_df["resolution"].values

    # Metric curves
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    for ax, col, color, title in zip(
        axes,
        ["silhouette", "ari_stability", "n_clusters"],
        ["#2196F3", "#4CAF50", "#FF9800"],
        ["Silhouette Score", "ARI Stability", "N Clusters"],
    ):
        ax.plot(res_vals, scores_df[col].values, "o-", color=color, linewidth=2, markersize=6)
        ax.axvline(best_res, color="red", linestyle="--", linewidth=1.5,
                   label=f"Best: {best_res}")
        ax.set_xlabel("Leiden Resolution"); ax.set_title(title)
        ax.legend(fontsize=9); ax.grid(alpha=0.3)
    plt.suptitle(f"Resolution Optimisation — {sample}", fontweight="bold", y=1.02)
    plt.tight_layout()
    _savefig(fig, os.path.join(fig_dir, f"resolution_metrics_{sample}.pdf"))

    # UMAP panel
    ncols = min(4, len(resolutions))
    nrows = int(np.ceil(len(resolutions) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4, nrows * 3.5))
    axes = np.array(axes).flatten()
    umap_xy = adata.obsm["X_umap"]
    for i, res in enumerate(resolutions):
        ax = axes[i]
        lbls = adata.obs[f"leiden_res{res}"].astype(int).values
        n_cl = len(np.unique(lbls))
        ax.scatter(umap_xy[:, 0], umap_xy[:, 1], c=lbls,
                   cmap=_get_cmap("tab20", n_cl), s=1, alpha=0.5, rasterized=True)
        is_best = np.isclose(res, best_res)
        ax.set_title(f"res={res} (n={n_cl})" + (" *" if is_best else ""),
                     fontsize=9, fontweight="bold" if is_best else "normal",
                     color="red" if is_best else "black")
        ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)
    plt.suptitle(f"UMAP at each resolution — {sample}", fontweight="bold")
    plt.tight_layout()
    _savefig(fig, os.path.join(fig_dir, f"resolution_umap_panel_{sample}.pdf"))

    # Composite score bar
    fig, ax = plt.subplots(figsize=(10, 4))
    bar_colors = ["red" if np.isclose(r, best_res) else "#90CAF9" for r in res_vals]
    ax.bar(res_vals.astype(str), scores_df["composite"].values,
           color=bar_colors, edgecolor="black", linewidth=0.5)
    ax.set_xlabel("Leiden Resolution")
    ax.set_ylabel("Composite Score")
    ax.set_title(f"Resolution Composite Score — {sample}\nBest: {best_res} (red)")
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    _savefig(fig, os.path.join(fig_dir, f"resolution_composite_{sample}.pdf"))

    return best_res, scores_df


# ─────────────────────────────────────────────────────────────────────────────
# Joint clustering
# ─────────────────────────────────────────────────────────────────────────────

def cluster_all(
    adata,
    batch_key,
    min_genes,
    min_cells,
    n_pcs,
    fig_dir,
    sample,
    optimize_resolution=True,
    resolutions=None,
    leiden_resolution=0.5,
):
    """
    Cluster all cells together after correcting for method only.

    Steps:
      1. Preprocess (filter, HVG, PCA)
      2. BBKNN — correct for 10X vs PIPseq only
      3. UMAP
      4. Leiden clustering (optimised or fixed resolution)

    Biological variation (cell line, timepoint, titer) is NOT corrected —
    it drives the clustering and is visible in the embedding.
    """
    print("\n" + "=" * 60)
    print("JOINT CLUSTERING — all conditions, method correction only")
    print("=" * 60)

    print(f"\nTotal cells: {adata.n_obs}")
    print("Condition breakdown:")
    print(adata.obs["bio_condition"].value_counts().to_string())
    print("\nMethod breakdown:")
    print(adata.obs["method"].value_counts().to_string())

    adata = preprocess(adata, min_genes, min_cells, n_pcs)
    print(f"\nAfter preprocessing: {adata.n_obs} cells, {adata.n_vars} genes")

    # Pre-correction UMAP (QC)
    tmp = adata.copy()
    sc.pp.neighbors(tmp, n_pcs=n_pcs)
    sc.tl.umap(tmp)
    sc.pl.umap(tmp, color=["method", "bio_condition"],
               save=f"_{sample}_before_correction.pdf", ncols=2,
               title=["Method (pre-correction)", "Bio condition (pre-correction)"])
    del tmp

    # BBKNN: correct for method only
    print("\nRunning BBKNN (batch_key='method') …")
    # Harmony: corrects PCA embedding for both method and replicate simultaneously.
    import harmonypy
    ho = harmonypy.run_harmony(
        adata.obsm["X_pca"][:, :n_pcs],
        adata.obs,
        vars_use=["method", "replicate"],
        max_iter_harmony=20,
        random_state=42,
    )
    adata.obsm["X_pca_harmony"] = ho.Z_corr.T
    sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_pcs=n_pcs)
    sc.tl.umap(adata)

    # Post-correction QC
    sc.pl.umap(adata, color=["method", "bio_condition", "cell_line", "timepoint_numeric"],
               save=f"_{sample}_after_correction.pdf", ncols=2,
               title=["Method (should overlap post-correction)",
                      "Bio condition (should still separate)",
                      "Cell line",
                      "Timepoint (numeric)"])

    # Resolution optimisation
    if optimize_resolution:
        final_resolution, _ = optimize_leiden_resolution(
            adata, resolutions=resolutions, n_pcs=n_pcs,
            fig_dir=fig_dir, sample=sample)
    else:
        final_resolution = leiden_resolution

    print(f"\nClustering at resolution={final_resolution} …")
    sc.tl.leiden(adata, resolution=final_resolution, key_added="leiden")

    n_clusters = adata.obs["leiden"].nunique()
    print(f"Clusters: {n_clusters}")
    print(adata.obs["leiden"].value_counts().sort_index().to_string())

    # Core UMAPs
    sc.pl.umap(adata, color=["leiden", "bio_condition", "method", "cell_line"],
               save=f"_{sample}_clusters.pdf", ncols=2, legend_loc="on data",
               title=["Leiden clusters", "Bio condition", "Method", "Cell line"])

    if "wolbachia_titer" in adata.obs.columns:
        vmax = float(adata.obs["wolbachia_titer"].quantile(0.95))
        sc.pl.umap(adata, color="wolbachia_titer", vmax=vmax, cmap="viridis",
                   save=f"_{sample}_titer.pdf",
                   title="Wolbachia titer")

    return adata, final_resolution


# ─────────────────────────────────────────────────────────────────────────────
# Condition enrichment analysis
# ─────────────────────────────────────────────────────────────────────────────

def analyze_condition_enrichment(adata, fig_dir, sample):
    """
    Which clusters are enriched at which conditions/timepoints?

    Produces:
      - Stacked bar: cluster composition per condition (% cells)
      - Stacked bar: condition composition per cluster (% cells)
      - Heatmap: fold-enrichment vs expected (chi-square residuals)
      - UMAP panels: one per condition, highlighting those cells
      - Chi-square test for cluster x condition association
    """
    print("\n" + "=" * 60)
    print("CONDITION ENRICHMENT ANALYSIS")
    print("=" * 60)

    colors = _leiden_colors(adata, key="leiden")
    clusters = sorted(adata.obs["leiden"].unique())

    # Define ordered conditions for consistent plotting
    condition_order = []
    for expected in ["JW18DOX-Ctrl", "JW18wMel-D7", "JW18wMel-D28",
                     "JW18wMel-D56", "JW18wMel-Ctrl"]:
        if expected in adata.obs["bio_condition"].values:
            condition_order.append(expected)
    # Add any remaining conditions not in expected list
    for cond in sorted(adata.obs["bio_condition"].unique()):
        if cond not in condition_order:
            condition_order.append(cond)

    # ── Chi-square ────────────────────────────────────────────────────────────
    contingency = pd.crosstab(adata.obs["leiden"], adata.obs["bio_condition"])
    chi2, p, dof, expected = chi2_contingency(contingency)
    n = contingency.sum().sum()
    cramers_v = np.sqrt(chi2 / (n * (min(contingency.shape) - 1)))
    print(f"\nChi-square: chi2={chi2:.2f}  p={p:.2e}  dof={dof}")
    print(f"Cramer's V: {cramers_v:.3f}")

    # ── Plot 1: % cluster composition per condition (stacked bar) ─────────────
    pct_by_cond = pd.crosstab(adata.obs["leiden"],
                               adata.obs["bio_condition"],
                               normalize="columns") * 100
    # Reorder columns
    pct_by_cond = pct_by_cond.reindex(
        columns=[c for c in condition_order if c in pct_by_cond.columns])

    fig, ax = plt.subplots(figsize=(max(10, len(condition_order) * 1.2), 6))
    pct_by_cond.T.plot(kind="bar", stacked=True, ax=ax,
                       color=colors, width=0.8, edgecolor="black", linewidth=0.2)
    ax.set_xlabel("Condition")
    ax.set_ylabel("% of cells")
    ax.set_title(f"Cluster composition per condition — {sample}\n"
                 f"chi2={chi2:.1f}  p={p:.2e}  Cramer's V={cramers_v:.3f}")
    ax.legend(title="Cluster", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    _savefig(fig, os.path.join(fig_dir, f"cluster_composition_by_condition_{sample}.pdf"))

    # ── Plot 2: % condition per cluster (stacked bar) ─────────────────────────
    pct_by_cluster = pd.crosstab(adata.obs["bio_condition"],
                                  adata.obs["leiden"],
                                  normalize="columns") * 100
    cond_colors = plt.cm.Set2(np.linspace(0, 1, len(condition_order)))

    fig, ax = plt.subplots(figsize=(max(10, len(clusters) * 1.0), 6))
    pct_by_cluster.reindex(
        index=[c for c in condition_order if c in pct_by_cluster.index]
    ).T.plot(kind="bar", stacked=True, ax=ax,
             color=cond_colors, width=0.8, edgecolor="black", linewidth=0.2)
    ax.set_xlabel("Leiden Cluster")
    ax.set_ylabel("% of cells")
    ax.set_title(f"Condition composition per cluster — {sample}")
    ax.legend(title="Condition", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.xticks(rotation=0)
    plt.tight_layout()
    _savefig(fig, os.path.join(fig_dir, f"condition_composition_by_cluster_{sample}.pdf"))

    # ── Plot 3: fold-enrichment heatmap ───────────────────────────────────────
    # Observed / expected (from chi-square) — highlights enrichment above baseline
    obs = contingency.values.astype(float)
    fold_enrich = np.log2((obs + 1) / (expected + 1))  # log2 fold enrichment
    fold_df = pd.DataFrame(fold_enrich,
                           index=contingency.index,
                           columns=contingency.columns)
    fold_df = fold_df.reindex(
        columns=[c for c in condition_order if c in fold_df.columns])

    fig, ax = plt.subplots(figsize=(max(10, len(condition_order) * 1.2), max(6, len(clusters) * 0.5)))
    sns.heatmap(fold_df, cmap="RdBu_r", center=0, ax=ax,
                cbar_kws={"label": "log2 fold enrichment vs expected"},
                linewidths=0.3, annot=(fold_df.shape[0] * fold_df.shape[1] <= 80),
                fmt=".1f")
    ax.set_xlabel("Condition")
    ax.set_ylabel("Leiden Cluster")
    ax.set_title(f"Cluster enrichment per condition (log2 fold) — {sample}\n"
                 "Positive = more cells than expected; negative = fewer")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    _savefig(fig, os.path.join(fig_dir, f"cluster_enrichment_heatmap_{sample}.pdf"))

    # ── Plot 4: UMAP panels per condition ─────────────────────────────────────
    conditions = [c for c in condition_order
                  if c in adata.obs["bio_condition"].values]
    ncols = min(3, len(conditions))
    nrows = int(np.ceil(len(conditions) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4.5, nrows * 4))
    axes = np.array(axes).flatten()
    umap_xy = adata.obsm["X_umap"]
    all_labels = adata.obs["leiden"].astype(int).values
    n_cl = len(np.unique(all_labels))

    for i, cond in enumerate(conditions):
        ax = axes[i]
        mask = adata.obs["bio_condition"].values == cond
        ax.scatter(umap_xy[~mask, 0], umap_xy[~mask, 1],
                   c="lightgrey", s=1, alpha=0.2, rasterized=True)
        ax.scatter(umap_xy[mask, 0], umap_xy[mask, 1],
                   c=all_labels[mask], cmap=_get_cmap("tab20", n_cl),
                   s=2, alpha=0.8, rasterized=True)
        ax.set_title(f"{cond}\n(n={mask.sum():,})", fontsize=9)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)
    plt.suptitle(f"Cells per condition on UMAP — {sample}", fontweight="bold")
    plt.tight_layout()
    _savefig(fig, os.path.join(fig_dir, f"umap_by_condition_{sample}.pdf"))

    # Save tables
    pct_by_cond.to_csv(os.path.join(fig_dir, f"cluster_pct_by_condition_{sample}.csv"))
    fold_df.to_csv(os.path.join(fig_dir, f"cluster_enrichment_by_condition_{sample}.csv"))
    pd.DataFrame({
        "chi2": [chi2], "p_value": [p], "dof": [dof], "cramers_v": [cramers_v]
    }).to_csv(os.path.join(fig_dir, f"enrichment_stats_{sample}.csv"), index=False)

    print(f"\n  Saved enrichment tables and plots to {fig_dir}/")
    return pct_by_cond, fold_df


# ─────────────────────────────────────────────────────────────────────────────
# Titer analysis
# ─────────────────────────────────────────────────────────────────────────────

def analyze_titer_by_cluster(adata, fig_dir, sample):
    """
    How does Wolbachia titer vary across clusters?

    Produces:
      - Boxplot + strip: titer per cluster
      - Violin: titer per cluster
      - UMAP: titer as continuous variable
      - UMAP: titer split by method (QC)
      - Bar: % infected cells (titer > 0) per cluster
      - Heatmap: mean titer per cluster × condition
      - Kruskal-Wallis test
    """
    if "wolbachia_titer" not in adata.obs.columns:
        print("  Skipping titer analysis: no 'wolbachia_titer' column")
        return

    print("\n" + "=" * 60)
    print("WOLBACHIA TITER BY CLUSTER")
    print("=" * 60)

    colors = _leiden_colors(adata, key="leiden")
    clusters = sorted(adata.obs["leiden"].unique())
    obs = adata.obs[["leiden", "wolbachia_titer", "bio_condition",
                     "method", "timepoint_numeric"]].copy()
    obs_titer = obs.dropna(subset=["wolbachia_titer"])

    # ── Kruskal-Wallis ────────────────────────────────────────────────────────
    groups = [obs_titer[obs_titer["leiden"] == c]["wolbachia_titer"].values
              for c in clusters if (obs_titer["leiden"] == c).sum() > 0]
    h, p_kw = kruskal(*groups)
    print(f"\nKruskal-Wallis: H={h:.2f}  p={p_kw:.2e}")
    print("Titer {'SIGNIFICANTLY' if p_kw < 0.05 else 'does NOT'} differ across clusters")

    # Summary stats
    stats = obs_titer.groupby("leiden")["wolbachia_titer"].agg(
        ["mean", "median", "std", "count"])
    print("\nTiter summary per cluster:")
    print(stats.to_string())
    stats.to_csv(os.path.join(fig_dir, f"titer_stats_by_cluster_{sample}.csv"))

    # ── Plot 1: boxplot + strip ───────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(max(10, len(clusters) * 0.8), 6))
    sns.stripplot(data=obs_titer, x="leiden", y="wolbachia_titer",
                  order=clusters, color="black", alpha=0.15, size=1.5, ax=ax)
    bp = ax.boxplot(
        [obs_titer[obs_titer["leiden"] == c]["wolbachia_titer"].values
         for c in clusters],
        positions=range(len(clusters)), widths=0.55, patch_artist=True,
        whiskerprops=dict(alpha=0.7), capprops=dict(alpha=0.7),
        medianprops=dict(color="black", linewidth=2),
        flierprops=dict(markersize=1),
    )
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color); patch.set_alpha(0.75)
    ax.set_xticklabels(clusters)
    ax.set_xlabel("Leiden Cluster")
    ax.set_ylabel("Wolbachia Titer")
    ax.set_title(f"Wolbachia titer by cluster — {sample}\n"
                 f"Kruskal-Wallis H={h:.2f}  p={p_kw:.2e}")
    plt.tight_layout()
    _savefig(fig, os.path.join(fig_dir, f"titer_boxplot_by_cluster_{sample}.pdf"))

    # ── Plot 2: violin ────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(max(10, len(clusters) * 0.8), 6))
    sc.pl.violin(adata, "wolbachia_titer", groupby="leiden",
                 ax=ax, show=False, rotation=0)
    ax.set_title(f"Wolbachia titer distribution by cluster — {sample}\n"
                 f"Kruskal-Wallis H={h:.2f}  p={p_kw:.2e}")
    plt.tight_layout()
    _savefig(fig, os.path.join(fig_dir, f"titer_violin_by_cluster_{sample}.pdf"))

    # ── Plot 3: % infected per cluster ────────────────────────────────────────
    pct_infected = obs_titer.groupby("leiden").apply(
        lambda x: (x["wolbachia_titer"] > 0).sum() / len(x) * 100
    )
    fig, ax = plt.subplots(figsize=(max(8, len(clusters) * 0.8), 5))
    ax.bar(range(len(pct_infected)), pct_infected.values,
           color=colors, alpha=0.85, edgecolor="black", linewidth=0.4)
    ax.set_xticks(range(len(pct_infected)))
    ax.set_xticklabels(pct_infected.index)
    ax.set_xlabel("Leiden Cluster")
    ax.set_ylabel("% Wolbachia-positive cells")
    ax.set_title(f"Infection rate by cluster — {sample}")
    ax.grid(axis="y", alpha=0.3)
    for i, v in enumerate(pct_infected.values):
        ax.text(i, v + 0.5, f"{v:.1f}%", ha="center", va="bottom", fontsize=7)
    plt.tight_layout()
    _savefig(fig, os.path.join(fig_dir, f"infection_pct_by_cluster_{sample}.pdf"))

    # ── Plot 4: mean titer heatmap — cluster × condition ──────────────────────
    mean_titer = obs_titer.pivot_table(
        values="wolbachia_titer", index="leiden",
        columns="bio_condition", aggfunc="mean")

    fig, ax = plt.subplots(figsize=(max(10, mean_titer.shape[1] * 1.2),
                                    max(5, mean_titer.shape[0] * 0.5)))
    sns.heatmap(mean_titer, cmap="viridis", ax=ax,
                cbar_kws={"label": "Mean Wolbachia Titer"},
                linewidths=0.3)
    ax.set_xlabel("Condition")
    ax.set_ylabel("Leiden Cluster")
    ax.set_title(f"Mean Wolbachia titer — cluster × condition — {sample}")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    _savefig(fig, os.path.join(fig_dir, f"titer_heatmap_cluster_condition_{sample}.pdf"))

    # ── Plot 5: titer by cluster × timepoint (infected cells only) ────────────
    infected = obs_titer[obs_titer["wolbachia_titer"] > 0]
    if len(infected) > 0:
        fig, ax = plt.subplots(figsize=(12, 5))
        sns.boxplot(data=infected, x="leiden", y="wolbachia_titer",
                    hue="timepoint_numeric", ax=ax,
                    flierprops=dict(markersize=1), order=clusters)
        ax.set_xlabel("Leiden Cluster")
        ax.set_ylabel("Wolbachia Titer (infected cells only)")
        ax.set_title(f"Titer by cluster and timepoint — {sample}")
        ax.legend(title="Timepoint", bbox_to_anchor=(1.05, 1), loc="upper left")
        plt.tight_layout()
        _savefig(fig, os.path.join(fig_dir,
                                   f"titer_by_cluster_timepoint_{sample}.pdf"))

    return stats


# ─────────────────────────────────────────────────────────────────────────────
# SCEPTIC export
# ─────────────────────────────────────────────────────────────────────────────

def export_sceptic(adata, fig_dir, sample, n_pcs=30):
    """
    Export SCEPTIC inputs from the jointly-clustered object.

    Uses wMel cells only (DOX→D7→D28→D56→wMel-Ctrl axis).
    DOX cells provide the t=0 anchor as the uninfected reference state.
    """
    print("\n" + "=" * 60)
    print("EXPORTING SCEPTIC INPUTS")
    print("=" * 60)

    os.makedirs(fig_dir, exist_ok=True)

    # Build the full pseudotime axis:
    # DOX-Ctrl (t=0) + wMel intermediates (D7/D28/D56) + wMel-Ctrl (t=999)
    sceptic_mask = (
        (adata.obs["cell_line"] == "JW18DOX") & (adata.obs["treatment"] == "Ctrl")
    ) | (
        adata.obs["cell_line"] == "JW18wMel"
    )
    adata_s = adata[sceptic_mask].copy()

    print(f"SCEPTIC cells: {adata_s.n_obs}")
    print("Timepoint breakdown:")
    print(adata_s.obs["timepoint_numeric"].value_counts().sort_index().to_string())

    # Feature matrix: PCA coords
    n_actual = min(n_pcs, adata_s.obsm["X_pca"].shape[1])
    if n_actual < n_pcs:
        print(f"  WARNING: Using {n_actual} PCs (requested {n_pcs})")

    data_mat     = adata_s.obsm["X_pca"][:, :n_actual]
    labels       = adata_s.obs["timepoint_numeric"].values.astype(int)
    label_list   = np.sort(np.unique(labels))
    cell_ids     = adata_s.obs_names.tolist()

    mat_path  = os.path.join(fig_dir, f"sceptic_matrix_{sample}.csv")
    lab_path  = os.path.join(fig_dir, f"sceptic_labels_{sample}.csv")
    ll_path   = os.path.join(fig_dir, f"sceptic_label_list_{sample}.csv")
    meta_path = os.path.join(fig_dir, f"sceptic_metadata_{sample}.csv")

    pd.DataFrame(data_mat, index=cell_ids,
                 columns=[f"PC{i+1}" for i in range(n_actual)]).to_csv(mat_path)
    pd.Series(labels, index=cell_ids, name="timepoint").to_csv(lab_path)
    pd.Series(label_list, name="timepoint").to_csv(ll_path, index=False)
    pd.DataFrame({
        "cell_id":   cell_ids,
        "timepoint": labels,
        "leiden":    adata_s.obs["leiden"].values,
        "method":    adata_s.obs["method"].values,
        "bio_condition": adata_s.obs["bio_condition"].values,
    }).to_csv(meta_path, index=False)

    print(f"  Matrix     -> {mat_path}  shape={data_mat.shape}")
    print(f"  Labels     -> {lab_path}")
    print(f"  Label list -> {ll_path}  values={label_list}")
    print(f"  Metadata   -> {meta_path}")


# ─────────────────────────────────────────────────────────────────────────────
# Main integration function
# ─────────────────────────────────────────────────────────────────────────────

def integrate(
    files,
    out_path,
    fig_dir,
    sample,
    batch_key,
    min_cells,
    min_genes,
    n_pcs=30,
    optimize_resolution=True,
    resolutions=None,
    leiden_resolution=0.5,
):
    os.makedirs(fig_dir, exist_ok=True)
    sc.settings.figdir = fig_dir

    # ── Load and concatenate ──────────────────────────────────────────────────
    print("Loading files …")
    adatas = []
    for fp in files:
        a = sc.read_h5ad(fp)
        a.obs[batch_key] = os.path.splitext(os.path.basename(fp))[0]
        adatas.append(a)

    adata = ad.concat(adatas, join="inner", merge="same", index_unique="-")
    adata = add_metadata(adata, batch_key)

    print(f"\nLoaded {adata.n_obs:,} cells from {len(files)} files")
    print("Condition breakdown:")
    print(adata.obs["bio_condition"].value_counts().to_string())

    # ── Joint clustering ──────────────────────────────────────────────────────
    adata, final_resolution = cluster_all(
        adata,
        batch_key=batch_key,
        min_genes=min_genes,
        min_cells=min_cells,
        n_pcs=n_pcs,
        fig_dir=fig_dir,
        sample=sample,
        optimize_resolution=optimize_resolution,
        resolutions=resolutions,
        leiden_resolution=leiden_resolution,
    )

    # ── Condition enrichment ──────────────────────────────────────────────────
    analyze_condition_enrichment(adata, fig_dir, sample)

    # ── Titer analysis ────────────────────────────────────────────────────────
    analyze_titer_by_cluster(adata, fig_dir, sample)

    # ── SCEPTIC export ────────────────────────────────────────────────────────
    export_sceptic(adata, fig_dir=fig_dir, sample=sample, n_pcs=n_pcs)

    # ── Save ──────────────────────────────────────────────────────────────────
    adata = _sanitize_obs(adata)
    adata.write(out_path)
    print(f"\nSaved: {out_path}")

    print("\n" + "=" * 60)
    print("COMPLETE")
    print("=" * 60)
    print(f"Integrated object -> {out_path}")
    print(f"Figures           -> {fig_dir}/")
    print(f"Cells: {adata.n_obs:,}  |  Clusters: {adata.obs['leiden'].nunique()}  "
          f"at resolution={final_resolution}")


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Joint clustering pipeline: method correction only, "
            "then analyse condition enrichment and Wolbachia titer per cluster."
        )
    )
    parser.add_argument("--files",      required=True, nargs="+",
                        help="h5ad files to integrate")
    parser.add_argument("--sample",     default="wolbachia_infection",
                        help="Label used in all output filenames")
    parser.add_argument("--batch_key",  default="batch",
                        help="obs column for sample/batch identity")
    parser.add_argument("--min_cells",  type=int, default=3)
    parser.add_argument("--min_genes",  type=int, default=200)
    parser.add_argument("--out_path",   default="integrated.h5ad")
    parser.add_argument("--fig_dir",    default="figures")
    parser.add_argument("--n_pcs",      type=int, default=30)
    parser.add_argument("--resolution", type=float, default=0.5,
                        help="Fixed Leiden resolution (ignored if --optimize_resolution)")
    parser.add_argument("--optimize_resolution", action="store_true", default=True)
    parser.add_argument("--no_optimize_resolution", dest="optimize_resolution",
                        action="store_false")
    parser.add_argument("--resolutions", type=float, nargs="+", default=None,
                        help="Custom resolution sweep values")

    args = parser.parse_args()

    integrate(
        files=args.files,
        out_path=args.out_path,
        fig_dir=args.fig_dir,
        sample=args.sample,
        batch_key=args.batch_key,
        min_cells=args.min_cells,
        min_genes=args.min_genes,
        n_pcs=args.n_pcs,
        optimize_resolution=args.optimize_resolution,
        resolutions=args.resolutions,
        leiden_resolution=args.resolution,
    )


if __name__ == "__main__":
    main()
