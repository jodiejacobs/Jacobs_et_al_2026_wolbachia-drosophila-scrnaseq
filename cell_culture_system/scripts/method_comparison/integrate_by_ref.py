"""
integrate_by_ref.py
============
Staged scRNA-seq integration pipeline for studying Wolbachia infection dynamics.

Strategy
--------
1. Correct for library-prep method (10X vs PIPseq) with BBKNN — preserves biology.
2. Build a reference embedding from uninfected control cells only.
3. Map (ingest) new-infection timepoint cells onto that reference.
4. Cluster on the reference; project cluster labels to query cells.
5. Export SCEPTIC-ready files for pseudotime inference on infected cells.

SCEPTIC inputs produced
-----------------------
- sceptic_matrix_{sample}.csv    : cells × PCA dims (or HVGs)
- sceptic_labels_{sample}.csv    : numeric timepoint label per cell (0 = uninfected)
- sceptic_label_list_{sample}.csv: unique ordered timepoints
- sceptic_metadata_{sample}.csv  : per-cell cluster, method, timepoint
"""

import os
import argparse
import warnings
import re

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.sparse
import anndata as ad
import scanpy as sc
import bbknn
from scipy.stats import chi2_contingency, mannwhitneyu
from sklearn.metrics import silhouette_score, adjusted_rand_score


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _leiden_colors(adata, key="leiden"):
    clusters = sorted(adata.obs[key].unique())
    cmap = matplotlib.colormaps["tab20"]
    return [cmap(i % 20) for i in range(len(clusters))]


def _extract_timepoint_numeric(row):
    """
    Assign a numeric position on the infection pseudotime axis.

        JW18DOX-Ctrl   = uninfected baseline        -> t=0
        D7 / D28 / D56 = infection intermediates    -> t=7, 28, 56
        JW18wMel-Ctrl  = stably infected endpoint   -> t=999
    """
    tp = row.get("timepoint", None)
    if tp is not None and pd.notna(tp) and str(tp) not in ("nan", "None", ""):
        m = re.search(r"(\d+)", str(tp))
        return int(m.group(1)) if m else 0
    elif str(row.get("cell_line", "")) == "JW18wMel":
        return 999
    else:
        return 0


# ─────────────────────────────────────────────────────────────────────────────
# Preprocessing
# ─────────────────────────────────────────────────────────────────────────────

def add_metadata(adata, batch_key):
    """Extract sample metadata from batch names."""
    adata.obs["cell_line"]  = adata.obs[batch_key].str.extract(r"(JW18DOX|JW18wMel)")[0]
    adata.obs["treatment"]  = adata.obs[batch_key].str.extract(r"-(Ctrl|SV)")[0]
    adata.obs["timepoint"]  = adata.obs[batch_key].str.extract(r"-(D\d+)-")[0]
    adata.obs["replicate"]  = adata.obs[batch_key].str.extract(r"-(\d+)_")[0]
    adata.obs["method"]     = adata.obs[batch_key].str.extract(r"_(10x|pipseq)$")[0]
    adata.obs["bio_condition"] = adata.obs.apply(
        lambda row: (f"{row['cell_line']}-{row['treatment']}-{row['timepoint']}"
                     if pd.notna(row["timepoint"])
                     else f"{row['cell_line']}-{row['treatment']}"),
        axis=1,
    )
    adata.obs["timepoint_numeric"] = adata.obs.apply(
        _extract_timepoint_numeric, axis=1).astype(int)
    return adata


def preprocess(adata, min_genes, min_cells, n_pcs, n_top_genes=2000):
    """
    Standard scanpy preprocessing pipeline.

    Steps (order matters):
      1.  Remove bacterial genes
      2.  eliminate_zeros() — remove explicitly stored zeros from sparse matrix
      3.  Filter cells (min_genes, min_counts) and genes (min_cells, min_counts)
      4.  Convert to dense float64, nan_to_num
      5.  Final zero-count cell/gene removal after dense conversion
      6.  HVG selection on raw counts (flavor='seurat', batch_key='method')
      7.  normalize_total + log1p + nan_to_num
      8.  Store adata.raw (normalized log counts, pre-scale)
      9.  Subset to HVGs, scale, PCA

    Without normalize + log1p + scale, PCA is driven by library size
    rather than biological variation, producing 1-2 clusters.
    """
    bacteria_genes = ["GQX67_00940", "GQX67_05945"] + [
        g for g in adata.var_names if g.startswith("16S_")]
    adata = adata[:, ~adata.var_names.isin(bacteria_genes)].copy()

    # Remove explicitly stored zeros — these pass min_counts filtering
    # but cause normalize_total to warn "Some cells have zero counts"
    if scipy.sparse.issparse(adata.X):
        adata.X.eliminate_zeros()

    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_cells(adata, min_counts=1)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    sc.pp.filter_genes(adata, min_counts=1)

    if scipy.sparse.issparse(adata.X):
        adata.X = adata.X.toarray()
    adata.X = np.nan_to_num(adata.X.astype(np.float64), nan=0.0, posinf=0.0, neginf=0.0)

    # Remove zero-count cells/genes created by nan_to_num filling NaN->0
    cell_sums = adata.X.sum(axis=1)
    gene_sums = adata.X.sum(axis=0)
    adata = adata[cell_sums > 0].copy()
    adata = adata[:, gene_sums > 0].copy()
    print(f"  After filtering: {adata.n_obs} cells, {adata.n_vars} genes")

    # Remove near-constant genes — genes with very low variance produce
    # zero or NaN dispersion in HVG selection (log(0) = NaN).
    # Use a small epsilon threshold rather than exactly zero to catch
    # genes that are effectively constant in float64 arithmetic.
    gene_var = adata.X.var(axis=0)
    gene_mean = adata.X.mean(axis=0)
    # Keep genes with variance > 1e-6 * mean (relative threshold) or variance > 1e-10 (absolute)
    keep = (gene_var > np.maximum(gene_mean * 1e-6, 1e-10))
    adata = adata[:, keep].copy()
    print(f"  After removing near-constant genes: {adata.n_vars} genes")

    # HVG on raw counts — flavor="seurat" uses mean/dispersion on raw counts.
    # batch_key="method" selects HVGs variable in BOTH 10X and PIPseq.
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=n_top_genes,
                                 batch_key="method")

    # Normalize + log1p
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    if scipy.sparse.issparse(adata.X):
        adata.X = adata.X.toarray()
    adata.X = np.nan_to_num(adata.X, nan=0.0, posinf=0.0, neginf=0.0)

    # Store normalized log counts before scaling (for DE, visualization)
    adata.raw = adata

    # Subset to HVGs, then scale
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)

    # 6. PCA on scaled HVGs
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
            ari_scores.append(adjusted_rand_score(labels, adata.obs["_tmp"].astype(int).values))
        adata.obs.drop(columns=["_tmp"], inplace=True)

        records.append(dict(resolution=res, n_clusters=n_clusters,
                            silhouette=sil, ari_stability=float(np.mean(ari_scores))))
        print(f"  res={res:.2f}  n_clusters={n_clusters:3d}  "
              f"silhouette={sil:.4f}  ari={records[-1]['ari_stability']:.4f}")

    scores_df = pd.DataFrame(records)
    scores_df.to_csv(os.path.join(fig_dir, f"resolution_scores_{sample}.csv"), index=False)

    def _norm(v):
        lo, hi = np.nanmin(v), np.nanmax(v)
        return (v - lo) / (hi - lo + 1e-12)

    scores_df["composite"] = (_norm(scores_df["silhouette"].values.astype(float)) +
                               _norm(scores_df["ari_stability"].values.astype(float))) / 2
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
        ax.axvline(best_res, color="red", linestyle="--", linewidth=1.5, label=f"Best: {best_res}")
        ax.set_xlabel("Leiden Resolution"); ax.set_title(title)
        ax.legend(fontsize=9); ax.grid(alpha=0.3)
    plt.suptitle(f"Resolution Optimisation — {sample}", fontweight="bold", y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f"resolution_metrics_{sample}.pdf"), bbox_inches="tight")
    plt.close()

    # UMAP panel
    ncols = min(4, len(resolutions))
    nrows = int(np.ceil(len(resolutions) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4, nrows * 3.5))
    axes = np.array(axes).flatten()
    for i, res in enumerate(resolutions):
        ax = axes[i]
        umap_xy = adata.obsm["X_umap"]
        lbls = adata.obs[f"leiden_res{res}"].astype(int).values
        n_cl = len(np.unique(lbls))
        ax.scatter(umap_xy[:, 0], umap_xy[:, 1], c=lbls,
                   cmap=matplotlib.colormaps["tab20"].resampled(n_cl),
                   s=1, alpha=0.5, rasterized=True)
        is_best = np.isclose(res, best_res)
        ax.set_title(f"res={res} (n={n_cl})" + (" *" if is_best else ""),
                     fontsize=9, fontweight="bold" if is_best else "normal",
                     color="red" if is_best else "black")
        ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)
    plt.suptitle(f"UMAP at each resolution — {sample}", fontweight="bold")
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f"resolution_umap_panel_{sample}.pdf"),
                bbox_inches="tight", dpi=150)
    plt.close()

    # Overlap heatmaps
    n_pairs = len(resolutions) - 1
    ncols_h = min(3, n_pairs)
    nrows_h = int(np.ceil(n_pairs / ncols_h))
    fig, axes = plt.subplots(nrows_h, ncols_h, figsize=(ncols_h * 4.5, nrows_h * 4))
    axes = np.array(axes).flatten()
    for i in range(n_pairs):
        res_lo, res_hi = resolutions[i], resolutions[i + 1]
        overlap = pd.crosstab(
            pd.Series(cluster_assignments[res_lo], name=f"res={res_lo}"),
            pd.Series(cluster_assignments[res_hi], name=f"res={res_hi}"),
            normalize="index",
        )
        sns.heatmap(overlap, ax=axes[i], cmap="YlOrRd", vmin=0, vmax=1,
                    linewidths=0.3, cbar_kws={"shrink": 0.7},
                    annot=(overlap.shape[0] * overlap.shape[1] <= 100), fmt=".2f")
        axes[i].set_title(f"{res_lo} -> {res_hi}", fontsize=9)
        axes[i].tick_params(labelsize=7)
    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)
    plt.suptitle(f"Cluster overlap — {sample}", fontweight="bold")
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f"resolution_overlap_heatmaps_{sample}.pdf"), bbox_inches="tight")
    plt.close()

    # Composite bar
    fig, ax = plt.subplots(figsize=(10, 4))
    bar_colors = ["red" if np.isclose(r, best_res) else "#90CAF9" for r in res_vals]
    ax.bar(res_vals.astype(str), scores_df["composite"].values,
           color=bar_colors, edgecolor="black", linewidth=0.5)
    ax.set_xlabel("Leiden Resolution")
    ax.set_ylabel("Composite Score (normalised silhouette + ARI)")
    ax.set_title(f"Resolution Composite Score — {sample}\nBest: {best_res} (red)")
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f"resolution_composite_{sample}.pdf"), bbox_inches="tight")
    plt.close()

    return best_res, scores_df


# ─────────────────────────────────────────────────────────────────────────────
# Stage 1 — Reference
# ─────────────────────────────────────────────────────────────────────────────

def build_reference(adata_ref, batch_key, min_genes, min_cells, n_pcs,
                    fig_dir, sample, optimize_resolution=True,
                    resolutions=None, leiden_resolution=0.5):
    """Build reference embedding from uninfected control cells (JW18DOX-Ctrl)."""
    print("\n" + "=" * 60)
    print("STAGE 1 — Building reference from uninfected controls")
    print("=" * 60)

    ref = adata_ref.copy()
    print(f"Reference cells: {ref.n_obs}")
    print(f"Methods: {ref.obs['method'].value_counts().to_dict()}")

    ref = preprocess(ref, min_genes, min_cells, n_pcs)

    # Pre-correction UMAP
    ref_tmp = ref.copy()
    sc.pp.neighbors(ref_tmp, n_pcs=n_pcs)
    sc.tl.umap(ref_tmp)
    sc.pl.umap(ref_tmp, color=["method", "bio_condition"],
               save=f"_{sample}_ref_before_correction.pdf", ncols=2,
               title=["Method (pre-correction)", "Bio condition (pre-correction)"])
    del ref_tmp

    # BBKNN: correct for library-prep method only.
    # batch_key="batch" (individual samples) is preferred over "method" (2 levels)
    # because correcting on only 2 batches forces every 10x cell to connect to
    # every pipseq cell, over-integrating and creating stringy cluster artifacts.
    # neighbors_within_batch=3 is conservative — default of 5 is too aggressive
    # when biological variation is subtle relative to batch effect.
    print("Correcting for library-prep method only (BBKNN key='batch') …")
    # Harmony: corrects PCA embedding for both method and replicate simultaneously.
    # Preferred over BBKNN when method separation dominates — Harmony regresses
    # covariates out of X_pca directly, producing a corrected embedding
    # (X_pca_harmony) that standard neighbor/UMAP/clustering then uses.
    import harmonypy
    ho = harmonypy.run_harmony(
        ref.obsm["X_pca"][:, :n_pcs],
        ref.obs,
        vars_use=["method", "replicate"],  # correct for both simultaneously
        max_iter_harmony=20,
        random_state=42,
    )
    ref.obsm["X_pca_harmony"] = ho.Z_corr.T
    sc.pp.neighbors(ref, use_rep="X_pca_harmony", n_pcs=n_pcs)
    sc.tl.umap(ref)

    sc.pl.umap(ref, color=["method", "bio_condition", "cell_line"],
               save=f"_{sample}_ref_after_method_correction.pdf", ncols=3,
               title=["Method (post-correction — should overlap)",
                      "Bio condition (should still separate)",
                      "Cell line"])

    if optimize_resolution:
        final_resolution, _ = optimize_leiden_resolution(
            ref, resolutions=resolutions, n_pcs=n_pcs,
            fig_dir=fig_dir, sample=f"{sample}_ref")
    else:
        final_resolution = leiden_resolution

    print(f"\nFinal clustering at resolution={final_resolution} …")
    sc.tl.leiden(ref, resolution=final_resolution, key_added="leiden")

    sc.pl.umap(ref, color=["leiden", "bio_condition", "method"],
               save=f"_{sample}_ref_final_clusters.pdf", ncols=3,
               legend_loc="on data",
               title=["Leiden clusters (reference)", "Bio condition", "Method"])

    print(f"Reference: {ref.n_obs} cells, {ref.obs['leiden'].nunique()} clusters")
    return ref, final_resolution


# ─────────────────────────────────────────────────────────────────────────────
# Stage 2 — Query: map onto reference via ingest
# ─────────────────────────────────────────────────────────────────────────────

def map_query_to_reference(adata_query, ref, batch_key, min_genes, n_pcs,
                            fig_dir, sample):
    """Project new-infection timepoint cells onto the reference using sc.tl.ingest."""
    print("\n" + "=" * 60)
    print("STAGE 2 — Mapping query onto reference")
    print("=" * 60)

    query = adata_query.copy()
    print(f"Query cells: {query.n_obs}")
    print(f"Timepoints:  {query.obs['timepoint'].value_counts().to_dict()}")
    print(f"Methods:     {query.obs['method'].value_counts().to_dict()}")

    # Preprocess query to match the reference exactly:
    #   normalize -> log1p -> subset to ref HVGs -> scale
    # This ensures query cells are in the same feature space as the reference
    # before sc.tl.ingest projects them into the reference PCA.

    # Remove bacterial genes
    bacteria_genes = ["GQX67_00940", "GQX67_05945"] + [
        g for g in query.var_names if g.startswith("16S_")]
    query = query[:, ~query.var_names.isin(bacteria_genes)].copy()

    # Subset to the full pre-HVG gene universe (ref.raw has pre-HVG-subset genes)
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

    # Normalize + log1p to match reference preprocessing
    sc.pp.normalize_total(query, target_sum=1e4)
    sc.pp.log1p(query)
    if scipy.sparse.issparse(query.X):
        query.X = query.X.toarray()
    query.X = np.nan_to_num(query.X, nan=0.0, posinf=0.0, neginf=0.0)

    # Store normalized log counts
    query.raw = query

    # Subset AND reorder query vars to exactly match reference var_names.
    # sc.tl.ingest requires identical var_names in the same order.
    shared_genes = ref.var_names.intersection(query.var_names)
    missing_genes = ref.var_names[~ref.var_names.isin(query.var_names)]
    if len(missing_genes) > 0:
        print(f"  WARNING: {len(missing_genes)} ref HVGs missing from query — filling with zeros")
        # Add missing genes as zero columns directly in numpy — avoids ad.concat
        # which drops obs metadata. query.X is already dense at this point.
        zero_cols = np.zeros((query.n_obs, len(missing_genes)), dtype=np.float64)
        combined_X = np.hstack([query.X, zero_cols])
        combined_var = pd.DataFrame(
            index=list(query.var_names) + list(missing_genes)
        )
        # Preserve all obs metadata
        query = ad.AnnData(X=combined_X, obs=query.obs.copy(), var=combined_var)
    # Reorder to exactly match reference var_names (same genes, same order)
    query = query[:, ref.var_names].copy()
    sc.pp.scale(query, max_value=10)
    query.var["highly_variable"] = True

    # ingest: project into reference PCA + UMAP, transfer leiden labels.
    # sc.tl.ingest overwrites query.obs with only the transferred column,
    # so save full obs beforehand and merge back after.
    obs_backup = query.obs.copy()
    print("Running sc.tl.ingest …")
    sc.tl.ingest(query, ref, obs="leiden", embedding_method="umap")
    # Restore all original obs columns, keeping the transferred leiden label
    leiden_transferred = query.obs["leiden"].copy()
    query.obs = obs_backup
    query.obs["leiden_ref"] = leiden_transferred.values

    sc.pl.umap(query, color=["leiden_ref", "timepoint", "method"],
               save=f"_{sample}_query_ingested.pdf", ncols=3,
               title=["Cluster (transferred from ref)", "Timepoint", "Method"])

    # Joint reference + query object for plotting
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
# Stage 3 — Export SCEPTIC inputs
# ─────────────────────────────────────────────────────────────────────────────

def export_sceptic(ref, query, fig_dir, sample, use_pca=True, n_pcs=30):
    """Export SCEPTIC inputs. Pseudotime axis: DOX(0) -> D7 -> D28 -> D56 -> wMel(999)."""
    print("\n" + "=" * 60)
    print("STAGE 3 — Exporting SCEPTIC inputs")
    print("=" * 60)

    os.makedirs(fig_dir, exist_ok=True)

    ref_sceptic   = ref
    query_sceptic = query

    print(f"SCEPTIC — DOX reference cells (t=0)  : {ref_sceptic.n_obs}")
    print(f"SCEPTIC — wMel query cells (all t)   : {query_sceptic.n_obs}")
    print(query_sceptic.obs["timepoint_numeric"].value_counts().sort_index().to_string())

    if use_pca:
        # Project query into reference PCA space if it doesn't have X_pca.
        # sc.tl.ingest transfers UMAP but not PCA — we project manually using
        # the reference PCA loadings (varm["PCs"]).
        if "X_pca" not in query_sceptic.obsm:
            print("  Projecting query into reference PCA space …")
            if "PCs" not in ref_sceptic.varm:
                raise ValueError("Reference has no PCA loadings (varm['PCs']). "
                                 "Run sc.pp.pca on ref before export.")
            pca_loadings = ref_sceptic.varm["PCs"]          # genes × n_pcs
            # query.X is scaled to match ref — project: cells × genes @ genes × pcs
            Q = query_sceptic.X
            if scipy.sparse.issparse(Q):
                Q = Q.toarray()
            query_sceptic.obsm["X_pca"] = Q @ pca_loadings  # cells × n_pcs

        ref_n_pcs   = ref_sceptic.obsm["X_pca"].shape[1]
        query_n_pcs = query_sceptic.obsm["X_pca"].shape[1]
        actual_pcs  = min(n_pcs, ref_n_pcs, query_n_pcs)
        if actual_pcs < n_pcs:
            print(f"  WARNING: Using {actual_pcs} PCs (ref={ref_n_pcs}, query={query_n_pcs})")
        ref_mat      = ref_sceptic.obsm.get("X_pca_harmony", ref_sceptic.obsm["X_pca"])[:, :actual_pcs]
        query_mat    = query_sceptic.obsm.get("X_pca_harmony", query_sceptic.obsm["X_pca"])[:, :actual_pcs]
        feature_cols = [f"PC{i+1}" for i in range(actual_pcs)]
    else:
        hvg = ref_sceptic.var_names[ref_sceptic.var["highly_variable"]]
        shared_hvg = hvg[hvg.isin(query_sceptic.var_names)]
        ref_mat   = ref_sceptic[:, shared_hvg].X
        query_mat = query_sceptic[:, shared_hvg].X
        if scipy.sparse.issparse(ref_mat):   ref_mat   = ref_mat.toarray()
        if scipy.sparse.issparse(query_mat): query_mat = query_mat.toarray()
        feature_cols = list(shared_hvg)

    data_concat = np.vstack([ref_mat, query_mat])
    ref_labels   = ref_sceptic.obs["timepoint_numeric"].values.astype(int)
    query_labels = query_sceptic.obs["timepoint_numeric"].values.astype(int)
    labels       = np.concatenate([ref_labels, query_labels])
    label_list   = np.sort(np.unique(labels))
    cell_ids     = np.concatenate([ref_sceptic.obs_names, query_sceptic.obs_names])

    leiden_labels = np.concatenate([
        ref_sceptic.obs["leiden"].values,
        query_sceptic.obs["leiden_ref"].values,
    ])
    method_labels = np.concatenate([
        ref_sceptic.obs["method"].values,
        query_sceptic.obs["method"].values,
    ])

    mat_path  = os.path.join(fig_dir, f"sceptic_matrix_{sample}.csv")
    lab_path  = os.path.join(fig_dir, f"sceptic_labels_{sample}.csv")
    ll_path   = os.path.join(fig_dir, f"sceptic_label_list_{sample}.csv")
    meta_path = os.path.join(fig_dir, f"sceptic_metadata_{sample}.csv")

    pd.DataFrame(data_concat, index=cell_ids, columns=feature_cols).to_csv(mat_path)
    pd.Series(labels, index=cell_ids, name="timepoint").to_csv(lab_path)
    pd.Series(label_list, name="timepoint").to_csv(ll_path, index=False)
    pd.DataFrame({
        "cell_id":   cell_ids,
        "timepoint": labels,
        "leiden":    leiden_labels,
        "method":    method_labels,
    }).to_csv(meta_path, index=False)

    print(f"SCEPTIC matrix     -> {mat_path}  shape={data_concat.shape}")
    print(f"SCEPTIC labels     -> {lab_path}")
    print(f"SCEPTIC label_list -> {ll_path}  values={label_list}")
    print(f"SCEPTIC metadata   -> {meta_path}")

    return data_concat, labels, label_list


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
# Main
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
    sceptic_use_pca=True,
    ref_condition=None,  # str or list of str
):
    os.makedirs(fig_dir, exist_ok=True)
    sc.settings.figdir = fig_dir

    adatas_all = []
    for fp in files:
        adata = sc.read_h5ad(fp)
        adata.obs[batch_key] = os.path.splitext(os.path.basename(fp))[0]
        adatas_all.append(adata)

    combined_tmp = ad.concat(adatas_all, join="inner", merge="same", index_unique="-")
    combined_tmp = add_metadata(combined_tmp, batch_key)

    available_cell_lines = combined_tmp.obs["cell_line"].dropna().unique().tolist()

    # --ref_conditions accepts one or more cell_line values.
    # Default: JW18DOX only. Pass "JW18DOX JW18wMel" to use both controls.
    # In all cases only treatment=="Ctrl" rows are included in the reference;
    # SV timepoint rows (e.g. JW18DOX-SV-D1) are always excluded.
    if ref_condition is None:
        ref_conditions = ["JW18DOX"]
        print(f"No --ref_condition specified. Defaulting to: {ref_conditions}")
    else:
        # ref_condition may be a single string or list depending on nargs
        ref_conditions = ref_condition if isinstance(ref_condition, list) else [ref_condition]

    invalid = [c for c in ref_conditions if c not in available_cell_lines]
    if invalid:
        raise ValueError(
            f"--ref_condition {invalid} not found in data. "
            f"Available cell lines: {available_cell_lines}"
        )

    # Reference = specified cell lines AND treatment == Ctrl only
    ref_mask = (
        combined_tmp.obs["cell_line"].isin(ref_conditions) &
        (combined_tmp.obs["treatment"] == "Ctrl")
    )
    # Query = any cell line NOT in ref_conditions AND treatment != Ctrl
    # (i.e. the intermediate SV timepoints)
    query_mask = (
        ~combined_tmp.obs["cell_line"].isin(ref_conditions) |
        (combined_tmp.obs["cell_line"].isin(ref_conditions) &
         (combined_tmp.obs["treatment"] != "Ctrl"))
    )
    # But exclude non-ref cell lines that are stable controls from the query too —
    # only include SV timepoint rows as query
    query_mask = combined_tmp.obs["treatment"] == "SV"

    # Report what's excluded
    excluded = (
        combined_tmp.obs["cell_line"].isin(ref_conditions) &
        (combined_tmp.obs["treatment"] != "Ctrl")
    )
    if excluded.sum() > 0:
        print(f"  NOTE: Excluded {excluded.sum()} ref cell_line cells with treatment != Ctrl:")
        print(combined_tmp.obs.loc[excluded, "bio_condition"].value_counts().to_string())

    print(f"Reference conditions : {ref_conditions} (treatment=Ctrl only)")
    print(f"Reference breakdown  :")
    print(combined_tmp.obs.loc[ref_mask, "bio_condition"].value_counts().to_string())
    print(f"Reference cells      : {ref_mask.sum()}")
    print(f"Query (SV timepoints):")
    print(combined_tmp.obs.loc[query_mask, "bio_condition"].value_counts().to_string())
    print(f"Query cells          : {query_mask.sum()}")

    if query_mask.sum() == 0:
        raise ValueError("No query cells found (expected treatment=SV timepoint rows).")
    if ref_mask.sum() == 0:
        raise ValueError(f"No reference cells found for conditions: {ref_conditions}.")

    adata_ref   = combined_tmp[ref_mask].copy()
    adata_query = combined_tmp[query_mask].copy()

    ref, final_resolution = build_reference(
        adata_ref, batch_key=batch_key,
        min_genes=min_genes, min_cells=min_cells, n_pcs=n_pcs,
        fig_dir=fig_dir, sample=sample,
        optimize_resolution=optimize_resolution,
        resolutions=resolutions,
        leiden_resolution=leiden_resolution,
    )

    query, combined = map_query_to_reference(
        adata_query, ref=ref, batch_key=batch_key,
        min_genes=min_genes, n_pcs=n_pcs,
        fig_dir=fig_dir, sample=sample,
    )

    export_sceptic(ref, query, fig_dir=fig_dir, sample=sample,
                   use_pca=sceptic_use_pca, n_pcs=n_pcs)

    def _sanitize_obs(adata):
        for col in adata.obs.columns:
            if adata.obs[col].dtype == object or str(adata.obs[col].dtype) == "category":
                adata.obs[col] = adata.obs[col].astype(str).replace("nan", "NA")
        return adata

    ref      = _sanitize_obs(ref)
    query    = _sanitize_obs(query)
    combined = _sanitize_obs(combined)

    ref.write(out_path.replace(".h5ad", "_reference.h5ad"))
    query.write(out_path.replace(".h5ad", "_query.h5ad"))
    combined.write(out_path.replace(".h5ad", "_combined.h5ad"))

    print("\n" + "=" * 60)
    print("COMPLETE")
    print("=" * 60)
    print(f"Reference  -> {out_path.replace('.h5ad','_reference.h5ad')}")
    print(f"Query      -> {out_path.replace('.h5ad','_query.h5ad')}")
    print(f"Combined   -> {out_path.replace('.h5ad','_combined.h5ad')}")
    print(f"Figures    -> {fig_dir}/")
    print(f"Reference: {ref.n_obs} cells, {ref.obs['leiden'].nunique()} clusters "
          f"at resolution={final_resolution}")
    print(f"Query:     {query.n_obs} cells across "
          f"{query.obs['timepoint'].nunique()} timepoints")


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Staged scRNA-seq integration: method-only correction -> ingest -> SCEPTIC export"
    )
    parser.add_argument("--files",      required=True, nargs="+")
    parser.add_argument("--sample",     default="wolbachia_infection")
    parser.add_argument("--batch_key",  default="batch")
    parser.add_argument("--min_cells",  type=int, default=3)
    parser.add_argument("--min_genes",  type=int, default=200)
    parser.add_argument("--out_path",   default="integrated.h5ad")
    parser.add_argument("--fig_dir",    default="figures")
    parser.add_argument("--n_pcs",      type=int, default=30)
    parser.add_argument("--resolution", type=float, default=0.5)
    parser.add_argument("--optimize_resolution", action="store_true", default=True)
    parser.add_argument("--no_optimize_resolution", dest="optimize_resolution",
                        action="store_false")
    parser.add_argument("--resolutions", type=float, nargs="+", default=None)
    parser.add_argument("--sceptic_raw_counts", action="store_true", default=False)
    parser.add_argument("--ref_condition", type=str, nargs="+", default=None,
                        help="One or more cell_line values to use as reference "
                             "(default: JW18DOX). Example: --ref_condition JW18DOX JW18wMel")

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
        sceptic_use_pca=not args.sceptic_raw_counts,
        ref_condition=args.ref_condition,
    )


if __name__ == "__main__":
    main()