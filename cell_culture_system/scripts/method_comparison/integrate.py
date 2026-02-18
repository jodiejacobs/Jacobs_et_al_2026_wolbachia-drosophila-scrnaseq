"""
integrate.py
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

 python scripts/method_comparison/integrate.py \
    --files results/filtered_h5ad/*.h5ad \
    --sample wolbachia_infection \
    --batch_key batch \
    --min_cells 3 \
    --min_genes 200 \
    --n_pcs 30 \
    --out_path results/integrated/integrated.h5ad \
    --fig_dir results/integrated/figures \
    --optimize_resolution \
    --resolutions 0.1 0.2 0.3 0.4 0.5 0.6 0.8 1.0 1.2 1.5
"""

import os
import argparse
import warnings

import numpy as np
import pandas as pd
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
    cmap = plt.cm.get_cmap("tab20")
    return [cmap(i % 20) for i in range(len(clusters))]


def _extract_timepoint_numeric(timepoint_str):
    """Convert 'D7' -> 7, 'D28' -> 28, NaN -> 0 (uninfected)."""
    if pd.isna(timepoint_str):
        return 0
    import re
    m = re.search(r"(\d+)", str(timepoint_str))
    return int(m.group(1)) if m else 0


# ─────────────────────────────────────────────────────────────────────────────
# Preprocessing (shared)
# ─────────────────────────────────────────────────────────────────────────────

def add_metadata(adata, batch_key):
    """Extract sample metadata from batch names."""
    adata.obs["cell_line"]  = adata.obs[batch_key].str.extract(r"(JW18DOX|JW18wMel)")[0]
    adata.obs["treatment"]  = adata.obs[batch_key].str.extract(r"-(Ctrl|SV)")[0]
    adata.obs["timepoint"]  = adata.obs[batch_key].str.extract(r"-(D\d+)-")[0]
    adata.obs["replicate"]  = adata.obs[batch_key].str.extract(r"-(Ctrl|SV)-([^_]+)")[1]
    adata.obs["method"]     = adata.obs[batch_key].str.extract(r"_(10x|pipseq)$")[0]
    adata.obs["bio_condition"] = adata.obs.apply(
        lambda row: (f"{row['cell_line']}-{row['treatment']}-{row['timepoint']}"
                     if pd.notna(row["timepoint"])
                     else f"{row['cell_line']}-{row['treatment']}"),
        axis=1,
    )
    adata.obs["timepoint_numeric"] = adata.obs["timepoint"].apply(
        _extract_timepoint_numeric).astype(int)
    return adata


def preprocess(adata, min_genes, min_cells, n_pcs, n_top_genes=2000):
    """Filter, HVG, PCA. Returns new adata."""
    bacteria_genes = ["GQX67_00940", "GQX67_05945"] + [
        g for g in adata.var_names if g.startswith("16S_")]
    adata = adata[:, ~adata.var_names.isin(bacteria_genes)].copy()
    if scipy.sparse.issparse(adata.X):
        adata.X = adata.X.toarray()
    adata.X = np.nan_to_num(adata.X, nan=0.0)
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=n_top_genes)
    sc.pp.pca(adata, n_comps=n_pcs)
    return adata


# ─────────────────────────────────────────────────────────────────────────────
# Resolution optimisation (runs on reference only)
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
    X_pca = adata.obsm["X_pca"][:, :n_pcs]
    rng = np.random.default_rng(random_state)
    sil_idx = rng.choice(adata.n_obs, min(silhouette_subsample, adata.n_obs), replace=False)

    records = []
    cluster_assignments = {}

    print(f"\nSweeping {len(resolutions)} resolutions on reference …")
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

    # Figure 1: metric curves
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
        ax.set_xlabel("Leiden Resolution")
        ax.set_title(title)
        ax.legend(fontsize=9); ax.grid(alpha=0.3)
    plt.suptitle(f"Resolution Optimisation — {sample}", fontweight="bold", y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f"resolution_metrics_{sample}.pdf"), bbox_inches="tight")
    plt.close()

    # Figure 2: UMAP panel
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
                   cmap=plt.cm.get_cmap("tab20", n_cl), s=1, alpha=0.5, rasterized=True)
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

    # Figure 3: overlap heatmaps between adjacent resolutions
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
    plt.savefig(os.path.join(fig_dir, f"resolution_overlap_heatmaps_{sample}.pdf"),
                bbox_inches="tight")
    plt.close()

    # Figure 4: composite bar chart
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
# Stage 1 — Reference: uninfected controls
# ─────────────────────────────────────────────────────────────────────────────

def build_reference(adata_ref, batch_key, min_genes, min_cells, n_pcs,
                    fig_dir, sample, optimize_resolution=True,
                    resolutions=None, leiden_resolution=0.5):
    """
    Build reference embedding from uninfected control cells.

    Corrects ONLY for library-prep method (BBKNN key='method').
    Biological variation (cell state, infection history) is NOT corrected —
    it remains visible in the embedding and drives clustering.
    """
    print("\n" + "=" * 60)
    print("STAGE 1 — Building reference from uninfected controls")
    print("=" * 60)

    ref = adata_ref.copy()
    print(f"Reference cells: {ref.n_obs}")
    print(f"Methods: {ref.obs['method'].value_counts().to_dict()}")

    ref = preprocess(ref, min_genes, min_cells, n_pcs)

    # Pre-correction UMAP for QC
    ref_tmp = ref.copy()
    sc.pp.neighbors(ref_tmp, n_pcs=n_pcs)
    sc.tl.umap(ref_tmp)
    sc.pl.umap(ref_tmp, color=["method", "bio_condition"],
               save=f"_{sample}_ref_before_correction.pdf", ncols=2,
               title=["Method (pre-correction)", "Bio condition (pre-correction)"])

    # Correct ONLY for library-prep method — NOT for biological condition
    print("Correcting for library-prep method only (BBKNN key='method') …")
    bbknn.bbknn(ref, batch_key="method", n_pcs=n_pcs, neighbors_within_batch=5)
    sc.tl.umap(ref)

    # Post-correction QC: methods should overlap, bio conditions should still separate
    sc.pl.umap(ref, color=["method", "bio_condition", "cell_line"],
               save=f"_{sample}_ref_after_method_correction.pdf", ncols=3,
               title=["Method (post-correction — should overlap)",
                      "Bio condition (should still separate)",
                      "Cell line"])

    # Optimise Leiden resolution on the reference
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
# Stage 2 — Query: map new-infection cells onto reference via ingest
# ─────────────────────────────────────────────────────────────────────────────

def map_query_to_reference(adata_query, ref, batch_key, min_genes, n_pcs,
                            fig_dir, sample):
    """
    Project new-infection timepoint cells onto the reference embedding using
    sc.tl.ingest.

    ingest:
      - Projects query cells into the reference PCA space
      - Transfers Leiden cluster labels from nearest reference neighbours
      - Does NOT re-cluster — cluster definitions remain anchored in the
        uninfected cell state landscape

    This means: if a D7-infected cell lands in cluster 3, it was in a cell
    state that exists in the uninfected reference; if it falls between
    clusters, it may be transitioning into a new infection-induced state.
    """
    print("\n" + "=" * 60)
    print("STAGE 2 — Mapping query (new infection timepoints) onto reference")
    print("=" * 60)

    query = adata_query.copy()
    print(f"Query cells: {query.n_obs}")
    print(f"Timepoints:  {query.obs['timepoint'].value_counts().to_dict()}")
    print(f"Methods:     {query.obs['method'].value_counts().to_dict()}")

    # Subset query to the same genes as reference, then filter cells
    query = query[:, query.var_names.isin(ref.var_names)].copy()
    if scipy.sparse.issparse(query.X):
        query.X = query.X.toarray()
    query.X = np.nan_to_num(query.X, nan=0.0)
    sc.pp.filter_cells(query, min_genes=min_genes)

    # Copy HVG flags from reference so ingest uses the same gene set
    query.var["highly_variable"] = query.var_names.isin(
        ref.var_names[ref.var["highly_variable"]])

    # ingest: project into reference PCA + UMAP, transfer leiden labels
    print("Running sc.tl.ingest …")
    sc.tl.ingest(query, ref, obs="leiden", embedding_method="umap")
    # ingest writes the transferred column under the same name; rename for clarity
    query.obs.rename(columns={"leiden": "leiden_ref"}, inplace=True)

    sc.pl.umap(query, color=["leiden_ref", "timepoint", "method"],
               save=f"_{sample}_query_ingested.pdf", ncols=3,
               title=["Cluster (transferred from ref)", "Timepoint", "Method"])

    # ── Joint reference + query object for plotting ───────────────────────────
    ref_plot = ref.copy()
    ref_plot.obs["dataset"]    = "uninfected (reference)"
    ref_plot.obs["leiden_ref"] = ref_plot.obs["leiden"]
    query.obs["dataset"]       = "new infection (" + query.obs["timepoint"].fillna("?") + ")"

    common_genes = ref_plot.var_names.intersection(query.var_names)
    combined = ad.concat(
        [ref_plot[:, common_genes], query[:, common_genes]],
        join="inner", index_unique="-",
    )

    sc.pl.umap(combined, color=["dataset", "leiden_ref", "timepoint", "method"],
               save=f"_{sample}_combined_ref_query.pdf", ncols=2,
               title=["Dataset (ref vs query)", "Cluster", "Timepoint", "Method"])

    # How does cluster composition change with infection?
    _plot_cluster_by_timepoint(combined, fig_dir, sample)

    # Validate method correction didn't bleed into biological signal
    _plot_method_validation(ref, query, fig_dir, sample)

    return query, combined


# ─────────────────────────────────────────────────────────────────────────────
# Stage 3 — Export SCEPTIC-ready files
# ─────────────────────────────────────────────────────────────────────────────

def export_sceptic(ref, query, fig_dir, sample, use_pca=True, n_pcs=30):
    """
    Prepare and save all inputs needed by SCEPTIC.

    SCEPTIC API:
        run_sceptic_and_evaluate(data, time_labels, label_list, method='xgboost')
        data         : cells x features  (numpy array)
        time_labels  : numeric timepoint per cell (0 = uninfected)
        label_list   : unique sorted timepoints

    We include both the uninfected reference cells (t=0) and all new-infection
    timepoint cells so SCEPTIC can learn the full time axis.

    Feature matrix options:
        use_pca=True  -> PCA coordinates (recommended: compact, de-noised,
                         and already in the shared method-corrected space)
        use_pca=False -> raw HVG counts (higher dimensional, noisier)
    """
    print("\n" + "=" * 60)
    print("STAGE 3 — Exporting SCEPTIC inputs")
    print("=" * 60)

    os.makedirs(fig_dir, exist_ok=True)

    # Use wMel cells only for pseudotime: uninfected wMel controls + new infections
    # Exclude JW18DOX cells — different perturbation, separate pseudotime analysis
    ref_wmel   = ref[ref.obs["cell_line"] == "JW18wMel"].copy()
    query_wmel = query[query.obs["cell_line"] == "JW18wMel"].copy()

    print(f"wMel reference cells  (t=0)  : {ref_wmel.n_obs}")
    print(f"wMel query cells (timepoints): {query_wmel.n_obs}")
    print(f"  timepoints: {query_wmel.obs['timepoint'].value_counts().to_dict()}")

    # Feature matrix
    if use_pca:
        # ingest projects query into the reference PCA space, so both share the same dims
        ref_mat   = ref_wmel.obsm["X_pca"][:, :n_pcs]
        query_mat = query_wmel.obsm["X_pca"][:, :n_pcs]
        feature_cols = [f"PC{i+1}" for i in range(n_pcs)]
    else:
        hvg = ref_wmel.var_names[ref_wmel.var["highly_variable"]]
        shared_hvg = hvg[hvg.isin(query_wmel.var_names)]
        ref_mat   = ref_wmel[:, shared_hvg].X
        query_mat = query_wmel[:, shared_hvg].X
        if scipy.sparse.issparse(ref_mat):   ref_mat   = ref_mat.toarray()
        if scipy.sparse.issparse(query_mat): query_mat = query_mat.toarray()
        feature_cols = list(shared_hvg)

    data_concat = np.vstack([ref_mat, query_mat])

    # Labels
    ref_labels   = ref_wmel.obs["timepoint_numeric"].values.astype(int)
    query_labels = query_wmel.obs["timepoint_numeric"].values.astype(int)
    labels       = np.concatenate([ref_labels, query_labels])
    label_list   = np.sort(np.unique(labels))
    cell_ids     = np.concatenate([ref_wmel.obs_names, query_wmel.obs_names])

    # Metadata: cluster + method per cell (useful for stratified SCEPTIC analysis)
    leiden_labels = np.concatenate([
        ref_wmel.obs["leiden"].values,
        query_wmel.obs["leiden_ref"].values,
    ])
    method_labels = np.concatenate([
        ref_wmel.obs["method"].values,
        query_wmel.obs["method"].values,
    ])

    # Save
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

    print(f"SCEPTIC matrix    -> {mat_path}  shape={data_concat.shape}")
    print(f"SCEPTIC labels    -> {lab_path}")
    print(f"SCEPTIC label_list -> {ll_path}  values={label_list}")
    print(f"SCEPTIC metadata  -> {meta_path}")
    print("\n" + "-" * 50)
    print("Example SCEPTIC usage:")
    print("  from sceptic import run_sceptic_and_evaluate")
    print("  import pandas as pd, numpy as np")
    print(f"  data       = pd.read_csv('{mat_path}', index_col=0).values")
    print(f"  labels     = pd.read_csv('{lab_path}', index_col=0)['timepoint'].values")
    print(f"  label_list = pd.read_csv('{ll_path}')['timepoint'].values")
    print("  cm, pred, pseudotime, prob = run_sceptic_and_evaluate(")
    print("      data, labels, label_list=label_list, method='xgboost')")
    print("  # For stratified analysis by cluster:")
    print("  meta = pd.read_csv('{meta_path}')")
    print("  from sceptic import plotting")
    print("  plotting.plot_pseudotime_by_group(pseudotime, labels, meta['leiden'].values,")
    print("      output_dir='pseudotime_by_cluster')")

    return data_concat, labels, label_list


# ─────────────────────────────────────────────────────────────────────────────
# Plotting helpers
# ─────────────────────────────────────────────────────────────────────────────

def _plot_cluster_by_timepoint(combined, fig_dir, sample):
    """Show how cluster composition shifts across infection timepoints."""
    cluster_col = "leiden_ref" if "leiden_ref" in combined.obs.columns else "leiden"
    colors = _leiden_colors(combined, key=cluster_col)
    timepoints = sorted(combined.obs["timepoint_numeric"].unique())

    # Stacked bar: % cluster composition per timepoint
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
    plt.savefig(os.path.join(fig_dir, f"cluster_by_timepoint_{sample}.pdf"),
                bbox_inches="tight")
    plt.close()

    # UMAP panels: one per timepoint, highlighting those cells
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
                   c=all_labels[mask], cmap=plt.cm.get_cmap("tab20", n_cl),
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

    # Titer by timepoint and cluster
    if "wolbachia_titer" in combined.obs.columns:
        fig, ax = plt.subplots(figsize=(12, 5))
        sns.boxplot(data=combined.obs, x="timepoint_numeric", y="wolbachia_titer",
                    hue=cluster_col, ax=ax, flierprops=dict(markersize=1))
        ax.set_xlabel("Timepoint (days post-infection)")
        ax.set_ylabel("Wolbachia Titer")
        ax.set_title("Wolbachia titer by timepoint and cluster")
        ax.legend(title="Cluster", bbox_to_anchor=(1.05, 1), loc="upper left")
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f"titer_by_timepoint_{sample}.pdf"),
                    bbox_inches="tight")
        plt.close()


def _plot_method_validation(ref, query, fig_dir, sample):
    """
    Key QC figure: confirm that method correction removed 10X/PIPseq technical
    variation without erasing biological signal.

    PASS criteria:
        - 10X and PIPseq cells overlap in UMAP space (method correction worked)
        - Biological conditions (uninfected vs timepoints) still separate
    """
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    cluster_col_map = {0: "leiden", 1: "leiden_ref"}

    for row, (adata, title_prefix) in enumerate([(ref, "Reference"), (query, "Query")]):
        umap_xy = adata.obsm["X_umap"]
        cluster_col = "leiden" if "leiden" in adata.obs.columns else "leiden_ref"

        # Method
        method_cmap = {"10x": "#1f77b4", "pipseq": "#ff7f0e"}
        c = [method_cmap.get(m, "grey") for m in adata.obs["method"].values]
        axes[row, 0].scatter(umap_xy[:, 0], umap_xy[:, 1], c=c, s=1, alpha=0.5, rasterized=True)
        axes[row, 0].set_title(f"{title_prefix}: Method (should overlap post-correction)")
        for m, col in method_cmap.items():
            axes[row, 0].scatter([], [], c=col, label=m)
        axes[row, 0].legend(); axes[row, 0].set_xticks([]); axes[row, 0].set_yticks([])

        # Bio condition
        conditions = adata.obs["bio_condition"].astype("category")
        cond_cmap = plt.cm.get_cmap("Set2", len(conditions.cat.categories))
        c = [cond_cmap(conditions.cat.codes.iloc[i]) for i in range(len(conditions))]
        axes[row, 1].scatter(umap_xy[:, 0], umap_xy[:, 1], c=c, s=1, alpha=0.5, rasterized=True)
        axes[row, 1].set_title(f"{title_prefix}: Bio condition (should still separate)")
        for i, cat in enumerate(conditions.cat.categories):
            axes[row, 1].scatter([], [], c=cond_cmap(i), label=cat)
        axes[row, 1].legend(fontsize=7); axes[row, 1].set_xticks([]); axes[row, 1].set_yticks([])

        # Clusters
        lbls = adata.obs[cluster_col].astype(int).values
        axes[row, 2].scatter(umap_xy[:, 0], umap_xy[:, 1], c=lbls,
                             cmap=plt.cm.get_cmap("tab20", len(np.unique(lbls))),
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
# Main entry point
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
):
    os.makedirs(fig_dir, exist_ok=True)
    sc.settings.figdir = fig_dir

    # Load all files
    adatas_all = []
    for fp in files:
        adata = sc.read_h5ad(fp)
        adata.obs[batch_key] = os.path.splitext(os.path.basename(fp))[0]
        adatas_all.append(adata)

    # Add metadata to a temporary combined object to classify cells
    combined_tmp = ad.concat(adatas_all, join="inner", merge="same", index_unique="-")
    combined_tmp = add_metadata(combined_tmp, batch_key)

    # Split: reference = uninfected (timepoint_numeric == 0)
    #        query     = new infection timepoints (timepoint_numeric > 0)
    ref_mask   = combined_tmp.obs["timepoint_numeric"] == 0
    query_mask = combined_tmp.obs["timepoint_numeric"] >  0

    print(f"Reference (uninfected) cells : {ref_mask.sum()}")
    print(f"Query (new infection) cells  : {query_mask.sum()}")

    if query_mask.sum() == 0:
        raise ValueError(
            "No query cells found (timepoint > 0). "
            "Check that infection timepoint labels (D7, D28, …) appear in filenames."
        )

    # Slice combined_tmp directly — avoids barcode mismatch caused by
    # index_unique="-" modifying barcodes during the initial concat
    adata_ref   = combined_tmp[ref_mask].copy()
    adata_query = combined_tmp[query_mask].copy()

    # Stage 1: build reference
    ref, final_resolution = build_reference(
        adata_ref, batch_key=batch_key,
        min_genes=min_genes, min_cells=min_cells, n_pcs=n_pcs,
        fig_dir=fig_dir, sample=sample,
        optimize_resolution=optimize_resolution,
        resolutions=resolutions,
        leiden_resolution=leiden_resolution,
    )

    # Stage 2: ingest query onto reference
    query, combined = map_query_to_reference(
        adata_query, ref=ref, batch_key=batch_key,
        min_genes=min_genes, n_pcs=n_pcs,
        fig_dir=fig_dir, sample=sample,
    )

    # Stage 3: export SCEPTIC inputs
    export_sceptic(ref, query, fig_dir=fig_dir, sample=sample,
                   use_pca=sceptic_use_pca, n_pcs=n_pcs)

    # Save objects
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
    print(f"\nReference: {ref.n_obs} cells, {ref.obs['leiden'].nunique()} clusters "
          f"at resolution={final_resolution}")
    print(f"Query:     {query.n_obs} cells across "
          f"{query.obs['timepoint'].nunique()} timepoints")


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Staged scRNA-seq integration: method-only correction -> "
            "ingest query onto reference -> SCEPTIC export"
        )
    )
    parser.add_argument("--files", required=True, nargs="+",
                        help="h5ad files to integrate")
    parser.add_argument("--sample", default="wolbachia_infection",
                        help="Label used in all output filenames")
    parser.add_argument("--batch_key", default="batch",
                        help="obs column for sample/batch identity")
    parser.add_argument("--min_cells", type=int, default=3)
    parser.add_argument("--min_genes",  type=int, default=200)
    parser.add_argument("--out_path",   default="integrated.h5ad",
                        help="Base path; _reference/_query/_combined suffixes added")
    parser.add_argument("--fig_dir",    default="figures")
    parser.add_argument("--n_pcs",      type=int, default=30)
    parser.add_argument("--resolution", type=float, default=0.5,
                        help="Fixed Leiden resolution (ignored if --optimize_resolution)")
    parser.add_argument("--optimize_resolution", default=True, action="store_true")
    parser.add_argument("--no_optimize_resolution", dest="optimize_resolution",
                        action="store_false",
                        help="Skip resolution sweep; use --resolution directly")
    parser.add_argument("--resolutions", type=float, nargs="+", default=None,
                        help="Custom resolution sweep values (default: 0.1 to 2.0)")
    parser.add_argument("--sceptic_raw_counts", action="store_true", default=False,
                        help="Export raw HVG counts for SCEPTIC instead of PCA coordinates")

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
    )


if __name__ == "__main__":
    main()
