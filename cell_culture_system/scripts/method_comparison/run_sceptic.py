"""
run_sceptic.py
==============
Run SCEPTIC pseudotime inference on Wolbachia infection time-series data
and analyse how Wolbachia titer relates to pseudotime.

Inputs (from integrate.py Stage 3 outputs)
-------------------------------------------
  sceptic_matrix_{sample}.csv    : cells x PCA features
  sceptic_labels_{sample}.csv    : numeric timepoint per cell (0 = uninfected)
  sceptic_label_list_{sample}.csv: unique sorted timepoints
  sceptic_metadata_{sample}.csv  : leiden cluster, method, timepoint per cell

The metadata file is joined with titer values from the original h5ad files
(reference + query) to enable titer-vs-pseudotime analysis.

Outputs
-------
  sceptic_results_{sample}.csv          : per-cell pseudotime + metadata
  confusion_matrix_{sample}.pdf         : SCEPTIC classification performance
  pseudotime_violin_{sample}.pdf        : pseudotime distribution by timepoint
  pseudotime_by_cluster_{sample}.pdf    : pseudotime stratified by leiden cluster
  titer_vs_pseudotime_{sample}.pdf      : scatter + regression titer ~ pseudotime
  titer_vs_pseudotime_by_cluster_{sample}.pdf : same, faceted by cluster
  titer_vs_pseudotime_by_timepoint_{sample}.pdf: same, faceted by timepoint
  titer_by_pseudotime_bin_{sample}.pdf  : titer boxplot across pseudotime bins
  pseudotime_umap_{sample}.pdf          : UMAP coloured by pseudotime (if h5ad provided)
  sceptic_stats_{sample}.csv            : Spearman/Pearson correlations + Kruskal-Wallis
"""

import os
import argparse
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns
import scanpy as sc
from scipy import stats
from scipy.stats import spearmanr, pearsonr, kruskal
from sklearn.metrics import confusion_matrix

try:
    from sceptic import run_sceptic_and_evaluate
    from sceptic import plotting as sceptic_plotting
    from sceptic import evaluation as sceptic_evaluation
    SCEPTIC_AVAILABLE = True
except ImportError:
    SCEPTIC_AVAILABLE = False
    warnings.warn("sceptic not installed. Install with: pip install sceptic")


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _savefig(fig, path):
    fig.savefig(path, bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  Saved: {path}")


def _bin_pseudotime(pseudotime, n_bins=10):
    """Divide pseudotime into n equal-width bins, return bin labels."""
    bins = np.linspace(pseudotime.min(), pseudotime.max(), n_bins + 1)
    labels = [f"{bins[i]:.2f}-{bins[i+1]:.2f}" for i in range(n_bins)]
    return pd.cut(pseudotime, bins=bins, labels=labels, include_lowest=True)


# ─────────────────────────────────────────────────────────────────────────────
# Load SCEPTIC inputs
# ─────────────────────────────────────────────────────────────────────────────
# ─────────────────────────────────────────────────────────────────────────────
# Load SCEPTIC inputs from h5ad
# ─────────────────────────────────────────────────────────────────────────────

def load_from_h5ad(h5ad_path, pca_key="X_pca_harmony", timepoint_col="timepoint"):
    """
    Extract everything SCEPTIC needs directly from an h5ad object.

    Returns
    -------
    data       : np.ndarray  (cells x PCs)
    labels     : np.ndarray  (numeric timepoint per cell)
    label_list : np.ndarray  (unique sorted timepoints)
    metadata   : pd.DataFrame
    cell_ids   : list[str]
    """
    import scanpy as sc
    print(f"Loading h5ad: {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)

    # ── PCA embedding ────────────────────────────────────────────────────────
    if pca_key in adata.obsm:
        data = adata.obsm[pca_key]
        print(f"  Using embedding: {pca_key}  shape={data.shape}")
    elif "X_pca" in adata.obsm:
        data = adata.obsm["X_pca"]
        print(f"  WARNING: '{pca_key}' not found, falling back to X_pca  shape={data.shape}")
    else:
        raise KeyError(f"No PCA embedding found in adata.obsm. "
                       f"Available: {list(adata.obsm.keys())}")

    # ── Timepoint labels ─────────────────────────────────────────────────────
    if timepoint_col not in adata.obs.columns:
        raise KeyError(f"Column '{timepoint_col}' not found in adata.obs. "
                       f"Available: {list(adata.obs.columns)}")

    # Convert timepoint strings like "D7" → numeric 7; "0" / NaN → 0
    def _parse_tp(row):
        """Parse a row of adata.obs to get a numeric timepoint."""
        bio = str(row.get("bio_condition", "")).strip()
        tp  = str(row.get("timepoint", "")).strip()

        # ── Named controls (identified via bio_condition) ─────────────────
        if "wMel-Ctrl" in bio:
            return 999          # persistently infected control → latest timepoint
        if "DOX-Ctrl" in bio:
            return 0            # cured/uninfected control → timepoint 0

        # ── Standard D1, D7, D28 … format ────────────────────────────────
        s = tp.lstrip("Dd")
        try:
            return int(s)
        except ValueError:
            return 0

    labels = adata.obs.apply(_parse_tp, axis=1).values.astype(int)
    label_list = np.array(sorted(np.unique(labels)))

    print(f"  Cells:        {adata.n_obs}")
    print(f"  Features:     {data.shape[1]}")
    print(f"  Timepoints:   {label_list}")
    print(f"  Label counts: {pd.Series(labels).value_counts().sort_index().to_dict()}")

    # ── Metadata ─────────────────────────────────────────────────────────────
    keep_cols = [c for c in ["leiden", "method", "timepoint", "bio_condition",
                              "cell_line", "treatment", "replicate",
                              "wolbachia_titer", "phase", "cyclum_theta"]
                 if c in adata.obs.columns]
    metadata = adata.obs[keep_cols].copy()

    return data, labels, label_list, metadata, adata.obs_names.tolist()

# ─────────────────────────────────────────────────────────────────────────────
# Add Wolbachia titer from h5ad
# ─────────────────────────────────────────────────────────────────────────────

def add_titer_from_h5ad(metadata, ref_h5ad_path, query_h5ad_path):
    """
    Pull wolbachia_titer from the reference and query h5ad objects and
    join onto the metadata dataframe by cell barcode.
    """
    print("\nLoading Wolbachia titer from h5ad files …")
    titer_frames = []

    for path, label in [(ref_h5ad_path, "reference"), (query_h5ad_path, "query")]:
        if path and os.path.exists(path):
            adata = sc.read_h5ad(path)
            if "wolbachia_titer" in adata.obs.columns:
                titer_frames.append(
                    adata.obs[["wolbachia_titer"]].rename_axis("cell_id")
                )
                print(f"  Loaded titer from {label}: {adata.n_obs} cells")
            else:
                print(f"  WARNING: 'wolbachia_titer' not found in {label} h5ad")
        else:
            print(f"  Skipping {label}: path not provided or file not found")

    if not titer_frames:
        print("  No titer data found — titer analysis will be skipped")
        metadata["wolbachia_titer"] = np.nan
        return metadata

    titer_df = pd.concat(titer_frames)
    # metadata index is cell barcodes (set in load_sceptic_inputs)
    metadata = metadata.join(titer_df, how="left")
    n_with_titer = metadata["wolbachia_titer"].notna().sum()
    print(f"  Cells with titer data: {n_with_titer}/{len(metadata)}")
    return metadata


# ─────────────────────────────────────────────────────────────────────────────
# Run SCEPTIC
# ─────────────────────────────────────────────────────────────────────────────

def run_sceptic(data, labels, label_list, method="xgboost"):
    if not SCEPTIC_AVAILABLE:
        raise RuntimeError("sceptic is not installed. Run: pip install sceptic")

    print(f"\nRunning SCEPTIC (method={method}) …")
    cm_result, label_predicted, pseudotime, prob = run_sceptic_and_evaluate(
        data, labels, label_list=label_list, method=method
    )
    print(f"  Pseudotime range: {pseudotime.min():.4f} – {pseudotime.max():.4f}")
    print(f"  Mean pseudotime:  {pseudotime.mean():.4f}")
    return cm_result, label_predicted, pseudotime, prob


# ─────────────────────────────────────────────────────────────────────────────
# Plotting
# ─────────────────────────────────────────────────────────────────────────────

def plot_confusion_matrix(cm_result, label_list, fig_dir, sample):
    """SCEPTIC confusion matrix heatmap."""
    fig, ax = plt.subplots(figsize=(max(6, len(label_list)), max(5, len(label_list))))
    labels_str = [str(l) for l in label_list]
    sns.heatmap(cm_result, annot=True, fmt=".2f", cmap="Blues",
                xticklabels=labels_str, yticklabels=labels_str,
                ax=ax, cbar_kws={"label": "Proportion"})
    ax.set_xlabel("Predicted Timepoint")
    ax.set_ylabel("True Timepoint")
    ax.set_title(f"SCEPTIC Confusion Matrix — {sample}\n"
                 f"(row-normalised; diagonal = correct classification)")
    plt.tight_layout()
    _savefig(fig, os.path.join(fig_dir, f"confusion_matrix_{sample}.pdf"))


def plot_pseudotime_violin(pseudotime, true_labels, label_list, fig_dir, sample):
    """Pseudotime distribution per true timepoint."""
    df = pd.DataFrame({
        "pseudotime": pseudotime,
        "timepoint":  [str(t) for t in true_labels],
    })
    order = [str(l) for l in label_list]

    fig, ax = plt.subplots(figsize=(max(8, len(label_list) * 1.5), 5))
    sns.violinplot(data=df, x="timepoint", y="pseudotime",
                   order=order, palette="viridis", inner="box", ax=ax)
    ax.set_xlabel("True Timepoint (days post-infection; 0 = uninfected)")
    ax.set_ylabel("SCEPTIC Pseudotime")
    ax.set_title(f"Pseudotime distribution by timepoint — {sample}")

    # Spearman correlation between pseudotime and true timepoint
    rho, p = spearmanr(pseudotime, true_labels)
    ax.text(0.02, 0.97, f"Spearman rho={rho:.3f}  p={p:.2e}",
            transform=ax.transAxes, va="top", fontsize=9,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.7))
    plt.tight_layout()
    _savefig(fig, os.path.join(fig_dir, f"pseudotime_violin_{sample}.pdf"))
    return rho, p


def plot_pseudotime_by_cluster(pseudotime, metadata, fig_dir, sample):
    """Pseudotime distribution per leiden cluster."""
    leiden_col = "leiden_ref" if "leiden_ref" in metadata.columns else "leiden"
    if leiden_col not in metadata.columns:
        print("  Skipping cluster pseudotime plot: no leiden column in metadata")
        return

    df = pd.DataFrame({
        "pseudotime": pseudotime,
        "cluster":    metadata[leiden_col].values,
    })
    clusters = sorted(df["cluster"].unique())
    cmap = plt.cm.get_cmap("tab20")
    palette = {c: cmap(i % 20) for i, c in enumerate(clusters)}

    fig, ax = plt.subplots(figsize=(max(8, len(clusters) * 1.2), 5))
    sns.violinplot(data=df, x="cluster", y="pseudotime",
                   order=clusters, palette=palette, inner="box", ax=ax)
    ax.set_xlabel("Leiden Cluster")
    ax.set_ylabel("SCEPTIC Pseudotime")
    ax.set_title(f"Pseudotime by cluster — {sample}")

    # Kruskal-Wallis across clusters
    groups = [df[df["cluster"] == c]["pseudotime"].values for c in clusters]
    if len(groups) >= 2:
        h, p = kruskal(*groups)
        ax.text(0.02, 0.97, f"Kruskal-Wallis H={h:.2f}  p={p:.2e}",
                transform=ax.transAxes, va="top", fontsize=9,
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.7))
    plt.tight_layout()
    _savefig(fig, os.path.join(fig_dir, f"pseudotime_by_cluster_{sample}.pdf"))


def plot_titer_vs_pseudotime(pseudotime, metadata, fig_dir, sample, n_bins=10):
    """
    Core analysis: Wolbachia titer as a function of SCEPTIC pseudotime.
    Produces three complementary views.
    """
    if "wolbachia_titer" not in metadata.columns:
        print("  Skipping titer plots: no wolbachia_titer column")
        return {}

    df = pd.DataFrame({
        "pseudotime":      pseudotime,
        "wolbachia_titer": metadata["wolbachia_titer"].values,
        "timepoint":       metadata["timepoint"].values.astype(str),
    })

    leiden_col = "leiden_ref" if "leiden_ref" in metadata.columns else "leiden"
    if leiden_col in metadata.columns:
        df["cluster"] = metadata[leiden_col].values

    # Drop cells without titer
    df_titer = df.dropna(subset=["wolbachia_titer"])
    print(f"\n  Cells with titer for correlation: {len(df_titer)}")

    if len(df_titer) < 10:
        print("  WARNING: Too few cells with titer data for reliable analysis")
        return {}

    # ── Stats ────────────────────────────────────────────────────────────────
    rho_sp, p_sp = spearmanr(df_titer["pseudotime"], df_titer["wolbachia_titer"])
    rho_pe, p_pe = pearsonr( df_titer["pseudotime"], df_titer["wolbachia_titer"])
    slope, intercept, r_value, p_lin, se = stats.linregress(
        df_titer["pseudotime"], df_titer["wolbachia_titer"])
    print(f"  Spearman rho={rho_sp:.4f}  p={p_sp:.2e}")
    print(f"  Pearson  r  ={rho_pe:.4f}  p={p_pe:.2e}")
    print(f"  Linear regression: slope={slope:.4f}  r²={r_value**2:.4f}  p={p_lin:.2e}")

    stat_results = {
        "spearman_rho": rho_sp, "spearman_p": p_sp,
        "pearson_r":    rho_pe, "pearson_p":  p_pe,
        "slope":        slope,  "r_squared":  r_value**2,
        "linreg_p":     p_lin,
    }

    # ── Plot 1: scatter + regression line ────────────────────────────────────
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(df_titer["pseudotime"], df_titer["wolbachia_titer"],
               alpha=0.3, s=8, color="#2196F3", rasterized=True)
    # Regression line
    x_line = np.linspace(df_titer["pseudotime"].min(), df_titer["pseudotime"].max(), 200)
    ax.plot(x_line, slope * x_line + intercept, color="red", linewidth=2,
            label=f"Linear fit (r²={r_value**2:.3f})")
    # LOWESS smoother
    from statsmodels.nonparametric.smoothers_lowess import lowess
    smoothed = lowess(df_titer["wolbachia_titer"].values,
                      df_titer["pseudotime"].values, frac=0.3)
    ax.plot(smoothed[:, 0], smoothed[:, 1], color="orange", linewidth=2,
            linestyle="--", label="LOWESS smoother")
    ax.set_xlabel("SCEPTIC Pseudotime")
    ax.set_ylabel("Wolbachia Titer")
    ax.set_title(f"Wolbachia titer vs pseudotime — {sample}")
    ax.text(0.02, 0.97,
            f"Spearman rho={rho_sp:.3f} p={p_sp:.2e}\n"
            f"Pearson r={rho_pe:.3f} p={p_pe:.2e}",
            transform=ax.transAxes, va="top", fontsize=9,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))
    ax.legend(fontsize=9)
    plt.tight_layout()
    _savefig(fig, os.path.join(fig_dir, f"titer_vs_pseudotime_{sample}.pdf"))

    # ── Plot 2: binned boxplot ────────────────────────────────────────────────
    df_titer = df_titer.copy()
    df_titer["pt_bin"] = _bin_pseudotime(df_titer["pseudotime"], n_bins=n_bins)

    fig, ax = plt.subplots(figsize=(12, 5))
    sns.boxplot(data=df_titer, x="pt_bin", y="wolbachia_titer",
                ax=ax, color="#90CAF9", flierprops=dict(markersize=2))
    ax.set_xlabel(f"Pseudotime bin (n={n_bins})")
    ax.set_ylabel("Wolbachia Titer")
    ax.set_title(f"Wolbachia titer across pseudotime bins — {sample}")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    _savefig(fig, os.path.join(fig_dir, f"titer_by_pseudotime_bin_{sample}.pdf"))

    # ── Plot 3: coloured by timepoint ─────────────────────────────────────────
    timepoints = sorted(df_titer["timepoint"].unique())
    tp_palette = dict(zip(timepoints, sns.color_palette("viridis", len(timepoints))))

    fig, ax = plt.subplots(figsize=(9, 6))
    for tp in timepoints:
        sub = df_titer[df_titer["timepoint"] == tp]
        ax.scatter(sub["pseudotime"], sub["wolbachia_titer"],
                   label=f"D{tp}" if tp != "0" else "uninfected",
                   color=tp_palette[tp], alpha=0.5, s=8, rasterized=True)
    ax.set_xlabel("SCEPTIC Pseudotime")
    ax.set_ylabel("Wolbachia Titer")
    ax.set_title(f"Wolbachia titer vs pseudotime — coloured by timepoint — {sample}")
    ax.legend(title="Timepoint", bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)
    plt.tight_layout()
    _savefig(fig, os.path.join(fig_dir, f"titer_vs_pseudotime_by_timepoint_{sample}.pdf"))

    # ── Plot 4: faceted by cluster (if available) ─────────────────────────────
    if "cluster" in df_titer.columns:
        clusters = sorted(df_titer["cluster"].unique())
        ncols = min(4, len(clusters))
        nrows = int(np.ceil(len(clusters) / ncols))
        fig, axes = plt.subplots(nrows, ncols,
                                 figsize=(ncols * 4, nrows * 3.5), sharey=True)
        axes = np.array(axes).flatten()
        cmap = plt.cm.get_cmap("tab20")

        for i, cluster in enumerate(clusters):
            ax = axes[i]
            sub = df_titer[df_titer["cluster"] == cluster]
            ax.scatter(sub["pseudotime"], sub["wolbachia_titer"],
                       color=cmap(i % 20), alpha=0.4, s=6, rasterized=True)
            if len(sub) >= 5:
                s, p_s = spearmanr(sub["pseudotime"], sub["wolbachia_titer"])
                # Mini regression line
                sl, ic, _, _, _ = stats.linregress(sub["pseudotime"],
                                                    sub["wolbachia_titer"])
                xl = np.linspace(sub["pseudotime"].min(), sub["pseudotime"].max(), 100)
                ax.plot(xl, sl * xl + ic, color="red", linewidth=1.5)
                ax.text(0.05, 0.95, f"rho={s:.2f} p={p_s:.1e}",
                        transform=ax.transAxes, va="top", fontsize=7,
                        bbox=dict(boxstyle="round", facecolor="white", alpha=0.7))
            ax.set_title(f"Cluster {cluster}", fontsize=9)
            ax.set_xlabel("Pseudotime", fontsize=8)
            if i % ncols == 0:
                ax.set_ylabel("Wolbachia Titer", fontsize=8)

        for j in range(i + 1, len(axes)):
            axes[j].set_visible(False)

        plt.suptitle(f"Titer vs pseudotime by cluster — {sample}", fontweight="bold")
        plt.tight_layout()
        _savefig(fig, os.path.join(fig_dir,
                                   f"titer_vs_pseudotime_by_cluster_{sample}.pdf"))

    return stat_results


def plot_pseudotime_on_umap(pseudotime, metadata, h5ad_path, fig_dir, sample):
    """Project pseudotime back onto the existing UMAP embedding."""
    if not h5ad_path or not os.path.exists(h5ad_path):
        print("  Skipping UMAP pseudotime plot: no h5ad path provided")
        return

    print(f"  Loading h5ad for UMAP: {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)

    if "X_umap" not in adata.obsm:
        print("  Skipping UMAP pseudotime plot: no X_umap in h5ad")
        return

    # Match pseudotime to cells present in this h5ad
    pt_series = pd.Series(pseudotime, index=metadata.index)
    common    = adata.obs_names.intersection(pt_series.index)
    if len(common) == 0:
        print("  Skipping UMAP pseudotime plot: no matching barcodes")
        return

    adata_sub = adata[common].copy()
    adata_sub.obs["sceptic_pseudotime"] = pt_series[common].values

    sc.settings.figdir = fig_dir
    sc.pl.umap(adata_sub, color="sceptic_pseudotime", cmap="viridis",
               save=f"_pseudotime_{sample}.pdf",
               title=f"SCEPTIC Pseudotime — {sample}")
    print(f"  Saved: {fig_dir}/umap_pseudotime_{sample}.pdf")


# ─────────────────────────────────────────────────────────────────────────────
# Save results
# ─────────────────────────────────────────────────────────────────────────────

def save_results(pseudotime, label_predicted, prob, labels, label_list,
                 metadata, cell_ids, stat_results, fig_dir, sample):
    """Write per-cell results and summary statistics to CSV."""
    results_df = pd.DataFrame({
        "cell_id":         cell_ids,
        "pseudotime":      pseudotime,
        "true_timepoint":  labels,
        "pred_timepoint":  label_predicted,
    })
    # Add probability columns
    for i, tp in enumerate(label_list):
        results_df[f"prob_t{tp}"] = prob[:, i]

    # Merge metadata
    metadata_reset = metadata.reset_index(drop=True)
    results_df = pd.concat([results_df, metadata_reset], axis=1)

    results_path = os.path.join(fig_dir, f"sceptic_results_{sample}.csv")
    results_df.to_csv(results_path, index=False)
    print(f"\n  Per-cell results -> {results_path}")

    if stat_results:
        stats_df = pd.DataFrame([stat_results])
        stats_path = os.path.join(fig_dir, f"sceptic_stats_{sample}.csv")
        stats_df.to_csv(stats_path, index=False)
        print(f"  Statistics       -> {stats_path}")

    return results_df


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Run SCEPTIC pseudotime and analyse Wolbachia titer",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--h5ad",         required=True,
                        help="Annotated h5ad (output of cell_cycle_analysis rule)")
    parser.add_argument("--sample",       default="wolbachia_infection")
    parser.add_argument("--fig_dir",      default="results/sceptic")
    parser.add_argument("--pca_key",      default="X_pca_harmony",
                        help="obsm key for PCA embedding (default: X_pca_harmony)")
    parser.add_argument("--timepoint_col", default="timepoint",
                        help="obs column with timepoint labels (default: timepoint)")
    parser.add_argument("--method",       default="xgboost",
                        choices=["xgboost", "svm"])
    parser.add_argument("--n_bins",       type=int, default=10)

    args = parser.parse_args()
    os.makedirs(args.fig_dir, exist_ok=True)

    # ── Load ─────────────────────────────────────────────────────────────────
    data, labels, label_list, metadata, cell_ids = load_from_h5ad(
        args.h5ad,
        pca_key=args.pca_key,
        timepoint_col=args.timepoint_col,
    )

    # titer is already in metadata if present in adata.obs — no separate h5ad needed

    # ── Run SCEPTIC ──────────────────────────────────────────────────────────
    cm_result, label_predicted, pseudotime, prob = run_sceptic(
        data, labels, label_list, method=args.method
    )

    # ── Plots ─────────────────────────────────────────────────────────────────
    print("\nGenerating plots …")
    plot_confusion_matrix(cm_result, label_list, args.fig_dir, args.sample)

    rho_sp, p_sp = plot_pseudotime_violin(
        pseudotime, labels, label_list, args.fig_dir, args.sample)

    plot_pseudotime_by_cluster(pseudotime, metadata, args.fig_dir, args.sample)

    stat_results = plot_titer_vs_pseudotime(
        pseudotime, metadata, args.fig_dir, args.sample, n_bins=args.n_bins)

    stat_results["pseudotime_timepoint_spearman_rho"] = rho_sp
    stat_results["pseudotime_timepoint_spearman_p"]   = p_sp

    # UMAP projection uses the same h5ad
    plot_pseudotime_on_umap(pseudotime, metadata, args.h5ad,
                            args.fig_dir, args.sample)

    # ── Save ──────────────────────────────────────────────────────────────────
    save_results(
        pseudotime, label_predicted, prob, labels, label_list,
        metadata, cell_ids, stat_results, args.fig_dir, args.sample
    )

    print("\n" + "=" * 60)
    print("SCEPTIC ANALYSIS COMPLETE")
    print("=" * 60)
    print(f"Pseudotime range : {pseudotime.min():.4f} – {pseudotime.max():.4f}")
    print(f"Pseudotime ~ true timepoint: Spearman rho={rho_sp:.3f}  p={p_sp:.2e}")
    if stat_results:
        print(f"Titer ~ pseudotime:          "
              f"Spearman rho={stat_results.get('spearman_rho', np.nan):.3f}  "
              f"p={stat_results.get('spearman_p', np.nan):.2e}")
    print(f"\nOutputs -> {args.fig_dir}/")


if __name__ == "__main__":
    main()