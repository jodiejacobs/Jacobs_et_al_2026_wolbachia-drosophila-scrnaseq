"""
run_sceptic.py
==============
Run SCEPTIC pseudotime inference on Wolbachia infection time-series data
and analyse how Wolbachia titer relates to pseudotime.

Inputs
------
  adata_with_programs.h5ad : h5ad containing NMF program scores (Program_0 …
                              Program_N), PCA embedding, timepoint_numeric, and
                              wolbachia_titer in adata.obs.

Outputs
-------
  sceptic_results_{sample}.csv              : per-cell pseudotime + metadata
  confusion_matrix_{sample}.pdf             : SCEPTIC classification performance
  pseudotime_violin_{sample}.pdf            : pseudotime distribution by timepoint
  pseudotime_by_program_{sample}.pdf        : pseudotime stratified by dominant NMF program
  titer_vs_pseudotime_{sample}.pdf          : scatter + regression titer ~ pseudotime
  titer_vs_pseudotime_by_program_{sample}.pdf : same, faceted by dominant NMF program
  titer_vs_pseudotime_by_timepoint_{sample}.pdf: same, coloured by timepoint
  titer_by_pseudotime_bin_{sample}.pdf      : titer boxplot across pseudotime bins
  program_vs_pseudotime_{sample}.pdf        : NMF program scores vs pseudotime (per program)
  program_pseudotime_corr_{sample}.pdf      : heatmap of Spearman rho (program x timepoint)
  pseudotime_umap_{sample}.pdf              : UMAP coloured by pseudotime
  sceptic_stats_{sample}.csv               : Spearman/Pearson correlations + Kruskal-Wallis
"""

import os
import argparse
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from scipy import stats
from scipy.stats import spearmanr, pearsonr, kruskal

try:
    from sceptic import run_sceptic_and_evaluate
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


def _dominant_program(metadata):
    """
    Return a Series of dominant NMF program labels (e.g. 'Program_3') per cell,
    computed as the argmax across all Program_* columns.
    """
    prog_cols = [c for c in metadata.columns if c.startswith("Program_")]
    if not prog_cols:
        return None
    return metadata[prog_cols].idxmax(axis=1).rename("dominant_program")


# ─────────────────────────────────────────────────────────────────────────────
# Load inputs from h5ad
# ─────────────────────────────────────────────────────────────────────────────

def load_from_h5ad(h5ad_path, pca_key="X_pca_harmony", timepoint_col="timepoint_numeric"):
    """
    Extract everything SCEPTIC needs directly from an h5ad object.

    Expects adata.obs to contain:
      - timepoint_numeric  : integer timepoint (0=uninfected, 999=persistent ctrl)
      - Program_0 … Program_N : NMF program usage scores
      - wolbachia_titer    : (optional) Wolbachia titer

    Returns
    -------
    data       : np.ndarray  (cells x PCs)
    labels     : np.ndarray  (numeric timepoint per cell)
    label_list : np.ndarray  (unique sorted timepoints)
    metadata   : pd.DataFrame  (includes Program_* columns)
    cell_ids   : list[str]
    """
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
    if timepoint_col in adata.obs.columns:
        labels = adata.obs[timepoint_col].values.astype(int)
        print(f"  Using timepoint column: {timepoint_col}")
    else:
        # Fallback: parse from bio_condition / timepoint string
        print(f"  WARNING: '{timepoint_col}' not found, parsing from bio_condition/timepoint")
        labels = adata.obs.apply(_parse_tp_fallback, axis=1).values.astype(int)

    label_list = np.array(sorted(np.unique(labels)))

    print(f"  Cells:        {adata.n_obs}")
    print(f"  Features:     {data.shape[1]}")
    print(f"  Timepoints:   {label_list}")
    print(f"  Label counts: {pd.Series(labels).value_counts().sort_index().to_dict()}")

    # ── Metadata ─────────────────────────────────────────────────────────────
    prog_cols = [c for c in adata.obs.columns if c.startswith("Program_")]
    base_cols = ["method", "timepoint", "timepoint_numeric", "bio_condition",
                 "cell_line", "treatment", "replicate", "wolbachia_titer",
                 "phase", "cyclum_theta"]
    keep_cols = [c for c in base_cols + prog_cols if c in adata.obs.columns]
    metadata = adata.obs[keep_cols].copy()

    print(f"  NMF programs found: {prog_cols}")

    return data, labels, label_list, metadata, adata.obs_names.tolist()


def _parse_tp_fallback(row):
    """Fallback timepoint parser for h5ads without timepoint_numeric."""
    bio = str(row.get("bio_condition", "")).strip()
    tp  = str(row.get("timepoint", "")).strip()
    if "wMel-Ctrl" in bio:
        return 999
    if "DOX-Ctrl" in bio:
        return 0
    s = tp.lstrip("Dd")
    try:
        return int(s)
    except ValueError:
        return 0


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

    rho, p = spearmanr(pseudotime, true_labels)
    ax.text(0.02, 0.97, f"Spearman rho={rho:.3f}  p={p:.2e}",
            transform=ax.transAxes, va="top", fontsize=9,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.7))
    plt.tight_layout()
    _savefig(fig, os.path.join(fig_dir, f"pseudotime_violin_{sample}.pdf"))
    return rho, p


def plot_pseudotime_by_program(pseudotime, metadata, fig_dir, sample):
    """Pseudotime distribution per dominant NMF program."""
    dom_prog = _dominant_program(metadata)
    if dom_prog is None:
        print("  Skipping program pseudotime plot: no Program_* columns in metadata")
        return

    df = pd.DataFrame({
        "pseudotime":       pseudotime,
        "dominant_program": dom_prog.values,
    })
    programs = sorted(df["dominant_program"].unique(),
                      key=lambda x: int(x.split("_")[1]))
    cmap = plt.cm.get_cmap("tab20")
    palette = {p: cmap(i % 20) for i, p in enumerate(programs)}

    fig, ax = plt.subplots(figsize=(max(10, len(programs) * 1.4), 5))
    sns.violinplot(data=df, x="dominant_program", y="pseudotime",
                   order=programs, palette=palette, inner="box", ax=ax)
    ax.set_xlabel("Dominant NMF Program")
    ax.set_ylabel("SCEPTIC Pseudotime")
    ax.set_title(f"Pseudotime by dominant NMF program — {sample}")
    plt.xticks(rotation=45, ha="right")

    groups = [df[df["dominant_program"] == p]["pseudotime"].values for p in programs]
    if len(groups) >= 2:
        h, p_kw = kruskal(*groups)
        ax.text(0.02, 0.97, f"Kruskal-Wallis H={h:.2f}  p={p_kw:.2e}",
                transform=ax.transAxes, va="top", fontsize=9,
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.7))
    plt.tight_layout()
    _savefig(fig, os.path.join(fig_dir, f"pseudotime_by_program_{sample}.pdf"))


def plot_program_vs_pseudotime(pseudotime, metadata, fig_dir, sample):
    """
    Two views of NMF program scores vs pseudotime:

    1. Per-program scatter (score vs pseudotime) with LOWESS smoother,
       one panel per program.
    2. Heatmap: Spearman rho between each program score and pseudotime,
       broken down by timepoint group.
    """
    prog_cols = [c for c in metadata.columns if c.startswith("Program_")]
    if not prog_cols:
        print("  Skipping program-vs-pseudotime plots: no Program_* columns")
        return

    from statsmodels.nonparametric.smoothers_lowess import lowess

    pt = pd.Series(pseudotime, index=metadata.index)

    # ── Plot 1: per-program scatter ──────────────────────────────────────────
    ncols = min(5, len(prog_cols))
    nrows = int(np.ceil(len(prog_cols) / ncols))
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(ncols * 3.5, nrows * 3), sharey=False)
    axes = np.array(axes).flatten()
    cmap = plt.cm.get_cmap("tab20")

    for i, prog in enumerate(prog_cols):
        ax = axes[i]
        scores = metadata[prog].values
        valid  = ~np.isnan(scores)
        x, y   = pt.values[valid], scores[valid]

        ax.scatter(x, y, alpha=0.2, s=4, color=cmap(i % 20), rasterized=True)

        if len(x) >= 10:
            smoothed = lowess(y, x, frac=0.3)
            ax.plot(smoothed[:, 0], smoothed[:, 1], color="black",
                    linewidth=1.5, label="LOWESS")
            rho, p_r = spearmanr(x, y)
            ax.text(0.05, 0.95, f"rho={rho:.2f}\np={p_r:.1e}",
                    transform=ax.transAxes, va="top", fontsize=7,
                    bbox=dict(boxstyle="round", facecolor="white", alpha=0.7))

        prog_num = prog.split("_")[1]
        ax.set_title(f"Program {prog_num}", fontsize=9)
        ax.set_xlabel("Pseudotime", fontsize=8)
        ax.set_ylabel("NMF score", fontsize=8)

    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    plt.suptitle(f"NMF program scores vs pseudotime — {sample}", fontweight="bold")
    plt.tight_layout()
    _savefig(fig, os.path.join(fig_dir, f"program_vs_pseudotime_{sample}.pdf"))

    # ── Plot 2: Spearman rho heatmap (programs × timepoint groups) ───────────
    tp_col = "timepoint_numeric" if "timepoint_numeric" in metadata.columns else "timepoint"
    if tp_col in metadata.columns:
        timepoints = sorted(metadata[tp_col].unique())
        rho_matrix = pd.DataFrame(index=prog_cols, columns=[str(t) for t in timepoints],
                                  dtype=float)
        for tp in timepoints:
            mask = metadata[tp_col] == tp
            for prog in prog_cols:
                x_tp = pt.values[mask]
                y_tp = metadata.loc[mask, prog].values
                valid = ~np.isnan(y_tp)
                if valid.sum() >= 5:
                    rho_matrix.loc[prog, str(tp)], _ = spearmanr(x_tp[valid], y_tp[valid])

        # Rename index for readability
        rho_matrix.index = [f"Prog {c.split('_')[1]}" for c in rho_matrix.index]

        fig, ax = plt.subplots(figsize=(max(6, len(timepoints) * 1.2),
                                        max(5, len(prog_cols) * 0.5)))
        sns.heatmap(rho_matrix.astype(float), cmap="RdBu_r", center=0,
                    vmin=-1, vmax=1, annot=True, fmt=".2f",
                    linewidths=0.5, ax=ax,
                    cbar_kws={"label": "Spearman rho (program ~ pseudotime)"})
        ax.set_xlabel("Timepoint")
        ax.set_ylabel("NMF Program")
        ax.set_title(f"Program–pseudotime correlation by timepoint — {sample}")
        plt.tight_layout()
        _savefig(fig, os.path.join(fig_dir, f"program_pseudotime_corr_{sample}.pdf"))


def plot_titer_vs_pseudotime(pseudotime, metadata, fig_dir, sample, n_bins=10):
    """
    Core analysis: Wolbachia titer as a function of SCEPTIC pseudotime.
    Produces scatter + regression, binned boxplot, timepoint-coloured scatter,
    and a per-dominant-program facet.
    """
    if "wolbachia_titer" not in metadata.columns:
        print("  Skipping titer plots: no wolbachia_titer column")
        return {}

    df = pd.DataFrame({
        "pseudotime":      pseudotime,
        "wolbachia_titer": metadata["wolbachia_titer"].values,
        "timepoint":       metadata["timepoint"].values.astype(str)
                           if "timepoint" in metadata.columns
                           else metadata.get("timepoint_numeric", pd.Series(0, index=metadata.index)).values.astype(str),
    }, index=metadata.index)

    # Add dominant program column
    dom_prog = _dominant_program(metadata)
    if dom_prog is not None:
        df["dominant_program"] = dom_prog.values

    df_titer = df.dropna(subset=["wolbachia_titer"])
    print(f"\n  Cells with titer for correlation: {len(df_titer)}")

    if len(df_titer) < 10:
        print("  WARNING: Too few cells with titer data for reliable analysis")
        return {}

    # ── Stats ────────────────────────────────────────────────────────────────
    rho_sp, p_sp = spearmanr(df_titer["pseudotime"], df_titer["wolbachia_titer"])
    rho_pe, p_pe = pearsonr( df_titer["pseudotime"], df_titer["wolbachia_titer"])
    slope, intercept, r_value, p_lin, _ = stats.linregress(
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
    from statsmodels.nonparametric.smoothers_lowess import lowess

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(df_titer["pseudotime"], df_titer["wolbachia_titer"],
               alpha=0.3, s=8, color="#2196F3", rasterized=True)
    x_line = np.linspace(df_titer["pseudotime"].min(), df_titer["pseudotime"].max(), 200)
    ax.plot(x_line, slope * x_line + intercept, color="red", linewidth=2,
            label=f"Linear fit (r²={r_value**2:.3f})")
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

    # ── Plot 3: coloured by timepoint ────────────────────────────────────────
    timepoints = sorted(df_titer["timepoint"].unique())
    tp_palette = dict(zip(timepoints, sns.color_palette("viridis", len(timepoints))))

    fig, ax = plt.subplots(figsize=(9, 6))
    for tp in timepoints:
        sub = df_titer[df_titer["timepoint"] == tp]
        label = "uninfected" if tp in ("0", "0.0") else f"D{tp}"
        ax.scatter(sub["pseudotime"], sub["wolbachia_titer"],
                   label=label, color=tp_palette[tp], alpha=0.5, s=8, rasterized=True)
    ax.set_xlabel("SCEPTIC Pseudotime")
    ax.set_ylabel("Wolbachia Titer")
    ax.set_title(f"Wolbachia titer vs pseudotime — coloured by timepoint — {sample}")
    ax.legend(title="Timepoint", bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)
    plt.tight_layout()
    _savefig(fig, os.path.join(fig_dir, f"titer_vs_pseudotime_by_timepoint_{sample}.pdf"))

    # ── Plot 4: faceted by dominant NMF program ───────────────────────────────
    if "dominant_program" in df_titer.columns:
        programs = sorted(df_titer["dominant_program"].unique(),
                          key=lambda x: int(x.split("_")[1]))
        ncols = min(4, len(programs))
        nrows = int(np.ceil(len(programs) / ncols))
        fig, axes = plt.subplots(nrows, ncols,
                                 figsize=(ncols * 4, nrows * 3.5), sharey=True)
        axes = np.array(axes).flatten()
        cmap = plt.cm.get_cmap("tab20")

        for i, prog in enumerate(programs):
            ax = axes[i]
            sub = df_titer[df_titer["dominant_program"] == prog]
            prog_num = int(prog.split("_")[1])
            ax.scatter(sub["pseudotime"], sub["wolbachia_titer"],
                       color=cmap(prog_num % 20), alpha=0.4, s=6, rasterized=True)
            if len(sub) >= 5:
                s, p_s = spearmanr(sub["pseudotime"], sub["wolbachia_titer"])
                sl, ic, _, _, _ = stats.linregress(sub["pseudotime"],
                                                    sub["wolbachia_titer"])
                xl = np.linspace(sub["pseudotime"].min(), sub["pseudotime"].max(), 100)
                ax.plot(xl, sl * xl + ic, color="red", linewidth=1.5)
                ax.text(0.05, 0.95, f"rho={s:.2f} p={p_s:.1e}",
                        transform=ax.transAxes, va="top", fontsize=7,
                        bbox=dict(boxstyle="round", facecolor="white", alpha=0.7))
            ax.set_title(f"Program {prog_num}", fontsize=9)
            ax.set_xlabel("Pseudotime", fontsize=8)
            if i % ncols == 0:
                ax.set_ylabel("Wolbachia Titer", fontsize=8)

        for j in range(i + 1, len(axes)):
            axes[j].set_visible(False)

        plt.suptitle(f"Titer vs pseudotime by dominant NMF program — {sample}",
                     fontweight="bold")
        plt.tight_layout()
        _savefig(fig, os.path.join(fig_dir,
                                   f"titer_vs_pseudotime_by_program_{sample}.pdf"))

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
        "cell_id":        cell_ids,
        "pseudotime":     pseudotime,
        "true_timepoint": labels,
        "pred_timepoint": label_predicted,
    })
    for i, tp in enumerate(label_list):
        results_df[f"prob_t{tp}"] = prob[:, i]

    # Add dominant program assignment
    dom_prog = _dominant_program(metadata)
    if dom_prog is not None:
        results_df["dominant_program"] = dom_prog.values

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
        description="Run SCEPTIC pseudotime and analyse Wolbachia titer with NMF programs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--h5ad",          required=True,
                        help="h5ad with NMF programs (adata_with_programs.h5ad)")
    parser.add_argument("--sample",        default="wolbachia_infection")
    parser.add_argument("--fig_dir",       default="results/sceptic")
    parser.add_argument("--pca_key",       default="X_pca_harmony",
                        help="obsm key for PCA embedding (default: X_pca_harmony)")
    parser.add_argument("--timepoint_col", default="timepoint_numeric",
                        help="obs column with numeric timepoints (default: timepoint_numeric)")
    parser.add_argument("--method",        default="xgboost",
                        choices=["xgboost", "svm"])
    parser.add_argument("--n_bins",        type=int, default=10)

    args = parser.parse_args()
    os.makedirs(args.fig_dir, exist_ok=True)

    # ── Load ─────────────────────────────────────────────────────────────────
    data, labels, label_list, metadata, cell_ids = load_from_h5ad(
        args.h5ad,
        pca_key=args.pca_key,
        timepoint_col=args.timepoint_col,
    )

    # ── Run SCEPTIC ──────────────────────────────────────────────────────────
    cm_result, label_predicted, pseudotime, prob = run_sceptic(
        data, labels, label_list, method=args.method
    )

    # ── Plots ─────────────────────────────────────────────────────────────────
    print("\nGenerating plots …")
    plot_confusion_matrix(cm_result, label_list, args.fig_dir, args.sample)

    rho_sp, p_sp = plot_pseudotime_violin(
        pseudotime, labels, label_list, args.fig_dir, args.sample)

    plot_pseudotime_by_program(pseudotime, metadata, args.fig_dir, args.sample)

    plot_program_vs_pseudotime(pseudotime, metadata, args.fig_dir, args.sample)

    stat_results = plot_titer_vs_pseudotime(
        pseudotime, metadata, args.fig_dir, args.sample, n_bins=args.n_bins)

    stat_results["pseudotime_timepoint_spearman_rho"] = rho_sp
    stat_results["pseudotime_timepoint_spearman_p"]   = p_sp

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
