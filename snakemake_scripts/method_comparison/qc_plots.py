"""
qc_plots.py
===========
QC and exploratory plots for the integrated scRNA-seq dataset.
Run after integrate.py using the _reference.h5ad, _query.h5ad, and _combined.h5ad outputs.

Generates
---------
Cluster composition:
  - cells_per_cluster_{sample}.pdf
  - cluster_composition_by_method_{sample}.pdf       (% cells per cluster, 10x vs pipseq)
  - cluster_composition_by_condition_{sample}.pdf    (% cells per cluster, by bio_condition)

Gene counts:
  - genes_per_cluster_bar_{sample}.pdf               (mean ± std)
  - genes_per_cluster_violin_{sample}.pdf
  - genes_per_cluster_by_method_{sample}.pdf         (10x vs pipseq per cluster)

Wolbachia titer:
  - titer_by_cluster_{sample}.pdf                    (boxplot + strip)
  - titer_by_method_{sample}.pdf                     (wMel and DOX separately)
  - titer_by_timepoint_cluster_{sample}.pdf          (query cells only)
  - infection_pct_by_cluster_{sample}.pdf            (% infected cells per cluster)

Method comparison (chi-square):
  - method_chi2_results_{sample}.txt

UMAPs:
  - umap_leiden_{sample}.pdf
  - umap_method_{sample}.pdf
  - umap_bio_condition_{sample}.pdf
  - umap_titer_{sample}.pdf
  - umap_split_by_method_{sample}.pdf                (side-by-side 10x / pipseq)

python qc_plots.py \
    --combined results/integrated_combined.h5ad \
    --ref      results/integrated_reference.h5ad \
    --query    results/integrated_query.h5ad \
    --sample   wolbachia_infection \
    --fig_dir  results/figures/qc
"""

import os
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from scipy.stats import chi2_contingency, mannwhitneyu


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _leiden_colors(adata, key="leiden_ref"):
    """Return one tab20 colour per cluster."""
    clusters = sorted(adata.obs[key].unique())
    cmap = plt.cm.get_cmap("tab20")
    return [cmap(i % 20) for i in range(len(clusters))]


def _cluster_col(adata):
    """Return whichever leiden column is present."""
    for col in ("leiden_ref", "leiden"):
        if col in adata.obs.columns:
            return col
    raise KeyError("No leiden column found in adata.obs")


def _savefig(fig, path):
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {path}")


# ─────────────────────────────────────────────────────────────────────────────
# Individual plot functions
# ─────────────────────────────────────────────────────────────────────────────

def plot_cells_per_cluster(adata, fig_dir, sample):
    col = _cluster_col(adata)
    colors = _leiden_colors(adata, key=col)
    counts = adata.obs[col].value_counts().sort_index()

    fig, ax = plt.subplots(figsize=(max(8, len(counts) * 0.6), 5))
    ax.bar(range(len(counts)), counts.values, color=colors, alpha=0.85, edgecolor="black", lw=0.4)
    ax.set_xticks(range(len(counts)))
    ax.set_xticklabels(counts.index)
    ax.set_xlabel("Leiden Cluster")
    ax.set_ylabel("Number of cells")
    ax.set_title(f"Cells per cluster — {sample}")
    ax.grid(axis="y", alpha=0.3)
    for i, v in enumerate(counts.values):
        ax.text(i, v + counts.max() * 0.01, str(v), ha="center", va="bottom", fontsize=7)
    plt.tight_layout()
    _savefig(fig, os.path.join(fig_dir, f"cells_per_cluster_{sample}.pdf"))
    return counts


def plot_genes_per_cluster(adata, fig_dir, sample):
    col = _cluster_col(adata)
    colors = _leiden_colors(adata, key=col)
    stats = adata.obs.groupby(col)["n_genes"].agg(["mean", "std"])
    clusters = stats.index.tolist()

    # Bar with error bars
    fig, ax = plt.subplots(figsize=(max(8, len(clusters) * 0.6), 5))
    ax.bar(range(len(stats)), stats["mean"], yerr=stats["std"],
           color=colors, alpha=0.85, capsize=4, edgecolor="black", lw=0.4)
    ax.set_xticks(range(len(stats)))
    ax.set_xticklabels(clusters)
    ax.set_xlabel("Leiden Cluster")
    ax.set_ylabel("Mean genes per cell")
    ax.set_title(f"Genes per cluster (mean ± SD) — {sample}")
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    _savefig(fig, os.path.join(fig_dir, f"genes_per_cluster_bar_{sample}.pdf"))

    # Violin
    fig, ax = plt.subplots(figsize=(max(10, len(clusters) * 0.7), 5))
    sc.pl.violin(adata, "n_genes", groupby=col, ax=ax, show=False, rotation=0)
    ax.set_title(f"Gene count distribution by cluster — {sample}")
    plt.tight_layout()
    _savefig(fig, os.path.join(fig_dir, f"genes_per_cluster_violin_{sample}.pdf"))

    return stats


def plot_genes_per_cluster_by_method(adata, fig_dir, sample):
    """Mean genes per cluster split by library prep method."""
    col = _cluster_col(adata)
    if "method" not in adata.obs.columns:
        print("  Skipping genes-by-method: no 'method' column")
        return

    pivot = adata.obs.groupby([col, "method"])["n_genes"].mean().unstack()
    fig, ax = plt.subplots(figsize=(max(10, len(pivot) * 0.7), 5))
    pivot.plot(kind="bar", ax=ax, color=["#1f77b4", "#ff7f0e"], edgecolor="black", lw=0.4)
    ax.set_xlabel("Leiden Cluster")
    ax.set_ylabel("Mean genes per cell")
    ax.set_title(f"Genes per cluster by library prep method — {sample}")
    ax.legend(title="Method")
    plt.xticks(rotation=0)
    plt.tight_layout()
    _savefig(fig, os.path.join(fig_dir, f"genes_per_cluster_by_method_{sample}.pdf"))


def plot_cluster_composition_by_method(adata, fig_dir, sample):
    """% cells per cluster for each library prep method + chi-square test."""
    col = _cluster_col(adata)
    if "method" not in adata.obs.columns:
        print("  Skipping cluster-by-method: no 'method' column")
        return

    pct = pd.crosstab(adata.obs[col], adata.obs["method"], normalize="columns") * 100
    fig, ax = plt.subplots(figsize=(max(10, len(pct) * 0.7), 5))
    pct.plot(kind="bar", ax=ax, color=["#1f77b4", "#ff7f0e"], edgecolor="black", lw=0.4)
    ax.set_xlabel("Leiden Cluster")
    ax.set_ylabel("% of cells")
    ax.set_title(f"Cluster composition by library prep method — {sample}")
    ax.legend(title="Method")
    plt.xticks(rotation=0)
    plt.tight_layout()
    _savefig(fig, os.path.join(fig_dir, f"cluster_composition_by_method_{sample}.pdf"))

    # Chi-square test
    contingency = pd.crosstab(adata.obs[col], adata.obs["method"])
    chi2, p, dof, _ = chi2_contingency(contingency)
    result = (
        f"Chi-square test — cluster distribution by method\n"
        f"chi2={chi2:.4f}  p={p:.2e}  dof={dof}\n"
        f"Methods show {'SIGNIFICANT' if p < 0.05 else 'NO'} difference in cluster distribution\n"
        f"\nContingency table:\n{contingency.to_string()}\n"
        f"\nColumn-normalised (%):\n{pct.to_string()}\n"
    )
    txt_path = os.path.join(fig_dir, f"method_chi2_results_{sample}.txt")
    with open(txt_path, "w") as fh:
        fh.write(result)
    print(f"  Saved: {txt_path}")
    print(f"  Chi-square: chi2={chi2:.4f}  p={p:.2e}  "
          f"({'SIGNIFICANT' if p < 0.05 else 'NOT significant'})")


def plot_cluster_composition_by_condition(adata, fig_dir, sample):
    """% cells per cluster for each biological condition."""
    col = _cluster_col(adata)
    if "bio_condition" not in adata.obs.columns:
        print("  Skipping cluster-by-condition: no 'bio_condition' column")
        return

    colors = _leiden_colors(adata, key=col)
    pct = pd.crosstab(adata.obs[col], adata.obs["bio_condition"], normalize="columns") * 100

    fig, ax = plt.subplots(figsize=(max(10, len(pct.columns) * 1.2), 6))
    pct.T.plot(kind="bar", stacked=True, ax=ax, color=colors, width=0.8, edgecolor="black", lw=0.2)
    ax.set_xlabel("Biological Condition")
    ax.set_ylabel("% of cells")
    ax.set_title(f"Cluster composition by biological condition — {sample}")
    ax.legend(title="Cluster", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    _savefig(fig, os.path.join(fig_dir, f"cluster_composition_by_condition_{sample}.pdf"))


def plot_titer_by_cluster(adata, fig_dir, sample):
    """Boxplot + strip of Wolbachia titer per cluster."""
    if "wolbachia_titer" not in adata.obs.columns:
        print("  Skipping titer-by-cluster: no 'wolbachia_titer' column")
        return

    col = _cluster_col(adata)
    colors = _leiden_colors(adata, key=col)
    clusters = sorted(adata.obs[col].unique())
    plot_data = adata.obs[[col, "wolbachia_titer"]].copy().sort_values(col)

    fig, ax = plt.subplots(figsize=(max(10, len(clusters) * 0.7), 6))
    sns.stripplot(data=plot_data, x=col, y="wolbachia_titer",
                  color="black", alpha=0.2, size=1.5, ax=ax, order=clusters)
    bp = ax.boxplot(
        [plot_data[plot_data[col] == c]["wolbachia_titer"].values for c in clusters],
        positions=range(len(clusters)), widths=0.55, patch_artist=True,
        whiskerprops=dict(alpha=0.7), capprops=dict(alpha=0.7),
        medianprops=dict(color="black", linewidth=2), flierprops=dict(markersize=1),
    )
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color); patch.set_alpha(0.75)
    ax.set_xlabel("Leiden Cluster")
    ax.set_ylabel("Wolbachia Titer")
    ax.set_xticklabels(clusters)
    ax.set_title(f"Wolbachia titer by cluster — {sample}")
    plt.tight_layout()
    _savefig(fig, os.path.join(fig_dir, f"titer_by_cluster_{sample}.pdf"))

    # Summary stats
    stats = adata.obs.groupby(col)["wolbachia_titer"].agg(["mean", "median", "std", "count"])
    print(stats.to_string())


def plot_infection_pct_by_cluster(adata, fig_dir, sample):
    """% Wolbachia-positive cells per cluster."""
    if "wolbachia_titer" not in adata.obs.columns:
        print("  Skipping infection-pct: no 'wolbachia_titer' column")
        return

    col = _cluster_col(adata)
    colors = _leiden_colors(adata, key=col)

    pct_infected = adata.obs.groupby(col).apply(
        lambda x: (x["wolbachia_titer"] > 0).sum() / len(x) * 100
    )

    fig, ax = plt.subplots(figsize=(max(8, len(pct_infected) * 0.6), 5))
    ax.bar(range(len(pct_infected)), pct_infected.values,
           color=colors, alpha=0.85, edgecolor="black", lw=0.4)
    ax.set_xticks(range(len(pct_infected)))
    ax.set_xticklabels(pct_infected.index)
    ax.set_xlabel("Leiden Cluster")
    ax.set_ylabel("% Wolbachia-positive cells")
    ax.set_title(f"Infection rate by cluster — {sample}")
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    _savefig(fig, os.path.join(fig_dir, f"infection_pct_by_cluster_{sample}.pdf"))


def plot_titer_by_method(adata, fig_dir, sample):
    """Titer comparison 10x vs PIPseq, split by cell line."""
    if "wolbachia_titer" not in adata.obs.columns or "method" not in adata.obs.columns:
        print("  Skipping titer-by-method: missing column(s)")
        return

    for cell_line, label in [("JW18wMel", "wMel"), ("JW18DOX", "DOX")]:
        if "cell_line" not in adata.obs.columns:
            break
        subset = adata[adata.obs["cell_line"] == cell_line]
        if subset.n_obs == 0:
            continue

        fig, ax = plt.subplots(figsize=(7, 5))
        sns.boxplot(data=subset.obs, x="method", y="wolbachia_titer",
                    ax=ax, palette=["#1f77b4", "#ff7f0e"], flierprops=dict(markersize=1))
        sns.stripplot(data=subset.obs, x="method", y="wolbachia_titer",
                      ax=ax, color="black", alpha=0.3, size=2)
        ax.set_xlabel("Library Prep Method")
        ax.set_ylabel("Wolbachia Titer")
        ax.set_title(f"{label} Ctrl — titer by method — {sample}")

        t10x = subset.obs[subset.obs["method"] == "10x"]["wolbachia_titer"].dropna()
        tpip = subset.obs[subset.obs["method"] == "pipseq"]["wolbachia_titer"].dropna()
        if len(t10x) > 0 and len(tpip) > 0:
            _, p = mannwhitneyu(t10x, tpip, alternative="two-sided")
            ax.text(0.5, 0.97, f"Mann-Whitney U p = {p:.2e}",
                    transform=ax.transAxes, ha="center", va="top", fontsize=9)

        plt.tight_layout()
        _savefig(fig, os.path.join(fig_dir, f"titer_by_method_{label}_{sample}.pdf"))


def plot_titer_by_timepoint(adata, fig_dir, sample):
    """Titer across infection timepoints, coloured by cluster (query cells only)."""
    if "wolbachia_titer" not in adata.obs.columns:
        print("  Skipping titer-by-timepoint: no 'wolbachia_titer' column")
        return
    if "timepoint_numeric" not in adata.obs.columns:
        print("  Skipping titer-by-timepoint: no 'timepoint_numeric' column")
        return

    col = _cluster_col(adata)
    fig, ax = plt.subplots(figsize=(12, 5))
    sns.boxplot(data=adata.obs, x="timepoint_numeric", y="wolbachia_titer",
                hue=col, ax=ax, flierprops=dict(markersize=1))
    ax.set_xlabel("Timepoint (days post-infection; 0 = uninfected)")
    ax.set_ylabel("Wolbachia Titer")
    ax.set_title(f"Wolbachia titer by timepoint and cluster — {sample}")
    ax.legend(title="Cluster", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    _savefig(fig, os.path.join(fig_dir, f"titer_by_timepoint_cluster_{sample}.pdf"))


def plot_umaps(adata, fig_dir, sample):
    """Core UMAP panels coloured by cluster, method, condition, titer."""
    col = _cluster_col(adata)

    sc.pl.umap(adata, color=col, legend_loc="on data",
               save=f"_leiden_{sample}.pdf",
               title=f"Leiden clusters — {sample}")

    if "method" in adata.obs.columns:
        sc.pl.umap(adata, color="method",
                   save=f"_method_{sample}.pdf",
                   title=f"Library prep method — {sample}")

        # Side-by-side split by method
        methods = adata.obs["method"].dropna().unique()
        if len(methods) == 2:
            fig, axes = plt.subplots(1, 2, figsize=(14, 5))
            for ax, meth in zip(axes, sorted(methods)):
                sub = adata[adata.obs["method"] == meth]
                sc.pl.umap(sub, color=col, ax=ax, show=False,
                           title=f"{meth}", frameon=False)
            plt.suptitle(f"Clusters split by method — {sample}", fontweight="bold")
            plt.tight_layout()
            _savefig(fig, os.path.join(fig_dir, f"umap_split_by_method_{sample}.pdf"))

    if "bio_condition" in adata.obs.columns:
        sc.pl.umap(adata, color="bio_condition",
                   save=f"_bio_condition_{sample}.pdf",
                   title=f"Biological condition — {sample}")

    if "wolbachia_titer" in adata.obs.columns:
        vmax = float(adata.obs["wolbachia_titer"].quantile(0.95))
        sc.pl.umap(adata, color="wolbachia_titer", vmax=vmax,
                   save=f"_titer_{sample}.pdf",
                   title=f"Wolbachia titer — {sample}")

        if "method" in adata.obs.columns:
            methods = adata.obs["method"].dropna().unique()
            if len(methods) == 2:
                fig, axes = plt.subplots(1, 2, figsize=(14, 5))
                for ax, meth in zip(axes, sorted(methods)):
                    sub = adata[adata.obs["method"] == meth]
                    sc.pl.umap(sub, color="wolbachia_titer", vmax=vmax,
                               ax=ax, show=False, title=f"{meth}", frameon=False)
                plt.suptitle(f"Wolbachia titer split by method — {sample}", fontweight="bold")
                plt.tight_layout()
                _savefig(fig, os.path.join(fig_dir, f"umap_titer_split_method_{sample}.pdf"))

    if "timepoint_numeric" in adata.obs.columns:
        sc.pl.umap(adata, color="timepoint_numeric",
                   save=f"_timepoint_{sample}.pdf",
                   title=f"Timepoint — {sample}")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def run_qc_plots(combined_path, ref_path, query_path, fig_dir, sample):
    os.makedirs(fig_dir, exist_ok=True)
    sc.settings.figdir = fig_dir

    print(f"\nLoading combined object: {combined_path}")
    combined = sc.read_h5ad(combined_path)

    print(f"Loading reference object: {ref_path}")
    ref = sc.read_h5ad(ref_path)

    print(f"Loading query object: {query_path}")
    query = sc.read_h5ad(query_path)

    print(f"\nCombined: {combined.n_obs} cells, {combined.n_vars} genes")
    print(f"Reference: {ref.n_obs} cells")
    print(f"Query:     {query.n_obs} cells")

    # ── Combined dataset plots ────────────────────────────────────────────────
    print("\n--- Combined dataset ---")
    plot_umaps(combined, fig_dir, f"{sample}_combined")
    plot_cells_per_cluster(combined, fig_dir, f"{sample}_combined")
    plot_genes_per_cluster(combined, fig_dir, f"{sample}_combined")
    plot_genes_per_cluster_by_method(combined, fig_dir, f"{sample}_combined")
    plot_cluster_composition_by_method(combined, fig_dir, f"{sample}_combined")
    plot_cluster_composition_by_condition(combined, fig_dir, f"{sample}_combined")
    plot_titer_by_cluster(combined, fig_dir, f"{sample}_combined")
    plot_infection_pct_by_cluster(combined, fig_dir, f"{sample}_combined")
    plot_titer_by_method(combined, fig_dir, f"{sample}_combined")
    plot_titer_by_timepoint(combined, fig_dir, f"{sample}_combined")

    # ── Reference-only plots ─────────────────────────────────────────────────
    print("\n--- Reference (uninfected controls) ---")
    plot_umaps(ref, fig_dir, f"{sample}_ref")
    plot_cells_per_cluster(ref, fig_dir, f"{sample}_ref")
    plot_genes_per_cluster(ref, fig_dir, f"{sample}_ref")
    plot_genes_per_cluster_by_method(ref, fig_dir, f"{sample}_ref")
    plot_cluster_composition_by_method(ref, fig_dir, f"{sample}_ref")
    plot_titer_by_cluster(ref, fig_dir, f"{sample}_ref")
    plot_titer_by_method(ref, fig_dir, f"{sample}_ref")

    # ── Query-only plots ─────────────────────────────────────────────────────
    print("\n--- Query (new infection timepoints) ---")
    plot_umaps(query, fig_dir, f"{sample}_query")
    plot_cells_per_cluster(query, fig_dir, f"{sample}_query")
    plot_genes_per_cluster(query, fig_dir, f"{sample}_query")
    plot_genes_per_cluster_by_method(query, fig_dir, f"{sample}_query")
    plot_cluster_composition_by_method(query, fig_dir, f"{sample}_query")
    plot_titer_by_cluster(query, fig_dir, f"{sample}_query")
    plot_titer_by_timepoint(query, fig_dir, f"{sample}_query")

    print(f"\nAll plots saved to {fig_dir}/")


def main():
    parser = argparse.ArgumentParser(
        description="QC plots for integrated scRNA-seq data. Run after integrate.py."
    )
    parser.add_argument("--combined", required=True,
                        help="Path to *_combined.h5ad from integrate.py")
    parser.add_argument("--ref", required=True,
                        help="Path to *_reference.h5ad from integrate.py")
    parser.add_argument("--query", required=True,
                        help="Path to *_query.h5ad from integrate.py")
    parser.add_argument("--sample", default="wolbachia_infection",
                        help="Label used in output filenames")
    parser.add_argument("--fig_dir", default="figures/qc",
                        help="Directory to save figures")

    args = parser.parse_args()

    run_qc_plots(
        combined_path=args.combined,
        ref_path=args.ref,
        query_path=args.query,
        fig_dir=args.fig_dir,
        sample=args.sample,
    )


if __name__ == "__main__":
    main()