#!/usr/bin/env python
"""
Randomization (permutation) test for association between cell cycle stage
(cc_stage) and Leiden cluster assignment.

For each Leiden cluster, computes the observed chi-squared contribution
(sum of squared standardized residuals across cc_stage categories), then
builds a null distribution by permuting cc_stage labels across all cells
many times (preserving cluster sizes and overall cc_stage proportions).
Empirical p-values are the fraction of permutations whose per-cluster
statistic meets or exceeds the observed statistic, BH-corrected across
clusters.

Also reports the global chi-squared statistic, its permutation p-value,
Cramer's V effect size, and per-cluster/per-stage standardized residuals
(to show which stage is over/under-represented in each cluster).

Usage:
    python cc_stage_leiden_permutation_test.py <path_to_h5ad> \
        [--cluster-col leiden] [--cc-col cc_stage] \
        [--n-perm 10000] [--seed 0] [--outdir cc_stage_leiden_test]
"""

import argparse
import os

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import chi2_contingency
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns


def benjamini_hochberg(pvals):
    """Return BH-adjusted p-values (FDR)."""
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order] * n / (np.arange(n) + 1)
    # enforce monotonicity from the largest p-value downward
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty(n)
    adjusted[order] = np.clip(ranked, 0, 1)
    return adjusted


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("h5ad_path", help="Path to input .h5ad file")
    parser.add_argument("--cluster-col", default="leiden", help="obs column with cluster labels (default: leiden)")
    parser.add_argument("--cc-col", default="cc_stage", help="obs column with cell cycle stage (default: cc_stage)")
    parser.add_argument("--n-perm", type=int, default=10000, help="number of permutations (default: 10000)")
    parser.add_argument("--seed", type=int, default=0, help="random seed (default: 0)")
    parser.add_argument("--outdir", default="cc_stage_leiden_test", help="output directory")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    print(f"Loading {args.h5ad_path}")
    adata = sc.read_h5ad(args.h5ad_path)

    for col in (args.cluster_col, args.cc_col):
        if col not in adata.obs.columns:
            raise ValueError(f"'{col}' not found in adata.obs")

    obs = adata.obs[[args.cluster_col, args.cc_col]].dropna().copy()
    n_dropped = adata.n_obs - len(obs)
    if n_dropped:
        print(f"Dropping {n_dropped} cells with NaN in '{args.cluster_col}' or '{args.cc_col}'")

    obs[args.cluster_col] = obs[args.cluster_col].astype("category")
    obs[args.cc_col] = obs[args.cc_col].astype("category")

    clusters = obs[args.cluster_col].cat.categories.tolist()
    stages = obs[args.cc_col].cat.categories.tolist()
    n_clusters = len(clusters)
    n_stages = len(stages)
    n = len(obs)

    cluster_codes = obs[args.cluster_col].cat.codes.to_numpy()
    stage_codes = obs[args.cc_col].cat.codes.to_numpy()

    # --- observed contingency table ---
    observed_table = pd.crosstab(obs[args.cluster_col], obs[args.cc_col])
    observed_table = observed_table.reindex(index=clusters, columns=stages)
    observed = observed_table.to_numpy(dtype=float)

    row_totals = observed.sum(axis=1)          # cluster sizes (fixed under permutation)
    col_totals = observed.sum(axis=0)          # global stage totals (fixed under permutation)
    expected = np.outer(row_totals, col_totals) / n

    # --- global chi-squared (asymptotic, for reference) ---
    chi2_stat, chi2_p_asymptotic, dof, _ = chi2_contingency(observed)
    cramers_v = np.sqrt(chi2_stat / (n * (min(n_clusters, n_stages) - 1)))

    # --- observed per-cluster chi-squared contribution ---
    def row_stats(table):
        return ((table - expected) ** 2 / expected).sum(axis=1)

    observed_row_stat = row_stats(observed)
    observed_global_stat = observed_row_stat.sum()  # equals chi2_stat

    # --- permutation null distribution ---
    print(f"Running {args.n_perm} permutations over {n} cells...")
    null_row_stats = np.zeros((args.n_perm, n_clusters))
    null_global_stats = np.zeros(args.n_perm)

    for p in range(args.n_perm):
        permuted_stage_codes = rng.permutation(stage_codes)
        idx = cluster_codes * n_stages + permuted_stage_codes
        counts = np.bincount(idx, minlength=n_clusters * n_stages).reshape(n_clusters, n_stages).astype(float)
        rs = row_stats(counts)
        null_row_stats[p] = rs
        null_global_stats[p] = rs.sum()
        if (p + 1) % 1000 == 0:
            print(f"  {p + 1}/{args.n_perm} permutations done")

    # --- empirical p-values ---
    global_perm_p = (1 + np.sum(null_global_stats >= observed_global_stat)) / (1 + args.n_perm)

    cluster_perm_p = np.array([
        (1 + np.sum(null_row_stats[:, i] >= observed_row_stat[i])) / (1 + args.n_perm)
        for i in range(n_clusters)
    ])
    cluster_perm_p_fdr = benjamini_hochberg(cluster_perm_p)

    # --- standardized residuals (which stage drives each cluster's signal) ---
    residuals = (observed - expected) / np.sqrt(expected)
    residuals_df = pd.DataFrame(residuals, index=clusters, columns=stages)

    proportions = observed_table.div(observed_table.sum(axis=1), axis=0) * 100

    # --- summary table ---
    summary = pd.DataFrame({
        "cluster": clusters,
        "n_cells": row_totals.astype(int),
        "chi2_contribution": observed_row_stat,
        "perm_pval": cluster_perm_p,
        "perm_pval_fdr": cluster_perm_p_fdr,
    })
    for stage in stages:
        summary[f"pct_{stage}"] = proportions[stage].values
        summary[f"residual_{stage}"] = residuals_df[stage].values

    summary = summary.sort_values("chi2_contribution", ascending=False)

    # --- save outputs ---
    summary_path = os.path.join(args.outdir, "cc_stage_leiden_permutation_summary.csv")
    summary.to_csv(summary_path, index=False)

    global_summary_path = os.path.join(args.outdir, "cc_stage_leiden_global_summary.txt")
    with open(global_summary_path, "w") as f:
        f.write(f"Cell cycle stage ({args.cc_col}) vs cluster ({args.cluster_col}) association\n")
        f.write(f"n cells analyzed: {n}\n")
        f.write(f"n clusters: {n_clusters}, n stages: {n_stages}\n\n")
        f.write(f"Global chi2 statistic: {chi2_stat:.2f} (dof={dof})\n")
        f.write(f"Asymptotic chi2 p-value: {chi2_p_asymptotic:.3e}\n")
        f.write(f"Permutation p-value ({args.n_perm} perms): {global_perm_p:.3e}\n")
        f.write(f"Cramer's V: {cramers_v:.4f}\n")

    print(open(global_summary_path).read())

    # --- residual heatmap ---
    fig, ax = plt.subplots(figsize=(max(6, n_stages * 1.2), max(6, n_clusters * 0.4)))
    sns.heatmap(
        residuals_df, cmap="RdBu_r", center=0, annot=True, fmt=".1f",
        cbar_kws={"label": "Standardized residual"}, ax=ax,
    )
    ax.set_title(f"{args.cc_col} vs {args.cluster_col}: standardized residuals")
    ax.set_xlabel(args.cc_col)
    ax.set_ylabel(args.cluster_col)
    fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, "residual_heatmap.png"), dpi=300)
    fig.savefig(os.path.join(args.outdir, "residual_heatmap.pdf"))
    plt.close(fig)

    print(f"\nSaved:\n  {summary_path}\n  {global_summary_path}\n  {os.path.join(args.outdir, 'residual_heatmap.png')}")


if __name__ == "__main__":
    main()
