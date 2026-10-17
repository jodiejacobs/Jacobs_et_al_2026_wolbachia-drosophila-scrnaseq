#!/usr/bin/env python3
"""
compare_velocity_pseudotime.py

Tests whether independently computed RNA velocity agrees with SCEPTIC
pseudotime. This is the only step where the two results meet; neither
input was computed using the other.

Tests
  1. Rank agreement: Spearman rho between sceptic_pseudotime and the
     velocity-derived orderings (latent_time, velocity_pseudotime),
     overall and per replicate.
  2. Directionality: for each cell, the expected SCEPTIC pseudotime of its
     velocity-predicted successors minus its own pseudotime,
         delta_i = sum_j T_ij * pt_j  -  pt_i
     where T is the scVelo forward transition matrix. delta > 0 means
     velocity points toward later SCEPTIC pseudotime. Significance comes
     from a permutation null that shuffles pseudotime across cells.

Usage:
    python compare_velocity_pseudotime.py \
        --velocity results/velocity/velocity_pipseq.h5ad \
        --sceptic  results/sceptic/wolbachia_infection/sceptic_wolbachia_infection.h5ad \
        --outdir   results/velocity/pseudotime_comparison
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats


def get_transition_matrix(vdata, backward=False):
    try:
        from scvelo.tools.transition_matrix import transition_matrix
    except ImportError:
        from scvelo.utils import get_transition_matrix as transition_matrix
    return transition_matrix(vdata, vkey="velocity", backward=backward).tocsr()


def expected_delta(T, pt, has_pt):
    """delta_i = E_T[pt_next] - pt_i, renormalised over successors with a pseudotime."""
    Tm = T[:, has_pt]
    w = np.asarray(Tm.sum(axis=1)).ravel()
    nxt = np.asarray(Tm @ pt[has_pt]).ravel()
    with np.errstate(invalid="ignore", divide="ignore"):
        nxt = nxt / w
    delta = nxt - pt
    delta[~has_pt | (w == 0)] = np.nan
    return delta


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--velocity", required=True)
    p.add_argument("--sceptic", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--pt-key", default="sceptic_pseudotime")
    p.add_argument("--cluster-key", default="leiden")
    p.add_argument("--n-perm", type=int, default=1000)
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    vdata = sc.read_h5ad(args.velocity)
    sobs = sc.read_h5ad(args.sceptic, backed="r").obs[[args.pt_key]].copy()

    pt = sobs[args.pt_key].reindex(vdata.obs_names).to_numpy(dtype=float)
    has_pt = ~np.isnan(pt)
    print(f"Velocity cells: {vdata.n_obs}; with {args.pt_key}: {has_pt.sum()}")
    if has_pt.sum() == 0:
        raise ValueError("No overlapping cell IDs between velocity and SCEPTIC objects")

    obs = vdata.obs.copy()
    obs[args.pt_key] = pt

    # ── 1. Rank agreement ────────────────────────────────────────────────────
    rows = []
    for vkey in ("latent_time", "velocity_pseudotime"):
        sub = obs.loc[has_pt]
        rho, pval = stats.spearmanr(sub[args.pt_key], sub[vkey])
        rows.append({"group": "all", "velocity_ordering": vkey, "n_cells": len(sub),
                     "spearman_rho": rho, "p_value": pval})
        if "replicate" in sub:
            for rep, g in sub.groupby("replicate", observed=True):
                r, pv = stats.spearmanr(g[args.pt_key], g[vkey])
                rows.append({"group": f"replicate_{rep}", "velocity_ordering": vkey,
                             "n_cells": len(g), "spearman_rho": r, "p_value": pv})
    rank_df = pd.DataFrame(rows)
    rank_df.to_csv(os.path.join(args.outdir, "rank_agreement.csv"), index=False)
    print(rank_df.to_string(index=False))

    # ── 2. Directionality ────────────────────────────────────────────────────
    T = get_transition_matrix(vdata)
    delta = expected_delta(T, pt, has_pt)
    obs["delta_pseudotime"] = delta
    obs_mean = np.nanmean(delta)
    frac_pos = np.nanmean(delta > 0)
    wil = stats.wilcoxon(delta[~np.isnan(delta)])

    null = np.empty(args.n_perm)
    idx = np.where(has_pt)[0]
    for i in range(args.n_perm):
        pt_perm = pt.copy()
        pt_perm[idx] = rng.permutation(pt[idx])
        null[i] = np.nanmean(expected_delta(T, pt_perm, has_pt))
    p_perm = (1 + np.sum(null >= obs_mean)) / (1 + args.n_perm)

    summary = pd.DataFrame([{
        "n_cells": int((~np.isnan(delta)).sum()),
        "mean_delta_pseudotime": obs_mean,
        "median_delta_pseudotime": np.nanmedian(delta),
        "fraction_cells_forward": frac_pos,
        "wilcoxon_stat": wil.statistic, "wilcoxon_p": wil.pvalue,
        "perm_null_mean": null.mean(), "perm_null_sd": null.std(),
        "perm_p_one_sided": p_perm, "n_perm": args.n_perm,
    }])
    summary.to_csv(os.path.join(args.outdir, "directionality_summary.csv"), index=False)
    print(summary.T.to_string(header=False))

    if args.cluster_key in obs:
        per_cl = (obs.loc[has_pt]
                  .groupby(args.cluster_key, observed=True)
                  .agg(n_cells=("delta_pseudotime", "size"),
                       median_pseudotime=(args.pt_key, "median"),
                       median_latent_time=("latent_time", "median"),
                       median_delta=("delta_pseudotime", "median"),
                       frac_forward=("delta_pseudotime", lambda x: np.nanmean(x > 0)))
                  .sort_values("median_pseudotime"))
        per_cl.to_csv(os.path.join(args.outdir, "directionality_by_cluster.csv"))
        print(per_cl.to_string())

    obs[[args.pt_key, "latent_time", "velocity_pseudotime", "velocity_confidence",
         "delta_pseudotime"]].to_csv(os.path.join(args.outdir, "per_cell.csv"))

    # ── Figures ──────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    sub = obs.loc[has_pt]
    for ax, vkey in zip(axes[:2], ("latent_time", "velocity_pseudotime")):
        ax.hexbin(sub[args.pt_key], sub[vkey], gridsize=40, mincnt=1, cmap="viridis")
        rho = rank_df.query("group == 'all' and velocity_ordering == @vkey")["spearman_rho"].iloc[0]
        ax.set(xlabel="SCEPTIC pseudotime", ylabel=vkey, title=f"Spearman rho = {rho:.2f}")
    axes[2].hist(null, bins=40, color="grey", label="permutation null")
    axes[2].axvline(obs_mean, color="crimson", label="observed")
    axes[2].set(xlabel="mean delta pseudotime", ylabel="count",
                title=f"Directionality, p = {p_perm:.3g}")
    axes[2].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, "velocity_vs_pseudotime.pdf"))

    if "X_umap" in vdata.obsm:
        vdata.obs["delta_pseudotime"] = obs["delta_pseudotime"]
        vdata.obs[args.pt_key] = obs[args.pt_key]
        sc.settings.figdir = args.outdir
        sc.pl.umap(vdata, color=[args.pt_key, "latent_time", "delta_pseudotime"],
                   cmap="coolwarm", show=False, save="_delta_pseudotime.pdf")


if __name__ == "__main__":
    main()
