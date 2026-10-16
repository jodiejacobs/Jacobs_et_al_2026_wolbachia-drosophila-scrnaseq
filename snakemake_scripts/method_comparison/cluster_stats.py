"""
cluster_stats.py
================
Report per-cluster statistics from an annotated h5ad file.

Usage
-----
python cluster_stats.py \
    --h5ad results/nmf_programs/adata_with_programs.h5ad \
    --leiden leiden_res0.5 \
    --out results/cluster_stats.csv
"""

import argparse
import numpy as np
import pandas as pd
import scanpy as sc


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad",    required=True)
    parser.add_argument("--leiden",  default="leiden",
                        help="Which leiden column to use (default: leiden)")
    parser.add_argument("--out",     default="cluster_stats.csv")
    args = parser.parse_args()

    print(f"Loading {args.h5ad} ...")
    adata = sc.read_h5ad(args.h5ad)
    obs   = adata.obs
    n_total = len(obs)

    leiden_col = args.leiden
    if leiden_col not in obs.columns:
        raise ValueError(f"'{leiden_col}' not found in obs. "
                         f"Available: {obs.columns.tolist()}")

    clusters = sorted(obs[leiden_col].unique())
    rows = []

    for c in clusters:
        mask = obs[leiden_col] == c
        cluster_obs = obs[mask]
        n_cells = mask.sum()

        # ── Majority bio_condition ────────────────────────────────────────────
        cond_counts   = cluster_obs["bio_condition"].value_counts()
        majority_cond = cond_counts.index[0]
        majority_pct  = round(cond_counts.iloc[0] / n_cells * 100, 2)

        # ── Wolbachia titer ───────────────────────────────────────────────────
        titer = cluster_obs["wolbachia_titer"].dropna()
        titer_median = np.median(titer) if len(titer) > 0 else np.nan
        titer_mean   = np.mean(titer)   if len(titer) > 0 else np.nan
        titer_n      = len(titer)

        rows.append({
            "cluster":               c,
            "n_cells":               n_cells,
            "pct_of_sample":         round(n_cells / n_total * 100, 2),
            "majority_bio_condition": majority_cond,
            "majority_pct":          majority_pct,
            "titer_n":               titer_n,
            "titer_median":          round(titer_median, 4) if not np.isnan(titer_median) else np.nan,
            "titer_mean":            round(titer_mean, 4)   if not np.isnan(titer_mean)   else np.nan,
        })

    df = pd.DataFrame(rows).sort_values("cluster").reset_index(drop=True)

    print(df.to_string(index=False))
    df.to_csv(args.out, index=False)
    print(f"\nSaved -> {args.out}")


if __name__ == "__main__":
    main()
