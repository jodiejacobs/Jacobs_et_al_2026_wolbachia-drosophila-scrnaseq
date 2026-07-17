"""
Test agreement between Leiden clusters and a categorical annotation
(e.g. cell cycle phase) in an AnnData object.

Gives you two complementary things:
  1. Chi-square test of independence on the contingency table
     -> tests whether phase is distributed non-randomly across clusters
     -> statistic + p-value, but p shrinks trivially with more cells,
        so also report Cramer's V as an effect size (0-1).
  2. Adjusted Rand Index (ARI) between leiden and phase labels
     -> measures direct label agreement, corrected for chance
     -> ARI has no closed-form null, so p-value comes from a
        label-permutation test (shuffle one column, recompute ARI,
        N times -> empirical p-value).

Usage:
    python cluster_vs_category_test.py adata.h5ad --cluster_key leiden --category_key phase
"""
import argparse
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score


def cramers_v(chi2_stat, contingency_table):
    n = contingency_table.to_numpy().sum()
    r, k = contingency_table.shape
    return np.sqrt((chi2_stat / n) / (min(r - 1, k - 1)))


def permutation_ari_pvalue(labels_a, labels_b, n_perm=1000, seed=0):
    rng = np.random.default_rng(seed)
    observed = adjusted_rand_score(labels_a, labels_b)
    labels_b = np.asarray(labels_b)
    null = np.empty(n_perm)
    for i in range(n_perm):
        null[i] = adjusted_rand_score(labels_a, rng.permutation(labels_b))
    p = (np.sum(null >= observed) + 1) / (n_perm + 1)
    return observed, p, null


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("h5ad", help="Path to AnnData .h5ad file")
    ap.add_argument("--cluster_key", default="leiden")
    ap.add_argument("--category_key", default="phase")
    ap.add_argument("--n_perm", type=int, default=1000)
    args = ap.parse_args()

    import anndata as ad
    adata = ad.read_h5ad(args.h5ad)
    obs = adata.obs
    clusters = obs[args.cluster_key].astype(str)
    category = obs[args.category_key].astype(str)

    # 1. Chi-square test of independence
    ct = pd.crosstab(clusters, category)
    chi2_stat, chi2_p, dof, expected = chi2_contingency(ct)
    v = cramers_v(chi2_stat, ct)

    print("Contingency table (clusters x category):")
    print(ct)
    print(f"\nChi-square test: chi2 = {chi2_stat:.2f}, dof = {dof}, p = {chi2_p:.3e}")
    print(f"Cramer's V (effect size, 0-1): {v:.3f}")

    # 2. Adjusted Rand Index + permutation p-value
    ari, p_perm, null_dist = permutation_ari_pvalue(
        clusters, category, n_perm=args.n_perm
    )
    nmi = normalized_mutual_info_score(clusters, category)
    print(f"\nAdjusted Rand Index: {ari:.3f}")
    print(f"Permutation p-value (n={args.n_perm}): {p_perm:.3e}")
    print(f"Normalized Mutual Information: {nmi:.3f}")


if __name__ == "__main__":
    main()
