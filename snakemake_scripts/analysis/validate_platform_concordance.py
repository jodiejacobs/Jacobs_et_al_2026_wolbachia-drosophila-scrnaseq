#!/usr/bin/env python
"""
Validate platform concordance (10x vs PIPseq) on shared samples.

Tests:
  1. Adjusted Rand Index (solo re-clustering vs joint labels)
  2. kNN label transfer (10x -> PIPseq)
  3. Marker gene Jaccard similarity per matched cluster
  4. Pseudobulk Spearman correlation per matched cluster

Usage:
    python validate_platform_concordance.py --h5ad /path/to/integrated.h5ad
"""

import argparse
import os
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import adjusted_rand_score, accuracy_score
from sklearn.neighbors import KNeighborsClassifier
from scipy.stats import spearmanr


# ── CLI args ───────────────────────────────────────────────────────────────

parser = argparse.ArgumentParser()
parser.add_argument('--h5ad', required=True, help='Path to integrated h5ad object')
parser.add_argument('--embedding', default='X_pca_harmony', help='obsm key for integrated embedding')
parser.add_argument('--cluster_key', default='leiden', help='obs key for joint cluster labels')
parser.add_argument('--platform_key', default='method', help='obs key for platform (10x/pipseq)')
parser.add_argument('--n_markers', type=int, default=50, help='top N marker genes per cluster')
parser.add_argument('--n_neighbors_knn', type=int, default=15, help='k for label transfer kNN')
parser.add_argument('--resolution', type=float, default=None,
                     help='leiden resolution for solo re-clustering; '
                          'if omitted, pulled from adata.uns[cluster_key]["params"]["resolution"]')
parser.add_argument('--counts_source', default='raw', choices=['raw', 'X', 'layer'],
                     help="where raw counts live: 'raw' = adata.raw.X, 'X' = adata.X, "
                          "'layer' = adata.layers[--counts_layer]")
parser.add_argument('--counts_layer', default='counts',
                     help="layer holding raw counts (only used if --counts_source layer)")
parser.add_argument('--outdir', default='.', help='directory to write output tables')
args = parser.parse_args()


os.makedirs(args.outdir, exist_ok=True)

SAMPLE_KEYS = ['bio_condition', 'replicate']
PLATFORM_VALUES = ('10x', 'pipseq')  # confirm these match adata.obs[platform_key].unique()


# ── Load ───────────────────────────────────────────────────────────────────

print(f"Loading {args.h5ad} ...")
adata = sc.read_h5ad(args.h5ad)
print(adata)

if adata.raw is not None:
    print(f"adata.raw: {adata.raw.shape[0]} obs x {adata.raw.shape[1]} vars")
else:
    print("adata.raw is None")

if args.counts_source == 'layer':
    assert args.counts_layer in adata.layers, (
        f"'{args.counts_layer}' not found in adata.layers. "
        f"Available layers: {list(adata.layers.keys())}"
    )
elif args.counts_source == 'raw':
    assert adata.raw is not None, (
        "adata.raw is None but --counts_source raw was requested. "
        "Use --counts_source X or --counts_source layer instead."
    )
# 'X' needs no check beyond adata.X existing, which it always does
assert args.embedding in adata.obsm, (
    f"'{args.embedding}' not found in adata.obsm. "
    f"Available obsm: {list(adata.obsm.keys())}"
)
assert args.cluster_key in adata.obs, (
    f"'{args.cluster_key}' not found in adata.obs."
)

platform_vals = set(adata.obs[args.platform_key].unique())
assert platform_vals == set(PLATFORM_VALUES), (
    f"Expected platform values {PLATFORM_VALUES}, found {platform_vals}. "
    f"Update PLATFORM_VALUES in the script."
)

resolution = args.resolution
if resolution is None:
    try:
        resolution = adata.uns[args.cluster_key]['params']['resolution']
        print(f"Using resolution from adata.uns['{args.cluster_key}']['params']: {resolution}")
    except KeyError:
        resolution = 1.0
        print(f"Could not find resolution in adata.uns; defaulting to {resolution}. "
              f"Pass --resolution explicitly to override.")


# ── 0. Filter to samples present in both platforms ────────────────────────

adata.obs['sample_id'] = (
    adata.obs[SAMPLE_KEYS]
    .astype(str)
    .agg('_'.join, axis=1)
)

sample_platform = (
    adata.obs[['sample_id', args.platform_key]]
    .drop_duplicates()
    .groupby('sample_id')[args.platform_key]
    .nunique()
)

shared_samples = sample_platform[sample_platform > 1].index.tolist()
print(f"\nSamples present in both platforms ({len(shared_samples)} of "
      f"{adata.obs['sample_id'].nunique()} total):")
for s in sorted(shared_samples):
    print(f"  {s}")

adata_shared = adata[adata.obs['sample_id'].isin(shared_samples)].copy()
print(f"\nCells before filter: {adata.n_obs}, after filter: {adata_shared.n_obs}")

print("\nFinal sample x platform breakdown:")
print(adata_shared.obs.groupby(['bio_condition', args.platform_key])['replicate'].unique())

tenx = adata_shared[adata_shared.obs[args.platform_key] == '10x'].copy()
pipseq = adata_shared[adata_shared.obs[args.platform_key] == 'pipseq'].copy()
print(f"\n10X: {tenx.n_obs} cells, PIPseq: {pipseq.n_obs} cells")


# ── 1. Adjusted Rand Index (solo re-clustering vs joint labels) ───────────

print("\n" + "=" * 70)
print("STEP 1: Adjusted Rand Index (solo re-clustering vs joint labels)")
print("=" * 70)

ari_results = {}
for platform_adata, name in [(tenx, '10x'), (pipseq, 'pipseq')]:
    sc.pp.neighbors(platform_adata, use_rep=args.embedding)
    sc.tl.leiden(platform_adata, resolution=resolution, key_added='leiden_solo')

    ari = adjusted_rand_score(
        platform_adata.obs[args.cluster_key],
        platform_adata.obs['leiden_solo']
    )
    ari_results[name] = ari
    print(f"ARI ({name}, solo vs joint): {ari:.3f}")


# ── 2. kNN Label Transfer (train on 10X, predict on PIPseq) ───────────────

print("\n" + "=" * 70)
print("STEP 2: kNN Label Transfer (10x -> PIPseq)")
print("=" * 70)

X_train = tenx.obsm[args.embedding]
y_train = tenx.obs[args.cluster_key].values

X_test = pipseq.obsm[args.embedding]
y_test = pipseq.obs[args.cluster_key].values

knn = KNeighborsClassifier(n_neighbors=args.n_neighbors_knn)
knn.fit(X_train, y_train)
y_pred = knn.predict(X_test)

accuracy = accuracy_score(y_test, y_pred)
ari_transfer = adjusted_rand_score(y_test, y_pred)
print(f"Label transfer accuracy: {accuracy:.3f}")
print(f"Label transfer ARI:      {ari_transfer:.3f}")

confusion = pd.crosstab(
    pd.Series(y_test, name='PIPseq_true'),
    pd.Series(y_pred, name='10X_predicted'),
    normalize='index'
)
confusion_path = f"{args.outdir}/label_transfer_confusion_matrix.csv"
confusion.to_csv(confusion_path)
print(f"Confusion matrix written to {confusion_path}")


# ── 3. Marker Gene Jaccard Similarity ──────────────────────────────────────

print("\n" + "=" * 70)
print("STEP 3: Marker Gene Jaccard Similarity")
print("=" * 70)

def get_top_markers(adata, groupby, n_genes):
    sc.tl.rank_genes_groups(adata, groupby=groupby, method='wilcoxon')
    markers = {}
    for cluster in adata.obs[groupby].cat.categories:
        genes = sc.get.rank_genes_groups_df(adata, group=cluster)
        genes = genes[genes['pvals_adj'] < 0.05]
        markers[cluster] = set(genes.head(n_genes)['names'])
    return markers

markers_10x = get_top_markers(tenx, args.cluster_key, args.n_markers)
markers_pipseq = get_top_markers(pipseq, args.cluster_key, args.n_markers)

def jaccard(set_a, set_b):
    if len(set_a | set_b) == 0:
        return np.nan
    return len(set_a & set_b) / len(set_a | set_b)

clusters = sorted(
    set(markers_10x.keys()) & set(markers_pipseq.keys()),
    key=lambda x: int(x)
)
jaccard_matrix = pd.DataFrame(index=clusters, columns=clusters, dtype=float)
for c1 in clusters:
    for c2 in clusters:
        jaccard_matrix.loc[c1, c2] = jaccard(markers_10x[c1], markers_pipseq[c2])

diag_vals = np.diag(jaccard_matrix.loc[clusters, clusters].values)
print(f"Median diagonal Jaccard (matched clusters): {np.nanmedian(diag_vals):.3f}")

jaccard_path = f"{args.outdir}/marker_gene_jaccard_matrix.csv"
jaccard_matrix.to_csv(jaccard_path)
print(f"Jaccard matrix written to {jaccard_path}")


# ── 4. Pseudobulk Spearman Correlation ─────────────────────────────────────

print("\n" + "=" * 70)
print("STEP 4: Pseudobulk Spearman Correlation")
print("=" * 70)

def pseudobulk(adata, cluster_col, counts_source, counts_layer=None):
    """
    Sum counts per cluster -> genes x clusters DataFrame.
    Uses adata.raw.var_names when counts_source == 'raw', since
    adata.raw may contain genes (e.g. pre-HVG-filtering) not present
    in adata.var_names.
    """
    if counts_source == 'raw':
        counts = adata.raw.X
        var_names = adata.raw.var_names
    elif counts_source == 'X':
        counts = adata.X
        var_names = adata.var_names
    elif counts_source == 'layer':
        counts = adata.layers[counts_layer]
        var_names = adata.var_names
    else:
        raise ValueError(f"Unknown counts_source: {counts_source}")

    pb = {}
    for cluster in adata.obs[cluster_col].cat.categories:
        mask = (adata.obs[cluster_col] == cluster).values
        sub = counts[mask]
        summed = np.asarray(sub.sum(axis=0)).flatten()
        pb[cluster] = summed
    return pd.DataFrame(pb, index=var_names)

pb_10x = pseudobulk(tenx, args.cluster_key, args.counts_source, args.counts_layer)
pb_pipseq = pseudobulk(pipseq, args.cluster_key, args.counts_source, args.counts_layer)

# Restrict to genes shared between the two pseudobulk matrices
# (relevant if --counts_source raw and raw var sets differ for any reason)
shared_genes = pb_10x.index.intersection(pb_pipseq.index)
if len(shared_genes) < len(pb_10x.index) or len(shared_genes) < len(pb_pipseq.index):
    print(f"Note: restricting pseudobulk correlation to {len(shared_genes)} genes "
          f"shared between platforms (10x: {len(pb_10x.index)}, pipseq: {len(pb_pipseq.index)})")
pb_10x = pb_10x.loc[shared_genes]
pb_pipseq = pb_pipseq.loc[shared_genes]

from mpmath import mp, mpf, sqrt, betainc, nstr
mp.dps = 50  # 50 decimal places of precision to avoid underflow

def spearman_pval(rho, n):
    """
    Two-sided Spearman p-value via the incomplete beta function at
    arbitrary precision (mpmath). scipy underflows to 0.0 for n > ~1000
    at typical rho values, so float64 is insufficient here.
    t = rho * sqrt((n-2) / (1 - rho^2))
    p = I(df/(df+t^2), df/2, 1/2)  [two-sided by symmetry]
    """
    rho = mpf(rho)
    n   = mpf(n)
    t_stat = rho * sqrt((n - 2) / (1 - rho ** 2))
    df = n - 2
    x  = df / (df + t_stat ** 2)
    pval = betainc(df / 2, mpf('0.5'), 0, x, regularized=True)
    return float(nstr(pval, 6))   # convert back to float-representable string then float

shared_clusters = sorted(
    set(pb_10x.columns) & set(pb_pipseq.columns),
    key=lambda x: int(x)
)
n_genes = len(shared_genes)
correlations = {}
for cluster in shared_clusters:
    rho, _ = spearmanr(pb_10x[cluster], pb_pipseq[cluster])
    pval = spearman_pval(rho, n_genes)
    correlations[cluster] = {'rho': rho, 'pval': pval}

corr_df = pd.DataFrame(correlations).T
print(f"Median pseudobulk Spearman rho: {corr_df['rho'].median():.3f}")
print(corr_df)

corr_path = f"{args.outdir}/pseudobulk_spearman_correlation.csv"
corr_df.to_csv(corr_path)
print(f"Correlation table written to {corr_path}")


# ── Summary ─────────────────────────────────────────────────────────────────

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"ARI (10x solo vs joint):      {ari_results['10x']:.3f}")
print(f"ARI (pipseq solo vs joint):   {ari_results['pipseq']:.3f}")
print(f"Label transfer accuracy:      {accuracy:.3f}")
print(f"Label transfer ARI:           {ari_transfer:.3f}")
print(f"Median marker gene Jaccard:   {np.nanmedian(diag_vals):.3f}")
print(f"Median pseudobulk Spearman:   {corr_df['rho'].median():.3f}")