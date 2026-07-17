#!/usr/bin/env python
"""
Validate that PIPseq captures the same biological signal as 10x (ground truth).
Framing: 10x is the reference; all tests ask whether PIPseq recapitulates it.

Tests:
  1. kNN label transfer (train on 10x, classify PIPseq)
     - Overall accuracy and ARI
     - Per-cluster recovery rate (what fraction of each 10x cluster is
       correctly recovered in PIPseq)
  2. Directed marker gene overlap
     - For each cluster: what fraction of 10x top markers appear in PIPseq markers?
     - Reported as recall (PIPseq recall of 10x markers), not symmetric Jaccard
  3. Pseudobulk Spearman correlation (10x vs PIPseq, marker genes only)
  4. UMAP: side-by-side 10x | PIPseq colored by joint cluster labels

Usage:
    python validate_platform_concordance.py --h5ad /path/to/integrated.h5ad
"""

import argparse
import os
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from sklearn.metrics import adjusted_rand_score, accuracy_score
from sklearn.neighbors import KNeighborsClassifier
from scipy.stats import spearmanr
from mpmath import mp, mpf, sqrt, betainc, nstr
mp.dps = 50


# ── CLI args ──────────────────────────────────────────────────────────────────

parser = argparse.ArgumentParser()
parser.add_argument('--h5ad',            required=True)
parser.add_argument('--embedding',       default='X_pca_harmony')
parser.add_argument('--cluster_key',     default='leiden')
parser.add_argument('--platform_key',    default='method')
parser.add_argument('--n_markers',       type=int, default=50)
parser.add_argument('--n_neighbors_knn', type=int, default=15)
parser.add_argument('--resolution',      type=float, default=None)
parser.add_argument('--counts_source',   default='raw', choices=['raw', 'X', 'layer'])
parser.add_argument('--counts_layer',    default='counts')
parser.add_argument('--outdir',          default='.')
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

SAMPLE_KEYS     = ['bio_condition', 'replicate']
PLATFORM_10X    = '10x'
PLATFORM_PIPSEQ = 'pipseq'


# ── Load ──────────────────────────────────────────────────────────────────────

print(f"Loading {args.h5ad} ...")
adata = sc.read_h5ad(args.h5ad)
print(adata)
print(f"adata.raw: {adata.raw.shape if adata.raw is not None else 'None'}")

assert adata.raw is not None, "adata.raw is None; set --counts_source X or layer"
assert args.embedding in adata.obsm
assert args.cluster_key in adata.obs

resolution = args.resolution
if resolution is None:
    try:
        resolution = adata.uns[args.cluster_key]['params']['resolution']
        print(f"Resolution from adata.uns: {resolution}")
    except KeyError:
        resolution = 1.0
        print(f"Resolution not found in uns; defaulting to {resolution}")


# ── 0. Filter to samples present in both platforms ────────────────────────────

adata.obs['sample_id'] = (
    adata.obs[SAMPLE_KEYS].astype(str).agg('_'.join, axis=1)
)
sample_platform = (
    adata.obs[['sample_id', args.platform_key]]
    .drop_duplicates()
    .groupby('sample_id')[args.platform_key]
    .nunique()
)
shared_samples  = sample_platform[sample_platform > 1].index.tolist()
adata_shared    = adata[adata.obs['sample_id'].isin(shared_samples)].copy()

print(f"\nShared samples ({len(shared_samples)}): {sorted(shared_samples)}")
print(f"Cells: {adata.n_obs} -> {adata_shared.n_obs} after filtering")
print(adata_shared.obs.groupby(['bio_condition', args.platform_key])['replicate'].unique())

tenx   = adata_shared[adata_shared.obs[args.platform_key] == PLATFORM_10X].copy()
pipseq = adata_shared[adata_shared.obs[args.platform_key] == PLATFORM_PIPSEQ].copy()
print(f"\n10x: {tenx.n_obs} cells | PIPseq: {pipseq.n_obs} cells")

clusters = sorted(
    adata_shared.obs[args.cluster_key].cat.categories.tolist(),
    key=lambda x: int(x)
)
palette = dict(zip(
    clusters,
    adata_shared.uns.get(f"{args.cluster_key}_colors", [None] * len(clusters))
))


# ── Helper: log-normalize from raw for marker testing ────────────────────────

def lognorm_from_raw(adata_sub, adata_ref):
    """
    Build a fresh AnnData with X = log-normalized counts from adata.raw,
    restricted to HVGs defined in adata_ref.var['highly_variable'].
    Used for rank_genes_groups (requires log-normalized, not scaled input).
    """
    hvg_names = adata_ref.var_names[adata_ref.var['highly_variable']]
    raw_var    = adata_sub.raw.var_names
    hvg_in_raw = raw_var.intersection(hvg_names)

    tmp = ad.AnnData(
        X   = adata_sub.raw[:, hvg_in_raw].X,
        obs = adata_sub.obs.copy(),
        var = pd.DataFrame(index=hvg_in_raw)
    )
    sc.pp.normalize_total(tmp, target_sum=1e4)
    sc.pp.log1p(tmp)
    tmp.obs[args.cluster_key] = (
        adata_sub.obs[args.cluster_key].values.astype('category')
    )
    return tmp


# ── Helper: pseudobulk ────────────────────────────────────────────────────────

def pseudobulk(adata_sub, cluster_col, counts_source, counts_layer=None):
    if counts_source == 'raw':
        counts   = adata_sub.raw.X
        var_names = adata_sub.raw.var_names
    elif counts_source == 'X':
        counts   = adata_sub.X
        var_names = adata_sub.var_names
    else:
        counts   = adata_sub.layers[counts_layer]
        var_names = adata_sub.var_names

    pb = {}
    for cluster in adata_sub.obs[cluster_col].cat.categories:
        mask = (adata_sub.obs[cluster_col] == cluster).values
        pb[cluster] = np.asarray(counts[mask].sum(axis=0)).flatten()
    return pd.DataFrame(pb, index=var_names)


# ── Helper: Spearman p-value (mpmath to avoid float underflow) ────────────────

def spearman_pval(rho, n):
    rho    = mpf(rho)
    n      = mpf(n)
    t_stat = rho * sqrt((n - 2) / (1 - rho ** 2))
    df     = n - 2
    x      = df / (df + t_stat ** 2)
    pval   = betainc(df / 2, mpf('0.5'), 0, x, regularized=True)
    return float(nstr(pval, 6))


# ══════════════════════════════════════════════════════════════════════════════
# STEP 1: kNN label transfer (10x -> PIPseq)
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("STEP 1: kNN Label Transfer (10x -> PIPseq)")
print("=" * 70)

X_train = tenx.obsm[args.embedding]
y_train = tenx.obs[args.cluster_key].values
X_test  = pipseq.obsm[args.embedding]
y_test  = pipseq.obs[args.cluster_key].values

knn = KNeighborsClassifier(n_neighbors=args.n_neighbors_knn)
knn.fit(X_train, y_train)
y_pred = knn.predict(X_test)

overall_accuracy = accuracy_score(y_test, y_pred)
overall_ari      = adjusted_rand_score(y_test, y_pred)
print(f"Overall accuracy: {overall_accuracy:.3f}")
print(f"Overall ARI:      {overall_ari:.3f}")

# Per-cluster recovery rate: for each 10x-defined cluster, what fraction
# of PIPseq cells assigned to that cluster are correctly identified?
# Uses precision (of predicted) and recall (of true) per cluster.
from sklearn.metrics import classification_report
report = classification_report(
    y_test, y_pred,
    labels=clusters,
    output_dict=True,
    zero_division=0
)
recovery_df = pd.DataFrame({
    c: {
        'precision': report[c]['precision'],
        'recall':    report[c]['recall'],
        'f1':        report[c]['f1-score'],
        'n_pipseq':  report[c]['support']
    }
    for c in clusters if c in report
}).T
recovery_df.index.name = 'cluster'
print("\nPer-cluster recovery (precision/recall/f1):")
print(recovery_df.round(3))

# Define well_sampled here so it's available for permutation test in step 2
# and summary section
well_sampled = recovery_df[recovery_df['n_pipseq'] >= 10].index.tolist()

# Confusion matrix (row = PIPseq true, col = 10x predicted, row-normalized)
confusion = pd.crosstab(
    pd.Series(y_test, name='PIPseq_true'),
    pd.Series(y_pred, name='10X_predicted'),
    normalize='index'
)
# Sort both axes numerically
num_order = sorted(confusion.index.tolist(),   key=lambda x: int(x))
col_order  = sorted(confusion.columns.tolist(), key=lambda x: int(x))
confusion  = confusion.loc[num_order, col_order]
confusion.to_csv(f"{args.outdir}/label_transfer_confusion_matrix.csv")
recovery_df.to_csv(f"{args.outdir}/per_cluster_recovery.csv")
print(f"Confusion matrix -> {args.outdir}/label_transfer_confusion_matrix.csv")
print(f"Recovery table   -> {args.outdir}/per_cluster_recovery.csv")

# Plot: confusion matrix heatmap
fig, ax = plt.subplots(figsize=(8, 6))
sns.heatmap(
    confusion, ax=ax, cmap='Blues', vmin=0, vmax=1,
    annot=True, fmt='.2f', annot_kws={'size': 7},
    linewidths=0.3, linecolor='lightgrey',
    cbar_kws={'label': 'Fraction of PIPseq cells', 'shrink': 0.8}
)
ax.set_xlabel("10x predicted cluster", fontsize=11, labelpad=8)
ax.set_ylabel("PIPseq true cluster",   fontsize=11, labelpad=8)
ax.set_title(
    f"Label transfer: 10x → PIPseq\n"
    f"accuracy={overall_accuracy:.3f}, ARI={overall_ari:.3f}",
    fontsize=11, pad=10
)
plt.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(f"{args.outdir}/label_transfer_confusion_matrix.{ext}",
                dpi=300, bbox_inches='tight')
plt.close()
print("Saved: label_transfer_confusion_matrix.pdf/.png")

# Plot: per-cluster recall bar chart
fig, ax = plt.subplots(figsize=(8, 4))
x = np.arange(len(recovery_df))
ax.bar(x, recovery_df['recall'],   width=0.4, label='Recall',    align='center')
ax.bar(x + 0.4, recovery_df['precision'], width=0.4, label='Precision', align='center')
ax.set_xticks(x + 0.2)
ax.set_xticklabels([f"C{c}" for c in recovery_df.index], fontsize=8)
ax.set_ylabel("Score", fontsize=11)
ax.set_xlabel("Cluster", fontsize=11)
ax.set_ylim(0, 1)
ax.axhline(overall_accuracy, color='grey', linestyle='--', linewidth=1,
           label=f"Overall accuracy ({overall_accuracy:.2f})")
ax.legend(fontsize=9, frameon=False)
ax.set_title("Per-cluster label transfer: PIPseq recall of 10x clusters", fontsize=11)
plt.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(f"{args.outdir}/per_cluster_recovery.{ext}",
                dpi=300, bbox_inches='tight')
plt.close()
print("Saved: per_cluster_recovery.pdf/.png")


# ══════════════════════════════════════════════════════════════════════════════
# STEP 2: Directed marker gene overlap (PIPseq recall of 10x markers)
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("STEP 2: Directed Marker Gene Overlap (PIPseq recall of 10x markers)")
print("=" * 70)

print("Building log-normalized HVG matrices for marker testing ...")
tenx_ln   = lognorm_from_raw(tenx,   adata_shared)
pipseq_ln = lognorm_from_raw(pipseq, adata_shared)
print(f"  10x:    {tenx_ln.n_obs} cells x {tenx_ln.n_vars} HVGs")
print(f"  PIPseq: {pipseq_ln.n_obs} cells x {pipseq_ln.n_vars} HVGs")

def get_top_markers(adata_ln, groupby, n_genes):
    sc.tl.rank_genes_groups(adata_ln, groupby=groupby, method='wilcoxon')
    markers = {}
    for cluster in adata_ln.obs[groupby].cat.categories:
        df = sc.get.rank_genes_groups_df(adata_ln, group=cluster)
        # No pvals_adj filter — taking top N by rank makes both platforms
        # comparable regardless of cell number difference (62k 10x vs 13k PIPseq).
        # pvals_adj cutoff is cell-number sensitive and would artificially reduce
        # the PIPseq marker list, inflating the apparent discordance.
        markers[cluster] = set(df.head(n_genes)['names'])
    return markers

markers_10x    = get_top_markers(tenx_ln,   args.cluster_key, args.n_markers)
markers_pipseq = get_top_markers(pipseq_ln, args.cluster_key, args.n_markers)

# Directed recall: for each cluster, what fraction of 10x top markers
# appear anywhere in the PIPseq marker list for that cluster?
overlap_results = {}
for cluster in clusters:
    m10x = markers_10x.get(cluster, set())
    mpip = markers_pipseq.get(cluster, set())
    if len(m10x) == 0:
        recall = np.nan
    else:
        recall = len(m10x & mpip) / len(m10x)
    overlap_results[cluster] = {
        'n_10x_markers':    len(m10x),
        'n_pipseq_markers': len(mpip),
        'n_shared':         len(m10x & mpip),
        'recall':           recall,
        'jaccard':          len(m10x & mpip) / len(m10x | mpip) if len(m10x | mpip) > 0 else np.nan
    }

overlap_df = pd.DataFrame(overlap_results).T
overlap_df.index.name = 'cluster'
print("\nDirected marker gene overlap (PIPseq recall of 10x markers):")
print(overlap_df.round(3))
print(f"\nMedian recall: {overlap_df['recall'].median():.3f}")

# ── Hypergeometric test for Jaccard significance ──────────────────────────
# For each cluster, test whether the observed overlap between 10x and PIPseq
# marker gene sets exceeds chance.
# Model: drawing n_pipseq_markers genes from N HVGs without replacement,
# given that n_10x_markers of those N genes are "successes" (10x markers).
# P-value = P(X >= k) where X ~ Hypergeometric(N, K, n)
#   N = total HVG space
#   K = number of 10x markers (successes in population)
#   n = number of PIPseq markers drawn
#   k = observed overlap
# Expected overlap by chance = n * K / N

from mpmath import mp, mpf, factorial, nstr as mpnstr
mp.dps = 50

def mpbinom(n, k):
    """Binomial coefficient at arbitrary precision."""
    n, k = int(n), int(k)
    if k < 0 or k > n:
        return mpf(0)
    return factorial(n) / (factorial(k) * factorial(n - k))

def hypergeom_sf_mpmath(k, N, K, n):
    """
    P(X >= k) for X ~ Hypergeometric(N, K, n) at arbitrary precision.
    Computed as sum of PMFs from k to min(K, n).
    PMF(i) = C(K,i)*C(N-K,n-i)/C(N,n)
    """
    N, K, n, k = int(N), int(K), int(n), int(k)
    denom = mpbinom(N, n)
    total = mpf(0)
    for i in range(k, min(K, n) + 1):
        total += mpbinom(K, i) * mpbinom(N - K, n - i)
    return total / denom

n_hvg = tenx_ln.n_vars
print(f"\nHypergeometric test (N={n_hvg} HVGs, K=n_10x=50, n=n_pipseq=50):")
print(f"  Expected overlap by chance = {args.n_markers}^2 / {n_hvg} = "
      f"{args.n_markers**2/n_hvg:.2f} genes")
print(f"  Expected Jaccard by chance = "
      f"{args.n_markers**2/n_hvg / (2*args.n_markers - args.n_markers**2/n_hvg):.4f}")

hyper_results = {}
for cluster in overlap_df.index:
    K = int(overlap_df.loc[cluster, 'n_10x_markers'])    # successes in population
    n = int(overlap_df.loc[cluster, 'n_pipseq_markers']) # draws
    k = int(overlap_df.loc[cluster, 'n_shared'])         # observed successes
    if K == 0 or n == 0:
        hyper_results[cluster] = {'expected': np.nan, 'pval': np.nan, 'pval_str': 'NA'}
        continue
    expected = n * K / n_hvg
    pval_mp  = hypergeom_sf_mpmath(k, n_hvg, K, n)
    pval_flt = float(pval_mp)   # may underflow to 0.0 for very small values
    pval_str = mpnstr(pval_mp, 4) if pval_flt == 0.0 else f"{pval_flt:.4e}"
    hyper_results[cluster] = {'expected': expected, 'pval': pval_flt,
                               'pval_str': pval_str}

hyper_df = pd.DataFrame(hyper_results).T
hyper_df.index.name = 'cluster'
overlap_df['expected_overlap'] = hyper_df['expected']
overlap_df['pval_hypergeom']   = hyper_df['pval']
overlap_df['pval_str']         = hyper_df['pval_str']

# Benjamini-Hochberg correction using mpmath to avoid float underflow.
# BH formula: padj[i] = min(1, p[i] * m / rank[i]) applied in rank order,
# then enforced to be monotone increasing from smallest to largest p.
valid_clusters = [c for c in overlap_df.index
                  if hyper_results[c].get('pval_str', 'NA') != 'NA'
                  and hyper_results[c]['pval_str'] != 'NA']

# Rebuild mpmath p-values from stored strings for BH
pval_mp_dict = {}
for cluster in valid_clusters:
    pval_mp_dict[cluster] = hypergeom_sf_mpmath(
        int(overlap_df.loc[cluster, 'n_shared']),
        n_hvg,
        int(overlap_df.loc[cluster, 'n_10x_markers']),
        int(overlap_df.loc[cluster, 'n_pipseq_markers'])
    )

m = len(pval_mp_dict)
sorted_clusters = sorted(pval_mp_dict, key=lambda c: pval_mp_dict[c])
padj_mp = {}
running_min = mpf(1)
for rank, cluster in enumerate(reversed(sorted_clusters), 1):
    padj_raw = pval_mp_dict[cluster] * mpf(m) / mpf(m - rank + 1)
    running_min = min(running_min, padj_raw)
    padj_mp[cluster] = running_min

overlap_df['padj_str'] = 'NA'
for cluster in valid_clusters:
    p = padj_mp[cluster]
    p_flt = float(p)
    overlap_df.loc[cluster, 'padj_str'] = (
        mpnstr(p, 4) if p_flt == 0.0 else f"{p_flt:.4e}"
    )

# Display table
display_df = overlap_df[['n_shared', 'expected_overlap', 'jaccard',
                          'pval_str', 'padj_str']].copy()
display_df.columns = ['n_shared', 'expected_overlap', 'jaccard',
                      'pval_hypergeom', 'padj_hypergeom']
print("\nHypergeometric test results per cluster:")
print(display_df)

well_jac    = overlap_df.loc[overlap_df.index.isin(well_sampled), 'jaccard'].dropna()
n_sig       = sum(
    1 for c in well_sampled
    if c in padj_mp and float(padj_mp[c]) < 0.05
)
exp_jaccard = args.n_markers**2 / n_hvg / (2*args.n_markers - args.n_markers**2/n_hvg)

print(f"\nSummary (well-sampled clusters, n={len(well_jac)}):")
print(f"  Median Jaccard:        {well_jac.median():.3f}")
print(f"  Expected by chance:    {exp_jaccard:.4f}")
print(f"  Significant (FDR<5%):  {n_sig}/{len(well_sampled)} clusters")
print(f"\n  Manuscript: median Jaccard index = {well_jac.median():.3f} vs. "
      f"{exp_jaccard:.3f} expected by chance for random gene sets of equivalent "
      f"size (n={args.n_markers}) drawn from the same HVG space (n={n_hvg} genes); "
      f"{n_sig}/{len(well_sampled)} clusters significant after "
      f"Benjamini-Hochberg correction (FDR < 5%), hypergeometric test)")


# Plot: Jaccard per cluster dot plot with cluster numbers as labels
fig, ax = plt.subplots(figsize=(9, 5))

jaccard_vals = overlap_df['jaccard'].astype(float)
null_line    = args.n_markers**2 / n_hvg / (2*args.n_markers - args.n_markers**2/n_hvg)

# Determine significance and sampling status per cluster
all_clusters_ordered = sorted(overlap_df.index.tolist(), key=lambda x: int(x))
x_pos = np.arange(len(all_clusters_ordered))

for i, cluster in enumerate(all_clusters_ordered):
    jac   = jaccard_vals.get(cluster, np.nan)
    n_pip = recovery_df.loc[cluster, 'n_pipseq'] if cluster in recovery_df.index else 0
    padj  = overlap_df.loc[cluster, 'padj_str'] if 'padj_str' in overlap_df.columns else 'NA'

    if np.isnan(jac) or n_pip == 0:
        color, zorder, alpha = '#d3d1c7', 1, 0.5
    elif n_pip < 10:
        color, zorder, alpha = '#d3d1c7', 2, 0.7
    elif padj == 'NA' or (padj not in ('NA',) and float(padj) >= 0.05):
        color, zorder, alpha = '#898781', 3, 0.9
    else:
        color, zorder, alpha = '#2a78d6', 4, 1.0

    # Draw circle
    ax.scatter(i, jac, s=520, color=color, zorder=zorder,
               alpha=alpha, linewidths=1.2,
               edgecolors='#185fa5' if color == '#2a78d6' else '#5f5e5a')
    # Cluster number inside dot
    ax.text(i, jac, str(cluster), ha='center', va='center',
            fontsize=7.5, fontweight='500', color='white' if color == '#2a78d6' else '#444441',
            zorder=zorder + 1)

# Null expectation line
ax.axhline(null_line, color='#e34948', linewidth=1.5,
           linestyle='--', label=f'Expected by chance ({null_line:.3f})', zorder=0)

# Median line (well-sampled clusters only)
well_jac_vals = [float(jaccard_vals[c]) for c in well_sampled
                 if c in jaccard_vals.index and not np.isnan(float(jaccard_vals[c]))]
obs_median = np.median(well_jac_vals)
ax.axhline(obs_median, color='#2a78d6', linewidth=1.5,
           linestyle=':', label=f'Median (well-sampled, {obs_median:.3f})', zorder=0)

ax.set_xticks([])
ax.set_ylabel("Jaccard index", fontsize=11)
ax.set_xlabel("Cluster (number shown in dot)", fontsize=11)
ax.set_xlim(-0.8, len(all_clusters_ordered) - 0.2)
ax.set_ylim(-0.03, max(jaccard_vals.dropna()) * 1.2)
ax.spines[['top', 'right']].set_visible(False)

# Legend
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#2a78d6',
           markersize=10, label='Significant (FDR < 5%)'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#898781',
           markersize=10, label='Not significant'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#d3d1c7',
           markersize=10, label=f'Low cell count (< 10 PIPseq cells)'),
    Line2D([0], [0], color='#e34948', linewidth=1.5, linestyle='--',
           label=f'Expected by chance ({null_line:.3f})'),
    Line2D([0], [0], color='#2a78d6', linewidth=1.5, linestyle=':',
           label=f'Median ({obs_median:.3f})'),
]
ax.legend(handles=legend_elements, fontsize=8, frameon=False,
          loc='upper left', bbox_to_anchor=(0, 1))

ax.set_title(f"Marker gene Jaccard index: 10x vs PIPseq\n"
             f"(top {args.n_markers} markers per cluster, {n_hvg} HVG space; "
             f"{len([c for c in well_sampled if overlap_df.loc[c,'padj_str'] not in ('NA',) and float(overlap_df.loc[c,'padj_str']) < 0.05])}/"
             f"{len(well_sampled)} well-sampled clusters FDR < 5%)",
             fontsize=10, pad=10)

plt.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(f"{args.outdir}/jaccard_per_cluster_dotplot.{ext}",
                dpi=300, bbox_inches='tight')
    print(f"Saved: {args.outdir}/jaccard_per_cluster_dotplot.{ext}")
plt.close()

# Also save display version with string p-values for supplement
display_df_out = overlap_df[['n_10x_markers','n_pipseq_markers','n_shared',
                              'recall','jaccard','expected_overlap',
                              'pval_str','padj_str']].copy()
display_df_out.columns = ['n_10x_markers','n_pipseq_markers','n_shared',
                           'recall','jaccard','expected_overlap',
                           'pval_hypergeom','padj_hypergeom']
display_df_out.to_csv(f"{args.outdir}/directed_marker_overlap_with_pvals.csv")
print(f"Overlap table (with string p-values) -> "
      f"{args.outdir}/directed_marker_overlap_with_pvals.csv")

# Plot: recall per cluster
fig, ax = plt.subplots(figsize=(8, 4))
x = np.arange(len(overlap_df))
bars = ax.bar(x, overlap_df['recall'], color='steelblue', width=0.6)
ax.set_xticks(x)
ax.set_xticklabels([f"C{c}" for c in overlap_df.index], fontsize=8)
ax.set_ylabel("Recall (fraction of 10x markers\nfound in PIPseq)", fontsize=10)
ax.set_xlabel("Cluster", fontsize=11)
ax.set_ylim(0, 1)
ax.axhline(overlap_df['recall'].median(), color='grey', linestyle='--',
           linewidth=1, label=f"Median ({overlap_df['recall'].median():.2f})")
ax.legend(fontsize=9, frameon=False)
ax.set_title(f"PIPseq recall of 10x marker genes (top {args.n_markers} per cluster)",
             fontsize=11)
plt.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(f"{args.outdir}/directed_marker_overlap.{ext}",
                dpi=300, bbox_inches='tight')
plt.close()
print("Saved: directed_marker_overlap.pdf/.png")


# ══════════════════════════════════════════════════════════════════════════════
# STEP 3: Pseudobulk Spearman correlation (marker genes only)
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("STEP 3: Pseudobulk Spearman Correlation (per-cluster marker genes)")
print("=" * 70)

# Use mean log-normalized expression per cluster (from the already-built
# lognorm matrices) rather than summed raw counts. Mean log-norm avoids
# highly expressed genes dominating the pseudobulk vector and better
# reflects per-cell transcriptional state.
#
# For each matched cluster pair, correlate only on that cluster's OWN
# 10x marker genes rather than the full union. This ensures the correlation
# reflects cluster-specific signal, not shared housekeeping expression.

def mean_expr_per_cluster(adata_ln, cluster_col):
    """Mean log-normalized expression per gene per cluster -> genes x clusters."""
    result = {}
    for cluster in adata_ln.obs[cluster_col].cat.categories:
        mask = (adata_ln.obs[cluster_col] == cluster).values
        if mask.sum() == 0:
            result[cluster] = np.zeros(adata_ln.n_vars)
        else:
            result[cluster] = np.asarray(adata_ln.X[mask].mean(axis=0)).flatten()
    return pd.DataFrame(result, index=adata_ln.var_names)

me_10x    = mean_expr_per_cluster(tenx_ln,   args.cluster_key)
me_pipseq = mean_expr_per_cluster(pipseq_ln, args.cluster_key)

shared_clusters = sorted(
    set(me_10x.columns) & set(me_pipseq.columns),
    key=lambda x: int(x)
)

# Full N×N matrix: for off-diagonal, use the row cluster's (10x) marker genes
# For diagonal (matched pair), use that cluster's own 10x marker genes
print("Computing full cluster x cluster Spearman matrix ...")
rho_matrix = pd.DataFrame(index=shared_clusters, columns=shared_clusters, dtype=float)
n_genes_diag = {}
for c_10x in shared_clusters:
    gene_set = list(markers_10x.get(c_10x, set()) & set(me_10x.index) & set(me_pipseq.index))
    if len(gene_set) < 2:
        for c_pip in shared_clusters:
            rho_matrix.loc[c_10x, c_pip] = np.nan
        n_genes_diag[c_10x] = 0
        continue
    n_genes_diag[c_10x] = len(gene_set)
    for c_pip in shared_clusters:
        rho, _ = spearmanr(me_10x.loc[gene_set, c_10x],
                           me_pipseq.loc[gene_set, c_pip])
        rho_matrix.loc[c_10x, c_pip] = rho

# Diagonal with p-values
correlations = {}
for cluster in shared_clusters:
    rho  = float(rho_matrix.loc[cluster, cluster])
    n    = n_genes_diag.get(cluster, 0)
    pval = spearman_pval(rho, n) if n >= 3 and not np.isnan(rho) else np.nan
    correlations[cluster] = {'rho': rho, 'pval': pval, 'n_genes': n}

corr_df = pd.DataFrame(correlations).T
valid   = corr_df['rho'].dropna()
print(f"\nMedian diagonal Spearman rho: {valid.median():.3f}")
print(f"Range: {valid.min():.3f} - {valid.max():.3f}")
print(f"Lowest cluster: {valid.idxmin()} (rho={valid.min():.3f})")
print(corr_df)

rho_matrix.to_csv(f"{args.outdir}/pseudobulk_spearman_matrix.csv")
corr_df.to_csv(f"{args.outdir}/pseudobulk_spearman_diagonal.csv")
print(f"Saved: pseudobulk_spearman_matrix.csv, pseudobulk_spearman_diagonal.csv")

# Plot: heatmap with diagonal highlighted
n    = len(shared_clusters)
data = rho_matrix.values.astype(float)

fig, ax = plt.subplots(figsize=(7, 6))
im = ax.imshow(data, cmap='RdYlBu_r', vmin=-1, vmax=1, aspect='equal')
cbar = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
cbar.set_label("Spearman ρ", fontsize=11)
cbar.ax.tick_params(labelsize=9)
ax.set_xticks(range(n))
ax.set_yticks(range(n))
ax.set_xticklabels([f"C{c}" for c in shared_clusters], rotation=45,
                   ha='right', fontsize=8)
ax.set_yticklabels([f"C{c}" for c in shared_clusters], fontsize=8)
ax.set_xlabel("PIPseq clusters", fontsize=11, labelpad=8)
ax.set_ylabel("10x clusters",    fontsize=11, labelpad=8)
for i in range(n):
    for j in range(n):
        val = data[i, j]
        ax.text(j, i, f"{val:.2f}", ha='center', va='center',
                fontsize=6, color='white' if abs(val) > 0.6 else 'black')
for k in range(n):
    ax.add_patch(mpatches.Rectangle(
        (k - 0.5, k - 0.5), 1, 1,
        fill=False, edgecolor='black', linewidth=1.5, zorder=3
    ))
ax.set_title(
    f"Pseudobulk Spearman ρ (marker genes only)\n"
    f"median diagonal ρ = {corr_df['rho'].median():.3f}, "
    f"range {corr_df['rho'].min():.3f}–{corr_df['rho'].max():.3f}",
    fontsize=10, pad=12
)
plt.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(f"{args.outdir}/pseudobulk_spearman_heatmap.{ext}",
                dpi=300, bbox_inches='tight')
plt.close()
print("Saved: pseudobulk_spearman_heatmap.pdf/.png")


# ══════════════════════════════════════════════════════════════════════════════
# STEP 4: UMAP side-by-side (10x | PIPseq)
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("STEP 4: UMAP by platform")
print("=" * 70)

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
for ax, (sub, name) in zip(axes, [(tenx, '10x'), (pipseq, 'PIPseq')]):
    coords = sub.obsm['X_umap']
    labels = sub.obs[args.cluster_key].astype(str).values
    for cluster in clusters:
        mask  = labels == str(cluster)
        color = palette.get(cluster)
        ax.scatter(
            coords[mask, 0], coords[mask, 1],
            s=1, alpha=0.4, c=color, rasterized=True
        )
    ax.set_title(f"{name}  (n={sub.n_obs:,})", fontsize=12)
    ax.set_xlabel("UMAP 1", fontsize=10)
    ax.set_ylabel("UMAP 2", fontsize=10)
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

handles = [
    plt.Line2D([0], [0], marker='o', color='w',
               markerfacecolor=palette.get(c, 'grey'),
               markersize=6, label=f"Cluster {c}")
    for c in clusters
]
fig.legend(handles=handles, title="Cluster", bbox_to_anchor=(1.01, 0.5),
           loc='center left', fontsize=8, title_fontsize=9, frameon=False)
fig.suptitle("UMAP colored by joint leiden clusters (shared samples only)",
             fontsize=12, y=1.01)
plt.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(f"{args.outdir}/umap_by_platform.{ext}",
                dpi=300, bbox_inches='tight')
plt.close()
print("Saved: umap_by_platform.pdf/.png")


# ══════════════════════════════════════════════════════════════════════════════
# MANUSCRIPT METRICS
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("MANUSCRIPT METRICS")
print("=" * 70)

n_clusters_10x    = tenx.obs[args.cluster_key].nunique()
n_clusters_pipseq = pipseq.obs[args.cluster_key].nunique()
print(f"Clusters (10x):    {n_clusters_10x}")
print(f"Clusters (PIPseq): {n_clusters_pipseq}")
print(f"  10x has {n_clusters_10x - n_clusters_pipseq} more clusters than PIPseq")

def genes_per_cluster(adata, cluster_key, gene_col='n_genes'):
    return adata.obs.groupby(cluster_key, observed=False)[gene_col].median()

gpc_10x    = genes_per_cluster(tenx,   args.cluster_key)
gpc_pipseq = genes_per_cluster(pipseq, args.cluster_key)

# Drop clusters with no cells in either platform
gpc_shared = pd.DataFrame({'10x': gpc_10x, 'pipseq': gpc_pipseq}).dropna()

from scipy.stats import mannwhitneyu
mwu_stat, mwu_pval = mannwhitneyu(
    gpc_shared['10x'], gpc_shared['pipseq'],
    alternative='two-sided'
)

print(f"\nMedian genes/cell per cluster:")
print(f"  10x:    median={gpc_shared['10x'].median():.0f}, "
      f"range={gpc_shared['10x'].min():.0f}-{gpc_shared['10x'].max():.0f}")
print(f"  PIPseq: median={gpc_shared['pipseq'].median():.0f}, "
      f"range={gpc_shared['pipseq'].min():.0f}-{gpc_shared['pipseq'].max():.0f}")
print(f"  Mann-Whitney U: U={mwu_stat:.0f}, p={mwu_pval:.4e} "
      f"(two-sided, n={len(gpc_shared)} shared clusters)")
print(f"\n  Manuscript: median {gpc_shared['10x'].median():.0f} genes/cell for 10x "
      f"vs. {gpc_shared['pipseq'].median():.0f} for PIPseq "
      f"(Mann-Whitney U={mwu_stat:.0f}, p={mwu_pval:.2e})")

gpc_df = gpc_shared.copy()
gpc_df.to_csv(f"{args.outdir}/genes_per_cluster.csv")
print(f"Genes per cluster -> {args.outdir}/genes_per_cluster.csv")

# Plot: boxplot of genes per cluster per platform
fig, ax = plt.subplots(figsize=(4, 5))
gpc_long = gpc_shared.reset_index().melt(
    id_vars=args.cluster_key, var_name='platform', value_name='median_genes'
)
platform_order = ['10x', 'pipseq']
colors = ['#4878CF', '#D65F5F']
for i, (platform, color) in enumerate(zip(platform_order, colors)):
    vals = gpc_long[gpc_long['platform'] == platform]['median_genes'].dropna()
    bp = ax.boxplot(
        vals, positions=[i], widths=0.4,
        patch_artist=True,
        boxprops=dict(facecolor=color, alpha=0.7),
        medianprops=dict(color='black', linewidth=2),
        whiskerprops=dict(color='black'),
        capprops=dict(color='black'),
        flierprops=dict(marker='o', markerfacecolor=color, markersize=5, alpha=0.7)
    )
    # Overlay individual cluster points
    ax.scatter(
        np.random.normal(i, 0.05, size=len(vals)),
        vals, color=color, alpha=0.6, s=20, zorder=3
    )

# Significance bracket
y_max = gpc_shared.values.max()
y_bracket = y_max * 1.08
ax.plot([0, 0, 1, 1], [y_bracket, y_bracket * 1.02,
                        y_bracket * 1.02, y_bracket], color='black', lw=1)
pval_label = (f"p = {mwu_pval:.2e}" if mwu_pval >= 1e-4
              else f"p = {mwu_pval:.2e}")
ax.text(0.5, y_bracket * 1.03, pval_label,
        ha='center', va='bottom', fontsize=9)

ax.set_xticks([0, 1])
ax.set_xticklabels(['10x Chromium', 'PIPseq T2'], fontsize=11)
ax.set_ylabel("Median genes detected per cell\n(per cluster)", fontsize=10)
ax.set_title("Genes per cluster by platform", fontsize=11, pad=10)
ax.spines[['top', 'right']].set_visible(False)

plt.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(f"{args.outdir}/genes_per_cluster_boxplot.{ext}",
                dpi=300, bbox_inches='tight')
    print(f"Saved: {args.outdir}/genes_per_cluster_boxplot.{ext}")
plt.close()

print(f"\nPseudo-bulk Spearman rho (diagonal, marker genes):")
print(f"  Median: {corr_df['rho'].median():.3f}")
print(f"  Range:  {corr_df['rho'].min():.3f} - {corr_df['rho'].max():.3f}")
print(f"  Lowest: cluster {corr_df['rho'].idxmin()} "
      f"(rho={corr_df['rho'].min():.3f})")

print(f"\nDirected marker gene overlap:")
print(f"  Median recall: {overlap_df['recall'].median():.3f}")
print(f"  Range:         {overlap_df['recall'].min():.3f} - "
      f"{overlap_df['recall'].max():.3f}")


# ══════════════════════════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

# Cell counts per cluster in PIPseq
pip_cluster_counts = pipseq.obs[args.cluster_key].value_counts()
pip_cluster_counts.index = pip_cluster_counts.index.astype(str)

# Identify outlier clusters: those with high precision but low recall
# in label transfer (systematically misassigned, not merely absent)
lt_recall_s    = recovery_df['recall'].astype(float)
lt_precision_s = recovery_df['precision'].astype(float)
outliers = recovery_df[
    (lt_precision_s >= 0.9) & (lt_recall_s < 0.3) &
    (recovery_df['n_pipseq'] > 50)
].index.tolist()

# well_sampled already defined after recovery_df above
well_no_out  = [c for c in well_sampled if c not in outliers]

lt_all   = lt_recall_s[well_sampled]
lt_clean = lt_recall_s[well_no_out]
mg_all   = overlap_df.loc[overlap_df.index.isin(well_sampled), 'recall'].dropna()
mg_clean = overlap_df.loc[overlap_df.index.isin(well_no_out),  'recall'].dropna()
jac_all  = overlap_df.loc[overlap_df.index.isin(well_sampled), 'jaccard'].dropna()
jac_clean= overlap_df.loc[overlap_df.index.isin(well_no_out),  'jaccard'].dropna()
pb_all   = corr_df.loc[corr_df.index.isin(well_sampled), 'rho'].astype(float)
pb_clean = corr_df.loc[corr_df.index.isin(well_no_out),  'rho'].astype(float)

print(f"\nOutlier clusters (high precision, low recall, n>=50 PIPseq cells): "
      f"{outliers}")
if outliers:
    for c in outliers:
        print(f"  Cluster {c}: precision={lt_precision_s[c]:.3f}, "
              f"recall={lt_recall_s[c]:.3f}, "
              f"n_pipseq={int(recovery_df.loc[c,'n_pipseq'])}")
    print(f"  -> These clusters likely reflect differential splitting between")
    print(f"     platforms, not failure to capture the population.")

print(f"\nAll well-sampled clusters (n_pipseq >= 10, n={len(well_sampled)}):")
print(f"  Label transfer recall:    "
      f"median={lt_all.median():.3f}, range={lt_all.min():.3f}-{lt_all.max():.3f}")
print(f"  Marker gene recall:       "
      f"median={mg_all.median():.3f}, range={mg_all.min():.3f}-{mg_all.max():.3f}")
print(f"  Marker gene Jaccard:      "
      f"median={jac_all.median():.3f}, range={jac_all.min():.3f}-{jac_all.max():.3f}")
print(f"  Pseudobulk Spearman rho:  "
      f"median={pb_all.median():.3f}, range={pb_all.min():.3f}-{pb_all.max():.3f}")

if outliers:
    print(f"\nExcluding outlier clusters {outliers} (n={len(well_no_out)} clusters):")
    print(f"  Label transfer recall:    "
          f"median={lt_clean.median():.3f}, range={lt_clean.min():.3f}-{lt_clean.max():.3f}")
    print(f"  Marker gene recall:       "
          f"median={mg_clean.median():.3f}, range={mg_clean.min():.3f}-{mg_clean.max():.3f}")
    print(f"  Marker gene Jaccard:      "
          f"median={jac_clean.median():.3f}, range={jac_clean.min():.3f}-{jac_clean.max():.3f}")
    print(f"  Pseudobulk Spearman rho:  "
          f"median={pb_clean.median():.3f}, range={pb_clean.min():.3f}-{pb_clean.max():.3f}")

print(f"\nOverall (all clusters, including empty):")
print(f"  Label transfer accuracy:  {overall_accuracy:.3f}")
print(f"  Label transfer ARI:       {overall_ari:.3f}")
print(f"\nAll outputs written to: {args.outdir}/")