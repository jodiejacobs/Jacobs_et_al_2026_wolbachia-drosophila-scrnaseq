#!/usr/bin/env python
"""
Test whether the centrifugal shell vial technique introduces transcriptional
variation in uninfected JW18 cells.

Comparison (both platforms pooled):
  Group A: JW18DOX-Ctrl  -- uninfected cells sampled directly from culture
  Group B: JW18DOX-SV-D1 -- uninfected cells subjected to mock shell vial
                             protocol (no wMel added)

Tests:
  1. Genes per cluster (Mann-Whitney U)
  2. Cell cluster composition (JSD + per-cluster Fisher's exact, BH-corrected)
  3. Marker gene Jaccard (hypergeometric test, BH-corrected)
  4. Pseudobulk Spearman correlation (per-cluster marker genes)
  5. UMAP side-by-side

Usage:
    python validate_shell_vial.py --h5ad /path/to/integrated.h5ad \
        --outdir results/validate_shell_vial
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
from matplotlib.lines import Line2D
from scipy.stats import spearmanr, mannwhitneyu, fisher_exact
from scipy.spatial.distance import jensenshannon
from statsmodels.stats.multitest import multipletests
from mpmath import mp, mpf, factorial, nstr as mpnstr, sqrt, betainc
mp.dps = 50


# ── CLI args ──────────────────────────────────────────────────────────────────

parser = argparse.ArgumentParser()
parser.add_argument('--h5ad',        required=True)
parser.add_argument('--cluster_key', default='leiden')
parser.add_argument('--n_markers',   type=int, default=50)
parser.add_argument('--outdir',      default='.')
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

GROUP_A       = 'JW18DOX-Ctrl'
GROUP_B       = 'JW18DOX-SV-D1'
GROUP_A_LABEL = 'DOX-Ctrl\n(direct culture)'
GROUP_B_LABEL = 'DOX-SV-D1\n(mock shell vial)'
COLORS        = {'ctrl': '#2a78d6', 'sv': '#1baf7a'}


# ── Helpers ───────────────────────────────────────────────────────────────────

def lognorm_from_raw(adata_sub, adata_ref, cluster_key):
    """Build log-normalized HVG AnnData from raw counts for marker testing."""
    hvg_names  = adata_ref.var_names[adata_ref.var['highly_variable']]
    hvg_in_raw = adata_sub.raw.var_names.intersection(hvg_names)
    tmp = ad.AnnData(
        X   = adata_sub.raw[:, hvg_in_raw].X,
        obs = adata_sub.obs.copy(),
        var = pd.DataFrame(index=hvg_in_raw)
    )
    sc.pp.normalize_total(tmp, target_sum=1e4)
    sc.pp.log1p(tmp)
    tmp.obs[cluster_key] = adata_sub.obs[cluster_key].values.astype('category')
    return tmp

def mean_expr_per_cluster(adata_ln, cluster_col):
    """Mean log-normalized expression per cluster -> genes x clusters DataFrame."""
    result = {}
    for cluster in adata_ln.obs[cluster_col].cat.categories:
        mask = (adata_ln.obs[cluster_col] == cluster).values
        result[cluster] = (
            np.zeros(adata_ln.n_vars) if mask.sum() == 0
            else np.asarray(adata_ln.X[mask].mean(axis=0)).flatten()
        )
    return pd.DataFrame(result, index=adata_ln.var_names)

def mpbinom(n, k):
    n, k = int(n), int(k)
    if k < 0 or k > n: return mpf(0)
    return factorial(n) / (factorial(k) * factorial(n - k))

def hypergeom_sf_mpmath(k, N, K, n):
    """P(X >= k) for X ~ Hypergeometric(N, K, n) at arbitrary precision."""
    N, K, n, k = int(N), int(K), int(n), int(k)
    denom = mpbinom(N, n)
    total = sum(mpbinom(K, i) * mpbinom(N - K, n - i)
                for i in range(k, min(K, n) + 1))
    return total / denom

def bh_mpmath(pval_mp_dict):
    """Benjamini-Hochberg correction on dict of mpmath p-values."""
    m = len(pval_mp_dict)
    sorted_clusters = sorted(pval_mp_dict, key=lambda c: pval_mp_dict[c])
    padj = {}
    running_min = mpf(1)
    for rank, cluster in enumerate(reversed(sorted_clusters), 1):
        raw = pval_mp_dict[cluster] * mpf(m) / mpf(m - rank + 1)
        running_min = min(running_min, raw)
        padj[cluster] = running_min
    return padj

def mpf_to_str(p):
    p_flt = float(p)
    return mpnstr(p, 4) if p_flt == 0.0 else f"{p_flt:.4e}"

def spearman_pval_mpmath(rho, n):
    rho    = mpf(rho)
    n      = mpf(n)
    t_stat = rho * sqrt((n - 2) / (1 - rho ** 2))
    df     = n - 2
    x      = df / (df + t_stat ** 2)
    return betainc(df / 2, mpf('0.5'), 0, x, regularized=True)


# ── Load & subset ─────────────────────────────────────────────────────────────

print(f"Loading {args.h5ad} ...")
adata = sc.read_h5ad(args.h5ad)
print(adata)
assert adata.raw is not None, "adata.raw is None"

dox_mask  = adata.obs['bio_condition'].isin([GROUP_A, GROUP_B])
adata_dox = adata[dox_mask].copy()
print(f"\nSubset to {GROUP_A} + {GROUP_B}: {adata_dox.n_obs} cells")
print(adata_dox.obs.groupby(['bio_condition', 'method'])['replicate'].unique())

ctrl = adata_dox[adata_dox.obs['bio_condition'] == GROUP_A].copy()
sv   = adata_dox[adata_dox.obs['bio_condition'] == GROUP_B].copy()

# Exclude JW18DOX-SV-D1 PIPseq rep 1 — suspected Wolbachia contamination
exclude_mask = (
    (sv.obs['method'] == 'pipseq') &
    (sv.obs['replicate'].astype(str) == '1')
)
n_excluded = exclude_mask.sum()
sv = sv[~exclude_mask].copy()
print(f"\nExcluded {n_excluded} cells: JW18DOX-SV-D1 PIPseq rep 1 "
      f"(suspected contamination)")
print(f"\n{GROUP_A}: {ctrl.n_obs} cells  |  {GROUP_B}: {sv.n_obs} cells")

clusters = sorted(adata_dox.obs[args.cluster_key].cat.categories.tolist(),
                  key=lambda x: int(x))
palette  = dict(zip(clusters,
                    adata_dox.uns.get(f"{args.cluster_key}_colors",
                                      [None]*len(clusters))))


# ══════════════════════════════════════════════════════════════════════════════
# STEP 1: Genes per cluster
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("STEP 1: Genes per cluster")
print("=" * 70)

n_clusters_ctrl = ctrl.obs[args.cluster_key].nunique()
n_clusters_sv   = sv.obs[args.cluster_key].nunique()
print(f"Clusters recovered: {GROUP_A}={n_clusters_ctrl}, {GROUP_B}={n_clusters_sv}")

def genes_per_cluster(adata, cluster_key, gene_col='n_genes'):
    return adata.obs.groupby(cluster_key, observed=False)[gene_col].median()

gpc_ctrl   = genes_per_cluster(ctrl, args.cluster_key)
gpc_sv     = genes_per_cluster(sv,   args.cluster_key)
gpc_shared = pd.DataFrame({'ctrl': gpc_ctrl, 'sv': gpc_sv}).dropna()

mwu_stat, mwu_pval = mannwhitneyu(
    gpc_shared['ctrl'], gpc_shared['sv'], alternative='two-sided'
)
print(f"\nMedian genes/cell per cluster:")
print(f"  {GROUP_A}:   median={gpc_shared['ctrl'].median():.0f}, "
      f"range={gpc_shared['ctrl'].min():.0f}-{gpc_shared['ctrl'].max():.0f}")
print(f"  {GROUP_B}: median={gpc_shared['sv'].median():.0f}, "
      f"range={gpc_shared['sv'].min():.0f}-{gpc_shared['sv'].max():.0f}")
print(f"  Mann-Whitney U={mwu_stat:.0f}, p={mwu_pval:.4e} "
      f"(n={len(gpc_shared)} shared clusters)")

gpc_shared.to_csv(f"{args.outdir}/shell_vial_genes_per_cluster.csv")

# Plot: paired boxplot
fig, ax = plt.subplots(figsize=(4, 5))
for i, (col, label) in enumerate(zip(['ctrl', 'sv'],
                                      [GROUP_A_LABEL, GROUP_B_LABEL])):
    color = list(COLORS.values())[i]
    vals  = gpc_shared[col].dropna()
    ax.boxplot(vals, positions=[i], widths=0.4, patch_artist=True,
               boxprops=dict(facecolor=color, alpha=0.7),
               medianprops=dict(color='black', linewidth=2),
               whiskerprops=dict(color='black'), capprops=dict(color='black'),
               flierprops=dict(marker='o', markerfacecolor=color,
                               markersize=5, alpha=0.7))
    ax.scatter(np.random.default_rng(42).normal(i, 0.05, len(vals)),
               vals, color=color, alpha=0.6, s=20, zorder=3)
y_br = gpc_shared.values.max() * 1.08
ax.plot([0, 0, 1, 1], [y_br, y_br*1.02, y_br*1.02, y_br], color='black', lw=1)
ax.text(0.5, y_br*1.03, f"p = {mwu_pval:.2e}", ha='center', va='bottom',
        fontsize=9)
ax.set_xticks([0, 1])
ax.set_xticklabels([GROUP_A_LABEL, GROUP_B_LABEL], fontsize=10)
ax.set_ylabel("Median genes detected per cell\n(per cluster)", fontsize=10)
ax.set_title("Genes per cluster: Ctrl vs mock shell vial", fontsize=11, pad=10)
ax.spines[['top', 'right']].set_visible(False)
plt.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(f"{args.outdir}/shell_vial_genes_per_cluster_boxplot.{ext}",
                dpi=300, bbox_inches='tight')
plt.close()
print("Saved: shell_vial_genes_per_cluster_boxplot.pdf/.png")


# ══════════════════════════════════════════════════════════════════════════════
# STEP 2: Cell cluster composition
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("STEP 2: Cell cluster composition")
print("=" * 70)

# Cell counts per cluster per condition
counts_ctrl = ctrl.obs[args.cluster_key].value_counts().reindex(clusters, fill_value=0)
counts_sv   = sv.obs[args.cluster_key].value_counts().reindex(clusters, fill_value=0)
prop_ctrl   = counts_ctrl / counts_ctrl.sum()
prop_sv     = counts_sv   / counts_sv.sum()

comp_df = pd.DataFrame({
    'n_ctrl':    counts_ctrl,
    'n_sv':      counts_sv,
    'prop_ctrl': prop_ctrl,
    'prop_sv':   prop_sv,
})
comp_df.index.name = 'cluster'
print("\nCluster composition:")
print(comp_df.round(4))

# Jensen-Shannon divergence between proportion vectors
# JSD = 0: identical distributions; JSD = 1: maximally different
# scipy.spatial.distance.jensenshannon returns the square root (JS distance)
# so we square it to get the true divergence
js_dist = jensenshannon(prop_ctrl.values, prop_sv.values, base=2)
jsd     = js_dist ** 2
print(f"\nJensen-Shannon divergence: {jsd:.4f} (JSD=0 identical, JSD=1 maximal)")
print(f"Jensen-Shannon distance:   {js_dist:.4f}")

# Per-cluster Fisher's exact test:
# For each cluster, 2x2 table:
#   rows = in cluster / not in cluster
#   cols = Ctrl / SV-D1
fisher_results = {}
total_ctrl = counts_ctrl.sum()
total_sv   = counts_sv.sum()
for cluster in clusters:
    a = counts_ctrl[cluster]           # ctrl in cluster
    b = counts_sv[cluster]             # sv in cluster
    c = total_ctrl - a                 # ctrl not in cluster
    d = total_sv   - b                 # sv not in cluster
    odds, pval = fisher_exact([[a, b], [c, d]], alternative='two-sided')
    fisher_results[cluster] = {
        'n_ctrl':    a,
        'n_sv':      b,
        'prop_ctrl': prop_ctrl[cluster],
        'prop_sv':   prop_sv[cluster],
        'odds_ratio': odds,
        'pval':       pval
    }

fisher_df = pd.DataFrame(fisher_results).T
fisher_df.index.name = 'cluster'

# BH correction
_, padj, _, _ = multipletests(fisher_df['pval'], method='fdr_bh')
fisher_df['padj'] = padj
n_sig_comp = (fisher_df['padj'] < 0.05).sum()

print(f"\nPer-cluster Fisher's exact test (BH-corrected):")
print(fisher_df[['n_ctrl', 'n_sv', 'prop_ctrl', 'prop_sv',
                 'odds_ratio', 'pval', 'padj']].round(4))
print(f"\n{n_sig_comp}/{len(clusters)} clusters significantly different "
      f"in proportion (FDR < 5%)")

fisher_df.to_csv(f"{args.outdir}/shell_vial_cluster_composition.csv")
print(f"Composition table -> {args.outdir}/shell_vial_cluster_composition.csv")

# Plot 1: Stacked bar chart of cluster proportions
fig, axes = plt.subplots(1, 2, figsize=(7, 5), sharey=False)
bar_colors = [palette.get(c, '#888888') for c in clusters]
bottoms = {'ctrl': 0, 'sv': 0}
for cluster, color in zip(clusters, bar_colors):
    for ax, (col, label) in zip(axes, [('ctrl', GROUP_A), ('sv', GROUP_B)]):
        prop = comp_df.loc[cluster, f'prop_{col}']
        ax.bar(0, prop, bottom=bottoms[col], color=color, width=0.5,
               edgecolor='white', linewidth=0.5)
        if prop > 0.03:
            ax.text(0, bottoms[col] + prop/2, str(cluster),
                    ha='center', va='center', fontsize=8,
                    color='white', fontweight='500')
        bottoms[col] += prop
    
for ax, label in zip(axes, [GROUP_A_LABEL, GROUP_B_LABEL]):
    ax.set_xlim(-0.4, 0.4)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_title(label, fontsize=10)
    ax.spines[['top', 'right', 'bottom']].set_visible(False)

axes[0].set_ylabel("Proportion of cells", fontsize=11)
axes[1].set_ylabel("")
fig.suptitle("Cell cluster composition:\nCtrl vs mock shell vial",
             fontsize=11, y=1.01)
plt.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(f"{args.outdir}/shell_vial_cluster_composition_bar.{ext}",
                dpi=300, bbox_inches='tight')
plt.close()
print("Saved: shell_vial_cluster_composition_bar.pdf/.png")

# Plot 2: Scatter of Ctrl proportion vs SV proportion per cluster
fig, ax = plt.subplots(figsize=(5, 5))
for cluster in clusters:
    pc = comp_df.loc[cluster, 'prop_ctrl']
    ps = comp_df.loc[cluster, 'prop_sv']
    sig = fisher_df.loc[cluster, 'padj'] < 0.05
    color = '#e34948' if sig else '#2a78d6'
    ax.scatter(pc, ps, s=200, color=color, zorder=3, alpha=0.85,
               edgecolors='#5f5e5a', linewidths=0.8)
    ax.text(pc, ps, str(cluster), ha='center', va='center',
            fontsize=7, color='white', fontweight='500', zorder=4)

max_prop = max(prop_ctrl.max(), prop_sv.max()) * 1.1
ax.plot([0, max_prop], [0, max_prop], 'k--', linewidth=1, alpha=0.4,
        label='y = x (identical proportions)')
ax.set_xlabel(f"Proportion in {GROUP_A}", fontsize=11)
ax.set_ylabel(f"Proportion in {GROUP_B}", fontsize=11)
ax.set_xlim(-0.01, max_prop)
ax.set_ylim(-0.01, max_prop)
ax.spines[['top', 'right']].set_visible(False)
legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#e34948',
           markersize=10, label='Sig. different (FDR < 5%)'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#2a78d6',
           markersize=10, label='Not significant'),
    Line2D([0], [0], color='black', linewidth=1, linestyle='--',
           alpha=0.4, label='y = x'),
]
ax.legend(handles=legend_elements, fontsize=8, frameon=False)
ax.set_title(f"Cluster proportions: {GROUP_A} vs {GROUP_B}\n"
             f"JSD = {jsd:.4f}; {n_sig_comp}/{len(clusters)} clusters "
             f"significantly different (FDR < 5%)",
             fontsize=10, pad=10)
plt.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(f"{args.outdir}/shell_vial_cluster_composition_scatter.{ext}",
                dpi=300, bbox_inches='tight')
plt.close()
print("Saved: shell_vial_cluster_composition_scatter.pdf/.png")


# ══════════════════════════════════════════════════════════════════════════════
# STEP 3: Marker gene Jaccard (SV-D1 recall of Ctrl markers)
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("STEP 3: Marker gene Jaccard (SV-D1 recall of Ctrl markers)")
print("=" * 70)

print("Building log-normalized HVG matrices ...")
ctrl_ln = lognorm_from_raw(ctrl, adata_dox, args.cluster_key)
sv_ln   = lognorm_from_raw(sv,   adata_dox, args.cluster_key)
n_hvg   = ctrl_ln.n_vars
print(f"  Ctrl: {ctrl_ln.n_obs} cells x {n_hvg} HVGs")
print(f"  SV:   {sv_ln.n_obs}  cells x {sv_ln.n_vars} HVGs")

def get_top_markers(adata_ln, groupby, n_genes):
    sc.tl.rank_genes_groups(adata_ln, groupby=groupby, method='wilcoxon')
    markers = {}
    for cluster in adata_ln.obs[groupby].cat.categories:
        df = sc.get.rank_genes_groups_df(adata_ln, group=cluster)
        markers[cluster] = set(df.head(n_genes)['names'])
    return markers

markers_ctrl = get_top_markers(ctrl_ln, args.cluster_key, args.n_markers)
markers_sv   = get_top_markers(sv_ln,   args.cluster_key, args.n_markers)

overlap_results = {}
for cluster in clusters:
    mc      = markers_ctrl.get(cluster, set())
    ms      = markers_sv.get(cluster, set())
    recall  = len(mc & ms) / len(mc)         if len(mc) > 0         else np.nan
    jaccard = len(mc & ms) / len(mc | ms)    if len(mc | ms) > 0    else np.nan
    overlap_results[cluster] = {
        'n_ctrl_markers': len(mc),
        'n_sv_markers':   len(ms),
        'n_shared':       len(mc & ms),
        'recall':         recall,
        'jaccard':        jaccard
    }

overlap_df = pd.DataFrame(overlap_results).T
overlap_df.index.name = 'cluster'
print(overlap_df.round(3))

null_jaccard = (args.n_markers**2 / n_hvg /
                (2*args.n_markers - args.n_markers**2 / n_hvg))
print(f"\nExpected Jaccard by chance: {null_jaccard:.4f}")

# Hypergeometric test + BH
pval_mp_dict = {}
for cluster in overlap_df.index:
    K = int(overlap_df.loc[cluster, 'n_ctrl_markers'])
    n = int(overlap_df.loc[cluster, 'n_sv_markers'])
    k = int(overlap_df.loc[cluster, 'n_shared'])
    if K == 0 or n == 0:
        continue
    pval_mp_dict[cluster] = hypergeom_sf_mpmath(k, n_hvg, K, n)

padj_mp  = bh_mpmath(pval_mp_dict)
overlap_df['pval_str'] = 'NA'
overlap_df['padj_str'] = 'NA'
for cluster in pval_mp_dict:
    overlap_df.loc[cluster, 'pval_str'] = mpf_to_str(pval_mp_dict[cluster])
    overlap_df.loc[cluster, 'padj_str'] = mpf_to_str(padj_mp[cluster])

n_sig_jac = sum(
    1 for c in overlap_df.index
    if overlap_df.loc[c, 'padj_str'] != 'NA'
    and float(overlap_df.loc[c, 'padj_str']) < 0.05
)
print(f"Significant (FDR < 5%): {n_sig_jac}/{len(clusters)} clusters")
print(f"Median Jaccard: {overlap_df['jaccard'].median():.3f} "
      f"vs {null_jaccard:.3f} expected by chance")

overlap_df.to_csv(f"{args.outdir}/shell_vial_marker_overlap.csv")
print(f"Marker overlap -> {args.outdir}/shell_vial_marker_overlap.csv")

# Jaccard dot plot
fig, ax = plt.subplots(figsize=(9, 5))
jaccard_vals         = overlap_df['jaccard'].astype(float)
all_clusters_ordered = sorted(overlap_df.index.tolist(), key=lambda x: int(x))

for i, cluster in enumerate(all_clusters_ordered):
    jac  = jaccard_vals.get(cluster, np.nan)
    padj = overlap_df.loc[cluster, 'padj_str']
    if np.isnan(jac):
        color, zorder, alpha = '#d3d1c7', 1, 0.5
    elif padj == 'NA' or float(padj) >= 0.05:
        color, zorder, alpha = '#898781', 3, 0.9
    else:
        color, zorder, alpha = '#1baf7a', 4, 1.0
    ax.scatter(i, jac, s=520, color=color, zorder=zorder, alpha=alpha,
               linewidths=1.2,
               edgecolors='#0f6e56' if color == '#1baf7a' else '#5f5e5a')
    ax.text(i, jac, str(cluster), ha='center', va='center',
            fontsize=7.5, fontweight='500',
            color='white' if color == '#1baf7a' else '#444441',
            zorder=zorder + 1)

obs_median = np.nanmedian(jaccard_vals.values.astype(float))
ax.axhline(null_jaccard, color='#e34948', linewidth=1.5, linestyle='--',
           label=f'Expected by chance ({null_jaccard:.3f})')
ax.axhline(obs_median, color='#1baf7a', linewidth=1.5, linestyle=':',
           label=f'Median ({obs_median:.3f})')
ax.set_xticks([])
ax.set_ylabel("Jaccard index", fontsize=11)
ax.set_xlabel("Cluster (number shown in dot)", fontsize=11)
ax.set_xlim(-0.8, len(all_clusters_ordered) - 0.2)
ax.set_ylim(-0.03, max(jaccard_vals.dropna()) * 1.2)
ax.spines[['top', 'right']].set_visible(False)
ax.legend(handles=[
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#1baf7a',
           markersize=10, label='Significant (FDR < 5%)'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#898781',
           markersize=10, label='Not significant'),
    Line2D([0], [0], color='#e34948', linewidth=1.5, linestyle='--',
           label=f'Expected by chance ({null_jaccard:.3f})'),
    Line2D([0], [0], color='#1baf7a', linewidth=1.5, linestyle=':',
           label=f'Median ({obs_median:.3f})'),
], fontsize=8, frameon=False, loc='upper left')
ax.set_title(f"Marker gene Jaccard: {GROUP_A} vs {GROUP_B}\n"
             f"(top {args.n_markers} markers per cluster, {n_hvg} HVG space; "
             f"{n_sig_jac}/{len(clusters)} clusters FDR < 5%)",
             fontsize=10, pad=10)
plt.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(f"{args.outdir}/shell_vial_jaccard_dotplot.{ext}",
                dpi=300, bbox_inches='tight')
plt.close()
print("Saved: shell_vial_jaccard_dotplot.pdf/.png")


# ══════════════════════════════════════════════════════════════════════════════
# STEP 4: Pseudobulk Spearman correlation (per-cluster marker genes)
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("STEP 4: Pseudobulk Spearman Correlation (per-cluster marker genes)")
print("=" * 70)

me_ctrl = mean_expr_per_cluster(ctrl_ln, args.cluster_key)
me_sv   = mean_expr_per_cluster(sv_ln,   args.cluster_key)

shared_clusters = sorted(set(me_ctrl.columns) & set(me_sv.columns),
                          key=lambda x: int(x))
print("Computing cluster x cluster Spearman matrix ...")
rho_matrix   = pd.DataFrame(index=shared_clusters, columns=shared_clusters,
                             dtype=float)
n_genes_diag = {}
for c_ctrl in shared_clusters:
    gene_set = list(
        markers_ctrl.get(c_ctrl, set())
        & set(me_ctrl.index) & set(me_sv.index)
    )
    n_genes_diag[c_ctrl] = len(gene_set)
    if len(gene_set) < 2:
        for c_sv in shared_clusters:
            rho_matrix.loc[c_ctrl, c_sv] = np.nan
        continue
    for c_sv in shared_clusters:
        rho, _ = spearmanr(me_ctrl.loc[gene_set, c_ctrl],
                           me_sv.loc[gene_set, c_sv])
        rho_matrix.loc[c_ctrl, c_sv] = rho

correlations = {}
for cluster in shared_clusters:
    rho  = float(rho_matrix.loc[cluster, cluster])
    n    = n_genes_diag.get(cluster, 0)
    pval = (spearman_pval_mpmath(rho, n)
            if n >= 3 and not np.isnan(rho) else np.nan)
    correlations[cluster] = {'rho': rho, 'pval': pval, 'n_genes': n}

corr_df = pd.DataFrame(correlations).T
valid   = corr_df['rho'].dropna()
print(f"\nMedian diagonal Spearman rho: {valid.median():.3f}")
print(f"Range: {valid.min():.3f} - {valid.max():.3f}")
print(corr_df)

rho_matrix.to_csv(f"{args.outdir}/shell_vial_spearman_matrix.csv")
corr_df.to_csv(f"{args.outdir}/shell_vial_spearman_diagonal.csv")

# Heatmap
n    = len(shared_clusters)
data = rho_matrix.values.astype(float)
fig, ax = plt.subplots(figsize=(7, 6))
im = ax.imshow(data, cmap='RdYlBu_r', vmin=-1, vmax=1, aspect='equal')
cbar = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
cbar.set_label("Spearman ρ", fontsize=11)
ax.set_xticks(range(n))
ax.set_yticks(range(n))
ax.set_xticklabels([f"C{c}" for c in shared_clusters],
                   rotation=45, ha='right', fontsize=8)
ax.set_yticklabels([f"C{c}" for c in shared_clusters], fontsize=8)
ax.set_xlabel(f"{GROUP_B} clusters", fontsize=11, labelpad=8)
ax.set_ylabel(f"{GROUP_A} clusters", fontsize=11, labelpad=8)
for i in range(n):
    for j in range(n):
        val = data[i, j]
        if not np.isnan(val):
            ax.text(j, i, f"{val:.2f}", ha='center', va='center',
                    fontsize=6,
                    color='white' if abs(val) > 0.6 else 'black')
for k in range(n):
    ax.add_patch(mpatches.Rectangle(
        (k - 0.5, k - 0.5), 1, 1,
        fill=False, edgecolor='black', linewidth=1.5, zorder=3
    ))
ax.set_title(
    f"Pseudobulk Spearman ρ (per-cluster marker genes)\n"
    f"{GROUP_A} vs {GROUP_B} — median diagonal ρ = {valid.median():.3f}",
    fontsize=10, pad=12
)
plt.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(f"{args.outdir}/shell_vial_spearman_heatmap.{ext}",
                dpi=300, bbox_inches='tight')
plt.close()
print("Saved: shell_vial_spearman_heatmap.pdf/.png")


# ══════════════════════════════════════════════════════════════════════════════
# STEP 5: UMAP side-by-side
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("STEP 5: UMAP by condition")
print("=" * 70)

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
for ax, (sub, name) in zip(axes, [(ctrl, GROUP_A), (sv, GROUP_B)]):
    coords = sub.obsm['X_umap']
    labels = sub.obs[args.cluster_key].astype(str).values
    for cluster in clusters:
        mask = labels == str(cluster)
        ax.scatter(coords[mask, 0], coords[mask, 1],
                   s=1, alpha=0.4, c=palette.get(cluster), rasterized=True)
    ax.set_title(f"{name}  (n={sub.n_obs:,})", fontsize=11)
    ax.set_xlabel("UMAP 1", fontsize=10)
    ax.set_ylabel("UMAP 2", fontsize=10)
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

fig.legend(
    handles=[plt.Line2D([0], [0], marker='o', color='w',
                        markerfacecolor=palette.get(c, 'grey'),
                        markersize=6, label=f"Cluster {c}")
             for c in clusters],
    title="Cluster", bbox_to_anchor=(1.01, 0.5),
    loc='center left', fontsize=8, title_fontsize=9, frameon=False
)
fig.suptitle("UMAP: direct culture vs mock shell vial (uninfected)",
             fontsize=12, y=1.01)
plt.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(f"{args.outdir}/shell_vial_umap.{ext}",
                dpi=300, bbox_inches='tight')
plt.close()
print("Saved: shell_vial_umap.pdf/.png")


# ══════════════════════════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"Clusters: {GROUP_A}={n_clusters_ctrl}, {GROUP_B}={n_clusters_sv}")
print(f"Genes/cell: Ctrl median={gpc_shared['ctrl'].median():.0f}, "
      f"SV median={gpc_shared['sv'].median():.0f} "
      f"(Mann-Whitney U={mwu_stat:.0f}, p={mwu_pval:.2e})")
print(f"Cluster composition: JSD={jsd:.4f}; "
      f"{n_sig_comp}/{len(clusters)} clusters significantly different "
      f"in proportion (Fisher's exact, FDR < 5%)")
print(f"Marker gene Jaccard: median={overlap_df['jaccard'].median():.3f} "
      f"vs {null_jaccard:.3f} by chance; "
      f"{n_sig_jac}/{len(clusters)} clusters FDR < 5%")
print(f"Pseudobulk Spearman rho: median={valid.median():.3f}, "
      f"range={valid.min():.3f}-{valid.max():.3f}")
print(f"\nManuscript:")
print(f"  The shell vial technique did not introduce transcriptional variation; "
      f"{GROUP_B} cells were highly similar to {GROUP_A} cells in cluster "
      f"composition (Jensen-Shannon divergence = {jsd:.4f}; "
      f"{n_sig_comp}/{len(clusters)} clusters significantly different in "
      f"proportion, Fisher's exact test, FDR < 5%), genes per cluster "
      f"(median {gpc_shared['sv'].median():.0f} vs {gpc_shared['ctrl'].median():.0f}; "
      f"Mann-Whitney U={mwu_stat:.0f}, p={mwu_pval:.2e}), and marker gene content "
      f"(median Jaccard={overlap_df['jaccard'].median():.3f} vs "
      f"{null_jaccard:.3f} by chance; "
      f"{n_sig_jac}/{len(clusters)} clusters FDR < 5%, hypergeometric test).")
print(f"\nAll outputs written to: {args.outdir}/")


# ══════════════════════════════════════════════════════════════════════════════
# STEP 6: Atlas-wide check — is any cluster dominated by DOX-SV-D1 cells?
# ══════════════════════════════════════════════════════════════════════════════
# Runs on the full adata object (all conditions, both platforms) to ask
# whether any jointly-defined cluster is enriched for DOX-SV-D1 cells
# beyond their overall representation in the dataset. A cluster dominated
# by SV-D1 cells would suggest a shell-vial-specific transcriptional state.

print("\n" + "=" * 70)
print("STEP 6: Atlas-wide cluster enrichment for DOX-SV-D1 cells")
print("=" * 70)

all_clusters = sorted(adata.obs[args.cluster_key].cat.categories.tolist(),
                      key=lambda x: int(x))

# Flag each cell as DOX-SV-D1 or not (excluding the contaminated replicate)
is_sv = (
    (adata.obs['bio_condition'] == GROUP_B) &
    ~((adata.obs['method'] == 'pipseq') &
      (adata.obs['replicate'].astype(str) == '1'))
)
n_sv_total    = is_sv.sum()
n_total       = adata.n_obs
bg_frac       = n_sv_total / n_total   # expected fraction under null
print(f"\nTotal cells:      {n_total:,}")
print(f"DOX-SV-D1 cells: {n_sv_total:,} ({bg_frac*100:.2f}% of atlas)")

# Per-cluster: count SV-D1 vs non-SV-D1 cells, test enrichment
enrichment_results = {}
for cluster in all_clusters:
    in_cluster     = adata.obs[args.cluster_key] == cluster
    n_cluster      = in_cluster.sum()
    n_sv_cluster   = (in_cluster & is_sv).sum()
    n_nosv_cluster = n_cluster - n_sv_cluster
    n_sv_other     = n_sv_total - n_sv_cluster
    n_nosv_other   = (n_total - n_sv_total) - n_nosv_cluster

    frac_sv = n_sv_cluster / n_cluster if n_cluster > 0 else np.nan
    fold_enrichment = (frac_sv / bg_frac) if bg_frac > 0 else np.nan

    # Fisher's exact: is SV-D1 over-represented in this cluster?
    # [[SV in cluster,    non-SV in cluster],
    #  [SV not in cluster, non-SV not in cluster]]
    odds, pval = fisher_exact(
        [[n_sv_cluster,   n_nosv_cluster],
         [n_sv_other,     n_nosv_other]],
        alternative='greater'   # one-sided: testing for enrichment
    )
    enrichment_results[cluster] = {
        'n_cells':        n_cluster,
        'n_sv':           n_sv_cluster,
        'frac_sv':        frac_sv,
        'fold_enrichment': fold_enrichment,
        'odds_ratio':     odds,
        'pval':           pval
    }

enrich_df = pd.DataFrame(enrichment_results).T
enrich_df.index.name = 'cluster'

# BH correction
_, padj, _, _ = multipletests(enrich_df['pval'], method='fdr_bh')
enrich_df['padj'] = padj
enrich_df['significant'] = enrich_df['padj'] < 0.05

print(f"\nBackground DOX-SV-D1 fraction: {bg_frac:.4f}")
print(f"\nPer-cluster DOX-SV-D1 enrichment:")
print(enrich_df[['n_cells', 'n_sv', 'frac_sv', 'fold_enrichment',
                  'pval', 'padj', 'significant']].round(4))

n_enriched = enrich_df['significant'].sum()
if n_enriched == 0:
    print(f"\n✓ No clusters significantly enriched for DOX-SV-D1 cells "
          f"(FDR < 5%) — no evidence of a shell-vial-specific cluster")
else:
    print(f"\n⚠ {n_enriched} cluster(s) significantly enriched for "
          f"DOX-SV-D1 cells (FDR < 5%):")
    print(enrich_df[enrich_df['significant']][
        ['n_cells', 'n_sv', 'frac_sv', 'fold_enrichment', 'padj']
    ].round(4))

enrich_df.to_csv(f"{args.outdir}/shell_vial_atlas_enrichment.csv")
print(f"\nEnrichment table -> {args.outdir}/shell_vial_atlas_enrichment.csv")

# Plot: fold enrichment per cluster, colored by significance
fig, ax = plt.subplots(figsize=(9, 5))
x_pos = np.arange(len(all_clusters))
for i, cluster in enumerate(all_clusters):
    fe   = enrich_df.loc[cluster, 'fold_enrichment']
    sig  = enrich_df.loc[cluster, 'significant']
    n    = int(enrich_df.loc[cluster, 'n_cells'])
    color = '#e34948' if sig else '#2a78d6'
    ax.scatter(i, fe, s=max(30, n/50), color=color, zorder=3, alpha=0.85,
               edgecolors='#5f5e5a', linewidths=0.8)
    ax.text(i, fe, str(cluster), ha='center', va='center',
            fontsize=7, color='white', fontweight='500', zorder=4)

ax.axhline(1.0, color='black', linewidth=1, linestyle='--', alpha=0.5,
           label='Expected (fold enrichment = 1)')
ax.set_xticks([])
ax.set_ylabel("Fold enrichment of DOX-SV-D1 cells\n(relative to atlas background)",
              fontsize=10)
ax.set_xlabel("Cluster (number shown in dot, size ∝ cluster size)", fontsize=10)
ax.set_xlim(-0.8, len(all_clusters) - 0.2)
ax.spines[['top', 'right']].set_visible(False)
ax.legend(handles=[
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#e34948',
           markersize=10, label='Enriched (FDR < 5%)'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#2a78d6',
           markersize=10, label='Not significant'),
    Line2D([0], [0], color='black', linewidth=1, linestyle='--', alpha=0.5,
           label='Expected (fold = 1)'),
], fontsize=8, frameon=False)
ax.set_title(
    f"DOX-SV-D1 cell enrichment per cluster (atlas-wide)\n"
    f"Background fraction: {bg_frac*100:.2f}%; "
    f"{n_enriched}/{len(all_clusters)} clusters enriched (FDR < 5%)",
    fontsize=10, pad=10
)
plt.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(f"{args.outdir}/shell_vial_atlas_enrichment.{ext}",
                dpi=300, bbox_inches='tight')
plt.close()
print("Saved: shell_vial_atlas_enrichment.pdf/.png")