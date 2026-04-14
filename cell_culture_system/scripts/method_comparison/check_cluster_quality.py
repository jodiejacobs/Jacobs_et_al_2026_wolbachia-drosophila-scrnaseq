# cluster_qc.py
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations

def flag_bad_clusters(
    adata,
    cluster_key='leiden',
    min_cell_frac=0.01,
    min_median_genes=200,
    corr_threshold=0.97,
    marker_score_threshold=2.0,
    de_method='wilcoxon',
    output_prefix='cluster_qc'
):
    """
    Flags clusters that may need merging based on:
      - Size (too few cells)
      - Complexity (too few genes)
      - Transcriptional similarity (high pairwise correlation)
      - Weak DE markers

    Returns a DataFrame of per-cluster flags and saves a Jaccard heatmap.
    """

    obs = adata.obs.copy()
    clusters = adata.obs[cluster_key].cat.categories.tolist()
    n_total = len(adata)

    # ── 1. Basic stats ─────────────────────────────────────────────────────────
    if 'n_genes' not in obs.columns:
        obs['n_genes'] = np.asarray((adata.X > 0).sum(axis=1)).flatten()

    stats = obs.groupby(cluster_key, observed=False)['n_genes'].agg(
        n_cells='count',
        median_genes='median'
    ).reset_index()
    stats['cell_frac'] = stats['n_cells'] / n_total
    stats['flag_small']      = stats['cell_frac'] < min_cell_frac
    stats['flag_low_genes']  = stats['median_genes'] < min_median_genes

    # ── 2. DE marker strength ──────────────────────────────────────────────────
    print("Running DE (Wilcoxon)...")
    sc.tl.rank_genes_groups(adata, groupby=cluster_key, method=de_method, pts=True)
    marker_df = sc.get.rank_genes_groups_df(adata, group=None)
    top_scores = marker_df.groupby('group')['scores'].max().reset_index()
    top_scores.columns = [cluster_key, 'max_marker_score']
    stats = stats.merge(top_scores, on=cluster_key, how='left')
    stats['flag_weak_markers'] = stats['max_marker_score'] < marker_score_threshold

    # ── 3. Pairwise Jaccard on top marker genes ────────────────────────────────
    top_genes = {}
    for cl in clusters:
        df_cl = marker_df[marker_df['group'] == cl].nlargest(50, 'scores')
        top_genes[cl] = set(df_cl['names'].tolist())

    jaccard_mat = pd.DataFrame(index=clusters, columns=clusters, dtype=float)
    for c1, c2 in combinations(clusters, 2):
        g1, g2 = top_genes[c1], top_genes[c2]
        j = len(g1 & g2) / len(g1 | g2) if len(g1 | g2) > 0 else 0.0
        jaccard_mat.loc[c1, c2] = j
        jaccard_mat.loc[c2, c1] = j
    for c in clusters:
        jaccard_mat.loc[c, c] = 1.0
    jaccard_mat = jaccard_mat.astype(float)

    high_jaccard_pairs = []
    for c1, c2 in combinations(clusters, 2):
        j = jaccard_mat.loc[c1, c2]
        if j > 0.3:
            high_jaccard_pairs.append((c1, c2, round(j, 3)))

    # ── 4. Mean expression correlation ────────────────────────────────────────
    print("Computing cluster mean expression...")
    cluster_means = pd.DataFrame(index=adata.var_names)
    for cl in clusters:
        mask = adata.obs[cluster_key] == cl
        cluster_means[cl] = np.asarray(adata[mask].X.mean(axis=0)).flatten()

    corr_mat = cluster_means.corr()

    high_corr_pairs = []
    for c1, c2 in combinations(clusters, 2):
        r = corr_mat.loc[c1, c2]
        if r > corr_threshold:
            high_corr_pairs.append((c1, c2, round(r, 4)))

    # ── 5. Summary flags ───────────────────────────────────────────────────────
    flagged_from_pairs = set()
    for c1, c2, _ in high_jaccard_pairs:
        flagged_from_pairs.update([c1, c2])
    for c1, c2, _ in high_corr_pairs:
        flagged_from_pairs.update([c1, c2])

    stats['flag_similar_to_another'] = stats[cluster_key].isin(flagged_from_pairs)
    stats['any_flag'] = stats[['flag_small', 'flag_low_genes',
                                'flag_weak_markers', 'flag_similar_to_another']].any(axis=1)

    # ── 6. Plots ───────────────────────────────────────────────────────────────

    # Jaccard + flag heatmaps
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    mask_diag = np.eye(len(clusters), dtype=bool)
    sns.heatmap(
        jaccard_mat, ax=axes[0],
        cmap='YlOrRd', vmin=0, vmax=1,
        annot=len(clusters) <= 20,
        fmt='.2f',
        linewidths=0.5,
        mask=mask_diag,
        square=True,
        cbar_kws={'label': 'Jaccard similarity'}
    )
    axes[0].set_title('Jaccard similarity (top 50 marker genes)')
    axes[0].set_xlabel('Cluster')
    axes[0].set_ylabel('Cluster')

    flag_cols = ['flag_small', 'flag_low_genes', 'flag_weak_markers', 'flag_similar_to_another']
    flag_display = stats.set_index(cluster_key)[flag_cols].astype(int)
    sns.heatmap(
        flag_display, ax=axes[1],
        cmap=['#f0f0f0', '#e74c3c'],
        vmin=0, vmax=1,
        linewidths=0.5,
        cbar=False,
        annot=True,
        fmt='d'
    )
    axes[1].set_title('Cluster flags (1 = flagged)')
    axes[1].set_xticklabels(
        ['Too small', 'Low genes', 'Weak markers', 'Similar to\nanother'],
        rotation=30, ha='right'
    )

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_flags.png', dpi=150, bbox_inches='tight')
    print(f"Saved: {output_prefix}_flags.png")
    plt.close()

    # UMAP panels with per-panel vmax
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    sc.pl.umap(adata, color='leiden',        ax=axes[0, 0], show=False)
    sc.pl.umap(adata, color='n_genes',       ax=axes[0, 1], show=False, vmax=1000)
    sc.pl.umap(adata, color='total_counts',  ax=axes[1, 0], show=False)
    sc.pl.umap(adata, color='pct_counts_mt', ax=axes[1, 1], show=False)
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_flag_umap.png', dpi=150, bbox_inches='tight')
    print(f"Saved: {output_prefix}_flag_umap.png")
    plt.close()

    # ── 7. Print report ────────────────────────────────────────────────────────
    print("\n── Cluster QC report ──────────────────────────────")
    print(stats[[cluster_key, 'n_cells', 'cell_frac', 'median_genes',
                 'max_marker_score', 'any_flag']].to_string(index=False))

    if high_jaccard_pairs:
        print("\nHigh Jaccard pairs (marker gene overlap >30%):")
        for c1, c2, j in sorted(high_jaccard_pairs, key=lambda x: -x[2]):
            print(f"  Clusters {c1} <-> {c2}  J={j}")

    if high_corr_pairs:
        print(f"\nHigh correlation pairs (r>{corr_threshold}):")
        for c1, c2, r in sorted(high_corr_pairs, key=lambda x: -x[2]):
            print(f"  Clusters {c1} <-> {c2}  r={r}")

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    # Histogram of n_genes across all cells
    axes[0].hist(adata.obs['n_genes_by_counts'], bins=200, color='steelblue', alpha=0.7)
    axes[0].axvline(500,  color='orange', linestyle='--', label='500')
    axes[0].axvline(1000, color='red',    linestyle='--', label='1000')
    axes[0].set_xlabel('n_genes per cell')
    axes[0].set_ylabel('n cells')
    axes[0].set_title('Gene count distribution')
    axes[0].legend()

    # Log scale to see the low end more clearly
    axes[1].hist(np.log1p(adata.obs['n_genes_by_counts']), bins=200, color='steelblue', alpha=0.7)
    axes[1].set_xlabel('log(n_genes + 1)')
    axes[1].set_title('Log gene count distribution')

    plt.tight_layout()
    plt.savefig('ngenes_distribution.png', dpi=150)

# ── 8. Flagged cluster inspection plots ───────────────────────────────────
    flagged_clusters = stats[stats['flag_similar_to_another']][cluster_key].tolist()
    small_clusters   = stats[stats['flag_small']][cluster_key].tolist()
    inspect_clusters = list(set(flagged_clusters + small_clusters))

    if inspect_clusters:
        # UMAP highlighting similar/small clusters vs rest
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        sc.pl.umap(adata, color=cluster_key, ax=axes[0], show=False,
                   groups=inspect_clusters, na_color='#e0e0e0',
                   title='Flagged clusters (similar/small)')

        sc.pl.umap(adata, color=cluster_key, ax=axes[1], show=False,
                   title='All clusters (labeled)')

        plt.tight_layout()
        plt.savefig(f'{output_prefix}_flagged_umap.png', dpi=150, bbox_inches='tight')
        print(f"Saved: {output_prefix}_flagged_umap.png")
        plt.close()

        # Dotplot of top 10 markers for flagged clusters
        sc.pl.rank_genes_groups_dotplot(
            adata,
            groups=inspect_clusters,
            n_genes=10,
            show=False
        )

        plt.savefig(f'{output_prefix}_flagged_dotplot.png', dpi=150, bbox_inches='tight')
        print(f"Saved: {output_prefix}_flagged_dotplot.png")
        plt.close()
        
        # Ensure group labels are strings to match rank_genes_groups keys
        inspect_clusters_str = [str(c) for c in inspect_clusters]
        sc.tl.dendrogram(adata, groupby=cluster_key)

        ax = sc.pl.rank_genes_groups_dotplot(
            adata,
            groups=inspect_clusters_str,
            n_genes=10,
            show=False,
            return_fig=True
        )
        ax.savefig(f'{output_prefix}_flagged_dotplot.png', dpi=150, bbox_inches='tight')
        plt.close()

    return stats, jaccard_mat, high_jaccard_pairs, high_corr_pairs


# ── Usage ──────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    adata = sc.read_h5ad('/private/groups/russelllab/jodie/scRNAseq/Jacobs_et_al_2026_wolbachia-drosophila-scrnaseq/cell_culture_system/results/integrated/integrated.h5ad')

    stats, jaccard_mat, jaccard_pairs, corr_pairs = flag_bad_clusters(
        adata,
        cluster_key='leiden',
        min_cell_frac=0.01,
        min_median_genes=200,
        corr_threshold=0.97,
        marker_score_threshold=2.0,
        output_prefix='results/qc/cluster_qc'
    )

    # Optional: merge flagged clusters
    # merge_map = {'5': '2', '8': '2'}  # based on output above
    # adata.obs['leiden_merged'] = adata.obs['leiden'].map(merge_map).fillna(adata.obs['leiden'])