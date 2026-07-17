'''
Cell cycle annotation using Scanpy with Drosophila-specific markers
and cluster-cell cycle association analysis

# Reference (leiden column, use --ctrl-only if you want)
python scanpy_cellcycle_analysis.py --input integrated_reference.h5ad \
    --output results/cellcycle/ref --sample ref --run-scoring --save-output

# Query (leiden_ref column, do NOT use --ctrl-only)
python scanpy_cellcycle_analysis.py --input integrated_query.h5ad \
    --output results/cellcycle/query --sample query --run-scoring

# Combined, ctrl only
python scanpy_cellcycle_analysis.py --input integrated_combined.h5ad \
    --output results/cellcycle/combined_ctrl --sample combined_ctrl \
    --ctrl-only --run-scoring
    
'''
import scanpy as sc 
import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chi2_contingency, kruskal

# Drosophila cell cycle genes (FlyBase IDs)
# Mapping from gene symbols to FlyBase IDs
FLYBASE_CELL_CYCLE_GENES = {
    # S phase genes
    'Pcna': 'FBgn0005655',
    'RPA1': 'FBgn0015806', 
    'RPA2': 'FBgn0034898',
    'pol-alpha1': 'FBgn0011230',
    'DNApol-alpha60': 'FBgn0015278',
    'DNApol-delta': 'FBgn0019624',
    'RnrL': 'FBgn0020369',
    'RnrS': 'FBgn0261933',
    'Mcm2': 'FBgn0020651',
    'Mcm3': 'FBgn0020652',
    'Mcm5': 'FBgn0015929',
    'Mcm6': 'FBgn0032435',
    'Mcm7': 'FBgn0015308',
    'E2f1': 'FBgn0011766',
    'E2f2': 'FBgn0262656',
    'CycE': 'FBgn0010382',
    'Cdk2': 'FBgn0010314',
    'Dp': 'FBgn0000499',
    'Rbf': 'FBgn0015799',
    'Rbf2': 'FBgn0028396',
    'Orc1': 'FBgn0015270',
    'Orc2': 'FBgn0015714',
    'Orc6': 'FBgn0025926',
    'Rrp1': 'FBgn0003257',
    # G2/M genes
    'CycA': 'FBgn0010114',
    'CycB': 'FBgn0010113',
    'CycB3': 'FBgn0011577',
    'Cdk1': 'FBgn0004107',
    'stg': 'FBgn0003525',
    'polo': 'FBgn0003124',
    'aurA': 'FBgn0025564',
    'aurB': 'FBgn0025948',
    'Nek2': 'FBgn0027548',
    'Pbl': 'FBgn0005619',
    'Wee1': 'FBgn0011739',
    'myt': 'FBgn0002863',
    'BubR1': 'FBgn0024822',
    'Mad2': 'FBgn0002610',
    'Cdc20': 'FBgn0010309',
    'APC2': 'FBgn0261823',
    'APC10': 'FBgn0036449',
}

# Define S and G2M gene sets by FlyBase ID
S_GENES_FBGN = [
    'FBgn0005655',  # Pcna
    'FBgn0015806',  # RPA1
    'FBgn0034898',  # RPA2
    'FBgn0011230',  # pol-alpha1
    'FBgn0015278',  # DNApol-alpha60
    'FBgn0019624',  # DNApol-delta
    'FBgn0020369',  # RnrL
    'FBgn0261933',  # RnrS
    'FBgn0020651',  # Mcm2
    'FBgn0020652',  # Mcm3
    'FBgn0015929',  # Mcm5
    'FBgn0032435',  # Mcm6
    'FBgn0015308',  # Mcm7
    'FBgn0011766',  # E2f1
    'FBgn0262656',  # E2f2
    'FBgn0010382',  # CycE
    'FBgn0010314',  # Cdk2
    'FBgn0000499',  # Dp
    'FBgn0015799',  # Rbf
    'FBgn0028396',  # Rbf2
    'FBgn0015270',  # Orc1
    'FBgn0015714',  # Orc2
    'FBgn0025926',  # Orc6
    'FBgn0003257',  # Rrp1
]

G2M_GENES_FBGN = [
    'FBgn0010114',  # CycA
    'FBgn0010113',  # CycB
    'FBgn0011577',  # CycB3
    'FBgn0004107',  # Cdk1
    'FBgn0003525',  # stg
    'FBgn0003124',  # polo
    'FBgn0025564',  # aurA
    'FBgn0025948',  # aurB
    'FBgn0027548',  # Nek2
    'FBgn0005619',  # Pbl
    'FBgn0011739',  # Wee1
    'FBgn0002863',  # myt
    'FBgn0024822',  # BubR1
    'FBgn0002610',  # Mad2
    'FBgn0010309',  # Cdc20
    'FBgn0261823',  # APC2
    'FBgn0036449',  # APC10
]

# Create reverse mapping for reporting
FBGN_TO_SYMBOL = {v: k for k, v in FLYBASE_CELL_CYCLE_GENES.items()}


def _get_leiden_col(adata):
    """
    Return the leiden column name present in adata.obs.
    Prefers 'leiden_ref' (query/combined from integrate.py) over 'leiden' (reference).
    Raises a clear error if neither is found.
    """
    for col in ('leiden_ref', 'leiden'):
        if col in adata.obs.columns:
            return col
    raise KeyError(
        "No leiden column found in adata.obs. "
        "Expected 'leiden_ref' (query/combined) or 'leiden' (reference). "
        f"Available columns: {list(adata.obs.columns)}"
    )


def _ensure_lognorm(adata):
    """
    sc.tl.score_genes_cell_cycle expects log-normalised counts.
    If the matrix looks like raw counts (max > 20 is a rough heuristic),
    normalise and log-transform in-place.
    """
    import scipy.sparse
    X = adata.X
    mat_max = X.max() if not scipy.sparse.issparse(X) else X.max()
    if hasattr(mat_max, 'A1'):      # sparse scalar
        mat_max = float(mat_max.A1[0])
    else:
        mat_max = float(mat_max)

    if mat_max > 20:
        print("WARNING: Data appears to be raw counts (max value = "
              f"{mat_max:.1f}). Running normalize_total + log1p before scoring.")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    else:
        print(f"Data appears log-normalised (max value = {mat_max:.3f}). Skipping normalisation.")


def check_gene_names(adata):
    """Check what format the gene names are in."""
    print("\nChecking gene naming format...")
    sample_genes = list(adata.var_names[:10])
    print(f"Sample gene names: {sample_genes}")

    fbgn_count   = sum(1 for g in adata.var_names if str(g).startswith('FBgn'))
    symbol_count = sum(1 for g in adata.var_names if not str(g).startswith('FBgn'))

    print(f"\nGenes starting with 'FBgn': {fbgn_count}")
    print(f"Genes not starting with 'FBgn': {symbol_count}")

    if fbgn_count > symbol_count:
        print("-> Detected FlyBase ID format")
        return 'flybase'
    else:
        print("-> Detected gene symbol format")
        return 'symbol'


def score_cell_cycle_scanpy(adata, output_dir, sample_name, s_genes=None, g2m_genes=None):
    """Score cell cycle using Scanpy with Drosophila genes."""
    print("\n" + "=" * 60)
    print("SCANPY CELL CYCLE SCORING (DROSOPHILA)")
    print("=" * 60)

    check_gene_names(adata)

    # Fix 1: ensure data is log-normalised before scoring
    _ensure_lognorm(adata)

    if s_genes is None:
        s_genes = S_GENES_FBGN
    if g2m_genes is None:
        g2m_genes = G2M_GENES_FBGN

    s_genes_present   = [g for g in s_genes   if g in adata.var_names]
    g2m_genes_present = [g for g in g2m_genes if g in adata.var_names]

    print(f"\nS phase genes: {len(s_genes_present)}/{len(s_genes)} found")
    if s_genes_present:
        print("  Present (showing symbols): " +
              ', '.join([FBGN_TO_SYMBOL.get(g, g) for g in s_genes_present[:10]]) +
              (f"... (+{len(s_genes_present)-10} more)" if len(s_genes_present) > 10 else ""))
    else:
        print("  None found")

    print(f"\nG2/M phase genes: {len(g2m_genes_present)}/{len(g2m_genes)} found")
    if g2m_genes_present:
        print("  Present (showing symbols): " +
              ', '.join([FBGN_TO_SYMBOL.get(g, g) for g in g2m_genes_present[:10]]) +
              (f"... (+{len(g2m_genes_present)-10} more)" if len(g2m_genes_present) > 10 else ""))
    else:
        print("  None found")

    if len(s_genes_present) < 3 or len(g2m_genes_present) < 3:
        print("\nWARNING: Very few marker genes found. Results may not be reliable.")
        missing_s   = [FBGN_TO_SYMBOL.get(g, g) for g in s_genes   if g not in adata.var_names]
        missing_g2m = [FBGN_TO_SYMBOL.get(g, g) for g in g2m_genes if g not in adata.var_names]
        print(f"Missing S genes:   {', '.join(missing_s[:20])}")
        print(f"Missing G2/M genes: {', '.join(missing_g2m[:20])}")
        if len(s_genes_present) == 0 and len(g2m_genes_present) == 0:
            print("\nERROR: No cell cycle genes found! Check gene annotations.")
            return None

    print("\nScoring cell cycle phases...")
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes_present, g2m_genes=g2m_genes_present)

    print("\nCell cycle phase distribution:")
    phase_counts = adata.obs['phase'].value_counts()
    for phase in ['G1', 'S', 'G2M']:
        count = phase_counts.get(phase, 0)
        print(f"  {phase}: {count} cells ({count / adata.n_obs * 100:.1f}%)")

    create_scanpy_plots(adata, output_dir, sample_name, s_genes_present, g2m_genes_present)
    return adata


def create_scanpy_plots(adata, output_dir, sample_name, s_genes, g2m_genes):
    """Create diagnostic plots for Scanpy cell cycle scoring."""
    os.makedirs(output_dir, exist_ok=True)

    colors = {'G1': '#FF6B6B', 'S': '#4ECDC4', 'G2M': '#45B7D1'}

    # 1. Phase distribution overview
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    phase_counts = adata.obs['phase'].value_counts()
    phase_counts.plot(kind='bar', ax=axes[0],
                      color=[colors.get(p, 'gray') for p in phase_counts.index])
    axes[0].set_xlabel('Cell Cycle Phase')
    axes[0].set_ylabel('Number of Cells')
    axes[0].set_title('Cell Cycle Phase Distribution')
    axes[0].set_xticklabels(axes[0].get_xticklabels(), rotation=45)

    for phase, color in colors.items():
        mask = adata.obs['phase'] == phase
        axes[1].scatter(adata.obs.loc[mask, 'S_score'],
                        adata.obs.loc[mask, 'G2M_score'],
                        c=color, label=phase, alpha=0.5, s=10)
    axes[1].set_xlabel('S Score')
    axes[1].set_ylabel('G2M Score')
    axes[1].set_title('Cell Cycle Scores')
    axes[1].legend()
    axes[1].axhline(y=0, color='k', linestyle='--', alpha=0.3)
    axes[1].axvline(x=0, color='k', linestyle='--', alpha=0.3)

    axes[2].hist(adata.obs['S_score'],   bins=50, alpha=0.5, label='S score',   color=colors['S'])
    axes[2].hist(adata.obs['G2M_score'], bins=50, alpha=0.5, label='G2M score', color=colors['G2M'])
    axes[2].set_xlabel('Score')
    axes[2].set_ylabel('Number of Cells')
    axes[2].set_title('Score Distributions')
    axes[2].legend()
    axes[2].axvline(x=0, color='k', linestyle='--', alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{sample_name}_cellcycle_overview.pdf'),
                dpi=300, bbox_inches='tight')
    plt.close()

    # 2. Marker gene expression heatmap
    import scipy.sparse
    all_markers = list(set(s_genes + g2m_genes))
    if all_markers:
        adata_sorted = adata[adata.obs.sort_values('phase').index, :]
        X = adata_sorted[:, all_markers].X
        if scipy.sparse.issparse(X):
            X = X.toarray()
        marker_expr = pd.DataFrame(X, index=adata_sorted.obs_names, columns=all_markers).T
        marker_expr_norm = marker_expr.apply(lambda x: (x - x.mean()) / (x.std() + 1e-10), axis=1)
        gene_labels = [FBGN_TO_SYMBOL.get(g, g) for g in all_markers]

        fig, ax = plt.subplots(figsize=(12, max(8, len(all_markers) * 0.3)))
        sns.heatmap(marker_expr_norm, cmap='RdBu_r', center=0,
                    cbar_kws={'label': 'Z-score'},
                    xticklabels=False, yticklabels=gene_labels,
                    vmin=-2, vmax=2, ax=ax)
        ax.set_xlabel('Cells (sorted by phase)')
        ax.set_ylabel('Cell Cycle Marker Genes')
        ax.set_title(f'Cell cycle marker expression ({len(all_markers)} genes)\n{sample_name}')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'{sample_name}_marker_heatmap.pdf'),
                    dpi=300, bbox_inches='tight')
        plt.close()

    print(f"\nCell cycle plots saved to {output_dir}")


def analyze_cluster_cellcycle_association(adata, fig_dir, sample):
    """Test and visualize association between Leiden clusters and cell cycle."""
    print("\n" + "=" * 60)
    print("CLUSTER - CELL CYCLE ASSOCIATION ANALYSIS")
    print("=" * 60)

    # Fix 1: handle both 'leiden' (reference) and 'leiden_ref' (query/combined)
    try:
        leiden_col = _get_leiden_col(adata)
    except KeyError as e:
        print(f"ERROR: {e}")
        return None

    print(f"Using leiden column: '{leiden_col}'")

    if 'phase' not in adata.obs.columns:
        print("ERROR: No 'phase' found in adata.obs. Run with --run-scoring first.")
        return None

    os.makedirs(fig_dir, exist_ok=True)
    # Set figdir before any sc.pl calls
    sc.settings.figdir = fig_dir

    clusters = sorted(adata.obs[leiden_col].unique())
    cmap = plt.cm.get_cmap('tab20')
    leiden_colors = [cmap(i % 20) for i in range(len(clusters))]

    # ── 1. Chi-square test ────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("1. CHI-SQUARE TEST: Cluster vs Cell Cycle Stage")
    print("=" * 60)

    contingency = pd.crosstab(adata.obs[leiden_col], adata.obs['phase'])
    chi2, p_value, dof, expected = chi2_contingency(contingency)

    print(f"chi2 = {chi2:.2f}")
    print(f"degrees of freedom = {dof}")
    print(f"p-value = {p_value:.2e}")
    print(f"\nConclusion: Clusters are "
          f"{'SIGNIFICANTLY' if p_value < 0.05 else 'NOT significantly'} "
          f"associated with cell cycle stage")

    n = contingency.sum().sum()
    cramers_v = np.sqrt(chi2 / (n * (min(contingency.shape) - 1)))
    print(f"Cramer's V = {cramers_v:.3f}")
    if cramers_v < 0.1:       effect = "negligible"
    elif cramers_v < 0.3:     effect = "weak"
    elif cramers_v < 0.5:     effect = "moderate"
    else:                     effect = "strong"
    print(f"Effect size: {effect}")

    # ── 2. Heatmap ────────────────────────────────────────────────────────────
    contingency_norm = contingency.div(contingency.sum(axis=1), axis=0) * 100
    print("\nCell cycle stage distribution by cluster (%):")
    print(contingency_norm.round(1))

    fig, ax = plt.subplots(figsize=(12, 8))
    sns.heatmap(contingency_norm, annot=True, fmt='.1f', cmap='YlOrRd', ax=ax,
                cbar_kws={'label': '% of cells in cluster'})
    ax.set_xlabel('Cell Cycle Phase')
    ax.set_ylabel(f'Leiden Cluster ({leiden_col})')
    ax.set_title(f"Cell cycle phase distribution by cluster\n"
                 f"chi2={chi2:.2f}, p={p_value:.2e}, Cramer's V={cramers_v:.3f}")
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'heatmap_cluster_cellcycle_{sample}.pdf'))
    plt.close()

    # ── 3. Stacked bar ────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(14, 7))
    contingency_norm.plot(kind='bar', stacked=True, ax=ax, width=0.8)
    ax.set_xlabel(f'Leiden Cluster ({leiden_col})', fontsize=12)
    ax.set_ylabel('Percentage of cells', fontsize=12)
    ax.set_title('Cell cycle phase composition by cluster', fontsize=14)
    ax.legend(title='Cell Cycle Phase', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks(rotation=0)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'barplot_cluster_cellcycle_{sample}.pdf'))
    plt.close()

    # ── 4. Enrichment summary ─────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("2. CELL CYCLE ENRICHED CLUSTERS")
    print("=" * 60)

    dominant_phase  = contingency_norm.idxmax(axis=1)
    max_percentage  = contingency_norm.max(axis=1)

    print(f"\n{'Cluster':<10} {'Dominant Phase':<15} {'Percentage':<12} Status")
    print("-" * 60)
    for cluster in clusters:
        phase = dominant_phase[cluster]
        pct   = max_percentage[cluster]
        enrichment = ("STRONGLY ENRICHED" if pct > 50
                      else "ENRICHED" if pct > 40
                      else "Mixed")
        print(f"{cluster:<10} {phase:<15} {pct:>6.1f}%      {enrichment}")

    # ── 5. Cell cycle scores by cluster ──────────────────────────────────────
    print("\n" + "=" * 60)
    print("3. CELL CYCLE SCORES BY CLUSTER")
    print("=" * 60)

    score_by_cluster = adata.obs.groupby(leiden_col)[['S_score', 'G2M_score']].agg(['mean', 'std'])
    print(score_by_cluster.round(3))

    groups_s   = [adata.obs[adata.obs[leiden_col] == c]['S_score'].dropna().values   for c in clusters]
    groups_g2m = [adata.obs[adata.obs[leiden_col] == c]['G2M_score'].dropna().values for c in clusters]
    h_stat_s,   p_value_s   = kruskal(*groups_s)
    h_stat_g2m, p_value_g2m = kruskal(*groups_g2m)

    print(f"\nKruskal-Wallis S_score:   H={h_stat_s:.2f},   p={p_value_s:.2e}")
    print(f"Kruskal-Wallis G2M_score: H={h_stat_g2m:.2f}, p={p_value_g2m:.2e}")

    # Violin plots
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    sc.pl.violin(adata, 'S_score',   groupby=leiden_col, ax=axes[0], show=False, rotation=0)
    axes[0].set_title(f'S Score by Cluster\nKruskal-Wallis H={h_stat_s:.2f}, p={p_value_s:.2e}')
    sc.pl.violin(adata, 'G2M_score', groupby=leiden_col, ax=axes[1], show=False, rotation=0)
    axes[1].set_title(f'G2M Score by Cluster\nKruskal-Wallis H={h_stat_g2m:.2f}, p={p_value_g2m:.2e}')
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'violin_scores_by_cluster_{sample}.pdf'))
    plt.close()

    # Bar plots for mean scores
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    for ax, score_key, title in zip(
        axes,
        ['S_score', 'G2M_score'],
        ['Mean S Score by Cluster', 'Mean G2M Score by Cluster'],
    ):
        means = score_by_cluster[score_key]['mean']
        stds  = score_by_cluster[score_key]['std']
        ax.bar(range(len(means)), means, yerr=stds,
               color=leiden_colors, alpha=0.7, capsize=5)
        ax.set_xlabel(f'Leiden Cluster ({leiden_col})')
        ax.set_ylabel(f'Mean {score_key}')
        ax.set_title(title)
        ax.set_xticks(range(len(means)))
        ax.set_xticklabels(means.index)
        ax.axhline(y=0, color='k', linestyle='--', alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'barplot_scores_by_cluster_{sample}.pdf'))
    plt.close()

    # ── 6. UMAPs ──────────────────────────────────────────────────────────────
    if 'X_umap' in adata.obsm:
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))
        sc.pl.umap(adata, color=leiden_col, ax=axes[0], show=False,
                   title=f'Leiden Clusters ({leiden_col})', frameon=False, legend_loc='on data')
        sc.pl.umap(adata, color='phase', ax=axes[1], show=False,
                   title='Cell Cycle Phase', frameon=False)
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f'umap_cluster_vs_cellcycle_{sample}.pdf'))
        plt.close()

        sc.pl.umap(adata, color=leiden_col,
                   save=f'_{sample}_leiden.pdf',
                   title=f'Leiden Clusters ({leiden_col})', legend_loc='on data')
        sc.pl.umap(adata, color='phase',
                   save=f'_{sample}_cellcycle_phase.pdf',
                   title='Cell Cycle Phase')
        sc.pl.umap(adata, color=['S_score', 'G2M_score'],
                   save=f'_{sample}_cellcycle_scores.pdf',
                   cmap='viridis')

    # ── 7. Summary table ──────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("4. COMPREHENSIVE SUMMARY TABLE")
    print("=" * 60)

    summary_df = pd.DataFrame({
        'Cluster':          clusters,
        'N_Cells':          [contingency.loc[c].sum()                             for c in clusters],
        'Dominant_Phase':   [dominant_phase[c]                                    for c in clusters],
        'Phase_Percentage': [max_percentage[c]                                    for c in clusters],
        'Mean_S_Score':     [score_by_cluster.loc[c, ('S_score',   'mean')]       for c in clusters],
        'Mean_G2M_Score':   [score_by_cluster.loc[c, ('G2M_score', 'mean')]       for c in clusters],
    })
    print(summary_df.to_string(index=False))

    summary_df.to_csv(os.path.join(fig_dir, f'cluster_cellcycle_summary_{sample}.csv'),   index=False)
    contingency.to_csv(os.path.join(fig_dir,      f'contingency_table_counts_{sample}.csv'))
    contingency_norm.to_csv(os.path.join(fig_dir, f'contingency_table_percentages_{sample}.csv'))

    stats_df = pd.DataFrame({
        'Test':       ['Chi-square', "Cramer's V",
                       'Kruskal-Wallis (S_score)', 'Kruskal-Wallis (G2M_score)'],
        'Statistic':  [f'{chi2:.6e}',     f'{cramers_v:.6e}',
                       f'{h_stat_s:.6e}', f'{h_stat_g2m:.6e}'],
        'P-value':    [f'{p_value:.6e}',  'NA',
                       f'{p_value_s:.6e}', f'{p_value_g2m:.6e}'],
        'Interpretation': [
            'Significant association'     if p_value   < 0.05 else 'No significant association',
            f'{effect.capitalize()} effect size',
            'S scores differ'             if p_value_s < 0.05 else 'No difference in S scores',
            'G2M scores differ'           if p_value_g2m < 0.05 else 'No difference in G2M scores',
        ],
    })
    stats_df.to_csv(os.path.join(fig_dir, f'statistical_tests_{sample}.csv'), index=False)

    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE")
    print("=" * 60)
    print(f"Output directory: {fig_dir}")

    return {
        'chi2':           chi2,
        'chi2_pvalue':    p_value,
        'cramers_v':      cramers_v,
        'kw_s_pvalue':    p_value_s,
        'kw_g2m_pvalue':  p_value_g2m,
        'summary':        summary_df,
    }


def main():
    parser = argparse.ArgumentParser(
        description='Cell cycle annotation using Scanpy with Drosophila genes',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Reference object (uninfected controls, has 'leiden' column)
  python scanpy_cellcycle_analysis.py \\
      --input results/integrated/integrated_reference.h5ad \\
      --output results/cellcycle/reference \\
      --sample reference --run-scoring --save-output

  # Query object (infected timepoints, has 'leiden_ref' column)
  python scanpy_cellcycle_analysis.py \\
      --input results/integrated/integrated_query.h5ad \\
      --output results/cellcycle/query \\
      --sample query --run-scoring

  # Combined object, ctrl cells only
  python scanpy_cellcycle_analysis.py \\
      --input results/integrated/integrated_combined.h5ad \\
      --output results/cellcycle/combined_ctrl \\
      --sample combined_ctrl --ctrl-only --run-scoring
        '''
    )

    parser.add_argument('--input',  '-i', required=True,
                        help='Path to h5ad file (_reference, _query, or _combined)')
    parser.add_argument('--output', '-o', default='cellcycle_analysis',
                        help='Output directory (default: cellcycle_analysis)')
    parser.add_argument('--sample', '-s', default='sample',
                        help='Sample name for output files (default: sample)')
    parser.add_argument('--run-scoring', action='store_true',
                        help='Run Scanpy cell cycle scoring (skip if already done)')
    parser.add_argument('--save-output', action='store_true',
                        help='Save updated h5ad with cell cycle annotations')
    # Fix 2: ctrl-only is now an opt-in flag, not hardcoded
    parser.add_argument('--ctrl-only', action='store_true',
                        help='Restrict analysis to treatment==Ctrl cells only '
                             '(do NOT use with _query.h5ad — query cells are '
                             'new infections and will be entirely dropped)')

    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    print(f"Loading data from {args.input}...")
    adata = sc.read_h5ad(args.input)

    # Fix 2: only filter to Ctrl if explicitly requested
    if args.ctrl_only:
        if 'treatment' not in adata.obs.columns:
            print("WARNING: --ctrl-only requested but 'treatment' column not found. Skipping filter.")
        else:
            before = adata.n_obs
            adata = adata[adata.obs['treatment'] == 'Ctrl'].copy()
            print(f"Filtered to Ctrl cells: {before} -> {adata.n_obs} cells")
            if adata.n_obs == 0:
                raise ValueError(
                    "No cells remaining after --ctrl-only filter. "
                    "If running on _query.h5ad, do not use --ctrl-only: "
                    "query cells are new infection timepoints with no 'Ctrl' treatment label."
                )

    print(f"\nLoaded AnnData: {adata.n_obs} cells, {adata.n_vars} genes")
    print(f"obs columns: {list(adata.obs.columns)}")

    has_scoring = all(c in adata.obs.columns for c in ['phase', 'S_score', 'G2M_score'])

    if args.run_scoring or not has_scoring:
        result = score_cell_cycle_scanpy(adata, args.output, args.sample)
        if result is None:
            print("\nERROR: Cell cycle scoring failed. Exiting.")
            return
        adata = result
        if args.save_output:
            out_path = args.input.replace('.h5ad', '_with_cellcycle.h5ad')
            print(f"\nSaving updated h5ad to {out_path}")
            adata.write(out_path)
    else:
        print("\nCell cycle scoring already present, skipping scoring step.")
        print(f"Phase distribution: {adata.obs['phase'].value_counts().to_dict()}")

    results = analyze_cluster_cellcycle_association(adata, args.output, args.sample)

    if results:
        print(f"\n{'=' * 60}")
        print("KEY FINDINGS")
        print(f"{'=' * 60}")
        print(f"Chi-square: chi2={results['chi2']:.2f}, p={results['chi2_pvalue']:.2e}")
        print(f"Cramer's V: {results['cramers_v']:.3f}")
        print(f"Kruskal-Wallis S score:   p={results['kw_s_pvalue']:.2e}")
        print(f"Kruskal-Wallis G2M score: p={results['kw_g2m_pvalue']:.2e}")


if __name__ == "__main__":
    main()