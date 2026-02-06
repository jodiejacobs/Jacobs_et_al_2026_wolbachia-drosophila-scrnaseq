'''
Cell cycle annotation using Scanpy with Drosophila-specific markers
and cluster-cell cycle association analysis
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


def check_gene_names(adata):
    """
    Check what format the gene names are in
    """
    print("\nChecking gene naming format...")
    sample_genes = list(adata.var_names[:10])
    print(f"Sample gene names: {sample_genes}")
    
    # Check if genes are FlyBase IDs
    fbgn_count = sum(1 for g in adata.var_names if str(g).startswith('FBgn'))
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
    """
    Score cell cycle using Scanpy with Drosophila genes
    """
    print("\n" + "="*60)
    print("SCANPY CELL CYCLE SCORING (DROSOPHILA)")
    print("="*60)
    
    # Check gene naming format
    gene_format = check_gene_names(adata)
    
    # Use FlyBase IDs
    if s_genes is None:
        s_genes = S_GENES_FBGN
    if g2m_genes is None:
        g2m_genes = G2M_GENES_FBGN
    
    # Check which genes are present
    s_genes_present = [g for g in s_genes if g in adata.var_names]
    g2m_genes_present = [g for g in g2m_genes if g in adata.var_names]
    
    print(f"\nS phase genes: {len(s_genes_present)}/{len(s_genes)} found")
    if s_genes_present:
        print(f"  Present (showing symbols): {', '.join([FBGN_TO_SYMBOL.get(g, g) for g in s_genes_present[:10]])}" + 
              (f"... (+{len(s_genes_present)-10} more)" if len(s_genes_present) > 10 else ""))
    else:
        print("  None found")
    
    print(f"\nG2/M phase genes: {len(g2m_genes_present)}/{len(g2m_genes)} found")
    if g2m_genes_present:
        print(f"  Present (showing symbols): {', '.join([FBGN_TO_SYMBOL.get(g, g) for g in g2m_genes_present[:10]])}" + 
              (f"... (+{len(g2m_genes_present)-10} more)" if len(g2m_genes_present) > 10 else ""))
    else:
        print("  None found")
    
    if len(s_genes_present) < 3 or len(g2m_genes_present) < 3:
        print("\nWARNING: Very few marker genes found. Results may not be reliable.")
        print("\nMissing S genes (showing symbols):")
        missing_s = [FBGN_TO_SYMBOL.get(g, g) for g in s_genes if g not in adata.var_names]
        print(f"  {', '.join(missing_s[:20])}")
        print("\nMissing G2/M genes (showing symbols):")
        missing_g2m = [FBGN_TO_SYMBOL.get(g, g) for g in g2m_genes if g not in adata.var_names]
        print(f"  {', '.join(missing_g2m[:20])}")
        
        if len(s_genes_present) == 0 and len(g2m_genes_present) == 0:
            print("\nERROR: No cell cycle genes found!")
            print("Please check your gene annotations.")
            return None
    
    # Score genes
    print("\nScoring cell cycle phases...")
    sc.tl.score_genes_cell_cycle(
        adata, 
        s_genes=s_genes_present, 
        g2m_genes=g2m_genes_present
    )
    
    # Print distribution
    print("\nCell cycle phase distribution:")
    phase_counts = adata.obs['phase'].value_counts()
    for phase in ['G1', 'S', 'G2M']:
        count = phase_counts.get(phase, 0)
        pct = (count / adata.n_obs) * 100
        print(f"  {phase}: {count} cells ({pct:.1f}%)")
    
    # Create plots
    create_scanpy_plots(adata, output_dir, sample_name, s_genes_present, g2m_genes_present)
    
    return adata


def create_scanpy_plots(adata, output_dir, sample_name, s_genes, g2m_genes):
    """
    Create diagnostic plots for Scanpy cell cycle scoring
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Phase distribution
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    # Bar plot
    ax = axes[0]
    phase_counts = adata.obs['phase'].value_counts()
    colors = {'G1': '#FF6B6B', 'S': '#4ECDC4', 'G2M': '#45B7D1'}
    phase_counts.plot(kind='bar', ax=ax, color=[colors.get(p, 'gray') for p in phase_counts.index])
    ax.set_xlabel('Cell Cycle Phase')
    ax.set_ylabel('Number of Cells')
    ax.set_title('Cell Cycle Phase Distribution')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
    
    # S score vs G2M score scatter
    ax = axes[1]
    for phase, color in colors.items():
        mask = adata.obs['phase'] == phase
        ax.scatter(adata.obs.loc[mask, 'S_score'], 
                  adata.obs.loc[mask, 'G2M_score'],
                  c=color, label=phase, alpha=0.5, s=10)
    ax.set_xlabel('S Score')
    ax.set_ylabel('G2M Score')
    ax.set_title('Cell Cycle Scores')
    ax.legend()
    ax.axhline(y=0, color='k', linestyle='--', alpha=0.3)
    ax.axvline(x=0, color='k', linestyle='--', alpha=0.3)
    
    # Score distributions
    ax = axes[2]
    ax.hist(adata.obs['S_score'], bins=50, alpha=0.5, label='S score', color=colors['S'])
    ax.hist(adata.obs['G2M_score'], bins=50, alpha=0.5, label='G2M score', color=colors['G2M'])
    ax.set_xlabel('Score')
    ax.set_ylabel('Number of Cells')
    ax.set_title('Score Distributions')
    ax.legend()
    ax.axvline(x=0, color='k', linestyle='--', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{sample_name}_cellcycle_overview.pdf'), 
                dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Marker gene expression heatmap
    all_markers = list(set(s_genes + g2m_genes))
    if len(all_markers) > 0:
        # Sort cells by phase
        adata_sorted = adata[adata.obs.sort_values('phase').index, :]
        
        fig, ax = plt.subplots(figsize=(12, max(8, len(all_markers) * 0.3)))
        
        # Get expression data
        marker_expr = pd.DataFrame(
            adata_sorted[:, all_markers].X.toarray() if hasattr(adata_sorted.X, 'toarray') else adata_sorted[:, all_markers].X,
            index=adata_sorted.obs_names,
            columns=all_markers
        ).T
        
        # Normalize by row (z-score)
        marker_expr_norm = marker_expr.apply(
            lambda x: (x - x.mean()) / (x.std() + 1e-10), 
            axis=1
        )
        
        # Use gene symbols for y-axis labels
        gene_labels = [FBGN_TO_SYMBOL.get(g, g) for g in all_markers]
        
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
    """
    Test and visualize association between Leiden clusters and cell cycle
    """
    print("\n" + "="*60)
    print("CLUSTER - CELL CYCLE ASSOCIATION ANALYSIS")
    print("="*60)
    
    # Check for required columns
    if 'leiden' not in adata.obs.columns:
        print("ERROR: No 'leiden' clustering found in adata.obs")
        return None
    
    if 'phase' not in adata.obs.columns:
        print("ERROR: No 'phase' found in adata.obs")
        print("Run cell cycle scoring first with --run-scoring flag")
        return None
    
    sc.settings.figdir = fig_dir
    
    # Get leiden colors
    leiden_colors = []
    clusters = sorted(adata.obs['leiden'].unique())
    cmap = plt.cm.get_cmap('tab20')
    for i, cluster in enumerate(clusters):
        leiden_colors.append(cmap(i % 20))
    
    # 1. Chi-square test
    print("\n" + "="*60)
    print("1. CHI-SQUARE TEST: Cluster vs Cell Cycle Stage")
    print("="*60)
    
    contingency = pd.crosstab(adata.obs['leiden'], adata.obs['phase'])
    chi2, p_value, dof, expected = chi2_contingency(contingency)
    
    print(f"χ² = {chi2:.2f}")
    print(f"degrees of freedom = {dof}")
    print(f"p-value = {p_value:.2e}")
    print(f"\nConclusion: Clusters are {'SIGNIFICANTLY' if p_value < 0.05 else 'NOT significantly'} associated with cell cycle stage")
    
    # Cramér's V
    n = contingency.sum().sum()
    cramers_v = np.sqrt(chi2 / (n * (min(contingency.shape) - 1)))
    print(f"Cramér's V = {cramers_v:.3f}")
    
    if cramers_v < 0.1:
        effect = "negligible"
    elif cramers_v < 0.3:
        effect = "weak"
    elif cramers_v < 0.5:
        effect = "moderate"
    else:
        effect = "strong"
    print(f"Effect size: {effect}")
    
    # 2. Heatmap
    contingency_norm = contingency.div(contingency.sum(axis=1), axis=0) * 100
    
    print("\nCell cycle stage distribution by cluster (%):")
    print(contingency_norm.round(1))
    
    fig, ax = plt.subplots(figsize=(12, 8))
    sns.heatmap(contingency_norm, annot=True, fmt='.1f', cmap='YlOrRd', 
                ax=ax, cbar_kws={'label': '% of cells in cluster'})
    ax.set_xlabel('Cell Cycle Phase')
    ax.set_ylabel('Leiden Cluster')
    ax.set_title(f'Cell cycle phase distribution by cluster\nχ² = {chi2:.2f}, p = {p_value:.2e}, Cramér\'s V = {cramers_v:.3f}')
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'heatmap_cluster_cellcycle_{sample}.pdf'))
    plt.close()
    
    # 3. Stacked bar
    fig, ax = plt.subplots(figsize=(14, 7))
    contingency_norm.plot(kind='bar', stacked=True, ax=ax, width=0.8)
    ax.set_xlabel('Leiden Cluster', fontsize=12)
    ax.set_ylabel('Percentage of cells', fontsize=12)
    ax.set_title(f'Cell cycle phase composition by cluster', fontsize=14)
    ax.legend(title='Cell Cycle Phase', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks(rotation=0)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'barplot_cluster_cellcycle_{sample}.pdf'))
    plt.close()
    
    # 4. Enrichment analysis
    print("\n" + "="*60)
    print("2. CELL CYCLE ENRICHED CLUSTERS")
    print("="*60)
    
    dominant_phase = contingency_norm.idxmax(axis=1)
    max_percentage = contingency_norm.max(axis=1)
    
    print("\nDominant cell cycle phase per cluster:")
    print(f"{'Cluster':<10} {'Dominant Phase':<15} {'Percentage':<12} {'Status'}")
    print("-" * 60)
    for cluster in clusters:
        phase = dominant_phase[cluster]
        pct = max_percentage[cluster]
        if pct > 50:
            enrichment = "STRONGLY ENRICHED"
        elif pct > 40:
            enrichment = "ENRICHED"
        else:
            enrichment = "Mixed"
        print(f"{cluster:<10} {phase:<15} {pct:>6.1f}%      {enrichment}")
    
    # 5. Cell cycle scores by cluster
    print("\n" + "="*60)
    print("3. CELL CYCLE SCORES BY CLUSTER")
    print("="*60)
    
    score_by_cluster = adata.obs.groupby('leiden')[['S_score', 'G2M_score']].agg(['mean', 'std'])
    print(score_by_cluster.round(3))
    
    # Kruskal-Wallis test for S_score
    groups_s = [adata.obs[adata.obs['leiden'] == cluster]['S_score'].dropna().values 
              for cluster in clusters]
    h_stat_s, p_value_s = kruskal(*groups_s)
    
    # Kruskal-Wallis test for G2M_score
    groups_g2m = [adata.obs[adata.obs['leiden'] == cluster]['G2M_score'].dropna().values 
                  for cluster in clusters]
    h_stat_g2m, p_value_g2m = kruskal(*groups_g2m)
    
    print(f"\nKruskal-Wallis test for S_score: H = {h_stat_s:.2f}, p = {p_value_s:.2e}")
    print(f"Kruskal-Wallis test for G2M_score: H = {h_stat_g2m:.2f}, p = {p_value_g2m:.2e}")
    
    # Violin plots for scores
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    sc.pl.violin(adata, 'S_score', groupby='leiden', 
                ax=axes[0], show=False, rotation=0)
    axes[0].set_title(f'S Score by Cluster\nKruskal-Wallis H = {h_stat_s:.2f}, p = {p_value_s:.2e}')
    
    sc.pl.violin(adata, 'G2M_score', groupby='leiden', 
                ax=axes[1], show=False, rotation=0)
    axes[1].set_title(f'G2M Score by Cluster\nKruskal-Wallis H = {h_stat_g2m:.2f}, p = {p_value_g2m:.2e}')
    
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'violin_scores_by_cluster_{sample}.pdf'))
    plt.close()
    
    # Bar plots for mean scores
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # S score
    s_means = score_by_cluster['S_score']['mean']
    s_stds = score_by_cluster['S_score']['std']
    axes[0].bar(range(len(s_means)), s_means, yerr=s_stds, 
               color=leiden_colors, alpha=0.7, capsize=5)
    axes[0].set_xlabel('Leiden Cluster')
    axes[0].set_ylabel('Mean S Score')
    axes[0].set_title('Mean S Score by Cluster')
    axes[0].set_xticks(range(len(s_means)))
    axes[0].set_xticklabels(s_means.index)
    axes[0].axhline(y=0, color='k', linestyle='--', alpha=0.3)
    
    # G2M score
    g2m_means = score_by_cluster['G2M_score']['mean']
    g2m_stds = score_by_cluster['G2M_score']['std']
    axes[1].bar(range(len(g2m_means)), g2m_means, yerr=g2m_stds,
               color=leiden_colors, alpha=0.7, capsize=5)
    axes[1].set_xlabel('Leiden Cluster')
    axes[1].set_ylabel('Mean G2M Score')
    axes[1].set_title('Mean G2M Score by Cluster')
    axes[1].set_xticks(range(len(g2m_means)))
    axes[1].set_xticklabels(g2m_means.index)
    axes[1].axhline(y=0, color='k', linestyle='--', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'barplot_scores_by_cluster_{sample}.pdf'))
    plt.close()
    
    # 6. UMAPs
    if 'X_umap' in adata.obsm:
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))
        sc.pl.umap(adata, color='leiden', ax=axes[0], show=False, 
                   title='Leiden Clusters', frameon=False, legend_loc='on data')
        sc.pl.umap(adata, color='phase', ax=axes[1], show=False,
                   title='Cell Cycle Phase', frameon=False)
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f'umap_cluster_vs_cellcycle_{sample}.pdf'))
        plt.close()
        
        # Individual UMAPs
        sc.pl.umap(adata, color='leiden', save=f'_{sample}_leiden.pdf',
                   title='Leiden Clusters', legend_loc='on data')
        sc.pl.umap(adata, color='phase', save=f'_{sample}_cellcycle_phase.pdf',
                   title='Cell Cycle Phase')
        sc.pl.umap(adata, color=['S_score', 'G2M_score'], 
                  save=f'_{sample}_cellcycle_scores.pdf',
                  cmap='viridis')
    
    # 7. Summary table
    print("\n" + "="*60)
    print("4. COMPREHENSIVE SUMMARY TABLE")
    print("="*60)
    
    summary_data = {
        'Cluster': clusters,
        'N_Cells': [contingency.loc[c].sum() for c in clusters],
        'Dominant_Phase': [dominant_phase[c] for c in clusters],
        'Phase_Percentage': [max_percentage[c] for c in clusters],
        'Mean_S_Score': [score_by_cluster.loc[c, ('S_score', 'mean')] for c in clusters],
        'Mean_G2M_Score': [score_by_cluster.loc[c, ('G2M_score', 'mean')] for c in clusters],
    }
    
    summary_df = pd.DataFrame(summary_data)
    print(summary_df.to_string(index=False))
    
    # Save outputs
    summary_df.to_csv(os.path.join(fig_dir, f'cluster_cellcycle_summary_{sample}.csv'), index=False)
    contingency.to_csv(os.path.join(fig_dir, f'contingency_table_counts_{sample}.csv'))
    contingency_norm.to_csv(os.path.join(fig_dir, f'contingency_table_percentages_{sample}.csv'))
    
    # Statistical results - FIXED with proper scientific notation
    stats_results = {
        'Test': ['Chi-square', 'Cramers V', 'Kruskal-Wallis (S_score)', 'Kruskal-Wallis (G2M_score)'],
        'Statistic': [f'{chi2:.6e}', f'{cramers_v:.6e}', f'{h_stat_s:.6e}', f'{h_stat_g2m:.6e}'],
        'P-value': [f'{p_value:.6e}', 'NA', f'{p_value_s:.6e}', f'{p_value_g2m:.6e}'],
        'Interpretation': [
            'Significant association' if p_value < 0.05 else 'No significant association',
            f'{effect.capitalize()} effect size',
            'S scores differ across clusters' if p_value_s < 0.05 else 'No difference in S scores',
            'G2M scores differ across clusters' if p_value_g2m < 0.05 else 'No difference in G2M scores'
        ]
    }
    
    stats_df = pd.DataFrame(stats_results)
    stats_df.to_csv(os.path.join(fig_dir, f'statistical_tests_{sample}.csv'), index=False)
    
    print("\n" + "="*60)
    print("ANALYSIS COMPLETE")
    print("="*60)
    print(f"Output directory: {fig_dir}")
    
    return {
        'chi2': chi2,
        'chi2_pvalue': p_value,
        'cramers_v': cramers_v,
        'kw_s_pvalue': p_value_s,
        'kw_g2m_pvalue': p_value_g2m,
        'summary': summary_df
    }

def main():
    parser = argparse.ArgumentParser(
        description='Cell cycle annotation using Scanpy with Drosophila genes',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Run Scanpy cell cycle scoring + analysis
  python scanpy_cellcycle_analysis.py \\
      --input integrated.h5ad \\
      --output cellcycle_analysis \\
      --sample all_conditions \\
      --run-scoring
  
  # Just run analysis (scoring already done)
  python scanpy_cellcycle_analysis.py \\
      --input integrated.h5ad \\
      --output cellcycle_analysis \\
      --sample all_conditions
  
  # Save updated h5ad
  python scanpy_cellcycle_analysis.py \\
      --input integrated.h5ad \\
      --output cellcycle_analysis \\
      --sample all_conditions \\
      --run-scoring \\
      --save-output
        '''
    )
    
    parser.add_argument('--input', '-i', required=True, type=str,
                        help='Path to h5ad file')
    parser.add_argument('--output', '-o', type=str, default='cellcycle_analysis',
                        help='Output directory (default: cellcycle_analysis)')
    parser.add_argument('--sample', '-s', type=str, default='sample',
                        help='Sample name for output files (default: sample)')
    parser.add_argument('--run-scoring', action='store_true',
                        help='Run Scanpy cell cycle scoring (skip if already done)')
    parser.add_argument('--save-output', action='store_true',
                        help='Save updated h5ad with cell cycle annotations')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    
    # Load data
    print(f"Loading data from {args.input}...")
    adata = sc.read_h5ad(args.input)

    # Filter for Control samples if needed 
    adata = adata[adata.obs['treatment'] == 'Ctrl', :]
    print(adata.obs)
    
    print(f"\nLoaded AnnData object:")
    print(f"  Cells: {adata.n_obs}")
    print(f"  Genes: {adata.n_vars}")
    print(f"  Observations: {list(adata.obs.columns)}")
    
    # Run Scanpy scoring if requested or if not present
    has_scoring = all(col in adata.obs.columns for col in ['phase', 'S_score', 'G2M_score'])
    
    if args.run_scoring or not has_scoring:
        result = score_cell_cycle_scanpy(adata, args.output, args.sample)
        if result is None:
            print("\nERROR: Cell cycle scoring failed. Exiting.")
            return
        adata = result
        
        # Save if requested
        if args.save_output:
            output_path = args.input.replace('.h5ad', '_with_cellcycle.h5ad')
            print(f"\nSaving updated h5ad to {output_path}")
            adata.write(output_path)
    else:
        print("\nCell cycle scoring already present, skipping...")
        print(f"Phase distribution: {adata.obs['phase'].value_counts().to_dict()}")
    
    # Run cluster-cell cycle association analysis
    results = analyze_cluster_cellcycle_association(adata, args.output, args.sample)
    
    if results:
        print(f"\n{'='*60}")
        print("KEY FINDINGS")
        print(f"{'='*60}")
        print(f"Chi-square test: χ² = {results['chi2']:.2f}, p = {results['chi2_pvalue']:.2e}")
        print(f"Effect size (Cramér's V): {results['cramers_v']:.3f}")
        print(f"Kruskal-Wallis (S score): p = {results['kw_s_pvalue']:.2e}")
        print(f"Kruskal-Wallis (G2M score): p = {results['kw_g2m_pvalue']:.2e}")


if __name__ == "__main__":
    main()