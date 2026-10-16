'''
Cell cycle annotation using Scanpy with Drosophila-specific markers
and analysis of Wolbachia titer associations with clusters and cell cycle
'''
import scanpy as sc 
import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chi2_contingency, kruskal, spearmanr, pearsonr, mannwhitneyu
from scipy.stats import f_oneway

# [Keep all the gene definitions from your original script]
FLYBASE_CELL_CYCLE_GENES = {
    'Pcna': 'FBgn0005655', 'RPA1': 'FBgn0015806', 'RPA2': 'FBgn0034898',
    'pol-alpha1': 'FBgn0011230', 'DNApol-alpha60': 'FBgn0015278',
    'DNApol-delta': 'FBgn0019624', 'RnrL': 'FBgn0020369', 'RnrS': 'FBgn0261933',
    'Mcm2': 'FBgn0020651', 'Mcm3': 'FBgn0020652', 'Mcm5': 'FBgn0015929',
    'Mcm6': 'FBgn0032435', 'Mcm7': 'FBgn0015308', 'E2f1': 'FBgn0011766',
    'E2f2': 'FBgn0262656', 'CycE': 'FBgn0010382', 'Cdk2': 'FBgn0010314',
    'Dp': 'FBgn0000499', 'Rbf': 'FBgn0015799', 'Rbf2': 'FBgn0028396',
    'Orc1': 'FBgn0015270', 'Orc2': 'FBgn0015714', 'Orc6': 'FBgn0025926',
    'Rrp1': 'FBgn0003257', 'CycA': 'FBgn0010114', 'CycB': 'FBgn0010113',
    'CycB3': 'FBgn0011577', 'Cdk1': 'FBgn0004107', 'stg': 'FBgn0003525',
    'polo': 'FBgn0003124', 'aurA': 'FBgn0025564', 'aurB': 'FBgn0025948',
    'Nek2': 'FBgn0027548', 'Pbl': 'FBgn0005619', 'Wee1': 'FBgn0011739',
    'myt': 'FBgn0002863', 'BubR1': 'FBgn0024822', 'Mad2': 'FBgn0002610',
    'Cdc20': 'FBgn0010309', 'APC2': 'FBgn0261823', 'APC10': 'FBgn0036449',
}

S_GENES_FBGN = ['FBgn0005655', 'FBgn0015806', 'FBgn0034898', 'FBgn0011230',
                'FBgn0015278', 'FBgn0019624', 'FBgn0020369', 'FBgn0261933',
                'FBgn0020651', 'FBgn0020652', 'FBgn0015929', 'FBgn0032435',
                'FBgn0015308', 'FBgn0011766', 'FBgn0262656', 'FBgn0010382',
                'FBgn0010314', 'FBgn0000499', 'FBgn0015799', 'FBgn0028396',
                'FBgn0015270', 'FBgn0015714', 'FBgn0025926', 'FBgn0003257']

G2M_GENES_FBGN = ['FBgn0010114', 'FBgn0010113', 'FBgn0011577', 'FBgn0004107',
                  'FBgn0003525', 'FBgn0003124', 'FBgn0025564', 'FBgn0025948',
                  'FBgn0027548', 'FBgn0005619', 'FBgn0011739', 'FBgn0002863',
                  'FBgn0024822', 'FBgn0002610', 'FBgn0010309', 'FBgn0261823',
                  'FBgn0036449']

FBGN_TO_SYMBOL = {v: k for k, v in FLYBASE_CELL_CYCLE_GENES.items()}


def check_gene_names(adata):
    """Check what format the gene names are in"""
    print("\nChecking gene naming format...")
    sample_genes = list(adata.var_names[:10])
    print(f"Sample gene names: {sample_genes}")
    
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


def detect_titer_column(adata):
    """
    Detect which column contains Wolbachia titer information
    """
    print("\n" + "="*60)
    print("DETECTING WOLBACHIA TITER COLUMN")
    print("="*60)
    
    # Common titer column names
    possible_names = [
        'wolbachia_titer', 'titer', 'Wolbachia_titer', 'Titer',
        'wolbachia_reads', 'wolbachia_counts', 'wol_titer',
        'wMel_titer', 'wRi_titer', 'wMelCS_titer'
    ]
    
    # Check for exact matches
    for name in possible_names:
        if name in adata.obs.columns:
            print(f"Found titer column: '{name}'")
            print(f"  Range: {adata.obs[name].min():.3f} to {adata.obs[name].max():.3f}")
            print(f"  Mean: {adata.obs[name].mean():.3f}")
            print(f"  Median: {adata.obs[name].median():.3f}")
            return name
    
    # Check for partial matches
    print("\nNo exact match found. Checking for partial matches...")
    for col in adata.obs.columns:
        if any(term in col.lower() for term in ['titer', 'wolbachia', 'wol', 'wmel', 'wri']):
            # Check if it's numeric
            if pd.api.types.is_numeric_dtype(adata.obs[col]):
                print(f"Potential titer column: '{col}'")
                print(f"  Range: {adata.obs[col].min():.3f} to {adata.obs[col].max():.3f}")
                print(f"  Mean: {adata.obs[col].mean():.3f}")
                response = input(f"Use '{col}' as titer column? (y/n): ")
                if response.lower() == 'y':
                    return col
    
    print("\nERROR: No titer column found!")
    print("Available columns:")
    for col in adata.obs.columns:
        if pd.api.types.is_numeric_dtype(adata.obs[col]):
            print(f"  - {col}")
    
    return None


def score_cell_cycle_scanpy(adata, output_dir, sample_name, s_genes=None, g2m_genes=None):
    """
    Score cell cycle using Scanpy with Drosophila genes
    """
    print("\n" + "="*60)
    print("SCANPY CELL CYCLE SCORING (DROSOPHILA)")
    print("="*60)
    
    gene_format = check_gene_names(adata)
    
    if s_genes is None:
        s_genes = S_GENES_FBGN
    if g2m_genes is None:
        g2m_genes = G2M_GENES_FBGN
    
    s_genes_present = [g for g in s_genes if g in adata.var_names]
    g2m_genes_present = [g for g in g2m_genes if g in adata.var_names]
    
    print(f"\nS phase genes: {len(s_genes_present)}/{len(s_genes)} found")
    if s_genes_present:
        print(f"  Present: {', '.join([FBGN_TO_SYMBOL.get(g, g) for g in s_genes_present[:10]])}" + 
              (f"... (+{len(s_genes_present)-10} more)" if len(s_genes_present) > 10 else ""))
    
    print(f"\nG2/M phase genes: {len(g2m_genes_present)}/{len(g2m_genes)} found")
    if g2m_genes_present:
        print(f"  Present: {', '.join([FBGN_TO_SYMBOL.get(g, g) for g in g2m_genes_present[:10]])}" + 
              (f"... (+{len(g2m_genes_present)-10} more)" if len(g2m_genes_present) > 10 else ""))
    
    if len(s_genes_present) == 0 and len(g2m_genes_present) == 0:
        print("\nERROR: No cell cycle genes found!")
        return None
    
    print("\nScoring cell cycle phases...")
    sc.tl.score_genes_cell_cycle(
        adata, 
        s_genes=s_genes_present, 
        g2m_genes=g2m_genes_present
    )
    
    print("\nCell cycle phase distribution:")
    phase_counts = adata.obs['phase'].value_counts()
    for phase in ['G1', 'S', 'G2M']:
        count = phase_counts.get(phase, 0)
        pct = (count / adata.n_obs) * 100
        print(f"  {phase}: {count} cells ({pct:.1f}%)")
    
    # Add continuous pseudotime position
    print("\nCalculating continuous cell cycle pseudotime...")
    adata = calculate_cellcycle_pseudotime(adata)
    
    create_scanpy_plots(adata, output_dir, sample_name, s_genes_present, g2m_genes_present)
    
    return adata


def calculate_cellcycle_pseudotime(adata):
    """
    Calculate a continuous pseudotime measure (0-1) through the cell cycle
    based on S and G2M scores and phase assignments
    """
    s_score = adata.obs['S_score'].values
    g2m_score = adata.obs['G2M_score'].values
    phase = adata.obs['phase'].values
    
    pseudotime = np.zeros(len(phase))
    
    for i, (p, s, g) in enumerate(zip(phase, s_score, g2m_score)):
        if p == 'G1':
            # G1: 0.0 to 0.33
            base = 0.165
            jitter = (s + g) * 0.05
            pseudotime[i] = np.clip(base + jitter, 0.0, 0.33)
        elif p == 'S':
            # S: 0.33 to 0.66
            base = 0.495
            s_in_s = s_score[phase == 'S']
            if len(s_in_s) > 0:
                s_norm = (s - s_in_s.min()) / (s_in_s.max() - s_in_s.min() + 1e-10)
                jitter = (s_norm - 0.5) * 0.25
            else:
                jitter = 0
            pseudotime[i] = np.clip(base + jitter, 0.33, 0.66)
        else:  # G2M
            # G2M: 0.66 to 1.0
            base = 0.83
            g_in_g2m = g2m_score[phase == 'G2M']
            if len(g_in_g2m) > 0:
                g_norm = (g - g_in_g2m.min()) / (g_in_g2m.max() - g_in_g2m.min() + 1e-10)
                jitter = (g_norm - 0.5) * 0.25
            else:
                jitter = 0
            pseudotime[i] = np.clip(base + jitter, 0.66, 1.0)
    
    adata.obs['cellcycle_pseudotime'] = pseudotime
    
    print(f"Pseudotime range: {pseudotime.min():.3f} to {pseudotime.max():.3f}")
    print(f"Mean pseudotime: {pseudotime.mean():.3f}")
    
    return adata


def create_scanpy_plots(adata, output_dir, sample_name, s_genes, g2m_genes):
    """
    Create diagnostic plots for Scanpy cell cycle scoring
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Phase distribution
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # Bar plot
    ax = axes[0, 0]
    phase_counts = adata.obs['phase'].value_counts()
    colors = {'G1': '#FF6B6B', 'S': '#4ECDC4', 'G2M': '#45B7D1'}
    phase_counts.plot(kind='bar', ax=ax, color=[colors.get(p, 'gray') for p in phase_counts.index])
    ax.set_xlabel('Cell Cycle Phase')
    ax.set_ylabel('Number of Cells')
    ax.set_title('Cell Cycle Phase Distribution')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
    
    # S score vs G2M score scatter
    ax = axes[0, 1]
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
    
    # Pseudotime distribution
    ax = axes[1, 0]
    ax.hist(adata.obs['cellcycle_pseudotime'], bins=50, alpha=0.7, edgecolor='black', color='purple')
    ax.set_xlabel('Cell Cycle Pseudotime (0-1)')
    ax.set_ylabel('Number of Cells')
    ax.set_title('Continuous Cell Cycle Pseudotime Distribution')
    ax.axvline(x=0.33, color='red', linestyle='--', alpha=0.5, label='G1/S')
    ax.axvline(x=0.66, color='green', linestyle='--', alpha=0.5, label='S/G2M')
    ax.legend()
    
    # Circular plot
    ax = axes[1, 1]
    ax = plt.subplot(2, 2, 4, projection='polar')
    theta = adata.obs['cellcycle_pseudotime'].values * 2 * np.pi
    phase_colors = [colors.get(p, 'gray') for p in adata.obs['phase']]
    ax.scatter(theta, np.ones(len(theta)), c=phase_colors, alpha=0.3, s=1)
    ax.set_ylim(0, 1.5)
    ax.set_yticks([])
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_title('Cell Cycle Position (Circular)', pad=20)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{sample_name}_cellcycle_overview.pdf'), 
                dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"\nCell cycle plots saved to {output_dir}")


def analyze_titer_cellcycle_association(adata, fig_dir, sample, titer_col):
    """
    Analyze association between Wolbachia titer and cell cycle
    """
    print("\n" + "="*60)
    print("WOLBACHIA TITER - CELL CYCLE ASSOCIATION ANALYSIS")
    print("="*60)
    
    if titer_col not in adata.obs.columns:
        print(f"ERROR: Titer column '{titer_col}' not found")
        return None
    
    if 'phase' not in adata.obs.columns:
        print("ERROR: Cell cycle phase annotations not found")
        return None
    
    if 'cellcycle_pseudotime' not in adata.obs.columns:
        print("ERROR: Cell cycle pseudotime not found")
        return None
    
    os.makedirs(fig_dir, exist_ok=True)
    
    # Remove cells with missing titer values
    valid_mask = ~adata.obs[titer_col].isna()
    n_removed = (~valid_mask).sum()
    if n_removed > 0:
        print(f"\nRemoving {n_removed} cells with missing titer values")
        adata_valid = adata[valid_mask, :].copy()
    else:
        adata_valid = adata
    
    print(f"Analyzing {adata_valid.n_obs} cells with valid titer data")
    
    # 1. Test association with discrete phases
    print("\n" + "="*60)
    print("1. TITER vs CELL CYCLE PHASE")
    print("="*60)
    
    # Kruskal-Wallis test (non-parametric ANOVA)
    phases = ['G1', 'S', 'G2M']
    groups = [adata_valid.obs[adata_valid.obs['phase'] == p][titer_col].values for p in phases]
    h_stat, p_value_kw = kruskal(*groups)
    
    print(f"Kruskal-Wallis test: H = {h_stat:.2f}, p = {p_value_kw:.2e}")
    print(f"Conclusion: Titer {'SIGNIFICANTLY differs' if p_value_kw < 0.05 else 'does NOT significantly differ'} across cell cycle phases")
    
    # Summary statistics by phase
    print("\nTiter statistics by phase:")
    titer_stats = adata_valid.obs.groupby('phase')[titer_col].agg(['mean', 'median', 'std', 'count'])
    print(titer_stats)
    
    # 2. Test correlation with continuous pseudotime
    print("\n" + "="*60)
    print("2. TITER vs CELL CYCLE PSEUDOTIME")
    print("="*60)
    
    pseudotime = adata_valid.obs['cellcycle_pseudotime'].values
    titer = adata_valid.obs[titer_col].values
    
    # Pearson correlation
    pearson_r, pearson_p = pearsonr(pseudotime, titer)
    print(f"Pearson correlation: r = {pearson_r:.4f}, p = {pearson_p:.2e}")
    
    # Spearman correlation (non-parametric)
    spearman_r, spearman_p = spearmanr(pseudotime, titer)
    print(f"Spearman correlation: ρ = {spearman_r:.4f}, p = {spearman_p:.2e}")
    
    if abs(spearman_r) < 0.1:
        effect = "negligible"
    elif abs(spearman_r) < 0.3:
        effect = "weak"
    elif abs(spearman_r) < 0.5:
        effect = "moderate"
    else:
        effect = "strong"
    
    print(f"Effect size: {effect}")
    print(f"Conclusion: Titer {'IS' if spearman_p < 0.05 else 'is NOT'} significantly correlated with cell cycle pseudotime")
    
    # 3. Create visualizations
    create_titer_cellcycle_plots(adata_valid, fig_dir, sample, titer_col, 
                                  h_stat, p_value_kw, spearman_r, spearman_p)
    
    # 4. Test by binned pseudotime (tertiles/quartiles)
    print("\n" + "="*60)
    print("3. TITER vs BINNED PSEUDOTIME")
    print("="*60)
    
    # Divide pseudotime into tertiles
    adata_valid.obs['pseudotime_tertile'] = pd.qcut(
        adata_valid.obs['cellcycle_pseudotime'], 
        q=3, 
        labels=['Early', 'Mid', 'Late']
    )
    
    # Test difference across tertiles
    tertile_groups = [
        adata_valid.obs[adata_valid.obs['pseudotime_tertile'] == t][titer_col].values 
        for t in ['Early', 'Mid', 'Late']
    ]
    h_stat_tertile, p_value_tertile = kruskal(*tertile_groups)
    
    print(f"Kruskal-Wallis (tertiles): H = {h_stat_tertile:.2f}, p = {p_value_tertile:.2e}")
    
    tertile_stats = adata_valid.obs.groupby('pseudotime_tertile')[titer_col].agg(['mean', 'median', 'std', 'count'])
    print("\nTiter by pseudotime tertile:")
    print(tertile_stats)
    
    # Save results
    results = {
        'kw_phase_H': h_stat,
        'kw_phase_p': p_value_kw,
        'pearson_r': pearson_r,
        'pearson_p': pearson_p,
        'spearman_r': spearman_r,
        'spearman_p': spearman_p,
        'kw_tertile_H': h_stat_tertile,
        'kw_tertile_p': p_value_tertile,
        'titer_stats_by_phase': titer_stats,
        'titer_stats_by_tertile': tertile_stats
    }
    
    # Save statistical results
    stats_results = {
        'Test': [
            'Kruskal-Wallis (Phase)', 
            'Kruskal-Wallis (Tertiles)',
            'Pearson correlation', 
            'Spearman correlation'
        ],
        'Statistic': [
            f'{h_stat:.6e}',
            f'{h_stat_tertile:.6e}',
            f'{pearson_r:.6e}',
            f'{spearman_r:.6e}'
        ],
        'P-value': [
            f'{p_value_kw:.6e}',
            f'{p_value_tertile:.6e}',
            f'{pearson_p:.6e}',
            f'{spearman_p:.6e}'
        ],
        'Interpretation': [
            'Titer differs across phases' if p_value_kw < 0.05 else 'No difference across phases',
            'Titer differs across tertiles' if p_value_tertile < 0.05 else 'No difference across tertiles',
            f'Pearson r = {pearson_r:.3f}',
            f'Spearman ρ = {spearman_r:.3f}, {effect} correlation'
        ]
    }
    
    stats_df = pd.DataFrame(stats_results)
    stats_df.to_csv(os.path.join(fig_dir, f'titer_cellcycle_stats_{sample}.csv'), index=False)
    
    # Save summary stats
    titer_stats.to_csv(os.path.join(fig_dir, f'titer_by_phase_{sample}.csv'))
    tertile_stats.to_csv(os.path.join(fig_dir, f'titer_by_tertile_{sample}.csv'))
    
    print("\n" + "="*60)
    print("TITER-CELL CYCLE ANALYSIS COMPLETE")
    print("="*60)
    
    return results


def create_titer_cellcycle_plots(adata, fig_dir, sample, titer_col, h_stat, p_kw, spearman_r, spearman_p):
    """
    Create comprehensive plots for titer-cell cycle associations
    """
    
    # 1. Violin plots by phase
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    ax = axes[0]
    phases = ['G1', 'S', 'G2M']
    colors = {'G1': '#FF6B6B', 'S': '#4ECDC4', 'G2M': '#45B7D1'}
    
    data_by_phase = [adata.obs[adata.obs['phase'] == p][titer_col].values for p in phases]
    parts = ax.violinplot(data_by_phase, positions=range(len(phases)), 
                          showmeans=True, showmedians=True)
    
    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(colors[phases[i]])
        pc.set_alpha(0.7)
    
    ax.set_xticks(range(len(phases)))
    ax.set_xticklabels(phases)
    ax.set_xlabel('Cell Cycle Phase', fontsize=12)
    ax.set_ylabel(f'Wolbachia Titer ({titer_col})', fontsize=12)
    ax.set_title(f'Titer by Cell Cycle Phase\nKruskal-Wallis H = {h_stat:.2f}, p = {p_kw:.2e}')
    ax.grid(axis='y', alpha=0.3)
    
    # Box plot overlay
    ax = axes[1]
    phase_data = [adata.obs[adata.obs['phase'] == p][titer_col] for p in phases]
    bp = ax.boxplot(phase_data, labels=phases, patch_artist=True, 
                    notch=True, showmeans=True)
    
    for patch, phase in zip(bp['boxes'], phases):
        patch.set_facecolor(colors[phase])
        patch.set_alpha(0.7)
    
    ax.set_xlabel('Cell Cycle Phase', fontsize=12)
    ax.set_ylabel(f'Wolbachia Titer ({titer_col})', fontsize=12)
    ax.set_title(f'Titer Distribution by Phase')
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'{sample}_titer_by_phase.pdf'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Scatter: titer vs pseudotime
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # Linear scatter
    ax = axes[0]
    colors_phase = [{'G1': '#FF6B6B', 'S': '#4ECDC4', 'G2M': '#45B7D1'}[p] 
                    for p in adata.obs['phase']]
    
    ax.scatter(adata.obs['cellcycle_pseudotime'], adata.obs[titer_col],
              c=colors_phase, alpha=0.3, s=10)
    
    # Add regression line
    z = np.polyfit(adata.obs['cellcycle_pseudotime'], adata.obs[titer_col], 1)
    p_fit = np.poly1d(z)
    x_line = np.linspace(0, 1, 100)
    ax.plot(x_line, p_fit(x_line), "r--", alpha=0.8, linewidth=2)
    
    ax.set_xlabel('Cell Cycle Pseudotime (0-1)', fontsize=12)
    ax.set_ylabel(f'Wolbachia Titer ({titer_col})', fontsize=12)
    ax.set_title(f'Titer vs Cell Cycle Pseudotime\nSpearman ρ = {spearman_r:.3f}, p = {spearman_p:.2e}')
    ax.axvline(x=0.33, color='gray', linestyle='--', alpha=0.3, label='G1/S')
    ax.axvline(x=0.66, color='gray', linestyle='--', alpha=0.3, label='S/G2M')
    ax.legend()
    ax.grid(alpha=0.3)
    
    # Circular scatter
    ax = axes[1]
    ax = plt.subplot(1, 2, 2, projection='polar')
    theta = adata.obs['cellcycle_pseudotime'].values * 2 * np.pi
    r = adata.obs[titer_col].values
    
    # Normalize r for better visualization
    r_norm = (r - r.min()) / (r.max() - r.min() + 1e-10)
    
    scatter = ax.scatter(theta, r_norm, c=r, cmap='viridis', alpha=0.5, s=20)
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_ylim(0, 1)
    ax.set_title(f'Circular View: Titer vs Pseudotime\n(color = titer)', pad=20)
    
    plt.colorbar(scatter, ax=ax, label=f'Titer ({titer_col})', pad=0.1)
    
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'{sample}_titer_vs_pseudotime.pdf'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. Hexbin density plot
    fig, ax = plt.subplots(figsize=(10, 8))
    
    hb = ax.hexbin(adata.obs['cellcycle_pseudotime'], adata.obs[titer_col],
                   gridsize=30, cmap='YlOrRd', mincnt=1)
    
    ax.set_xlabel('Cell Cycle Pseudotime (0-1)', fontsize=12)
    ax.set_ylabel(f'Wolbachia Titer ({titer_col})', fontsize=12)
    ax.set_title(f'Density: Titer vs Pseudotime\nSpearman ρ = {spearman_r:.3f}, p = {spearman_p:.2e}')
    ax.axvline(x=0.33, color='blue', linestyle='--', alpha=0.5, label='G1/S')
    ax.axvline(x=0.66, color='green', linestyle='--', alpha=0.5, label='S/G2M')
    ax.legend()
    
    plt.colorbar(hb, ax=ax, label='Cell count')
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'{sample}_titer_pseudotime_density.pdf'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # 4. Binned analysis
    if 'pseudotime_tertile' in adata.obs.columns:
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))
        
        # Violin by tertile
        ax = axes[0]
        tertiles = ['Early', 'Mid', 'Late']
        tertile_data = [adata.obs[adata.obs['pseudotime_tertile'] == t][titer_col].values for t in tertiles]
        
        parts = ax.violinplot(tertile_data, positions=range(len(tertiles)), 
                              showmeans=True, showmedians=True)
        
        ax.set_xticks(range(len(tertiles)))
        ax.set_xticklabels(tertiles)
        ax.set_xlabel('Pseudotime Tertile', fontsize=12)
        ax.set_ylabel(f'Wolbachia Titer ({titer_col})', fontsize=12)
        ax.set_title('Titer by Pseudotime Tertile')
        ax.grid(axis='y', alpha=0.3)
        
        # Bar plot with error bars
        ax = axes[1]
        tertile_stats = adata.obs.groupby('pseudotime_tertile')[titer_col].agg(['mean', 'std'])
        
        x_pos = range(len(tertiles))
        ax.bar(x_pos, tertile_stats['mean'], yerr=tertile_stats['std'],
              alpha=0.7, capsize=5, color=['#FF6B6B', '#4ECDC4', '#45B7D1'])
        
        ax.set_xticks(x_pos)
        ax.set_xticklabels(tertiles)
        ax.set_xlabel('Pseudotime Tertile', fontsize=12)
        ax.set_ylabel(f'Mean Titer ({titer_col})', fontsize=12)
        ax.set_title('Mean Titer by Pseudotime Tertile')
        ax.grid(axis='y', alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f'{sample}_titer_by_tertile.pdf'), dpi=300, bbox_inches='tight')
        plt.close()
    
    # 5. UMAP with titer
    if 'X_umap' in adata.obsm:
        fig, axes = plt.subplots(1, 3, figsize=(21, 6))
        
        # Phase
        sc.pl.umap(adata, color='phase', ax=axes[0], show=False, title='Cell Cycle Phase')
        
        # Pseudotime
        sc.pl.umap(adata, color='cellcycle_pseudotime', ax=axes[1], show=False, 
                  title='Cell Cycle Pseudotime', cmap='twilight')
        
        # Titer
        sc.pl.umap(adata, color=titer_col, ax=axes[2], show=False,
                  title=f'Wolbachia Titer', cmap='viridis')
        
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f'{sample}_umap_phase_pseudotime_titer.pdf'), 
                   dpi=300, bbox_inches='tight')
        plt.close()
    
    print(f"\nTiter-cell cycle plots saved to {fig_dir}")


def analyze_titer_cluster_association(adata, fig_dir, sample, titer_col):
    """
    Analyze association between Wolbachia titer and Leiden clusters
    """
    print("\n" + "="*60)
    print("WOLBACHIA TITER - CLUSTER ASSOCIATION ANALYSIS")
    print("="*60)
    
    if 'leiden' not in adata.obs.columns:
        print("ERROR: No Leiden clustering found")
        return None
    
    if titer_col not in adata.obs.columns:
        print(f"ERROR: Titer column '{titer_col}' not found")
        return None
    
    # Remove cells with missing titer
    valid_mask = ~adata.obs[titer_col].isna()
    adata_valid = adata[valid_mask, :].copy()
    
    clusters = sorted(adata_valid.obs['leiden'].unique())
    
    # Kruskal-Wallis test
    groups = [adata_valid.obs[adata_valid.obs['leiden'] == c][titer_col].values for c in clusters]
    h_stat, p_value = kruskal(*groups)
    
    print(f"Kruskal-Wallis test: H = {h_stat:.2f}, p = {p_value:.2e}")
    print(f"Conclusion: Titer {'SIGNIFICANTLY differs' if p_value < 0.05 else 'does NOT significantly differ'} across clusters")
    
    # Summary statistics
    titer_by_cluster = adata_valid.obs.groupby('leiden')[titer_col].agg(['mean', 'median', 'std', 'count'])
    print("\nTiter by cluster:")
    print(titer_by_cluster)
    
    # Create plots
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # Violin plot
    ax = axes[0]
    cluster_data = [adata_valid.obs[adata_valid.obs['leiden'] == c][titer_col].values for c in clusters]
    parts = ax.violinplot(cluster_data, positions=range(len(clusters)), 
                          showmeans=True, showmedians=True)
    
    ax.set_xticks(range(len(clusters)))
    ax.set_xticklabels(clusters)
    ax.set_xlabel('Leiden Cluster', fontsize=12)
    ax.set_ylabel(f'Wolbachia Titer ({titer_col})', fontsize=12)
    ax.set_title(f'Titer by Cluster\nKruskal-Wallis H = {h_stat:.2f}, p = {p_value:.2e}')
    ax.grid(axis='y', alpha=0.3)
    
    # Bar plot with error bars
    ax = axes[1]
    x_pos = range(len(clusters))
    ax.bar(x_pos, titer_by_cluster['mean'], yerr=titer_by_cluster['std'],
          alpha=0.7, capsize=5)
    
    ax.set_xticks(x_pos)
    ax.set_xticklabels(clusters)
    ax.set_xlabel('Leiden Cluster', fontsize=12)
    ax.set_ylabel(f'Mean Titer ({titer_col})', fontsize=12)
    ax.set_title('Mean Titer by Cluster')
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'{sample}_titer_by_cluster.pdf'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # Save results
    titer_by_cluster.to_csv(os.path.join(fig_dir, f'titer_by_cluster_{sample}.csv'))
    
    stats_results = {
        'Test': ['Kruskal-Wallis (Clusters)'],
        'Statistic': [f'{h_stat:.6e}'],
        'P-value': [f'{p_value:.6e}'],
        'Interpretation': ['Titer differs across clusters' if p_value < 0.05 else 'No difference across clusters']
    }
    
    stats_df = pd.DataFrame(stats_results)
    stats_df.to_csv(os.path.join(fig_dir, f'titer_cluster_stats_{sample}.csv'), index=False)
    
    print("\n" + "="*60)
    print("TITER-CLUSTER ANALYSIS COMPLETE")
    print("="*60)
    
    return {'kw_H': h_stat, 'kw_p': p_value, 'titer_by_cluster': titer_by_cluster}


def analyze_cluster_cellcycle_association(adata, fig_dir, sample):
    """
    Test and visualize association between Leiden clusters and cell cycle
    (Keeping your original function with minor formatting updates)
    """
    print("\n" + "="*60)
    print("CLUSTER - CELL CYCLE ASSOCIATION ANALYSIS")
    print("="*60)
    
    if 'leiden' not in adata.obs.columns:
        print("ERROR: No 'leiden' clustering found in adata.obs")
        return None
    
    if 'phase' not in adata.obs.columns:
        print("ERROR: No 'phase' found in adata.obs")
        return None
    
    sc.settings.figdir = fig_dir
    
    leiden_colors = []
    clusters = sorted(adata.obs['leiden'].unique())
    cmap = plt.cm.get_cmap('tab20')
    for i, cluster in enumerate(clusters):
        leiden_colors.append(cmap(i % 20))
    
    # Chi-square test
    contingency = pd.crosstab(adata.obs['leiden'], adata.obs['phase'])
    chi2, p_value, dof, expected = chi2_contingency(contingency)
    
    print(f"\nχ² = {chi2:.2f}")
    print(f"p-value = {p_value:.2e}")
    
    n = contingency.sum().sum()
    cramers_v = np.sqrt(chi2 / (n * (min(contingency.shape) - 1)))
    print(f"Cramér's V = {cramers_v:.3f}")
    
    # Heatmap
    contingency_norm = contingency.div(contingency.sum(axis=1), axis=0) * 100
    
    fig, ax = plt.subplots(figsize=(12, 8))
    sns.heatmap(contingency_norm, annot=True, fmt='.1f', cmap='YlOrRd', ax=ax)
    ax.set_xlabel('Cell Cycle Phase')
    ax.set_ylabel('Leiden Cluster')
    ax.set_title(f'Cell cycle phase by cluster\nχ² = {chi2:.2f}, p = {p_value:.2e}')
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'heatmap_cluster_cellcycle_{sample}.pdf'))
    plt.close()
    
    # Save results
    contingency.to_csv(os.path.join(fig_dir, f'contingency_counts_{sample}.csv'))
    contingency_norm.to_csv(os.path.join(fig_dir, f'contingency_percentages_{sample}.csv'))
    
    stats_results = {
        'Test': ['Chi-square', 'Cramers V'],
        'Statistic': [f'{chi2:.6e}', f'{cramers_v:.6e}'],
        'P-value': [f'{p_value:.6e}', 'NA'],
        'Interpretation': [
            'Significant association' if p_value < 0.05 else 'No significant association',
            f'Effect size'
        ]
    }
    
    stats_df = pd.DataFrame(stats_results)
    stats_df.to_csv(os.path.join(fig_dir, f'cluster_cellcycle_stats_{sample}.csv'), index=False)
    
    return {'chi2': chi2, 'chi2_pvalue': p_value, 'cramers_v': cramers_v}


def main():
    parser = argparse.ArgumentParser(
        description='Cell cycle annotation with Wolbachia titer association analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Full analysis with auto-detected titer column
  python scanpy_titer_cellcycle.py \\
      --input integrated.h5ad \\
      --output titer_cellcycle_analysis \\
      --sample all_conditions \\
      --run-scoring
  
  # Specify titer column
  python scanpy_titer_cellcycle.py \\
      --input integrated.h5ad \\
      --output titer_cellcycle_analysis \\
      --sample all_conditions \\
      --titer-column wolbachia_titer
  
  # Just titer analysis (scoring already done)
  python scanpy_titer_cellcycle.py \\
      --input integrated.h5ad \\
      --output titer_cellcycle_analysis \\
      --sample all_conditions \\
      --titer-column wolbachia_titer
        '''
    )
    
    parser.add_argument('--input', '-i', required=True, help='Path to h5ad file')
    parser.add_argument('--output', '-o', default='titer_cellcycle_analysis', help='Output directory')
    parser.add_argument('--sample', '-s', default='sample', help='Sample name')
    parser.add_argument('--run-scoring', action='store_true', help='Run cell cycle scoring')
    parser.add_argument('--titer-column', type=str, default=None,
                        help='Name of Wolbachia titer column (auto-detect if not specified)')
    parser.add_argument('--save-output', action='store_true', help='Save updated h5ad')
    
    args = parser.parse_args()
    
    os.makedirs(args.output, exist_ok=True)
    
    print(f"Loading data from {args.input}...")
    adata = sc.read_h5ad(args.input)
    
    print(f"\nLoaded: {adata.n_obs} cells × {adata.n_vars} genes")
    
    # Run cell cycle scoring
    has_scoring = all(col in adata.obs.columns for col in ['phase', 'S_score', 'G2M_score'])
    
    if args.run_scoring or not has_scoring:
        result = score_cell_cycle_scanpy(adata, args.output, args.sample)
        if result is None:
            print("\nERROR: Cell cycle scoring failed")
            return
        adata = result
        
        if args.save_output:
            output_path = args.input.replace('.h5ad', '_with_cellcycle.h5ad')
            print(f"\nSaving to {output_path}")
            adata.write(output_path)
    else:
        # Still need pseudotime if not present
        if 'cellcycle_pseudotime' not in adata.obs.columns:
            print("\nCalculating cell cycle pseudotime...")
            adata = calculate_cellcycle_pseudotime(adata)
    
    # Detect or use specified titer column
    if args.titer_column:
        titer_col = args.titer_column
        if titer_col not in adata.obs.columns:
            print(f"\nERROR: Specified titer column '{titer_col}' not found!")
            return
    else:
        titer_col = detect_titer_column(adata)
        if titer_col is None:
            print("\nERROR: Could not detect titer column")
            return
    
    # Run analyses
    print("\n" + "="*60)
    print("RUNNING ASSOCIATION ANALYSES")
    print("="*60)
    
    # 1. Titer vs cell cycle
    titer_cc_results = analyze_titer_cellcycle_association(adata, args.output, args.sample, titer_col)
    
    # 2. Titer vs clusters
    titer_cluster_results = analyze_titer_cluster_association(adata, args.output, args.sample, titer_col)
    
    # 3. Clusters vs cell cycle (original analysis)
    cluster_cc_results = analyze_cluster_cellcycle_association(adata, args.output, args.sample)
    
    # Summary
    print("\n" + "="*60)
    print("SUMMARY OF KEY FINDINGS")
    print("="*60)
    
    if titer_cc_results:
        print("\n1. TITER vs CELL CYCLE:")
        print(f"   Kruskal-Wallis (phase): p = {titer_cc_results['kw_phase_p']:.2e}")
        print(f"   Spearman correlation: ρ = {titer_cc_results['spearman_r']:.3f}, p = {titer_cc_results['spearman_p']:.2e}")
    
    if titer_cluster_results:
        print("\n2. TITER vs CLUSTERS:")
        print(f"   Kruskal-Wallis: p = {titer_cluster_results['kw_p']:.2e}")
    
    if cluster_cc_results:
        print("\n3. CLUSTERS vs CELL CYCLE:")
        print(f"   Chi-square: p = {cluster_cc_results['chi2_pvalue']:.2e}")
        print(f"   Cramér's V = {cluster_cc_results['cramers_v']:.3f}")
    
    print("\n" + "="*60)
    print(f"All results saved to: {args.output}")
    print("="*60)


if __name__ == "__main__":
    main()