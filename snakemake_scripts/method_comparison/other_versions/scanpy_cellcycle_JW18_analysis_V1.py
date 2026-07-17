'''
Cell cycle annotation using Scanpy with Drosophila-specific markers
Filtered for JW18wMel-Ctrl samples with comparison to JW18DOX-Ctrl
Analysis of Wolbachia titer associations and infection effects on cell cycle
'''
import scanpy as sc 
import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chi2_contingency, kruskal, spearmanr, pearsonr, mannwhitneyu
import scikit_posthocs as sp

# Drosophila cell cycle genes (FlyBase IDs)
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


def filter_samples(adata, wmel_filter='JW18wMel-Ctrl', dox_filter='JW18DOX-Ctrl'):
    """
    Filter AnnData for specific samples
    Returns two filtered datasets: wMel-infected and Uninfected
    """
    print("\n" + "="*60)
    print("FILTERING DATA FOR ANALYSIS")
    print("="*60)
    
    print(f"Original data: {adata.n_obs} cells × {adata.n_vars} genes")
    
    # Use bio_condition column
    sample_col = 'bio_condition'
    
    if sample_col not in adata.obs.columns:
        print(f"\nERROR: Column '{sample_col}' not found!")
        print("Available columns:")
        print(adata.obs.columns.tolist())
        return None, None
    
    print(f"\nUsing column: '{sample_col}'")
    print(f"Unique values: {adata.obs[sample_col].unique()}")
    
    # Filter for wMel-infected samples (JW18wMel-Ctrl)
    print(f"\n{'='*60}")
    print(f"FILTERING FOR wMel-INFECTED: {wmel_filter}")
    print(f"{'='*60}")
    
    wmel_mask = adata.obs[sample_col] == wmel_filter
    
    if wmel_mask.sum() == 0:
        print(f"ERROR: No samples found matching '{wmel_filter}'")
        print(f"Available values in '{sample_col}': {adata.obs[sample_col].unique()}")
        return None, None
    
    adata_wmel = adata[wmel_mask, :].copy()
    print(f"wMel-infected samples: {adata_wmel.n_obs} cells")
    print(f"bio_condition values: {adata_wmel.obs[sample_col].unique()}")
    
    # Filter for uninfected samples (JW18DOX-Ctrl)
    print(f"\n{'='*60}")
    print(f"FILTERING FOR UNINFECTED: {dox_filter}")
    print(f"{'='*60}")
    
    dox_mask = adata.obs[sample_col] == dox_filter
    
    if dox_mask.sum() == 0:
        print(f"WARNING: No samples found matching '{dox_filter}'")
        print(f"Available values in '{sample_col}': {adata.obs[sample_col].unique()}")
        print("Uninfected comparison will be skipped.")
        adata_uninfected = None
    else:
        adata_uninfected = adata[dox_mask, :].copy()
        print(f"Uninfected samples: {adata_uninfected.n_obs} cells")
        print(f"bio_condition values: {adata_uninfected.obs[sample_col].unique()}")
    
    return adata_wmel, adata_uninfected

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
            if pd.api.types.is_numeric_dtype(adata.obs[col]):
                print(f"Potential titer column: '{col}'")
                print(f"  Range: {adata.obs[col].min():.3f} to {adata.obs[col].max():.3f}")
                print(f"  Mean: {adata.obs[col].mean():.3f}")
                return col
    
    print("\nWARNING: No titer column found!")
    return None


def score_cell_cycle_scanpy(adata, output_dir, sample_name, s_genes=None, g2m_genes=None):
    """
    Score cell cycle using Scanpy with Drosophila genes
    """
    print("\n" + "="*60)
    print(f"SCANPY CELL CYCLE SCORING: {sample_name}")
    print("="*60)
    
    if s_genes is None:
        s_genes = S_GENES_FBGN
    if g2m_genes is None:
        g2m_genes = G2M_GENES_FBGN
    
    s_genes_present = [g for g in s_genes if g in adata.var_names]
    g2m_genes_present = [g for g in g2m_genes if g in adata.var_names]
    
    print(f"\nS phase genes: {len(s_genes_present)}/{len(s_genes)} found")
    print(f"G2/M phase genes: {len(g2m_genes_present)}/{len(g2m_genes)} found")
    
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
    
    # Add continuous pseudotime
    print("\nCalculating continuous cell cycle pseudotime...")
    adata = calculate_cellcycle_pseudotime(adata)
    
    return adata


def calculate_cellcycle_pseudotime(adata):
    """
    Calculate a continuous pseudotime measure (0-1) through the cell cycle
    """
    s_score = adata.obs['S_score'].values
    g2m_score = adata.obs['G2M_score'].values
    phase = adata.obs['phase'].values
    
    pseudotime = np.zeros(len(phase))
    
    for i, (p, s, g) in enumerate(zip(phase, s_score, g2m_score)):
        if p == 'G1':
            base = 0.165
            jitter = (s + g) * 0.05
            pseudotime[i] = np.clip(base + jitter, 0.0, 0.33)
        elif p == 'S':
            base = 0.495
            s_in_s = s_score[phase == 'S']
            if len(s_in_s) > 0:
                s_norm = (s - s_in_s.min()) / (s_in_s.max() - s_in_s.min() + 1e-10)
                jitter = (s_norm - 0.5) * 0.25
            else:
                jitter = 0
            pseudotime[i] = np.clip(base + jitter, 0.33, 0.66)
        else:  # G2M
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
    
    return adata


def compare_cellcycle_distributions(adata_wmel, adata_uninfected, fig_dir, sample_name):
    """
    Compare cell cycle phase distributions between wMel-infected and Uninfected
    """
    print("\n" + "="*60)
    print("COMPARING CELL CYCLE DISTRIBUTIONS: wMel-infected vs Uninfected")
    print("="*60)
    
    # Create combined dataframe
    wmel_df = pd.DataFrame({
        'phase': adata_wmel.obs['phase'],
        'condition': 'wMel-infected'
    })
    
    uninfected_df = pd.DataFrame({
        'phase': adata_uninfected.obs['phase'],
        'condition': 'Uninfected'
    })
    
    combined_df = pd.concat([wmel_df, uninfected_df])
    
    # 1. Chi-square test
    print("\n1. CHI-SQUARE TEST")
    contingency = pd.crosstab(combined_df['condition'], combined_df['phase'])
    chi2, p_value, dof, expected = chi2_contingency(contingency)
    
    print(f"χ² = {chi2:.2f}, p = {p_value:.2e}")
    print(f"Conclusion: Cell cycle distributions {'SIGNIFICANTLY differ' if p_value < 0.05 else 'do NOT differ'} between wMel-infected and Uninfected")
    
    # 2. Phase-specific comparisons
    print("\n2. PHASE-SPECIFIC COMPARISONS")
    phase_tests = {}
    for phase in ['G1', 'S', 'G2M']:
        wmel_count = (adata_wmel.obs['phase'] == phase).sum()
        uninfected_count = (adata_uninfected.obs['phase'] == phase).sum()
        
        wmel_pct = (wmel_count / len(adata_wmel)) * 100
        uninfected_pct = (uninfected_count / len(adata_uninfected)) * 100
        
        print(f"\n{phase} phase:")
        print(f"  wMel-infected: {wmel_count} cells ({wmel_pct:.1f}%)")
        print(f"  Uninfected: {uninfected_count} cells ({uninfected_pct:.1f}%)")
        print(f"  Difference: {wmel_pct - uninfected_pct:+.1f}%")
        
        phase_tests[phase] = {
            'wmel_pct': wmel_pct,
            'uninfected_pct': uninfected_pct,
            'diff': wmel_pct - uninfected_pct
        }
    
    # 3. Create comparison plots
    create_infection_comparison_plots(adata_wmel, adata_uninfected, fig_dir, sample_name,
                                      chi2, p_value, phase_tests)
    
    # Save results
    contingency_norm = contingency.div(contingency.sum(axis=1), axis=0) * 100
    contingency.to_csv(os.path.join(fig_dir, f'{sample_name}_infection_comparison_counts.csv'))
    contingency_norm.to_csv(os.path.join(fig_dir, f'{sample_name}_infection_comparison_percentages.csv'))
    
    stats_results = {
        'Test': ['Chi-square'],
        'Statistic': [f'{chi2:.6e}'],
        'P-value': [f'{p_value:.6e}'],
        'Interpretation': ['Distributions differ' if p_value < 0.05 else 'No difference']
    }
    
    stats_df = pd.DataFrame(stats_results)
    stats_df.to_csv(os.path.join(fig_dir, f'{sample_name}_infection_comparison_stats.csv'), index=False)
    
    return {'chi2': chi2, 'p_value': p_value, 'phase_tests': phase_tests}


def create_infection_comparison_plots(adata_wmel, adata_uninfected, fig_dir, sample, chi2, p_val, phase_tests):
    """
    Create plots comparing cell cycle between wMel-infected and Uninfected
    """
    colors = {'G1': '#FF6B6B', 'S': '#4ECDC4', 'G2M': '#45B7D1'}
    
    # 1. Side-by-side bar plots
    fig, axes = plt.subplots(1, 3, figsize=(21, 6))
    
    # Bar plot comparison
    ax = axes[0]
    phases = ['G1', 'S', 'G2M']
    wmel_pcts = [phase_tests[p]['wmel_pct'] for p in phases]
    uninfected_pcts = [phase_tests[p]['uninfected_pct'] for p in phases]
    
    x = np.arange(len(phases))
    width = 0.35
    
    ax.bar(x - width/2, wmel_pcts, width, label='wMel-infected', color='#E74C3C', alpha=0.8)
    ax.bar(x + width/2, uninfected_pcts, width, label='Uninfected', color='#3498DB', alpha=0.8)
    
    ax.set_xlabel('Cell Cycle Phase', fontsize=12)
    ax.set_ylabel('Percentage of Cells', fontsize=12)
    ax.set_title(f'Cell Cycle Distribution\nχ² = {chi2:.2f}, p = {p_val:.2e}')
    ax.set_xticks(x)
    ax.set_xticklabels(phases)
    ax.legend()
    ax.grid(axis='y', alpha=0.3)
    
    # Difference plot
    ax = axes[1]
    diffs = [phase_tests[p]['diff'] for p in phases]
    colors_diff = ['red' if d < 0 else 'green' for d in diffs]
    
    ax.bar(phases, diffs, color=colors_diff, alpha=0.7)
    ax.axhline(y=0, color='k', linestyle='-', linewidth=1)
    ax.set_xlabel('Cell Cycle Phase', fontsize=12)
    ax.set_ylabel('Difference (wMel-infected - Uninfected) %', fontsize=12)
    ax.set_title('Cell Cycle Phase Enrichment\n(Positive = more in wMel-infected)')
    ax.grid(axis='y', alpha=0.3)
    
    # Stacked bar plot
    ax = axes[2]
    
    wmel_counts = [adata_wmel.obs['phase'].value_counts().get(p, 0) for p in phases]
    uninfected_counts = [adata_uninfected.obs['phase'].value_counts().get(p, 0) for p in phases]
    
    wmel_norm = np.array(wmel_counts) / sum(wmel_counts) * 100
    uninfected_norm = np.array(uninfected_counts) / sum(uninfected_counts) * 100
    
    ax.bar(['wMel-infected'], [wmel_norm[0]], label='G1', color=colors['G1'])
    ax.bar(['wMel-infected'], [wmel_norm[1]], bottom=wmel_norm[0], label='S', color=colors['S'])
    ax.bar(['wMel-infected'], [wmel_norm[2]], bottom=wmel_norm[0]+wmel_norm[1], label='G2M', color=colors['G2M'])
    
    ax.bar(['Uninfected'], [uninfected_norm[0]], color=colors['G1'])
    ax.bar(['Uninfected'], [uninfected_norm[1]], bottom=uninfected_norm[0], color=colors['S'])
    ax.bar(['Uninfected'], [uninfected_norm[2]], bottom=uninfected_norm[0]+uninfected_norm[1], color=colors['G2M'])
    
    ax.set_ylabel('Percentage', fontsize=12)
    ax.set_title('Cell Cycle Composition')
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'{sample}_infection_comparison.pdf'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Circular plots
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), subplot_kw=dict(projection='polar'))
    
    # wMel-infected
    ax = axes[0]
    theta_wmel = adata_wmel.obs['cellcycle_pseudotime'].values * 2 * np.pi
    phase_colors_wmel = [colors[p] for p in adata_wmel.obs['phase']]
    ax.scatter(theta_wmel, np.ones(len(theta_wmel)), c=phase_colors_wmel, alpha=0.3, s=5)
    ax.set_ylim(0, 1.5)
    ax.set_yticks([])
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_title('wMel-infected (JW18wMel-Ctrl)', pad=20)
    
    # Uninfected
    ax = axes[1]
    theta_uninfected = adata_uninfected.obs['cellcycle_pseudotime'].values * 2 * np.pi
    phase_colors_uninfected = [colors[p] for p in adata_uninfected.obs['phase']]
    ax.scatter(theta_uninfected, np.ones(len(theta_uninfected)), c=phase_colors_uninfected, alpha=0.3, s=5)
    ax.set_ylim(0, 1.5)
    ax.set_yticks([])
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_title('Uninfected (JW18DOX-Ctrl)', pad=20)
    
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'{sample}_circular_comparison.pdf'), dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"\nInfection comparison plots saved to {fig_dir}")


def analyze_titer_cellcycle_association(adata, fig_dir, sample, titer_col):
    """
    Analyze association between Wolbachia titer and cell cycle
    """
    print("\n" + "="*60)
    print("WOLBACHIA TITER - CELL CYCLE ASSOCIATION ANALYSIS")
    print("="*60)
    
    if titer_col is None or titer_col not in adata.obs.columns:
        print(f"WARNING: Titer column not found, skipping titer analysis")
        return None
    
    # Remove cells with missing titer
    valid_mask = ~adata.obs[titer_col].isna()
    n_removed = (~valid_mask).sum()
    if n_removed > 0:
        print(f"\nRemoving {n_removed} cells with missing titer values")
        adata_valid = adata[valid_mask, :].copy()
    else:
        adata_valid = adata
    
    print(f"Analyzing {adata_valid.n_obs} cells with valid titer data")
    
    # Test correlation with pseudotime
    print("\nTITER vs CELL CYCLE PSEUDOTIME:")
    pseudotime = adata_valid.obs['cellcycle_pseudotime'].values
    titer = adata_valid.obs[titer_col].values
    
    pearson_r, pearson_p = pearsonr(pseudotime, titer)
    spearman_r, spearman_p = spearmanr(pseudotime, titer)
    
    print(f"Pearson correlation: r = {pearson_r:.4f}, p = {pearson_p:.2e}")
    print(f"Spearman correlation: ρ = {spearman_r:.4f}, p = {spearman_p:.2e}")
    
    # Test by phase
    print("\nTITER vs CELL CYCLE PHASE:")
    phases = ['G1', 'S', 'G2M']
    groups = [adata_valid.obs[adata_valid.obs['phase'] == p][titer_col].values for p in phases]
    h_stat, p_kw = kruskal(*groups)
    
    print(f"Kruskal-Wallis test: H = {h_stat:.2f}, p = {p_kw:.2e}")
    
    titer_stats = adata_valid.obs.groupby('phase')[titer_col].agg(['mean', 'median', 'std'])
    print("\nTiter by phase:")
    print(titer_stats)
    
    # Create plots
    create_titer_plots(adata_valid, fig_dir, sample, titer_col, spearman_r, spearman_p, h_stat, p_kw)
    
    # Save results
    stats_results = {
        'Test': ['Kruskal-Wallis (Phase)', 'Pearson correlation', 'Spearman correlation'],
        'Statistic': [f'{h_stat:.6e}', f'{pearson_r:.6e}', f'{spearman_r:.6e}'],
        'P-value': [f'{p_kw:.6e}', f'{pearson_p:.6e}', f'{spearman_p:.6e}'],
        'Interpretation': [
            'Titer differs across phases' if p_kw < 0.05 else 'No difference',
            f'Pearson r = {pearson_r:.3f}',
            f'Spearman ρ = {spearman_r:.3f}'
        ]
    }
    
    stats_df = pd.DataFrame(stats_results)
    stats_df.to_csv(os.path.join(fig_dir, f'{sample}_titer_cellcycle_stats.csv'), index=False)
    titer_stats.to_csv(os.path.join(fig_dir, f'{sample}_titer_by_phase.csv'))
    
    return {
        'spearman_r': spearman_r,
        'spearman_p': spearman_p,
        'kw_p': p_kw
    }


def create_titer_plots(adata, fig_dir, sample, titer_col, spearman_r, spearman_p, h_stat, p_kw):
    """
    Create titer-cell cycle plots
    """
    colors = {'G1': '#FF6B6B', 'S': '#4ECDC4', 'G2M': '#45B7D1'}
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # Violin plot by phase
    ax = axes[0, 0]
    phases = ['G1', 'S', 'G2M']
    data_by_phase = [adata.obs[adata.obs['phase'] == p][titer_col].values for p in phases]
    parts = ax.violinplot(data_by_phase, positions=range(len(phases)), showmeans=True, showmedians=True)
    
    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(colors[phases[i]])
        pc.set_alpha(0.7)
    
    ax.set_xticks(range(len(phases)))
    ax.set_xticklabels(phases)
    ax.set_xlabel('Cell Cycle Phase')
    ax.set_ylabel(f'Wolbachia Titer')
    ax.set_title(f'Titer by Phase\nKruskal-Wallis H = {h_stat:.2f}, p = {p_kw:.2e}')
    ax.grid(axis='y', alpha=0.3)
    
    # Scatter: titer vs pseudotime
    ax = axes[0, 1]
    phase_colors_list = [colors[p] for p in adata.obs['phase']]
    ax.scatter(adata.obs['cellcycle_pseudotime'], adata.obs[titer_col],
              c=phase_colors_list, alpha=0.3, s=10)
    
    # Regression line
    z = np.polyfit(adata.obs['cellcycle_pseudotime'], adata.obs[titer_col], 1)
    p_fit = np.poly1d(z)
    x_line = np.linspace(0, 1, 100)
    ax.plot(x_line, p_fit(x_line), "r--", alpha=0.8, linewidth=2)
    
    ax.set_xlabel('Cell Cycle Pseudotime (0-1)')
    ax.set_ylabel(f'Wolbachia Titer')
    ax.set_title(f'Titer vs Pseudotime\nSpearman ρ = {spearman_r:.3f}, p = {spearman_p:.2e}')
    ax.axvline(x=0.33, color='gray', linestyle='--', alpha=0.3, label='G1/S')
    ax.axvline(x=0.66, color='gray', linestyle='--', alpha=0.3, label='S/G2M')
    ax.legend()
    ax.grid(alpha=0.3)
    
    # Hexbin density
    ax = axes[1, 0]
    hb = ax.hexbin(adata.obs['cellcycle_pseudotime'], adata.obs[titer_col],
                   gridsize=30, cmap='YlOrRd', mincnt=1)
    ax.set_xlabel('Cell Cycle Pseudotime (0-1)')
    ax.set_ylabel(f'Wolbachia Titer')
    ax.set_title('Density: Titer vs Pseudotime')
    ax.axvline(x=0.33, color='blue', linestyle='--', alpha=0.5)
    ax.axvline(x=0.66, color='green', linestyle='--', alpha=0.5)
    plt.colorbar(hb, ax=ax, label='Cell count')
    
    # Circular plot
    ax = axes[1, 1]
    ax = plt.subplot(2, 2, 4, projection='polar')
    theta = adata.obs['cellcycle_pseudotime'].values * 2 * np.pi
    r = adata.obs[titer_col].values
    r_norm = (r - r.min()) / (r.max() - r.min() + 1e-10)
    
    scatter = ax.scatter(theta, r_norm, c=r, cmap='viridis', alpha=0.5, s=20)
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_ylim(0, 1)
    ax.set_title('Circular: Titer vs Pseudotime', pad=20)
    plt.colorbar(scatter, ax=ax, label='Titer', pad=0.1)
    
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'{sample}_titer_cellcycle.pdf'), dpi=300, bbox_inches='tight')
    plt.close()

def analyze_titer_by_cluster(adata, fig_dir, sample, titer_col, cluster_col='leiden'):
    """
    Test if Wolbachia titer is associated with clusters and perform 
    Dunn's post-hoc test for pairwise differences.
    """
    print("\n" + "="*60)
    print(f"WOLBACHIA TITER - CLUSTER ASSOCIATION ({cluster_col})")
    print("="*60)

    if cluster_col not in adata.obs.columns:
        print(f"WARNING: Cluster column '{cluster_col}' not found.")
        return None

    # Clean data
    df = adata.obs[[cluster_col, titer_col]].dropna()
    df[cluster_col] = df[cluster_col].astype(str) # Ensure categories are strings
    
    # 1. Kruskal-Wallis Global Test
    cluster_groups = [group[titer_col].values for name, group in df.groupby(cluster_col)]
    h_stat, p_kw = kruskal(*cluster_groups)
    print(f"Kruskal-Wallis p-value: {p_kw:.2e}")

    # 2. Dunn's Post-hoc Test (Pairwise)
    # This identifies exactly which clusters are different from each other
    p_values_dunn = sp.posthoc_dunn(df, val_col=titer_col, group_col=cluster_col, p_adjust='bonferroni')
    
    # 3. Visualization: Heatmap of Pairwise Significance
    plt.figure(figsize=(10, 8))
    # We use a log scale for the heatmap to highlight small p-values (< 0.05)
    sns.heatmap(p_values_dunn, annot=False, cmap='rocket_r', vmin=0, vmax=0.05)
    plt.title(f"Dunn's Test: Pairwise Titer Significance (p < 0.05)\n{sample}")
    plt.savefig(os.path.join(fig_dir, f'{sample}_titer_cluster_dunn_heatmap.pdf'))
    plt.close()

    # 4. Save Significance Matrix
    p_values_dunn.to_csv(os.path.join(fig_dir, f'{sample}_titer_cluster_dunn_pvals.csv'))

    # Summary Stats
    cluster_stats = df.groupby(cluster_col)[titer_col].agg(['count', 'mean', 'median'])
    cluster_stats = cluster_stats.sort_values(by='median', ascending=False)
    
    print("\nTop 5 Clusters by Titer:")
    print(cluster_stats.head())

    return {'p_kw': p_kw, 'dunn_p': p_values_dunn, 'stats': cluster_stats}

def find_titer_associated_genes(adata, cluster_results, cluster_col='leiden'):
    """
    Identifies genes differentially expressed between clusters with 
    the highest and lowest Wolbachia titers.
    """
    print("\n" + "="*60)
    print("DIFFERENTIAL EXPRESSION: HIGH VS LOW TITER CLUSTERS")
    print("="*60)
    
    # 1. Identify High and Low groups (Top 2 and Bottom 2 clusters)
    stats = cluster_results['stats']
    high_clusters = stats.index[:2].tolist()
    low_clusters = stats.index[-2:].tolist()
    
    print(f"High-titer clusters: {high_clusters}")
    print(f"Low-titer clusters: {low_clusters}")

    # Create a temporary group for DE
    adata.obs['titer_group'] = 'Medium'
    adata.obs.loc[adata.obs[cluster_col].isin(high_clusters), 'titer_group'] = 'High'
    adata.obs.loc[adata.obs[cluster_col].isin(low_clusters), 'titer_group'] = 'Low'

    # 2. Run Rank Genes Groups (Wilcoxon)
    # Comparing High vs Low specifically
    sc.tl.rank_genes_groups(adata, 'titer_group', groups=['High'], reference='Low', method='wilcoxon')

    # 3. Extract results
    result = sc.get.rank_genes_groups_df(adata, group="High")
    
    # Filter for significance and logfoldchange
    sig_genes = result[(result['pvals_adj'] < 0.05) & (abs(result['logfoldchanges']) > 0.5)]
    
    print(f"\nFound {len(sig_genes)} genes differentially expressed between High and Low titer clusters.")
    print("\nTop genes enriched in High-titer clusters:")
    print(sig_genes.head(10)[['names', 'logfoldchanges', 'pvals_adj']])
    
    return sig_genes

def analyze_titer_phase_intersection(adata, fig_dir, sample, cluster_col='leiden'):
    """
    Visualizes the relationship between Leiden clusters, Wolbachia titer, 
    and Cell Cycle phase.
    """
    print("\n" + "="*60)
    print("ANALYZING TITER-PHASE-CLUSTER INTERSECTION")
    print("="*60)

    # 1. Calculate phase proportions per cluster
    # This shows if certain clusters are 'enriched' for specific phases
    phase_proportions = pd.crosstab(adata.obs[cluster_col], adata.obs['phase'], normalize='index') * 100

    # 2. Get median titer per cluster
    titer_col = detect_titer_column(adata)
    cluster_titers = adata.obs.groupby(cluster_col)[titer_col].median()

    # 3. Create a composite visualization
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6), gridspec_kw={'width_ratios': [2, 1]})

    # Stacked Bar: Phase composition per cluster
    phase_proportions.plot(kind='bar', stacked=True, 
                           color={'G1': '#FF6B6B', 'S': '#4ECDC4', 'G2M': '#45B7D1'}, 
                           ax=ax1)
    ax1.set_title("Cell Cycle Phase Distribution per Cluster")
    ax1.set_ylabel("Percentage of Cells")
    ax1.legend(title="Phase", bbox_to_anchor=(1.05, 1), loc='upper left')

    # Scatter: Median Titer vs G2M Percentage
    # High-titer clusters that also have high G2M suggest a cycle delay
    ax2.scatter(phase_proportions['G2M'], cluster_titers, s=100, color='purple', alpha=0.6)
    
    # Label the points with cluster IDs
    for i, txt in enumerate(phase_proportions.index):
        ax2.annotate(txt, (phase_proportions['G2M'].iloc[i], cluster_titers.iloc[i]), fontsize=9)

    ax2.set_xlabel("Percentage of cells in G2/M")
    ax2.set_ylabel("Median Wolbachia Titer")
    ax2.set_title("Titer vs. Mitotic Enrichment")
    ax2.grid(True, linestyle='--', alpha=0.5)

    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'{sample}_titer_phase_intersection.pdf'))
    plt.close()

    print(f"Intersection analysis saved to {fig_dir}")

def main():
    parser = argparse.ArgumentParser(
        description='Cell cycle analysis: JW18wMel-Ctrl vs JW18DOX-Ctrl comparison',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Example:
  python scanpy_cellcycle_JW18_analysis.py \\
      --input integrated.h5ad \\
      --output JW18_cellcycle_analysis \\
      --wmel-filter JW18wMel-Ctrl \\
      --dox-filter JW18DOX-Ctrl \\
      --run-scoring \\
      --save-output
        '''
    )
    
    parser.add_argument('--input', '-i', required=True, help='Path to h5ad file')
    parser.add_argument('--output', '-o', default='JW18_cellcycle_analysis', help='Output directory')
    parser.add_argument('--wmel-filter', type=str, default='JW18wMel-Ctrl',
                        help='Filter string for wMel-infected samples (default: JW18wMel-Ctrl)')
    parser.add_argument('--dox-filter', type=str, default='JW18DOX-Ctrl',
                        help='Filter string for uninfected samples (default: JW18DOX-Ctrl)')
    parser.add_argument('--run-scoring', action='store_true', help='Run cell cycle scoring')
    parser.add_argument('--titer-column', type=str, default=None,
                        help='Name of Wolbachia titer column (auto-detect if not specified)')
    parser.add_argument('--save-output', action='store_true', help='Save filtered h5ad files')
    
    args = parser.parse_args()
    
    os.makedirs(args.output, exist_ok=True)
    
    print(f"Loading data from {args.input}...")
    adata = sc.read_h5ad(args.input)
    print(f"Loaded: {adata.n_obs} cells × {adata.n_vars} genes")
    
    # Filter samples
    adata_wmel, adata_uninfected = filter_samples(adata, args.wmel_filter, args.dox_filter)
    
    if adata_wmel is None:
        print("\nERROR: Could not filter wMel-infected samples. Exiting.")
        return
    
    # Detect titer column (for wMel-infected only)
    if args.titer_column:
        titer_col = args.titer_column
    else:
        titer_col = detect_titer_column(adata_wmel)
    
    # Score cell cycle for wMel-infected samples
    print("\n" + "="*60)
    print("SCORING wMel-INFECTED SAMPLES")
    print("="*60)
    
    has_scoring = all(col in adata_wmel.obs.columns for col in ['phase', 'S_score', 'G2M_score'])
    
    if args.run_scoring or not has_scoring:
        adata_wmel = score_cell_cycle_scanpy(adata_wmel, args.output, 'wMel-infected')
        if adata_wmel is None:
            print("\nERROR: Cell cycle scoring failed for wMel-infected samples")
            return
    else:
        if 'cellcycle_pseudotime' not in adata_wmel.obs.columns:
            adata_wmel = calculate_cellcycle_pseudotime(adata_wmel)
    
    # Score cell cycle for uninfected samples (if available)
    if adata_uninfected is not None:
        print("\n" + "="*60)
        print("SCORING UNINFECTED SAMPLES")
        print("="*60)
        
        has_scoring_uninf = all(col in adata_uninfected.obs.columns for col in ['phase', 'S_score', 'G2M_score'])
        
        if args.run_scoring or not has_scoring_uninf:
            adata_uninfected = score_cell_cycle_scanpy(adata_uninfected, args.output, 'Uninfected')
            if adata_uninfected is None:
                print("\nWARNING: Cell cycle scoring failed for uninfected samples")
                adata_uninfected = None
        else:
            if 'cellcycle_pseudotime' not in adata_uninfected.obs.columns:
                adata_uninfected = calculate_cellcycle_pseudotime(adata_uninfected)
    
    # Save filtered datasets if requested
    if args.save_output:
        wmel_path = os.path.join(args.output, 'JW18wMel-Ctrl_filtered.h5ad')
        print(f"\nSaving wMel-infected samples to {wmel_path}")
        adata_wmel.write(wmel_path)
        
        if adata_uninfected is not None:
            uninf_path = os.path.join(args.output, 'JW18DOX-Ctrl_filtered.h5ad')
            print(f"Saving uninfected samples to {uninf_path}")
            adata_uninfected.write(uninf_path)
    
    # Run analyses
    print("\n" + "="*60)
    print("RUNNING ANALYSES")
    print("="*60)
    
    # 1. Titer-cell cycle analysis (wMel-infected only)
    titer_results = analyze_titer_cellcycle_association(adata_wmel, args.output, 
                                                        'wMel-infected', titer_col)
    
    # 2. wMel-infected vs Uninfected comparison
    if adata_uninfected is not None:
        comparison_results = compare_cellcycle_distributions(adata_wmel, adata_uninfected, 
                                                            args.output, 'JW18_comparison')
    # Inside main(), under "RUN ANALYSES" section:
    
    # 3. Titer-Cluster Association
    cluster_results = analyze_titer_by_cluster(
        adata_wmel, 
        args.output, 
        'wMel-infected', 
        titer_col, 
        cluster_col='leiden' # Change this if your clusters are named differently
    )

    if cluster_results and cluster_results['p_kw'] < 0.05:
        de_genes = find_titer_associated_genes(adata_wmel, cluster_results)
        
        # Optional: Save the gene list to CSV
        de_genes.to_csv(os.path.join(args.output, 'titer_associated_DE_genes.csv'), index=False)
        
        # Optional: Plot the top genes
        sc.pl.rank_genes_groups(adata_wmel, n_genes=20, sharey=False, show=False)
        plt.savefig(os.path.join(args.output, 'titer_group_de_genes_plot.pdf'))

        top_genes = de_genes['names'].head(5).tolist()
        sc.pl.dotplot(adata_wmel, top_genes, groupby='leiden', standard_scale='var', save=os.path.join(args.output, 'titer_associated_top_genes_dotplot.pdf'))

    analyze_titer_phase_intersection(adata_wmel, args.output, 'JW18wMel-Ctrl')

    # Summary
    print("\n" + "="*60)
    print("ANALYSIS SUMMARY")
    print("="*60)
    print(f"\nwMel-infected samples (JW18wMel-Ctrl): {adata_wmel.n_obs} cells")
    
    if adata_uninfected is not None:
        print(f"Uninfected samples (JW18DOX-Ctrl): {adata_uninfected.n_obs} cells")
        if comparison_results:
            print(f"\nwMel-infected vs Uninfected comparison:")
            print(f"  Chi-square: p = {comparison_results['p_value']:.2e}")
            for phase, data in comparison_results['phase_tests'].items():
                print(f"  {phase}: {data['diff']:+.1f}% difference")
    
    if titer_results:
        print(f"\nTiter-cell cycle association:")
        print(f"  Spearman ρ = {titer_results['spearman_r']:.3f}, p = {titer_results['spearman_p']:.2e}")
        print(f"  Kruskal-Wallis: p = {titer_results['kw_p']:.2e}")
    
        # Add to the Summary section at the end of main():
    if cluster_results:
        print(f"\nTiter-Cluster association:")
        print(f"  Kruskal-Wallis p = {cluster_results['p_value']:.2e}")
        top_cluster = cluster_results['stats'].index[0]
        print(f"  Highest titer cluster: {top_cluster} (median: {cluster_results['stats'].loc[top_cluster, 'median']:.3f})")

    print(f"\n{'='*60}")
    print(f"All results saved to: {args.output}")
    print("="*60)


if __name__ == "__main__":
    main()