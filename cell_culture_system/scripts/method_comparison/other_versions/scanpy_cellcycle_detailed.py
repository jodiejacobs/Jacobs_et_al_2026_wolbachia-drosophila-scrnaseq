'''
Cell cycle annotation using Scanpy with Drosophila-specific markers
Includes continuous position and marker-based detailed sub-phases
'''
import scanpy as sc 
import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chi2_contingency, kruskal
import warnings
warnings.filterwarnings('ignore')

# [Keep all the FlyBase gene definitions from before - S_GENES_FBGN, G2M_GENES_FBGN, etc.]
# ... [Copy the entire gene definition section from the previous script]

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

SUBPHASE_MARKERS = {
    'G1_early': ['FBgn0010382', 'FBgn0011766'],
    'G1_late': ['FBgn0005655', 'FBgn0020651', 'FBgn0015929'],
    'S_early': ['FBgn0015806', 'FBgn0020369', 'FBgn0034898'],
    'S_late': ['FBgn0011230', 'FBgn0019624'],
    'G2_early': ['FBgn0010114'],
    'G2_late': ['FBgn0010113', 'FBgn0011739'],
    'M_early': ['FBgn0003525', 'FBgn0003124', 'FBgn0025564'],
    'M_late': ['FBgn0010309', 'FBgn0002610', 'FBgn0025948'],
}

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


def assign_continuous_cellcycle_position(adata):
    """
    Assign continuous cell cycle position using the basic 3-phase assignments
    This maps G1 -> 0-0.33, S -> 0.33-0.66, G2M -> 0.66-1.0 with jitter
    """
    print("\n" + "="*60)
    print("CALCULATING CONTINUOUS CELL CYCLE POSITION")
    print("="*60)
    
    s_score = adata.obs['S_score'].values
    g2m_score = adata.obs['G2M_score'].values
    phase = adata.obs['phase'].values
    
    # Create position based on phase and scores
    positions = np.zeros(len(phase))
    
    for i, (p, s, g) in enumerate(zip(phase, s_score, g2m_score)):
        if p == 'G1':
            # G1 occupies 0.0 to 0.33
            # Use scores to position within G1
            # Lower scores = earlier in G1
            base_pos = 0.165  # middle of G1
            jitter = (s + g) * 0.05  # small adjustment
            positions[i] = np.clip(base_pos + jitter, 0.0, 0.33)
        
        elif p == 'S':
            # S occupies 0.33 to 0.66
            # Use S score to position within S
            base_pos = 0.495  # middle of S
            # Normalize S score within S phase cells
            s_in_s = s_score[phase == 'S']
            if len(s_in_s) > 0:
                s_norm = (s - s_in_s.min()) / (s_in_s.max() - s_in_s.min() + 1e-10)
                jitter = (s_norm - 0.5) * 0.25  # spread across S phase
            else:
                jitter = 0
            positions[i] = np.clip(base_pos + jitter, 0.33, 0.66)
        
        else:  # G2M
            # G2M occupies 0.66 to 1.0
            # Use G2M score to position within G2M
            base_pos = 0.83  # middle of G2M
            # Normalize G2M score within G2M phase cells
            g_in_g2m = g2m_score[phase == 'G2M']
            if len(g_in_g2m) > 0:
                g_norm = (g - g_in_g2m.min()) / (g_in_g2m.max() - g_in_g2m.min() + 1e-10)
                jitter = (g_norm - 0.5) * 0.25  # spread across G2M phase
            else:
                jitter = 0
            positions[i] = np.clip(base_pos + jitter, 0.66, 1.0)
    
    adata.obs['cellcycle_position'] = positions
    
    print(f"Cell cycle position range: {positions.min():.3f} to {positions.max():.3f}")
    print(f"Mean position: {positions.mean():.3f}")
    print(f"Positions by phase:")
    for p in ['G1', 'S', 'G2M']:
        phase_pos = positions[phase == p]
        if len(phase_pos) > 0:
            print(f"  {p}: {phase_pos.min():.3f} to {phase_pos.max():.3f} (mean: {phase_pos.mean():.3f})")
    
    return adata


def assign_subphases_by_markers(adata):
    """
    Assign detailed sub-phases using specific marker genes
    Fixed version that handles zero expression properly
    """
    print("\n" + "="*60)
    print("ASSIGNING SUB-PHASES USING MARKER GENES")
    print("="*60)
    
    # Calculate expression scores for each sub-phase
    subphase_scores = {}
    genes_found = {}
    
    for phase, genes in SUBPHASE_MARKERS.items():
        present_genes = [g for g in genes if g in adata.var_names]
        genes_found[phase] = len(present_genes)
        
        if present_genes:
            expr_data = adata[:, present_genes].X
            if hasattr(expr_data, 'toarray'):
                expr_data = expr_data.toarray()
            # Use mean expression across marker genes
            subphase_scores[phase] = expr_data.mean(axis=1)
        else:
            # If no genes found, use very low scores
            subphase_scores[phase] = np.full(adata.n_obs, -999.0)
    
    print("\nMarker genes found per sub-phase:")
    for phase in ['G1_early', 'G1_late', 'S_early', 'S_late',
                  'G2_early', 'G2_late', 'M_early', 'M_late']:
        n_genes = genes_found.get(phase, 0)
        total = len(SUBPHASE_MARKERS[phase])
        if n_genes > 0:
            symbols = [FBGN_TO_SYMBOL.get(g, g) for g in SUBPHASE_MARKERS[phase] if g in adata.var_names]
            print(f"  {phase}: {n_genes}/{total} genes - {', '.join(symbols)}")
        else:
            print(f"  {phase}: {n_genes}/{total} genes - NONE FOUND")
    
    # Assign phase based on highest score
    subphase_df = pd.DataFrame(subphase_scores)
    
    # Get the phase with max score
    adata.obs['phase_marker'] = subphase_df.idxmax(axis=1)
    adata.obs['phase_marker_score'] = subphase_df.max(axis=1)
    
    # Store individual scores
    for phase, scores in subphase_scores.items():
        adata.obs[f'{phase}_score'] = scores
    
    # Check if we got valid assignments
    valid_assignments = (adata.obs['phase_marker_score'] > -900).sum()
    print(f"\nCells with valid marker-based assignments: {valid_assignments}/{adata.n_obs}")
    
    if valid_assignments == 0:
        print("\nWARNING: No valid marker-based assignments!")
        print("This might be because too few marker genes were found.")
        print("Falling back to position-based phases only.")
        # Copy position-based to marker-based
        if 'phase_position' in adata.obs.columns:
            adata.obs['phase_marker'] = adata.obs['phase_position']
            adata.obs['phase_marker_score'] = np.ones(adata.n_obs)
    
    print("\nMarker-based sub-phase distribution:")
    phase_order = ['G1_early', 'G1_late', 'S_early', 'S_late',
                   'G2_early', 'G2_late', 'M_early', 'M_late']
    phase_counts = adata.obs['phase_marker'].value_counts()
    for phase in phase_order:
        count = phase_counts.get(phase, 0)
        pct = (count / adata.n_obs) * 100
        print(f"  {phase}: {count} cells ({pct:.1f}%)")
    
    return adata


def assign_position_based_phases(adata):
    """
    Assign 8-stage phases based on continuous position
    """
    print("\n" + "="*60)
    print("ASSIGNING POSITION-BASED SUB-PHASES")
    print("="*60)
    
    if 'cellcycle_position' not in adata.obs.columns:
        print("ERROR: cellcycle_position not found.")
        return adata
    
    cc_position = adata.obs['cellcycle_position'].values
    
    # 8-stage assignment
    phases_position = []
    for pos in cc_position:
        if pos < 0.165:
            phase = 'G1_early'
        elif pos < 0.33:
            phase = 'G1_late'
        elif pos < 0.415:
            phase = 'S_early'
        elif pos < 0.66:
            phase = 'S_late'
        elif pos < 0.745:
            phase = 'G2_early'
        elif pos < 0.83:
            phase = 'G2_late'
        elif pos < 0.915:
            phase = 'M_early'
        else:
            phase = 'M_late'
        phases_position.append(phase)
    
    adata.obs['phase_position'] = phases_position
    
    print("\nPosition-based sub-phase distribution:")
    phase_order = ['G1_early', 'G1_late', 'S_early', 'S_late',
                   'G2_early', 'G2_late', 'M_early', 'M_late']
    phase_counts = adata.obs['phase_position'].value_counts()
    for phase in phase_order:
        count = phase_counts.get(phase, 0)
        pct = (count / adata.n_obs) * 100
        print(f"  {phase}: {count} cells ({pct:.1f}%)")
    
    return adata


def create_consensus_phases(adata):
    """
    Create consensus phases combining position and marker-based assignments
    """
    print("\n" + "="*60)
    print("CREATING CONSENSUS SUB-PHASES")
    print("="*60)
    
    if 'phase_position' not in adata.obs.columns or 'phase_marker' not in adata.obs.columns:
        print("WARNING: Missing required columns for consensus")
        return adata
    
    # Check if marker-based is valid
    has_valid_markers = (adata.obs['phase_marker_score'] > -900).sum() > 0
    
    if not has_valid_markers:
        print("Marker-based phases not valid, using position-based only")
        adata.obs['phase_consensus'] = adata.obs['phase_position']
        return adata
    
    # Create consensus
    consensus = []
    agreement_count = 0
    
    for idx, row in adata.obs.iterrows():
        pos_phase = row['phase_position']
        marker_phase = row['phase_marker']
        marker_score = row['phase_marker_score']
        
        # If marker score is invalid, use position
        if marker_score < -900:
            consensus.append(pos_phase)
        elif pos_phase == marker_phase:
            consensus.append(pos_phase)
            agreement_count += 1
        else:
            # Use whichever has more confidence
            # Higher marker score = more confident
            if marker_score > 0.5:
                consensus.append(marker_phase)
            else:
                consensus.append(pos_phase)
    
    adata.obs['phase_consensus'] = consensus
    
    agreement_pct = (agreement_count / adata.n_obs) * 100
    print(f"\nAgreement between methods: {agreement_count}/{adata.n_obs} cells ({agreement_pct:.1f}%)")
    
    print("\nConsensus sub-phase distribution:")
    phase_order = ['G1_early', 'G1_late', 'S_early', 'S_late',
                   'G2_early', 'G2_late', 'M_early', 'M_late']
    phase_counts = adata.obs['phase_consensus'].value_counts()
    for phase in phase_order:
        count = phase_counts.get(phase, 0)
        pct = (count / adata.n_obs) * 100
        print(f"  {phase}: {count} cells ({pct:.1f}%)")
    
    return adata


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
        print(f"  Present: {', '.join([FBGN_TO_SYMBOL.get(g, g) for g in s_genes_present[:10]])}" + 
              (f"... (+{len(s_genes_present)-10} more)" if len(s_genes_present) > 10 else ""))
    
    print(f"\nG2/M phase genes: {len(g2m_genes_present)}/{len(g2m_genes)} found")
    if g2m_genes_present:
        print(f"  Present: {', '.join([FBGN_TO_SYMBOL.get(g, g) for g in g2m_genes_present[:10]])}" + 
              (f"... (+{len(g2m_genes_present)-10} more)" if len(g2m_genes_present) > 10 else ""))
    
    if len(s_genes_present) == 0 and len(g2m_genes_present) == 0:
        print("\nERROR: No cell cycle genes found!")
        return None
    
    # Score genes
    print("\nScoring cell cycle phases...")
    sc.tl.score_genes_cell_cycle(
        adata, 
        s_genes=s_genes_present, 
        g2m_genes=g2m_genes_present
    )
    
    # Print distribution
    print("\nBasic cell cycle phase distribution:")
    phase_counts = adata.obs['phase'].value_counts()
    for phase in ['G1', 'S', 'G2M']:
        count = phase_counts.get(phase, 0)
        pct = (count / adata.n_obs) * 100
        print(f"  {phase}: {count} cells ({pct:.1f}%)")
    
    # Assign continuous position (FIXED VERSION)
    adata = assign_continuous_cellcycle_position(adata)
    
    # Assign position-based phases
    adata = assign_position_based_phases(adata)
    
    # Assign marker-based phases
    adata = assign_subphases_by_markers(adata)
    
    # Create consensus phases
    adata = create_consensus_phases(adata)
    
    # Create plots (with error handling)
    try:
        create_basic_plots(adata, output_dir, sample_name, s_genes_present, g2m_genes_present)
    except Exception as e:
        print(f"\nWarning: Could not create all plots: {e}")
        print("Continuing anyway...")
    
    return adata


def create_basic_plots(adata, output_dir, sample_name, s_genes, g2m_genes):
    """
    Create basic diagnostic plots only (simplified version)
    """
    os.makedirs(output_dir, exist_ok=True)
    
    colors_basic = {'G1': '#FF6B6B', 'S': '#4ECDC4', 'G2M': '#45B7D1'}
    colors_8stage = {
        'G1_early': '#E74C3C', 'G1_late': '#EC7063',
        'S_early': '#3498DB', 'S_late': '#5DADE2',
        'G2_early': '#F39C12', 'G2_late': '#F8C471',
        'M_early': '#9B59B6', 'M_late': '#BB8FCE'
    }
    
    # 1. Basic overview
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Basic phase
    ax = axes[0, 0]
    phase_counts = adata.obs['phase'].value_counts()
    phase_counts.plot(kind='bar', ax=ax, color=[colors_basic.get(p, 'gray') for p in phase_counts.index])
    ax.set_title('Basic Cell Cycle Phases (3 stages)')
    ax.set_ylabel('Number of Cells')
    
    # S vs G2M scores
    ax = axes[0, 1]
    for phase, color in colors_basic.items():
        mask = adata.obs['phase'] == phase
        ax.scatter(adata.obs.loc[mask, 'S_score'], 
                  adata.obs.loc[mask, 'G2M_score'],
                  c=color, label=phase, alpha=0.5, s=5)
    ax.set_xlabel('S Score')
    ax.set_ylabel('G2M Score')
    ax.set_title('Cell Cycle Scores')
    ax.legend()
    ax.axhline(0, color='k', linestyle='--', alpha=0.3)
    ax.axvline(0, color='k', linestyle='--', alpha=0.3)
    
    # Position distribution
    ax = axes[1, 0]
    ax.hist(adata.obs['cellcycle_position'], bins=50, alpha=0.7, edgecolor='black')
    ax.set_xlabel('Cell Cycle Position')
    ax.set_ylabel('Number of Cells')
    ax.set_title('Continuous Position Distribution')
    
    # Detailed phases
    ax = axes[1, 1]
    if 'phase_consensus' in adata.obs.columns:
        phase_col = 'phase_consensus'
    else:
        phase_col = 'phase_position'
    
    phase_counts_det = adata.obs[phase_col].value_counts()
    phase_order = ['G1_early', 'G1_late', 'S_early', 'S_late',
                   'G2_early', 'G2_late', 'M_early', 'M_late']
    phase_counts_det = phase_counts_det.reindex([p for p in phase_order if p in phase_counts_det.index])
    phase_counts_det.plot(kind='bar', ax=ax,
                         color=[colors_8stage.get(p, 'gray') for p in phase_counts_det.index])
    ax.set_title('Detailed Sub-Phases (8 stages)')
    ax.set_ylabel('Number of Cells')
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{sample_name}_cellcycle_overview.pdf'),
                dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Circular plot
    fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(projection='polar'))
    
    theta = adata.obs['cellcycle_position'].values * 2 * np.pi
    phase_colors = [colors_8stage.get(p, 'gray') for p in adata.obs[phase_col]]
    
    ax.scatter(theta, np.ones(len(theta)), c=phase_colors, alpha=0.3, s=1)
    ax.set_ylim(0, 1.5)
    ax.set_yticks([])
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_title('Cell Cycle Position (Circular)', pad=20, fontsize=14)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{sample_name}_circular_cellcycle.pdf'),
                dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"\nCell cycle plots saved to {output_dir}")


# [Keep the analyze_cluster_cellcycle_association function from before - it's fine]
# [Keep the main() function from before - it's fine]

def analyze_cluster_cellcycle_association(adata, fig_dir, sample, phase_column='phase_consensus'):
    """Test and visualize association between Leiden clusters and cell cycle"""
    
    print("\n" + "="*60)
    print(f"CLUSTER - CELL CYCLE ASSOCIATION ANALYSIS ({phase_column})")
    print("="*60)
    
    if 'leiden' not in adata.obs.columns:
        print("ERROR: No 'leiden' clustering found")
        return None
    
    if phase_column not in adata.obs.columns:
        print(f"ERROR: '{phase_column}' not found")
        # Try alternatives
        if 'phase_position' in adata.obs.columns:
            phase_column = 'phase_position'
            print(f"Using '{phase_column}' instead")
        else:
            return None
    
    sc.settings.figdir = fig_dir
    
    clusters = sorted(adata.obs['leiden'].unique())
    contingency = pd.crosstab(adata.obs['leiden'], adata.obs[phase_column])
    chi2, p_value, dof, expected = chi2_contingency(contingency)
    
    n = contingency.sum().sum()
    cramers_v = np.sqrt(chi2 / (n * (min(contingency.shape) - 1)))
    
    print(f"\nχ² = {chi2:.2f}, p = {p_value:.2e}")
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
    
    # Summary
    summary_df = pd.DataFrame({
        'Cluster': clusters,
        'N_Cells': [contingency.loc[c].sum() for c in clusters],
        'Dominant_Phase': contingency_norm.idxmax(axis=1),
        'Phase_Percentage': contingency_norm.max(axis=1),
    })
    
    summary_df.to_csv(os.path.join(fig_dir, f'cluster_cellcycle_summary_{sample}.csv'), index=False)
    
    print("\n" + "="*60)
    print("ANALYSIS COMPLETE")
    print("="*60)
    
    return {'chi2': chi2, 'chi2_pvalue': p_value, 'cramers_v': cramers_v, 'summary': summary_df}


def main():
    parser = argparse.ArgumentParser(
        description='Detailed cell cycle annotation',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('--input', '-i', required=True, help='Path to h5ad file')
    parser.add_argument('--output', '-o', default='cellcycle_analysis', help='Output directory')
    parser.add_argument('--sample', '-s', default='sample', help='Sample name')
    parser.add_argument('--run-scoring', action='store_true', help='Run cell cycle scoring')
    parser.add_argument('--save-output', action='store_true', help='Save updated h5ad')
    parser.add_argument('--phase-column', default='phase_consensus',
                        choices=['phase', 'phase_position', 'phase_marker', 'phase_consensus'],
                        help='Phase column for cluster analysis')
    
    args = parser.parse_args()
    
    os.makedirs(args.output, exist_ok=True)
    
    print(f"Loading data from {args.input}...")
    adata = sc.read_h5ad(args.input)
    
    print(f"\nLoaded: {adata.n_obs} cells × {adata.n_vars} genes")
    
    has_scoring = 'phase' in adata.obs.columns
    has_detailed = 'phase_consensus' in adata.obs.columns or 'phase_position' in adata.obs.columns
    
    if args.run_scoring or not (has_scoring and has_detailed):
        result = score_cell_cycle_scanpy(adata, args.output, args.sample)
        if result is None:
            print("\nERROR: Scoring failed")
            return
        adata = result
        
        if args.save_output:
            output_path = args.input.replace('.h5ad', '_with_detailed_cellcycle.h5ad')
            print(f"\nSaving to {output_path}")
            adata.write(output_path)
    else:
        print("\nUsing existing annotations")
    
    results = analyze_cluster_cellcycle_association(adata, args.output, args.sample, args.phase_column)
    
    if results:
        print(f"\n{'='*60}")
        print("KEY FINDINGS")
        print(f"{'='*60}")
        print(f"χ² = {results['chi2']:.2f}, p = {results['chi2_pvalue']:.2e}")
        print(f"Cramér's V = {results['cramers_v']:.3f}")


if __name__ == "__main__":
    main()