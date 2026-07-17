'''
Combined Cyclum annotation and cluster-cell cycle association analysis
'''
import cyclum 
import cyclum.models
import cyclum.tuning 
import cyclum.illustration
import scanpy as sc 
import argparse
import os
import sklearn
from sklearn.neighbors import NearestNeighbors
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chi2_contingency, kruskal

def assign_cell_cycle_stage_simple(pseudotime_flat):
    """
    Simple, gap-free cell cycle assignment
    """
    print("Step 1: Converting to circular coordinates...")
    
    # Normalize pseudotime to [0, 2π] 
    if pseudotime_flat.max() <= 1:
        angles = pseudotime_flat * 2 * np.pi
    else:
        # Normalize to [0, 2π] range
        angles = ((pseudotime_flat - pseudotime_flat.min()) / 
                 (pseudotime_flat.max() - pseudotime_flat.min())) * 2 * np.pi
    
    print("Step 2: Direct angle-based assignment...")
    
    # Simple direct assignment based on angle ranges
    phases = []
    for angle in angles:
        normalized_angle = angle % (2 * np.pi)
        
        if normalized_angle < (2 * np.pi / 3):  # 0 to 2π/3
            phases.append('g0/g1')
        elif normalized_angle < (4 * np.pi / 3):  # 2π/3 to 4π/3
            phases.append('s')
        else:  # 4π/3 to 2π
            phases.append('g2/m')
    
    print("Step 3: Light smoothing at boundaries...")
    
    # Apply light smoothing only near the boundaries
    boundary1 = 2 * np.pi / 3
    boundary2 = 4 * np.pi / 3
    boundary_width = np.pi / 12
    
    n_cells = len(angles)
    if n_cells > 10:
        nn = NearestNeighbors(n_neighbors=min(10, n_cells//10))
        circular_coords = np.column_stack([np.cos(angles), np.sin(angles)])
        nn.fit(circular_coords)
        
        smoothed_phases = phases.copy()
        changes_made = 0
        
        for i, angle in enumerate(angles):
            normalized_angle = angle % (2 * np.pi)
            
            near_boundary = (abs(normalized_angle - boundary1) < boundary_width or 
                           abs(normalized_angle - boundary2) < boundary_width or
                           abs(normalized_angle - 0) < boundary_width or
                           abs(normalized_angle - 2*np.pi) < boundary_width)
            
            if near_boundary:
                distances, indices = nn.kneighbors([circular_coords[i]])
                neighbor_indices = indices[0][1:]
                neighbor_phases = [phases[j] for j in neighbor_indices]
                
                current_phase = phases[i]
                phase_counts = {}
                for phase in neighbor_phases:
                    phase_counts[phase] = phase_counts.get(phase, 0) + 1
                
                if phase_counts:
                    most_common = max(phase_counts, key=phase_counts.get)
                    if (phase_counts[most_common] > len(neighbor_phases) * 0.7 and 
                        most_common != current_phase):
                        smoothed_phases[i] = most_common
                        changes_made += 1
        
        phases = smoothed_phases
        print(f"Boundary smoothing changed {changes_made} assignments")
    
    print("Step 4: Calculate confidence scores...")
    
    confidence_scores = np.ones(len(angles))
    
    for i, angle in enumerate(angles):
        normalized_angle = angle % (2 * np.pi)
        
        dist_to_b1 = min(abs(normalized_angle - boundary1), 2*np.pi - abs(normalized_angle - boundary1))
        dist_to_b2 = min(abs(normalized_angle - boundary2), 2*np.pi - abs(normalized_angle - boundary2))
        dist_to_start = min(normalized_angle, 2*np.pi - normalized_angle)
        
        min_dist_to_boundary = min(dist_to_b1, dist_to_b2, dist_to_start)
        confidence_scores[i] = min(1.0, min_dist_to_boundary / (boundary_width * 2))
    
    # Final validation
    phase_counts = pd.Series(phases).value_counts()
    total_cells = len(phases)
    
    print("Final phase distribution:")
    for phase in ['g0/g1', 's', 'g2/m']:
        count = phase_counts.get(phase, 0)
        percentage = (count / total_cells) * 100
        print(f"  {phase}: {count} cells ({percentage:.1f}%)")
    
    unassigned = sum(1 for phase in phases if phase not in ['g0/g1', 's', 'g2/m'])
    print(f"Unassigned cells: {unassigned} (should be 0)")
    
    return phases, confidence_scores


def run_cyclum_analysis(adata, output_dir, sample_name, force_retrain=False):
    """
    Run Cyclum on integrated data or load existing annotations
    """
    print("\n" + "="*60)
    print("CYCLUM CELL CYCLE ANNOTATION")
    print("="*60)
    
    # Check if Cyclum has already been run
    has_cyclum = all(col in adata.obs.columns for col in ['cyclum_stage', 'cyclum_pseudotime', 'cyclum_confidence'])
    
    if has_cyclum and not force_retrain:
        print("Cyclum annotations already present in adata.obs")
        print("Columns found: cyclum_stage, cyclum_pseudotime, cyclum_confidence")
        print(f"Cell cycle stages: {adata.obs['cyclum_stage'].value_counts().to_dict()}")
        
        # Just create the plots
        create_cyclum_plots(adata, output_dir, sample_name)
        return adata
    
    print("Running Cyclum on integrated dataset...")
    print(f"Dataset: {adata.n_obs} cells × {adata.n_vars} genes")
    
    # Get expression matrix
    # Use raw counts if available, otherwise use X
    if hasattr(adata, 'raw') and adata.raw is not None:
        print("Using raw counts for Cyclum")
        mtx = adata.raw.X
    else:
        print("Using X matrix for Cyclum")
        mtx = adata.X
    
    # Train model 
    print("\nTraining Cyclum model...")
    model = cyclum.tuning.CyclumAutoTune(mtx)
    model.train(mtx, epochs=800, verbose=100, rate=2e-4)
    
    # Extract pseudotime
    pseudotime = model.predict_pseudotime(mtx)
    pseudotime_flat = pseudotime.flatten()
    
    print(f"\nPseudotime shape: {pseudotime.shape}")
    print(f"Pseudotime range: {pseudotime.min():.3f} to {pseudotime.max():.3f}")
    
    # Assign cell cycle stages
    stages, confidence_scores = assign_cell_cycle_stage_simple(pseudotime_flat)
    
    # Add to adata
    adata.obs['cyclum_stage'] = stages
    adata.obs['cyclum_pseudotime'] = pseudotime_flat
    adata.obs['cyclum_confidence'] = confidence_scores
    
    print("\nAdded to adata.obs:")
    print(f"  cyclum_stage: {len(adata.obs['cyclum_stage'])} cells")
    print(f"  cyclum_pseudotime: {len(adata.obs['cyclum_pseudotime'])} cells")
    print(f"  cyclum_confidence: {len(adata.obs['cyclum_confidence'])} cells")
    
    # Create plots
    create_cyclum_plots(adata, output_dir, sample_name, model=model)
    
    return adata


def create_cyclum_plots(adata, output_dir, sample_name, model=None):
    """
    Create Cyclum diagnostic plots
    """
    os.makedirs(output_dir, exist_ok=True)
    
    pseudotime_flat = adata.obs['cyclum_pseudotime'].values
    stages = adata.obs['cyclum_stage'].values
    confidence_scores = adata.obs['cyclum_confidence'].values
    
    # Color map
    color_map = {'stage': {"g0/g1": "red", "s": "green", "g2/m": "blue"}}
    
    # 1. Circular cell cycle plot
    fig = cyclum.illustration.plot_round_distr_color(pseudotime_flat, stages, color_map['stage'])
    plt.savefig(os.path.join(output_dir, f'{sample_name}_cyclum_cell_cycle.pdf'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Elbow and bar plots (if model available)
    if model is not None:
        elbow_fig = model.show_elbow()
        plt.savefig(os.path.join(output_dir, f'{sample_name}_cyclum_elbow.pdf'), dpi=300, bbox_inches='tight')
        plt.close()
        
        bar_fig = model.show_bar()
        plt.savefig(os.path.join(output_dir, f'{sample_name}_cyclum_bar.pdf'), dpi=300, bbox_inches='tight')
        plt.close()
    
    # 3. Confidence plots
    plt.figure(figsize=(16, 4))
    
    plt.subplot(1, 4, 1)
    plt.hist(confidence_scores, bins=30, alpha=0.7, edgecolor='black')
    plt.xlabel('Confidence Score')
    plt.ylabel('Number of Cells')
    plt.title('Assignment Confidence')
    
    plt.subplot(1, 4, 2)
    colors = [color_map['stage'][stage] for stage in stages]
    plt.scatter(pseudotime_flat, confidence_scores, c=colors, alpha=0.6, s=10)
    plt.xlabel('Pseudotime')
    plt.ylabel('Confidence Score')
    plt.title('Confidence vs Pseudotime')
    
    plt.subplot(1, 4, 3)
    angles = pseudotime_flat * 2 * np.pi if pseudotime_flat.max() <= 1 else ((pseudotime_flat - pseudotime_flat.min()) / (pseudotime_flat.max() - pseudotime_flat.min())) * 2 * np.pi
    plt.hist(angles, bins=60, alpha=0.7, edgecolor='black')
    plt.xlabel('Angle (radians)')
    plt.ylabel('Number of Cells')
    plt.title('Cell Distribution Around Circle')
    plt.axvline(2*np.pi/3, color='red', linestyle='--', alpha=0.7, label='G1/S boundary')
    plt.axvline(4*np.pi/3, color='green', linestyle='--', alpha=0.7, label='S/G2M boundary')
    plt.legend()
    
    plt.subplot(1, 4, 4)
    stage_counts = pd.Series(stages).value_counts()
    stage_counts.plot(kind='bar', color=[color_map['stage'][s] for s in stage_counts.index], alpha=0.7)
    plt.xlabel('Cell Cycle Stage')
    plt.ylabel('Number of Cells')
    plt.title('Stage Distribution')
    plt.xticks(rotation=45)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{sample_name}_cyclum_confidence.pdf'), dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"\nCyclum plots saved to {output_dir}")


def analyze_cluster_cellcycle_association(adata, fig_dir, sample):
    """Test and visualize association between Leiden clusters and cell cycle"""
    
    print("\n" + "="*60)
    print("CLUSTER - CELL CYCLE ASSOCIATION ANALYSIS")
    print("="*60)
    
    # Check for required columns
    if 'leiden' not in adata.obs.columns:
        print("ERROR: No 'leiden' clustering found in adata.obs")
        return None
    
    if 'cyclum_stage' not in adata.obs.columns:
        print("ERROR: No 'cyclum_stage' found in adata.obs")
        print("Run Cyclum first with --run-cyclum flag")
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
    
    contingency = pd.crosstab(adata.obs['leiden'], adata.obs['cyclum_stage'])
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
    ax.set_xlabel('Cell Cycle Stage')
    ax.set_ylabel('Leiden Cluster')
    ax.set_title(f'Cell cycle stage distribution by cluster\nχ² = {chi2:.2f}, p = {p_value:.2e}, Cramér\'s V = {cramers_v:.3f}')
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'heatmap_cluster_cellcycle_{sample}.pdf'))
    plt.close()
    
    # 3. Stacked bar
    fig, ax = plt.subplots(figsize=(14, 7))
    contingency_norm.plot(kind='bar', stacked=True, ax=ax, width=0.8)
    ax.set_xlabel('Leiden Cluster', fontsize=12)
    ax.set_ylabel('Percentage of cells', fontsize=12)
    ax.set_title(f'Cell cycle stage composition by cluster', fontsize=14)
    ax.legend(title='Cell Cycle Stage', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks(rotation=0)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'barplot_cluster_cellcycle_{sample}.pdf'))
    plt.close()
    
    # 4. Enrichment analysis
    print("\n" + "="*60)
    print("2. CELL CYCLE ENRICHED CLUSTERS")
    print("="*60)
    
    dominant_stage = contingency_norm.idxmax(axis=1)
    max_percentage = contingency_norm.max(axis=1)
    
    print("\nDominant cell cycle stage per cluster:")
    print(f"{'Cluster':<10} {'Dominant Stage':<15} {'Percentage':<12} {'Status'}")
    print("-" * 60)
    for cluster in clusters:
        stage = dominant_stage[cluster]
        pct = max_percentage[cluster]
        if pct > 50:
            enrichment = "STRONGLY ENRICHED"
        elif pct > 40:
            enrichment = "ENRICHED"
        else:
            enrichment = "Mixed"
        print(f"{cluster:<10} {stage:<15} {pct:>6.1f}%      {enrichment}")
    
    # 5. Kruskal-Wallis test
    if 'cyclum_pseudotime' in adata.obs.columns:
        print("\n" + "="*60)
        print("3. KRUSKAL-WALLIS TEST: Pseudotime across clusters")
        print("="*60)
        
        groups = [adata.obs[adata.obs['leiden'] == cluster]['cyclum_pseudotime'].dropna().values 
                  for cluster in clusters]
        
        h_stat, p_value_kw = kruskal(*groups)
        print(f"H-statistic = {h_stat:.2f}")
        print(f"p-value = {p_value_kw:.2e}")
        print(f"\nConclusion: Pseudotime {'SIGNIFICANTLY' if p_value_kw < 0.05 else 'DOES NOT significantly'} differ across clusters")
        
        # Violin plot
        fig, ax = plt.subplots(figsize=(14, 6))
        sc.pl.violin(adata, 'cyclum_pseudotime', groupby='leiden', 
                    ax=ax, show=False, rotation=0)
        ax.set_title(f'Cell cycle pseudotime by cluster\nKruskal-Wallis H = {h_stat:.2f}, p = {p_value_kw:.2e}')
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f'violin_pseudotime_by_cluster_{sample}.pdf'))
        plt.close()
        
        # Bar plot
        pseudotime_by_cluster = adata.obs.groupby('leiden')['cyclum_pseudotime'].agg(['mean', 'std', 'median'])
        
        fig, ax = plt.subplots(figsize=(12, 6))
        ax.bar(range(len(pseudotime_by_cluster)), pseudotime_by_cluster['mean'],
               yerr=pseudotime_by_cluster['std'], color=leiden_colors, 
               alpha=0.7, capsize=5)
        ax.set_xlabel('Leiden Cluster', fontsize=12)
        ax.set_ylabel('Mean Cell Cycle Pseudotime', fontsize=12)
        ax.set_title('Mean pseudotime by cluster (with std dev)')
        ax.set_xticks(range(len(pseudotime_by_cluster)))
        ax.set_xticklabels(pseudotime_by_cluster.index)
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f'barplot_pseudotime_by_cluster_{sample}.pdf'))
        plt.close()
    else:
        p_value_kw = None
        pseudotime_by_cluster = None
    
    # 6. Confidence by cluster
    if 'cyclum_confidence' in adata.obs.columns:
        print("\n" + "="*60)
        print("4. CELL CYCLE CONFIDENCE BY CLUSTER")
        print("="*60)
        
        confidence_by_cluster = adata.obs.groupby('leiden')['cyclum_confidence'].agg(['mean', 'std', 'median'])
        print(confidence_by_cluster.round(3))
        
        fig, ax = plt.subplots(figsize=(14, 6))
        sc.pl.violin(adata, 'cyclum_confidence', groupby='leiden', 
                    ax=ax, show=False, rotation=0)
        ax.set_title('Cell cycle confidence by cluster')
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f'violin_confidence_by_cluster_{sample}.pdf'))
        plt.close()
    else:
        confidence_by_cluster = None
    
    # 7. UMAPs
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    sc.pl.umap(adata, color='leiden', ax=axes[0], show=False, 
               title='Leiden Clusters', frameon=False, legend_loc='on data')
    sc.pl.umap(adata, color='cyclum_stage', ax=axes[1], show=False,
               title='Cell Cycle Stage', frameon=False)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'umap_cluster_vs_cellcycle_{sample}.pdf'))
    plt.close()
    
    # Individual UMAPs
    sc.pl.umap(adata, color='leiden', save=f'_{sample}_leiden.pdf',
               title='Leiden Clusters', legend_loc='on data')
    sc.pl.umap(adata, color='cyclum_stage', save=f'_{sample}_cyclum_stage.pdf',
               title='Cell Cycle Stage')
    
    if 'cyclum_pseudotime' in adata.obs.columns:
        sc.pl.umap(adata, color='cyclum_pseudotime', save=f'_{sample}_cyclum_pseudotime.pdf',
                   title='Cell Cycle Pseudotime', cmap='twilight')
    
    # 8. Polar plot
    if 'cyclum_pseudotime' in adata.obs.columns:
        fig, ax = plt.subplots(figsize=(12, 12), subplot_kw=dict(projection='polar'))
        
        for i, cluster in enumerate(clusters):
            cluster_data = adata.obs[adata.obs['leiden'] == cluster]['cyclum_pseudotime']
            theta = cluster_data.values * 2 * np.pi
            r = np.ones(len(theta)) * (i + 1)
            ax.scatter(theta, r, c=[leiden_colors[i]], alpha=0.3, s=1, label=f'Cluster {cluster}')
        
        ax.set_ylim(0, len(clusters) + 1)
        ax.set_yticks(range(1, len(clusters) + 1))
        ax.set_yticklabels(clusters)
        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)
        ax.set_title('Cell cycle pseudotime distribution by cluster\n(radial = cluster, angle = pseudotime)', 
                     pad=20, fontsize=14)
        ax.legend(bbox_to_anchor=(1.3, 1.0), fontsize=10)
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f'polar_pseudotime_by_cluster_{sample}.pdf'))
        plt.close()
    
    # 9. Summary table
    print("\n" + "="*60)
    print("5. COMPREHENSIVE SUMMARY TABLE")
    print("="*60)
    
    summary_data = {
        'Cluster': clusters,
        'N_Cells': [contingency.loc[c].sum() for c in clusters],
        'Dominant_Stage': [dominant_stage[c] for c in clusters],
        'Stage_Percentage': [max_percentage[c] for c in clusters],
    }
    
    if pseudotime_by_cluster is not None:
        summary_data['Mean_Pseudotime'] = [pseudotime_by_cluster.loc[c, 'mean'] for c in clusters]
        summary_data['Median_Pseudotime'] = [pseudotime_by_cluster.loc[c, 'median'] for c in clusters]
        summary_data['Std_Pseudotime'] = [pseudotime_by_cluster.loc[c, 'std'] for c in clusters]
    
    if confidence_by_cluster is not None:
        summary_data['Mean_Confidence'] = [confidence_by_cluster.loc[c, 'mean'] for c in clusters]
        summary_data['Median_Confidence'] = [confidence_by_cluster.loc[c, 'median'] for c in clusters]
    
    summary_df = pd.DataFrame(summary_data)
    print(summary_df.to_string(index=False))
    
    # Save outputs
    summary_df.to_csv(os.path.join(fig_dir, f'cluster_cellcycle_summary_{sample}.csv'), index=False)
    contingency.to_csv(os.path.join(fig_dir, f'contingency_table_counts_{sample}.csv'))
    contingency_norm.to_csv(os.path.join(fig_dir, f'contingency_table_percentages_{sample}.csv'))
    
    # Statistical results
    stats_results = {
        'Test': ['Chi-square', 'Cramers V'],
        'Statistic': [chi2, cramers_v],
        'P-value': [p_value, np.nan],
        'Interpretation': [
            'Significant association' if p_value < 0.05 else 'No significant association',
            f'{effect.capitalize()} effect size'
        ]
    }
    
    if p_value_kw is not None:
        stats_results['Test'].append('Kruskal-Wallis')
        stats_results['Statistic'].append(h_stat)
        stats_results['P-value'].append(p_value_kw)
        stats_results['Interpretation'].append(
            'Pseudotime differs across clusters' if p_value_kw < 0.05 else 'No difference in pseudotime'
        )
    
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
        'kw_pvalue': p_value_kw,
        'summary': summary_df
    }


def main():
    parser = argparse.ArgumentParser(
        description='Cyclum annotation and cluster-cell cycle association analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Run Cyclum + analysis
  python cyclum_cluster_analysis.py \\
      --input integrated.h5ad \\
      --output cluster_cellcycle_analysis \\
      --sample all_conditions \\
      --run-cyclum
  
  # Just run analysis (Cyclum already done)
  python cyclum_cluster_analysis.py \\
      --input integrated.h5ad \\
      --output cluster_cellcycle_analysis \\
      --sample all_conditions
  
  # Force re-run Cyclum
  python cyclum_cluster_analysis.py \\
      --input integrated.h5ad \\
      --output cluster_cellcycle_analysis \\
      --sample all_conditions \\
      --run-cyclum \\
      --force-retrain
        '''
    )
    
    parser.add_argument('--input', '-i', required=True, type=str,
                        help='Path to integrated h5ad file')
    parser.add_argument('--output', '-o', type=str, default='cluster_cellcycle_analysis',
                        help='Output directory (default: cluster_cellcycle_analysis)')
    parser.add_argument('--sample', '-s', type=str, default='sample',
                        help='Sample name for output files (default: sample)')
    parser.add_argument('--run-cyclum', action='store_true',
                        help='Run Cyclum annotation (skip if already in adata)')
    parser.add_argument('--force-retrain', action='store_true',
                        help='Force re-train Cyclum even if annotations exist')
    parser.add_argument('--save-output', action='store_true',
                        help='Save updated h5ad with Cyclum annotations')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    
    # Load data
    print(f"Loading data from {args.input}...")
    adata = sc.read_h5ad(args.input)
    
    print(f"\nLoaded AnnData object:")
    print(f"  Cells: {adata.n_obs}")
    print(f"  Genes: {adata.n_vars}")
    print(f"  Observations: {list(adata.obs.columns)}")
    
    # Run or load Cyclum
    if args.run_cyclum or args.force_retrain:
        adata = run_cyclum_analysis(adata, args.output, args.sample, force_retrain=args.force_retrain)
        
        # Save if requested
        if args.save_output:
            output_path = args.input.replace('.h5ad', '_with_cyclum.h5ad')
            print(f"\nSaving updated h5ad to {output_path}")
            adata.write(output_path)
    
    # Run cluster-cell cycle association analysis
    results = analyze_cluster_cellcycle_association(adata, args.output, args.sample)
    
    if results:
        print(f"\n{'='*60}")
        print("KEY FINDINGS")
        print(f"{'='*60}")
        print(f"Chi-square test: χ² = {results['chi2']:.2f}, p = {results['chi2_pvalue']:.2e}")
        print(f"Effect size (Cramér's V): {results['cramers_v']:.3f}")
        if results['kw_pvalue'] is not None:
            print(f"Kruskal-Wallis test: p = {results['kw_pvalue']:.2e}")


if __name__ == "__main__":
    main()