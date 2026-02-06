'''
Validate Cyclum cell cycle predictions using known Drosophila cell cycle genes
Reference: https://www.sdbonline.org/sites/fly/aignfam/cellcycl.htm
'''
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os
from scipy.stats import mannwhitneyu, kruskal

# Known Drosophila cell cycle genes by phase
CELL_CYCLE_GENES = {
    'G1/S': [
        'E2f1', 'E2f2',  # E2F transcription factors
        'CycE',  # Cyclin E
        'Cdk2',  # Cyclin-dependent kinase 2
        'Dp',  # DP transcription factor
        'Rbf', 'Rbf2',  # Retinoblastoma-family proteins
        'Mcm2', 'Mcm3', 'Mcm5', 'Mcm6', 'Mcm7',  # MCM complex (DNA replication)
        'Orc1', 'Orc2', 'Orc6',  # Origin recognition complex
        'Pcna',  # PCNA (replication)
        'RnrL', 'RnrS',  # Ribonucleotide reductase
        'Rrp1',  # Replication protein
    ],
    'S': [
        'Pcna',  # PCNA continues in S
        'RPA1', 'RPA2',  # Replication protein A
        'pol-alpha1', 'DNApol-alpha60',  # DNA polymerase alpha
        'DNApol-delta',  # DNA polymerase delta
        'RnrL', 'RnrS',  # Continue in S phase
        'Mcm2', 'Mcm3', 'Mcm5', 'Mcm6', 'Mcm7',  # MCM continues
    ],
    'G2': [
        'CycA',  # Cyclin A (peaks in G2)
        'Cdk1',  # Cdk1 builds up in G2
        'Wee1',  # Wee1 kinase (G2 checkpoint)
        'myt',  # Myt1 kinase
    ],
    'G2/M': [
        'CycA',  # Cyclin A
        'CycB', 'CycB3',  # Cyclin B
        'Cdk1',  # Cdk1/Cdc2
        'stg',  # String/Cdc25 (triggers mitosis)
        'polo',  # Polo kinase
        'aurA', 'aurB',  # Aurora kinases
        'Nek2',  # NIMA-related kinase
        'PBl',  # Polo-binding protein
    ],
    'M': [
        'CycB', 'CycB3',  # Cyclin B peaks in M
        'Cdk1',  # Active Cdk1
        'polo',  # Polo kinase
        'aurA', 'aurB',  # Aurora A and B
        'BubR1', 'Mad2',  # Spindle checkpoint
        'Cdc20',  # APC/C activator
        'APC2', 'APC10',  # Anaphase-promoting complex
    ],
    'M/G1': [
        'fzr',  # Fizzy-related (APC/C-Cdh1, degrades mitotic cyclins)
        'Cdc20',  # Transitions M to G1
        'APC2', 'APC10',  # APC/C activity continues
    ]
}

# Consolidate into major phases
MAJOR_PHASE_GENES = {
    'g0/g1': [
        'E2f1', 'E2f2', 'CycE', 'Cdk2', 'Dp', 'Rbf', 'Rbf2',
        'Mcm2', 'Mcm3', 'Mcm5', 'Mcm6', 'Mcm7',
        'Orc1', 'Orc2', 'Orc6', 'Pcna', 'RnrL', 'RnrS', 'Rrp1',
        'fzr'  # Also active at G1
    ],
    's': [
        'Pcna', 'RPA1', 'RPA2', 'pol-alpha1', 'DNApol-alpha60',
        'DNApol-delta', 'RnrL', 'RnrS',
        'Mcm2', 'Mcm3', 'Mcm5', 'Mcm6', 'Mcm7'
    ],
    'g2/m': [
        'CycA', 'CycB', 'CycB3', 'Cdk1', 'stg', 'polo',
        'aurA', 'aurB', 'Nek2', 'PBl', 'Wee1', 'myt',
        'BubR1', 'Mad2', 'Cdc20', 'APC2', 'APC10'
    ]
}


def validate_cyclum_with_markers(adata, output_dir, sample_name):
    """
    Validate Cyclum predictions using known cell cycle marker genes
    """
    print("\n" + "="*70)
    print("VALIDATING CYCLUM PREDICTIONS WITH CELL CYCLE MARKER GENES")
    print("="*70)
    
    os.makedirs(output_dir, exist_ok=True)
    sc.settings.figdir = output_dir
    
    # Check required columns
    if 'cyclum_stage' not in adata.obs.columns:
        print("ERROR: cyclum_stage not found in adata.obs")
        return None
    
    # Find which marker genes are present
    all_markers = set()
    for phase_genes in MAJOR_PHASE_GENES.values():
        all_markers.update(phase_genes)
    
    present_markers = [g for g in all_markers if g in adata.var_names]
    missing_markers = all_markers - set(present_markers)
    
    print(f"\nMarker genes present: {len(present_markers)}/{len(all_markers)}")
    print(f"Missing markers: {len(missing_markers)}")
    
    if len(present_markers) < 5:
        print("WARNING: Too few marker genes found. Results may not be reliable.")
    
    # Organize present markers by phase
    phase_markers_present = {}
    for phase, genes in MAJOR_PHASE_GENES.items():
        present = [g for g in genes if g in adata.var_names]
        if present:
            phase_markers_present[phase] = present
    
    print("\nMarker genes found per phase:")
    for phase, genes in phase_markers_present.items():
        print(f"  {phase}: {len(genes)} genes - {', '.join(genes[:5])}" + 
              (f"... (+{len(genes)-5} more)" if len(genes) > 5 else ""))
    
    if missing_markers:
        print(f"\nMissing markers (sample): {list(missing_markers)[:10]}")
    
    # 1. Calculate mean expression per cyclum stage
    print("\n" + "="*70)
    print("CALCULATING EXPRESSION BY CYCLUM STAGE")
    print("="*70)
    
    stage_expression = {}
    for phase, genes in phase_markers_present.items():
        phase_expr = pd.DataFrame(index=adata.obs['cyclum_stage'].unique())
        
        for gene in genes:
            if gene in adata.var_names:
                gene_expr = pd.DataFrame({
                    'stage': adata.obs['cyclum_stage'],
                    'expression': adata[:, gene].X.toarray().flatten() if hasattr(adata.X, 'toarray') else adata[:, gene].X.flatten()
                })
                mean_expr = gene_expr.groupby('stage')['expression'].mean()
                phase_expr[gene] = mean_expr
        
        stage_expression[phase] = phase_expr
    
    # 2. Heatmap of all marker genes by cyclum stage
    print("\nCreating heatmap of marker expression by stage...")
    
    # Prepare data for heatmap
    expr_matrix = []
    gene_labels = []
    phase_labels = []
    
    for phase, genes in phase_markers_present.items():
        for gene in genes:
            if gene in adata.var_names:
                gene_expr = pd.DataFrame({
                    'stage': adata.obs['cyclum_stage'],
                    'expression': adata[:, gene].X.toarray().flatten() if hasattr(adata.X, 'toarray') else adata[:, gene].X.flatten()
                })
                mean_expr = gene_expr.groupby('stage')['expression'].mean()
                expr_matrix.append(mean_expr.values)
                gene_labels.append(gene)
                phase_labels.append(phase)
    
    expr_df = pd.DataFrame(expr_matrix, 
                           columns=sorted(adata.obs['cyclum_stage'].unique()),
                           index=gene_labels)
    
    # Normalize by row (z-score)
    expr_df_norm = expr_df.apply(lambda x: (x - x.mean()) / (x.std() + 1e-10), axis=1)
    
    # Create heatmap with phase annotations
    fig, ax = plt.subplots(figsize=(8, max(12, len(gene_labels) * 0.3)))
    
    # Create color map for phases
    phase_colors = {'g0/g1': '#FF6B6B', 's': '#4ECDC4', 'g2/m': '#45B7D1'}
    row_colors = [phase_colors[p] for p in phase_labels]
    
    sns.heatmap(expr_df_norm, cmap='RdBu_r', center=0, 
                cbar_kws={'label': 'Z-score normalized expression'},
                yticklabels=gene_labels, xticklabels=expr_df.columns,
                ax=ax, vmin=-2, vmax=2)
    
    ax.set_xlabel('Cyclum Stage', fontsize=12)
    ax.set_ylabel('Cell Cycle Marker Genes', fontsize=12)
    ax.set_title(f'Cell cycle marker expression by Cyclum stage\n{sample_name}', 
                 fontsize=14, pad=20)
    
    # Add phase color bar on the left
    for i, (gene, phase_color) in enumerate(zip(gene_labels, row_colors)):
        ax.add_patch(plt.Rectangle((-0.5, i), 0.3, 1, 
                                   color=phase_color, clip_on=False, 
                                   transform=ax.get_yaxis_transform()))
    
    # Add legend for phases
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=color, label=phase.upper()) 
                      for phase, color in phase_colors.items()]
    ax.legend(handles=legend_elements, loc='upper left', 
             bbox_to_anchor=(1.15, 1), title='Expected Phase')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'heatmap_markers_by_stage_{sample_name}.pdf'), 
                bbox_inches='tight', dpi=300)
    plt.close()
    
    # 3. Calculate validation scores
    print("\n" + "="*70)
    print("VALIDATION SCORING")
    print("="*70)
    
    validation_scores = {}
    
    for expected_phase, genes in phase_markers_present.items():
        if not genes:
            continue
            
        # Calculate mean expression of these genes in each stage
        phase_means = {}
        for stage in adata.obs['cyclum_stage'].unique():
            stage_cells = adata.obs['cyclum_stage'] == stage
            stage_expr = []
            
            for gene in genes:
                if gene in adata.var_names:
                    gene_vals = adata[stage_cells, gene].X.toarray().flatten() if hasattr(adata.X, 'toarray') else adata[stage_cells, gene].X.flatten()
                    stage_expr.append(gene_vals.mean())
            
            phase_means[stage] = np.mean(stage_expr) if stage_expr else 0
        
        # Check if expression is highest in expected stage
        max_stage = max(phase_means, key=phase_means.get)
        is_correct = (max_stage == expected_phase)
        
        validation_scores[expected_phase] = {
            'mean_expression_by_stage': phase_means,
            'peak_stage': max_stage,
            'expected_stage': expected_phase,
            'is_correct': is_correct,
            'n_genes': len(genes)
        }
        
        print(f"\n{expected_phase.upper()} markers ({len(genes)} genes):")
        print(f"  Peak expression in: {max_stage}")
        print(f"  Expected peak in: {expected_phase}")
        print(f"  Match: {'✓' if is_correct else '✗'}")
        print(f"  Expression by stage: {phase_means}")
    
    # 4. Violin plots for key marker genes
    print("\nCreating violin plots for key markers...")
    
    key_markers = {
        'g0/g1': ['CycE', 'E2f1', 'Pcna'],
        's': ['Pcna', 'RnrL', 'Mcm2'],
        'g2/m': ['CycB', 'CycA', 'polo', 'aurA', 'stg']
    }
    
    n_plots = sum(len([g for g in genes if g in adata.var_names]) 
                  for genes in key_markers.values())
    
    if n_plots > 0:
        fig, axes = plt.subplots(3, max(3, n_plots//3 + 1), 
                                figsize=(15, 10))
        axes = axes.flatten()
        
        plot_idx = 0
        for phase, genes in key_markers.items():
            for gene in genes:
                if gene in adata.var_names and plot_idx < len(axes):
                    ax = axes[plot_idx]
                    
                    # Prepare data
                    gene_data = pd.DataFrame({
                        'stage': adata.obs['cyclum_stage'],
                        'expression': adata[:, gene].X.toarray().flatten() if hasattr(adata.X, 'toarray') else adata[:, gene].X.flatten()
                    })
                    
                    sns.violinplot(data=gene_data, x='stage', y='expression', 
                                 ax=ax, palette='Set2')
                    ax.set_title(f'{gene} (expected: {phase})', fontsize=10, fontweight='bold')
                    ax.set_xlabel('')
                    ax.set_ylabel('Expression')
                    
                    plot_idx += 1
        
        # Hide unused subplots
        for idx in range(plot_idx, len(axes)):
            axes[idx].axis('off')
        
        plt.suptitle(f'Key cell cycle marker expression by Cyclum stage\n{sample_name}', 
                    fontsize=14, y=0.995)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'violin_key_markers_{sample_name}.pdf'), 
                   bbox_inches='tight', dpi=300)
        plt.close()
    
    # 5. Expression across pseudotime
    if 'cyclum_pseudotime' in adata.obs.columns:
        print("\nCreating pseudotime expression plots...")
        
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))
        
        for idx, (phase, genes) in enumerate(phase_markers_present.items()):
            ax = axes[idx]
            
            # Get pseudotime
            pseudotime = adata.obs['cyclum_pseudotime'].values
            
            # Calculate mean expression for this phase's genes
            phase_expr = []
            for gene in genes:
                if gene in adata.var_names:
                    gene_vals = adata[:, gene].X.toarray().flatten() if hasattr(adata.X, 'toarray') else adata[:, gene].X.flatten()
                    phase_expr.append(gene_vals)
            
            if phase_expr:
                mean_phase_expr = np.mean(phase_expr, axis=0)
                
                # Sort by pseudotime for smooth line
                sort_idx = np.argsort(pseudotime)
                sorted_pseudo = pseudotime[sort_idx]
                sorted_expr = mean_phase_expr[sort_idx]
                
                # Smooth with rolling average
                window = len(sorted_pseudo) // 50
                if window > 10:
                    import pandas as pd
                    smoothed = pd.Series(sorted_expr).rolling(window=window, center=True).mean()
                else:
                    smoothed = sorted_expr
                
                ax.scatter(sorted_pseudo, sorted_expr, alpha=0.1, s=1, c='gray')
                ax.plot(sorted_pseudo, smoothed, linewidth=3, 
                       label=f'{phase.upper()} markers (n={len(genes)})',
                       color=phase_colors.get(phase, 'black'))
                
                ax.set_xlabel('Cyclum Pseudotime', fontsize=11)
                ax.set_ylabel('Mean Expression', fontsize=11)
                ax.set_title(f'{phase.upper()} markers across cell cycle', fontsize=12)
                ax.legend()
                ax.grid(True, alpha=0.3)
        
        plt.suptitle(f'Cell cycle marker expression across pseudotime\n{sample_name}', 
                    fontsize=14)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'pseudotime_expression_{sample_name}.pdf'), 
                   bbox_inches='tight', dpi=300)
        plt.close()
    
    # 6. Statistical validation
    print("\n" + "="*70)
    print("STATISTICAL VALIDATION")
    print("="*70)
    
    stats_results = []
    
    for expected_phase, genes in phase_markers_present.items():
        if not genes:
            continue
        
        # Compare expression in expected stage vs other stages
        expected_stage_cells = adata.obs['cyclum_stage'] == expected_phase
        other_stage_cells = adata.obs['cyclum_stage'] != expected_phase
        
        # Get mean expression across all genes for this phase
        expected_expr = []
        other_expr = []
        
        for gene in genes:
            if gene in adata.var_names:
                gene_vals = adata[:, gene].X.toarray().flatten() if hasattr(adata.X, 'toarray') else adata[:, gene].X.flatten()
                expected_expr.extend(gene_vals[expected_stage_cells])
                other_expr.extend(gene_vals[other_stage_cells])
        
        if expected_expr and other_expr:
            # Mann-Whitney U test
            u_stat, p_val = mannwhitneyu(expected_expr, other_expr, alternative='greater')
            
            mean_expected = np.mean(expected_expr)
            mean_other = np.mean(other_expr)
            fold_change = mean_expected / (mean_other + 1e-10)
            
            stats_results.append({
                'phase': expected_phase,
                'n_genes': len(genes),
                'mean_in_phase': mean_expected,
                'mean_other': mean_other,
                'fold_change': fold_change,
                'p_value': p_val,
                'significant': p_val < 0.05
            })
            
            print(f"\n{expected_phase.upper()}:")
            print(f"  Mean expr in {expected_phase}: {mean_expected:.3f}")
            print(f"  Mean expr in other stages: {mean_other:.3f}")
            print(f"  Fold change: {fold_change:.2f}x")
            print(f"  P-value: {p_val:.2e} {'***' if p_val < 0.001 else '**' if p_val < 0.01 else '*' if p_val < 0.05 else 'ns'}")
    
    stats_df = pd.DataFrame(stats_results)
    
    # Save validation results
    stats_df.to_csv(os.path.join(output_dir, f'validation_statistics_{sample_name}.csv'), 
                    index=False)
    
    # 7. Summary report
    print("\n" + "="*70)
    print("VALIDATION SUMMARY")
    print("="*70)
    
    n_correct = sum(1 for score in validation_scores.values() if score['is_correct'])
    n_total = len(validation_scores)
    accuracy = n_correct / n_total if n_total > 0 else 0
    
    n_significant = sum(1 for result in stats_results if result['significant'])
    
    print(f"\nPhase assignment accuracy: {n_correct}/{n_total} ({accuracy*100:.1f}%)")
    print(f"Statistically significant enrichments: {n_significant}/{len(stats_results)}")
    print(f"\nTotal marker genes used: {len(present_markers)}")
    print(f"Missing marker genes: {len(missing_markers)}")
    
    summary = {
        'sample': sample_name,
        'accuracy': accuracy,
        'n_correct_phases': n_correct,
        'n_total_phases': n_total,
        'n_significant': n_significant,
        'n_markers_present': len(present_markers),
        'n_markers_missing': len(missing_markers),
        'validation_scores': validation_scores,
        'statistics': stats_df
    }
    
    # Save summary
    with open(os.path.join(output_dir, f'validation_summary_{sample_name}.txt'), 'w') as f:
        f.write("="*70 + "\n")
        f.write("CYCLUM VALIDATION SUMMARY\n")
        f.write("="*70 + "\n\n")
        f.write(f"Sample: {sample_name}\n")
        f.write(f"Total cells: {adata.n_obs}\n\n")
        f.write(f"Phase assignment accuracy: {n_correct}/{n_total} ({accuracy*100:.1f}%)\n")
        f.write(f"Statistically significant enrichments: {n_significant}/{len(stats_results)}\n")
        f.write(f"Marker genes present: {len(present_markers)}/{len(all_markers)}\n\n")
        
        f.write("Phase-specific results:\n")
        f.write("-"*70 + "\n")
        for phase, score in validation_scores.items():
            f.write(f"\n{phase.upper()}:\n")
            f.write(f"  Genes tested: {score['n_genes']}\n")
            f.write(f"  Peak stage: {score['peak_stage']}\n")
            f.write(f"  Expected: {score['expected_stage']}\n")
            f.write(f"  Match: {'YES' if score['is_correct'] else 'NO'}\n")
            f.write(f"  Expression by stage: {score['mean_expression_by_stage']}\n")
    
    print(f"\nValidation complete! Results saved to {output_dir}")
    
    return summary


def main():
    parser = argparse.ArgumentParser(
        description='Validate Cyclum cell cycle predictions with marker genes',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  python validate_cyclum_markers.py \\
      --input results/method_comparison/all_conditions_all_methods_with_cyclum.h5ad \\
      --output results/cellcycle_validation \\
      --sample all_conditions
        '''
    )
    
    parser.add_argument('--input', '-i', required=True, type=str,
                       help='Path to h5ad file with Cyclum annotations')
    parser.add_argument('--output', '-o', type=str, default='cellcycle_validation',
                       help='Output directory (default: cellcycle_validation)')
    parser.add_argument('--sample', '-s', type=str, default='sample',
                       help='Sample name for output files (default: sample)')
    
    args = parser.parse_args()
    
    # Load data
    print(f"Loading data from {args.input}...")
    adata = sc.read_h5ad(args.input)
    
    print(f"\nLoaded AnnData:")
    print(f"  Cells: {adata.n_obs}")
    print(f"  Genes: {adata.n_vars}")
    print(f"  Cyclum stages: {adata.obs['cyclum_stage'].value_counts().to_dict() if 'cyclum_stage' in adata.obs.columns else 'NOT FOUND'}")
    
    # Run validation
    summary = validate_cyclum_with_markers(adata, args.output, args.sample)
    
    if summary:
        print("\n" + "="*70)
        print("VALIDATION COMPLETE")
        print("="*70)
        print(f"\nOverall accuracy: {summary['accuracy']*100:.1f}%")
        print(f"Check {args.output} for detailed results")


if __name__ == "__main__":
    main()