#!/usr/bin/env python3
"""
Script to plot scRNA-seq QC metrics as boxplots grouped by platform
Author: Jodie Jacobs
Date: 2026-01-29

Example usage:
python scrna_qc_boxplots.py --filtered filtered_dataset.csv --fastqc fastqc_data.csv --output qc_plots

"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
import argparse
from scipy import stats

def load_data(filtered_csv, fastqc_csv):
    """Load and preprocess the CSV files"""
    # Load filtered dataset (comma-separated)
    filtered_df = pd.read_csv(filtered_csv, sep=',')
    
    # Load FastQC data (comma-separated)
    fastqc_df = pd.read_csv(fastqc_csv, sep=',')
    
    # Clean the data - remove rows with missing platform information
    filtered_df = filtered_df.dropna(subset=['platform'])
    fastqc_df = fastqc_df.dropna(subset=['platform'])
    
    # Process FastQC data to get both R1 and R2 scores
    fastqc_r1 = fastqc_df[fastqc_df['ReadNumber'] == 1].copy()
    fastqc_r2 = fastqc_df[fastqc_df['ReadNumber'] == 2].copy()
    
    # Additional cleaning for numeric columns
    numeric_columns = ['n_cells', 'genes_per_cell_median', 'umis_per_cell_median', 
                      'n_genes_total', 'transcriptome_coverage_median', 'doublet_rate']
    
    for col in numeric_columns:
        if col in filtered_df.columns:
            filtered_df[col] = pd.to_numeric(filtered_df[col], errors='coerce')
    
    # Remove rows where critical numeric columns are NaN
    filtered_df = filtered_df.dropna(subset=['n_cells', 'genes_per_cell_median', 'umis_per_cell_median'])
    
    print(f"Loaded {len(filtered_df)} filtered samples, {len(fastqc_r1)} FastQC R1 records, and {len(fastqc_r2)} FastQC R2 records")
    print(f"Platforms found: {filtered_df['platform'].unique()}")
    
    return filtered_df, fastqc_r1, fastqc_r2

def map_sample_names(sample_name):
    """Map sample names to more readable format"""
    # Remove 'kallisto_' prefix
    clean_name = sample_name.replace('kallisto_', '')
    
    # Map specific names
    if 'JW18DOX' in clean_name:
        return clean_name.replace('JW18DOX', 'uninfected')
    elif 'JW18wMel' in clean_name:
        return clean_name.replace('JW18wMel', 'wMel-infected')
    else:
        return clean_name

def perform_statistical_test(data, metric_col):
    """
    Perform Mann-Whitney U test to compare two platforms
    Returns U statistic and p-value
    """
    platform1_data = data[data['platform'] == '10x'][metric_col].values
    platform2_data = data[data['platform'] == 'pipseq'][metric_col].values
    
    # Perform Mann-Whitney U test (non-parametric)
    statistic, p_value = stats.mannwhitneyu(platform1_data, platform2_data, alternative='two-sided')
    
    return statistic, p_value

def format_pvalue(p_value):
    """Format p-value for display"""
    if p_value < 0.001:
        return "p < 0.001"
    elif p_value < 0.01:
        return f"p = {p_value:.3f}"
    elif p_value < 0.05:
        return f"p = {p_value:.3f}"
    else:
        return f"p = {p_value:.3f}"

def get_significance_stars(p_value):
    """Convert p-value to significance stars"""
    if p_value < 0.001:
        return "***"
    elif p_value < 0.01:
        return "**"
    elif p_value < 0.05:
        return "*"
    else:
        return "ns"

def create_boxplot(ax, data, metric_col, ylabel, platform_colors, log_scale=False):
    """Helper function to create a single boxplot with overlayed points and statistics"""
    
    # Perform statistical test
    statistic, p_value = perform_statistical_test(data, metric_col)
    sig_stars = get_significance_stars(p_value)
    
    # Create boxplot
    box_parts = ax.boxplot([data[data['platform'] == platform][metric_col].values 
                            for platform in ['10x', 'pipseq']],
                          positions=[0, 1],
                          widths=0.5,
                          patch_artist=True,
                          showfliers=False,  # Don't show outliers, we'll plot all points
                          medianprops=dict(color='black', linewidth=1.5),
                          boxprops=dict(linewidth=1),
                          whiskerprops=dict(linewidth=1),
                          capprops=dict(linewidth=1))
    
    # Color the boxes
    for patch, platform in zip(box_parts['boxes'], ['10x', 'pipseq']):
        patch.set_facecolor(platform_colors[platform])
        patch.set_alpha(0.6)
    
    # Overlay individual points with slight jitter
    for i, platform in enumerate(['10x', 'pipseq']):
        platform_data = data[data['platform'] == platform][metric_col].values
        x_jitter = np.random.normal(i, 0.04, size=len(platform_data))
        ax.scatter(x_jitter, platform_data, 
                  color=platform_colors[platform], 
                  alpha=0.8, 
                  s=30, 
                  edgecolors='black', 
                  linewidths=0.5,
                  zorder=3)
    
    # Add significance annotation
    y_max = data[metric_col].max()
    y_min = data[metric_col].min()
    y_range = y_max - y_min
    
    if log_scale:
        # For log scale, work in log space
        log_max = np.log10(y_max)
        log_min = np.log10(y_min)
        log_range = log_max - log_min
        y_sig = 10 ** (log_max + 0.1 * log_range)
    else:
        y_sig = y_max + 0.15 * y_range
    
    # Draw significance bracket
    ax.plot([0, 0, 1, 1], [y_sig, y_sig * 1.02, y_sig * 1.02, y_sig], 
            'k-', linewidth=1)
    ax.text(0.5, y_sig * 1.03, sig_stars, ha='center', va='bottom', fontsize=8)
    
    # Set labels and formatting
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['10x', 'PIPseq'], fontsize=8)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.tick_params(axis='both', labelsize=8)
    
    if log_scale:
        ax.set_yscale('log')
    
    # Show all spines (full box)
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(True)
    
    return p_value

def create_boxplot_plots(filtered_df, fastqc_r1, fastqc_r2, output_dir="plots"):
    """Create boxplot plots for all specified metrics"""
    
    # Create output directory
    Path(output_dir).mkdir(exist_ok=True)
    
    # Set up clean plotting style with Arial font
    plt.style.use('default')
    plt.rcParams['figure.facecolor'] = 'white'
    plt.rcParams['axes.facecolor'] = 'white'
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 8
    
    # Define colors for platforms
    platform_colors = {'10x': '#4682B4', 'pipseq': '#FF8C00'}
    
    # Define figure size in inches
    fig_size = (2, 2.25)
    
    # Clean sample names
    filtered_df['sample_name_clean'] = filtered_df['sample_name'].apply(map_sample_names)
    fastqc_r1['sample_name_clean'] = fastqc_r1['SampleID'].apply(map_sample_names)
    fastqc_r2['sample_name_clean'] = fastqc_r2['SampleID'].apply(map_sample_names)
    
    # Store p-values for summary
    p_values = {}
    
    # 1. Raw data quality (FastQC score R1)
    fig1, ax1 = plt.subplots(figsize=fig_size)
    p_values['FastQC R1'] = create_boxplot(ax1, fastqc_r1, 'Mean Quality Score (PF)', 
                                           'Mean Quality Score (R1)', platform_colors)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/fastqc_quality_r1_boxplot.svg', bbox_inches='tight')
    plt.savefig(f'{output_dir}/fastqc_quality_r1_boxplot.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Raw data quality (FastQC score R2)
    fig2, ax2 = plt.subplots(figsize=fig_size)
    p_values['FastQC R2'] = create_boxplot(ax2, fastqc_r2, 'Mean Quality Score (PF)', 
                                           'Mean Quality Score (R2)', platform_colors)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/fastqc_quality_r2_boxplot.svg', bbox_inches='tight')
    plt.savefig(f'{output_dir}/fastqc_quality_r2_boxplot.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. Number of cells (filtered) - with log scale
    fig3, ax3 = plt.subplots(figsize=fig_size)
    p_values['Number of Cells'] = create_boxplot(ax3, filtered_df, 'n_cells', 
                                                  'Number of Cells', platform_colors, log_scale=True)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/n_cells_filtered_boxplot.svg', bbox_inches='tight')
    plt.savefig(f'{output_dir}/n_cells_filtered_boxplot.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 4. Genes per cell (filtered)
    fig4, ax4 = plt.subplots(figsize=fig_size)
    p_values['Genes per Cell'] = create_boxplot(ax4, filtered_df, 'genes_per_cell_median', 
                                                 'Genes per Cell (Median)', platform_colors)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/genes_per_cell_filtered_boxplot.svg', bbox_inches='tight')
    plt.savefig(f'{output_dir}/genes_per_cell_filtered_boxplot.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 5. UMIs per cell (filtered)
    fig5, ax5 = plt.subplots(figsize=fig_size)
    p_values['UMIs per Cell'] = create_boxplot(ax5, filtered_df, 'umis_per_cell_median', 
                                                'UMIs per Cell (Median)', platform_colors)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/umis_per_cell_filtered_boxplot.svg', bbox_inches='tight')
    plt.savefig(f'{output_dir}/umis_per_cell_filtered_boxplot.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 6. Number of genes total (filtered)
    fig6, ax6 = plt.subplots(figsize=fig_size)
    p_values['Total Genes'] = create_boxplot(ax6, filtered_df, 'n_genes_total', 
                                              'Total Number of Genes', platform_colors)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/n_genes_total_filtered_boxplot.svg', bbox_inches='tight')
    plt.savefig(f'{output_dir}/n_genes_total_filtered_boxplot.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 7. Transcriptome coverage (median)
    fig7, ax7 = plt.subplots(figsize=fig_size)
    p_values['Transcriptome Coverage'] = create_boxplot(ax7, filtered_df, 'transcriptome_coverage_median', 
                                                         'Transcriptome Coverage (Median)', platform_colors)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/transcriptome_coverage_median_boxplot.svg', bbox_inches='tight')
    plt.savefig(f'{output_dir}/transcriptome_coverage_median_boxplot.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 8. Doublet rate
    fig8, ax8 = plt.subplots(figsize=fig_size)
    p_values['Doublet Rate'] = create_boxplot(ax8, filtered_df, 'doublet_rate', 
                                               'Doublet Rate', platform_colors)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/doublet_rate_boxplot.svg', bbox_inches='tight')
    plt.savefig(f'{output_dir}/doublet_rate_boxplot.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create combined figure with all 8 plots
    create_combined_figure(fastqc_r1, fastqc_r2, filtered_df, platform_colors, output_dir)
    
    return p_values

def create_combined_figure(fastqc_r1, fastqc_r2, filtered_df, platform_colors, output_dir):
    """Create a combined figure with all 8 boxplots"""
    
    # Set up the combined figure (6.5 x 9 inches, 4 rows x 2 columns)
    fig, axes = plt.subplots(4, 2, figsize=(6.5, 9))
    fig.subplots_adjust(hspace=0.5, wspace=0.4)
    
    # Flatten axes for easier iteration
    axes_flat = axes.flatten()
    
    # Define all metrics to plot
    metrics = [
        (fastqc_r1, 'Mean Quality Score (PF)', 'Mean Quality Score (R1)', False),
        (fastqc_r2, 'Mean Quality Score (PF)', 'Mean Quality Score (R2)', False),
        (filtered_df, 'n_cells', 'Number of Cells', True),
        (filtered_df, 'genes_per_cell_median', 'Genes per Cell (Median)', False),
        (filtered_df, 'umis_per_cell_median', 'UMIs per Cell (Median)', False),
        (filtered_df, 'n_genes_total', 'Total Number of Genes', False),
        (filtered_df, 'transcriptome_coverage_median', 'Transcriptome Coverage (Median)', False),
        (filtered_df, 'doublet_rate', 'Doublet Rate', False),
    ]
    
    # Create each boxplot
    for ax, (data, metric_col, ylabel, log_scale) in zip(axes_flat, metrics):
        create_boxplot(ax, data, metric_col, ylabel, platform_colors, log_scale)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/combined_qc_metrics_boxplot.svg', bbox_inches='tight')
    plt.savefig(f'{output_dir}/combined_qc_metrics_boxplot.pdf', dpi=300, bbox_inches='tight')
    plt.close()

def create_summary_stats(filtered_df, fastqc_r1, fastqc_r2, p_values):
    """Generate summary statistics for all metrics, grouped by platform"""
    
    print("="*60)
    print("SUMMARY STATISTICS FOR scRNA-seq QC METRICS BY PLATFORM")
    print("="*60)
    
    for platform in ['10x', 'pipseq']:
        print(f"\n{'='*60}")
        print(f"PLATFORM: {platform.upper()}")
        print(f"{'='*60}")
        
        # Filter data for this platform
        fastqc_r1_plat = fastqc_r1[fastqc_r1['platform'] == platform]
        fastqc_r2_plat = fastqc_r2[fastqc_r2['platform'] == platform]
        filtered_plat = filtered_df[filtered_df['platform'] == platform]
        
        # FastQC Quality (R1)
        print(f"\n1. FastQC Quality Scores (R1):")
        print(f"   Mean: {fastqc_r1_plat['Mean Quality Score (PF)'].mean():.2f}")
        print(f"   Median: {fastqc_r1_plat['Mean Quality Score (PF)'].median():.2f}")
        print(f"   Range: {fastqc_r1_plat['Mean Quality Score (PF)'].min():.2f} - {fastqc_r1_plat['Mean Quality Score (PF)'].max():.2f}")
        
        # FastQC Quality (R2)
        print(f"\n2. FastQC Quality Scores (R2):")
        print(f"   Mean: {fastqc_r2_plat['Mean Quality Score (PF)'].mean():.2f}")
        print(f"   Median: {fastqc_r2_plat['Mean Quality Score (PF)'].median():.2f}")
        print(f"   Range: {fastqc_r2_plat['Mean Quality Score (PF)'].min():.2f} - {fastqc_r2_plat['Mean Quality Score (PF)'].max():.2f}")
        
        # Cell numbers
        print(f"\n3. Number of Cells (Filtered):")
        print(f"   Mean: {filtered_plat['n_cells'].mean():.0f}")
        print(f"   Median: {filtered_plat['n_cells'].median():.0f}")
        print(f"   Range: {filtered_plat['n_cells'].min()} - {filtered_plat['n_cells'].max()}")
        
        # Genes per cell
        print(f"\n4. Genes per Cell (Median):")
        print(f"   Mean: {filtered_plat['genes_per_cell_median'].mean():.0f}")
        print(f"   Median: {filtered_plat['genes_per_cell_median'].median():.0f}")
        print(f"   Range: {filtered_plat['genes_per_cell_median'].min()} - {filtered_plat['genes_per_cell_median'].max()}")
        
        # UMIs per cell
        print(f"\n5. UMIs per Cell (Median):")
        print(f"   Mean: {filtered_plat['umis_per_cell_median'].mean():.0f}")
        print(f"   Median: {filtered_plat['umis_per_cell_median'].median():.0f}")
        print(f"   Range: {filtered_plat['umis_per_cell_median'].min()} - {filtered_plat['umis_per_cell_median'].max()}")
        
        # Total genes
        print(f"\n6. Total Genes (Filtered):")
        print(f"   Mean: {filtered_plat['n_genes_total'].mean():.0f}")
        print(f"   Median: {filtered_plat['n_genes_total'].median():.0f}")
        print(f"   Range: {filtered_plat['n_genes_total'].min()} - {filtered_plat['n_genes_total'].max()}")
        
        # Transcriptome coverage
        print(f"\n7. Transcriptome Coverage (Median):")
        print(f"   Mean: {filtered_plat['transcriptome_coverage_median'].mean():.4f}")
        print(f"   Median: {filtered_plat['transcriptome_coverage_median'].median():.4f}")
        print(f"   Range: {filtered_plat['transcriptome_coverage_median'].min():.4f} - {filtered_plat['transcriptome_coverage_median'].max():.4f}")
        
        # Doublet rate
        print(f"\n8. Doublet Rate:")
        print(f"   Mean: {filtered_plat['doublet_rate'].mean():.4f}")
        print(f"   Median: {filtered_plat['doublet_rate'].median():.4f}")
        print(f"   Range: {filtered_plat['doublet_rate'].min():.4f} - {filtered_plat['doublet_rate'].max():.4f}")
        
        print(f"\n   N samples: {len(filtered_plat)}")
    
    # Print statistical test results
    print(f"\n{'='*60}")
    print("STATISTICAL COMPARISONS (Mann-Whitney U Test)")
    print(f"{'='*60}")
    print("\nMetric                          p-value      Significance")
    print("-" * 60)
    
    for metric_name, p_val in p_values.items():
        sig = get_significance_stars(p_val)
        print(f"{metric_name:<30} {p_val:>10.4f}      {sig}")
    
    print("\nSignificance levels: *** p<0.001, ** p<0.01, * p<0.05, ns = not significant")

def main():
    """Main function to run the analysis"""
    parser = argparse.ArgumentParser(description='Generate boxplots for scRNA-seq QC metrics grouped by platform')
    parser.add_argument('--filtered', '-f', required=True, 
                       help='Path to filtered dataset CSV file')
    parser.add_argument('--fastqc', '-q', required=True,
                       help='Path to FastQC CSV file')
    parser.add_argument('--output', '-o', default='plots',
                       help='Output directory for plots (default: plots)')
    
    args = parser.parse_args()
    
    # Load data
    print("Loading data...")
    filtered_df, fastqc_r1, fastqc_r2 = load_data(args.filtered, args.fastqc)
    
    # Create plots and get p-values
    print(f"\nGenerating boxplot plots...")
    p_values = create_boxplot_plots(filtered_df, fastqc_r1, fastqc_r2, args.output)
    
    # Generate summary statistics with p-values
    create_summary_stats(filtered_df, fastqc_r1, fastqc_r2, p_values)
    
    print(f"\nPlots saved to: {args.output}/")
    print("Analysis complete!")

if __name__ == "__main__":
    main()

