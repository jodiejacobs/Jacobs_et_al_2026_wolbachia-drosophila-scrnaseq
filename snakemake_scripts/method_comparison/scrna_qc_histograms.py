#!/usr/bin/env python3
"""
Script to plot scRNA-seq QC metrics as histograms
Author: Bioinformatics analysis script
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
import argparse

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

def create_histogram_plots(filtered_df, fastqc_r1, fastqc_r2, output_dir="plots"):
    """Create histogram plots for all specified metrics"""
    
    # Create output directory
    Path(output_dir).mkdir(exist_ok=True)
    
    # Set up clean plotting style
    plt.style.use('default')
    plt.rcParams['figure.facecolor'] = 'white'
    plt.rcParams['axes.facecolor'] = 'white'
    plt.rcParams['font.size'] = 5  # Set minimum font size
    
    # Define colors for platforms
    platform_colors = {'10x': '#4682B4', 'pipseq': '#FF8C00'}  # Medium blue and orange
    
    # Define figure size in inches
    fig_size = (3.25, 2.25)
    
    # Clean sample names
    filtered_df['sample_name_clean'] = filtered_df['sample_name'].apply(map_sample_names)
    fastqc_r1['sample_name_clean'] = fastqc_r1['SampleID'].apply(map_sample_names)
    fastqc_r2['sample_name_clean'] = fastqc_r2['SampleID'].apply(map_sample_names)
    
    # Store all figure objects for the combined plot
    figures = []
    
    # 1. Raw data quality (FastQC score R1) - Bar plot
    fig1, ax1 = plt.subplots(figsize=fig_size)
    fastqc_r1_sorted = fastqc_r1.sort_values(['platform', 'sample_name_clean'])
    x_pos = range(len(fastqc_r1_sorted))
    
    bars = ax1.bar(x_pos, fastqc_r1_sorted['Mean Quality Score (PF)'], 
                   color=[platform_colors[platform] for platform in fastqc_r1_sorted['platform']],
                   alpha=0.8, edgecolor='none')
    
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels(fastqc_r1_sorted['sample_name_clean'], rotation=45, ha='right')
    ax1.set_ylabel('Mean Quality Score')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/fastqc_quality_r1.svg', bbox_inches='tight')
    plt.savefig(f'{output_dir}/fastqc_quality_r1.pdf', dpi=300, bbox_inches='tight')
    figures.append((fig1, ax1, 'FastQC R1'))
    
    # 2. Raw data quality (FastQC score R2) - Bar plot  
    fig2, ax2 = plt.subplots(figsize=fig_size)
    fastqc_r2_sorted = fastqc_r2.sort_values(['platform', 'sample_name_clean'])
    x_pos = range(len(fastqc_r2_sorted))
    
    bars = ax2.bar(x_pos, fastqc_r2_sorted['Mean Quality Score (PF)'], 
                   color=[platform_colors[platform] for platform in fastqc_r2_sorted['platform']],
                   alpha=0.8, edgecolor='none')
    
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(fastqc_r2_sorted['sample_name_clean'], rotation=45, ha='right')
    ax2.set_ylabel('Mean Quality Score')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/fastqc_quality_r2.svg', bbox_inches='tight')
    plt.savefig(f'{output_dir}/fastqc_quality_r2.pdf', dpi=300, bbox_inches='tight')
    figures.append((fig2, ax2, 'FastQC R2'))
    
    # 3. Number of cells (filtered) - Bar plot with log y-axis
    fig3, ax3 = plt.subplots(figsize=fig_size)
    filtered_sorted = filtered_df.sort_values(['platform', 'sample_name_clean'])
    x_pos = range(len(filtered_sorted))
    
    bars = ax3.bar(x_pos, filtered_sorted['n_cells'], 
                   color=[platform_colors[platform] for platform in filtered_sorted['platform']],
                   alpha=0.8, edgecolor='none')
    
    ax3.set_xticks(x_pos)
    ax3.set_xticklabels(filtered_sorted['sample_name_clean'], rotation=45, ha='right')
    ax3.set_yscale('log')
    ax3.set_ylabel('Number of Cells (log scale)')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/n_cells_filtered.svg', bbox_inches='tight')
    plt.savefig(f'{output_dir}/n_cells_filtered.pdf', dpi=300, bbox_inches='tight')
    figures.append((fig3, ax3, 'Number of Cells'))
    
    # 4. Genes per cell (filtered) - Bar plot using median
    fig4, ax4 = plt.subplots(figsize=fig_size)
    bars = ax4.bar(x_pos, filtered_sorted['genes_per_cell_median'], 
                   color=[platform_colors[platform] for platform in filtered_sorted['platform']],
                   alpha=0.8, edgecolor='none')
    
    ax4.set_xticks(x_pos)
    ax4.set_xticklabels(filtered_sorted['sample_name_clean'], rotation=45, ha='right')
    ax4.set_ylabel('Genes per Cell (Median)')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/genes_per_cell_filtered.svg', bbox_inches='tight')
    plt.savefig(f'{output_dir}/genes_per_cell_filtered.pdf', dpi=300, bbox_inches='tight')
    figures.append((fig4, ax4, 'Genes per Cell'))
    
    # 5. UMIs per cell (filtered) - Bar plot using median
    fig5, ax5 = plt.subplots(figsize=fig_size)
    bars = ax5.bar(x_pos, filtered_sorted['umis_per_cell_median'], 
                   color=[platform_colors[platform] for platform in filtered_sorted['platform']],
                   alpha=0.8, edgecolor='none')
    
    ax5.set_xticks(x_pos)
    ax5.set_xticklabels(filtered_sorted['sample_name_clean'], rotation=45, ha='right')
    ax5.set_ylabel('UMIs per Cell (Median)')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/umis_per_cell_filtered.svg', bbox_inches='tight')
    plt.savefig(f'{output_dir}/umis_per_cell_filtered.pdf', dpi=300, bbox_inches='tight')
    figures.append((fig5, ax5, 'UMIs per Cell'))
    
    # 6. Number of genes total (filtered) - Bar plot
    fig6, ax6 = plt.subplots(figsize=fig_size)
    bars = ax6.bar(x_pos, filtered_sorted['n_genes_total'], 
                   color=[platform_colors[platform] for platform in filtered_sorted['platform']],
                   alpha=0.8, edgecolor='none')
    
    ax6.set_xticks(x_pos)
    ax6.set_xticklabels(filtered_sorted['sample_name_clean'], rotation=45, ha='right')
    ax6.set_ylabel('Total Number of Genes')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/n_genes_total_filtered.svg', bbox_inches='tight')
    plt.savefig(f'{output_dir}/n_genes_total_filtered.pdf', dpi=300, bbox_inches='tight')
    figures.append((fig6, ax6, 'Total Genes'))
    
    # 7. Transcriptome coverage (median) - Bar plot
    fig7, ax7 = plt.subplots(figsize=fig_size)
    bars = ax7.bar(x_pos, filtered_sorted['transcriptome_coverage_median'], 
                   color=[platform_colors[platform] for platform in filtered_sorted['platform']],
                   alpha=0.8, edgecolor='none')
    
    ax7.set_xticks(x_pos)
    ax7.set_xticklabels(filtered_sorted['sample_name_clean'], rotation=45, ha='right')
    ax7.set_ylabel('Transcriptome Coverage (Median)')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/transcriptome_coverage_median.svg', bbox_inches='tight')
    plt.savefig(f'{output_dir}/transcriptome_coverage_median.pdf', dpi=300, bbox_inches='tight')
    figures.append((fig7, ax7, 'Transcriptome Coverage'))
    
    # 8. Doublet rate - Bar plot
    fig8, ax8 = plt.subplots(figsize=fig_size)
    bars = ax8.bar(x_pos, filtered_sorted['doublet_rate'], 
                   color=[platform_colors[platform] for platform in filtered_sorted['platform']],
                   alpha=0.8, edgecolor='none')
    
    ax8.set_xticks(x_pos)
    ax8.set_xticklabels(filtered_sorted['sample_name_clean'], rotation=45, ha='right')
    ax8.set_ylabel('Doublet Rate')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/doublet_rate.svg', bbox_inches='tight')
    plt.savefig(f'{output_dir}/doublet_rate.pdf', dpi=300, bbox_inches='tight')
    figures.append((fig8, ax8, 'Doublet Rate'))
    
    # Create combined figure with all 8 plots
    create_combined_figure(fastqc_r1, fastqc_r2, filtered_df, platform_colors, output_dir)

def create_combined_figure(fastqc_r1, fastqc_r2, filtered_df, platform_colors, output_dir):
    """Create a combined figure with all 8 plots"""
    
    # Set up the combined figure (6.5 x 9 inches, 4 rows x 2 columns)
    fig, axes = plt.subplots(4, 2, figsize=(6.5, 9))
    fig.subplots_adjust(hspace=0.4, wspace=0.3)
    
    # Sort data
    fastqc_r1_sorted = fastqc_r1.sort_values(['platform', 'sample_name_clean'])
    fastqc_r2_sorted = fastqc_r2.sort_values(['platform', 'sample_name_clean'])
    filtered_sorted = filtered_df.sort_values(['platform', 'sample_name_clean'])
    
    # Plot 1: FastQC R1
    ax = axes[0, 0]
    x_pos = range(len(fastqc_r1_sorted))
    ax.bar(x_pos, fastqc_r1_sorted['Mean Quality Score (PF)'], 
           color=[platform_colors[platform] for platform in fastqc_r1_sorted['platform']],
           alpha=0.8, edgecolor='none')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(fastqc_r1_sorted['sample_name_clean'], rotation=45, ha='right', fontsize=5)
    ax.set_ylabel('Mean Quality Score', fontsize=5)
    ax.tick_params(axis='y', labelsize=5)
    
    # Plot 2: FastQC R2
    ax = axes[0, 1]
    x_pos = range(len(fastqc_r2_sorted))
    ax.bar(x_pos, fastqc_r2_sorted['Mean Quality Score (PF)'], 
           color=[platform_colors[platform] for platform in fastqc_r2_sorted['platform']],
           alpha=0.8, edgecolor='none')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(fastqc_r2_sorted['sample_name_clean'], rotation=45, ha='right', fontsize=5)
    ax.set_ylabel('Mean Quality Score', fontsize=5)
    ax.tick_params(axis='y', labelsize=5)
    
    # Plot 3: Number of cells
    ax = axes[1, 0]
    x_pos = range(len(filtered_sorted))
    ax.bar(x_pos, filtered_sorted['n_cells'], 
           color=[platform_colors[platform] for platform in filtered_sorted['platform']],
           alpha=0.8, edgecolor='none')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(filtered_sorted['sample_name_clean'], rotation=45, ha='right', fontsize=5)
    ax.set_yscale('log')
    ax.set_ylabel('Number of Cells (log scale)', fontsize=5)
    ax.tick_params(axis='y', labelsize=5)
    
    # Plot 4: Genes per cell
    ax = axes[1, 1]
    ax.bar(x_pos, filtered_sorted['genes_per_cell_median'], 
           color=[platform_colors[platform] for platform in filtered_sorted['platform']],
           alpha=0.8, edgecolor='none')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(filtered_sorted['sample_name_clean'], rotation=45, ha='right', fontsize=5)
    ax.set_ylabel('Genes per Cell (Median)', fontsize=5)
    ax.tick_params(axis='y', labelsize=5)
    
    # Plot 5: UMIs per cell
    ax = axes[2, 0]
    ax.bar(x_pos, filtered_sorted['umis_per_cell_median'], 
           color=[platform_colors[platform] for platform in filtered_sorted['platform']],
           alpha=0.8, edgecolor='none')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(filtered_sorted['sample_name_clean'], rotation=45, ha='right', fontsize=5)
    ax.set_ylabel('UMIs per Cell (Median)', fontsize=5)
    ax.tick_params(axis='y', labelsize=5)
    
    # Plot 6: Total genes
    ax = axes[2, 1]
    ax.bar(x_pos, filtered_sorted['n_genes_total'], 
           color=[platform_colors[platform] for platform in filtered_sorted['platform']],
           alpha=0.8, edgecolor='none')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(filtered_sorted['sample_name_clean'], rotation=45, ha='right', fontsize=5)
    ax.set_ylabel('Total Number of Genes', fontsize=5)
    ax.tick_params(axis='y', labelsize=5)
    
    # Plot 7: Transcriptome coverage
    ax = axes[3, 0]
    ax.bar(x_pos, filtered_sorted['transcriptome_coverage_median'], 
           color=[platform_colors[platform] for platform in filtered_sorted['platform']],
           alpha=0.8, edgecolor='none')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(filtered_sorted['sample_name_clean'], rotation=45, ha='right', fontsize=5)
    ax.set_ylabel('Transcriptome Coverage (Median)', fontsize=5)
    ax.tick_params(axis='y', labelsize=5)
    
    # Plot 8: Doublet rate
    ax = axes[3, 1]
    ax.bar(x_pos, filtered_sorted['doublet_rate'], 
           color=[platform_colors[platform] for platform in filtered_sorted['platform']],
           alpha=0.8, edgecolor='none')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(filtered_sorted['sample_name_clean'], rotation=45, ha='right', fontsize=5)
    ax.set_ylabel('Doublet Rate', fontsize=5)
    ax.tick_params(axis='y', labelsize=5)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/combined_qc_metrics.svg', bbox_inches='tight')
    plt.savefig(f'{output_dir}/combined_qc_metrics.pdf', dpi=300, bbox_inches='tight')

def create_summary_stats(filtered_df, fastqc_r1, fastqc_r2):
    """Generate summary statistics for all metrics"""
    
    print("="*60)
    print("SUMMARY STATISTICS FOR scRNA-seq QC METRICS")
    print("="*60)
    
    # FastQC Quality (R1)
    print(f"\n1. FastQC Quality Scores (R1):")
    print(f"   Mean: {fastqc_r1['Mean Quality Score (PF)'].mean():.2f}")
    print(f"   Median: {fastqc_r1['Mean Quality Score (PF)'].median():.2f}")
    print(f"   Range: {fastqc_r1['Mean Quality Score (PF)'].min():.2f} - {fastqc_r1['Mean Quality Score (PF)'].max():.2f}")
    
    # FastQC Quality (R2)
    print(f"\n2. FastQC Quality Scores (R2):")
    print(f"   Mean: {fastqc_r2['Mean Quality Score (PF)'].mean():.2f}")
    print(f"   Median: {fastqc_r2['Mean Quality Score (PF)'].median():.2f}")
    print(f"   Range: {fastqc_r2['Mean Quality Score (PF)'].min():.2f} - {fastqc_r2['Mean Quality Score (PF)'].max():.2f}")
    
    # Cell numbers
    print(f"\n3. Number of Cells (Filtered):")
    print(f"   Mean: {filtered_df['n_cells'].mean():.0f}")
    print(f"   Median: {filtered_df['n_cells'].median():.0f}")
    print(f"   Range: {filtered_df['n_cells'].min()} - {filtered_df['n_cells'].max()}")
    
    # Genes per cell
    print(f"\n4. Genes per Cell (Median):")
    print(f"   Mean: {filtered_df['genes_per_cell_median'].mean():.0f}")
    print(f"   Median: {filtered_df['genes_per_cell_median'].median():.0f}")
    print(f"   Range: {filtered_df['genes_per_cell_median'].min()} - {filtered_df['genes_per_cell_median'].max()}")
    
    # UMIs per cell
    print(f"\n5. UMIs per Cell (Median):")
    print(f"   Mean: {filtered_df['umis_per_cell_median'].mean():.0f}")
    print(f"   Median: {filtered_df['umis_per_cell_median'].median():.0f}")
    print(f"   Range: {filtered_df['umis_per_cell_median'].min()} - {filtered_df['umis_per_cell_median'].max()}")
    
    # Total genes
    print(f"\n6. Total Genes (Filtered):")
    print(f"   Mean: {filtered_df['n_genes_total'].mean():.0f}")
    print(f"   Median: {filtered_df['n_genes_total'].median():.0f}")
    print(f"   Range: {filtered_df['n_genes_total'].min()} - {filtered_df['n_genes_total'].max()}")
    
    # Transcriptome coverage
    print(f"\n7. Transcriptome Coverage (Median):")
    print(f"   Mean: {filtered_df['transcriptome_coverage_median'].mean():.4f}")
    print(f"   Median: {filtered_df['transcriptome_coverage_median'].median():.4f}")
    print(f"   Range: {filtered_df['transcriptome_coverage_median'].min():.4f} - {filtered_df['transcriptome_coverage_median'].max():.4f}")
    
    # Doublet rate
    print(f"\n8. Doublet Rate:")
    print(f"   Mean: {filtered_df['doublet_rate'].mean():.4f}")
    print(f"   Median: {filtered_df['doublet_rate'].median():.4f}")
    print(f"   Range: {filtered_df['doublet_rate'].min():.4f} - {filtered_df['doublet_rate'].max():.4f}")
    
    # Platform breakdown
    print(f"\n9. Platform Distribution:")
    platform_counts = filtered_df['platform'].value_counts()
    for platform, count in platform_counts.items():
        print(f"   {platform}: {count} samples")

def main():
    """Main function to run the analysis"""
    parser = argparse.ArgumentParser(description='Generate histograms for scRNA-seq QC metrics')
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
    
    # Generate summary statistics
    create_summary_stats(filtered_df, fastqc_r1, fastqc_r2)
    
    # Create plots
    print(f"\nGenerating histogram plots...")
    create_histogram_plots(filtered_df, fastqc_r1, fastqc_r2, args.output)
    
    print(f"\nPlots saved to: {args.output}/")
    print("Analysis complete!")

if __name__ == "__main__":
    main()

# Example usage:
# python scrna_qc_histograms.py --filtered filtered_dataset.csv --fastqc fastqc_data.csv --output qc_plots