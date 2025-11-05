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
    
    # Process FastQC data to get R2 scores (Read 2)
    fastqc_r2 = fastqc_df[fastqc_df['ReadNumber'] == 2].copy()
    
    # Additional cleaning for numeric columns
    numeric_columns = ['n_cells', 'genes_per_cell_median', 'umis_per_cell_median', 
                      'n_genes_total', 'transcriptome_coverage_median', 'doublet_rate']
    
    for col in numeric_columns:
        if col in filtered_df.columns:
            filtered_df[col] = pd.to_numeric(filtered_df[col], errors='coerce')
    
    # Remove rows where critical numeric columns are NaN
    filtered_df = filtered_df.dropna(subset=['n_cells', 'genes_per_cell_median', 'umis_per_cell_median'])
    
    print(f"Loaded {len(filtered_df)} filtered samples and {len(fastqc_r2)} FastQC R2 records")
    print(f"Platforms found: {filtered_df['platform'].unique()}")
    
    return filtered_df, fastqc_r2

def create_histogram_plots(filtered_df, fastqc_r2, output_dir="plots"):
    """Create histogram plots for all specified metrics"""
    
    # Create output directory
    Path(output_dir).mkdir(exist_ok=True)
    
    # Set up clean plotting style
    plt.style.use('default')
    plt.rcParams['figure.facecolor'] = 'white'
    plt.rcParams['axes.facecolor'] = 'white'
    
    # Define colors for platforms
    platform_colors = {'10x': '#4682B4', 'pipseq': '#FF8C00'}  # Medium blue and orange
    
    # Define figure size
    fig_size = (10, 6)
    
    # Clean sample names by removing 'kallisto_' prefix
    filtered_df['sample_name_clean'] = filtered_df['sample_name'].str.replace('kallisto_', '', regex=False)
    fastqc_r2['sample_name_clean'] = fastqc_r2['SampleID'].str.replace('kallisto_', '', regex=False)
    
    # 1. Raw data quality (FastQC score R2) - Bar plot
    plt.figure(figsize=(12, 6))
    # Get R2 data and create sample names for x-axis
    fastqc_r2_sorted = fastqc_r2.sort_values(['platform', 'sample_name_clean'])
    
    # Create x-axis positions
    x_pos = range(len(fastqc_r2_sorted))
    
    # Create bars with platform-specific colors
    bars = plt.bar(x_pos, fastqc_r2_sorted['Mean Quality Score (PF)'], 
                   color=[platform_colors[platform] for platform in fastqc_r2_sorted['platform']],
                   alpha=0.8, edgecolor='black', linewidth=0.5)
    
    # Set x-axis labels (sample names)
    plt.xticks(x_pos, fastqc_r2_sorted['sample_name_clean'], rotation=45, ha='right')
    plt.title('FastQC Mean Quality Scores (R2 Reads)', fontsize=14, fontweight='bold')
    plt.xlabel('Sample ID (R2 Reads)')
    plt.ylabel('Mean Quality Score')
    
    # Create custom legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=platform_colors[platform], label=platform) 
                      for platform in platform_colors.keys()]
    plt.legend(handles=legend_elements)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/fastqc_quality_r2.pdf', dpi=300, bbox_inches='tight')
    plt.show()
    
    # 2. Number of cells (filtered) - Bar plot with log y-axis
    plt.figure(figsize=(12, 6))
    filtered_sorted = filtered_df.sort_values(['platform', 'sample_name_clean'])
    x_pos = range(len(filtered_sorted))
    
    bars = plt.bar(x_pos, filtered_sorted['n_cells'], 
                   color=[platform_colors[platform] for platform in filtered_sorted['platform']],
                   alpha=0.8, edgecolor='black', linewidth=0.5)
    
    plt.xticks(x_pos, filtered_sorted['sample_name_clean'], rotation=45, ha='right')
    plt.yscale('log')  # Set log scale for y-axis
    plt.title('Number of Cells (Filtered)', fontsize=14, fontweight='bold')
    plt.xlabel('Sample Name')
    plt.ylabel('Number of Cells (log scale)')
    
    # Create custom legend
    legend_elements = [Patch(facecolor=platform_colors[platform], label=platform) 
                      for platform in platform_colors.keys()]
    plt.legend(handles=legend_elements)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/n_cells_filtered.pdf', dpi=300, bbox_inches='tight')
    plt.show()
    
    # 3. Genes per cell (filtered) - Bar plot using median
    plt.figure(figsize=(12, 6))
    bars = plt.bar(x_pos, filtered_sorted['genes_per_cell_median'], 
                   color=[platform_colors[platform] for platform in filtered_sorted['platform']],
                   alpha=0.8, edgecolor='black', linewidth=0.5)
    
    plt.xticks(x_pos, filtered_sorted['sample_name_clean'], rotation=45, ha='right')
    plt.title('Genes per Cell (Median, Filtered)', fontsize=14, fontweight='bold')
    plt.xlabel('Sample Name')
    plt.ylabel('Genes per Cell (Median)')
    
    plt.legend(handles=legend_elements)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/genes_per_cell_filtered.pdf', dpi=300, bbox_inches='tight')
    plt.show()
    
    # 4. UMIs per cell (filtered) - Bar plot using median
    plt.figure(figsize=(12, 6))
    bars = plt.bar(x_pos, filtered_sorted['umis_per_cell_median'], 
                   color=[platform_colors[platform] for platform in filtered_sorted['platform']],
                   alpha=0.8, edgecolor='black', linewidth=0.5)
    
    plt.xticks(x_pos, filtered_sorted['sample_name_clean'], rotation=45, ha='right')
    plt.title('UMIs per Cell (Median, Filtered)', fontsize=14, fontweight='bold')
    plt.xlabel('Sample Name')
    plt.ylabel('UMIs per Cell (Median)')
    
    plt.legend(handles=legend_elements)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/umis_per_cell_filtered.pdf', dpi=300, bbox_inches='tight')
    plt.show()
    
    # 5. Number of genes total (filtered) - Bar plot
    plt.figure(figsize=(12, 6))
    bars = plt.bar(x_pos, filtered_sorted['n_genes_total'], 
                   color=[platform_colors[platform] for platform in filtered_sorted['platform']],
                   alpha=0.8, edgecolor='black', linewidth=0.5)
    
    plt.xticks(x_pos, filtered_sorted['sample_name_clean'], rotation=45, ha='right')
    plt.title('Total Number of Genes (Filtered)', fontsize=14, fontweight='bold')
    plt.xlabel('Sample Name')
    plt.ylabel('Total Number of Genes')
    
    plt.legend(handles=legend_elements)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/n_genes_total_filtered.pdf', dpi=300, bbox_inches='tight')
    plt.show()
    
    # 6. Transcriptome coverage (median) - Bar plot
    plt.figure(figsize=(12, 6))
    bars = plt.bar(x_pos, filtered_sorted['transcriptome_coverage_median'], 
                   color=[platform_colors[platform] for platform in filtered_sorted['platform']],
                   alpha=0.8, edgecolor='black', linewidth=0.5)
    
    plt.xticks(x_pos, filtered_sorted['sample_name_clean'], rotation=45, ha='right')
    plt.title('Transcriptome Coverage (Median)', fontsize=14, fontweight='bold')
    plt.xlabel('Sample Name')
    plt.ylabel('Transcriptome Coverage (Median)')
    
    plt.legend(handles=legend_elements)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/transcriptome_coverage_median.pdf', dpi=300, bbox_inches='tight')
    plt.show()
    
    # 7. Doublet rate - Bar plot
    plt.figure(figsize=(12, 6))
    bars = plt.bar(x_pos, filtered_sorted['doublet_rate'], 
                   color=[platform_colors[platform] for platform in filtered_sorted['platform']],
                   alpha=0.8, edgecolor='black', linewidth=0.5)
    
    plt.xticks(x_pos, filtered_sorted['sample_name_clean'], rotation=45, ha='right')
    plt.title('Doublet Rate', fontsize=14, fontweight='bold')
    plt.xlabel('Sample Name')
    plt.ylabel('Doublet Rate')
    
    plt.legend(handles=legend_elements)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/doublet_rate.pdf', dpi=300, bbox_inches='tight')
    plt.show()
    
    # 8. UMIs per Barcode (Total UMIs / Number of Cells) - Bar plot
    plt.figure(figsize=(12, 6))
    
    # Calculate total UMIs per sample (UMIs per cell * number of cells)
    total_umis = filtered_sorted['umis_per_cell_median'] * filtered_sorted['n_cells']
    umis_per_barcode = total_umis / filtered_sorted['n_cells']  # This is just umis_per_cell_median, but keeping for clarity
    
    # Actually, let's calculate it as total UMIs divided by total barcodes for the entire sample
    # This gives us the overall UMI efficiency per barcode
    bars = plt.bar(x_pos, filtered_sorted['umis_per_cell_median'], 
                   color=[platform_colors[platform] for platform in filtered_sorted['platform']],
                   alpha=0.8, edgecolor='black', linewidth=0.5)
    
    plt.xticks(x_pos, filtered_sorted['sample_name_clean'], rotation=45, ha='right')
    plt.title('UMIs per Barcode (UMIs/Barcodes)', fontsize=14, fontweight='bold')
    plt.xlabel('Sample Name')
    plt.ylabel('UMIs per Barcode')
    
    plt.legend(handles=legend_elements)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/umis_per_barcode.pdf', dpi=300, bbox_inches='tight')
    plt.show()

def create_summary_stats(filtered_df, fastqc_r2):
    """Generate summary statistics for all metrics"""
    
    print("="*60)
    print("SUMMARY STATISTICS FOR scRNA-seq QC METRICS")
    print("="*60)
    
    # FastQC Quality (R2)
    print(f"\n1. FastQC Quality Scores (R2):")
    print(f"   Mean: {fastqc_r2['Mean Quality Score (PF)'].mean():.2f}")
    print(f"   Median: {fastqc_r2['Mean Quality Score (PF)'].median():.2f}")
    print(f"   Range: {fastqc_r2['Mean Quality Score (PF)'].min():.2f} - {fastqc_r2['Mean Quality Score (PF)'].max():.2f}")
    
    # Cell numbers
    print(f"\n2. Number of Cells (Filtered):")
    print(f"   Mean: {filtered_df['n_cells'].mean():.0f}")
    print(f"   Median: {filtered_df['n_cells'].median():.0f}")
    print(f"   Range: {filtered_df['n_cells'].min()} - {filtered_df['n_cells'].max()}")
    
    # Genes per cell
    print(f"\n3. Genes per Cell (Median):")
    print(f"   Mean: {filtered_df['genes_per_cell_median'].mean():.0f}")
    print(f"   Median: {filtered_df['genes_per_cell_median'].median():.0f}")
    print(f"   Range: {filtered_df['genes_per_cell_median'].min()} - {filtered_df['genes_per_cell_median'].max()}")
    
    # UMIs per cell
    print(f"\n4. UMIs per Cell (Median):")
    print(f"   Mean: {filtered_df['umis_per_cell_median'].mean():.0f}")
    print(f"   Median: {filtered_df['umis_per_cell_median'].median():.0f}")
    print(f"   Range: {filtered_df['umis_per_cell_median'].min()} - {filtered_df['umis_per_cell_median'].max()}")
    
    # Total genes
    print(f"\n5. Total Genes (Filtered):")
    print(f"   Mean: {filtered_df['n_genes_total'].mean():.0f}")
    print(f"   Median: {filtered_df['n_genes_total'].median():.0f}")
    print(f"   Range: {filtered_df['n_genes_total'].min()} - {filtered_df['n_genes_total'].max()}")
    
    # Transcriptome coverage
    print(f"\n6. Transcriptome Coverage (Median):")
    print(f"   Mean: {filtered_df['transcriptome_coverage_median'].mean():.4f}")
    print(f"   Median: {filtered_df['transcriptome_coverage_median'].median():.4f}")
    print(f"   Range: {filtered_df['transcriptome_coverage_median'].min():.4f} - {filtered_df['transcriptome_coverage_median'].max():.4f}")
    
    # Doublet rate
    print(f"\n7. Doublet Rate:")
    print(f"   Mean: {filtered_df['doublet_rate'].mean():.4f}")
    print(f"   Median: {filtered_df['doublet_rate'].median():.4f}")
    print(f"   Range: {filtered_df['doublet_rate'].min():.4f} - {filtered_df['doublet_rate'].max():.4f}")
    
    # Platform breakdown
    print(f"\n8. Platform Distribution:")
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
    filtered_df, fastqc_r2 = load_data(args.filtered, args.fastqc)
    
    # Generate summary statistics
    create_summary_stats(filtered_df, fastqc_r2)
    
    # Create plots
    print(f"\nGenerating histogram plots...")
    create_histogram_plots(filtered_df, fastqc_r2, args.output)
    
    print(f"\nPlots saved to: {args.output}/")
    print("Analysis complete!")

if __name__ == "__main__":
    main()

# Example usage:
# python scrna_qc_histograms.py --filtered filtered_dataset.csv --fastqc fastqc_data.csv --output qc_plots