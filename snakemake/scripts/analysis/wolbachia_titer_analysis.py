#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Wolbachia Titer Analysis Script
-------------------------------
This script analyzes the relationship between samples, batches and Wolbachia titer
in a single-cell RNA-seq dataset.

Usage:
    python wolbachia_analysis.py <path_to_adata> <output_directory>

Example:
    python wolbachia_analysis.py /path/to/adata.h5ad /path/to/output_dir

Author: Claude
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc

def main():
    # Check arguments
    if len(sys.argv) != 3:
        print(f"Usage: python {sys.argv[0]} <path_to_adata> <output_directory>")
        sys.exit(1)
    
    adata_path = sys.argv[1]
    output_dir = sys.argv[2]
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Load the AnnData object
    print(f"Loading AnnData from {adata_path}...")
    adata = sc.read(adata_path)
    print(f"Loaded AnnData with {adata.n_obs} cells and {adata.n_vars} genes")
    
    # 1. Batch to Sample mapping
    print("\n--- Batch to Sample Mapping ---")
    analyze_batch_sample_mapping(adata, output_dir)
    
    # 2. Wolbachia titer analysis
    print("\n--- Wolbachia Titer Analysis ---")
    analyze_wolbachia_titer(adata, output_dir)
    
    # 3. PCA colored by log Wolbachia titer
    print("\n--- Creating PCA Visualization ---")
    create_pca_with_titer(adata, output_dir)
    
    print(f"\nAll analyses completed. Results saved to {output_dir}")

def analyze_batch_sample_mapping(adata, output_dir):
    """Analyze the mapping between batches and samples"""
    
    # Create a DataFrame from the obs annotation
    batch_sample_df = pd.DataFrame({
        'batch': adata.obs['batch'],
        'Sample': adata.obs['Sample']
    })
    
    # Get unique batch-to-sample mappings
    batch_sample_mapping = batch_sample_df.drop_duplicates()
    
    # Save to file
    mapping_file = os.path.join(output_dir, "batch_sample_mapping.csv")
    batch_sample_mapping.to_csv(mapping_file, index=False)
    print(f"Batch to Sample mapping saved to {mapping_file}")
    
    # Print summary
    print("\nSummary of samples in each batch:")
    batch_summary = {}
    for batch in sorted(adata.obs['batch'].unique()):
        samples = adata.obs.loc[adata.obs['batch'] == batch, 'Sample'].unique()
        batch_summary[batch] = samples
        print(f"Batch {batch}: {', '.join(samples)}")
    
    # Count cells per batch and sample
    cell_counts = adata.obs.groupby(['batch', 'Sample']).size().reset_index(name='n_cells')
    counts_file = os.path.join(output_dir, "cell_counts_by_batch_sample.csv")
    cell_counts.to_csv(counts_file, index=False)
    print(f"Cell counts saved to {counts_file}")
    
    # Create a bar plot of cells per batch colored by sample
    plt.figure(figsize=(12, 6))
    
    # Get unique samples
    unique_samples = adata.obs['Sample'].unique()
    
    # Get a colormap with enough colors
    cmap = plt.cm.get_cmap('tab20', len(unique_samples))
    
    # Create a dictionary mapping samples to colors
    sample_colors = {sample: cmap(i) for i, sample in enumerate(unique_samples)}
    
    # Plot each sample as a different color
    sample_handles = []
    for sample in unique_samples:
        subset = cell_counts[cell_counts['Sample'] == sample]
        bars = plt.bar(subset['batch'].astype(str), subset['n_cells'], 
                        label=sample, color=sample_colors[sample], 
                        alpha=0.8)
        sample_handles.append(bars[0])
    
    plt.xlabel('Batch')
    plt.ylabel('Number of cells')
    plt.title('Number of cells per batch by sample')
    plt.legend(handles=sample_handles, title='Sample')
    plt.tight_layout()
    
    # Save as SVG
    cell_counts_plot_svg = os.path.join(output_dir, "cell_counts_by_batch_sample.svg")
    plt.savefig(cell_counts_plot_svg, format='svg')
    
    # Save as PDF
    cell_counts_plot_pdf = os.path.join(output_dir, "cell_counts_by_batch_sample.pdf")
    plt.savefig(cell_counts_plot_pdf, format='pdf')
    print(f"Cell counts plot saved to {cell_counts_plot_svg} and {cell_counts_plot_pdf}")
    plt.close()

def analyze_wolbachia_titer(adata, output_dir):
    """Analyze Wolbachia titer by sample and batch"""
    
    # Calculate average Wolbachia titer per sample
    titer_by_sample = adata.obs.groupby('Sample')['wolbachia_titer'].agg(['mean', 'std', 'count']).reset_index()
    titer_by_sample.columns = ['Sample', 'Mean_Wolbachia_Titer', 'Std_Wolbachia_Titer', 'Cell_Count']
    
    # Get batch for each sample
    sample_to_batch = adata.obs.groupby('Sample')['batch'].first().reset_index()
    titer_by_sample = pd.merge(titer_by_sample, sample_to_batch, on='Sample')
    
    # Sort by batch and then by mean titer
    titer_by_sample = titer_by_sample.sort_values(['batch', 'Mean_Wolbachia_Titer'], ascending=[True, False])
    
    # Save to file
    titer_file = os.path.join(output_dir, "wolbachia_titer_by_sample.csv")
    titer_by_sample.to_csv(titer_file, index=False)
    print(f"Wolbachia titer statistics saved to {titer_file}")
    
    # Print results
    print("\nAverage Wolbachia titer by sample (sorted by batch):")
    print(titer_by_sample[['Sample', 'batch', 'Mean_Wolbachia_Titer', 'Std_Wolbachia_Titer', 'Cell_Count']])
    
    # Create a bar plot
    plt.figure(figsize=(12, 6))
    
    # Create a custom colormap for the batches
    batch_values = sorted(adata.obs['batch'].unique())
    batch_cmap = plt.cm.get_cmap('Dark2', len(batch_values))
    batch_colors = {batch: batch_cmap(i) for i, batch in enumerate(batch_values)}
    
    # Plot bars with colors by batch
    for i, (_, row) in enumerate(titer_by_sample.iterrows()):
        plt.bar(i, row['Mean_Wolbachia_Titer'], 
                yerr=row['Std_Wolbachia_Titer'] / np.sqrt(row['Cell_Count']),  # Standard error
                color=batch_colors[row['batch']], 
                alpha=0.8)
    
    # Create batch legend
    batch_handles = [plt.Rectangle((0,0), 1, 1, color=batch_colors[batch]) for batch in batch_values]
    plt.legend(batch_handles, [f'Batch {batch}' for batch in batch_values], 
               title='Batch', loc='upper right')
    
    # Customize x-axis with sample names
    plt.xticks(range(len(titer_by_sample)), titer_by_sample['Sample'], rotation=45, ha='right')
    plt.xlabel('Sample')
    plt.ylabel('Mean Wolbachia Titer')
    plt.title('Average Wolbachia Titer by Sample')
    plt.tight_layout()
    
    # Save as SVG
    titer_bar_plot_svg = os.path.join(output_dir, "wolbachia_titer_barplot.svg")
    plt.savefig(titer_bar_plot_svg, format='svg')
    
    # Save as PDF
    titer_bar_plot_pdf = os.path.join(output_dir, "wolbachia_titer_barplot.pdf")
    plt.savefig(titer_bar_plot_pdf, format='pdf')
    print(f"Wolbachia titer bar plot saved to {titer_bar_plot_svg} and {titer_bar_plot_pdf}")
    plt.close()
    
    # Create a box plot showing the distribution within each sample
    plt.figure(figsize=(14, 8))
    
    # Convert batch to string for better visualization
    plot_data = adata.obs.copy()
    plot_data['batch'] = plot_data['batch'].astype(str)
    
    # Create box plot
    sns.boxplot(x='Sample', y='wolbachia_titer', hue='batch', data=plot_data)
    
    plt.xticks(rotation=45, ha='right')
    plt.title('Distribution of Wolbachia Titer by Sample')
    plt.xlabel('Sample')
    plt.ylabel('Wolbachia Titer')
    plt.tight_layout()
    
    # Save as SVG
    titer_box_plot_svg = os.path.join(output_dir, "wolbachia_titer_boxplot.svg")
    plt.savefig(titer_box_plot_svg, format='svg')
    
    # Save as PDF
    titer_box_plot_pdf = os.path.join(output_dir, "wolbachia_titer_boxplot.pdf")
    plt.savefig(titer_box_plot_pdf, format='pdf')
    print(f"Wolbachia titer box plot saved to {titer_box_plot_svg} and {titer_box_plot_pdf}")
    plt.close()

def create_pca_with_titer(adata, output_dir):
    """Create PCA visualization colored by log Wolbachia titer"""
    
    # Check if PCA exists in the object
    if 'X_pca' not in adata.obsm:
        print("Warning: PCA not found in the AnnData object. Computing PCA...")
        sc.pp.pca(adata)
    
    # Check if log1p_wolbachia_titer exists, otherwise create it
    if 'log1p_wolbachia_titer' not in adata.obs:
        print("Adding log1p transformed Wolbachia titer...")
        adata.obs['log1p_wolbachia_titer'] = np.log1p(adata.obs['wolbachia_titer'])
    
    # Create PCA plot colored by log1p Wolbachia titer
    plt.figure(figsize=(10, 8))
    
    # Get the PCA coordinates
    pca_coords = adata.obsm['X_pca']
    
    # Create a scatter plot colored by log1p_wolbachia_titer
    scatter = plt.scatter(pca_coords[:, 0], pca_coords[:, 1], 
                         c=adata.obs['log1p_wolbachia_titer'], 
                         cmap='viridis', 
                         alpha=0.7,
                         s=5)  # Point size
    
    # Add a colorbar
    cbar = plt.colorbar(scatter)
    cbar.set_label('log1p(Wolbachia Titer)')
    
    # Set labels and title
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.title('PCA colored by log1p(Wolbachia Titer)')
    
    # Save the plot
    plt.tight_layout()
    
    # Save as SVG
    pca_plot_svg = os.path.join(output_dir, "pca_wolbachia_titer.svg")
    plt.savefig(pca_plot_svg, format='svg')
    
    # Save as PDF
    pca_plot_pdf = os.path.join(output_dir, "pca_wolbachia_titer.pdf")
    plt.savefig(pca_plot_pdf, format='pdf')
    print(f"PCA plot saved to {pca_plot_svg} and {pca_plot_pdf}")
    plt.close()
    
    # Create an additional PCA plot split by batch
    batch_values = sorted(adata.obs['batch'].unique())
    
    # Determine grid layout (roughly square)
    n_batches = len(batch_values)
    n_cols = int(np.ceil(np.sqrt(n_batches)))
    n_rows = int(np.ceil(n_batches / n_cols))
    
    # Create the figure
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
    axes = axes.flatten() if n_batches > 1 else [axes]
    
    # Plot each batch separately
    for i, batch in enumerate(batch_values):
        # Get data for this batch
        mask = adata.obs['batch'] == batch
        
        # Plot
        scatter = axes[i].scatter(
            pca_coords[mask, 0], pca_coords[mask, 1],
            c=adata.obs.loc[mask, 'log1p_wolbachia_titer'],
            cmap='viridis',
            alpha=0.7,
            s=5
        )
        
        # Add a title
        axes[i].set_title(f'Batch {batch}')
        
        # Add labels
        axes[i].set_xlabel('PC1')
        axes[i].set_ylabel('PC2')
        
        # Add colorbar
        fig.colorbar(scatter, ax=axes[i], label='log1p(Wolbachia Titer)')
    
    # Hide any unused subplots
    for j in range(i+1, len(axes)):
        axes[j].axis('off')
    
    # Add a title to the overall figure
    fig.suptitle('PCA by Batch Colored by log1p(Wolbachia Titer)', fontsize=16)
    fig.tight_layout(rect=[0, 0, 1, 0.97])  # Adjust for the suptitle
    
    # Save the plot
    # Save as SVG
    pca_batch_plot_svg = os.path.join(output_dir, "pca_wolbachia_titer_by_batch.svg")
    plt.savefig(pca_batch_plot_svg, format='svg')
    
    # Save as PDF
    pca_batch_plot_pdf = os.path.join(output_dir, "pca_wolbachia_titer_by_batch.pdf")
    plt.savefig(pca_batch_plot_pdf, format='pdf')
    print(f"PCA by batch plot saved to {pca_batch_plot_svg} and {pca_batch_plot_pdf}")
    plt.close()

if __name__ == "__main__":
    main()
