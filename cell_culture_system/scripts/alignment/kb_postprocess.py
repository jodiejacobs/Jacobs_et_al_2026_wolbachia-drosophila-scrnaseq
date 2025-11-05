#!/usr/bin/env python3
"""
Post-processing script for kallisto bustools output
Based on: https://www.kallistobus.tools/tutorials/kb_getting_started/python/kb_intro_2_python/

This script takes an adata.h5ad file from kb count and performs standard scRNA-seq processing
including filtering, normalization, dimensionality reduction, and clustering.

Usage:
    python kb_postprocess.py --input adata.h5ad --output processed_data.h5ad --sample sample_name
"""

import argparse
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Set scanpy settings
sc.settings.verbosity = 3  # verbosity level
sc.settings.set_figure_params(dpi=80, facecolor='white')

def load_and_inspect_data(input_file, sample_name):
    """Load the h5ad file and perform initial inspection."""
    print(f"Loading data from {input_file}")
    adata = sc.read_h5ad(input_file)
    
    # Add sample information
    adata.obs['sample'] = sample_name
    
    print(f"Initial data shape: {adata.shape}")
    print(f"Variables (genes): {adata.n_vars}")
    print(f"Observations (cells): {adata.n_obs}")
    
    return adata

def calculate_qc_metrics(adata):
    """Calculate quality control metrics."""
    print("Calculating QC metrics...")
    
    # Make variable names unique (good practice)
    adata.var_names_unique()
    
    # Mitochondrial genes
    adata.var['mt'] = adata.var_names.str.startswith('mt:') | adata.var_names.str.startswith('MT-')
    
    # Ribosomal genes  
    adata.var['ribo'] = adata.var_names.str.startswith('Rp') | adata.var_names.str.startswith('RP')
    
    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    
    # Add mitochondrial and ribosomal gene percentages
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    sc.pp.calculate_qc_metrics(adata, qc_vars=['ribo'], percent_top=None, log1p=False, inplace=True)
    
    print(f"Number of mitochondrial genes: {adata.var['mt'].sum()}")
    print(f"Number of ribosomal genes: {adata.var['ribo'].sum()}")
    
    return adata

def plot_qc_metrics(adata, output_dir):
    """Plot quality control metrics."""
    print("Plotting QC metrics...")
    
    # Create QC plots
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Number of genes per cell
    sc.pl.violin(adata, ['n_genes_by_counts'], 
                 jitter=0.4, multi_panel=True, ax=axes[0,0])
    axes[0,0].set_title('Number of genes per cell')
    
    # Total counts per cell
    sc.pl.violin(adata, ['total_counts'], 
                 jitter=0.4, multi_panel=True, ax=axes[0,1])
    axes[0,1].set_title('Total counts per cell')
    
    # Mitochondrial gene percentage
    sc.pl.violin(adata, ['pct_counts_mt'], 
                 jitter=0.4, multi_panel=True, ax=axes[1,0])
    axes[1,0].set_title('Mitochondrial gene %')
    
    # Ribosomal gene percentage
    sc.pl.violin(adata, ['pct_counts_ribo'], 
                 jitter=0.4, multi_panel=True, ax=axes[1,1])
    axes[1,1].set_title('Ribosomal gene %')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'qc_metrics.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Scatter plots
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', ax=axes[0])
    axes[0].set_title('Total counts vs MT%')
    
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', ax=axes[1])
    axes[1].set_title('Total counts vs N genes')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'qc_scatter.png', dpi=300, bbox_inches='tight')
    plt.close()

def filter_cells_and_genes(adata, min_genes=200, min_cells=3, max_genes=5000, max_mt_percent=20):
    """Filter cells and genes based on QC metrics."""
    print("Filtering cells and genes...")
    
    print(f"Initial: {adata.n_obs} cells, {adata.n_vars} genes")
    
    # Filter genes that are detected in at least min_cells cells
    sc.pp.filter_genes(adata, min_cells=min_cells)
    print(f"After gene filtering: {adata.n_obs} cells, {adata.n_vars} genes")
    
    # Filter cells with too few or too many genes
    sc.pp.filter_cells(adata, min_genes=min_genes)
    print(f"After min gene filtering: {adata.n_obs} cells, {adata.n_vars} genes")
    
    # Filter cells with too many genes (potential doublets)
    adata = adata[adata.obs.n_genes_by_counts < max_genes, :]
    print(f"After max gene filtering: {adata.n_obs} cells, {adata.n_vars} genes")
    
    # Filter cells with high mitochondrial gene expression
    adata = adata[adata.obs.pct_counts_mt < max_mt_percent, :]
    print(f"After MT filtering: {adata.n_obs} cells, {adata.n_vars} genes")
    
    return adata

def normalize_and_log_transform(adata):
    """Normalize and log transform the data."""
    print("Normalizing and log transforming...")
    
    # Save raw counts
    adata.raw = adata
    
    # Normalize to 10,000 reads per cell
    sc.pp.normalize_total(adata, target_sum=1e4)
    
    # Log transform
    sc.pp.log1p(adata)
    
    return adata

def find_highly_variable_genes(adata, n_top_genes=2000):
    """Find highly variable genes."""
    print("Finding highly variable genes...")
    
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=n_top_genes)
    
    # Plot highly variable genes
    sc.pl.highly_variable_genes(adata)
    plt.savefig('highly_variable_genes.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Number of highly variable genes: {adata.var.highly_variable.sum()}")
    
    return adata

def scale_data(adata, max_value=10):
    """Scale data to unit variance."""
    print("Scaling data...")
    
    # Keep only highly variable genes for downstream analysis
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    
    # Scale to unit variance
    sc.pp.scale(adata, max_value=max_value)
    
    return adata

def perform_pca(adata, n_comps=50):
    """Perform Principal Component Analysis."""
    print("Performing PCA...")
    
    sc.tl.pca(adata, svd_solver='arpack', n_comps=n_comps)
    
    # Plot PCA
    sc.pl.pca_variance_ratio(adata, log=True, n_pcs=n_comps)
    plt.savefig('pca_variance.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    return adata

def compute_neighborhood_graph(adata, n_neighbors=10, n_pcs=40):
    """Compute the neighborhood graph."""
    print("Computing neighborhood graph...")
    
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    
    return adata

def perform_umap(adata):
    """Perform UMAP embedding."""
    print("Performing UMAP...")
    
    sc.tl.umap(adata)
    
    return adata

def perform_clustering(adata, resolution=0.5):
    """Perform Leiden clustering."""
    print("Performing clustering...")
    
    sc.tl.leiden(adata, resolution=resolution)
    
    print(f"Number of clusters: {adata.obs['leiden'].nunique()}")
    
    return adata

def plot_results(adata, output_dir, sample_name):
    """Plot final results."""
    print("Creating final plots...")
    
    # UMAP with clusters
    sc.pl.umap(adata, color=['leiden'], legend_loc='on data', 
               title=f'{sample_name} - Leiden Clusters', frameon=False, save=f'_{sample_name}_clusters.pdf')
    
    # UMAP with QC metrics
    sc.pl.umap(adata, color=['total_counts', 'n_genes_by_counts', 'pct_counts_mt'], 
               ncols=3, save=f'_{sample_name}_qc_metrics.pdf')
    
    # Move plots to output directory
    import shutil
    for file in Path('.').glob('figures/*'):
        shutil.move(str(file), output_dir / file.name)

def save_results(adata, output_file):
    """Save the processed data."""
    print(f"Saving results to {output_file}")
    adata.write(output_file)

def main():
    parser = argparse.ArgumentParser(description='Process kallisto bustools scRNA-seq data')
    parser.add_argument('--input', required=True, help='Input h5ad file from kb count')
    parser.add_argument('--output', required=True, help='Output h5ad file')
    parser.add_argument('--sample', required=True, help='Sample name')
    parser.add_argument('--min-genes', type=int, default=200, help='Minimum genes per cell')
    parser.add_argument('--min-cells', type=int, default=3, help='Minimum cells per gene')
    parser.add_argument('--max-genes', type=int, default=5000, help='Maximum genes per cell')
    parser.add_argument('--max-mt-percent', type=float, default=20, help='Maximum mitochondrial percentage')
    parser.add_argument('--n-top-genes', type=int, default=2000, help='Number of highly variable genes')
    parser.add_argument('--resolution', type=float, default=0.5, help='Clustering resolution')
    parser.add_argument('--n-neighbors', type=int, default=10, help='Number of neighbors for graph')
    parser.add_argument('--n-pcs', type=int, default=40, help='Number of PCs for neighbors')
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Set scanpy figure directory
    sc.settings.figdir = output_dir / 'figures'
    sc.settings.figdir.mkdir(exist_ok=True)
    
    # Load and process data
    adata = load_and_inspect_data(args.input, args.sample)
    adata = calculate_qc_metrics(adata)
    
    # Plot initial QC metrics
    plot_qc_metrics(adata, output_dir)
    
    # Filter data
    adata = filter_cells_and_genes(adata, args.min_genes, args.min_cells, 
                                  args.max_genes, args.max_mt_percent)
    
    # Normalize and transform
    adata = normalize_and_log_transform(adata)
    
    # Find highly variable genes
    adata = find_highly_variable_genes(adata, args.n_top_genes)
    
    # Scale data
    adata = scale_data(adata)
    
    # Dimensionality reduction and clustering
    adata = perform_pca(adata)
    adata = compute_neighborhood_graph(adata, args.n_neighbors, args.n_pcs)
    adata = perform_umap(adata)
    adata = perform_clustering(adata, args.resolution)
    
    # Plot results
    plot_results(adata, output_dir, args.sample)
    
    # Save results
    save_results(adata, args.output)
    
    print("Processing complete!")
    print(f"Final data shape: {adata.shape}")
    print(f"Number of clusters: {adata.obs['leiden'].nunique()}")

if __name__ == "__main__":
    main()
