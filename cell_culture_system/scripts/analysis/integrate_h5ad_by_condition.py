import os
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from glob import glob
import argparse
import scipy.sparse
import re
import anndata as ad
import bbknn
import harmonypy
import warnings

def integrate(files, out_path, fig_dir, sample, batch_key, min_cells, min_genes, calculate_titer=True, n_pcs=30):
    """
    files: list of h5ad file paths
    """
    # Make output directory if it doesn't exist
    os.makedirs(fig_dir, exist_ok=True)

    # Load all files
    adatas = []
    for file_path in files:
        adata = sc.read_h5ad(file_path)
        # Add batch information based on filename
        batch_name = os.path.splitext(os.path.basename(file_path))[0]
        adata.obs[batch_key] = batch_name
        adatas.append(adata)

    # Set output directory for figures
    sc.settings.figdir = fig_dir

    combined = ad.concat(adatas, join='inner', merge='same', index_unique='-')

    # # Calculate Wolbachia titer if requested
    # if calculate_titer:
    #     combined = calculate_wolbachia_titer(combined)

    print(f"Combined data shape for {sample}: {combined.shape}")
    print(combined)
    print(f"Data range: {combined.X.min():.3f} to {combined.X.max():.3f}")
    print(combined.X)
    
    # Filter out bacterial genes
    bacteria_genes = ['GQX67_00940', 'GQX67_05945'] + [gene for gene in combined.var_names if gene.startswith('16S_')]
    
    bacteria_mask = combined.var_names.isin(bacteria_genes)

    combined = combined[:, ~bacteria_mask]

    # Basic preprocessing
    print("Performing basic preprocessing...")
    combined.X = np.nan_to_num(combined.X, nan=0.0)  # Replace NaN with 0

    sc.pp.filter_cells(combined, min_genes=min_genes)
    sc.pp.filter_genes(combined, min_cells=min_cells)
    
    # Calculate QC metrics
    # sc.pp.calculate_qc_metrics(combined, inplace=True)
    
    # # Normalize the data
    # print("Normalizing data...")
    # sc.pp.normalize_total(combined, target_sum=1e4)
    # sc.pp.log1p(combined)

    # print("Data after normalization:")
    # print(combined)
    # print(f"Data range: {combined.X.min():.3f} to {combined.X.max():.3f}")
          
    # Find highly variable genes
    print("Finding highly variable genes...")
    print(f"Data shape before HVG: {combined.shape}")
    print(f"Data range: {combined.X.min():.3f} to {combined.X.max():.3f}")

    # Use seurat method which is more robust
    sc.pp.highly_variable_genes(combined, flavor='seurat', n_top_genes=2000)
    combined = combined[:, combined.var.highly_variable]
        
    # Run PCA
    print("Running PCA...")
    sc.pp.pca(combined, n_comps=n_pcs)
    
    # Save a copy of the unintegrated data for comparison
    combined_unintegrated = combined.copy()
    sc.pp.neighbors(combined_unintegrated, n_pcs=n_pcs)
    sc.tl.umap(combined_unintegrated)
    sc.pl.umap(combined_unintegrated, color=batch_key, save=f'_{sample}_before_batch_correction.pdf')

    bbknn.bbknn(combined, batch_key=batch_key, n_pcs=n_pcs, neighbors_within_batch=5)
            
    # Run UMAP and clustering
    print("Running UMAP and Leiden clustering...")
    sc.tl.umap(combined)
    sc.tl.leiden(combined, resolution=0.8)
    
    # Save the integrated object
    print(f"Saving integrated object for {sample} to {out_path}")
    combined.write(out_path)
    
    # Generate diagnostic plots
    print("Generating diagnostic plots...")
    sc.pl.umap(combined, color=batch_key, save=f'_{sample}_bbknn.pdf')
    sc.pl.umap(combined, color='leiden', save=f'_{sample}_bbknn_leiden.pdf')
    
    # If wolbachia_titer exists, plot it too
    if 'wolbachia_titer' in combined.obs.columns:
        sc.pl.umap(combined, color='wolbachia_titer', save=f'_{sample}_bbknn_titer.pdf')
        sc.pl.umap(combined, color='log1p_wolbachia_titer', save=f'_{sample}_log1p_bbknn_titer.pdf')
        
        # Create a violin plot of titer by batch
        sc.pl.violin(combined, 'wolbachia_titer', groupby=batch_key, save=f'_{sample}_wolbachia_titer_by_rep.pdf')
        
        # Create a violin plot of titer by cluster
        sc.pl.violin(combined, 'wolbachia_titer', groupby='leiden', save=f'_{sample}_wolbachia_titer_by_cluster.pdf')
    
    print(f"Integration complete for sample type {sample}!")
    
    # Print summary for this sample type
    print(f"Summary of integrated data for {sample}:")
    print(f"Number of cells: {combined.n_obs}")
    print(f"Number of genes: {combined.n_vars}")
    print(f"Number of batches: {combined.obs[batch_key].nunique()}")
    print(f"Number of clusters: {combined.obs['leiden'].nunique()}")
    
    if 'wolbachia_titer' in combined.obs.columns:
        # Calculate percentage of infected cells
        n_infected = np.sum(combined.obs['wolbachia_titer'] > 0)
        print(f"Number of cells with Wolbachia: {n_infected} ({n_infected/combined.n_obs*100:.2f}%)")
        
        # Calculate average titer
        mean_titer = np.nanmean(combined.obs['wolbachia_titer'])
        median_titer = np.nanmedian(combined.obs['wolbachia_titer'])
        print(f"Average Wolbachia titer: mean={mean_titer:.4f}, median={median_titer:.4f}")
        
        # Calculate titer by batch
        for batch in combined.obs[batch_key].unique():
            batch_cells = combined[combined.obs[batch_key] == batch]
            n_batch_infected = np.sum(batch_cells.obs['wolbachia_titer'] > 0)
            mean_batch_titer = np.nanmean(batch_cells.obs['wolbachia_titer'])
            print(f"  {batch}: {n_batch_infected}/{batch_cells.n_obs} cells infected ({n_batch_infected/batch_cells.n_obs*100:.2f}%), mean titer={mean_batch_titer:.4f}")


def main():
    parser = argparse.ArgumentParser(description='Integrate h5ad files by sample type with batch correction')

    parser.add_argument('--files', required=True, nargs='+', type=str,
                        help='List of h5ad files to integrate')
    parser.add_argument('--sample', type=str, default='NA',
                        help='Sample type (e.g., Infected, Uninfected)')
    parser.add_argument('--batch_key', type=str, default='batch',
                        help='Key in .obs to use for batch information')
    parser.add_argument('--min_cells', type=int, default=3,
                        help='Minimum cells per gene for filtering')
    parser.add_argument('--min_genes', type=int, default=200,
                        help='Minimum genes per cell for filtering')  
    parser.add_argument('--out_path', type=str, default='test_integrated.h5ad',
                        help='Path to save the integrated h5ad file')      
    parser.add_argument('--fig_dir', type=str, default='figures',
                        help='Directory to save figures')    

    args = parser.parse_args()
    
    # Run the integration with list of files
    integrate(
        files=args.files,  # Pass the list of files
        out_path=args.out_path,
        fig_dir=args.fig_dir,
        sample=args.sample, 
        batch_key=args.batch_key,
        min_cells=args.min_cells,
        min_genes=args.min_genes,
    )

if __name__ == "__main__":
    main()