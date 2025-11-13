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

# rRNA gene dictionaries with lengths calculated from transcripts_to_genes.txt
wMel_rRNA={
    "GQX67_00940": 2772,
    "GQX67_00945": 107,
    "GQX67_05945": 1505
}

# Drosophila rRNA genes with accurate lengths from your reference
all_rRNA={
    # Mitochondrial rRNAs (keeping your original entries)
    "FBgn0013686": 1324, # Dmel mtrRNA
    "FBgn0013688": 786,  # Dmel mtrRNA
    
    # # 2S rRNA genes (all 30 bp)
    "FBgn0267496": 30,   # 2SrRNA:CR45836
    "FBgn0267500": 30,   # 2SrRNA:CR45840
    "FBgn0267503": 30,   # 2SrRNA:CR45843
    "FBgn0085765": 30,   # 2SrRNA-Psi:CR40677
    "FBgn0267518": 30,   # 2SrRNA-Psi:CR45858
    "FBgn0267524": 30,   # 2SrRNA:CR45864
    
    # # 5.8S rRNA genes (all 123 bp)
    "FBgn0267509": 123,  # 5.8SrRNA-Psi:CR45849
    "FBgn0267499": 123,  # 5.8SrRNA:CR45839
    "FBgn0267502": 123,  # 5.8SrRNA:CR45842
    "FBgn0267512": 123,  # 5.8SrRNA:CR45852
    "FBgn0267517": 123,  # 5.8SrRNA-Psi:CR45857
    "FBgn0267523": 123,  # 5.8SrRNA-Psi:CR45863
    "FBgn0250731": 123,  # 5.8SrRNA:CR40454
    "FBgn0267514": 123,  # 5.8SrRNA-Psi:CR45854
    
    # # 18S rRNA genes 
    "FBgn0085802": 1995, # 18SrRNA:CR41548
    "FBgn0267498": 1995, # 18SrRNA:CR45838
    "FBgn0267501": 1995, # 18SrRNA:CR45841
    "FBgn0267521": 1934, # 18SrRNA-Psi:CR45861
    "FBgn0085813": 1975, # 18SrRNA-Psi:CR41602 #This is the one in the dataset
    
    # # 28S rRNA genes (variable lengths)
    "FBgn0267504": 3970, # 28SrRNA:CR45844
    "FBgn0267508": 821,  # 28SrRNA-Psi:CR45848
    "FBgn0267511": 2800, # 28SrRNA-Psi:CR45851
    "FBgn0085753": 6005, # 28SrRNA-Psi:CR40596
    "FBgn0267497": 2715, # 28SrRNA:CR45837
    "FBgn0267522": 2004, # 28SrRNA-Psi:CR45862
    "FBgn0085771": 1258, # 28SrRNA-Psi:CR40741
    "FBgn0267519": 2689, # 28SrRNA-Psi:CR45859
    "FBgn0085819": 895,  # 28SrRNA-Psi:CR41609
    "FBgn0267513": 255,  # 28SrRNA-Psi:CR45853
    "FBgn0267520": 357,  # 28SrRNA-Psi:CR45860
    "FBgn0267515": 704   # 28SrRNA-Psi:CR45855
}

# Set plotting settings
sc.settings.set_figure_params(dpi=100, frameon=False)

def calculate_wolbachia_titer(adata):
    '''
    Calculate Wolbachia titer for each cell in the AnnData object using normalized data.
    The titer is calculated using relative expression levels of wMel vs Dmel rRNA genes.
    '''
    print("Calculating Wolbachia titer using normalized expression data...")
    
    print(f"Using {len(all_rRNA)} Dmel rRNA genes: {list(all_rRNA.keys())}")
    print(f"Using {len(wMel_rRNA)} wMel rRNA genes: {list(wMel_rRNA.keys())}")
    
    # Check if we need to look in var_names or gene_ids column
    if 'gene_ids' in adata.var.columns:
        # Create a mapping from gene names/indices to gene_ids
        gene_id_map = dict(zip(adata.var.index, adata.var['gene_ids']))
        
        # Find which genes are present in our datasets
        wMel_genes_present = []
        for gene_id in wMel_rRNA.keys():
            if gene_id in adata.var['gene_ids'].values:
                wMel_genes_present.append(gene_id)
                
        dmel_genes_present = []
        for gene_id in all_rRNA.keys():
            if gene_id in adata.var['gene_ids'].values:
                dmel_genes_present.append(gene_id)
        
        # Create masks for the genes
        wMel_mask = [gene_id_map.get(idx) in wMel_genes_present for idx in adata.var.index]
        dmel_mask = [gene_id_map.get(idx) in dmel_genes_present for idx in adata.var.index]
    else:
        # Use original approach with var_names
        wMel_genes_present = [gene for gene in wMel_rRNA.keys() if gene in adata.var_names]
        dmel_genes_present = [gene for gene in all_rRNA.keys() if gene in adata.var_names]
        
        # Create masks based on var_names
        wMel_mask = [gene in wMel_genes_present for gene in adata.var_names]
        dmel_mask = [gene in dmel_genes_present for gene in adata.var_names]
    
    print(f"Found {len(wMel_genes_present)} wMel rRNA genes and {len(dmel_genes_present)} Dmel rRNA genes in dataset")
    
    # Get gene indices from the masks
    wMel_indices = np.where(wMel_mask)[0]
    dmel_indices = np.where(dmel_mask)[0]
    
    # Convert sparse matrix to dense if necessary
    is_sparse = scipy.sparse.issparse(adata.X)
    
    # For normalized data, we'll use mean expression across rRNA genes
    if len(wMel_indices) > 0:
        if is_sparse:
            wMel_mean_expr = np.array(adata.X[:, wMel_indices].mean(axis=1)).flatten()
        else:
            wMel_mean_expr = np.mean(adata.X[:, wMel_indices], axis=1)
    else:
        wMel_mean_expr = np.full(adata.n_obs, -np.inf)  # Very low expression if no genes
        print("Warning: No wMel rRNA genes found in dataset")
    
    if len(dmel_indices) > 0:
        if is_sparse:
            dmel_mean_expr = np.array(adata.X[:, dmel_indices].mean(axis=1)).flatten()
        else:
            dmel_mean_expr = np.mean(adata.X[:, dmel_indices], axis=1)
    else:
        dmel_mean_expr = np.full(adata.n_obs, -np.inf)  # Very low expression if no genes
        print("Warning: No Dmel rRNA genes found in dataset")
    
    # Convert normalized expression to relative expression (softmax-like approach)
    # Add small constant to avoid issues with negative values
    min_expr = min(wMel_mean_expr.min(), dmel_mean_expr.min()) - 1
    wMel_shifted = wMel_mean_expr - min_expr
    dmel_shifted = dmel_mean_expr - min_expr
    
    # Calculate titer using shifted values
    # Only use NA if both expressions are very low (below -2 in normalized space)
    very_low_threshold = -2.0
    has_detectable_rRNA = (wMel_mean_expr > very_low_threshold) | (dmel_mean_expr > very_low_threshold)
    
    titer = np.full(adata.n_obs, np.nan)
    
    if np.any(has_detectable_rRNA):
        # For cells with detectable rRNA, calculate relative expression titer
        total_shifted = wMel_shifted + dmel_shifted
        titer[has_detectable_rRNA] = np.where(
            total_shifted[has_detectable_rRNA] > 0,
            wMel_shifted[has_detectable_rRNA] / total_shifted[has_detectable_rRNA],
            np.nan
        )
    
    # Add the titer to the AnnData object
    adata.obs['wolbachia_titer'] = titer
    adata.obs['log1p_wolbachia_titer'] = np.log1p(np.where(np.isfinite(titer), titer, 0))
    
    # Add raw expression values for reference
    adata.obs['wMel_mean_expr'] = wMel_mean_expr
    adata.obs['dmel_mean_expr'] = dmel_mean_expr
    
    # Alternative simpler titer: just use the difference in expression
    # adata.obs['wolbachia_expr_diff'] = wMel_mean_expr - dmel_mean_expr
    
    # Count cells with different expression states
    n_cells_detectable = np.sum(has_detectable_rRNA)
    n_high_wolbachia = np.sum(wMel_mean_expr > 0)  # Above average Wolbachia
    n_no_detection = np.sum(~has_detectable_rRNA)
    
    print(f"Cells with detectable rRNA: {n_cells_detectable} out of {adata.n_obs} ({n_cells_detectable/adata.n_obs*100:.2f}%)")
    print(f"Cells with above-average Wolbachia: {n_high_wolbachia} out of {adata.n_obs} ({n_high_wolbachia/adata.n_obs*100:.2f}%)")
    print(f"Cells with low rRNA detection: {n_no_detection} out of {adata.n_obs} ({n_no_detection/adata.n_obs*100:.2f}%)")
    
    print(f"Wolbachia expression range: {wMel_mean_expr.min():.3f} to {wMel_mean_expr.max():.3f}")
    print(f"Drosophila rRNA expression range: {dmel_mean_expr.min():.3f} to {dmel_mean_expr.max():.3f}")
    
    return adata

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

    # Calculate Wolbachia titer if requested
    if calculate_titer:
        combined = calculate_wolbachia_titer(combined)

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