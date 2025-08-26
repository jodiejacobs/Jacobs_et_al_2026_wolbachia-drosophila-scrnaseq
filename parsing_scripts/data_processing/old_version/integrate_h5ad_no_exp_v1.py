import os
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from glob import glob
import argparse
import scipy.sparse

# rRNA gene dictionaries with lengths
wMel_rRNA={
    "GQX67_00940": 2772,
    "GQX67_00945": 107,
    "GQX67_05945": 1505
}

Dmel_rRNA={
    "FBgn0085802": 1995,
    "FBgn0267504": 3970,
    "FBgn0267512": 123,
    "FBgn0267500": 30,
    "FBgn0267508": 821,
    "FBgn0267509": 123,
    "FBgn0267510": 30,
    "FBgn0267498": 1995,
    "FBgn0267502": 123,
    "FBgn0267503": 30,
    "FBgn0267505": 6863,
    "FBgn0267506": 4191,
    "FBgn0267511": 1721,
    "FBgn0085753": 2775,
    "FBgn0267507": 8118,
    "FBgn0267501": 1995,
    "FBgn0267496": 30,
    "FBgn0267497": 2715,
    "FBgn0267499": 123
}

Dmel_mt={
    "FBgn0013686": 1324,
    "FBgn0013688": 786
}

# Set plotting settings
sc.settings.set_figure_params(dpi=100, frameon=False)
sc.settings.verbosity = 1

def calculate_wolbachia_titer(adata, min_total_reads=100, min_wolb_reads=3, depth_normalize=True):
    """
    Calculate Wolbachia titer for each cell in the AnnData object with depth bias correction.
    
    Parameters:
    -----------
    adata : AnnData
        AnnData object containing gene expression data
    min_total_reads : int
        Minimum total reads (sum of all genes) required to consider a cell for titer calculation
    min_wolb_reads : int
        Minimum Wolbachia rRNA reads required to consider a cell infected
    depth_normalize : bool
        Whether to normalize titer by total sequencing depth
        
    Returns:
    --------
    adata : AnnData
        The input AnnData object with additional columns for Wolbachia titer
    """
    print("Calculating Wolbachia titer with depth bias correction...")
    
    # Check if we need to look in var_names or gene_ids column
    if 'gene_ids' in adata.var.columns:
        # Create a mapping from gene names/indices to gene_ids
        gene_id_map = dict(zip(adata.var.index, adata.var['gene_ids']))
        
        # Find which genes are present in our dictionaries
        wMel_genes_present = []
        for gene_id in wMel_rRNA.keys():
            # Check if this gene ID is in the gene_ids column
            if gene_id in adata.var['gene_ids'].values:
                wMel_genes_present.append(gene_id)
                
        Dmel_genes_present = []
        for gene_id in Dmel_rRNA.keys():
            # Check if this gene ID is in the gene_ids column
            if gene_id in adata.var['gene_ids'].values:
                Dmel_genes_present.append(gene_id)
        
        # Create masks for the genes
        wMel_mask = [gene_id_map.get(idx, "") in wMel_genes_present for idx in adata.var.index]
        Dmel_mask = [gene_id_map.get(idx, "") in Dmel_genes_present for idx in adata.var.index]
    else:
        # Use original approach with var_names
        wMel_genes_present = [gene for gene in wMel_rRNA.keys() if gene in adata.var_names]
        Dmel_genes_present = [gene for gene in Dmel_rRNA.keys() if gene in adata.var_names]
        
        # Create masks based on var_names
        wMel_mask = [gene in wMel_genes_present for gene in adata.var_names]
        Dmel_mask = [gene in Dmel_genes_present for gene in adata.var_names]
    
    if not wMel_genes_present:
        print("Warning: No wMel rRNA genes found in the dataset!")
        adata.obs['wolbachia_titer'] = np.nan
        return adata
        
    if not Dmel_genes_present:
        print("Warning: No Dmel rRNA genes found in the dataset!")
        adata.obs['wolbachia_titer'] = np.nan
        return adata
    
    print(f"Found {len(wMel_genes_present)} wMel rRNA genes and {len(Dmel_genes_present)} Dmel rRNA genes")
    
    # Get gene indices from the masks
    wMel_indices = np.where(wMel_mask)[0]
    Dmel_indices = np.where(Dmel_mask)[0]
    
    # Convert sparse matrix to dense if necessary
    is_sparse = scipy.sparse.issparse(adata.X)
    
    # Calculate total reads per cell (for sequencing depth normalization)
    if is_sparse:
        total_reads_per_cell = np.array(adata.X.sum(axis=1)).flatten()
    else:
        total_reads_per_cell = np.sum(adata.X, axis=1)
    
    # Store the total reads in the AnnData object
    adata.obs['total_reads'] = total_reads_per_cell
    
    # Calculate counts per length for each gene
    wMel_counts_per_length = np.zeros((adata.n_obs, len(wMel_indices)))
    Dmel_counts_per_length = np.zeros((adata.n_obs, len(Dmel_indices)))
    
    # Calculate raw counts for Wolbachia and Drosophila (for thresholding)
    wMel_raw_counts = np.zeros(adata.n_obs)
    Dmel_raw_counts = np.zeros(adata.n_obs)
    
    # Calculate counts per length for wMel genes
    for i, idx in enumerate(wMel_indices):
        gene_idx = adata.var.index[idx]
        gene_id = gene_id_map.get(gene_idx, gene_idx) if 'gene_ids' in adata.var.columns else gene_idx
        length = wMel_rRNA[gene_id]
        
        if is_sparse:
            counts = adata.X[:, idx].toarray().flatten()
        else:
            counts = adata.X[:, idx]
            
        wMel_counts_per_length[:, i] = counts / length
        wMel_raw_counts += counts
    
    # Calculate counts per length for Dmel genes
    for i, idx in enumerate(Dmel_indices):
        gene_idx = adata.var.index[idx]
        gene_id = gene_id_map.get(gene_idx, gene_idx) if 'gene_ids' in adata.var.columns else gene_idx
        length = Dmel_rRNA[gene_id]
        
        if is_sparse:
            counts = adata.X[:, idx].toarray().flatten()
        else:
            counts = adata.X[:, idx]
            
        Dmel_counts_per_length[:, i] = counts / length
        Dmel_raw_counts += counts
    
    # Store raw counts in the AnnData object
    adata.obs['wMel_raw_counts'] = wMel_raw_counts
    adata.obs['Dmel_raw_counts'] = Dmel_raw_counts
    
    # Calculate the mean counts per length for each organism and cell
    wMel_mean_per_cell = np.mean(wMel_counts_per_length, axis=1)
    Dmel_mean_per_cell = np.mean(Dmel_counts_per_length, axis=1)
    
    # Create a mask for cells with sufficient sequencing depth and Wolbachia reads
    depth_mask = total_reads_per_cell >= min_total_reads
    wolb_mask = wMel_raw_counts >= min_wolb_reads
    
    # Calculate the titer with depth normalization if requested
    with np.errstate(divide='ignore', invalid='ignore'):
        if depth_normalize:
            # Calculate titer as ratio of proportion of reads
            wMel_proportion = wMel_mean_per_cell / total_reads_per_cell
            Dmel_proportion = Dmel_mean_per_cell / total_reads_per_cell
            titer = np.where(Dmel_proportion > 0, 
                            wMel_proportion / Dmel_proportion, 
                            np.nan)
        else:
            # Calculate raw titer (original method)
            titer = np.where(Dmel_mean_per_cell > 0, 
                            wMel_mean_per_cell / Dmel_mean_per_cell, 
                            np.nan)
    
    # Apply the depth and Wolbachia read thresholds
    # Set titer to 0 for cells that don't meet criteria (instead of NaN)
    titer_with_threshold = titer.copy()
    
    # First handle NaN values (where Dmel_mean_per_cell is 0)
    titer_with_threshold[np.isnan(titer_with_threshold)] = 0
    
    # Then apply the thresholds - if cell doesn't meet the criteria, set titer to 0
    titer_with_threshold[~(depth_mask & wolb_mask)] = 0
    
    # Add the metrics to the AnnData object
    adata.obs['wolbachia_titer_raw'] = titer  # Original titer without thresholds (may contain NaN)
    adata.obs['wolbachia_titer'] = titer_with_threshold  # Titer with thresholds applied and NaN set to 0
    adata.obs['log1p_wolbachia_titer'] = np.log1p(titer_with_threshold)  # This works with zeros
    adata.obs['wolbachia_detected'] = wolb_mask & depth_mask  # Boolean indicator of infection
    
    # Calculate additional metrics for QC
    adata.obs['wMel_to_total_ratio'] = wMel_raw_counts / total_reads_per_cell
    
    # Count cells with Wolbachia
    n_infected = np.sum(titer_with_threshold > 0)
    n_filtered = np.sum(~depth_mask)
    print(f"Filtered out {n_filtered} cells with < {min_total_reads} total reads")
    print(f"Detected Wolbachia in {n_infected} out of {adata.n_obs} cells ({n_infected/adata.n_obs*100:.2f}%)")
    
    # Plot detection metrics
    try:
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm
        import seaborn as sns
        
        # Create a color palette based on unique samples if 'Sample' exists in obs
        if 'Sample' in adata.obs.columns:
            samples = adata.obs['Sample'].unique()
            sample_colors = sns.color_palette("husl", len(samples))
            sample_color_dict = dict(zip(samples, sample_colors))
            
            # Create array of colors for each cell based on its sample
            cell_colors = [sample_color_dict[s] for s in adata.obs['Sample']]
            
            # Create plots with sample-based coloring
            fig, axs = plt.subplots(1, 3, figsize=(18, 5))
            
            # Plot 1: Wolbachia counts vs. total counts
            for sample in samples:
                sample_mask = adata.obs['Sample'] == sample
                axs[0].scatter(
                    total_reads_per_cell[sample_mask], 
                    wMel_raw_counts[sample_mask],
                    label=sample,
                    alpha=0.6, s=15
                )
                
            axs[0].set_xscale('log')
            axs[0].set_yscale('log')
            axs[0].set_xlabel('Total reads per cell')
            axs[0].set_ylabel('Wolbachia rRNA reads')
            axs[0].axvline(x=min_total_reads, color='r', linestyle='--')
            axs[0].axhline(y=min_wolb_reads, color='r', linestyle='--')
            axs[0].set_title('Wolbachia detection by Sample')
            axs[0].legend(title='Sample', bbox_to_anchor=(1.05, 1), loc='upper left')
            
            # Plot 2: Wolbachia ratio vs. total counts
            for sample in samples:
                sample_mask = adata.obs['Sample'] == sample
                axs[1].scatter(
                    total_reads_per_cell[sample_mask], 
                    wMel_raw_counts[sample_mask] / total_reads_per_cell[sample_mask],
                    label=sample,
                    alpha=0.6, s=15
                )
                
            axs[1].set_xscale('log')
            axs[1].set_yscale('log')
            axs[1].set_xlabel('Total reads per cell')
            axs[1].set_ylabel('Wolbachia reads / Total reads')
            axs[1].axvline(x=min_total_reads, color='r', linestyle='--')
            axs[1].set_title('Wolbachia proportion by Sample')
            # Don't repeat the legend on the second plot
            
            # Plot 3: Histogram of titer values by sample
            for i, sample in enumerate(samples):
                sample_mask = adata.obs['Sample'] == sample
                sample_titer = titer_with_threshold[sample_mask]
                sample_titer_nonzero = sample_titer[sample_titer > 0]  # Only plot non-zero values
                
                if len(sample_titer_nonzero) > 0:  # Only plot if there are valid titer values
                    axs[2].hist(sample_titer_nonzero, bins=30, alpha=0.5, 
                               label=f"{sample} ({len(sample_titer_nonzero)}/{sum(sample_mask)} cells)", 
                               color=sample_colors[i])
                    
            axs[2].set_xlabel('Wolbachia titer (only showing values > 0)')
            axs[2].set_ylabel('Number of cells')
            axs[2].set_title('Titer distribution by Sample')
            axs[2].legend(title='Sample', bbox_to_anchor=(1.05, 1), loc='upper left')
            
        else:
            # Original plotting code if 'Sample' column doesn't exist
            fig, axs = plt.subplots(1, 3, figsize=(15, 5))
            
            # Plot 1: Wolbachia counts vs. total counts
            sc = axs[0].scatter(total_reads_per_cell, wMel_raw_counts, 
                               c=wolb_mask & depth_mask, cmap='viridis', 
                               alpha=0.6, s=10)
            axs[0].set_xscale('log')
            axs[0].set_yscale('log')
            axs[0].set_xlabel('Total reads per cell')
            axs[0].set_ylabel('Wolbachia rRNA reads')
            axs[0].axvline(x=min_total_reads, color='r', linestyle='--')
            axs[0].axhline(y=min_wolb_reads, color='r', linestyle='--')
            axs[0].set_title('Wolbachia detection')
            
            # Plot 2: Wolbachia ratio vs. total counts
            sc = axs[1].scatter(total_reads_per_cell, wMel_raw_counts / total_reads_per_cell, 
                               c=wolb_mask & depth_mask, cmap='viridis', 
                               alpha=0.6, s=10)
            axs[1].set_xscale('log')
            axs[1].set_yscale('log')
            axs[1].set_xlabel('Total reads per cell')
            axs[1].set_ylabel('Wolbachia reads / Total reads')
            axs[1].axvline(x=min_total_reads, color='r', linestyle='--')
            axs[1].set_title('Wolbachia proportion')
            
            # Plot 3: Histogram of titer values
            titer_nonzero = titer_with_threshold[titer_with_threshold > 0]
            if len(titer_nonzero) > 0:
                axs[2].hist(titer_nonzero, bins=50, alpha=0.7)
                axs[2].set_xlabel('Wolbachia titer (only showing values > 0)')
                axs[2].set_ylabel('Number of cells')
                axs[2].set_title(f'Titer distribution ({len(titer_nonzero)}/{len(titer_with_threshold)} cells)')
            else:
                axs[2].text(0.5, 0.5, "No cells with Wolbachia detected", 
                           horizontalalignment='center', verticalalignment='center',
                           transform=axs[2].transAxes, fontsize=12)
        
        plt.tight_layout()
        plt.savefig('wolbachia_detection_metrics.pdf')
        plt.close()
        
    except Exception as e:
        print(f"Could not generate diagnostic plots: {e}")
    
    return adata

def integrate_h5ad_files(directory_path, output_path, batch_key='batch', 
                         min_cells=3, min_genes=200, n_pcs=30, n_neighbors=15,
                         method='bbknn', calculate_titer=True,
                         min_total_reads=100, min_wolb_reads=3, depth_normalize=True):
    """
    Integrate multiple h5ad files with batch correction
    
    Parameters:
    -----------
    directory_path : str
        Path to directory containing h5ad files
    output_path : str
        Path to save integrated h5ad file
    batch_key : str
        Name of column to use for batch correction
    min_cells : int
        Minimum number of cells expressing a gene
    min_genes : int
        Minimum number of genes expressed in a cell
    n_pcs : int
        Number of principal components to use
    n_neighbors : int
        Number of neighbors for neighborhood graph
    method : str
        Batch correction method: 'bbknn', 'harmony', or 'both'
    calculate_titer : bool
        Whether to calculate Wolbachia titer before integration
    min_total_reads : int
        Minimum total reads required to consider a cell for Wolbachia titer calculation
    min_wolb_reads : int
        Minimum Wolbachia rRNA reads required to consider a cell infected
    depth_normalize : bool
        Whether to normalize titer by total sequencing depth
    """
    # Set the figure directory to the directory containing the output file
    output_dir = os.path.dirname(output_path)
    figure_dir = os.path.join(output_dir, "figures")
    # Make sure the directory exists
    os.makedirs(figure_dir, exist_ok=True)
    # Set Scanpy figure directory
    sc.settings.figdir = figure_dir
    # Get all h5ad files in the directory
    
    h5ad_files = glob(os.path.join(directory_path, "*.h5ad"))
    
    if len(h5ad_files) == 0:
        print(f"No h5ad files found in {directory_path}")
        return None
    
    print(f"Found {len(h5ad_files)} h5ad files")
    
    # List to store individual datasets
    adatas = []
    
    # Load each dataset and add batch information
    for i, file_path in enumerate(h5ad_files):
        file_name = os.path.basename(file_path)
        batch_id = file_name.split('.')[0]  # Use filename without extension as batch ID
        print(f"Processing file {i+1}/{len(h5ad_files)}: {file_name}")
        
        # Load the dataset
        adata = sc.read_h5ad(file_path)
        adata.obs['Sample'] = os.path.splitext(os.path.basename(file_path))[0]
        
        # Check if var index contains tab characters and fix if needed
        if any('\t' in idx for idx in adata.var_names):
            print("Fixing var index with tab characters...")
            # Extract gene IDs and gene names from the index
            gene_ids = []
            gene_names = []
            
            for idx in adata.var_names:
                if '\t' in idx:
                    parts = idx.split('\t')
                    gene_id = parts[0]
                    gene_name = parts[1]
                    gene_ids.append(gene_id)
                    gene_names.append(gene_name)
                else:
                    # If there's no tab, use the same value for both id and name
                    gene_ids.append(idx)
                    gene_names.append(idx)
            
            # Create a new var DataFrame with gene names as index
            new_var = pd.DataFrame(index=gene_names)
            new_var['gene_ids'] = gene_ids
            new_var['feature_types'] = 'Gene Expression'
            
            # Copy over other columns from original var
            for col in adata.var.columns:
                new_var[col] = adata.var[col].values
            
            # Create a new AnnData object with fixed var
            # We need to get the X matrix and other components
            import anndata as ad
            new_adata = ad.AnnData(
                X=adata.X,
                obs=adata.obs,
                var=new_var,
                uns=adata.uns,
                obsm=adata.obsm if hasattr(adata, 'obsm') else None,
                varm=adata.varm if hasattr(adata, 'varm') else None,
                obsp=adata.obsp if hasattr(adata, 'obsp') else None,
                varp=adata.varp if hasattr(adata, 'varp') else None
            )
            adata = new_adata
        
        # Add batch information based on filename
        batch_id = file_name.split('.')[0]  # Use filename without extension as batch ID
        adata.obs[batch_key] = batch_id
        
        # Calculate Wolbachia titer if requested
        if calculate_titer:
            adata = calculate_wolbachia_titer(adata, 
                                            min_total_reads=min_total_reads, 
                                            min_wolb_reads=min_wolb_reads, 
                                            depth_normalize=depth_normalize)
            
        adatas.append(adata)
        
    # Concatenate all datasets
    print("Concatenating datasets...")
    try:
        # Try using anndata.concat as recommended by the FutureWarning
        import anndata as ad
        combined = ad.concat(adatas, join='outer', merge='same', label=batch_key, index_unique='-')
    except:
        # Fall back to concatenate method if concat fails
        combined = adatas[0].concatenate(adatas[1:], join='outer', batch_key=batch_key)

    print(f"Combined data shape: {combined.shape}")
    
    # Basic preprocessing
    print("Performing basic preprocessing...")
    sc.pp.filter_cells(combined, min_genes=min_genes)
    sc.pp.filter_genes(combined, min_cells=min_cells)
    
    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(combined, inplace=True)
    
    # Normalize the data
    print("Normalizing data...")
    sc.pp.normalize_total(combined, target_sum=1e4)
    sc.pp.log1p(combined)
    
    # Find highly variable genes
    print("Finding highly variable genes...")
    sc.pp.highly_variable_genes(combined, batch_key=batch_key)
    combined = combined[:, combined.var.highly_variable]
    
    # Run PCA
    print("Running PCA...")
    sc.pp.pca(combined, n_comps=n_pcs)
    
    # Save a copy of the unintegrated data for comparison
    combined_unintegrated = combined.copy()
    sc.pp.neighbors(combined_unintegrated, n_pcs=n_pcs)
    sc.tl.umap(combined_unintegrated)
    sc.pl.umap(combined_unintegrated, color=batch_key, save='_before_batch_correction.pdf')
    
    # Optional: Also plot Wolbachia titer on unintegrated data
    if 'wolbachia_titer' in combined_unintegrated.obs.columns:
        sc.pl.umap(combined_unintegrated, color='wolbachia_titer', save='_unintegrated_wolbachia_titer.pdf')
    
    # Run batch correction using the specified method
    if method == 'bbknn' or method == 'both':
        print("Performing batch correction with BBKNN...")
        try:
            import bbknn
            bbknn.bbknn(combined, batch_key=batch_key, n_pcs=n_pcs, neighbors_within_batch=5)
            
            # Save BBKNN-corrected object
            if method == 'both':
                combined_bbknn = combined.copy()
        except ImportError:
            print("BBKNN package not available. Please install it with: pip install bbknn")
            if method == 'bbknn':
                print("Falling back to standard PCA without batch correction.")
                sc.pp.neighbors(combined, n_pcs=n_pcs)
    
    if method == 'harmony' or method == 'both':
        print("Performing batch correction with Harmony...")
        # If we're using both methods, we need to start from scratch for harmony
        if method == 'both':
            combined = combined_unintegrated.copy()
        
        try:
            sc.external.pp.harmony_integrate(combined, batch_key)
            # Use harmony embedding for downstream analysis
            combined.obsm['X_pca_harmony'] = combined.obsm['X_pca_harmony']
            sc.pp.neighbors(combined, use_rep='X_pca_harmony')
        except ImportError:
            print("Harmony package not available. Please install it with: pip install harmonypy")
            if method == 'harmony':
                print("Falling back to standard PCA without batch correction.")
                sc.pp.neighbors(combined, n_pcs=n_pcs)
    
    # Run UMAP 
    print("Running UMAP and Leiden clustering...")
    sc.tl.umap(combined)
    sc.tl.leiden(combined, resolution=0.8)
    
    # Save the integrated object
    print(f"Saving integrated object to {output_path}")
    combined.write(output_path)
    
    # Generate diagnostic plots
    print("Generating diagnostic plots...")
    sc.pl.umap(combined, color=batch_key, save='_batch.pdf')
    sc.pl.umap(combined, color='leiden', save='_leiden.pdf')
    
    # If wolbachia_titer exists, plot it too
    if 'wolbachia_titer' in combined.obs.columns:
        sc.pl.umap(combined, color='wolbachia_titer', save='_wolbachia_titer.pdf')
        sc.pl.umap(combined, color='log1p_wolbachia_titer', save='_log1p_wolbachia_titer.pdf')
        
        # Create a violin plot of titer by batch
        sc.pl.violin(combined, 'wolbachia_titer', groupby=batch_key, save='_wolbachia_titer_by_batch.pdf')
        
        # Create a violin plot of titer by cluster
        sc.pl.violin(combined, 'wolbachia_titer', groupby='leiden', save='_wolbachia_titer_by_cluster.pdf')
    
    # If both methods are used, compare them
    if method == 'both':
        # Run UMAP on BBKNN object for comparison
        sc.tl.umap(combined_bbknn)
        sc.tl.leiden(combined_bbknn, resolution=0.8)
        
        # Generate comparison plots
        sc.pl.umap(combined_bbknn, color=batch_key, save='_bbknn_batch.pdf')
        sc.pl.umap(combined_bbknn, color='leiden', save='_bbknn_leiden.pdf')
        sc.tl.rank_genes_groups(combined_bbknn, 'leiden', method='wilcoxon', n_genes=50)
        sc.pl.rank_genes_groups(combined_bbknn, n_genes=50, save='_bbknn_rank_genes.pdf')
        #Output .txt file of rank_genes_groups
        rank_genes = combined_bbknn.uns['rank_genes_groups']
        
        
        if 'wolbachia_titer' in combined_bbknn.obs.columns:
            sc.pl.umap(combined_bbknn, color='wolbachia_titer', save='_bbknn_wolbachia_titer.pdf')
            sc.pl.umap(combined_bbknn, color='log1p_wolbachia_titer', save='_bbknn_log1p_wolbachia_titer.pdf')
        
        # Save the BBKNN object as well
        combined_bbknn.write(output_path.replace('.h5ad', '_bbknn.h5ad'))
        
        print("Saved both Harmony and BBKNN integrated objects for comparison")
    
    print("Integration complete!")
    return combined
def main():
    parser = argparse.ArgumentParser(description='Integrate multiple h5ad files with batch correction')
    
    # Required arguments
    parser.add_argument('--input_dir', type=str, required=True, 
                        help='Directory containing h5ad files to integrate')
    parser.add_argument('--output_file', type=str, required=True, 
                        help='Path to save the integrated h5ad file')
    
    # Optional arguments
    parser.add_argument('--batch_key', type=str, default='batch',
                        help='Name of column to use for batch correction')
    parser.add_argument('--min_cells', type=int, default=3,
                        help='Minimum number of cells expressing a gene')
    parser.add_argument('--min_genes', type=int, default=200,
                        help='Minimum number of genes expressed in a cell')
    parser.add_argument('--n_pcs', type=int, default=30,
                        help='Number of principal components to use')
    parser.add_argument('--n_neighbors', type=int, default=15,
                        help='Number of neighbors for neighborhood graph')
    parser.add_argument('--method', type=str, default='bbknn', choices=['bbknn', 'harmony', 'both'],
                        help='Batch correction method to use')
    parser.add_argument('--calculate_titer', action='store_true', 
                        help='Calculate Wolbachia titer for each cell before integration')
    
    # New arguments for Wolbachia titer calculation
    parser.add_argument('--min_total_reads', type=int, default=100,
                        help='Minimum total reads required to consider a cell for titer calculation')
    parser.add_argument('--min_wolb_reads', type=int, default=3,
                        help='Minimum Wolbachia rRNA reads required to consider a cell infected')
    parser.add_argument('--depth_normalize', action='store_true', default=True,
                        help='Normalize titer by sequencing depth')
    parser.add_argument('--no_depth_normalize', action='store_false', dest='depth_normalize',
                        help='Do not normalize titer by sequencing depth')
    
    args = parser.parse_args()
    
    # Run the integration
    integrated_adata = integrate_h5ad_files(
        directory_path=args.input_dir,
        output_path=args.output_file,
        batch_key=args.batch_key,
        min_cells=args.min_cells,
        min_genes=args.min_genes,
        n_pcs=args.n_pcs,
        n_neighbors=args.n_neighbors,
        method=args.method,
        calculate_titer=args.calculate_titer,
        min_total_reads=args.min_total_reads,
        min_wolb_reads=args.min_wolb_reads,
        depth_normalize=args.depth_normalize
    )
    
    # Print summary
    if integrated_adata is not None:
        print("Summary of integrated data:")
        print(f"Number of cells: {integrated_adata.n_obs}")
        print(f"Number of genes: {integrated_adata.n_vars}")
        print(f"Number of batches: {integrated_adata.obs[args.batch_key].nunique()}")
        print(f"Number of clusters: {integrated_adata.obs['leiden'].nunique()}")
        
        if 'wolbachia_titer' in integrated_adata.obs.columns:
            # Calculate percentage of infected cells (now using > 0 instead of ~np.isnan)
            n_infected = np.sum(integrated_adata.obs['wolbachia_titer'] > 0)
            print(f"Number of cells with Wolbachia: {n_infected} ({n_infected/integrated_adata.n_obs*100:.2f}%)")
            
            # Calculate average titer among infected cells only
            infected_mask = integrated_adata.obs['wolbachia_titer'] > 0
            if np.sum(infected_mask) > 0:
                mean_titer = np.mean(integrated_adata.obs['wolbachia_titer'][infected_mask])
                median_titer = np.median(integrated_adata.obs['wolbachia_titer'][infected_mask])
                print(f"Average Wolbachia titer (infected cells only): mean={mean_titer:.4f}, median={median_titer:.4f}")
            else:
                print("No infected cells detected.")
            
            # Calculate titer by Sample instead of batch
            if 'Sample' in integrated_adata.obs.columns:
                print("\nWolbachia statistics by Sample:")
                for sample in integrated_adata.obs['Sample'].unique():
                    sample_cells = integrated_adata[integrated_adata.obs['Sample'] == sample]
                    n_sample_infected = np.sum(sample_cells.obs['wolbachia_titer'] > 0)
                    n_sample_total = sample_cells.n_obs
                    if n_sample_infected > 0:
                        sample_infected_mask = sample_cells.obs['wolbachia_titer'] > 0
                        mean_sample_titer = np.mean(sample_cells.obs['wolbachia_titer'][sample_infected_mask])
                        print(f"  {sample}: {n_sample_infected}/{n_sample_total} cells infected ({n_sample_infected/n_sample_total*100:.2f}%), mean titer={mean_sample_titer:.4f}")
                    else:
                        print(f"  {sample}: 0/{n_sample_total} cells infected (0.00%)")
            else:
                # Fallback to batch if Sample doesn't exist
                print("\nWolbachia statistics by batch:")
                for batch in integrated_adata.obs[args.batch_key].unique():
                    batch_cells = integrated_adata[integrated_adata.obs[args.batch_key] == batch]
                    n_batch_infected = np.sum(batch_cells.obs['wolbachia_titer'] > 0)
                    n_batch_total = batch_cells.n_obs
                    if n_batch_infected > 0:
                        batch_infected_mask = batch_cells.obs['wolbachia_titer'] > 0
                        mean_batch_titer = np.mean(batch_cells.obs['wolbachia_titer'][batch_infected_mask])
                        print(f"  {batch}: {n_batch_infected}/{n_batch_total} cells infected ({n_batch_infected/n_batch_total*100:.2f}%), mean titer={mean_batch_titer:.4f}")
                    else:
                        print(f"  {batch}: 0/{n_batch_total} cells infected (0.00%)")
            
            # Print threshold information
            print(f"\nTiter calculation parameters:")
            print(f"  Minimum total reads threshold: {args.min_total_reads}")
            print(f"  Minimum Wolbachia reads threshold: {args.min_wolb_reads}")
            print(f"  Depth normalization: {'Yes' if args.depth_normalize else 'No'}")
            
            # If we have wolbachia_titer_raw, calculate how many cells were filtered by thresholds
            if 'wolbachia_titer_raw' in integrated_adata.obs.columns and 'wolbachia_detected' in integrated_adata.obs.columns:
                # Count cells that have any signal in raw counts but were filtered out
                n_any_signal = np.sum(integrated_adata.obs['wMel_raw_counts'] > 0)
                n_passed_thresholds = np.sum(integrated_adata.obs['wolbachia_detected'])
                n_filtered = n_any_signal - n_passed_thresholds
                if n_any_signal > 0:
                    print(f"  Cells filtered by thresholds: {n_filtered} ({n_filtered/n_any_signal*100:.2f}% of cells with any Wolbachia signal)")

if __name__ == "__main__":
    main()