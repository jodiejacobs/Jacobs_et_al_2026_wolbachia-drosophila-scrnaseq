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

# rRNA gene dictionaries with lengths
wMel_rRNA={
    "GQX67_00940": 2772,
    "GQX67_00945": 107,
    "GQX67_05945": 1505
}

all_rRNA={
    "FBgn0013686": 1324, # Dmel mtrRNA
    "FBgn0013688": 786,  # Dmel mtrRNA
    "FBgn0053353": 135,
    "FBgn0053354": 135,
    "FBgn0053355": 135,
    "FBgn0053356": 143,
    "FBgn0053357": 135,
    "FBgn0053358": 135,
    "FBgn0053359": 135,
    "FBgn0053360": 135,
    "FBgn0053361": 135,
    "FBgn005336": 135,
    "FBgn0053363": 144,
    "FBgn0053364": 135,
    "FBgn005336": 135,
    "FBgn0053366": 135,
    "FBgn0053367": 135,
    "FBgn0053368": 135,
    "FBgn0053369": 135,
    "FBgn0053370": 135,
    "FBgn0053371": 137,
    "FBgn0053372": 135,
    "FBgn0053373": 135,
    "FBgn0053374": 135,
    "FBgn0053375": 135,
    "FBgn0053376": 135,
    "FBgn0053377": 135,
    "FBgn0053378": 135,
    "FBgn0053379": 135,
    "FBgn0053380": 135,
    "FBgn0053381": 135,
    "FBgn0053382": 135,
    "FBgn0053383": 135,
    "FBgn0053384": 135,
    "FBgn0053385": 135,
    "FBgn0053386": 135,
    "FBgn0053387": 135,
    "FBgn0053388": 135,
    "FBgn0053389": 135,
    "FBgn0053390": 135,
    "FBgn0053391": 134,
    "FBgn0053392": 135,
    "FBgn0053393": 135,
    "FBgn0053394": 135,
    "FBgn0053395": 135,
    "FBgn0053396": 135,
    "FBgn0053397": 135,
    "FBgn0053398": 135,
    "FBgn0053399": 135,
    "FBgn0053400": 135,
    "FBgn0053401": 135,
    "FBgn0053402": 135,
    "FBgn0053403": 135,
    "FBgn0053404": 135,
    "FBgn0053405": 134,
    "FBgn0053406": 135,
    "FBgn0053407": 135,
    "FBgn0053408": 135,
    "FBgn0053409": 135,
    "FBgn0053410": 135,
    "FBgn0053411": 135,
    "FBgn0053412": 135,
    "FBgn0053413": 135,
    "FBgn0053414": 135,
    "FBgn0053415": 135,
    "FBgn0053416": 143,
    "FBgn0053417": 135,
    "FBgn0053418": 135,
    "FBgn0053419": 135,
    "FBgn0053420": 135,
    "FBgn0053421": 135,
    "FBgn0053422": 135,
    "FBgn0053423": 135,
    "FBgn0053424": 135,
    "FBgn0053425": 135,
    "FBgn0053426": 135,
    "FBgn0053427": 135,
    "FBgn0053428": 135,
    "FBgn0053429": 135,
    "FBgn0053430": 135,
    "FBgn0053431": 135,
    "FBgn0053432": 135,
    "FBgn0053433": 135,
    "FBgn0053434": 135,
    "FBgn0053435": 135,
    "FBgn0053436": 135,
    "FBgn0053437": 135,
    "FBgn0053438": 135,
    "FBgn0053439": 135,
    "FBgn0053440": 135,
    "FBgn0053441": 135,
    "FBgn0053442": 135,
    "FBgn0053443": 135,
    "FBgn0053444": 135,
    "FBgn0053445": 135,
    "FBgn0053446": 135,
    "FBgn0053447": 135,
    "FBgn0053448": 135,
    "FBgn0053449": 135,
    "FBgn0053450": 135,
    "FBgn0053451": 135,
    "FBgn0053452": 135,
    "FBgn0085753": 6005,
    "FBgn0085765": 30,
    "FBgn0085771": 1258,
    "FBgn0085802": 1995,
    "FBgn0085813": 1975,
    "FBgn0085819": 895,
    "FBgn0250731": 123,
    "FBgn0267496": 30,
    "FBgn0267497": 2715,
    "FBgn0267498": 1995,
    "FBgn0267499": 123,
    "FBgn0267500": 30,
    "FBgn0267501": 1995,
    "FBgn0267502": 123,
    "FBgn0267503": 30,
    "FBgn0267504": 3970,
    "FBgn0267505": 6863,
    "FBgn0267506": 4191,
    "FBgn0267507": 8118,
    "FBgn0267508": 821,
    "FBgn0267509": 123,
    "FBgn0267510": 30,
    "FBgn0267511": 2800,
    "FBgn0267512": 123,
    "FBgn0267513": 255,
    "FBgn0267514": 123,
    "FBgn0267515": 704,
    "FBgn0267516": 5725,
    "FBgn0267517": 123,
    "FBgn0267518": 30,
    "FBgn0267519": 2689,
    "FBgn0267520": 357,
    "FBgn0267521": 1934,
    "FBgn0267522": 2004,
    "FBgn0267523": 123,
    "FBgn0267524": 30
}

# Set plotting settings
sc.settings.set_figure_params(dpi=100, frameon=False)
sc.settings.verbosity = 1

# Import bacteria gene list
def import_bacteria_genes(filepath):
    with open(filepath, 'r') as f:
        bacteria_genes = {line.strip() for line in f if line.strip()} # Makes a set of bacteria genes
    return bacteria_genes

bacteria_genes = import_bacteria_genes('/private/groups/russelllab/jodie/scRNAseq/reference/gene_list_bacteria.txt')

def calculate_wolbachia_titer(adata):
    '''
    Calculate Wolbachia titer for each cell in the AnnData object.
    The titer is calculated as the ratio of wMel rRNA counts per length to Dmel rRNA counts per length.
    '''
    print("Calculating Wolbachia titer...")
    
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
                
        all_genes_present = []
        for gene_id in all_rRNA.keys():
            # Check if this gene ID is in the gene_ids column
            if gene_id in adata.var['gene_ids'].values:
                all_genes_present.append(gene_id)
        
        # Create masks for the genes
        wMel_mask = [gene_id_map.get(idx) in wMel_genes_present for idx in adata.var.index]
        all_mask = [gene_id_map.get(idx) in all_genes_present for idx in adata.var.index]
    else:
        # Use original approach with var_names
        wMel_genes_present = [gene for gene in wMel_rRNA.keys() if gene in adata.var_names]
        all_genes_present = [gene for gene in all_rRNA.keys() if gene in adata.var_names]
        
        # Create masks based on var_names
        wMel_mask = [gene in wMel_genes_present for gene in adata.var_names]
        all_mask = [gene in all_genes_present for gene in adata.var_names]
    
    if not wMel_genes_present:
        print("Warning: No wMel rRNA genes found in the dataset!")
        adata.obs['wolbachia_titer'] = np.nan
        return adata
        
    if not all_genes_present:
        print("Warning: No Dmel rRNA genes found in the dataset!")
        adata.obs['wolbachia_titer'] = np.nan
        return adata
    
    print(f"Found {len(wMel_genes_present)} wMel rRNA genes and {len(all_genes_present)} Dmel rRNA genes")
    
    # Get gene indices from the masks
    wMel_indices = np.where(wMel_mask)[0]
    all_indices = np.where(all_mask)[0]
    
    # Convert sparse matrix to dense if necessary
    is_sparse = scipy.sparse.issparse(adata.X)
    
    # Initialize arrays properly - FIXED BUG HERE
    wMel_counts_per_length = np.zeros((adata.n_obs, len(wMel_indices)))
    all_counts_per_length = np.zeros((adata.n_obs, len(all_indices)))
    
    # Calculate counts per length for wMel genes
    for i, idx in enumerate(wMel_indices):
        gene_idx = adata.var.index[idx]
        gene_id = gene_id_map.get(gene_idx, gene_idx) if 'gene_ids' in adata.var.columns else gene_idx

        if is_sparse:
            counts = adata.X[:, idx].toarray().flatten()
        else:
            counts = adata.X[:, idx]
            
        wMel_counts_per_length[:, i] = counts
    
    # Calculate counts per length for Dmel genes
    for i, idx in enumerate(all_indices):
        gene_idx = adata.var.index[idx]
        gene_id = gene_id_map.get(gene_idx, gene_idx) if 'gene_ids' in adata.var.columns else gene_idx

        if is_sparse:
            counts = adata.X[:, idx].toarray().flatten()
        else:
            counts = adata.X[:, idx]
            
        all_counts_per_length[:, i] = counts 
    
    # Calculate the mean counts per length for each organism and cell
    wMel_mean_per_cell = np.mean(wMel_counts_per_length, axis=1)
    all_mean_per_cell = np.mean(all_counts_per_length, axis=1)
    
    # Calculate the titer (ratio of wMel to Dmel rRNA counts per length)
    with np.errstate(divide='ignore', invalid='ignore'):
        titer = np.where(all_mean_per_cell > 0, 
                         wMel_mean_per_cell / all_mean_per_cell, 
                         np.nan)
    
    # Add the titer to the AnnData object
    adata.obs['wolbachia_titer'] = titer
    adata.obs['log1p_wolbachia_titer'] = np.log1p(titer)
    
    # Count cells with Wolbachia
    n_infected = np.sum(titer > 0)
    print(f"Detected Wolbachia in {n_infected} out of {adata.n_obs} cells ({n_infected/adata.n_obs*100:.2f}%)")
    
    return adata

def extract_sample_type(filename, sample_type_pattern=None):
    """
    Extract sample type from filename based on a pattern or delimiter
    
    Parameters:
    -----------
    filename : str
        Filename to extract sample type from
    sample_type_pattern : str or None
        Regex pattern to extract sample type. If None, assumes filename format: "SampleType_OtherInfo.h5ad"
        If pattern has multiple capture groups, they'll be combined with '_' delimiter.
        
    Returns:
    --------
    sample_type : str
        Extracted sample type
    """
    basename = os.path.basename(filename)
    name_without_ext = os.path.splitext(basename)[0]
    
    if sample_type_pattern:
        # Use provided regex pattern to extract sample type
        match = re.search(sample_type_pattern, name_without_ext)
        if match:
            # Check if we have multiple capture groups
            if match.lastindex and match.lastindex > 1:
                # Combine all captured groups using '_' delimiter
                parts = [match.group(i) for i in range(1, match.lastindex + 1)]
                return '_'.join(parts)
            else:
                return match.group(1)
        else:
            # If pattern doesn't match, use the whole name as the sample type
            print(f"Warning: Pattern didn't match for file {basename}. Using whole name as sample type.")
            return name_without_ext
    else:
        # Default behavior: assume SampleType_OtherInfo format with underscore delimiter
        parts = name_without_ext.split('_')
        if len(parts) > 0:
            return parts[0]
        else:
            return name_without_ext

def integrate_h5ad_files_by_sample_type(directory_path, output_dir, sample_type_pattern=None, batch_key='batch', 
                                       min_cells=3, min_genes=200, n_pcs=30, n_neighbors=15,
                                       method='bbknn', calculate_titer=True, prefix=None):
    """
    Group h5ad files by sample type and integrate each group separately
    
    Parameters:
    -----------
    directory_path : str
        Path to directory containing h5ad files
    output_dir : str
        Directory to save integrated h5ad files
    sample_type_pattern : str or None
        Regex pattern to extract sample type from filenames. If None, uses SampleType_OtherInfo assumption
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
    calculate_titer : bool
        Whether to calculate Wolbachia titer before integration
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Create figures directory
    figure_dir = os.path.join(output_dir, "figures")
    os.makedirs(figure_dir, exist_ok=True)
    sc.settings.figdir = figure_dir
    
    # Get all h5ad files in the directory
    h5ad_files = glob(os.path.join(directory_path, "*.h5ad"))
    
    if len(h5ad_files) == 0:
        print(f"No h5ad files found in {directory_path}")
        return
    
    print(f"Found {len(h5ad_files)} h5ad files")
    
    # Group files by sample type
    sample_type_to_files = {}
    
    for file_path in h5ad_files:
        sample_type = extract_sample_type(file_path, sample_type_pattern)
        if sample_type not in sample_type_to_files:
            sample_type_to_files[sample_type] = []
        sample_type_to_files[sample_type].append(file_path)
    
    print(f"Grouped files into {len(sample_type_to_files)} sample types:")
    for sample_type, files in sample_type_to_files.items():
        print(f"  {sample_type}: {len(files)} files")
    
    # Process each sample type
    for sample_type, files in sample_type_to_files.items():
        if len(files) < 2:
            print(f"Skipping sample type '{sample_type}' as it only has {len(files)} file(s)")
            continue
            
        print(f"\nProcessing sample type: {sample_type} ({len(files)} files)")
        output_path = os.path.join(output_dir, f"{prefix}_integrated.h5ad")
        
        # List to store individual datasets for this sample type
        adatas = []
        
        # Load each dataset and add batch information
        for i, file_path in enumerate(files):
            file_name = os.path.basename(file_path)
            batch_id = os.path.splitext(file_name)[0]  # Use filename without extension as batch ID
            print(f"Processing file {i+1}/{len(files)}: {file_name}")
            
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
            adata.obs[batch_key] = batch_id
            
            # Calculate Wolbachia titer if requested
            if calculate_titer:
                adata = calculate_wolbachia_titer(adata)
            
            # Filter out bacterial genes
            bacteria_mask = adata.var_names.isin(bacteria_genes)
            adata = adata[:, ~bacteria_mask]

            adatas.append(adata)
            
        # Concatenate all datasets for this sample type
        print(f"Concatenating {len(adatas)} datasets for sample type {sample_type}...")
        try:
            # Try using anndata.concat as recommended by the FutureWarning
            combined = ad.concat(adatas, join='outer', merge='same', label=batch_key, index_unique='-')
        except:
            # Fall back to concatenate method if concat fails
            combined = adatas[0].concatenate(adatas[1:], join='outer', batch_key=batch_key)
        
        print(f"Combined data shape for {sample_type}: {combined.shape}")
        
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
        sc.pl.umap(combined_unintegrated, color=batch_key, save=f'_{prefix}_before_batch_correction.pdf')

        bbknn.bbknn(combined, batch_key=batch_key, n_pcs=n_pcs, neighbors_within_batch=5)
                
        # Run UMAP and clustering
        print("Running UMAP and Leiden clustering...")
        sc.tl.umap(combined)
        sc.tl.leiden(combined, resolution=0.8)
        
        # Save the integrated object
        print(f"Saving integrated object for {prefix} to {output_path}")
        combined.write(output_path)
        
        # Generate diagnostic plots
        print("Generating diagnostic plots...")
        sc.pl.umap(combined, color=batch_key, save=f'_{prefix}_bbknn.pdf')
        sc.pl.umap(combined, color='leiden', save=f'_{prefix}_bbknn_leiden.pdf')
        
        # If wolbachia_titer exists, plot it too
        if 'wolbachia_titer' in combined.obs.columns:
            sc.pl.umap(combined, color='wolbachia_titer', save=f'_{prefix}_bbknn_titer.pdf')
            sc.pl.umap(combined, color='log1p_wolbachia_titer', save=f'_{prefix}_log1p_bbknn_titer.pdf')
            
            # Create a violin plot of titer by batch
            sc.pl.violin(combined, 'wolbachia_titer', groupby=batch_key, save=f'_{prefix}_wolbachia_titer_by_rep.pdf')
            
            # Create a violin plot of titer by cluster
            sc.pl.violin(combined, 'wolbachia_titer', groupby='leiden', save=f'_{prefix}_wolbachia_titer_by_cluster.pdf')
        
        print(f"Integration complete for sample type {prefix}!")
        
        # Print summary for this sample type
        print(f"Summary of integrated data for {prefix}:")
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
    
    print("\nAll sample types processed!")

def main():
    parser = argparse.ArgumentParser(description='Integrate h5ad files by sample type with batch correction')
    
    # Required arguments
    parser.add_argument('--input_dir', type=str, required=True, 
                        help='Directory containing h5ad files to integrate')
    parser.add_argument('--output_dir', type=str, required=True, 
                        help='Directory to save the integrated h5ad files')
    
    # Optional arguments
    parser.add_argument('--sample_type_pattern', type=str, default=None,
                        help='Regex pattern to extract sample type from filename (e.g., "^([^_]+)_")')
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
    parser.add_argument('--prefix', type=str, default=None,
                        help='Prefix to use for output filename')

    args = parser.parse_args()
    
    # Run the integration by sample type
    integrate_h5ad_files_by_sample_type(
        directory_path=args.input_dir,
        output_dir=args.output_dir,
        sample_type_pattern=args.sample_type_pattern,
        batch_key=args.batch_key,
        min_cells=args.min_cells,
        min_genes=args.min_genes,
        n_pcs=args.n_pcs,
        n_neighbors=args.n_neighbors,
        calculate_titer=args.calculate_titer,
        prefix=args.prefix
    )

if __name__ == "__main__":
    main()