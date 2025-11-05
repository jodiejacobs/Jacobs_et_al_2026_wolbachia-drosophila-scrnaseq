'''
Filtering script for kallisto bustools aligned single cell data. 
    Input: raw data object 
    Output: QC plots and filtered adata object
'''

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import scrublet as scr
from scipy.io import mmread
from scipy.interpolate import interpn
from sklearn.decomposition import TruncatedSVD
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.sparse import csr_matrix
import os
from pybiomart import Dataset
import sys
import gzip
import shutil
import gc #Garbage Collection to remove old data 
import argparse

# Set up the Argument Parser
parser = argparse.ArgumentParser(description="Kallisto Bustools QC Filter")
parser.add_argument("--input", type=str, required=False, help="Input h5ad file")
parser.add_argument("--output", type=str, required=False, help="Output h5ad file")

args = parser.parse_args()
input = args.input
output = args.output    

# Get the basename for the h5ad input (remove the .h5ad):
sample_name = os.path.basename(input).replace(".h5ad", "")

# Get the output directory
output_dir = os.path.dirname(output)
fig_dir = f"{output_dir}/{sample_name}"

# Make output directories
os.makedirs(output_dir, exist_ok=True)
os.makedirs(fig_dir, exist_ok=True)

# Set figdir for scanpy:
sc.settings.autosave = True
sc.settings.figdir = fig_dir

# Set global matplotlib parameters
matplotlib.rcParams.update({
    'font.size': 6,
    'figure.figsize': [2, 2],
    'axes.titlesize': 6,
    'axes.labelsize': 6,
    'xtick.labelsize': 6,
    'ytick.labelsize': 6,
    'legend.fontsize': 6,
    'figure.titlesize': 6
})

# Set scanpy figure parameters
sc.settings.set_figure_params(
    dpi=300,  # High resolution
    dpi_save=300,
    figsize=(2, 2),  # 2x2 inches
    fontsize=6
)

dataset = Dataset(name='dmelanogaster_gene_ensembl', host='http://www.ensembl.org')

# Function to get mitochondrial genes from biomart, properly filtered
mito_genes = dataset.query(
    attributes=['ensembl_gene_id', 'external_gene_name'],
    filters={'chromosome_name': ['mitochondrion_genome']}
)
mito_genes = mito_genes['Gene stable ID'].values.flatten()

# density display for PCA plot

def identify_doublets(adata, fig_dir):
    '''
    Doublet identification with scrublet 
    '''
    
    print("Starting scrublet doublet detection:")
    
    # Initialize Scrublet with adata object
    scrub = scr.Scrublet(adata.X, expected_doublet_rate=0.1) # requires an estimated rate as a prior 
    # Run the doublet detection
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)
    # scrub.call_doublets(threshold=0.25)
    scrub.call_doublets()
    
    print("Plotting scrublet histogram")
    # Plot doublet score histogram
    scrub.plot_histogram()
    plt.savefig(f"{fig_dir}/doublet_histogram.pdf", bbox_inches='tight', pad_inches=0.1)
    plt.close()

    # Plot UMAP:
    print('Running UMAP...')
    scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
    scrub.plot_embedding('UMAP', order_points=True)
    plt.savefig(f"{fig_dir}/doublet_umap.pdf", bbox_inches='tight', pad_inches=0.1)
    plt.close()
    
    print("Saving scrublet data to adata")
    
    # Add scrublet results back to adata object FIRST
    adata.obs['doublet_score'] = doublet_scores
    adata.obs['predicted_doublet'] = predicted_doublets
    
    # THEN create the categorical version
    adata.obs['predicted_doublet_cat'] = adata.obs['predicted_doublet'].astype(str).astype('category')
    
    return adata

def calculate_qc_metrics(adata, sample_name="", stage=""):
    """
    Memory-optimized QC metrics calculation - avoids converting entire matrix to dense
    """
    import numpy as np
    import pandas as pd
    from scipy import sparse
    
    # Basic cell and gene counts
    n_cells = adata.n_obs
    n_genes_total = adata.n_vars
    
    # Work with sparse matrix directly - don't convert to dense!
    X = adata.X
    
    # Genes per cell (number of detected genes per cell) - sparse friendly
    if sparse.issparse(X):
        genes_per_cell = np.asarray((X > 0).sum(axis=1)).flatten()
        umis_per_cell = np.asarray(X.sum(axis=1)).flatten()
    else:
        genes_per_cell = np.sum(X > 0, axis=1)
        umis_per_cell = np.sum(X, axis=1)
    
    # Transcriptome coverage per cell
    transcriptome_coverage_per_cell = genes_per_cell / n_genes_total
    
    # Doublet rate (if available)
    doublet_rate = None
    n_doublets = None
    if 'predicted_doublet' in adata.obs.columns:
        doublet_rate = adata.obs['predicted_doublet'].sum() / len(adata.obs)
        n_doublets = adata.obs['predicted_doublet'].sum()
    
    # Compile metrics
    metrics = {
        'sample_name': sample_name,
        'stage': stage,
        'n_cells': n_cells,
        'n_genes_total': n_genes_total,
        
        # Genes per cell statistics
        'genes_per_cell_mean': float(np.mean(genes_per_cell)),
        'genes_per_cell_median': float(np.median(genes_per_cell)),
        'genes_per_cell_std': float(np.std(genes_per_cell)),
        'genes_per_cell_min': int(np.min(genes_per_cell)),
        'genes_per_cell_max': int(np.max(genes_per_cell)),
        
        # UMIs per cell statistics
        'umis_per_cell_mean': float(np.mean(umis_per_cell)),
        'umis_per_cell_median': float(np.median(umis_per_cell)),
        'umis_per_cell_std': float(np.std(umis_per_cell)),
        'umis_per_cell_min': float(np.min(umis_per_cell)),
        'umis_per_cell_max': float(np.max(umis_per_cell)),
        
        # Transcriptome coverage statistics
        'transcriptome_coverage_mean': float(np.mean(transcriptome_coverage_per_cell)),
        'transcriptome_coverage_median': float(np.median(transcriptome_coverage_per_cell)),
        'transcriptome_coverage_std': float(np.std(transcriptome_coverage_per_cell)),
        'transcriptome_coverage_min': float(np.min(transcriptome_coverage_per_cell)),
        'transcriptome_coverage_max': float(np.max(transcriptome_coverage_per_cell)),
        
        # Doublet information
        'doublet_rate': doublet_rate,
        'n_doublets': n_doublets,
        
        # Additional useful metrics
        'mitochondrial_percent_mean': float(adata.obs['percent_mito'].mean()) if 'percent_mito' in adata.obs.columns else None,
        'mitochondrial_percent_median': float(adata.obs['percent_mito'].median()) if 'percent_mito' in adata.obs.columns else None,
    }
    
    # Clean up temporary arrays
    del genes_per_cell, umis_per_cell, transcriptome_coverage_per_cell
    gc.collect()
    
    return metrics

def print_qc_summary(metrics):
    """
    Print a formatted summary of QC metrics.
    
    Parameters:
    -----------
    metrics : dict
        Dictionary of QC metrics from calculate_qc_metrics()
    """
    print(f"\n{'='*50}")
    print(f"QC METRICS SUMMARY")
    if metrics['sample_name']:
        print(f"Sample: {metrics['sample_name']}")
    if metrics['stage']:
        print(f"Stage: {metrics['stage']}")
    print(f"{'='*50}")
    
    print(f"Dataset Overview:")
    print(f"  Total cells: {metrics['n_cells']:,}")
    print(f"  Total genes: {metrics['n_genes_total']:,}")
    
    print(f"\nGenes per cell:")
    print(f"  Mean: {metrics['genes_per_cell_mean']:.1f}")
    print(f"  Median: {metrics['genes_per_cell_median']:.1f}")
    print(f"  Range: {metrics['genes_per_cell_min']:.0f} - {metrics['genes_per_cell_max']:.0f}")
    
    print(f"\nUMIs per cell:")
    print(f"  Mean: {metrics['umis_per_cell_mean']:.1f}")
    print(f"  Median: {metrics['umis_per_cell_median']:.1f}")
    print(f"  Range: {metrics['umis_per_cell_min']:.0f} - {metrics['umis_per_cell_max']:.0f}")
    
    print(f"\nTranscriptome coverage per cell:")
    print(f"  Mean: {metrics['transcriptome_coverage_mean']:.3f} ({metrics['transcriptome_coverage_mean']*100:.1f}%)")
    print(f"  Median: {metrics['transcriptome_coverage_median']:.3f} ({metrics['transcriptome_coverage_median']*100:.1f}%)")
    print(f"  Range: {metrics['transcriptome_coverage_min']:.3f} - {metrics['transcriptome_coverage_max']:.3f}")
    
    if metrics['doublet_rate'] is not None:
        print(f"\nDoublet information:")
        print(f"  Doublet rate: {metrics['doublet_rate']:.3f} ({metrics['doublet_rate']*100:.1f}%)")
        print(f"  Number of doublets: {metrics['n_doublets']}")
    
    if metrics['mitochondrial_percent_mean'] is not None:
        print(f"\nMitochondrial gene expression:")
        print(f"  Mean: {metrics['mitochondrial_percent_mean']:.3f} ({metrics['mitochondrial_percent_mean']*100:.1f}%)")
        print(f"  Median: {metrics['mitochondrial_percent_median']:.3f} ({metrics['mitochondrial_percent_median']*100:.1f}%)")
    
    print(f"{'='*50}\n")

def save_metrics_to_csv(metrics_list, output_path):
    """
    Save a list of metrics dictionaries to CSV file.
    
    Parameters:
    -----------
    metrics_list : list
        List of metrics dictionaries from calculate_qc_metrics()
    output_path : str
        Path to save the CSV file
    """
    import pandas as pd
    
    df = pd.DataFrame(metrics_list)
    df.to_csv(output_path, index=False)
    print(f"Metrics saved to {output_path}")
    

def density_scatter( x , y, ax = None, sort = True, bins = 20, **kwargs )   :
    """
    Scatter plot colored by 2d histogram
    """
    if ax is None :
        fig , ax = plt.subplots()
    data , x_e, y_e = np.histogram2d( x, y, bins = bins)
    z = interpn( ( 0.5*(x_e[1:] + x_e[:-1]) , 0.5*(y_e[1:]+y_e[:-1]) ) , data , np.vstack([x,y]).T , method = "splinef2d", bounds_error = False )

    # Sort the points by density, so that the densest points are plotted last
    if sort :
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]

    scatter_obj = ax.scatter( x, y, c=z, **kwargs )
    return scatter_obj


### Basic QC, Check the cells in 2D space
def qc_plots(adata, fig_dir):
    """
    Perform basic QC and visualization for a count matrix and its PCA projection.
    Saves plots to disk.
    """
    # Perform SVD
    tsvd = TruncatedSVD(n_components=2)
    tsvd.fit(adata.X)
    X = tsvd.transform(adata.X)

    # Plot the cells in the 2D PCA projection
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.scatter(X[:,0], X[:,1], alpha=0.5, c="green")
    plt.axis('off')
    plt.savefig(f"{fig_dir}/pca_plot.pdf", bbox_inches='tight', pad_inches=0.1)
    plt.close(fig)

    # density display for PCA plot
    fig, ax = plt.subplots(figsize=(7,7))
    x = X[:,0]
    y = X[:,1]
    scatter_obj = density_scatter(x, y, ax=ax, cmap="Greens")
    fig.colorbar(scatter_obj, ax=ax)
    plt.axis('off')
    # plt.show()
    plt.savefig(f"{fig_dir}/density_scatter.pdf", bbox_inches='tight', pad_inches=0.1)
    plt.close(fig)

    # Test for library saturation:
    # Create a plot showing genes detected as a function of UMI counts.
    fig, ax = plt.subplots(figsize=(10, 7))
    x = np.asarray(adata.X.sum(axis=1))[:,0]
    y = np.asarray(np.sum(adata.X>0, axis=1))[:,0]
    ax.scatter(x, y, color="green", alpha=0.25)
    ax.set_xlabel("UMI Counts")
    ax.set_ylabel("Genes Detected")
    ax.set_xscale('log')
    ax.set_yscale('log')
    # ax.set_xlim((0.5, 4500))
    # ax.set_ylim((0.5,2000))
    plt.savefig(f"{fig_dir}/genes_vs_umi.pdf", bbox_inches='tight', pad_inches=0.1)
    plt.close(fig)

    # Knee plot for cutoffs
    #Threshold cells according to knee plot { run: "auto", vertical-output: true }
    cutoff =  100
    knee = np.sort((np.array(adata.X.sum(axis=1))).flatten())[::-1]
    cell_set = np.arange(len(knee))
    num_cells = cell_set[knee > cutoff][::-1][0]
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.loglog(knee, cell_set, linewidth=5, color="g")
    ax.axvline(x=cutoff, linewidth=3, color="k")
    ax.axhline(y=num_cells, linewidth=3, color="k")
    ax.set_xlabel("UMI Counts")
    ax.set_ylabel("Set of Barcodes")
    plt.grid(True, which="both")
    plt.savefig(f"{fig_dir}/knee_plot.pdf", bbox_inches='tight', pad_inches=0.1)
    plt.close(fig)

    print(f"{fig_dir}\n{num_cells:,.0f} cells passed the {cutoff} UMI threshold")

    # Filter the empty droplets
    print("Raw adata matrix:")
    print(adata)
    
    # Filter the cells according to the knee plot threshold
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_cells(adata, min_counts=knee[num_cells])
    print(f'After knee plot thresholding:')
    print(adata)

    # Annotate mitochondrial genes
    print("Annotating mitochondrial genes")
    # Get mitochondrial genes that are actually present in the data
    valid_mito_genes = [gene for gene in mito_genes if gene in adata.var_names]

    # Calculate percent mitochondrial
    if len(valid_mito_genes) > 0:
        print(f"Using {len(valid_mito_genes)} mitochondrial genes for filtering:")
        for gene in valid_mito_genes:
            print(f"  {gene}")
            
        # Calculate mitochondrial percentage
        mito_counts = adata[:, valid_mito_genes].X.sum(axis=1)
        total_counts = adata.X.sum(axis=1)
        
        # Handle sparse matrices properly
        if hasattr(mito_counts, 'A1'):
            mito_counts = mito_counts.A1
        if hasattr(total_counts, 'A1'):
            total_counts = total_counts.A1
            
        adata.obs['percent_mito'] = mito_counts / total_counts
    else:
        print("No mitochondrial genes found in dataset, setting percent_mito to 0")
        adata.obs['percent_mito'] = 0.0

    # add the total counts per cell as observations-annotation to adata
    adata.obs['n_counts'] = adata.X.sum(axis=1)
    if hasattr(adata.obs['n_counts'], 'A1'):
        adata.obs['n_counts'] = adata.obs['n_counts'].A1

    # Create the scatter plot and properly close it
    fig, ax = plt.subplots(figsize=(10, 7))
    sc.pl.scatter(adata, x='n_counts', y='percent_mito', ax=ax, show=False)
    plt.savefig(f"{fig_dir}/mito_scatter.pdf", bbox_inches='tight', pad_inches=0.1)
    plt.close(fig)

    print(f"Mitochondrial percentage stats:")
    print(f"Mean: {adata.obs['percent_mito'].mean():.3f}")
    print(f"Median: {adata.obs['percent_mito'].median():.3f}")
    print(f"95th percentile: {adata.obs['percent_mito'].quantile(0.95):.3f}")
    print(f"Max: {adata.obs['percent_mito'].max():.3f}")
    # print("After filtering for high mitochondrial counts:")
    # print(adata)

    # Filter out genes that are not present in any cells:
    sc.pp.filter_genes(adata, min_cells=3)
    print("After filtering unexpressed genes:")
    print(adata)

    # Create violin plot and save it
    sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'], jitter=0.4, multi_panel=True, show=False)
    plt.savefig(f"{fig_dir}/violin_plot.pdf", bbox_inches='tight', pad_inches=0.1)
    plt.close()

# Add a function here for Wolbachia titer:


def analyze_filtered_adata(adata, output_prefix=""):

    # Filter out the doublets:
    print(f"Cells before doublet removal: {adata.n_obs}")
    adata = adata[~adata.obs['predicted_doublet']]
    print(f"Cells after doublet removal: {adata.n_obs}")

    #Normalize the filtered counts:
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    sc.pp.log1p(adata)

    # Identify highly variable genes:
    # flavor="cell_ranger" is consistent with Seurat and flavor="seurat" is not consistent with Seurat
    sc.pp.highly_variable_genes(adata, min_mean=0.01, max_mean=10, min_disp=.1, n_top_genes=1000, n_bins=20, flavor="seurat")
    sc.pl.highly_variable_genes(adata, show=False)
    # plt.savefig(f"{output_prefix}/highly_variable_genes.pdf", bbox_inches='tight', pad_inches=0.1)
    plt.close()

    # Check how many HVGs were found
    n_hvg = adata.var['highly_variable'].sum()
    print(f"Found {n_hvg} highly variable genes")
    
    sc.pp.scale(adata, max_value=10)

    # Cluster the cells using leiden clustering
    sc.tl.pca(adata, svd_solver='arpack', use_highly_variable=True, n_comps=10)
    sc.pp.neighbors(adata, n_neighbors=30, n_pcs=10, knn=True)
    sc.tl.leiden(adata)

    # PCA projection:
    sc.pl.pca(adata, color='leiden', show=False)
    plt.close()

    sc.tl.umap(adata)
    sc.pl.umap(adata, color='leiden', show=False)
    # plt.savefig(f"{output_prefix}/umap_leiden.pdf", bbox_inches='tight', pad_inches=0.1)
    plt.close()

    # Compute and plot the variance explained by the PC subspaces.
    sc.pl.pca_variance_ratio(adata, show=False) # The variance explained by each principal component is a measure of how well a projection to that component represents the data. 
    plt.close()

    print('Final Data')
    print(adata)
    return adata 

def process_data_with_metrics(key, matrix, log_to_file=True):
    # Load the matrix as a sparse matrix
    adata = ad.read_h5ad(matrix)

    # Set up logging if requested
    if log_to_file:
        log_file = open(f"{output_dir}/{key}_stats.txt", 'w')
        original_stdout = sys.stdout
        sys.stdout = log_file
    
    try:
        # Set the output directory for all scanpy objects
        # sc.settings.figdir = output_dir
        
        # Store metrics at each stage
        all_metrics = []
        
        # Raw data metrics
        raw_metrics = calculate_qc_metrics(adata, sample_name=key, stage="raw")
        all_metrics.append(raw_metrics)
        print_qc_summary(raw_metrics)
        
        # Apply QC filtering
        filtered_adata = qc_plots(adata, fig_dir)
        # Identify doublets:
        filtered_adata = identify_doublets(adata, fig_dir)

        # Post-filtering metrics
        filtered_metrics = calculate_qc_metrics(filtered_adata, sample_name=key, stage="filtered")
        all_metrics.append(filtered_metrics)
        print_qc_summary(filtered_metrics)
        
        # Apply normalization and analysis
        normalized_adata = analyze_filtered_adata(filtered_adata, output_dir)
        
        # Save metrics to CSV
        metrics_csv_path = f"{fig_dir}/{key}_qc_metrics.csv"
        save_metrics_to_csv(all_metrics, metrics_csv_path)
        
        h5ad_path = f"{output.rstrip('/')}/{key}.h5ad"  # Fix the output path
        # Save the filtered adata object
        normalized_adata.write(output)
        
        # # Gzip the saved .h5ad file
        # gzipped_path = f"{h5ad_path}.gz"
        # with open(h5ad_path, 'rb') as f_in:
        #     with gzip.open(gzipped_path, 'wb') as f_out:
        #         shutil.copyfileobj(f_in, f_out)

        # Remove the uncompressed version to save space
        # os.remove(h5ad_path)
        # print(f"Saved compressed file: {gzipped_path}")
        
    finally:
        # Always restore stdout
        if log_to_file:
            sys.stdout = original_stdout
            log_file.close()
    
    return all_metrics

print(f"Analyzing {sample_name}, {input}")
sample_metrics = process_data_with_metrics(sample_name, input)
print(f"Completed {sample_name}")

