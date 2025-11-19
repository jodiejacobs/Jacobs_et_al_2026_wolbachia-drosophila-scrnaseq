import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np

# Install packages for analysis and plotting
from scipy.io import mmread
from sklearn.decomposition import TruncatedSVD
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.sparse import csr_matrix
import os
matplotlib.rcParams.update({'font.size': 22})


# File with all the matrix paths, open as dictionary
# matrix_paths = pd.read_csv('/Users/jodiejacobs/Downloads/matrix_paths.csv').to_dict(orient='records')

# Matix file paths:
matrix = {
    # # Cellranger 10X
    # 'cellranger_JW18DOX-Ctrl-1_X': 'scRNAseq_pilot_study/snakemake_pipeline/results_all_cellranger/10x/JW18DOX-Ctrl-1_10x/outs/filtered_feature_bc_matrix/matrix.mtx.gz',
    # 'cellranger_JW18DOX-Ctrl-2_X': 'scRNAseq_pilot_study/snakemake_pipeline/results_all_cellranger/10x/JW18DOX-Ctrl-2_10x/outs/filtered_feature_bc_matrix/matrix.mtx.gz',
    # 'cellranger_JW18DOX-SV-1_X': 'scRNAseq_pilot_study/snakemake_pipeline/results_all_cellranger/10x/JW18DOX-SV-1_10x/outs/filtered_feature_bc_matrix/matrix.mtx.gz',
    # 'cellranger_JW18DOX-SV-2_X': 'scRNAseq_pilot_study/snakemake_pipeline/results_all_cellranger/10x/JW18DOX-SV-2_10x/outs/filtered_feature_bc_matrix/matrix.mtx.gz',
    # 'cellranger_JW18wMel-Ctrl-1_X': 'scRNAseq_pilot_study/snakemake_pipeline/results_all_cellranger/10x/JW18wMel-Ctrl-1_10x/outs/filtered_feature_bc_matrix/matrix.mtx.gz',
    # 'cellranger_JW18wMel-Ctrl-2_X': 'scRNAseq_pilot_study/snakemake_pipeline/results_all_cellranger/10x/JW18wMel-Ctrl-2_10x/outs/filtered_feature_bc_matrix/matrix.mtx.gz',
    # 'cellranger_JW18wMel-SV-1_X': 'scRNAseq_pilot_study/snakemake_pipeline/results_all_cellranger/10x/JW18wMel-SV-1_10x/outs/filtered_feature_bc_matrix/matrix.mtx.gz',
    # 'cellranger_JW18wMel-SV-2_X': 'scRNAseq_pilot_study/snakemake_pipeline/results_all_cellranger/10x/JW18wMel-SV-2_10x/outs/filtered_feature_bc_matrix/matrix.mtx.gz',

    # # Cellranger pipseq
    # 'cellranger_JW18DOX-Ctrl-1_P': 'scRNAseq_pilot_study/snakemake_pipeline/results_all_cellranger/pipseq/JW18DOX-Ctrl-1_pipseq/outs/filtered_feature_bc_matrix/matrix.mtx.gz',
    # 'cellranger_JW18DOX-Ctrl-2_P': 'scRNAseq_pilot_study/snakemake_pipeline/results_all_cellranger/pipseq/JW18DOX-Ctrl-2_pipseq/outs/filtered_feature_bc_matrix/matrix.mtx.gz',
    # 'cellranger_JW18DOX-SV-1_P': 'scRNAseq_pilot_study/snakemake_pipeline/results_all_cellranger/pipseq/JW18DOX-SV-1_pipseq/outs/filtered_feature_bc_matrix/matrix.mtx.gz',
    # 'cellranger_JW18DOX-SV-2_P': 'scRNAseq_pilot_study/snakemake_pipeline/results_all_cellranger/pipseq/JW18DOX-SV-2_pipseq/outs/filtered_feature_bc_matrix/matrix.mtx.gz',
    # 'cellranger_JW18wMel-Ctrl-1_P': 'scRNAseq_pilot_study/snakemake_pipeline/results_all_cellranger/pipseq/JW18wMel-Ctrl-1_pipseq/outs/filtered_feature_bc_matrix/matrix.mtx.gz',
    # 'cellranger_JW18wMel-Ctrl-2_P': 'scRNAseq_pilot_study/snakemake_pipeline/results_all_cellranger/pipseq/JW18wMel-Ctrl-2_pipseq/outs/filtered_feature_bc_matrix/matrix.mtx.gz',
    # 'cellranger_JW18wMel-SV-1_P': 'scRNAseq_pilot_study/snakemake_pipeline/results_all_cellranger/pipseq/JW18wMel-SV-1_pipseq/outs/filtered_feature_bc_matrix/matrix.mtx.gz',
    # 'cellranger_JW18wMel-SV-2_P': 'scRNAseq_pilot_study/snakemake_pipeline/results_all_cellranger/pipseq/JW18wMel-SV-2_pipseq/outs/filtered_feature_bc_matrix/matrix.mtx.gz',

    # Kallisto bustools 10X
    'kallisto_JW18DOX-Ctrl-1_X': 'scRNAseq_pilot_study/snakemake_pipeline/results_kallisto_bustools/10x/JW18DOX-Ctrl-1_10x/counts_unfiltered/cells_x_genes.mtx',
    'kallisto_JW18DOX-Ctrl-2_X': 'scRNAseq_pilot_study/snakemake_pipeline/results_kallisto_bustools/10x/JW18DOX-Ctrl-2_10x/counts_unfiltered/cells_x_genes.mtx',
    'kallisto_JW18DOX-SV-1_X': 'scRNAseq_pilot_study/snakemake_pipeline/results_kallisto_bustools/10x/JW18DOX-SV-1_10x/counts_unfiltered/cells_x_genes.mtx',
    'kallisto_JW18DOX-SV-2_X': 'scRNAseq_pilot_study/snakemake_pipeline/results_kallisto_bustools/10x/JW18DOX-SV-2_10x/counts_unfiltered/cells_x_genes.mtx',
    'kallisto_JW18wMel-Ctrl-1_X': 'scRNAseq_pilot_study/snakemake_pipeline/results_kallisto_bustools/10x/JW18wMel-Ctrl-1_10x/counts_unfiltered/cells_x_genes.mtx',
    'kallisto_JW18wMel-Ctrl-2_X': 'scRNAseq_pilot_study/snakemake_pipeline/results_kallisto_bustools/10x/JW18wMel-Ctrl-2_10x/counts_unfiltered/cells_x_genes.mtx',
    'kallisto_JW18wMel-SV-1_X': 'scRNAseq_pilot_study/snakemake_pipeline/results_kallisto_bustools/10x/JW18wMel-SV-1_10x/counts_unfiltered/cells_x_genes.mtx',
    'kallisto_JW18wMel-SV-2_X': 'scRNAseq_pilot_study/snakemake_pipeline/results_kallisto_bustools/10x/JW18wMel-SV-2_10x/counts_unfiltered/cells_x_genes.mtx',
    # Kallisto bustools pipseq
    'kallisto_JW18DOX-Ctrl-1_P': 'scRNAseq_pilot_study/snakemake_pipeline/results_kallisto_bustools/pipseq/JW18DOX-Ctrl-1_pipseq/counts_unfiltered/cells_x_genes.mtx',
    'kallisto_JW18DOX-Ctrl-2_P': 'scRNAseq_pilot_study/snakemake_pipeline/results_kallisto_bustools/pipseq/JW18DOX-Ctrl-2_pipseq/counts_unfiltered/cells_x_genes.mtx',
    'kallisto_JW18DOX-SV-1_P': 'scRNAseq_pilot_study/snakemake_pipeline/results_kallisto_bustools/pipseq/JW18DOX-SV-1_pipseq/counts_unfiltered/cells_x_genes.mtx',
    'kallisto_JW18DOX-SV-2_P': 'scRNAseq_pilot_study/snakemake_pipeline/results_kallisto_bustools/pipseq/JW18DOX-SV-2_pipseq/counts_unfiltered/cells_x_genes.mtx',
    'kallisto_JW18wMel-Ctrl-1_P': 'scRNAseq_pilot_study/snakemake_pipeline/results_kallisto_bustools/pipseq/JW18wMel-Ctrl-1_pipseq/counts_unfiltered/cells_x_genes.mtx',
    'kallisto_JW18wMel-Ctrl-2_P': 'scRNAseq_pilot_study/snakemake_pipeline/results_kallisto_bustools/pipseq/JW18wMel-Ctrl-2_pipseq/counts_unfiltered/cells_x_genes.mtx',
    'kallisto_JW18wMel-SV-1_P': 'scRNAseq_pilot_study/snakemake_pipeline/results_kallisto_bustools/pipseq/JW18wMel-SV-1_pipseq/counts_unfiltered/cells_x_genes.mtx',
    'kallisto_JW18wMel-SV-2_P': 'scRNAseq_pilot_study/snakemake_pipeline/results_kallisto_bustools/pipseq/JW18wMel-SV-2_pipseq/counts_unfiltered/cells_x_genes.mtx',
    'pipseeker_JW18DOX-Ctrl-1_P': 'scRNAseq_pilot_study/snakemake_pipeline/results_pip_and_cellranger/pipseeker/JW18DOX-Ctrl-1_pipseq/raw_matrix/matrix.mtx.gz',
    'pipseeker_JW18DOX-Ctrl-2_P': 'scRNAseq_pilot_study/snakemake_pipeline/results_pip_and_cellranger/pipseeker/JW18DOX-Ctrl-2_pipseq/raw_matrix/matrix.mtx.gz',
    'pipseeker_JW18DOX-SV-1_P': 'scRNAseq_pilot_study/snakemake_pipeline/results_pip_and_cellranger/pipseeker/JW18DOX-SV-1_pipseq/raw_matrix/matrix.mtx.gz',
    'pipseeker_JW18DOX-SV-2_P': 'scRNAseq_pilot_study/snakemake_pipeline/results_pip_and_cellranger/pipseeker/JW18DOX-SV-2_pipseq/raw_matrix/matrix.mtx.gz',
    'pipseeker_JW18wMel-Ctrl-1_P': 'scRNAseq_pilot_study/snakemake_pipeline/results_pip_and_cellranger/pipseeker/JW18wMel-Ctrl-1_pipseq/raw_matrix/matrix.mtx.gz',
    'pipseeker_JW18wMel-Ctrl-2_P': 'scRNAseq_pilot_study/snakemake_pipeline/results_pip_and_cellranger/pipseeker/JW18wMel-Ctrl-2_pipseq/raw_matrix/matrix.mtx.gz',
    'pipseeker_JW18wMel-SV-1_P': 'scRNAseq_pilot_study/snakemake_pipeline/results_pip_and_cellranger/pipseeker/JW18wMel-SV-1_pipseq/raw_matrix/matrix.mtx.gz',
    'pipseeker_JW18wMel-SV-2_P': 'scRNAseq_pilot_study/snakemake_pipeline/results_pip_and_cellranger/pipseeker/JW18wMel-SV-2_pipseq/raw_matrix/matrix.mtx.gz'
}

# density display for PCA plot
from scipy.interpolate import interpn

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

    sc = ax.scatter( x, y, c=z, **kwargs )
    return sc


### Basic QC, Check the cells in 2D space
def qc_plots(mtx_sparse, X, output_prefix="/content"):
    """
    Perform basic QC and visualization for a count matrix and its PCA projection.
    Saves plots to disk.
    """
    # Plot the cells in the 2D PCA projection
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.scatter(X[:,0], X[:,1], alpha=0.5, c="green")
    plt.axis('off')
    # plt.show()
    plt.savefig(f"{output_prefix}pca_plot.pdf", bbox_inches='tight', pad_inches=0.1)
    plt.close(fig)

    # density display for PCA plot
    fig, ax = plt.subplots(figsize=(7,7))
    x = X[:,0]
    y = X[:,1]
    sc = density_scatter(x, y, ax=ax, cmap="Greens")
    fig.colorbar(sc, ax=ax)
    plt.axis('off')
    # plt.show()
    plt.savefig(f"{output_prefix}density_scatter.pdf", bbox_inches='tight', pad_inches=0.1)
    plt.close(fig)

    # Test for library saturation:
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.scatter(np.asarray(mtx_sparse.sum(axis=1))[:,0], np.asarray(np.sum(mtx_sparse>0, axis=1))[:,0], color="green", alpha=0.01)
    ax.set_xlabel("UMI Counts")
    ax.set_ylabel("Genes Detected")
    ax.set_xscale('log')
    ax.set_yscale('log') # Removed nonposy='clip'
    # ax.set_xlim((0.5, 4500))
    # ax.set_ylim((0.5,2000))
    # plt.show()
    plt.savefig(f"{output_prefix}genes_vs_umi.pdf", bbox_inches='tight', pad_inches=0.1)
    plt.close(fig)

    # Knee plot for library saturation
    knee = np.sort((np.array(mtx_sparse.sum(axis=1))).flatten())[::-1]
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.loglog(knee, range(len(knee)),linewidth=5, color="g")
    ax.set_xlabel("UMI Counts")
    ax.set_ylabel("Set of Barcodes")
    plt.grid(True, which="both")
    # plt.show()
    plt.savefig(f"{output_prefix}knee_plot.pdf", bbox_inches='tight', pad_inches=0.1)
    plt.close(fig)



# JW18DOX_10X_KB='/content/10x/cells_x_genes.mtx'
# matrix = JW18DOX_10X_KB
# # Load the matrix as a sparse matrix
# mtx_sparse = mmread(matrix).T

# # Perform PCA on the sparse matrix
# svd = TruncatedSVD(n_components=2)
# X = svd.fit_transform(mtx_sparse)

# # Run QC plots (Note: qc_plots might need adjustments to handle sparse input if not already)
# qc_plots(mtx_sparse, X, '/content/JW18DOX_10X_KB')

# # Process pipseek matrix:

# JW18DOX_pip_KB='scRNAseq_pilot_study/snakemake_pipeline/results_kallisto_bustools/pipseq/JW18DOX-Ctrl-1_pipseq/counts_unfiltered/cells_x_genes.mtx'
# matrix = JW18DOX_pip_KB
# # Load the matrix as a sparse matrix
# mtx_sparse = mmread(matrix).T

# # Perform PCA on the sparse matrix
# svd = TruncatedSVD(n_components=2)
# X = svd.fit_transform(mtx_sparse)

# qc_plots(mtx_sparse, X, '/content/JW18DOX_pip_KB')

def load_matrix(key, dictionary):
  matrix = dictionary[key]
  # Load the matrix as a sparse matrix
  mtx_sparse = mmread(matrix).T

  # Perform PCA on the sparse matrix
  svd = TruncatedSVD(n_components=2)
  X = svd.fit_transform(mtx_sparse)

  output_dir = f'scRNAseq_pilot_study/snakemake_pipeline/qc_matrix/results/{key}/'
  os.makedirs(output_dir, exist_ok=True)

  qc_plots(mtx_sparse, X, output_dir)

keys = list(matrix.keys())
for key in keys:
  load_matrix(key, matrix)

