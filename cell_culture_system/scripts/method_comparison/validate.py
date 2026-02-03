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
import seaborn as sns
from scipy.stats import chi2_contingency, mannwhitneyu

def open_adata(file_path):
    """Open an h5ad file and return the AnnData object"""
    adata = sc.read_h5ad(file_path)
    print(f"Opened {file_path} with {adata.n_obs} cells and {adata.n_vars} genes.")
    print(adata)
    return adata

def filter_adata(adata, conditions):
    """Filter AnnData based on given conditions"""
    for key, value in conditions.items():
        adata = adata[adata.obs[key] == value]
    return adata


def plot_condition_comparison(combined, fig_dir, sample):
    """Additional plots when multiple biological conditions are present"""
    
    # Stacked bar: cluster composition by method and condition
    fig, ax = plt.subplots(figsize=(12, 6))
    
    comp_data = pd.crosstab([combined.obs['bio_condition'], combined.obs['leiden']], 
                            combined.obs['method'], normalize='index') * 100
    
    comp_data.plot(kind='bar', stacked=True, ax=ax, color=['#1f77b4', '#ff7f0e'])
    ax.set_xlabel('Biological Condition - Cluster')
    ax.set_ylabel('Percentage of cells')
    ax.set_title(f'Cluster composition across conditions - {sample}')
    ax.legend(title='Method')
    ax.axhline(y=50, color='gray', linestyle='--', alpha=0.5)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'condition_comparison_{sample}.pdf'))
    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Compare library prep methods with detailed analysis')

    parser.add_argument('--files', required=True, nargs='+', type=str,
                        help='List of h5ad files to integrate')
    parser.add_argument('--sample', type=str, default='method_comparison',
                        help='Sample type label')
    parser.add_argument('--bio_condition', type=str, default=None,
                        help='Filter to specific biological condition (e.g., JW18DOX-Ctrl)')
    parser.add_argument('--batch_key', type=str, default='batch',
                        help='Key in .obs to use for batch information')
    parser.add_argument('--min_cells', type=int, default=3,
                        help='Minimum cells per gene for filtering')
    parser.add_argument('--min_genes', type=int, default=200,
                        help='Minimum genes per cell for filtering')  
    parser.add_argument('--out_path', type=str, default='integrated.h5ad',
                        help='Path to save the integrated h5ad file')      
    parser.add_argument('--fig_dir', type=str, default='figures',
                        help='Directory to save figures')    
    parser.add_argument('--n_pcs', type=int, default=30,
                        help='Number of principal components to use')

    args = parser.parse_args()
    
    integrate(
        files=args.files,
        out_path=args.out_path,
        fig_dir=args.fig_dir,
        sample=args.sample, 
        batch_key=args.batch_key,
        min_cells=args.min_cells,
        min_genes=args.min_genes,
        n_pcs=args.n_pcs,
        bio_condition=args.bio_condition,
    )

if __name__ == "__main__":
    main()
