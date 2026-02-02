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

def integrate(files, out_path, fig_dir, sample, batch_key, min_cells, min_genes, 
              calculate_titer=True, n_pcs=30, bio_condition=None):
    """
    files: list of h5ad file paths
    bio_condition: optional filter like 'JW18DOX-Ctrl' to compare only those samples
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

    # Extract metadata from sample names
    combined.obs['cell_line'] = combined.obs[batch_key].str.extract(r'(JW18DOX|JW18wMel)')[0]
    combined.obs['treatment'] = combined.obs[batch_key].str.extract(r'-(Ctrl|SV)')[0]
    combined.obs['replicate'] = combined.obs[batch_key].str.extract(r'-(Ctrl|SV)-(\d+|D\d+)')[1]
    combined.obs['method'] = combined.obs[batch_key].str.extract(r'_(10x|pipseq)$')[0]
    
    # Create biological condition column
    combined.obs['bio_condition'] = combined.obs['cell_line'] + '-' + combined.obs['treatment']
    
    # Filter to specific biological condition if requested
    if bio_condition:
        print(f"Filtering to biological condition: {bio_condition}")
        combined = combined[combined.obs['bio_condition'] == bio_condition].copy()
        sample = f"{sample}_{bio_condition}"
    
    print(f"Combined data shape for {sample}: {combined.shape}")
    print(f"Samples included: {combined.obs[batch_key].unique()}")
    print(f"Methods: {combined.obs['method'].value_counts()}")
    print(f"Data range: {combined.X.min():.3f} to {combined.X.max():.3f}")
    
    # Filter out bacterial genes
    bacteria_genes = ['GQX67_00940', 'GQX67_05945'] + [gene for gene in combined.var_names if gene.startswith('16S_')]
    bacteria_mask = combined.var_names.isin(bacteria_genes)
    combined = combined[:, ~bacteria_mask]

    # Basic preprocessing
    print("Performing basic preprocessing...")
    combined.X = np.nan_to_num(combined.X, nan=0.0)

    sc.pp.filter_cells(combined, min_genes=min_genes)
    sc.pp.filter_genes(combined, min_cells=min_cells)
          
    # Find highly variable genes
    print("Finding highly variable genes...")
    sc.pp.highly_variable_genes(combined, flavor='seurat', n_top_genes=2000)
    combined = combined[:, combined.var.highly_variable]
        
    # Run PCA
    print("Running PCA...")
    sc.pp.pca(combined, n_comps=n_pcs)
    
    # Save unintegrated version
    combined_unintegrated = combined.copy()
    sc.pp.neighbors(combined_unintegrated, n_pcs=n_pcs)
    sc.tl.umap(combined_unintegrated)
    sc.pl.umap(combined_unintegrated, color=['method', batch_key], 
               save=f'_{sample}_before_integration.pdf', ncols=2)

    # Batch correction
    print(f"Running BBKNN batch correction on {batch_key}...")
    bbknn.bbknn(combined, batch_key=batch_key, n_pcs=n_pcs, neighbors_within_batch=5)
            
    # Run UMAP and clustering
    print("Running UMAP and Leiden clustering...")
    sc.tl.umap(combined)
    sc.tl.leiden(combined, resolution=0.8)
    
    # Save the integrated object
    print(f"Saving integrated object for {sample} to {out_path}")
    combined.write(out_path)
    
    # Generate comprehensive comparison plots
    print("Generating method comparison analysis...")
    
    # 1. Side-by-side UMAPs with same coordinates
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    sc.pl.umap(combined, color='method', ax=axes[0], show=False, 
               title='Library Prep Method', frameon=False)
    sc.pl.umap(combined, color='leiden', ax=axes[1], show=False, 
               title='Leiden Clusters', frameon=False, legend_loc='on data')
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'umap_{sample}_method_overview.pdf'))
    plt.close()
    
    # 2. Split UMAPs by method - same shape
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    sc.pl.umap(combined[combined.obs['method'] == '10x'], 
               color='leiden', ax=axes[0], show=False, 
               title='10x Genomics', frameon=False)
    sc.pl.umap(combined[combined.obs['method'] == 'pipseq'], 
               color='leiden', ax=axes[1], show=False, 
               title='PIPseq', frameon=False)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'umap_{sample}_split_by_method.pdf'))
    plt.close()
    
    # 3. Cluster composition analysis
    cluster_comp = pd.crosstab(combined.obs['leiden'], combined.obs['method'], 
                               normalize='index') * 100
    
    # Perform chi-square test for cluster distribution
    contingency = pd.crosstab(combined.obs['leiden'], combined.obs['method'])
    chi2, p_value, dof, expected = chi2_contingency(contingency)
    
    print("\n" + "="*60)
    print("CLUSTER COMPOSITION BY METHOD")
    print("="*60)
    print(cluster_comp)
    print(f"\nChi-square test: χ² = {chi2:.2f}, p = {p_value:.2e}")
    print(f"Methods show {'SIGNIFICANT' if p_value < 0.05 else 'NO'} difference in cluster distribution")
    
    # Plot cluster composition as stacked bar
    fig, ax = plt.subplots(figsize=(10, 6))
    cluster_comp.plot(kind='bar', stacked=True, ax=ax, 
                      color=['#1f77b4', '#ff7f0e'])
    ax.set_xlabel('Leiden Cluster')
    ax.set_ylabel('Percentage of cells')
    ax.set_title(f'Cluster composition by method - {sample}')
    ax.legend(title='Method')
    ax.axhline(y=50, color='gray', linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'cluster_composition_{sample}.pdf'))
    plt.close()
    
    # 4. Cell count comparison
    method_counts = combined.obs.groupby(['method', 'leiden']).size().unstack(fill_value=0)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    method_counts.T.plot(kind='bar', ax=ax, color=['#1f77b4', '#ff7f0e'])
    ax.set_xlabel('Leiden Cluster')
    ax.set_ylabel('Number of cells')
    ax.set_title(f'Cell counts by method and cluster - {sample}')
    ax.legend(title='Method')
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'cell_counts_{sample}.pdf'))
    plt.close()
    
    # 5. If multiple biological conditions, compare across them
    if combined.obs['bio_condition'].nunique() > 1:
        plot_condition_comparison(combined, fig_dir, sample)
    
    # 6. Wolbachia titer analysis
    if 'wolbachia_titer' in combined.obs.columns:
        print("\n" + "="*60)
        print("WOLBACHIA TITER ANALYSIS")
        print("="*60)
        
        # Titer by method
        titer_by_method = combined.obs.groupby('method')['wolbachia_titer'].agg(['mean', 'median', 'std'])
        print("\nTiter by method:")
        print(titer_by_method)
        
        # Statistical test
        titer_10x = combined.obs[combined.obs['method'] == '10x']['wolbachia_titer'].dropna()
        titer_pipseq = combined.obs[combined.obs['method'] == 'pipseq']['wolbachia_titer'].dropna()
        u_stat, p_val = mannwhitneyu(titer_10x, titer_pipseq, alternative='two-sided')
        print(f"\nMann-Whitney U test: U = {u_stat:.2f}, p = {p_val:.2e}")
        print(f"Titer shows {'SIGNIFICANT' if p_val < 0.05 else 'NO'} difference between methods")
        
        # Violin plot: titer by method
        fig, ax = plt.subplots(figsize=(8, 6))
        sc.pl.violin(combined, 'wolbachia_titer', groupby='method', ax=ax, show=False)
        ax.set_title(f'Wolbachia titer by method - {sample}')
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f'titer_by_method_{sample}.pdf'))
        plt.close()
        
        # Titer by cluster and method
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))
        
        # Split violin plots
        for i, method in enumerate(['10x', 'pipseq']):
            subset = combined[combined.obs['method'] == method]
            sc.pl.violin(subset, 'wolbachia_titer', groupby='leiden', 
                        ax=axes[i], show=False)
            axes[i].set_title(f'{method} - Titer by cluster')
        
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f'titer_by_cluster_split_{sample}.pdf'))
        plt.close()
        
        # Heatmap: mean titer by cluster and method
        titer_heatmap = combined.obs.groupby(['leiden', 'method'])['wolbachia_titer'].mean().unstack()
        
        fig, ax = plt.subplots(figsize=(8, 10))
        sns.heatmap(titer_heatmap, annot=True, fmt='.3f', cmap='viridis', ax=ax)
        ax.set_title(f'Mean Wolbachia titer by cluster and method - {sample}')
        ax.set_xlabel('Method')
        ax.set_ylabel('Leiden Cluster')
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f'titer_heatmap_{sample}.pdf'))
        plt.close()
        
        # Correlation plot: titer between methods per cluster
        titer_corr_data = []
        for cluster in combined.obs['leiden'].unique():
            cluster_cells = combined[combined.obs['leiden'] == cluster]
            titer_10x = cluster_cells[cluster_cells.obs['method'] == '10x'].obs['wolbachia_titer'].mean()
            titer_pipseq = cluster_cells[cluster_cells.obs['method'] == 'pipseq'].obs['wolbachia_titer'].mean()
            titer_corr_data.append({'cluster': cluster, '10x': titer_10x, 'pipseq': titer_pipseq})
        
        titer_corr_df = pd.DataFrame(titer_corr_data)
        
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.scatter(titer_corr_df['10x'], titer_corr_df['pipseq'], s=100, alpha=0.6)
        for idx, row in titer_corr_df.iterrows():
            ax.annotate(row['cluster'], (row['10x'], row['pipseq']), 
                       xytext=(5, 5), textcoords='offset points', fontsize=8)
        
        # Add diagonal line
        max_val = max(titer_corr_df['10x'].max(), titer_corr_df['pipseq'].max())
        ax.plot([0, max_val], [0, max_val], 'r--', alpha=0.5, label='y=x')
        
        ax.set_xlabel('Mean titer - 10x')
        ax.set_ylabel('Mean titer - PIPseq')
        ax.set_title(f'Titer correlation between methods by cluster - {sample}')
        ax.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f'titer_correlation_{sample}.pdf'))
        plt.close()
        
        # UMAP split by method with titer
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        sc.pl.umap(combined[combined.obs['method'] == '10x'], 
                   color='wolbachia_titer', ax=axes[0], show=False, 
                   title='10x Genomics - Wolbachia titer', frameon=False, vmax=np.percentile(combined.obs['wolbachia_titer'], 95))
        sc.pl.umap(combined[combined.obs['method'] == 'pipseq'], 
                   color='wolbachia_titer', ax=axes[1], show=False, 
                   title='PIPseq - Wolbachia titer', frameon=False, vmax=np.percentile(combined.obs['wolbachia_titer'], 95))
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f'umap_titer_split_{sample}.pdf'))
        plt.close()
    
    # Standard diagnostic plots
    sc.pl.umap(combined, color=batch_key, save=f'_{sample}_by_sample.pdf')
    sc.pl.umap(combined, color='leiden', save=f'_{sample}_clusters.pdf')
    
    # Summary statistics
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"Total cells: {combined.n_obs}")
    print(f"Total genes: {combined.n_vars}")
    print(f"\nCells per method:")
    print(combined.obs['method'].value_counts())
    print(f"\nClusters: {combined.obs['leiden'].nunique()}")
    print(f"\nCells per cluster:")
    print(combined.obs['leiden'].value_counts().sort_index())
    
    if 'wolbachia_titer' in combined.obs.columns:
        n_infected = np.sum(combined.obs['wolbachia_titer'] > 0)
        print(f"\nInfected cells: {n_infected} ({n_infected/combined.n_obs*100:.2f}%)")


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
