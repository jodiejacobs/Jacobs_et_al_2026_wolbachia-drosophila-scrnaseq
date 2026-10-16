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
              calculate_titer=True, n_pcs=30, bio_condition=None, leiden_resolution=0.5):
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
    
    # Concatenate all datasets
    combined = ad.concat(adatas, join='inner', merge='same', index_unique='-')

    # Extract metadata from sample names
    combined.obs['cell_line'] = combined.obs[batch_key].str.extract(r'(JW18DOX|JW18wMel)')[0]
    combined.obs['treatment'] = combined.obs[batch_key].str.extract(r'-(Ctrl|SV)')[0]

    # Extract timepoint if present (D7, D28, D56)
    combined.obs['timepoint'] = combined.obs[batch_key].str.extract(r'-(D\d+)-')[0]

    # Extract replicate - everything between treatment and method
    # This captures: 1, 2, 3, D7-1, D7-2, D28-1, etc.
    combined.obs['replicate'] = combined.obs[batch_key].str.extract(r'-(Ctrl|SV)-([^_]+)')[1]

    # Extract method
    combined.obs['method'] = combined.obs[batch_key].str.extract(r'_(10x|pipseq)$')[0]

    # Create biological condition column
    # Include timepoint if it exists, otherwise just cell_line-treatment
    combined.obs['bio_condition'] = combined.obs.apply(
        lambda row: f"{row['cell_line']}-{row['treatment']}-{row['timepoint']}" 
                    if pd.notna(row['timepoint']) 
                    else f"{row['cell_line']}-{row['treatment']}", 
        axis=1
    )    
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
    # combined = combined[:, combined.var.highly_variable]
        
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
    sc.tl.leiden(combined, resolution=leiden_resolution) # Manually adjust leiden resolution if needed based on cluster numbers
    
    # Save the integrated object
    print(f"Saving integrated object for {sample} to {out_path}")
    combined.write(out_path)
    
    # Generate Ctrl-only comparison plots
    plot_ctrl_comparison(combined, fig_dir, sample)
    
    # Generate full dataset analysis
    plot_full_dataset_analysis(combined, fig_dir, sample)
    
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


def plot_ctrl_comparison(combined, fig_dir, sample):
    """Compare Ctrl samples between 10x and PIPseq methods"""
    
    print("\n" + "="*60)
    print("CTRL SAMPLES - METHOD COMPARISON")
    print("="*60)
    
    # Filter to only Ctrl samples
    combined_ctrl = combined[combined.obs['treatment'] == 'Ctrl'].copy()
    
    if combined_ctrl.n_obs == 0:
        print("No Ctrl samples found, skipping comparison")
        return
    
    print(f"Ctrl samples: {combined_ctrl.n_obs} cells")
    print(f"Methods in Ctrl: {combined_ctrl.obs['method'].value_counts()}")
    
    # Get leiden colors for consistent coloring
    leiden_colors = []
    clusters = sorted(combined_ctrl.obs['leiden'].unique())
    cmap = plt.cm.get_cmap('tab20')
    for i, cluster in enumerate(clusters):
        leiden_colors.append(cmap(i % 20))
    
    # 1. UMAP colored by sample type (bio_condition)
    sc.pl.umap(combined_ctrl, color='bio_condition', 
               save=f'_{sample}_ctrl_by_bio_condition.pdf',
               title='Ctrl samples by biological condition')
    
    # 2. UMAP by titer
    if 'wolbachia_titer' in combined_ctrl.obs.columns:
        sc.pl.umap(combined_ctrl, color='wolbachia_titer', 
                   save=f'_{sample}_ctrl_by_titer.pdf',
                   title='Ctrl samples - Wolbachia titer',
                   vmax=np.percentile(combined_ctrl.obs['wolbachia_titer'].dropna(), 95))
    
    # 3. UMAP by leiden
    sc.pl.umap(combined_ctrl, color='leiden', 
               save=f'_{sample}_ctrl_by_leiden.pdf',
               title='Ctrl samples - Leiden clusters',
               legend_loc='on data')
    
    # 4. UMAP by pipseq vs 10x
    sc.pl.umap(combined_ctrl, color='method', 
               save=f'_{sample}_ctrl_by_method.pdf',
               title='Ctrl samples - Library prep method')
    
    # 5. Split UMAPs by method - same shape
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    sc.pl.umap(combined_ctrl[combined_ctrl.obs['method'] == '10x'], 
               color='leiden', ax=axes[0], show=False, 
               title='10x Genomics - Ctrl', frameon=False)
    sc.pl.umap(combined_ctrl[combined_ctrl.obs['method'] == 'pipseq'], 
               color='leiden', ax=axes[1], show=False, 
               title='PIPseq - Ctrl', frameon=False)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'umap_{sample}_ctrl_split_by_method.pdf'))
    plt.close()
    
    # 6. UMAP by cell cycle stage (if present)
    if 'phase' in combined_ctrl.obs.columns:
        sc.pl.umap(combined_ctrl, color='phase', 
                   save=f'_{sample}_ctrl_by_cellcycle.pdf',
                   title='Ctrl samples - Cell cycle phase')
    elif 'cyclum_theta' in combined_ctrl.obs.columns:
        sc.pl.umap(combined_ctrl, color='cyclum_theta', 
                   save=f'_{sample}_ctrl_by_cyclum_theta.pdf',
                   title='Ctrl samples - Cyclum theta',
                   cmap='twilight')
    
    # 7. Wolbachia titer comparison - wMel Ctrl only (10x vs pipseq)
    if 'wolbachia_titer' in combined_ctrl.obs.columns:
        wmel_ctrl = combined_ctrl[combined_ctrl.obs['cell_line'] == 'JW18wMel']
        
        if wmel_ctrl.n_obs > 0:
            fig, ax = plt.subplots(figsize=(8, 6))
            
            # Box plot with strip plot
            sns.boxplot(data=wmel_ctrl.obs, x='method', y='wolbachia_titer', 
                       ax=ax, palette=['#1f77b4', '#ff7f0e'])
            sns.stripplot(data=wmel_ctrl.obs, x='method', y='wolbachia_titer',
                         ax=ax, color='black', alpha=0.3, size=2)
            
            ax.set_xlabel('Library Prep Method')
            ax.set_ylabel('Wolbachia Titer')
            ax.set_title('wMel Ctrl - Wolbachia titer comparison')
            
            # Add stats
            titer_10x = wmel_ctrl.obs[wmel_ctrl.obs['method'] == '10x']['wolbachia_titer'].dropna()
            titer_pipseq = wmel_ctrl.obs[wmel_ctrl.obs['method'] == 'pipseq']['wolbachia_titer'].dropna()
            if len(titer_10x) > 0 and len(titer_pipseq) > 0:
                u_stat, p_val = mannwhitneyu(titer_10x, titer_pipseq, alternative='two-sided')
                ax.text(0.5, 0.95, f'Mann-Whitney U: p = {p_val:.2e}',
                       transform=ax.transAxes, ha='center', va='top')
            
            plt.tight_layout()
            plt.savefig(os.path.join(fig_dir, f'titer_comparison_{sample}_wMel_ctrl.pdf'))
            plt.close()
    
    # 8. DOX titer comparison (10x vs pipseq)
    if 'wolbachia_titer' in combined_ctrl.obs.columns:
        dox_ctrl = combined_ctrl[combined_ctrl.obs['cell_line'] == 'JW18DOX']
        
        if dox_ctrl.n_obs > 0:
            fig, ax = plt.subplots(figsize=(8, 6))
            
            sns.boxplot(data=dox_ctrl.obs, x='method', y='wolbachia_titer', 
                       ax=ax, palette=['#1f77b4', '#ff7f0e'])
            sns.stripplot(data=dox_ctrl.obs, x='method', y='wolbachia_titer',
                         ax=ax, color='black', alpha=0.3, size=2)
            
            ax.set_xlabel('Library Prep Method')
            ax.set_ylabel('Wolbachia Titer')
            ax.set_title('DOX Ctrl - Wolbachia titer comparison')
            
            # Add stats
            titer_10x = dox_ctrl.obs[dox_ctrl.obs['method'] == '10x']['wolbachia_titer'].dropna()
            titer_pipseq = dox_ctrl.obs[dox_ctrl.obs['method'] == 'pipseq']['wolbachia_titer'].dropna()
            if len(titer_10x) > 0 and len(titer_pipseq) > 0:
                u_stat, p_val = mannwhitneyu(titer_10x, titer_pipseq, alternative='two-sided')
                ax.text(0.5, 0.95, f'Mann-Whitney U: p = {p_val:.2e}',
                       transform=ax.transAxes, ha='center', va='top')
            
            plt.tight_layout()
            plt.savefig(os.path.join(fig_dir, f'titer_comparison_{sample}_DOX_ctrl.pdf'))
            plt.close()
    
    # 9. Titer boxplot by cluster with leiden colors
    if 'wolbachia_titer' in combined_ctrl.obs.columns:
        fig, ax = plt.subplots(figsize=(12, 6))
        
        # Prepare data
        plot_data = combined_ctrl.obs[['leiden', 'wolbachia_titer']].copy()
        plot_data = plot_data.sort_values('leiden')
        
        clusters = sorted(plot_data['leiden'].unique())
        
        # Plot individual points first (so they appear under the boxes)
        sns.stripplot(data=plot_data, x='leiden', y='wolbachia_titer', 
                     color='black', alpha=0.3, size=2, ax=ax)
        
        # Plot box plot with leiden colors
        box_parts = ax.boxplot([plot_data[plot_data['leiden'] == cluster]['wolbachia_titer'].values 
                                for cluster in clusters],
                                positions=range(len(clusters)),
                                widths=0.6,
                                patch_artist=True,
                                whiskerprops=dict(alpha=0.7),
                                capprops=dict(alpha=0.7),
                                medianprops=dict(color='black', linewidth=2))
        
        # Color each box with its corresponding leiden color
        for patch, color in zip(box_parts['boxes'], leiden_colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        
        ax.set_xlabel('Leiden Cluster')
        ax.set_ylabel('Wolbachia Titer')
        ax.set_xticklabels(clusters)
        ax.set_title('Ctrl samples - Wolbachia titer by cluster')
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f'boxplot_{sample}_ctrl_titer_by_cluster.pdf'))
        plt.close()
    
    # 10. Number of genes per cluster (pipseq vs 10x)
    genes_per_cluster = combined_ctrl.obs.groupby(['leiden', 'method'])['n_genes'].mean().unstack()
    
    fig, ax = plt.subplots(figsize=(10, 6))
    genes_per_cluster.plot(kind='bar', ax=ax, color=['#1f77b4', '#ff7f0e'])
    ax.set_xlabel('Leiden Cluster')
    ax.set_ylabel('Mean number of genes per cell')
    ax.set_title('Ctrl samples - Genes per cluster by method')
    ax.legend(title='Method')
    plt.xticks(rotation=0)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'genes_per_cluster_{sample}_ctrl.pdf'))
    plt.close()
    
    # 11. Percentage of cells per cluster (pipseq vs 10x)
    cluster_comp = pd.crosstab(combined_ctrl.obs['leiden'], combined_ctrl.obs['method'], 
                               normalize='columns') * 100
    
    fig, ax = plt.subplots(figsize=(10, 6))
    cluster_comp.plot(kind='bar', ax=ax, color=['#1f77b4', '#ff7f0e'])
    ax.set_xlabel('Leiden Cluster')
    ax.set_ylabel('Percentage of cells')
    ax.set_title('Ctrl samples - Cell distribution by method')
    ax.legend(title='Method')
    plt.xticks(rotation=0)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'cell_percentage_per_cluster_{sample}_ctrl.pdf'))
    plt.close()
    
    # Chi-square test for cluster distribution
    contingency = pd.crosstab(combined_ctrl.obs['leiden'], combined_ctrl.obs['method'])
    chi2, p_value, dof, expected = chi2_contingency(contingency)
    
    print("\nCTRL SAMPLES - CLUSTER COMPOSITION BY METHOD")
    print(cluster_comp)
    print(f"\nChi-square test: χ² = {chi2:.2f}, p = {p_value:.2e}")
    print(f"Methods show {'SIGNIFICANT' if p_value < 0.05 else 'NO'} difference in cluster distribution")
    
    # 12. Side-by-side titer split by method
    if 'wolbachia_titer' in combined_ctrl.obs.columns:
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        
        vmax = np.percentile(combined_ctrl.obs['wolbachia_titer'].dropna(), 95)
        
        sc.pl.umap(combined_ctrl[combined_ctrl.obs['method'] == '10x'], 
                   color='wolbachia_titer', ax=axes[0], show=False, 
                   title='10x Genomics - Wolbachia titer', frameon=False, vmax=vmax)
        sc.pl.umap(combined_ctrl[combined_ctrl.obs['method'] == 'pipseq'], 
                   color='wolbachia_titer', ax=axes[1], show=False, 
                   title='PIPseq - Wolbachia titer', frameon=False, vmax=vmax)
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f'umap_{sample}_ctrl_titer_split.pdf'))
        plt.close()


def plot_full_dataset_analysis(combined, fig_dir, sample):
    """Analysis of full dataset: cell cycle, clustering, and Wolbachia infection"""
    
    print("\n" + "="*60)
    print("FULL DATASET ANALYSIS")
    print("="*60)
    
    # Get leiden colors for consistent coloring
    leiden_colors = []
    clusters = sorted(combined.obs['leiden'].unique())
    cmap = plt.cm.get_cmap('tab20')
    for i, cluster in enumerate(clusters):
        leiden_colors.append(cmap(i % 20))
    
    # 1. Leiden cluster plots
    sc.pl.umap(combined, color='leiden', 
               save=f'_{sample}_full_leiden_clusters.pdf',
               title='All samples - Leiden clusters',
               legend_loc='on data')
    
    # 2. UMAP by cell cycle stage
    if 'phase' in combined.obs.columns:
        sc.pl.umap(combined, color='phase', 
                   save=f'_{sample}_full_by_cellcycle_phase.pdf',
                   title='All samples - Cell cycle phase')
        
        # Cell cycle distribution per cluster
        fig, ax = plt.subplots(figsize=(12, 6))
        phase_comp = pd.crosstab(combined.obs['leiden'], combined.obs['phase'], 
                                normalize='index') * 100
        phase_comp.plot(kind='bar', stacked=True, ax=ax)
        ax.set_xlabel('Leiden Cluster')
        ax.set_ylabel('Percentage of cells')
        ax.set_title('Cell cycle phase distribution by cluster')
        ax.legend(title='Cell cycle phase')
        plt.xticks(rotation=0)
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f'cellcycle_by_cluster_{sample}_full.pdf'))
        plt.close()
    
    if 'cyclum_theta' in combined.obs.columns:
        sc.pl.umap(combined, color='cyclum_theta', 
                   save=f'_{sample}_full_by_cyclum_theta.pdf',
                   title='All samples - Cyclum theta',
                   cmap='twilight')
        
        # Cyclum theta violin plot by cluster
        fig, ax = plt.subplots(figsize=(14, 6))
        sc.pl.violin(combined, 'cyclum_theta', groupby='leiden', 
                    ax=ax, show=False, rotation=0)
        ax.set_title('Cyclum theta distribution by cluster')
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f'cyclum_theta_by_cluster_{sample}_full.pdf'))
        plt.close()
    
    # 3. Genes per cluster
    genes_per_cluster = combined.obs.groupby('leiden')['n_genes'].agg(['mean', 'std'])
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Bar plot with error bars
    x_pos = range(len(genes_per_cluster))
    ax.bar(x_pos, genes_per_cluster['mean'], 
           yerr=genes_per_cluster['std'],
           color=leiden_colors, alpha=0.7, capsize=5)
    
    ax.set_xlabel('Leiden Cluster')
    ax.set_ylabel('Mean number of genes per cell')
    ax.set_title('Gene count by cluster (all samples)')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(genes_per_cluster.index)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'genes_per_cluster_{sample}_full.pdf'))
    plt.close()
    
    print("\nGenes per cluster:")
    print(genes_per_cluster)
    
    # Violin plot version
    fig, ax = plt.subplots(figsize=(14, 6))
    sc.pl.violin(combined, 'n_genes', groupby='leiden', 
                ax=ax, show=False, rotation=0)
    ax.set_title('Gene count distribution by cluster')
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'genes_per_cluster_violin_{sample}_full.pdf'))
    plt.close()
    
    # 4. Titer by cluster - boxplot with leiden colors
    if 'wolbachia_titer' in combined.obs.columns:
        # UMAP colored by titer
        sc.pl.umap(combined, color='wolbachia_titer', 
                   save=f'_{sample}_full_wolbachia_titer.pdf',
                   title='All samples - Wolbachia titer',
                   vmax=np.percentile(combined.obs['wolbachia_titer'].dropna(), 95))
        
        # Boxplot with stripplot
        fig, ax = plt.subplots(figsize=(14, 6))
        
        plot_data = combined.obs[['leiden', 'wolbachia_titer']].copy()
        plot_data = plot_data.sort_values('leiden')
        
        clusters = sorted(plot_data['leiden'].unique())
        
        # Plot individual points first
        sns.stripplot(data=plot_data, x='leiden', y='wolbachia_titer', 
                     color='black', alpha=0.2, size=1, ax=ax)
        
        # Box plot with leiden colors
        box_parts = ax.boxplot([plot_data[plot_data['leiden'] == cluster]['wolbachia_titer'].values 
                                for cluster in clusters],
                                positions=range(len(clusters)),
                                widths=0.6,
                                patch_artist=True,
                                whiskerprops=dict(alpha=0.7),
                                capprops=dict(alpha=0.7),
                                medianprops=dict(color='black', linewidth=2))
        
        # Color each box
        for patch, color in zip(box_parts['boxes'], leiden_colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        
        ax.set_xlabel('Leiden Cluster')
        ax.set_ylabel('Wolbachia Titer')
        ax.set_title('Wolbachia titer by cluster (all samples)')
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f'boxplot_{sample}_full_titer_by_cluster.pdf'))
        plt.close()
        
        # Summary stats
        titer_stats = combined.obs.groupby('leiden')['wolbachia_titer'].agg(['mean', 'median', 'std'])
        print("\nWolbachia titer by cluster:")
        print(titer_stats)
        
        # Percentage of infected cells per cluster
        infected_per_cluster = combined.obs.groupby('leiden').apply(
            lambda x: (x['wolbachia_titer'] > 0).sum() / len(x) * 100
        )
        
        fig, ax = plt.subplots(figsize=(12, 6))
        ax.bar(range(len(infected_per_cluster)), infected_per_cluster.values,
               color=leiden_colors, alpha=0.7)
        ax.set_xlabel('Leiden Cluster')
        ax.set_ylabel('% Infected cells')
        ax.set_title('Percentage of Wolbachia-infected cells by cluster')
        ax.set_xticks(range(len(infected_per_cluster)))
        ax.set_xticklabels(infected_per_cluster.index)
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f'infection_percentage_by_cluster_{sample}_full.pdf'))
        plt.close()
        
        print("\nPercentage of infected cells per cluster:")
        print(infected_per_cluster)
    
    # 5. Cluster by treatment and cell line
    sc.pl.umap(combined, color=['cell_line', 'treatment'], 
               save=f'_{sample}_full_by_condition.pdf',
               ncols=2)
    
    # 6. Cluster composition by biological condition
    fig, ax = plt.subplots(figsize=(14, 8))
    
    comp_by_condition = pd.crosstab(combined.obs['leiden'], 
                                    combined.obs['bio_condition'], 
                                    normalize='columns') * 100
    
    comp_by_condition.T.plot(kind='bar', stacked=True, ax=ax, 
                             color=leiden_colors, width=0.8)
    ax.set_xlabel('Biological Condition')
    ax.set_ylabel('Percentage of cells')
    ax.set_title('Cluster composition by biological condition')
    ax.legend(title='Leiden Cluster', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'cluster_composition_by_condition_{sample}_full.pdf'))
    plt.close()
    
    # 7. Cell cycle vs Wolbachia infection
    if 'wolbachia_titer' in combined.obs.columns and 'phase' in combined.obs.columns:
        fig, ax = plt.subplots(figsize=(10, 6))
        
        sns.boxplot(data=combined.obs, x='phase', y='wolbachia_titer', ax=ax)
        sns.stripplot(data=combined.obs, x='phase', y='wolbachia_titer',
                     ax=ax, color='black', alpha=0.2, size=1)
        
        ax.set_xlabel('Cell Cycle Phase')
        ax.set_ylabel('Wolbachia Titer')
        ax.set_title('Wolbachia titer by cell cycle phase')
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f'titer_by_cellcycle_{sample}_full.pdf'))
        plt.close()
        
        # Stats by phase
        titer_by_phase = combined.obs.groupby('phase')['wolbachia_titer'].agg(['mean', 'median', 'std', 'count'])
        print("\nWolbachia titer by cell cycle phase:")
        print(titer_by_phase)
    
    # 8. Cells per cluster
    cells_per_cluster = combined.obs['leiden'].value_counts().sort_index()
    
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.bar(range(len(cells_per_cluster)), cells_per_cluster.values,
           color=leiden_colors, alpha=0.7)
    ax.set_xlabel('Leiden Cluster')
    ax.set_ylabel('Number of cells')
    ax.set_title('Cell count by cluster (all samples)')
    ax.set_xticks(range(len(cells_per_cluster)))
    ax.set_xticklabels(cells_per_cluster.index)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'cells_per_cluster_{sample}_full.pdf'))
    plt.close()
    
    print("\nCells per cluster:")
    print(cells_per_cluster)


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
    parser.add_argument('--resolution', type=float, default=0.5,
                        help='Leiden clustering resolution')

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
        leiden_resolution=args.resolution,
    )

if __name__ == "__main__":
    main()