#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Wolbachia Analysis Script
-------------------------
This script performs analysis on single-cell RNA-seq data with Wolbachia infection:
1. Generates UMAP visualizations of log1p_wolbachia_titer
2. Performs differential gene expression analysis between leiden clusters
3. Calculates and visualizes gene coverage statistics
4. Analyzes genes per cell and other quality metrics

Usage:
    python /private/groups/russelllab/jodie/scRNAseq/scRNAseq_methods_comparison/analysis/summary_statistics.py <input_h5ad_path> <output_directory>

Author: Claude
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import scipy.sparse

def main():
    # Check arguments
    if len(sys.argv) != 3:
        print(f"Usage: python {sys.argv[0]} <input_h5ad_path> <output_directory>")
        sys.exit(1)
    
    input_path = sys.argv[1]
    output_dir = sys.argv[2]
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Extract filename without extension for naming outputs
    filename_base = os.path.basename(input_path).split('.')[0]
    
    # Load the AnnData object
    print(f"Loading AnnData from {input_path}...")
    adata = sc.read(input_path)
    print(f"Loaded AnnData with {adata.n_obs} cells and {adata.n_vars} genes")
    
    # Ensure log1p_wolbachia_titer is present
    if 'wolbachia_titer' in adata.obs and 'log1p_wolbachia_titer' not in adata.obs:
        print("Adding log1p transformed Wolbachia titer...")
        adata.obs['log1p_wolbachia_titer'] = np.log1p(adata.obs['wolbachia_titer'])
    
    # 1. Generate UMAP of log1p_wolbachia_titer
    print("\n--- Generating UMAP visualizations ---")
    generate_umap_visualizations(adata, output_dir, filename_base)
    
    # 2. Perform differential gene expression analysis between leiden clusters
    print("\n--- Performing differential gene expression analysis ---")
    perform_differential_expression(adata, output_dir, filename_base)

    # 3. After perform_differential_expression function call
    print("\n--- Calculating cluster similarity ---")
    calculate_cluster_similarity(adata, output_dir, filename_base)
    
    # 4. Analyze genes per cell and quality metrics
    print("\n--- Analyzing genes per cell ---")
    analyze_genes_per_cell(adata, output_dir, filename_base)
    
    # 5. Calculate gene coverage across cells
    print("\n--- Calculating gene coverage across cells ---")
    analyze_gene_coverage(adata, output_dir, filename_base)
    
    print(f"\nAll analyses completed. Results saved to {output_dir}")

def generate_umap_visualizations(adata, output_dir, filename_base):
    """Generate UMAP visualizations colored by log1p_wolbachia_titer"""
    
    # Check if UMAP exists in the object
    if 'X_umap' not in adata.obsm:
        print("Warning: UMAP not found in the AnnData object. Computing UMAP...")
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)
    
    # Check if log1p_wolbachia_titer exists
    if 'log1p_wolbachia_titer' not in adata.obs:
        print("Warning: log1p_wolbachia_titer not found. Please ensure Wolbachia titer data is available.")
        return
    
    # Generate UMAP colored by log1p_wolbachia_titer
    print("Generating UMAP colored by log1p_wolbachia_titer...")
    
    # Using scanpy's plotting function
    sc.settings.figdir = output_dir
    sc.pl.umap(adata, color='log1p_wolbachia_titer', 
               save=f'_{filename_base}_log1p_wolbachia_titer.pdf',
               show=False)
    plt.close()
    print(f"UMAP visualizations saved to {output_dir}")

def perform_differential_expression(adata, output_dir, filename_base):
    """Perform differential gene expression analysis between leiden clusters"""
    
    # Check if leiden clustering exists
    if 'leiden' not in adata.obs:
        print("Warning: leiden clustering not found. Computing leiden clusters...")
        if 'neighbors' not in adata.uns:
            sc.pp.neighbors(adata)
        sc.tl.leiden(adata)
    
    # Perform rank_genes_groups analysis
    print("Performing differential expression analysis using Wilcoxon rank-sum test...")
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon', n_genes=50)
    
    # Generate plot
    sc.settings.figdir = output_dir
    sc.pl.rank_genes_groups(adata, n_genes=50, 
                           save=f'_{filename_base}_rank_genes.pdf',
                           show=False)
    
    # Export results to CSV
    print("Exporting differential expression results to CSV...")
    
    # Extract the data from the rank_genes dictionary
    rank_genes = adata.uns['rank_genes_groups']
    groups = adata.obs['leiden'].cat.categories.tolist()  # Get cluster names
    
    # Create a DataFrame for storing results
    df = pd.DataFrame()
    
    # Add columns for each statistic
    for stat in ['names', 'scores', 'pvals', 'pvals_adj', 'logfoldchanges']:
        if stat in rank_genes:
            stat_df = pd.DataFrame(rank_genes[stat])
            stat_df.columns = [f'{group}_{stat}' for group in groups]
            df = pd.concat([df, stat_df], axis=1)
    
    # Save to CSV
    de_results_file = os.path.join(output_dir, f"{filename_base}_ranked_genes.csv")
    df.to_csv(de_results_file, index=False)
    print(f"Differential expression results saved to {de_results_file}")

def calculate_cluster_similarity(adata, output_dir, filename_base):
    """
    Calculate Jaccard similarity index between leiden clusters based on marker genes
    """
    print("Calculating Jaccard similarity index between clusters...")
    
    # Check if leiden and rank_genes_groups are available
    if 'leiden' not in adata.obs or 'rank_genes_groups' not in adata.uns:
        print("Warning: leiden clusters or rank_genes_groups not found. Skipping Jaccard similarity calculation.")
        return
    
    # Get cluster information
    clusters = adata.obs['leiden'].cat.categories.tolist()
    n_clusters = len(clusters)
    
    # Get the top marker genes for each cluster
    rank_genes = adata.uns['rank_genes_groups']
    n_top_genes = 50  # Number of top genes to consider for each cluster (or all available)
    
    # Extract the top genes for each cluster
    markers_by_cluster = {}
    
    # Looking at your specific error, rank_genes['names'] is 1D with shape (50,)
    # This likely means it's a recarray where each element contains a tuple of genes for each cluster
    if 'names' in rank_genes:
        print(f"Shape of rank_genes['names']: {rank_genes['names'].shape}")
        print(f"Type of rank_genes['names']: {type(rank_genes['names'])}")
        
        # For structured arrays/recarrays
        if hasattr(rank_genes['names'], 'dtype') and rank_genes['names'].dtype.names is not None:
            print("Detected recarray structure")
            # Each row in the array contains values for all clusters
            for i, row in enumerate(rank_genes['names']):
                if i >= n_top_genes:
                    break
                # Each element of row contains the gene for a specific cluster
                for j, cluster in enumerate(clusters):
                    if j >= len(row):
                        continue
                    if cluster not in markers_by_cluster:
                        markers_by_cluster[cluster] = set()
                    gene_name = row[j]
                    if gene_name:  # Ensure it's not empty
                        markers_by_cluster[cluster].add(gene_name)
        else:
            print("Unknown structure for rank_genes['names']. Using alternative approach.")
            # Last resort
            for cluster in clusters:
                markers_by_cluster[cluster] = set()
            
            # Try to find other ways to get marker genes
            if 'scores' in rank_genes and hasattr(rank_genes['scores'], 'dtype') and rank_genes['scores'].dtype.names is not None:
                print("Using scores to identify top genes")
                # Use the score array which might have the same structure
                for i, row in enumerate(rank_genes['scores']):
                    if i >= n_top_genes:
                        break
                    for j, cluster in enumerate(clusters):
                        if j >= len(row):
                            continue
                        if rank_genes['names'][i] and rank_genes['scores'][i, j] > 0:
                            markers_by_cluster[cluster].add(rank_genes['names'][i])
    
    # If no markers could be extracted, use a different approach
    if not any(markers_by_cluster.values()):
        print("Warning: Could not extract marker genes using standard methods. Falling back to scanpy's get_rank_genes().")
        # Fallback: Use scanpy to get rank genes directly for each group
        for cluster in clusters:
            # Get top genes for this cluster if available
            try:
                sc_markers = sc.get.rank_genes_groups_df(adata, group=cluster).head(n_top_genes)
                markers_by_cluster[cluster] = set(sc_markers['names'].values)
            except:
                print(f"Could not get markers for cluster {cluster}. Using empty set.")
                markers_by_cluster[cluster] = set()
    
    # If still no markers could be extracted, exit
    if not any(markers_by_cluster.values()):
        print("Warning: Could not extract any marker genes. Skipping Jaccard similarity calculation.")
        return
    
    # Print the number of marker genes found for each cluster
    print("Number of marker genes found per cluster:")
    for cluster, markers in markers_by_cluster.items():
        print(f"  Cluster {cluster}: {len(markers)} genes")
    
    # Calculate Jaccard similarity matrix
    similarity_matrix = np.zeros((n_clusters, n_clusters))
    
    for i, cluster1 in enumerate(clusters):
        for j, cluster2 in enumerate(clusters):
            genes1 = markers_by_cluster[cluster1]
            genes2 = markers_by_cluster[cluster2]
            
            # Calculate Jaccard similarity: size of intersection / size of union
            intersection = len(genes1.intersection(genes2))
            union = len(genes1.union(genes2))
            
            if union > 0:
                jaccard = intersection / union
            else:
                jaccard = 0
                
            similarity_matrix[i, j] = jaccard
    
    # Create a DataFrame for easy visualization
    similarity_df = pd.DataFrame(similarity_matrix, 
                                index=clusters, 
                                columns=clusters)
    
    # Save to CSV
    similarity_file = os.path.join(output_dir, f"{filename_base}_cluster_jaccard_similarity.csv")
    similarity_df.to_csv(similarity_file)
    
    # Create a heatmap visualization
    plt.figure(figsize=(10, 8))
    sns.heatmap(similarity_df, annot=True, cmap="YlGnBu", vmin=0, vmax=1)
    plt.title('Jaccard Similarity Between Clusters Based on Marker Genes')
    plt.tight_layout()
    
    # Save the heatmap
    heatmap_file = os.path.join(output_dir, f"{filename_base}_cluster_jaccard_similarity_heatmap.pdf")
    plt.savefig(heatmap_file)
    
    # Also save as PNG for easy viewing
    heatmap_file_png = os.path.join(output_dir, f"{filename_base}_cluster_jaccard_similarity_heatmap.png")
    plt.savefig(heatmap_file_png)
    plt.close()
    
    print(f"Jaccard similarity matrix saved to {similarity_file}")
    print(f"Jaccard similarity heatmap saved to {heatmap_file}")
    
    # Additionally, create a cluster dendrogram based on Jaccard distances
    from scipy.cluster import hierarchy
    from scipy.spatial.distance import squareform
    
    # Convert similarity to distance (1 - similarity)
    jaccard_distance = 1 - similarity_matrix
    
    # Perform hierarchical clustering
    linkage = hierarchy.linkage(squareform(jaccard_distance), method='average')
    
    # Create dendrogram plot
    plt.figure(figsize=(12, 6))
    hierarchy.dendrogram(linkage, labels=clusters, leaf_rotation=90)
    plt.title('Cluster Dendrogram Based on Jaccard Distance of Marker Genes')
    plt.xlabel('Cluster')
    plt.ylabel('Jaccard Distance')
    plt.tight_layout()
    
    # Save the dendrogram
    dendrogram_file = os.path.join(output_dir, f"{filename_base}_cluster_dendrogram.pdf")
    plt.savefig(dendrogram_file)
    plt.close()
    
    print(f"Cluster dendrogram saved to {dendrogram_file}")

def analyze_genes_per_cell(adata, output_dir, filename_base):
    """Analyze the number of genes expressed per cell"""
    
    print("Calculating genes per cell...")
    
    # Count number of expressed genes per cell (non-zero counts)
    if 'n_genes' not in adata.obs:
        genes_per_cell = (adata.X > 0).sum(axis=1)
        
        # If adata.X is a sparse matrix, convert the result to a dense array
        if scipy.sparse.issparse(adata.X):
            genes_per_cell = genes_per_cell.A1  # Convert to 1D array
        
        # Add this information to adata.obs
        adata.obs['n_genes'] = genes_per_cell
    
    # Get summary statistics
    stats = {
        'Mean genes per cell': adata.obs['n_genes'].mean(),
        'Median genes per cell': adata.obs['n_genes'].median(),
        'Min genes per cell': adata.obs['n_genes'].min(),
        'Max genes per cell': adata.obs['n_genes'].max(),
        'Total cells': len(adata.obs)
    }
    
    # Save statistics to file
    stats_file = os.path.join(output_dir, f"{filename_base}_genes_per_cell_stats.csv")
    stats_df = pd.DataFrame([stats])
    stats_df.to_csv(stats_file, index=False)
    
    # Print statistics
    for key, value in stats.items():
        if key in ['Mean genes per cell', 'Median genes per cell']:
            print(f"{key}: {value:.2f}")
        else:
            print(f"{key}: {value:.0f}")
    
    # Visualize distribution
    plt.figure(figsize=(10, 6))
    sns.histplot(adata.obs['n_genes'], bins=50)
    plt.xlabel('Number of expressed genes')
    plt.ylabel('Number of cells')
    plt.title('Distribution of gene expression per cell')
    plt.grid(alpha=0.3)
    plt.tight_layout()
    
    # Save the plot
    genes_per_cell_plot = os.path.join(output_dir, f"{filename_base}_genes_per_cell_distribution.pdf")
    plt.savefig(genes_per_cell_plot)
    
    # Also save as PNG for easy viewing
    genes_per_cell_plot_png = os.path.join(output_dir, f"{filename_base}_genes_per_cell_distribution.png")
    plt.savefig(genes_per_cell_plot_png)
    plt.close()
    
    print(f"Genes per cell analysis saved to {stats_file} and {genes_per_cell_plot}")
    
    # Create a violin plot of genes per cell grouped by leiden cluster
    if 'leiden' in adata.obs:
        plt.figure(figsize=(12, 6))
        sns.violinplot(x='leiden', y='n_genes', data=adata.obs)
        plt.xlabel('Leiden Cluster')
        plt.ylabel('Number of expressed genes')
        plt.title('Genes per cell by cluster')
        plt.xticks(rotation=45)
        plt.tight_layout()
        
        # Save the plot
        genes_by_cluster_plot = os.path.join(output_dir, f"{filename_base}_genes_per_cell_by_cluster.pdf")
        plt.savefig(genes_by_cluster_plot)
        plt.close()
        print(f"Genes per cell by cluster plot saved to {genes_by_cluster_plot}")

def analyze_gene_coverage(adata, output_dir, filename_base):
    """Analyze gene coverage across cells"""
    
    print("Calculating gene coverage across cells...")
    
    # Calculate the percentage of cells expressing each gene
    # Get binary matrix of expression (1 if gene expressed, 0 if not)
    if scipy.sparse.issparse(adata.X):
        binary_matrix = (adata.X > 0).astype(int)
        # Keep it sparse for efficiency
        gene_coverage = binary_matrix.mean(axis=0).A1 * 100  # Get percentage
    else:
        binary_matrix = (adata.X > 0).astype(int)
        gene_coverage = binary_matrix.mean(axis=0) * 100
    
    # Create a DataFrame with gene names and coverage percentages
    coverage_df = pd.DataFrame({
        'gene': adata.var_names,
        'pct_cells': gene_coverage
    })
    
    # Calculate mean expression of each gene across cells
    if scipy.sparse.issparse(adata.X):
        mean_expression = adata.X.mean(axis=0).A1
    else:
        mean_expression = adata.X.mean(axis=0)
    
    # Add mean expression to the DataFrame
    coverage_df['mean_expression'] = mean_expression
    
    # Sort by coverage
    coverage_df = coverage_df.sort_values('pct_cells', ascending=False)
    
    # Add this information to adata.var if not already present
    if 'pct_cells' not in adata.var:
        adata.var['pct_cells'] = gene_coverage
    
    if 'mean_expression' not in adata.var:
        adata.var['mean_expression'] = mean_expression
    
    # Save the data to CSV
    coverage_file = os.path.join(output_dir, f"{filename_base}_gene_coverage.csv")
    coverage_df.to_csv(coverage_file, index=False)
    
    # Save basic stats
    top_genes_file = os.path.join(output_dir, f"{filename_base}_top30_expressed_genes.csv")
    coverage_df.head(30).to_csv(top_genes_file, index=False)
    
    # Visualize the distribution
    plt.figure(figsize=(10, 6))
    sns.histplot(coverage_df['pct_cells'], bins=50)
    plt.xlabel('Percentage of cells expressing gene')
    plt.ylabel('Number of genes')
    plt.title('Distribution of gene coverage across cells')
    plt.grid(alpha=0.3)
    plt.tight_layout()
    
    # Save the plot
    coverage_hist_plot = os.path.join(output_dir, f"{filename_base}_gene_coverage_distribution.pdf")
    plt.savefig(coverage_hist_plot)
    plt.close()
    
    # Plot the top expressed genes
    plt.figure(figsize=(12, 8))
    top_genes = coverage_df.head(30)
    sns.barplot(x='pct_cells', y='gene', data=top_genes)
    plt.xlabel('Percentage of cells expressing gene')
    plt.ylabel('Gene')
    plt.title('Top 30 genes by expression across cells')
    plt.tight_layout()
    
    # Save the plot
    top_genes_plot = os.path.join(output_dir, f"{filename_base}_top30_genes.pdf")
    plt.savefig(top_genes_plot)
    
    # Also save as PNG for easy viewing
    top_genes_plot_png = os.path.join(output_dir, f"{filename_base}_top30_genes.png")
    plt.savefig(top_genes_plot_png)
    plt.close()
    
    print(f"Gene coverage analysis saved to {coverage_file} and {top_genes_file}")
    print(f"Gene coverage plots saved to {coverage_hist_plot} and {top_genes_plot}")
    
    # Create a scatter plot of gene coverage vs mean expression
    plt.figure(figsize=(10, 8))
    plt.scatter(coverage_df['pct_cells'], 
               coverage_df['mean_expression'], 
               alpha=0.5, 
               s=3)
    
    # Annotate some top genes
    for _, row in coverage_df.head(15).iterrows():
        plt.annotate(row['gene'], 
                    (row['pct_cells'], row['mean_expression']),
                    fontsize=8)
    
    plt.xlabel('Percentage of cells expressing gene')
    plt.ylabel('Mean expression level')
    plt.title('Gene coverage vs expression level')
    plt.grid(alpha=0.3)
    plt.tight_layout()
    
    # Save the plot
    scatter_plot = os.path.join(output_dir, f"{filename_base}_gene_coverage_vs_expression.pdf")
    plt.savefig(scatter_plot)
    plt.close()
    print(f"Gene coverage vs expression plot saved to {scatter_plot}")

if __name__ == "__main__":
    main()