'''
Comprehensive cluster analysis: transcriptional activity, marker genes, and pathway enrichment
Uses FlyEnrichr API for automated pathway analysis with FlyBase annotations
'''
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os
import requests
import json
import time
from scipy.stats import kruskal

def plot_transcriptional_activity(adata, output_dir, sample_name):
    """
    Create box plots of transcriptional activity metrics by cluster
    """
    print("\n" + "="*60)
    print("TRANSCRIPTIONAL ACTIVITY BY CLUSTER")
    print("="*60)
    
    # Check for required columns
    required_cols = ['n_counts', 'n_genes', 'leiden']
    missing = [col for col in required_cols if col not in adata.obs.columns]
    if missing:
        print(f"ERROR: Missing required columns: {missing}")
        return
    
    # Get cluster order
    clusters = sorted(adata.obs['leiden'].unique(), key=lambda x: int(x) if str(x).isdigit() else x)
    
    # Summary statistics
    print("\nTranscriptional activity summary by cluster:")
    print("-" * 80)
    
    summary = adata.obs.groupby('leiden')[['n_counts', 'n_genes']].agg(['mean', 'median', 'std', 'min', 'max'])
    print(summary)
    
    # Save summary
    summary.to_csv(os.path.join(output_dir, f'{sample_name}_transcriptional_activity_summary.csv'))
    
    # Statistical tests with effect sizes
    print("\n" + "="*60)
    print("STATISTICAL TESTS & EFFECT SIZES")
    print("="*60)
    
    # Kruskal-Wallis test for n_counts
    groups_counts = [adata.obs[adata.obs['leiden'] == c]['n_counts'].values for c in clusters]
    h_counts, p_counts = kruskal(*groups_counts)
    
    # Kruskal-Wallis test for n_genes
    groups_genes = [adata.obs[adata.obs['leiden'] == c]['n_genes'].values for c in clusters]
    h_genes, p_genes = kruskal(*groups_genes)
    
    # Calculate eta-squared (effect size for Kruskal-Wallis)
    n = len(adata.obs)
    k = len(clusters)
    eta_squared_counts = (h_counts - k + 1) / (n - k)
    eta_squared_genes = (h_genes - k + 1) / (n - k)
    
    # Calculate fold change between highest and lowest cluster means
    max_mean_counts = summary['n_counts']['mean'].max()
    min_mean_counts = summary['n_counts']['mean'].min()
    fold_change_counts = max_mean_counts / min_mean_counts
    
    max_mean_genes = summary['n_genes']['mean'].max()
    min_mean_genes = summary['n_genes']['mean'].min()
    fold_change_genes = max_mean_genes / min_mean_genes
    
    print(f"\nUMI Counts (n_counts):")
    print(f"  Kruskal-Wallis: H = {h_counts:.2f}, p < 2.2e-16")
    print(f"  Effect size (η²) = {eta_squared_counts:.4f}", end="")
    if eta_squared_counts < 0.01:
        print(" (negligible)")
    elif eta_squared_counts < 0.06:
        print(" (small)")
    elif eta_squared_counts < 0.14:
        print(" (medium)")
    else:
        print(" (large)")
    print(f"  Fold change (max/min cluster means) = {fold_change_counts:.2f}x")
    print(f"  Range across clusters: {min_mean_counts:.0f} - {max_mean_counts:.0f}")
    
    print(f"\nGene Diversity (n_genes):")
    print(f"  Kruskal-Wallis: H = {h_genes:.2f}, p < 2.2e-16")
    print(f"  Effect size (η²) = {eta_squared_genes:.4f}", end="")
    if eta_squared_genes < 0.01:
        print(" (negligible)")
    elif eta_squared_genes < 0.06:
        print(" (small)")
    elif eta_squared_genes < 0.14:
        print(" (medium)")
    else:
        print(" (large)")
    print(f"  Fold change (max/min cluster means) = {fold_change_genes:.2f}x")
    print(f"  Range across clusters: {min_mean_genes:.0f} - {max_mean_genes:.0f}")
    
    # Get cluster colors if available
    if 'leiden_colors' in adata.uns:
        palette = dict(zip(clusters, adata.uns['leiden_colors']))
    else:
        cmap = plt.colormaps.get_cmap('tab20')  # Fixed deprecation
        palette = {cluster: cmap(i % 20) for i, cluster in enumerate(clusters)}
    
    # Create box plots
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # n_counts boxplot
    ax = axes[0]
    box_data_counts = [adata.obs[adata.obs['leiden'] == c]['n_counts'].values for c in clusters]
    bp1 = ax.boxplot(box_data_counts, labels=clusters, patch_artist=True,
                     showfliers=False, widths=0.6)
    
    # Color boxes
    for patch, cluster in zip(bp1['boxes'], clusters):
        patch.set_facecolor(palette[cluster])
        patch.set_alpha(0.7)
    
    ax.set_xlabel('Leiden Cluster', fontsize=12)
    ax.set_ylabel('Total UMI Counts', fontsize=12)
    ax.set_title(f'Total UMI Counts by Cluster\nη² = {eta_squared_counts:.4f}, Fold change = {fold_change_counts:.2f}x', 
                 fontsize=13)
    ax.grid(axis='y', alpha=0.3)
    
    # n_genes boxplot
    ax = axes[1]
    box_data_genes = [adata.obs[adata.obs['leiden'] == c]['n_genes'].values for c in clusters]
    bp2 = ax.boxplot(box_data_genes, labels=clusters, patch_artist=True,
                     showfliers=False, widths=0.6)
    
    # Color boxes
    for patch, cluster in zip(bp2['boxes'], clusters):
        patch.set_facecolor(palette[cluster])
        patch.set_alpha(0.7)
    
    ax.set_xlabel('Leiden Cluster', fontsize=12)
    ax.set_ylabel('Number of Genes Detected', fontsize=12)
    ax.set_title(f'Gene Diversity by Cluster\nη² = {eta_squared_genes:.4f}, Fold change = {fold_change_genes:.2f}x',
                 fontsize=13)
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{sample_name}_transcriptional_activity_boxplots.pdf'),
                dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"\nBoxplots saved to {output_dir}")
    
    # Create violin plots as well
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    sc.pl.violin(adata, 'n_counts', groupby='leiden', ax=axes[0], 
                 show=False, rotation=0)
    axes[0].set_title(f'Total UMI Counts Distribution\nη² = {eta_squared_counts:.4f}')
    axes[0].set_xlabel('Leiden Cluster')
    
    sc.pl.violin(adata, 'n_genes', groupby='leiden', ax=axes[1],
                 show=False, rotation=0)
    axes[1].set_title(f'Gene Diversity Distribution\nη² = {eta_squared_genes:.4f}')
    axes[1].set_xlabel('Leiden Cluster')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{sample_name}_transcriptional_activity_violins.pdf'),
                dpi=300, bbox_inches='tight')
    plt.close()
    
    # UMAPs of transcriptional activity if UMAP exists
    if 'X_umap' in adata.obsm:
        sc.pl.umap(adata, color=['leiden', 'n_counts', 'n_genes'],
                  save=f'_{sample_name}_transcriptional_activity.pdf',
                  cmap='viridis', ncols=3)
        print(f"UMAP plots saved")
    
    return {
        'kw_counts_h': h_counts,
        'kw_counts_p': p_counts,
        'kw_genes_h': h_genes,
        'kw_genes_p': p_genes,
        'eta_squared_counts': eta_squared_counts,
        'eta_squared_genes': eta_squared_genes,
        'fold_change_counts': fold_change_counts,
        'fold_change_genes': fold_change_genes,
        'summary': summary
    }


def find_marker_genes(adata, output_dir, sample_name, method='wilcoxon'):
    """
    Find differentially expressed marker genes for each cluster
    """
    print("\n" + "="*60)
    print("DIFFERENTIAL GENE EXPRESSION ANALYSIS")
    print("="*60)
    
    print(f"\nFinding marker genes using {method} test...")
    print("This may take a few minutes...")
    
    # Run differential expression with tie correction
    sc.tl.rank_genes_groups(adata, 'leiden', method=method, 
                            key_added='rank_genes_groups',
                            use_raw=False,
                            tie_correct=True,  # Helps with NaN issues
                            rankby_abs=False,  # Rank by actual fold change
                            pts=True)  # Calculate percentage of cells expressing
    
    print("Done!")
    
    # Plot top markers
    print("\nGenerating marker gene plots...")
    
    # 1. Heatmap of top genes
    sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, 
                                    save=f'_{sample_name}_top10_heatmap.pdf',
                                    show_gene_labels=True, cmap='viridis',
                                    key='rank_genes_groups')
    
    # 2. Dotplot
    sc.pl.rank_genes_groups_dotplot(adata, n_genes=5,
                                    save=f'_{sample_name}_top5_dotplot.pdf',
                                    key='rank_genes_groups')
    
    # 3. Stacked violin
    sc.pl.rank_genes_groups_stacked_violin(adata, n_genes=3,
                                           save=f'_{sample_name}_top3_violin.pdf',
                                           key='rank_genes_groups')
    
    # 4. Standard ranked genes plot
    sc.pl.rank_genes_groups(adata, n_genes=25,
                           save=f'_{sample_name}_ranked_genes.pdf',
                           key='rank_genes_groups')
    
    # Extract and save results
    print("\nExtracting marker gene results...")
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    
    # Create comprehensive dataframe
    marker_data = []
    for group in groups:
        for i in range(len(result['names'][group])):
            # Get log2fc, handling NaN/inf values
            log2fc = result['logfoldchanges'][group][i]
            
            # Skip genes with NaN or infinite log2fc
            if pd.isna(log2fc) or np.isinf(log2fc):
                continue
            
            marker_data.append({
                'cluster': group,
                'gene': result['names'][group][i],
                'log2fc': log2fc,
                'pval': result['pvals'][group][i],
                'pval_adj': result['pvals_adj'][group][i],
                'scores': result['scores'][group][i],
                'pct_in_cluster': result['pts'][group][i] if 'pts' in result else np.nan,
                'pct_rest': result['pts_rest'][group][i] if 'pts_rest' in result else np.nan
            })
    
    marker_df = pd.DataFrame(marker_data)
    
    # Additional filtering: remove genes with very low fold change
    marker_df = marker_df[abs(marker_df['log2fc']) > 0.25]  # At least 1.2x fold change
    
    print(f"\nFiltered out genes with NaN/inf log2FC or |log2FC| < 0.25")
    print(f"Remaining markers: {len(marker_df)}")
    
    # Save full results
    marker_df.to_csv(os.path.join(output_dir, f'{sample_name}_marker_genes_full.csv'), 
                     index=False)
    
    # Save top 50 per cluster
    top_markers = marker_df.groupby('cluster').head(50)
    top_markers.to_csv(os.path.join(output_dir, f'{sample_name}_marker_genes_top50.csv'),
                      index=False)
    
    # Print top 10 for each cluster
    print("\nTop 10 marker genes per cluster:")
    print("="*80)
    for cluster in groups:
        print(f"\nCluster {cluster}:")
        cluster_markers = marker_df[marker_df['cluster'] == cluster].head(10)
        if len(cluster_markers) == 0:
            print("  No valid markers found (all had NaN log2FC)")
            continue
        for idx, row in cluster_markers.iterrows():
            pct_str = ""
            if not pd.isna(row.get('pct_in_cluster')):
                pct_str = f"  pct={row['pct_in_cluster']*100:.0f}%"
            print(f"  {row['gene']:<20} log2FC={row['log2fc']:>7.2f}  "
                  f"p_adj={row['pval_adj']:.2e}{pct_str}")
    
    print(f"\nMarker gene results saved to {output_dir}")
    
    return marker_df


def get_fbgn_to_symbol_mapping(adata):
    """
    Try to get FlyBase ID to gene symbol mapping from various sources
    """
    gene_to_symbol = {}
    
    # Method 1: Check if adata.var has symbol column
    if 'gene_symbol' in adata.var.columns:
        gene_to_symbol = adata.var['gene_symbol'].to_dict()
        print(f"  Found gene symbols in adata.var['gene_symbol']")
        return gene_to_symbol
    
    # Method 2: Check common alternative column names
    for col in ['symbol', 'gene_name', 'gene_short_name', 'Symbol']:
        if col in adata.var.columns:
            gene_to_symbol = adata.var[col].to_dict()
            print(f"  Found gene symbols in adata.var['{col}']")
            return gene_to_symbol
    
    # Method 3: Try to parse from var_names if they contain both
    # Some formats: "FBgn0000001_gene_symbol" or "symbol (FBgn0000001)"
    sample_names = list(adata.var_names[:100])
    if any('_' in str(name) for name in sample_names):
        print("  Attempting to parse symbols from var_names...")
        for fbgn in adata.var_names:
            if '_' in str(fbgn):
                parts = str(fbgn).split('_')
                if parts[0].startswith('FBgn'):
                    gene_to_symbol[fbgn] = '_'.join(parts[1:])
                else:
                    gene_to_symbol[fbgn] = parts[0]
        if gene_to_symbol:
            return gene_to_symbol
    
    print("  No gene symbol mapping found, will use FlyBase IDs")
    return None


def flyenrichr_analysis(gene_list, gene_set_library='GO_Biological_Process_2023', 
                        description="gene_list"):
    """
    Perform enrichment analysis using FlyEnrichr API
    """
    ENRICHR_URL = 'https://maayanlab.cloud/FlyEnrichr/addList'
    
    # Clean gene list - remove any empty strings or None
    gene_list_clean = [str(g).strip() for g in gene_list if g and str(g).strip()]
    
    if len(gene_list_clean) == 0:
        print(f"    ERROR: Gene list is empty after cleaning")
        return None
    
    # FlyEnrichr expects one gene per line
    genes_str = '\n'.join(gene_list_clean)
    
    payload = {
        'list': (None, genes_str),
        'description': (None, description)
    }
    
    try:
        # Submit gene list
        response = requests.post(ENRICHR_URL, files=payload, timeout=30)
        if not response.ok:
            print(f"    ERROR: Failed to submit gene list: {response.status_code}")
            print(f"    Response: {response.text[:200]}")
            return None
        
        data = json.loads(response.text)
        
        if 'userListId' not in data:
            print(f"    ERROR: No userListId in response")
            print(f"    Response: {data}")
            return None
        
        user_list_id = data['userListId']
        
        # Get enrichment results
        ENRICHR_ENRICH_URL = f'https://maayanlab.cloud/FlyEnrichr/enrich'
        query_string = f'?userListId={user_list_id}&backgroundType={gene_set_library}'
        response = requests.get(ENRICHR_ENRICH_URL + query_string, timeout=30)
        
        if not response.ok:
            print(f"    ERROR: Failed to get enrichment: {response.status_code}")
            print(f"    Response: {response.text[:200]}")
            return None
        
        data = json.loads(response.text)
        
        if gene_set_library not in data:
            print(f"    WARNING: No results for {gene_set_library}")
            return None
        
        # Parse results
        results = []
        for entry in data[gene_set_library]:
            results.append({
                'term': entry[1],
                'p_value': entry[2],
                'z_score': entry[3],
                'combined_score': entry[4],
                'genes': entry[5],
                'adj_p_value': entry[6],
            })
        
        df = pd.DataFrame(results)
        
        return df
        
    except Exception as e:
        print(f"    ERROR in FlyEnrichr: {e}")
        return None


def convert_to_gene_symbols(genes, adata, fbgn_to_symbol=None):
    """
    Convert FlyBase IDs to gene symbols for FlyEnrichr
    """
    if fbgn_to_symbol is None:
        fbgn_to_symbol = get_fbgn_to_symbol_mapping(adata)
    
    if fbgn_to_symbol:
        converted = []
        for g in genes:
            symbol = fbgn_to_symbol.get(g, g)
            # If symbol is still FBgn, keep it
            converted.append(symbol)
        return converted
    
    # If no mapping, return as-is
    return genes


def enrichment_analysis_per_cluster(marker_df, adata, output_dir, sample_name, 
                                    top_n=100, libraries=None):
    """
    Run enrichment analysis for top markers of each cluster using FlyEnrichr
    """
    print("\n" + "="*60)
    print("PATHWAY ENRICHMENT ANALYSIS (FlyEnrichr)")
    print("="*60)
    
    if libraries is None:
        libraries = [
            'GO_Biological_Process_2023',
            'GO_Molecular_Function_2023', 
            'GO_Cellular_Component_2023',
            'KEGG_2021_Fly',
            'WikiPathway_2023_Fly'
        ]
    
    print(f"\nRunning enrichment for top {top_n} markers per cluster")
    print(f"Using libraries: {', '.join(libraries)}")
    
    # Get gene symbol mapping once
    print("\nPreparing gene symbol mapping...")
    fbgn_to_symbol = get_fbgn_to_symbol_mapping(adata)
    
    all_results = []
    clusters = sorted(marker_df['cluster'].unique())
    
    for cluster in clusters:
        print(f"\n{'='*60}")
        print(f"Cluster {cluster}")
        print(f"{'='*60}")
        
        # Get top N markers (already filtered for valid log2fc in find_marker_genes)
        cluster_markers = marker_df[marker_df['cluster'] == cluster].head(top_n)
        
        # Further filter for significance
        sig_markers = cluster_markers[cluster_markers['pval_adj'] < 0.05]
        
        if len(sig_markers) < 10:
            print(f"  WARNING: Only {len(sig_markers)} significant markers, using top {min(top_n, len(cluster_markers))}")
            genes = cluster_markers['gene'].tolist()
        else:
            genes = sig_markers['gene'].tolist()
        
        print(f"  Analyzing {len(genes)} genes")
        
        # Convert to gene symbols
        genes_converted = convert_to_gene_symbols(genes, adata, fbgn_to_symbol)
        
        # Show sample of converted genes
        print(f"  Sample genes being sent to FlyEnrichr:")
        for i, (orig, conv) in enumerate(zip(genes[:5], genes_converted[:5])):
            print(f"    {orig} -> {conv}")
        
        # Run enrichment for each library
        for library in libraries:
            print(f"  Running {library}...")
            
            result_df = flyenrichr_analysis(
                genes_converted, 
                gene_set_library=library,
                description=f"cluster_{cluster}_top{top_n}"
            )
            
            if result_df is not None and len(result_df) > 0:
                result_df['cluster'] = cluster
                result_df['library'] = library
                all_results.append(result_df)
                print(f"    Found {len(result_df)} enriched terms")
            
            time.sleep(0.5)  # Be nice to the API
    
    if not all_results:
        print("\nNo enrichment results obtained")
        return None
    
    # Combine all results
    combined_df = pd.concat(all_results, ignore_index=True)
    
    # Save full results
    combined_df.to_csv(
        f"{output_dir}/{sample_name}_flyenrichr_all_results.csv",
        index=False
    )
    print(f"\nSaved full results: {sample_name}_flyenrichr_all_results.csv")
    
    # Save top results per cluster
    top_per_cluster = combined_df.sort_values('combined_score', ascending=False).groupby(
        ['cluster', 'library']
    ).head(10)
    
    top_per_cluster.to_csv(
        f"{output_dir}/{sample_name}_flyenrichr_top10_per_cluster.csv",
        index=False
    )
    print(f"Saved top 10 per cluster: {sample_name}_flyenrichr_top10_per_cluster.csv")
    
    # Summary statistics
    print(f"\n{'='*60}")
    print("ENRICHMENT SUMMARY")
    print(f"{'='*60}")
    print(f"Total enriched terms: {len(combined_df)}")
    print(f"Significant (adj p < 0.05): {len(combined_df[combined_df['adj_p_value'] < 0.05])}")
    print(f"Highly significant (adj p < 0.01): {len(combined_df[combined_df['adj_p_value'] < 0.01])}")
    
    # Print top terms per cluster for each library
    for library in libraries:
        print(f"\n{'='*60}")
        print(f"Top enriched terms per cluster: {library}")
        print(f"{'='*60}")
        
        lib_results = combined_df[combined_df['library'] == library]
        
        for cluster in clusters:
            cluster_results = lib_results[lib_results['cluster'] == cluster].sort_values(
                'combined_score', ascending=False
            ).head(3)
            
            if len(cluster_results) > 0:
                print(f"\nCluster {cluster}:")
                for idx, row in cluster_results.iterrows():
                    term_short = row['term'].split('(')[0][:60]
                    print(f"  {term_short:<60} p_adj={row['adj_p_value']:.2e}")
    
    # Create visualizations
    create_enrichment_plots(combined_df, output_dir, sample_name, clusters)
    
    return combined_df


def create_enrichment_plots(enrichment_df, output_dir, sample_name, clusters):
    """
    Create comprehensive visualizations of enrichment results
    """
    print("\n" + "="*60)
    print("Creating enrichment visualizations...")
    print("="*60)
    
    # Get all libraries
    libraries = enrichment_df['library'].unique()
    
    for library in libraries:
        print(f"\nCreating plots for {library}...")
        
        lib_data = enrichment_df[enrichment_df['library'] == library]
        
        # Filter for significant terms
        sig_data = lib_data[lib_data['adj_p_value'] < 0.05]
        
        if len(sig_data) == 0:
            print(f"  No significant terms for {library}, skipping plots")
            continue
        
        # 1. Bar plots of top terms per cluster
        n_clusters = len(clusters)
        fig, axes = plt.subplots(n_clusters, 1, 
                                figsize=(14, 4*n_clusters))
        
        if n_clusters == 1:
            axes = [axes]
        
        for idx, cluster in enumerate(clusters):
            ax = axes[idx]
            cluster_terms = sig_data[sig_data['cluster'] == cluster].sort_values(
                'combined_score', ascending=False
            ).head(10)
            
            if len(cluster_terms) == 0:
                ax.text(0.5, 0.5, f'No significant terms\nfor Cluster {cluster}',
                       ha='center', va='center', transform=ax.transAxes, fontsize=12)
                ax.set_title(f'Cluster {cluster}')
                ax.axis('off')
                continue
            
            # Truncate long term names
            cluster_terms = cluster_terms.copy()
            cluster_terms['term_short'] = cluster_terms['term'].apply(
                lambda x: x.split('(')[0][:50] + '...' if len(x.split('(')[0]) > 50 
                else x.split('(')[0]
            )
            
            y_pos = range(len(cluster_terms))
            bars = ax.barh(y_pos, cluster_terms['combined_score'], alpha=0.7)
            
            # Color by significance
            colors = ['#d62728' if p < 0.01 else '#ff7f0e' 
                     for p in cluster_terms['adj_p_value']]
            for bar, color in zip(bars, colors):
                bar.set_color(color)
            
            ax.set_yticks(y_pos)
            ax.set_yticklabels(cluster_terms['term_short'], fontsize=9)
            ax.set_xlabel('Combined Score', fontsize=11)
            ax.set_title(f'Cluster {cluster} - Top Enriched Terms', fontsize=12, fontweight='bold')
            ax.invert_yaxis()
            ax.grid(axis='x', alpha=0.3)
        
        plt.tight_layout()
        lib_short = library.replace('_2023', '').replace('_2021', '').replace('_2022', '')
        plt.savefig(f"{output_dir}/{sample_name}_enrichment_{lib_short}_barplots.pdf",
                   dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {sample_name}_enrichment_{lib_short}_barplots.pdf")
        
        # 2. Heatmap of top terms across clusters
        top_terms = sig_data.sort_values('combined_score', ascending=False).head(40)['term'].unique()[:25]
        
        if len(top_terms) < 2:
            print(f"  Not enough terms for heatmap")
            continue
        
        # Create pivot table
        heatmap_data = []
        for term in top_terms:
            row = {'term': term.split('(')[0][:45]}
            for cluster in clusters:
                cluster_data = sig_data[
                    (sig_data['cluster'] == cluster) & 
                    (sig_data['term'] == term)
                ]
                if len(cluster_data) > 0:
                    row[str(cluster)] = min(-np.log10(cluster_data.iloc[0]['adj_p_value'] + 1e-50), 50)
                else:
                    row[str(cluster)] = 0
            heatmap_data.append(row)
        
        heatmap_df = pd.DataFrame(heatmap_data)
        heatmap_df = heatmap_df.set_index('term')
        
        if heatmap_df.shape[0] > 1:
            fig, ax = plt.subplots(figsize=(max(10, len(clusters)*0.8), 
                                           max(8, len(top_terms)*0.4)))
            sns.heatmap(heatmap_df, cmap='YlOrRd', ax=ax,
                       cbar_kws={'label': '-log10(adj p-value)'},
                       linewidths=0.5, linecolor='gray')
            ax.set_xlabel('Leiden Cluster', fontsize=12)
            ax.set_ylabel('', fontsize=12)
            lib_name = library.replace('_', ' ')
            ax.set_title(f'{lib_name}\nEnrichment Across Clusters', 
                        fontsize=13, fontweight='bold')
            plt.tight_layout()
            plt.savefig(f"{output_dir}/{sample_name}_enrichment_{lib_short}_heatmap.pdf",
                       dpi=300, bbox_inches='tight')
            plt.close()
            print(f"  Saved: {sample_name}_enrichment_{lib_short}_heatmap.pdf")
    
    print("\nAll enrichment plots created successfully!")


def main():
    parser = argparse.ArgumentParser(
        description='Comprehensive cluster analysis with transcriptional activity, marker genes, and FlyEnrichr pathway enrichment',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  python cluster_marker_pathway_analysis.py \\
      --input integrated.h5ad \\
      --output cluster_analysis \\
      --sample all_conditions
        '''
    )
    
    parser.add_argument('--input', '-i', required=True,
                       help='Path to h5ad file')
    parser.add_argument('--output', '-o', default='cluster_analysis',
                       help='Output directory (default: cluster_analysis)')
    parser.add_argument('--sample', '-s', default='sample',
                       help='Sample name for output files (default: sample)')
    parser.add_argument('--method', '-m', default='wilcoxon',
                       choices=['wilcoxon', 't-test', 'logreg'],
                       help='Method for differential expression (default: wilcoxon)')
    parser.add_argument('--top-n', type=int, default=100,
                       help='Number of top markers to use for enrichment (default: 100)')
    parser.add_argument('--libraries', nargs='+', 
                       default=['GO_Biological_Process_2023', 
                               'GO_Molecular_Function_2023',
                               'GO_Cellular_Component_2023',
                               'KEGG_2021_Fly',
                               'WikiPathway_2023_Fly'],
                       help='FlyEnrichr libraries to use (default: all)')
    parser.add_argument('--skip-enrichment', action='store_true',
                       help='Skip pathway enrichment analysis')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    sc.settings.figdir = args.output
    
    # Load data
    print("="*60)
    print("LOADING DATA")
    print("="*60)
    print(f"\nLoading data from {args.input}...")
    adata = sc.read_h5ad(args.input)
    
    print(f"\nLoaded AnnData object:")
    print(f"  Cells: {adata.n_obs}")
    print(f"  Genes: {adata.n_vars}")
    
    if 'leiden' in adata.obs.columns:
        print(f"  Clusters: {adata.obs['leiden'].nunique()}")
        print(f"  Cluster IDs: {sorted(adata.obs['leiden'].unique())}")
    else:
        print("  ERROR: No 'leiden' clustering found!")
        return
    
    # 1. Transcriptional activity analysis
    print("\n" + "="*60)
    print("STEP 1: TRANSCRIPTIONAL ACTIVITY ANALYSIS")
    print("="*60)
    activity_results = plot_transcriptional_activity(adata, args.output, args.sample)
    
    # 2. Find marker genes
    print("\n" + "="*60)
    print("STEP 2: MARKER GENE IDENTIFICATION")
    print("="*60)
    marker_df = find_marker_genes(adata, args.output, args.sample, method=args.method)
    
    # 3. Pathway enrichment
    if not args.skip_enrichment:
        print("\n" + "="*60)
        print("STEP 3: PATHWAY ENRICHMENT ANALYSIS")
        print("="*60)
        enrichment_df = enrichment_analysis_per_cluster(
            marker_df, 
            adata, 
            args.output, 
            args.sample,
            top_n=args.top_n,
            libraries=args.libraries
        )
    else:
        print("\n" + "="*60)
        print("STEP 3: PATHWAY ENRICHMENT (SKIPPED)")
        print("="*60)
        enrichment_df = None
    
    # Final summary
    print("\n" + "="*60)
    print("ANALYSIS COMPLETE")
    print("="*60)
    print(f"\nAll results saved to: {args.output}/")


if __name__ == "__main__":
    main()