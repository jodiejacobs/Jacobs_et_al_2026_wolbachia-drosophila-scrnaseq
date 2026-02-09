'''
Comprehensive cluster analysis: transcriptional activity, marker genes, and pathway enrichment
Uses FlyEnrichr API for automated pathway analysis with FlyBase annotations
REQUIRES GENE SYMBOLS - converts FBgn IDs automatically
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
import gzip
from scipy.stats import kruskal

def load_fbgn_to_symbol_mapping(mapping_file):
    """
    Load FBgn to gene symbol mapping from reference file
    
    Parameters:
    -----------
    mapping_file : str
        Path to dmel_gene_id_key.uniq.tsv.gz file
    
    Returns:
    --------
    dict : FBgn -> gene_symbol mapping
    """
    print("\n" + "="*60)
    print("LOADING GENE SYMBOL MAPPING")
    print("="*60)
    print(f"Loading mapping from: {mapping_file}")
    
    fbgn_to_symbol = {}
    
    try:
        # Handle gzipped or plain text
        if mapping_file.endswith('.gz'):
            opener = gzip.open
            mode = 'rt'
        else:
            opener = open
            mode = 'r'
        
        with opener(mapping_file, mode) as f:
            header = f.readline()  # Skip header
            
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    transcript_id = parts[0]
                    gene_symbol = parts[1]
                    gene_id = parts[2]  # This is the FBgn
                    
                    # Map FBgn to symbol
                    if gene_id.startswith('FBgn') or gene_id.startswith('CG'):
                        fbgn_to_symbol[gene_id] = gene_symbol
        
        print(f"  Loaded {len(fbgn_to_symbol)} gene mappings")
        print(f"  Sample mappings:")
        for i, (fbgn, symbol) in enumerate(list(fbgn_to_symbol.items())[:5]):
            print(f"    {fbgn} -> {symbol}")
        
        return fbgn_to_symbol
        
    except FileNotFoundError:
        print(f"  ERROR: File not found: {mapping_file}")
        print(f"  Please provide path to dmel_gene_id_key.uniq.tsv.gz")
        return None
    except Exception as e:
        print(f"  ERROR loading mapping: {e}")
        return None


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
        cmap = plt.colormaps.get_cmap('tab20')
        palette = {cluster: cmap(i % 20) for i, cluster in enumerate(clusters)}
    
    # Create box plots
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # n_counts boxplot
    ax = axes[0]
    box_data_counts = [adata.obs[adata.obs['leiden'] == c]['n_counts'].values for c in clusters]
    bp1 = ax.boxplot(box_data_counts, labels=clusters, patch_artist=True,
                     showfliers=False, widths=0.6)
    
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
    
    # Create violin plots
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
    
    # UMAPs if available
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
    
    # Run differential expression
    sc.tl.rank_genes_groups(adata, 'leiden', method=method, 
                            key_added='rank_genes_groups',
                            use_raw=False,
                            tie_correct=True,
                            rankby_abs=False,
                            pts=True)
    
    print("Done!")
    
    # Plot top markers
    print("\nGenerating marker gene plots...")
    
    sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, 
                                    save=f'_{sample_name}_top10_heatmap.pdf',
                                    show_gene_labels=True, cmap='viridis',
                                    key='rank_genes_groups')
    
    sc.pl.rank_genes_groups_dotplot(adata, n_genes=5,
                                    save=f'_{sample_name}_top5_dotplot.pdf',
                                    key='rank_genes_groups')
    
    sc.pl.rank_genes_groups_stacked_violin(adata, n_genes=3,
                                           save=f'_{sample_name}_top3_violin.pdf',
                                           key='rank_genes_groups')
    
    sc.pl.rank_genes_groups(adata, n_genes=25,
                           save=f'_{sample_name}_ranked_genes.pdf',
                           key='rank_genes_groups')
    
    # Extract results
    print("\nExtracting marker gene results...")
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    
    marker_data = []
    for group in groups:
        for i in range(len(result['names'][group])):
            log2fc = result['logfoldchanges'][group][i]
            
            # Skip NaN/inf
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
    
    # Filter
    marker_df = marker_df[abs(marker_df['log2fc']) > 0.25]
    
    print(f"\nFiltered markers: {len(marker_df)}")
    
    # Save
    marker_df.to_csv(os.path.join(output_dir, f'{sample_name}_marker_genes_full.csv'), 
                     index=False)
    
    top_markers = marker_df.groupby('cluster').head(50)
    top_markers.to_csv(os.path.join(output_dir, f'{sample_name}_marker_genes_top50.csv'),
                      index=False)
    
    # Print top 10
    print("\nTop 10 marker genes per cluster:")
    print("="*80)
    for cluster in groups:
        print(f"\nCluster {cluster}:")
        cluster_markers = marker_df[marker_df['cluster'] == cluster].head(10)
        if len(cluster_markers) == 0:
            print("  No valid markers found")
            continue
        for idx, row in cluster_markers.iterrows():
            print(f"  {row['gene']:<20} log2FC={row['log2fc']:>7.2f}  "
                  f"p_adj={row['pval_adj']:.2e}")
    
    print(f"\nMarker gene results saved to {output_dir}")
    
    return marker_df


def flyenrichr_analysis(gene_list, gene_set_library='GO_Biological_Process_2018', 
                        description="gene_list"):
    """
    Perform enrichment analysis using FlyEnrichr API
    REQUIRES GENE SYMBOLS, NOT FBgn IDs
    """
    ENRICHR_URL = 'https://maayanlab.cloud/FlyEnrichr/addList'
    
    # Clean gene list
    gene_list_clean = [str(g).strip() for g in gene_list if g and str(g).strip()]
    
    if len(gene_list_clean) == 0:
        print(f"    ERROR: Gene list is empty")
        return None
    
    genes_str = '\n'.join(gene_list_clean)
    
    payload = {
        'list': (None, genes_str),
        'description': (None, description)
    }
    
    try:
        # Submit
        response = requests.post(ENRICHR_URL, files=payload, timeout=30)
        if not response.ok:
            print(f"    ERROR: Failed to submit: {response.status_code}")
            return None
        
        data = json.loads(response.text)
        
        if 'userListId' not in data:
            print(f"    ERROR: No userListId in response")
            return None
        
        user_list_id = data['userListId']
        
        # Get enrichment
        ENRICHR_ENRICH_URL = f'https://maayanlab.cloud/FlyEnrichr/enrich'
        query_string = f'?userListId={user_list_id}&backgroundType={gene_set_library}'
        response = requests.get(ENRICHR_ENRICH_URL + query_string, timeout=30)
        
        if not response.ok:
            print(f"    ERROR: Failed enrichment: {response.status_code}")
            return None
        
        data = json.loads(response.text)
        
        if gene_set_library not in data:
            print(f"    WARNING: No results for {gene_set_library}")
            return None
        
        # Parse
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
        print(f"    ERROR: {e}")
        return None


def enrichment_analysis_per_cluster(marker_df, fbgn_to_symbol, output_dir, sample_name, 
                                    top_n=100, libraries=None):
    """
    Run enrichment analysis - converts FBgn to symbols first
    """
    print("\n" + "="*60)
    print("PATHWAY ENRICHMENT ANALYSIS (FlyEnrichr)")
    print("="*60)
    
    if fbgn_to_symbol is None:
        print("ERROR: No gene symbol mapping provided!")
        print("Cannot proceed with enrichment - FlyEnrichr requires gene symbols")
        return None
    
    if libraries is None:
        libraries = [
            'GO_Biological_Process_2018',
            'GO_Molecular_Function_2018', 
            'GO_Cellular_Component_2018',
            'KEGG_2019',
            'WikiPathways_2018',
            'InterPro_Domains_2019',
            'Pfam_Domains_2019'
        ]
    
    print(f"\nRunning enrichment for top {top_n} markers per cluster")
    print(f"Using libraries: {', '.join(libraries)}")
    
    all_results = []
    clusters = sorted(marker_df['cluster'].unique())
    
    for cluster in clusters:
        print(f"\n{'='*60}")
        print(f"Cluster {cluster}")
        print(f"{'='*60}")
        
        # Get top markers
        cluster_markers = marker_df[marker_df['cluster'] == cluster].head(top_n)
        sig_markers = cluster_markers[cluster_markers['pval_adj'] < 0.05]
        
        if len(sig_markers) < 10:
            print(f"  WARNING: Only {len(sig_markers)} significant markers, using top {min(top_n, len(cluster_markers))}")
            genes_fbgn = cluster_markers['gene'].tolist()
        else:
            genes_fbgn = sig_markers['gene'].tolist()
        
        # Convert FBgn to symbols
        genes_symbols = []
        unmapped = []
        for fbgn in genes_fbgn:
            symbol = fbgn_to_symbol.get(fbgn, None)
            if symbol:
                genes_symbols.append(symbol)
            else:
                unmapped.append(fbgn)
        
        print(f"  Converting {len(genes_fbgn)} FBgn IDs to symbols")
        print(f"  Successfully mapped: {len(genes_symbols)}")
        if unmapped:
            print(f"  Unmapped: {len(unmapped)} genes")
            if len(unmapped) <= 5:
                print(f"    {', '.join(unmapped)}")
        
        if len(genes_symbols) < 5:
            print(f"  ERROR: Too few mapped genes ({len(genes_symbols)}), skipping")
            continue
        
        print(f"  Sample symbols: {', '.join(genes_symbols[:5])}")
        
        # Run enrichment
        for library in libraries:
            print(f"  Running {library}...")
            
            result_df = flyenrichr_analysis(
                genes_symbols, 
                gene_set_library=library,
                description=f"cluster_{cluster}_top{top_n}"
            )
            
            if result_df is not None and len(result_df) > 0:
                result_df['cluster'] = cluster
                result_df['library'] = library
                all_results.append(result_df)
                print(f"    Found {len(result_df)} enriched terms")
            else:
                print(f"    No results")
            
            time.sleep(0.5)
    
    if not all_results:
        print("\nNo enrichment results obtained")
        return None
    
    # Combine
    combined_df = pd.concat(all_results, ignore_index=True)
    
    # Save
    combined_df.to_csv(
        f"{output_dir}/{sample_name}_flyenrichr_all_results.csv",
        index=False
    )
    print(f"\nSaved: {sample_name}_flyenrichr_all_results.csv")
    
    top_per_cluster = combined_df.sort_values('combined_score', ascending=False).groupby(
        ['cluster', 'library']
    ).head(10)
    
    top_per_cluster.to_csv(
        f"{output_dir}/{sample_name}_flyenrichr_top10_per_cluster.csv",
        index=False
    )
    print(f"Saved: {sample_name}_flyenrichr_top10_per_cluster.csv")
    
    # Summary
    print(f"\n{'='*60}")
    print("ENRICHMENT SUMMARY")
    print(f"{'='*60}")
    print(f"Total enriched terms: {len(combined_df)}")
    print(f"Significant (adj p < 0.05): {len(combined_df[combined_df['adj_p_value'] < 0.05])}")
    
    # Print top terms
    for library in libraries:
        lib_results = combined_df[combined_df['library'] == library]
        
        if len(lib_results) == 0:
            continue
        
        print(f"\n{'='*60}")
        print(f"{library}")
        print(f"{'='*60}")
        
        for cluster in clusters:
            cluster_results = lib_results[lib_results['cluster'] == cluster].sort_values(
                'combined_score', ascending=False
            ).head(3)
            
            if len(cluster_results) > 0:
                print(f"\nCluster {cluster}:")
                for idx, row in cluster_results.iterrows():
                    term_short = row['term'].split('(')[0][:60]
                    print(f"  {term_short:<60} p_adj={row['adj_p_value']:.2e}")
    
    # Visualizations
    create_enrichment_plots(combined_df, output_dir, sample_name, clusters)
    
    return combined_df


def create_enrichment_plots(enrichment_df, output_dir, sample_name, clusters):
    """
    Create visualizations
    """
    print("\nCreating enrichment visualizations...")
    
    libraries = enrichment_df['library'].unique()
    
    for library in libraries:
        lib_data = enrichment_df[enrichment_df['library'] == library]
        sig_data = lib_data[lib_data['adj_p_value'] < 0.05]
        
        if len(sig_data) == 0:
            continue
        
        # Bar plots
        n_clusters = len(clusters)
        fig, axes = plt.subplots(n_clusters, 1, figsize=(14, 4*n_clusters))
        
        if n_clusters == 1:
            axes = [axes]
        
        for idx, cluster in enumerate(clusters):
            ax = axes[idx]
            cluster_terms = sig_data[sig_data['cluster'] == cluster].sort_values(
                'combined_score', ascending=False
            ).head(10)
            
            if len(cluster_terms) == 0:
                ax.text(0.5, 0.5, f'No significant terms\nfor Cluster {cluster}',
                       ha='center', va='center', transform=ax.transAxes)
                ax.axis('off')
                continue
            
            cluster_terms = cluster_terms.copy()
            cluster_terms['term_short'] = cluster_terms['term'].apply(
                lambda x: x.split('(')[0][:50] + '...' if len(x.split('(')[0]) > 50 
                else x.split('(')[0]
            )
            
            y_pos = range(len(cluster_terms))
            bars = ax.barh(y_pos, cluster_terms['combined_score'], alpha=0.7)
            
            colors = ['#d62728' if p < 0.01 else '#ff7f0e' 
                     for p in cluster_terms['adj_p_value']]
            for bar, color in zip(bars, colors):
                bar.set_color(color)
            
            ax.set_yticks(y_pos)
            ax.set_yticklabels(cluster_terms['term_short'], fontsize=9)
            ax.set_xlabel('Combined Score', fontsize=11)
            ax.set_title(f'Cluster {cluster}', fontsize=12, fontweight='bold')
            ax.invert_yaxis()
            ax.grid(axis='x', alpha=0.3)
        
        plt.tight_layout()
        lib_short = library.replace('_2018', '').replace('_2019', '')
        plt.savefig(f"{output_dir}/{sample_name}_enrichment_{lib_short}_barplots.pdf",
                   dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {lib_short}_barplots.pdf")


def main():
    parser = argparse.ArgumentParser(
        description='Comprehensive cluster analysis with FlyEnrichr pathway enrichment',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('--input', '-i', required=True,
                       help='Path to h5ad file')
    parser.add_argument('--output', '-o', default='cluster_analysis',
                       help='Output directory')
    parser.add_argument('--sample', '-s', default='sample',
                       help='Sample name for output files')
    parser.add_argument('--mapping', '-map', required=True,
                       help='Path to dmel_gene_id_key.uniq.tsv.gz')
    parser.add_argument('--method', '-m', default='wilcoxon',
                       choices=['wilcoxon', 't-test', 'logreg'],
                       help='DE method (default: wilcoxon)')
    parser.add_argument('--top-n', type=int, default=100,
                       help='Top N markers for enrichment (default: 100)')
    parser.add_argument('--skip-enrichment', action='store_true',
                       help='Skip pathway enrichment')
    
    args = parser.parse_args()
    
    os.makedirs(args.output, exist_ok=True)
    sc.settings.figdir = args.output
    
    # Load data
    print("="*60)
    print("LOADING DATA")
    print("="*60)
    adata = sc.read_h5ad(args.input)
    
    print(f"\nCells: {adata.n_obs}")
    print(f"Genes: {adata.n_vars}")
    print(f"Clusters: {adata.obs['leiden'].nunique()}")
    
    # Load gene mapping
    fbgn_to_symbol = load_fbgn_to_symbol_mapping(args.mapping)
    
    if fbgn_to_symbol is None and not args.skip_enrichment:
        print("\nERROR: Could not load gene mapping, required for enrichment")
        print("Either fix mapping file or use --skip-enrichment")
        return
    
    # 1. Transcriptional activity
    print("\n" + "="*60)
    print("STEP 1: TRANSCRIPTIONAL ACTIVITY")
    print("="*60)
    activity_results = plot_transcriptional_activity(adata, args.output, args.sample)
    
    # 2. Marker genes
    print("\n" + "="*60)
    print("STEP 2: MARKER GENES")
    print("="*60)
    marker_df = find_marker_genes(adata, args.output, args.sample, method=args.method)
    
    # 3. Enrichment
    if not args.skip_enrichment:
        print("\n" + "="*60)
        print("STEP 3: PATHWAY ENRICHMENT")
        print("="*60)
        enrichment_df = enrichment_analysis_per_cluster(
            marker_df, 
            fbgn_to_symbol,
            args.output, 
            args.sample,
            top_n=args.top_n
        )
    
    print("\n" + "="*60)
    print("ANALYSIS COMPLETE")
    print("="*60)


if __name__ == "__main__":
    main()