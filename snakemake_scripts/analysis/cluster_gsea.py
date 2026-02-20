#!/usr/bin/env python3
"""
Cluster GSEA analysis for Snakemake workflow
Performs GSEA on each cluster vs rest
"""

import argparse
import scanpy as sc
import pandas as pd
import numpy as np
import gseapy as gp
import os
import warnings
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings('ignore')

def parse_args():
    parser = argparse.ArgumentParser(description='Cluster GSEA analysis')
    parser.add_argument('--input', required=True, help='Input h5ad file')
    parser.add_argument('--output_dir', required=True, help='Output directory')
    parser.add_argument('--organism', default='Fly', help='Organism for gene sets')
    parser.add_argument('--gene_id_type', default='flybase', 
                       help='Gene ID type: flybase, symbol, ensembl')
    parser.add_argument('--flybase_annotation', default=None,
                       help='Path to FlyBase annotation file')
    parser.add_argument('--cluster_column', default=None,
                       help='Cluster column name (auto-detect if None)')
    parser.add_argument('--threads', type=int, default=4, help='Number of threads')
    return parser.parse_args()

def detect_cluster_column(adata):
    """Auto-detect cluster column"""
    possible = ['leiden', 'louvain', 'clusters', 'cluster', 'cell_type', 'celltype']
    for col in possible:
        if col in adata.obs.columns:
            print(f"  ✓ Detected cluster column: '{col}'")
            return col
    raise ValueError(f"No cluster column found. Available: {adata.obs.columns.tolist()}")

def load_gene_mapping(flybase_path):
    """Load FlyBase ID to symbol mapping"""
    import gzip
    from io import StringIO
    
    with gzip.open(flybase_path, 'rt') as f:
        lines = [line for line in f if not line.startswith('#')]
    
    mapping_df = pd.read_csv(StringIO(''.join(lines)), sep='\t', header=None)
    gene_map = dict(zip(mapping_df[2], mapping_df[0]))
    gene_map = {k: v for k, v in gene_map.items() if pd.notna(v) and v != ''}
    
    return gene_map

def calculate_ranking_metric(adata, cluster_id, cluster_col):
    """Calculate ranking metric for GSEA"""
    cluster_mask = adata.obs[cluster_col] == cluster_id
    
    if hasattr(adata.X, 'toarray'):
        expr = adata.X.toarray()
    else:
        expr = adata.X
    
    results = []
    
    for gene_idx in range(adata.n_vars):
        gene_id = adata.var_names[gene_idx]
        expr_in = expr[cluster_mask, gene_idx]
        expr_out = expr[~cluster_mask, gene_idx]
        
        if expr_in.std() == 0 and expr_out.std() == 0:
            continue
        
        try:
            statistic, pval = stats.mannwhitneyu(expr_in, expr_out, alternative='two-sided')
        except:
            continue
        
        mean_in = expr_in.mean()
        mean_out = expr_out.mean()
        log2fc = np.log2(mean_in / mean_out) if mean_out > 0 else np.nan
        
        # Rank-biserial correlation
        n1, n2 = len(expr_in), len(expr_out)
        rbc = 1 - (2 * statistic) / (n1 * n2)
        
        results.append({
            'gene_id': gene_id,
            'mean_in': mean_in,
            'mean_out': mean_out,
            'log2fc': log2fc,
            'pval': pval,
            'effect_size': rbc
        })
    
    results_df = pd.DataFrame(results)
    from scipy.stats import false_discovery_control
    results_df['pval_adj'] = false_discovery_control(results_df['pval'].values)
    
    return results_df

def main():
    args = parse_args()
    
    print("="*70)
    print("CLUSTER GSEA ANALYSIS")
    print("="*70)
    
    # Load data
    print("\n1. Loading data...")
    adata = sc.read_h5ad(args.input)
    print(f"   Loaded: {adata.n_obs} cells × {adata.n_vars} genes")
    
    # Detect cluster column
    cluster_col = args.cluster_column or detect_cluster_column(adata)
    cluster_ids = sorted(adata.obs[cluster_col].unique())
    print(f"   Found {len(cluster_ids)} clusters")
    
    # Setup output
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load gene mapping
    print("\n2. Loading gene mapping...")
    if args.gene_id_type == 'flybase' and args.flybase_annotation:
        gene_map = load_gene_mapping(args.flybase_annotation)
        print(f"   Loaded {len(gene_map)} gene mappings")
    else:
        gene_map = {g: g for g in adata.var_names}
        print("   Using gene names as-is")
    
    # Load gene sets
    print("\n3. Loading gene sets...")
    gene_sets_dir = os.path.join(args.output_dir, 'gene_sets')
    os.makedirs(gene_sets_dir, exist_ok=True)
    
    libraries = [
        'GO_Biological_Process_2018',
        'GO_Cellular_Component_2018',
        'GO_Molecular_Function_2018',
        'KEGG_2019'
    ]
    
    gene_sets = {}
    for library in libraries:
        gmt_path = os.path.join(gene_sets_dir, f'{library}.gmt')
        if not os.path.exists(gmt_path):
            print(f"   Downloading {library}...")
            try:
                gs = gp.get_library(name=library, organism=args.organism)
                with open(gmt_path, 'w') as f:
                    for term, genes in gs.items():
                        f.write(f"{term}\tNA\t" + "\t".join(genes) + "\n")
                gene_sets.update(gs)
            except:
                print(f"   WARNING: Could not download {library}")
        else:
            with open(gmt_path, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        gene_sets[parts[0]] = parts[2:]
    
    print(f"   Loaded {len(gene_sets)} gene sets")
    
    # Process each cluster
    print("\n4. Processing clusters...")
    summary_data = []
    
    for cluster_id in cluster_ids:
        print(f"\n   Cluster {cluster_id}")
        cluster_size = (adata.obs[cluster_col] == cluster_id).sum()
        
        # Calculate DE
        de_results = calculate_ranking_metric(adata, cluster_id, cluster_col)
        de_results['gene_symbol'] = de_results['gene_id'].map(gene_map)
        de_results = de_results.dropna(subset=['gene_symbol'])
        
        # Ranking metric
        de_results['ranking_metric'] = (
            np.sign(de_results['log2fc']) * 
            (-np.log10(de_results['pval_adj'])) * 
            np.abs(de_results['effect_size'])
        )
        de_results = de_results[np.isfinite(de_results['ranking_metric'])]
        
        if len(de_results) < 100:
            print(f"   WARNING: Only {len(de_results)} genes, skipping")
            continue
        
        # Create ranked list
        de_results = de_results.sort_values('ranking_metric', ascending=False)
        rnk = de_results[['gene_symbol', 'ranking_metric']].drop_duplicates(
            subset='gene_symbol', keep='first').set_index('gene_symbol')['ranking_metric']
        
        # Save
        cluster_dir = os.path.join(args.output_dir, f'cluster_{cluster_id}')
        os.makedirs(cluster_dir, exist_ok=True)
        de_results.to_csv(os.path.join(cluster_dir, 'DE_results.csv'), index=False)
        
        # Run GSEA
        try:
            pre_res = gp.prerank(
                rnk=rnk,
                gene_sets=gene_sets,
                processes=args.threads,
                permutation_num=1000,
                outdir=cluster_dir,
                seed=42,
                min_size=15,
                max_size=500,
                verbose=False
            )
            
            if pre_res.res2d is not None and not pre_res.res2d.empty:
                pre_res.res2d.to_csv(os.path.join(cluster_dir, 'gsea_results.csv'), index=False)
                sig = pre_res.res2d[pre_res.res2d['FDR q-val'] < 0.05]
                n_sig = len(sig)
                print(f"   Found {n_sig} significant pathways")
                
                summary_data.append({
                    'Cluster': cluster_id,
                    'N_cells': cluster_size,
                    'N_DE_genes': len(de_results[de_results['pval_adj'] < 0.05]),
                    'N_pathways': n_sig
                })
        except Exception as e:
            print(f"   ERROR: {e}")
    
    # Save summary
    if summary_data:
        summary_df = pd.DataFrame(summary_data)
        summary_df.to_csv(os.path.join(args.output_dir, 'cluster_summary.csv'), index=False)
        print("\n" + "="*70)
        print("SUMMARY")
        print("="*70)
        print(summary_df.to_string(index=False))
    
    print("\n✓ Analysis complete!")
    print(f"Results in: {args.output_dir}")

if __name__ == '__main__':
    main()
