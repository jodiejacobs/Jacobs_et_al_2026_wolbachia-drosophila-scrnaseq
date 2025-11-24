import scanpy as sc
import pandas as pd
import numpy as np
import gseapy as gp
import os
import gzip
import warnings

# Suppress the specific log2 warning from scanpy (it's expected)
warnings.filterwarnings('ignore', message='invalid value encountered in log2')

# Test data:
adata = sc.read_h5ad("/private/groups/russelllab/jodie/scRNAseq/Jacobs_et_al_2026_wolbachia-drosophila-scrnaseq/in_vivo_translation/results/combined/mei_P26_wMel_pipseq.h5ad")

OUTPUT_BASE_DIR = 'test_meiP26wMel_gseapy_prerank_results_all_clusters'
os.makedirs(OUTPUT_BASE_DIR, exist_ok=True)

cluster_ids = adata.obs['leiden'].unique().categories

print(f"Starting GSEA for {len(cluster_ids)} clusters...")

# Load FlyBase annotation mapping
print("Loading FlyBase annotation mapping...")
flybase_annot_path = "/private/groups/russelllab/jodie/scRNAseq/Jacobs_et_al_2026_wolbachia-drosophila-scrnaseq/reference/fbgn_annotation_ID_fb_2025_04.tsv.gz"

with gzip.open(flybase_annot_path, 'rt') as f:
    lines = [line for line in f if not line.startswith('#')]
    
from io import StringIO
mapping_df = pd.read_csv(StringIO(''.join(lines)), sep='\t', header=None)

fbgn_to_symbol = dict(zip(mapping_df[2], mapping_df[0]))
fbgn_to_symbol = {k: v for k, v in fbgn_to_symbol.items() if pd.notna(v) and v != '' and pd.notna(k)}

print(f"Created mapping for {len(fbgn_to_symbol)} genes")

# Check gene sets
print("\nChecking gene sets...")
gene_sets_dir = os.path.join(OUTPUT_BASE_DIR, 'gene_sets')
os.makedirs(gene_sets_dir, exist_ok=True)

FLY_GENE_SETS = []
gene_set_libraries = [
    'GO_Biological_Process_2018',
    'GO_Cellular_Component_2018', 
    'GO_Molecular_Function_2018',
    'KEGG_2019'
]

for library in gene_set_libraries:
    gmt_path = os.path.join(gene_sets_dir, f'{library}.gmt')
    if not os.path.exists(gmt_path):
        print(f"  Downloading {library}...")
        gene_set = gp.get_library(name=library, organism='Fly')
        with open(gmt_path, 'w') as f:
            for term, genes in gene_set.items():
                f.write(f"{term}\tNA\t" + "\t".join(genes) + "\n")
    else:
        print(f"  Using existing {library}")
    FLY_GENE_SETS.append(gmt_path)

print(f"\nSuccessfully loaded {len(FLY_GENE_SETS)} gene set libraries")

# Create summary dataframe
summary_data = []

# Loop through each cluster ID
for cluster_id in cluster_ids:
    print(f"\n{'='*60}")
    print(f"Processing Cluster {cluster_id}")
    print('='*60)
    
    cluster_size = (adata.obs['leiden'] == cluster_id).sum()
    print(f"Cluster size: {cluster_size} cells")
    
    # Run Differential Expression
    DE_KEY = f'rank_genes_{cluster_id}_vs_rest'
    sc.tl.rank_genes_groups(adata, 
                            groupby='leiden', 
                            groups=[cluster_id],        
                            reference='rest', 
                            method='wilcoxon',          
                            key_added=DE_KEY,
                            n_genes=adata.shape[1]
                           )
    
    # Extract Results
    de_results_df = sc.get.rank_genes_groups_df(adata, group=cluster_id, key=DE_KEY)
    
    print(f"\nInitial DE results: {len(de_results_df)} genes")
    print(f"Genes with valid logFC: {np.isfinite(de_results_df['logfoldchanges']).sum()}")
    print(f"Genes with pvals_adj < 0.05: {(de_results_df['pvals_adj'] < 0.05).sum()}")
    
    # Map to gene symbols
    de_results_df['gene_symbol'] = de_results_df['names'].map(fbgn_to_symbol)
    de_results_df = de_results_df.dropna(subset=['gene_symbol'])
    de_results_df = de_results_df[de_results_df['gene_symbol'] != '']
    
    print(f"After gene symbol mapping: {len(de_results_df)} genes")
    
    # Keep only genes with finite logfoldchanges
    de_results_valid = de_results_df[np.isfinite(de_results_df['logfoldchanges'])].copy()
    print(f"After removing NaN/Inf logFC: {len(de_results_valid)} genes")
    
    # Count significant genes with different thresholds
    sig_genes_fdr05 = de_results_valid[de_results_valid['pvals_adj'] < 0.05]
    sig_genes_fdr05_lfc05 = sig_genes_fdr05[np.abs(sig_genes_fdr05['logfoldchanges']) > 0.5]
    sig_genes_fdr05_lfc1 = sig_genes_fdr05[np.abs(sig_genes_fdr05['logfoldchanges']) > 1.0]
    
    print(f"\nSignificant genes (FDR<0.05): {len(sig_genes_fdr05)}")
    print(f"  With |logFC| > 0.5: {len(sig_genes_fdr05_lfc05)}")
    print(f"  With |logFC| > 1.0: {len(sig_genes_fdr05_lfc1)}")
    
    # Use a more lenient threshold for GSEA (just FDR < 0.05)
    sig_genes = sig_genes_fdr05
    sig_up = sig_genes[sig_genes['logfoldchanges'] > 0]
    sig_down = sig_genes[sig_genes['logfoldchanges'] < 0]
    
    print(f"\nUsing FDR<0.05 threshold (no logFC filter for GSEA):")
    print(f"  Total significantly DE genes: {len(sig_genes)}")
    print(f"    Upregulated: {len(sig_up)}")
    print(f"    Downregulated: {len(sig_down)}")
    print(f"  Max |log2FC|: {de_results_valid['logfoldchanges'].abs().max():.2f}")
    
    # Prepare data for GSEA using ALL genes (not just significant ones)
    # GSEA works better with the full ranked list
    
    # Handle 0 adjusted p-values
    min_nonzero_pval = de_results_valid.loc[de_results_valid['pvals_adj'] > 0, 'pvals_adj'].min()
    if pd.isna(min_nonzero_pval) or min_nonzero_pval == 0:
        min_nonzero_pval = 1e-300
    de_results_valid.loc[de_results_valid['pvals_adj'] == 0, 'pvals_adj'] = min_nonzero_pval

    # Calculate ranking metric
    de_results_valid['rank_score'] = np.sign(de_results_valid['logfoldchanges']) * \
                                      (-np.log10(de_results_valid['pvals_adj']))
    
    # Add noise to break ties
    np.random.seed(42)
    noise = np.random.uniform(-0.001, 0.001, size=len(de_results_valid))
    de_results_valid['rank_score'] = de_results_valid['rank_score'] + noise

    # Sort and format for GSEApy
    rnk_df = de_results_valid[['gene_symbol', 'rank_score', 'names', 'pvals_adj', 'logfoldchanges']].sort_values(by='rank_score', ascending=False)
    
    # Save full DE results
    de_output = os.path.join(OUTPUT_BASE_DIR, f'cluster_{cluster_id}_DE_results.csv')
    rnk_df.to_csv(de_output, index=False)
    
    ranked_gene_series = pd.Series(rnk_df['rank_score'].values, index=rnk_df['gene_symbol'].values)
    ranked_gene_series = ranked_gene_series.loc[~ranked_gene_series.index.duplicated(keep='first')]
    
    print(f"\nRanked {len(ranked_gene_series)} unique genes for GSEA")
    if len(ranked_gene_series) > 0:
        print(f"Top 5 upregulated: {ranked_gene_series.head().index.tolist()}")
        print(f"  Scores: {ranked_gene_series.head().values}")
        print(f"Top 5 downregulated: {ranked_gene_series.tail().index.tolist()}")
        print(f"  Scores: {ranked_gene_series.tail().values}")
    
    # Save ranked gene list
    rank_file = os.path.join(OUTPUT_BASE_DIR, f'cluster_{cluster_id}_ranked_genes.rnk')
    ranked_gene_series.to_csv(rank_file, sep='\t', header=False)

    # Run GSEApy prerank
    cluster_outdir = os.path.join(OUTPUT_BASE_DIR, f'cluster_{cluster_id}_results')
    
    print(f"\nRunning GSEApy prerank...")
    
    n_sig_pathways = 0
    try:
        prerank_results = gp.prerank(
            rnk=ranked_gene_series,
            gene_sets=FLY_GENE_SETS,
            outdir=cluster_outdir,
            seed=42,
            threads=4,
            min_size=5,
            max_size=500,
            permutation_num=1000,
            figsize=(6, 5),
            verbose=True
        )
        
        if prerank_results.res2d is not None and not prerank_results.res2d.empty:
            results_path = os.path.join(cluster_outdir, f'GSEApy_results_summary_cluster_{cluster_id}.csv')
            prerank_results.res2d.to_csv(results_path, index=False)
            
            sig_results = prerank_results.res2d[prerank_results.res2d['FDR q-val'] < 0.05]
            n_sig_pathways = len(sig_results)
            
            print(f"\n✓ Results saved to {results_path}")
            print(f"✓ Found {n_sig_pathways} significant pathways (FDR < 0.05)")
            
            if n_sig_pathways > 0:
                print("\nTop 10 significant pathways:")
                top_pathways = sig_results.sort_values('FDR q-val').head(10)
                for idx, row in top_pathways.iterrows():
                    print(f"  {row['Term'][:60]:<60} NES={row['NES']:6.2f}  FDR={row['FDR q-val']:.3e}")
            else:
                print("\nTop 10 pathways by nominal p-value:")
                top_nominal = prerank_results.res2d.sort_values('NOM p-val').head(10)
                for idx, row in top_nominal.iterrows():
                    print(f"  {row['Term'][:60]:<60} NES={row['NES']:6.2f}  Nom p={row['NOM p-val']:.3e}")
                    
    except Exception as e:
        print(f"\nERROR in cluster {cluster_id}: {e}")
        import traceback
        traceback.print_exc()
    
    # Add to summary
    summary_data.append({
        'Cluster': cluster_id,
        'N_cells': cluster_size,
        'N_DE_genes_FDR05': len(sig_genes),
        'N_upregulated': len(sig_up),
        'N_downregulated': len(sig_down),
        'N_DE_logFC1': len(sig_genes_fdr05_lfc1),
        'Max_abs_logFC': de_results_valid['logfoldchanges'].abs().max(),
        'N_genes_ranked': len(ranked_gene_series),
        'N_sig_pathways_FDR05': n_sig_pathways
    })

# Save summary
summary_df = pd.DataFrame(summary_data)
summary_path = os.path.join(OUTPUT_BASE_DIR, 'cluster_summary_detailed.csv')
summary_df.to_csv(summary_path, index=False)

print("\n" + "="*60)
print("SUMMARY TABLE")
print("="*60)
print(summary_df.to_string(index=False))
print(f"\nSummary saved to {summary_path}")

print("\n" + "="*60)
print("All GSEA analyses complete!")
print("="*60)