import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

# Quick functional check for your titer-responsive programs
print("="*70)
print("QUICK FUNCTIONAL CHECK OF TITER-RESPONSIVE PROGRAMS")
print("="*70)

titer_responsive_programs = [0, 3, 6, 9, 11, 13]

import gzip
from io import StringIO

# Load mapping
flybase_annot_path = "/private/groups/russelllab/jodie/scRNAseq/Jacobs_et_al_2026_wolbachia-drosophila-scrnaseq/reference/fbgn_annotation_ID_fb_2025_04.tsv.gz"
with gzip.open(flybase_annot_path, 'rt') as f:
    lines = [line for line in f if not line.startswith('#')]
mapping_df = pd.read_csv(StringIO(''.join(lines)), sep='\t', header=None)
fbgn_to_symbol = dict(zip(mapping_df[2], mapping_df[0]))

# Functional gene sets
functional_keywords = {
    'Ribosomal': ['RpL', 'RpS', 'mrpl', 'mrps'],
    'Mitochondrial': ['mt:', 'ND', 'COX', 'ATP5', 'NADH', 'Cyc', 'Cyt-c'],
    'Immune': ['Rel', 'Dif', 'Toll', 'imd', 'Def', 'Dro', 'Att', 'Cec', 'PGRP'],
    'Metabolism': ['Idh', 'Mdh', 'Gdh', 'G6P', 'Pgk', 'Eno', 'Gapdh', 'Ldh'],
    'Stress_response': ['Hsp', 'Hsc', 'JNK', 'p38'],
    'Proteasome': ['Pros', 'Rpn', 'Rpt', 'Psm'],
    'Translation': ['eIF', 'eEF', 'eRF', 'Taf'],
    'Cell_cycle': ['CycA', 'CycB', 'CycE', 'cdc', 'aurA', 'polo', 'stg']
}

for prog_num in titer_responsive_programs:
    print(f"\n{'='*70}")
    
    # Load correlation info
    corr_df = pd.read_csv('nmf_gene_programs/titer_analysis/program_titer_correlations.csv')
    prog_corr = corr_df[corr_df['Program'] == f'Program_{prog_num}'].iloc[0]
    
    direction = "↑ INDUCED" if prog_corr['Spearman_r'] > 0 else "↓ SUPPRESSED"
    print(f"Program {prog_num} {direction} (r={prog_corr['Spearman_r']:.3f})")
    print(f"{'='*70}")
    
    # Load genes
    gene_file = f'nmf_gene_programs/Program_{prog_num}_genes.csv'
    genes_df = pd.read_csv(gene_file)
    genes_df['symbol'] = genes_df['gene'].map(fbgn_to_symbol).fillna(genes_df['gene'])
    
    top_genes = genes_df['symbol'].head(50).tolist()
    
    # Check functional enrichment
    print("\nFunctional composition (top 50 genes):")
    found_any = False
    
    for func_name, keywords in functional_keywords.items():
        matches = []
        for gene in top_genes:
            gene_upper = str(gene).upper()
            for keyword in keywords:
                if keyword.upper() in gene_upper:
                    matches.append(gene)
                    break
        
        if len(matches) > 0:
            pct = 100 * len(matches) / len(top_genes)
            print(f"  {func_name:20s}: {len(matches):2d}/50 ({pct:4.0f}%) - {', '.join(matches[:5])}")
            found_any = True
    
    if not found_any:
        print("  No obvious functional enrichment detected")
        print(f"\n  Top 15 genes: {', '.join(top_genes[:15])}")
    
    # Load GSEA if available
    gsea_file = f'nmf_gene_programs/gsea_results/Program_{prog_num}/gsea_results_Program_{prog_num}.csv'
    if os.path.exists(gsea_file):
        gsea_df = pd.read_csv(gsea_file)
        sig_gsea = gsea_df[gsea_df['FDR q-val'] < 0.05]
        
        if len(sig_gsea) > 0:
            print(f"\n  GSEA (FDR < 0.05): {len(sig_gsea)} pathways")
            print(f"    Top: {sig_gsea.iloc[0]['Term'][:60]}")

print("\n" + "="*70)
print("SUMMARY")
print("="*70)
print("\nBased on titer correlations:")
print("  - Program 6 (r=0.186): STRONGEST Wolbachia response")
print("  - Program 0 (r=-0.117): SUPPRESSED by Wolbachia")
print("  - Programs 9, 11, 3, 13: Moderate responses")
print("\nCheck the gene lists above to understand what each program represents!")
print("="*70)