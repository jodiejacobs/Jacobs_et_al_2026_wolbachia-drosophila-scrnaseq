import pandas as pd
import os

# Load all the results
titer_dir = 'nmf_gene_programs/titer_analysis'
prog_dir = 'nmf_gene_programs/program_interpretation'

print("="*70)
print("WOLBACHIA TITER-RESPONSIVE PROGRAMS - BIOLOGICAL INTERPRETATION")
print("="*70)

# Load correlation results
corr_df = pd.read_csv(os.path.join(titer_dir, 'program_titer_correlations.csv'))
corr_df = corr_df.sort_values('Abs_Spearman_r', ascending=False)

# Load program annotations (if you ran the interpretation script)
if os.path.exists(os.path.join(prog_dir, 'PROGRAM_SUMMARY.csv')):
    summary_df = pd.read_csv(os.path.join(prog_dir, 'PROGRAM_SUMMARY.csv'))
else:
    summary_df = None

# Get titer-responsive programs
responsive = corr_df[corr_df['Spearman_FDR'] < 0.05].copy()

print("\n" + "="*70)
print("TITER-RESPONSIVE PROGRAMS (FDR < 0.05)")
print("="*70)

for idx, row in responsive.iterrows():
    program_num = int(row['Program'].split('_')[1])
    
    print(f"\n{'='*70}")
    print(f"{row['Program']} - {'INDUCED' if row['Spearman_r'] > 0 else 'SUPPRESSED'} by Wolbachia")
    print(f"{'='*70}")
    
    print(f"\nStatistics:")
    print(f"  Spearman correlation: {row['Spearman_r']:.3f}")
    print(f"  FDR: {row['Spearman_FDR']:.2e}")
    print(f"  Pearson correlation: {row['Pearson_r']:.3f}")
    
    # Load gene list
    gene_file = f'nmf_gene_programs/Program_{program_num}_genes.csv'
    if os.path.exists(gene_file):
        genes_df = pd.read_csv(gene_file)
        
        # Map to symbols if we have them
        import gzip
        from io import StringIO
        
        flybase_annot_path = "/private/groups/russelllab/jodie/scRNAseq/Jacobs_et_al_2026_wolbachia-drosophila-scrnaseq/reference/fbgn_annotation_ID_fb_2025_04.tsv.gz"
        with gzip.open(flybase_annot_path, 'rt') as f:
            lines = [line for line in f if not line.startswith('#')]
        mapping_df = pd.read_csv(StringIO(''.join(lines)), sep='\t', header=None)
        fbgn_to_symbol = dict(zip(mapping_df[2], mapping_df[0]))
        
        genes_df['symbol'] = genes_df['gene'].map(fbgn_to_symbol)
        genes_df['symbol'] = genes_df['symbol'].fillna(genes_df['gene'])
        
        print(f"\nTop 30 genes:")
        for i, gene_row in genes_df.head(30).iterrows():
            print(f"  {i+1:2d}. {gene_row['symbol']:20s} (weight: {gene_row['weight']:.3f})")
    
    # Check GSEA results
    gsea_file = f'nmf_gene_programs/gsea_results/Program_{program_num}/gsea_results_Program_{program_num}.csv'
    if os.path.exists(gsea_file):
        gsea_df = pd.read_csv(gsea_file)
        sig_gsea = gsea_df[gsea_df['FDR q-val'] < 0.05]
        
        if len(sig_gsea) > 0:
            print(f"\nEnriched pathways (FDR < 0.05):")
            for i, gsea_row in sig_gsea.head(10).iterrows():
                print(f"  - {gsea_row['Term'][:65]}")
                print(f"    NES: {gsea_row['NES']:.2f}, FDR: {gsea_row['FDR q-val']:.2e}")
        else:
            print(f"\nNo significantly enriched pathways (FDR < 0.05)")
            print(f"Top nominal pathways:")
            for i, gsea_row in gsea_df.sort_values('NOM p-val').head(5).iterrows():
                print(f"  - {gsea_row['Term'][:65]}")
                print(f"    NES: {gsea_row['NES']:.2f}, Nom p: {gsea_row['NOM p-val']:.2e}")
    
    # Add functional annotation if available
    if summary_df is not None:
        prog_summary = summary_df[summary_df['Program'] == row['Program']]
        if len(prog_summary) > 0:
            print(f"\nAnnotation: {prog_summary.iloc[0]['Annotation']}")

print("\n" + "="*70)
print("BIOLOGICAL INTERPRETATION GUIDE")
print("="*70)

print("\n1. PROGRAM 6 (r=0.186, STRONGEST POSITIVE):")
print("   - Most strongly induced by Wolbachia")
print("   - Look for: ribosomal genes, translation machinery")
print("   - Interpretation: Wolbachia hijacking host translation?")

print("\n2. PROGRAM 0 (r=-0.117, STRONGEST NEGATIVE):")
print("   - Suppressed by Wolbachia infection")
print("   - Look for: mitochondrial genes, metabolism")
print("   - Interpretation: Metabolic changes or mitochondrial interference?")

print("\n3. PROGRAM 9 (r=0.136, POSITIVE):")
print("   - Induced by Wolbachia")
print("   - Check top genes for functional clues")

print("\n4. PROGRAMS 3, 11, 13 (weak but significant):")
print("   - Subtle titer-dependent responses")
print("   - May represent secondary effects or specific cell types")

print("\n" + "="*70)
print("NEXT STEPS")
print("="*70)
print("\n1. Look at the PDFs in nmf_gene_programs/titer_analysis/:")
print("   - umap_titer_and_programs.pdf: Are these programs spatially clustered?")
print("   - top_programs_scatter.pdf: Linear or threshold response?")
print("   - programs_by_titer_quartile.pdf: Dose-dependent patterns?")

print("\n2. Check if titer-responsive programs match specific clusters:")
print("   - Load: nmf_gene_programs/program_cluster_enrichment.csv")

print("\n3. Compare to your cluster GSEA:")
print("   - Do titer-responsive programs match pathways from infected clusters?")

print("\n4. Biological questions to ask:")
print("   - Program 6: Is this ribosomal? (Wolbachia needs host ribosomes)")
print("   - Program 0: Is this mitochondrial? (Wolbachia-mitochondria competition)")
print("   - Program 9: What cellular process does this represent?")

print("\n" + "="*70)