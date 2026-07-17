import scanpy as sc
import pandas as pd
import numpy as np
from sklearn.decomposition import NMF
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy as gp
import os
from scipy.stats import spearmanr

# Load your data
adata = sc.read_h5ad("/private/groups/russelllab/jodie/scRNAseq/Jacobs_et_al_2026_wolbachia-drosophila-scrnaseq/in_vivo_translation/results/combined/mei_P26_wMel_pipseq.h5ad")

OUTPUT_DIR = 'nmf_gene_programs'
os.makedirs(OUTPUT_DIR, exist_ok=True)

print("="*60)
print("NMF GENE PROGRAM DISCOVERY")
print("="*60)

# ============================================================================
# STEP 1: Prepare Data
# ============================================================================
print("\n1. Preparing data...")

# Use highly variable genes only (reduces noise and computation)
if 'highly_variable' not in adata.var.columns:
    print("   Computing highly variable genes...")
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat_v3')

adata_hvg = adata[:, adata.var.highly_variable].copy()
print(f"   Using {adata_hvg.n_vars} highly variable genes")
print(f"   Across {adata_hvg.n_obs} cells")

# NMF requires non-negative input - use raw counts or normalized data
# Make sure data is dense (not sparse) and non-negative
X = adata_hvg.X.toarray() if hasattr(adata_hvg.X, 'toarray') else adata_hvg.X
X = np.maximum(X, 0)  # Ensure non-negative

print(f"   Data matrix shape: {X.shape}")

# ============================================================================
# STEP 2: Choose Number of Programs
# ============================================================================
print("\n2. Determining optimal number of programs...")

# Try different numbers of programs and evaluate
n_programs_range = range(5, 21, 5)
reconstruction_errors = []

for n in n_programs_range:
    model = NMF(n_components=n, init='nndsvda', random_state=42, max_iter=500)
    W = model.fit_transform(X)
    reconstruction_errors.append(model.reconstruction_err_)
    print(f"   n={n}: reconstruction error = {model.reconstruction_err_:.2f}")

# Plot elbow curve
plt.figure(figsize=(8, 5))
plt.plot(n_programs_range, reconstruction_errors, 'bo-')
plt.xlabel('Number of Programs')
plt.ylabel('Reconstruction Error')
plt.title('NMF Elbow Plot')
plt.grid(True)
plt.savefig(os.path.join(OUTPUT_DIR, 'nmf_elbow_plot.pdf'))
plt.close()

print("   Elbow plot saved!")

# ============================================================================
# STEP 3: Run NMF with Chosen Number of Programs
# ============================================================================
n_programs = 15  # Adjust based on elbow plot or biological knowledge
print(f"\n3. Running NMF with {n_programs} programs...")

model = NMF(
    n_components=n_programs,
    init='nndsvda',  # Non-negative double SVD initialization (usually best)
    random_state=42,
    max_iter=1000,
    tol=1e-4,
    verbose=1
)

W = model.fit_transform(X)  # Cell × Program (usage scores)
H = model.components_       # Program × Gene (gene weights)

print(f"   Final reconstruction error: {model.reconstruction_err_:.2f}")
print(f"   W shape (Cell × Program): {W.shape}")
print(f"   H shape (Program × Gene): {H.shape}")

# ============================================================================
# STEP 4: Extract and Characterize Gene Programs
# ============================================================================
print("\n4. Extracting gene programs...")

programs = {}
program_stats = []

for i in range(n_programs):
    # Get top genes for this program
    gene_weights = H[i, :]
    top_gene_idx = np.argsort(-gene_weights)[:200]  # Top 200 genes
    top_genes = adata_hvg.var_names[top_gene_idx].tolist()
    top_weights = gene_weights[top_gene_idx]
    
    programs[f'Program_{i}'] = {
        'genes': top_genes,
        'weights': top_weights
    }
    
    # Calculate program statistics
    program_usage = W[:, i]
    program_stats.append({
        'Program': f'Program_{i}',
        'N_cells_high': (program_usage > np.percentile(program_usage, 75)).sum(),
        'Mean_usage': program_usage.mean(),
        'Max_usage': program_usage.max(),
        'Top_gene': top_genes[0],
        'Top_5_genes': ', '.join(top_genes[:5])
    })
    
    print(f"\n   === Program {i} ===")
    print(f"   Top 10 genes: {', '.join(top_genes[:10])}")
    print(f"   Mean usage: {program_usage.mean():.3f}")
    print(f"   Max usage: {program_usage.max():.3f}")
    
    # Save full gene list for this program
    program_df = pd.DataFrame({
        'gene': top_genes[:200],
        'weight': top_weights[:200],
        'rank': range(1, 201)
    })
    program_df.to_csv(
        os.path.join(OUTPUT_DIR, f'Program_{i}_genes.csv'),
        index=False
    )

# Save program statistics
stats_df = pd.DataFrame(program_stats)
stats_df.to_csv(os.path.join(OUTPUT_DIR, 'program_statistics.csv'), index=False)
print("\n   Program statistics saved!")

# ============================================================================
# STEP 5: Add Program Scores to AnnData
# ============================================================================
print("\n5. Adding program scores to AnnData...")

for i in range(n_programs):
    adata.obs[f'Program_{i}'] = W[:, i]

# Save updated adata
adata.write(os.path.join(OUTPUT_DIR, 'adata_with_programs.h5ad'))
print("   Updated AnnData saved!")

# ============================================================================
# STEP 6: Visualize Programs on UMAP
# ============================================================================
print("\n6. Creating visualizations...")

# Compute UMAP if not already present
if 'X_umap' not in adata.obsm:
    print("   Computing UMAP...")
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

# Plot all programs
fig, axes = plt.subplots(3, 5, figsize=(20, 12))
axes = axes.flatten()

for i in range(n_programs):
    sc.pl.umap(adata, 
               color=f'Program_{i}',
               ax=axes[i],
               show=False,
               title=f'Program {i}',
               cmap='viridis')

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'all_programs_umap.pdf'))
plt.close()
print("   UMAP plots saved!")

# Plot program usage by cluster
if 'leiden' in adata.obs.columns:
    print("   Creating program × cluster heatmap...")
    
    program_by_cluster = pd.DataFrame()
    for cluster in sorted(adata.obs['leiden'].unique()):
        cluster_cells = adata.obs['leiden'] == cluster
        for i in range(n_programs):
            program_by_cluster.loc[f'Program_{i}', f'Cluster_{cluster}'] = \
                adata.obs.loc[cluster_cells, f'Program_{i}'].mean()
    
    plt.figure(figsize=(10, 8))
    sns.heatmap(program_by_cluster, cmap='viridis', annot=True, fmt='.2f')
    plt.title('Mean Program Usage by Cluster')
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'programs_by_cluster_heatmap.pdf'))
    plt.close()

# ============================================================================
# STEP 7: Run GSEA on Each Program
# ============================================================================
print("\n7. Running GSEA on gene programs...")

# Load FlyBase mapping (reusing your code)
import gzip
from io import StringIO

flybase_annot_path = "/private/groups/russelllab/jodie/scRNAseq/Jacobs_et_al_2026_wolbachia-drosophila-scrnaseq/reference/fbgn_annotation_ID_fb_2025_04.tsv.gz"

with gzip.open(flybase_annot_path, 'rt') as f:
    lines = [line for line in f if not line.startswith('#')]
mapping_df = pd.read_csv(StringIO(''.join(lines)), sep='\t', header=None)
fbgn_to_symbol = dict(zip(mapping_df[2], mapping_df[0]))
fbgn_to_symbol = {k: v for k, v in fbgn_to_symbol.items() 
                  if pd.notna(v) and v != '' and pd.notna(k)}

# Load gene sets
gene_sets_dir = 'test_meiP26wMel_gseapy_prerank_results_all_clusters/gene_sets'
FLY_GENE_SETS = [
    os.path.join(gene_sets_dir, 'GO_Biological_Process_2018.gmt'),
    os.path.join(gene_sets_dir, 'GO_Cellular_Component_2018.gmt'),
    os.path.join(gene_sets_dir, 'GO_Molecular_Function_2018.gmt'),
    os.path.join(gene_sets_dir, 'KEGG_2019.gmt')
]

gsea_results_dir = os.path.join(OUTPUT_DIR, 'gsea_results')
os.makedirs(gsea_results_dir, exist_ok=True)

all_gsea_results = []

for i in range(n_programs):
    print(f"\n   Running GSEA for Program {i}...")
    
    # Get genes and weights for this program
    program_genes = programs[f'Program_{i}']['genes'][:500]  # Top 500
    program_weights = programs[f'Program_{i}']['weights'][:500]
    
    # Map to gene symbols
    gene_symbols = [fbgn_to_symbol.get(g, None) for g in program_genes]
    
    # Create ranked list
    ranked_data = []
    for gene, symbol, weight in zip(program_genes, gene_symbols, program_weights):
        if symbol is not None and symbol != '':
            ranked_data.append((symbol, weight))
    
    if len(ranked_data) < 10:
        print(f"      WARNING: Only {len(ranked_data)} genes mapped. Skipping.")
        continue
    
    # Create Series for GSEA
    ranked_df = pd.DataFrame(ranked_data, columns=['gene', 'weight'])
    ranked_df = ranked_df.sort_values('weight', ascending=False)
    ranked_series = pd.Series(
        ranked_df['weight'].values,
        index=ranked_df['gene'].values
    )
    
    # Remove duplicates
    ranked_series = ranked_series[~ranked_series.index.duplicated(keep='first')]
    
    print(f"      {len(ranked_series)} genes for GSEA")
    
    # Run GSEA
    try:
        program_outdir = os.path.join(gsea_results_dir, f'Program_{i}')
        
        prerank_results = gp.prerank(
            rnk=ranked_series,
            gene_sets=FLY_GENE_SETS,
            outdir=program_outdir,
            seed=42,
            threads=4,
            min_size=5,
            max_size=500,
            permutation_num=1000,
            verbose=False
        )
        
        if prerank_results.res2d is not None and not prerank_results.res2d.empty:
            # Get significant results
            sig_results = prerank_results.res2d[prerank_results.res2d['FDR q-val'] < 0.05]
            
            print(f"      Found {len(sig_results)} significant pathways")
            
            if len(sig_results) > 0:
                print(f"      Top pathway: {sig_results.iloc[0]['Term'][:60]}")
                
                # Add program info to results
                sig_results['Program'] = f'Program_{i}'
                all_gsea_results.append(sig_results)
            
            # Save full results
            prerank_results.res2d.to_csv(
                os.path.join(program_outdir, f'gsea_results_Program_{i}.csv'),
                index=False
            )
    
    except Exception as e:
        print(f"      ERROR: {e}")
        continue

# Combine all GSEA results
if len(all_gsea_results) > 0:
    combined_gsea = pd.concat(all_gsea_results, ignore_index=True)
    combined_gsea = combined_gsea.sort_values(['Program', 'FDR q-val'])
    combined_gsea.to_csv(
        os.path.join(gsea_results_dir, 'all_programs_significant_pathways.csv'),
        index=False
    )
    print(f"\n   Combined GSEA results saved ({len(combined_gsea)} significant pathways)")

# ============================================================================
# STEP 8: Correlate Programs with Metadata
# ============================================================================
print("\n8. Correlating programs with metadata...")

if 'leiden' in adata.obs.columns:
    # Which clusters are enriched for each program?
    enrichment_results = []
    
    for i in range(n_programs):
        program_scores = adata.obs[f'Program_{i}']
        
        for cluster in sorted(adata.obs['leiden'].unique()):
            in_cluster = adata.obs['leiden'] == cluster
            out_cluster = ~in_cluster
            
            mean_in = program_scores[in_cluster].mean()
            mean_out = program_scores[out_cluster].mean()
            fold_enrichment = mean_in / mean_out if mean_out > 0 else float('inf')
            
            enrichment_results.append({
                'Program': f'Program_{i}',
                'Cluster': cluster,
                'Mean_in_cluster': mean_in,
                'Mean_out_cluster': mean_out,
                'Fold_enrichment': fold_enrichment
            })
    
    enrichment_df = pd.DataFrame(enrichment_results)
    enrichment_df.to_csv(
        os.path.join(OUTPUT_DIR, 'program_cluster_enrichment.csv'),
        index=False
    )
    print("   Program-cluster enrichment saved!")

# ============================================================================
# STEP 9: Create Summary Report
# ============================================================================
print("\n9. Creating summary report...")

with open(os.path.join(OUTPUT_DIR, 'SUMMARY.txt'), 'w') as f:
    f.write("="*60 + "\n")
    f.write("NMF GENE PROGRAM DISCOVERY - SUMMARY\n")
    f.write("="*60 + "\n\n")
    
    f.write(f"Dataset: {adata.n_obs} cells × {adata.n_vars} genes\n")
    f.write(f"Highly variable genes used: {adata_hvg.n_vars}\n")
    f.write(f"Number of programs: {n_programs}\n")
    f.write(f"Reconstruction error: {model.reconstruction_err_:.2f}\n\n")
    
    f.write("PROGRAM SUMMARIES:\n")
    f.write("-"*60 + "\n\n")
    
    for i in range(n_programs):
        f.write(f"Program {i}:\n")
        f.write(f"  Top 10 genes: {', '.join(programs[f'Program_{i}']['genes'][:10])}\n")
        f.write(f"  Mean usage: {W[:, i].mean():.3f}\n")
        f.write(f"  Cells with high usage (>75th percentile): {(W[:, i] > np.percentile(W[:, i], 75)).sum()}\n")
        
        # Add top enriched pathway if available
        if len(all_gsea_results) > 0:
            program_pathways = combined_gsea[combined_gsea['Program'] == f'Program_{i}']
            if len(program_pathways) > 0:
                top_pathway = program_pathways.iloc[0]
                f.write(f"  Top enriched pathway: {top_pathway['Term']}\n")
                f.write(f"    NES: {top_pathway['NES']:.2f}, FDR: {top_pathway['FDR q-val']:.3e}\n")
        f.write("\n")

print("\n" + "="*60)
print("NMF ANALYSIS COMPLETE!")
print("="*60)
print(f"\nResults saved to: {OUTPUT_DIR}/")
print("\nKey outputs:")
print(f"  - Program gene lists: Program_X_genes.csv")
print(f"  - GSEA results: gsea_results/")
print(f"  - Visualizations: *.pdf")
print(f"  - Updated AnnData: adata_with_programs.h5ad")
print(f"  - Summary: SUMMARY.txt")
print("="*60)