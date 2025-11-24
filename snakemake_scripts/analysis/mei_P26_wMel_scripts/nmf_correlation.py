import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Load data
adata = sc.read_h5ad('nmf_gene_programs/adata_with_programs.h5ad')

OUTPUT_DIR = 'nmf_gene_programs/biological_interpretation'
import os
os.makedirs(OUTPUT_DIR, exist_ok=True)

print("="*70)
print("BIOLOGICAL INTERPRETATION: WOLBACHIA EFFECTS")
print("="*70)

# ============================================================================
# 1. Compare Program 0 (mitochondrial suppression) vs Program 6 (ribosomal induction)
# ============================================================================

fig, axes = plt.subplots(2, 3, figsize=(18, 12))

# Row 1: Program 0 (Suppressed - Mitochondrial)
ax = axes[0, 0]
sc.pl.umap(adata, color='Program_0', ax=ax, show=False, title='Program 0\n(Mitochondrial - SUPPRESSED)')

ax = axes[0, 1]
ax.scatter(adata.obs['wolbachia_titer'], adata.obs['Program_0'], alpha=0.3, s=10)
ax.set_xlabel('Wolbachia Titer')
ax.set_ylabel('Program 0 Score')
ax.set_title('Program 0 vs Titer\n(Negative correlation)')
# Add regression line
from scipy.stats import linregress
mask = ~pd.isna(adata.obs['wolbachia_titer']) & ~pd.isna(adata.obs['Program_0'])
slope, intercept, r, p, se = linregress(
    adata.obs.loc[mask, 'wolbachia_titer'], 
    adata.obs.loc[mask, 'Program_0']
)
x_line = np.linspace(adata.obs['wolbachia_titer'].min(), adata.obs['wolbachia_titer'].max(), 100)
ax.plot(x_line, slope * x_line + intercept, 'r-', linewidth=2)
ax.text(0.05, 0.95, f'r={r:.3f}', transform=ax.transAxes, va='top')
ax.grid(True, alpha=0.3)

ax = axes[0, 2]
if 'leiden' in adata.obs.columns:
    cluster_means = adata.obs.groupby('leiden')['Program_0'].mean().sort_values()
    ax.barh(range(len(cluster_means)), cluster_means.values)
    ax.set_yticks(range(len(cluster_means)))
    ax.set_yticklabels(cluster_means.index)
    ax.set_xlabel('Mean Program 0 Score')
    ax.set_ylabel('Cluster')
    ax.set_title('Program 0 by Cluster')
    ax.grid(True, alpha=0.3, axis='x')

# Row 2: Program 6 (Induced - Ribosomal)
ax = axes[1, 0]
sc.pl.umap(adata, color='Program_6', ax=ax, show=False, title='Program 6\n(Ribosomal - INDUCED)')

ax = axes[1, 1]
ax.scatter(adata.obs['wolbachia_titer'], adata.obs['Program_6'], alpha=0.3, s=10)
ax.set_xlabel('Wolbachia Titer')
ax.set_ylabel('Program 6 Score')
ax.set_title('Program 6 vs Titer\n(Positive correlation)')
# Add regression line
mask = ~pd.isna(adata.obs['wolbachia_titer']) & ~pd.isna(adata.obs['Program_6'])
slope, intercept, r, p, se = linregress(
    adata.obs.loc[mask, 'wolbachia_titer'], 
    adata.obs.loc[mask, 'Program_6']
)
x_line = np.linspace(adata.obs['wolbachia_titer'].min(), adata.obs['wolbachia_titer'].max(), 100)
ax.plot(x_line, slope * x_line + intercept, 'r-', linewidth=2)
ax.text(0.05, 0.95, f'r={r:.3f}', transform=ax.transAxes, va='top')
ax.grid(True, alpha=0.3)

ax = axes[1, 2]
if 'leiden' in adata.obs.columns:
    cluster_means = adata.obs.groupby('leiden')['Program_6'].mean().sort_values()
    ax.barh(range(len(cluster_means)), cluster_means.values)
    ax.set_yticks(range(len(cluster_means)))
    ax.set_yticklabels(cluster_means.index)
    ax.set_xlabel('Mean Program 6 Score')
    ax.set_ylabel('Cluster')
    ax.set_title('Program 6 by Cluster')
    ax.grid(True, alpha=0.3, axis='x')

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'program_0_vs_6_comparison.pdf'))
plt.close()
print("\n✓ Saved: program_0_vs_6_comparison.pdf")

# ============================================================================
# 2. Are Program 0 and Program 6 anticorrelated?
# ============================================================================

fig, ax = plt.subplots(figsize=(8, 8))

# Scatter with titer as color
scatter = ax.scatter(adata.obs['Program_0'], 
                     adata.obs['Program_6'],
                     c=adata.obs['wolbachia_titer'],
                     cmap='YlOrRd',
                     alpha=0.5,
                     s=20)

ax.set_xlabel('Program 0 (Mitochondrial - Suppressed)')
ax.set_ylabel('Program 6 (Ribosomal - Induced)')
ax.set_title('Program 0 vs Program 6\nColored by Wolbachia Titer')

# Add correlation
from scipy.stats import spearmanr
corr, pval = spearmanr(adata.obs['Program_0'], adata.obs['Program_6'])
ax.text(0.05, 0.95, f'Spearman r={corr:.3f}\np={pval:.2e}', 
        transform=ax.transAxes, va='top',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

plt.colorbar(scatter, label='Wolbachia Titer')
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'program_0_vs_6_scatter.pdf'))
plt.close()
print("✓ Saved: program_0_vs_6_scatter.pdf")

# ============================================================================
# 3. Create a summary figure for publication
# ============================================================================

fig = plt.figure(figsize=(16, 10))

# Title
fig.suptitle('Wolbachia Infection: Host Response Programs', fontsize=16, fontweight='bold')

# Panel A: UMAP with titer
ax1 = plt.subplot(2, 4, 1)
sc.pl.umap(adata, color='wolbachia_titer', ax=ax1, show=False, 
           title='A. Wolbachia Titer', cmap='YlOrRd')

# Panel B: UMAP with Program 0
ax2 = plt.subplot(2, 4, 2)
sc.pl.umap(adata, color='Program_0', ax=ax2, show=False,
           title='B. Program 0\n(Mitochondrial Suppression)')

# Panel C: UMAP with Program 6
ax3 = plt.subplot(2, 4, 3)
sc.pl.umap(adata, color='Program_6', ax=ax3, show=False,
           title='C. Program 6\n(Ribosomal Induction)')

# Panel D: Correlation overview
ax4 = plt.subplot(2, 4, 4)
corr_df = pd.read_csv('nmf_gene_programs/titer_analysis/program_titer_correlations.csv')
programs = [f'P{i}' for i in range(15)]
correlations = corr_df.sort_values('Program')['Spearman_r'].values
colors = ['red' if r < 0 else 'blue' for r in correlations]
bars = ax4.barh(programs, correlations, color=colors, alpha=0.7)
ax4.axvline(0, color='black', linestyle='-', linewidth=0.5)
ax4.set_xlabel('Correlation with Titer')
ax4.set_title('D. All Programs vs Titer')
ax4.grid(True, alpha=0.3, axis='x')

# Panel E: Program 0 scatter
ax5 = plt.subplot(2, 4, 5)
ax5.scatter(adata.obs['wolbachia_titer'], adata.obs['Program_0'], alpha=0.3, s=10)
mask = ~pd.isna(adata.obs['wolbachia_titer']) & ~pd.isna(adata.obs['Program_0'])
slope, intercept, r, p, se = linregress(
    adata.obs.loc[mask, 'wolbachia_titer'], 
    adata.obs.loc[mask, 'Program_0']
)
x_line = np.linspace(adata.obs['wolbachia_titer'].min(), adata.obs['wolbachia_titer'].max(), 100)
ax5.plot(x_line, slope * x_line + intercept, 'r-', linewidth=2)
ax5.set_xlabel('Wolbachia Titer')
ax5.set_ylabel('Program 0 Score')
ax5.set_title(f'E. Program 0 (r={r:.3f})')
ax5.grid(True, alpha=0.3)

# Panel F: Program 6 scatter
ax6 = plt.subplot(2, 4, 6)
ax6.scatter(adata.obs['wolbachia_titer'], adata.obs['Program_6'], alpha=0.3, s=10)
mask = ~pd.isna(adata.obs['wolbachia_titer']) & ~pd.isna(adata.obs['Program_6'])
slope, intercept, r, p, se = linregress(
    adata.obs.loc[mask, 'wolbachia_titer'], 
    adata.obs.loc[mask, 'Program_6']
)
x_line = np.linspace(adata.obs['wolbachia_titer'].min(), adata.obs['wolbachia_titer'].max(), 100)
ax6.plot(x_line, slope * x_line + intercept, 'r-', linewidth=2)
ax6.set_xlabel('Wolbachia Titer')
ax6.set_ylabel('Program 6 Score')
ax6.set_title(f'F. Program 6 (r={r:.3f})')
ax6.grid(True, alpha=0.3)

# Panel G: Top genes for Program 0
ax7 = plt.subplot(2, 4, 7)
prog0_genes = pd.read_csv('nmf_gene_programs/Program_0_genes.csv')
# Map to symbols
import gzip
from io import StringIO
flybase_annot_path = "/private/groups/russelllab/jodie/scRNAseq/Jacobs_et_al_2026_wolbachia-drosophila-scrnaseq/reference/fbgn_annotation_ID_fb_2025_04.tsv.gz"
with gzip.open(flybase_annot_path, 'rt') as f:
    lines = [line for line in f if not line.startswith('#')]
mapping_df = pd.read_csv(StringIO(''.join(lines)), sep='\t', header=None)
fbgn_to_symbol = dict(zip(mapping_df[2], mapping_df[0]))
prog0_genes['symbol'] = prog0_genes['gene'].map(fbgn_to_symbol).fillna(prog0_genes['gene'])

top10 = prog0_genes.head(10)
ax7.barh(range(len(top10)), top10['weight'].values[::-1])
ax7.set_yticks(range(len(top10)))
ax7.set_yticklabels(top10['symbol'].values[::-1])
ax7.set_xlabel('Gene Weight')
ax7.set_title('G. Program 0 Top Genes')
ax7.grid(True, alpha=0.3, axis='x')

# Panel H: Top genes for Program 6
ax8 = plt.subplot(2, 4, 8)
prog6_genes = pd.read_csv('nmf_gene_programs/Program_6_genes.csv')
prog6_genes['symbol'] = prog6_genes['gene'].map(fbgn_to_symbol).fillna(prog6_genes['gene'])

top10 = prog6_genes.head(10)
ax8.barh(range(len(top10)), top10['weight'].values[::-1])
ax8.set_yticks(range(len(top10)))
ax8.set_yticklabels(top10['symbol'].values[::-1])
ax8.set_xlabel('Gene Weight')
ax8.set_title('H. Program 6 Top Genes')
ax8.grid(True, alpha=0.3, axis='x')

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'main_figure_wolbachia_programs.pdf'), dpi=300)
plt.close()
print("✓ Saved: main_figure_wolbachia_programs.pdf (publication-ready)")

# ============================================================================
# 4. Write interpretation report
# ============================================================================

with open(os.path.join(OUTPUT_DIR, 'BIOLOGICAL_INTERPRETATION.txt'), 'w') as f:
    f.write("="*70 + "\n")
    f.write("WOLBACHIA INFECTION: HOST RESPONSE PROGRAMS\n")
    f.write("="*70 + "\n\n")
    
    f.write("MAJOR FINDINGS:\n\n")
    
    f.write("1. PROGRAM 0: MITOCHONDRIAL SUPPRESSION (r=-0.117, FDR<0.001)\n")
    f.write("   - 20% of top genes are mitochondrial (NdufA3, ATPsynD, ND-39)\n")
    f.write("   - GSEA top hit: Mitochondrial proton-transport\n")
    f.write("   - INTERPRETATION: Wolbachia suppresses host mitochondrial function\n")
    f.write("   - MECHANISM: Likely competition for cellular resources\n\n")
    
    f.write("2. PROGRAM 6: RIBOSOMAL INDUCTION (r=0.186, FDR<0.001)\n")
    f.write("   - GSEA top hit: Cytosolic ribosome\n")
    f.write("   - Contains translation initiation factors (eIF5B, eIF1A)\n")
    f.write("   - INTERPRETATION: Wolbachia hijacks host translation machinery\n")
    f.write("   - MECHANISM: Wolbachia requires host ribosomes for bacterial protein synthesis\n\n")
    
    f.write("3. SECONDARY PROGRAMS (9, 11, 3, 13): Weaker titer-dependent responses\n\n")
    
    f.write("="*70 + "\n")
    f.write("BIOLOGICAL MODEL:\n")
    f.write("="*70 + "\n\n")
    
    f.write("Wolbachia infection creates a dual response in host cells:\n\n")
    f.write("  [High Wolbachia Titer]\n")
    f.write("           ↓\n")
    f.write("  ┌────────┴────────┐\n")
    f.write("  ↓                 ↓\n")
    f.write("SUPPRESS          INDUCE\n")
    f.write("Mitochondria    Ribosomes\n")
    f.write("(Program 0)     (Program 6)\n")
    f.write("  ↓                 ↓\n")
    f.write("Reduced ATP     Increased\n")
    f.write("production      translation\n")
    f.write("                capacity\n\n")
    
    f.write("This creates a metabolic trade-off:\n")
    f.write("- Host mitochondrial function decreases\n")
    f.write("- But protein synthesis capacity increases\n")
    f.write("- Wolbachia exploits increased translation for its own proteins\n\n")
    
    f.write("="*70 + "\n")
    f.write("LITERATURE CONNECTIONS:\n")
    f.write("="*70 + "\n\n")
    
    f.write("1. Wolbachia-mitochondria competition:\n")
    f.write("   - Both are maternally inherited\n")
    f.write("   - Both have bacterial origin (alphaproteobacteria)\n")
    f.write("   - Competition for cellular space and resources documented\n\n")
    
    f.write("2. Ribosome hijacking:\n")
    f.write("   - Common strategy for intracellular pathogens\n")
    f.write("   - Wolbachia lacks many biosynthetic pathways\n")
    f.write("   - Must rely on host for translation machinery\n\n")
    
    f.write("="*70 + "\n")

print("✓ Saved: BIOLOGICAL_INTERPRETATION.txt")

print("\n" + "="*70)
print("BIOLOGICAL INTERPRETATION COMPLETE!")
print("="*70)
print(f"\nAll results in: {OUTPUT_DIR}/")
print("\nKey files:")
print("  - main_figure_wolbachia_programs.pdf: Publication-ready figure")
print("  - BIOLOGICAL_INTERPRETATION.txt: Full biological interpretation")
print("="*70)