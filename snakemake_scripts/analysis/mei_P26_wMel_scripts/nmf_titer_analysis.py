import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr, pearsonr
from scipy.stats import false_discovery_control
import os

# Load the results
adata = sc.read_h5ad('nmf_gene_programs/adata_with_programs.h5ad')

OUTPUT_DIR = 'nmf_gene_programs/titer_analysis'
os.makedirs(OUTPUT_DIR, exist_ok=True)

print("="*60)
print("WOLBACHIA TITER-DEPENDENT PROGRAM ANALYSIS")
print("="*60)

# Check titer column
titer_col = 'wolbachia_titer'
if titer_col not in adata.obs.columns:
    print(f"\nERROR: '{titer_col}' not found in adata.obs")
    print(f"Available columns: {adata.obs.columns.tolist()}")
    exit(1)

# Get all program columns
program_cols = [col for col in adata.obs.columns if col.startswith('Program_')]
n_programs = len(program_cols)

print(f"\nAnalyzing {n_programs} programs")
print(f"Titer range: {adata.obs[titer_col].min():.2f} - {adata.obs[titer_col].max():.2f}")
print(f"Titer mean: {adata.obs[titer_col].mean():.2f}")
print(f"N cells: {len(adata)}")

# ============================================================================
# 1. Correlate each program with Wolbachia titer
# ============================================================================
print("\n1. Computing correlations with Wolbachia titer...")

correlation_results = []

for program in program_cols:
    program_scores = adata.obs[program]
    titer_values = adata.obs[titer_col]
    
    # Remove any NaN values
    valid_mask = ~(pd.isna(program_scores) | pd.isna(titer_values))
    
    if valid_mask.sum() < 10:
        print(f"   WARNING: Only {valid_mask.sum()} valid values for {program}")
        continue
    
    program_clean = program_scores[valid_mask]
    titer_clean = titer_values[valid_mask]
    
    # Calculate both Pearson and Spearman correlations
    pearson_r, pearson_p = pearsonr(titer_clean, program_clean)
    spearman_r, spearman_p = spearmanr(titer_clean, program_clean)
    
    correlation_results.append({
        'Program': program,
        'Pearson_r': pearson_r,
        'Pearson_pval': pearson_p,
        'Spearman_r': spearman_r,
        'Spearman_pval': spearman_p,
        'N_cells': valid_mask.sum()
    })

corr_df = pd.DataFrame(correlation_results)

# FDR correction
corr_df['Pearson_FDR'] = false_discovery_control(corr_df['Pearson_pval'])
corr_df['Spearman_FDR'] = false_discovery_control(corr_df['Spearman_pval'])

# Sort by absolute Spearman correlation
corr_df['Abs_Spearman_r'] = corr_df['Spearman_r'].abs()
corr_df = corr_df.sort_values('Abs_Spearman_r', ascending=False)

# Save results
corr_df.to_csv(os.path.join(OUTPUT_DIR, 'program_titer_correlations.csv'), index=False)

print("\n   Top positively correlated programs (Spearman):")
top_pos = corr_df[corr_df['Spearman_r'] > 0].head(5)
for idx, row in top_pos.iterrows():
    sig_star = '***' if row['Spearman_FDR'] < 0.001 else '**' if row['Spearman_FDR'] < 0.01 else '*' if row['Spearman_FDR'] < 0.05 else ''
    print(f"   {row['Program']}: r={row['Spearman_r']:.3f}, FDR={row['Spearman_FDR']:.3e} {sig_star}")

print("\n   Top negatively correlated programs (Spearman):")
top_neg = corr_df[corr_df['Spearman_r'] < 0].head(5)
for idx, row in top_neg.iterrows():
    sig_star = '***' if row['Spearman_FDR'] < 0.001 else '**' if row['Spearman_FDR'] < 0.01 else '*' if row['Spearman_FDR'] < 0.05 else ''
    print(f"   {row['Program']}: r={row['Spearman_r']:.3f}, FDR={row['Spearman_FDR']:.3e} {sig_star}")

# ============================================================================
# 2. Visualization: Correlation bar plot
# ============================================================================
print("\n2. Creating correlation visualizations...")

# Bar plot of correlations
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

# Spearman correlations
colors = ['red' if r < 0 else 'blue' for r in corr_df['Spearman_r']]
ax1.barh(range(len(corr_df)), corr_df['Spearman_r'], color=colors, alpha=0.7)
ax1.set_yticks(range(len(corr_df)))
ax1.set_yticklabels(corr_df['Program'])
ax1.set_xlabel('Spearman Correlation with Wolbachia Titer')
ax1.set_title('Program-Titer Correlations (Spearman)')
ax1.axvline(0, color='black', linestyle='-', linewidth=0.5)
ax1.grid(True, alpha=0.3)

# Add significance markers
for i, (idx, row) in enumerate(corr_df.iterrows()):
    if row['Spearman_FDR'] < 0.05:
        ax1.text(row['Spearman_r'] + 0.02 * np.sign(row['Spearman_r']), 
                i, '*', fontsize=12, va='center')

# Pearson correlations
colors = ['red' if r < 0 else 'blue' for r in corr_df['Pearson_r']]
ax2.barh(range(len(corr_df)), corr_df['Pearson_r'], color=colors, alpha=0.7)
ax2.set_yticks(range(len(corr_df)))
ax2.set_yticklabels(corr_df['Program'])
ax2.set_xlabel('Pearson Correlation with Wolbachia Titer')
ax2.set_title('Program-Titer Correlations (Pearson)')
ax2.axvline(0, color='black', linestyle='-', linewidth=0.5)
ax2.grid(True, alpha=0.3)

for i, (idx, row) in enumerate(corr_df.iterrows()):
    if row['Pearson_FDR'] < 0.05:
        ax2.text(row['Pearson_r'] + 0.02 * np.sign(row['Pearson_r']), 
                i, '*', fontsize=12, va='center')

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'correlation_barplot.pdf'))
plt.close()
print("   Correlation bar plots saved!")

# ============================================================================
# 3. Scatter plots for top correlated programs
# ============================================================================
print("\n3. Creating scatter plots for top correlated programs...")

# Get top 6 most strongly correlated (by absolute value)
top_6_programs = corr_df.head(6)['Program'].tolist()

fig, axes = plt.subplots(2, 3, figsize=(18, 12))
axes = axes.flatten()

for idx, program in enumerate(top_6_programs):
    ax = axes[idx]
    
    x = adata.obs[titer_col]
    y = adata.obs[program]
    
    # Scatter plot
    ax.scatter(x, y, alpha=0.5, s=20, edgecolors='none')
    
    # Add regression line
    valid_mask = ~(pd.isna(x) | pd.isna(y))
    if valid_mask.sum() > 2:
        z = np.polyfit(x[valid_mask], y[valid_mask], 1)
        p = np.poly1d(z)
        x_line = np.linspace(x[valid_mask].min(), x[valid_mask].max(), 100)
        ax.plot(x_line, p(x_line), "r-", linewidth=2, alpha=0.8)
    
    # Get correlation stats
    corr_row = corr_df[corr_df['Program'] == program].iloc[0]
    
    ax.set_xlabel('Wolbachia Titer')
    ax.set_ylabel('Program Score')
    ax.set_title(f"{program}\nr={corr_row['Spearman_r']:.3f}, FDR={corr_row['Spearman_FDR']:.2e}")
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'top_programs_scatter.pdf'))
plt.close()
print("   Scatter plots saved!")

# ============================================================================
# 4. UMAP with titer and top programs
# ============================================================================
print("\n4. Creating UMAP visualizations...")

if 'X_umap' in adata.obsm:
    # Create figure with UMAP colored by titer and top 4 programs
    top_4_programs = corr_df.head(4)['Program'].tolist()
    
    fig = plt.figure(figsize=(20, 4))
    
    # UMAP colored by titer
    ax1 = plt.subplot(1, 5, 1)
    sc1 = ax1.scatter(adata.obsm['X_umap'][:, 0],
                     adata.obsm['X_umap'][:, 1],
                     c=adata.obs[titer_col],
                     cmap='YlOrRd',
                     s=10,
                     alpha=0.7)
    ax1.set_title('Wolbachia Titer')
    ax1.set_xlabel('UMAP1')
    ax1.set_ylabel('UMAP2')
    plt.colorbar(sc1, ax=ax1)
    
    # UMAPs for top 4 correlated programs
    for idx, program in enumerate(top_4_programs):
        ax = plt.subplot(1, 5, idx+2)
        corr_row = corr_df[corr_df['Program'] == program].iloc[0]
        
        sc = ax.scatter(adata.obsm['X_umap'][:, 0],
                       adata.obsm['X_umap'][:, 1],
                       c=adata.obs[program],
                       cmap='viridis',
                       s=10,
                       alpha=0.7)
        ax.set_title(f"{program}\nr={corr_row['Spearman_r']:.2f}")
        ax.set_xlabel('UMAP1')
        ax.set_ylabel('UMAP2')
        plt.colorbar(sc, ax=ax)
    
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'umap_titer_and_programs.pdf'))
    plt.close()
    print("   UMAP plots saved!")

# ============================================================================
# 5. Bin titer and compare programs
# ============================================================================
print("\n5. Binning titer and comparing programs...")

# Create titer bins (quartiles)
adata.obs['titer_quartile'] = pd.qcut(adata.obs[titer_col], 
                                       q=4, 
                                       labels=['Q1_Low', 'Q2', 'Q3', 'Q4_High'],
                                       duplicates='drop')

# Compare program scores across quartiles
quartile_comparison = []

for program in program_cols:
    for quartile in ['Q1_Low', 'Q2', 'Q3', 'Q4_High']:
        mask = adata.obs['titer_quartile'] == quartile
        if mask.sum() > 0:
            quartile_comparison.append({
                'Program': program,
                'Quartile': quartile,
                'Mean_score': adata.obs.loc[mask, program].mean(),
                'Std_score': adata.obs.loc[mask, program].std(),
                'N_cells': mask.sum()
            })

quartile_df = pd.DataFrame(quartile_comparison)
quartile_df.to_csv(os.path.join(OUTPUT_DIR, 'program_by_titer_quartile.csv'), index=False)

# Violin plots for top programs across quartiles
top_6_for_violin = corr_df.head(6)['Program'].tolist()

fig, axes = plt.subplots(2, 3, figsize=(18, 12))
axes = axes.flatten()

for idx, program in enumerate(top_6_for_violin):
    ax = axes[idx]
    
    # Prepare data for violin plot
    data_by_quartile = []
    quartiles = ['Q1_Low', 'Q2', 'Q3', 'Q4_High']
    
    for q in quartiles:
        mask = adata.obs['titer_quartile'] == q
        if mask.sum() > 0:
            data_by_quartile.append(adata.obs.loc[mask, program].values)
        else:
            data_by_quartile.append([])
    
    # Create violin plot
    parts = ax.violinplot(data_by_quartile, 
                          positions=range(len(quartiles)),
                          showmeans=True,
                          showmedians=True)
    
    # Add box plot overlay
    ax.boxplot(data_by_quartile,
               positions=range(len(quartiles)),
               widths=0.1,
               showfliers=False)
    
    ax.set_xticks(range(len(quartiles)))
    ax.set_xticklabels(quartiles, rotation=45, ha='right')
    ax.set_ylabel('Program Score')
    
    corr_row = corr_df[corr_df['Program'] == program].iloc[0]
    ax.set_title(f"{program}\nr={corr_row['Spearman_r']:.3f}")
    ax.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'programs_by_titer_quartile.pdf'))
plt.close()
print("   Quartile violin plots saved!")

# ============================================================================
# 6. High vs Low titer comparison
# ============================================================================
print("\n6. Comparing high vs low titer cells...")

# Define high and low titer (top and bottom 25%)
titer_25 = adata.obs[titer_col].quantile(0.25)
titer_75 = adata.obs[titer_col].quantile(0.75)

high_titer = adata.obs[titer_col] >= titer_75
low_titer = adata.obs[titer_col] <= titer_25

print(f"   High titer (≥75th percentile): {high_titer.sum()} cells, titer ≥ {titer_75:.2f}")
print(f"   Low titer (≤25th percentile): {low_titer.sum()} cells, titer ≤ {titer_25:.2f}")

from scipy.stats import mannwhitneyu

high_low_comparison = []

for program in program_cols:
    high_scores = adata.obs.loc[high_titer, program]
    low_scores = adata.obs.loc[low_titer, program]
    
    stat, pval = mannwhitneyu(high_scores, low_scores, alternative='two-sided')
    
    fold_change = high_scores.mean() / low_scores.mean() if low_scores.mean() > 0 else float('inf')
    
    high_low_comparison.append({
        'Program': program,
        'High_titer_mean': high_scores.mean(),
        'Low_titer_mean': low_scores.mean(),
        'Fold_change': fold_change,
        'Log2_FC': np.log2(fold_change) if fold_change > 0 else np.nan,
        'P_value': pval
    })

high_low_df = pd.DataFrame(high_low_comparison)
high_low_df['FDR'] = false_discovery_control(high_low_df['P_value'])
high_low_df = high_low_df.sort_values('FDR')

high_low_df.to_csv(os.path.join(OUTPUT_DIR, 'high_vs_low_titer_programs.csv'), index=False)

print("\n   Programs significantly different between high/low titer (FDR < 0.05):")
sig_high_low = high_low_df[high_low_df['FDR'] < 0.05]
for idx, row in sig_high_low.head(10).iterrows():
    direction = "UP" if row['Fold_change'] > 1 else "DOWN"
    print(f"   {row['Program']}: {direction} {row['Fold_change']:.2f}x, FDR={row['FDR']:.3e}")

# ============================================================================
# 7. Identify titer-responsive vs titer-independent programs
# ============================================================================
print("\n7. Categorizing programs by titer responsiveness...")

# Classify programs
titer_responsive = corr_df[corr_df['Spearman_FDR'] < 0.05]['Program'].tolist()
titer_independent = corr_df[corr_df['Spearman_FDR'] >= 0.05]['Program'].tolist()

print(f"\n   Titer-responsive programs (FDR < 0.05): {len(titer_responsive)}")
print(f"   Titer-independent programs: {len(titer_independent)}")

classification = pd.DataFrame({
    'Program': program_cols,
    'Category': ['Titer-responsive' if p in titer_responsive else 'Titer-independent' 
                 for p in program_cols]
})
classification.to_csv(os.path.join(OUTPUT_DIR, 'program_classification.csv'), index=False)

# ============================================================================
# 8. Summary Report
# ============================================================================
print("\n8. Creating summary report...")

with open(os.path.join(OUTPUT_DIR, 'TITER_ANALYSIS_SUMMARY.txt'), 'w') as f:
    f.write("="*60 + "\n")
    f.write("WOLBACHIA TITER-DEPENDENT PROGRAM ANALYSIS\n")
    f.write("="*60 + "\n\n")
    
    f.write(f"Dataset: {len(adata)} cells\n")
    f.write(f"Titer range: {adata.obs[titer_col].min():.2f} - {adata.obs[titer_col].max():.2f}\n")
    f.write(f"Titer mean ± std: {adata.obs[titer_col].mean():.2f} ± {adata.obs[titer_col].std():.2f}\n")
    f.write(f"Number of programs analyzed: {n_programs}\n\n")
    
    f.write("TITER-RESPONSIVE PROGRAMS (FDR < 0.05):\n")
    f.write("-"*60 + "\n\n")
    
    for idx, row in corr_df[corr_df['Spearman_FDR'] < 0.05].iterrows():
        direction = "Positive" if row['Spearman_r'] > 0 else "Negative"
        f.write(f"{row['Program']}:\n")
        f.write(f"  Direction: {direction} correlation\n")
        f.write(f"  Spearman r: {row['Spearman_r']:.3f}\n")
        f.write(f"  FDR: {row['Spearman_FDR']:.3e}\n\n")
    
    f.write("\nTITER-INDEPENDENT PROGRAMS:\n")
    f.write("-"*60 + "\n\n")
    
    for idx, row in corr_df[corr_df['Spearman_FDR'] >= 0.05].iterrows():
        f.write(f"{row['Program']}: r={row['Spearman_r']:.3f}, FDR={row['Spearman_FDR']:.3e}\n")

print("\n" + "="*60)
print("TITER ANALYSIS COMPLETE!")
print("="*60)
print(f"\nResults saved to: {OUTPUT_DIR}/")
print("\nKey outputs:")
print("  - program_titer_correlations.csv: Correlation statistics")
print("  - high_vs_low_titer_programs.csv: Differential program usage")
print("  - correlation_barplot.pdf: Overview of all correlations")
print("  - top_programs_scatter.pdf: Detailed scatter plots")
print("  - umap_titer_and_programs.pdf: Spatial patterns")
print("  - programs_by_titer_quartile.pdf: Dose-response patterns")
print("  - TITER_ANALYSIS_SUMMARY.txt: Summary report")
print("="*60)