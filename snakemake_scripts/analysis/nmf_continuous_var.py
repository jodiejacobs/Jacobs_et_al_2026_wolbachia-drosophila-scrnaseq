#!/usr/bin/env python3
"""
Correlate NMF programs with continuous variable
"""

import argparse
import scanpy as sc
import pandas as pd
import numpy as np
from scipy.stats import spearmanr, pearsonr, false_discovery_control
import matplotlib.pyplot as plt
import seaborn as sns
import os

def parse_args():
    parser = argparse.ArgumentParser(description='NMF continuous variable analysis')
    parser.add_argument('--input', required=True, help='Input h5ad with programs')
    parser.add_argument('--output_dir', required=True, help='Output directory')
    parser.add_argument('--continuous_var', default=None, help='Continuous variable column')
    parser.add_argument('--flybase_annotation', default=None, help='FlyBase annotation')
    return parser.parse_args()

def detect_continuous_var(adata):
    """Auto-detect continuous variable"""
    possible = ['wolbachia_titer', 'titer', 'infection_level', 'pseudotime', 'dpt_pseudotime']
    for col in possible:
        if col in adata.obs.columns:
            print(f"  ✓ Detected continuous variable: '{col}'")
            return col
    return None

def main():
    args = parse_args()
    
    print("="*70)
    print("NMF CONTINUOUS VARIABLE ANALYSIS")
    print("="*70)
    
    # Load data
    print("\n1. Loading data...")
    adata = sc.read_h5ad(args.input)
    
    # Detect continuous variable
    cont_var = args.continuous_var or detect_continuous_var(adata)
    if cont_var is None:
        numeric_cols = adata.obs.select_dtypes(include=[np.number]).columns.tolist()
        print(f"\nERROR: No continuous variable found!")
        print(f"Available numeric columns: {numeric_cols}")
        return
    
    print(f"\n   Using: '{cont_var}'")
    print(f"   Range: {adata.obs[cont_var].min():.2f} - {adata.obs[cont_var].max():.2f}")
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Get programs
    program_cols = [c for c in adata.obs.columns if c.startswith('Program_')]
    print(f"   Found {len(program_cols)} programs")
    
    # Correlate
    print("\n2. Computing correlations...")
    results = []
    
    for prog in program_cols:
        valid = ~(pd.isna(adata.obs[prog]) | pd.isna(adata.obs[cont_var]))
        if valid.sum() < 10:
            continue
        
        r_s, p_s = spearmanr(adata.obs.loc[valid, cont_var], adata.obs.loc[valid, prog])
        r_p, p_p = pearsonr(adata.obs.loc[valid, cont_var], adata.obs.loc[valid, prog])
        
        results.append({
            'Program': prog,
            'Spearman_r': r_s,
            'Spearman_p': p_s,
            'Pearson_r': r_p,
            'Pearson_p': p_p,
            'N_cells': valid.sum()
        })
    
    corr_df = pd.DataFrame(results)
    corr_df['Spearman_FDR'] = false_discovery_control(corr_df['Spearman_p'])
    corr_df['Pearson_FDR'] = false_discovery_control(corr_df['Pearson_p'])
    corr_df = corr_df.sort_values('Spearman_r', key=abs, ascending=False)
    
    # Save
    corr_df.to_csv(os.path.join(args.output_dir, 'program_correlations.csv'), index=False)
    
    # Report
    print("\n3. Significant correlations (FDR < 0.05):")
    sig = corr_df[corr_df['Spearman_FDR'] < 0.05]
    for _, row in sig.iterrows():
        direction = "↑" if row['Spearman_r'] > 0 else "↓"
        print(f"   {row['Program']}: {direction} r={row['Spearman_r']:.3f}, FDR={row['Spearman_FDR']:.2e}")
    
    # Plot
    fig, ax = plt.subplots(figsize=(10, 8))
    colors = ['red' if r < 0 else 'blue' for r in corr_df['Spearman_r']]
    ax.barh(range(len(corr_df)), corr_df['Spearman_r'], color=colors, alpha=0.7)
    ax.set_yticks(range(len(corr_df)))
    ax.set_yticklabels(corr_df['Program'])
    ax.set_xlabel(f'Correlation with {cont_var}')
    ax.axvline(0, color='black', linewidth=1)
    
    for i, row in corr_df.iterrows():
        if row['Spearman_FDR'] < 0.05:
            ax.text(row['Spearman_r'], i, ' *', fontsize=12, va='center')
    
    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, 'correlations.pdf'), dpi=300)
    
    # Summary
    with open(os.path.join(args.output_dir, 'SUMMARY.txt'), 'w') as f:
        f.write(f"CONTINUOUS VARIABLE ANALYSIS\n")
        f.write(f"{'='*70}\n\n")
        f.write(f"Variable: {cont_var}\n")
        f.write(f"Range: {adata.obs[cont_var].min():.2f} - {adata.obs[cont_var].max():.2f}\n")
        f.write(f"Programs analyzed: {len(program_cols)}\n")
        f.write(f"Significant (FDR<0.05): {len(sig)}\n\n")
        for _, row in sig.iterrows():
            f.write(f"{row['Program']}: r={row['Spearman_r']:.3f}, FDR={row['Spearman_FDR']:.2e}\n")
    
    print(f"\n✓ Analysis complete!")
    print(f"Results in: {args.output_dir}")

if __name__ == '__main__':
    main()
