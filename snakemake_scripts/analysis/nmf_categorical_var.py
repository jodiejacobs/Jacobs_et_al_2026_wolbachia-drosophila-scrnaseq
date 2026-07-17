#!/usr/bin/env python3
"""
Compare NMF programs across categorical groups
"""

import argparse
import scanpy as sc
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu, kruskal, false_discovery_control
import matplotlib.pyplot as plt
import os

def parse_args():
    parser = argparse.ArgumentParser(description='NMF categorical variable analysis')
    parser.add_argument('--input', required=True, help='Input h5ad with programs')
    parser.add_argument('--output_dir', required=True, help='Output directory')
    parser.add_argument('--categorical_var', default=None, help='Categorical variable column')
    return parser.parse_args()

def detect_categorical_var(adata):
    possible = ['treatment', 'condition', 'strain', 'genotype', 'cell_line', 'bio_condition']
    for col in possible:
        if col in adata.obs.columns:
            print(f"  ✓ Detected categorical variable: '{col}'")
            return col
    return None

def main():
    args = parse_args()
    
    print("="*70)
    print("NMF CATEGORICAL VARIABLE ANALYSIS")
    print("="*70)
    
    adata = sc.read_h5ad(args.input)
    
    cat_var = args.categorical_var or detect_categorical_var(adata)
    if cat_var is None:
        print("\nNo categorical variable found")
        return
    
    print(f"\n   Using: '{cat_var}'")
    categories = adata.obs[cat_var].unique()
    print(f"   Categories: {categories}")
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    program_cols = [c for c in adata.obs.columns if c.startswith('Program_')]
    
    results = []
    
    for prog in program_cols:
        groups = [adata.obs.loc[adata.obs[cat_var] == cat, prog].values 
                 for cat in categories]
        
        if len(categories) == 2:
            stat, p = mannwhitneyu(groups[0], groups[1])
            test = 'Mann-Whitney'
        else:
            stat, p = kruskal(*groups)
            test = 'Kruskal-Wallis'
        
        results.append({
            'Program': prog,
            'Test': test,
            'Statistic': stat,
            'P_value': p
        })
    
    results_df = pd.DataFrame(results)
    results_df['FDR'] = false_discovery_control(results_df['P_value'])
    results_df = results_df.sort_values('FDR')
    
    results_df.to_csv(os.path.join(args.output_dir, 'program_comparison.csv'), index=False)
    
    print("\nSignificant programs (FDR < 0.05):")
    for _, row in results_df[results_df['FDR'] < 0.05].iterrows():
        print(f"   {row['Program']}: p={row['P_value']:.2e}, FDR={row['FDR']:.2e}")
    
    print(f"\n✓ Analysis complete!")

if __name__ == '__main__':
    main()
