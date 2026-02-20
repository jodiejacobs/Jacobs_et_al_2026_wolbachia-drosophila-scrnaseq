#!/usr/bin/env python3
"""
NMF program discovery for Snakemake workflow
"""

import argparse
import scanpy as sc
import pandas as pd
import numpy as np
import gseapy as gp
from sklearn.decomposition import NMF
import matplotlib.pyplot as plt
import os

def parse_args():
    parser = argparse.ArgumentParser(description='NMF program discovery')
    parser.add_argument('--input', required=True, help='Input h5ad file')
    parser.add_argument('--output_dir', required=True, help='Output directory')
    parser.add_argument('--n_programs', type=int, default=15, help='Number of programs')
    parser.add_argument('--n_top_genes', type=int, default=2000, help='Number of HVGs')
    parser.add_argument('--organism', default='Fly', help='Organism for gene sets')
    parser.add_argument('--gene_id_type', default='flybase', help='Gene ID type')
    parser.add_argument('--flybase_annotation', default=None, help='FlyBase annotation')
    return parser.parse_args()

def load_gene_mapping(flybase_path):
    import gzip
    from io import StringIO
    
    with gzip.open(flybase_path, 'rt') as f:
        lines = [line for line in f if not line.startswith('#')]
    
    mapping_df = pd.read_csv(StringIO(''.join(lines)), sep='\t', header=None)
    return dict(zip(mapping_df[2], mapping_df[0]))

def main():
    args = parse_args()
    
    print("="*70)
    print("NMF GENE PROGRAM DISCOVERY")
    print("="*70)
    
    # Load data
    print("\n1. Loading data...")
    adata = sc.read_h5ad(args.input)
    print(f"   {adata.n_obs} cells × {adata.n_vars} genes")
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Prepare HVGs
    print("\n2. Preparing highly variable genes...")
    if 'highly_variable' not in adata.var.columns:
        sc.pp.highly_variable_genes(adata, n_top_genes=args.n_top_genes, flavor='seurat_v3')
    
    adata_hvg = adata[:, adata.var.highly_variable].copy()
    print(f"   Using {adata_hvg.n_vars} HVGs")
    
    X = adata_hvg.X.toarray() if hasattr(adata_hvg.X, 'toarray') else adata_hvg.X
    X = np.maximum(X, 0)
    
    # Run NMF
    print(f"\n3. Running NMF with {args.n_programs} programs...")
    model = NMF(n_components=args.n_programs, init='nndsvda', random_state=42, 
                max_iter=1000, verbose=1)
    W = model.fit_transform(X)
    H = model.components_
    
    print(f"   Reconstruction error: {model.reconstruction_err_:.2f}")
    
    # Extract programs
    print("\n4. Extracting programs...")
    for i in range(args.n_programs):
        top_idx = np.argsort(-H[i, :])[:200]
        pd.DataFrame({
            'gene': adata_hvg.var_names[top_idx],
            'weight': H[i, top_idx],
            'rank': range(1, 201)
        }).to_csv(os.path.join(args.output_dir, f'Program_{i}_genes.csv'), index=False)
    
    # Add to adata
    for i in range(args.n_programs):
        adata.obs[f'Program_{i}'] = W[:, i]
    
    # Save
    adata.write(os.path.join(args.output_dir, 'adata_with_programs.h5ad'))
    
    # Summary
    with open(os.path.join(args.output_dir, 'SUMMARY.txt'), 'w') as f:
        f.write(f"NMF PROGRAM DISCOVERY\n")
        f.write(f"{'='*70}\n\n")
        f.write(f"Dataset: {adata.n_obs} cells × {adata.n_vars} genes\n")
        f.write(f"HVGs used: {adata_hvg.n_vars}\n")
        f.write(f"Programs: {args.n_programs}\n")
        f.write(f"Reconstruction error: {model.reconstruction_err_:.2f}\n")
    
    print("\n✓ NMF complete!")
    print(f"Results in: {args.output_dir}")

if __name__ == '__main__':
    main()
