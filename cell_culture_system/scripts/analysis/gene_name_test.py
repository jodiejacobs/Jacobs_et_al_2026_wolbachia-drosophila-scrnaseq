import scanpy as sc
import pandas as pd

# Quick script to examine gene names in your dataset
adata = sc.read_h5ad('snakemake/results/filtered_h5ad/JW18DOX-Ctrl-2_10x.h5ad')

print(f"Dataset: {adata.n_obs} cells, {adata.n_vars} genes")
print(f"\nFirst 20 gene names:")
for i, gene in enumerate(adata.var_names[:20]):
    print(f"{i+1:2d}. {gene}")

print(f"\nLast 20 gene names:")
for i, gene in enumerate(adata.var_names[-20:]):
    print(f"{len(adata.var_names)-19+i:2d}. {gene}")

# Look for patterns
print(f"\nGene name patterns:")
print(f"Genes starting with 'CG': {sum(1 for g in adata.var_names if g.startswith('CG'))}")
print(f"Genes with numbers: {sum(1 for g in adata.var_names if any(c.isdigit() for c in g))}")
print(f"Genes with hyphens: {sum(1 for g in adata.var_names if '-' in g)}")
print(f"Genes with underscores: {sum(1 for g in adata.var_names if '_' in g)}")

# Search for known cell cycle genes
known_cc = ['CycA', 'CycB', 'CycE', 'pcna', 'Mcm2', 'polo', 'string', 'stg']
found_cc = [g for g in known_cc if g in adata.var_names]
print(f"\nKnown cell cycle genes found: {found_cc}")

# Search for ribosomal patterns
ribo_patterns = ['RpL', 'RpS', 'rpl', 'rps']
ribo_found = []
for pattern in ribo_patterns:
    matches = [g for g in adata.var_names if pattern in g]
    if matches:
        ribo_found.extend(matches[:3])  # Show first 3 matches
print(f"Ribosomal genes found: {ribo_found[:10]}")  # Show first 10