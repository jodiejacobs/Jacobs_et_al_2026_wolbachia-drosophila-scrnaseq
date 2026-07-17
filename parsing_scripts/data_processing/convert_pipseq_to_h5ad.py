#!/usr/bin/env python3

import argparse
import pandas as pd
import h5py
import scanpy as sc
import numpy as np
import anndata
from pathlib import Path
import gc
from tqdm import tqdm
import scipy.sparse as sp

def parse_args():
    parser = argparse.ArgumentParser(description='Convert PIPseq molecule_info.h5 to h5ad format')
    parser.add_argument('--input', '-i', required=True, help='Input molecule_info.h5 file from PIPseq')
    parser.add_argument('--output', '-o', required=True, help='Output .h5ad file')
    parser.add_argument('--min_cells', type=int, default=3, 
                      help='Minimum number of cells expressing a gene for it to be included')
    parser.add_argument('--min_genes', type=int, default=200, 
                      help='Minimum number of genes detected in a cell for it to be included')
    parser.add_argument('--compression', default='gzip', 
                      choices=['gzip', 'lzf', None], 
                      help='Compression for h5ad output')
    parser.add_argument('--compression_opts', type=int, default=4,
                      help='Compression level (1-9, higher=slower but smaller file)')
    parser.add_argument('--verbose', action='store_true', help='Print detailed information')
    return parser.parse_args()

def safe_print_h5_keys(h5_file):
    """Safely print h5 file keys without causing errors"""
    try:
        print(f"Keys in h5 file: {list(h5_file.keys())}")
        for key in h5_file.keys():
            try:
                item = h5_file[key]
                if isinstance(item, h5py.Dataset):
                    shape = item.shape
                    dtype = item.dtype
                    shape_str = "scalar" if len(shape) == 0 else str(shape)
                    print(f"  {key}: Dataset {shape_str}, {dtype}")
                elif isinstance(item, h5py.Group):
                    print(f"  {key}: Group with keys {list(item.keys())}")
            except Exception as e:
                print(f"  Error accessing {key}: {str(e)}")
    except Exception as e:
        print(f"Error exploring h5 file: {str(e)}")

def convert_pipseq_to_h5ad(input_file, output_file, min_cells=3, min_genes=200, 
                          compression='gzip', compression_opts=4, verbose=False):
    """
    Convert PIPseq molecule_info.h5 to h5ad format
    This handles the specific case of molecule-level data where each entry 
    represents a single RNA molecule (UMI)
    """
    print(f"Converting {input_file} to {output_file}")
    
    # Process the file
    with h5py.File(input_file, 'r') as f:
        # First, safely explore what's in the file
        safe_print_h5_keys(f)
        
        # Based on the previous error, we know some keys exist:
        # - barcode_list
        # - barcodes
        # - counts
        
        print("\nReading molecule data...")
        
        # Load barcodes and gene/feature information
        if 'barcode_list' in f:
            try:
                barcode_list = f['barcode_list'][:]
                if isinstance(barcode_list[0], bytes):
                    barcode_list = [b.decode('utf-8') for b in barcode_list]
                else:
                    barcode_list = [str(b) for b in barcode_list]
                print(f"Loaded {len(barcode_list)} unique cell barcodes")
            except Exception as e:
                print(f"Error loading barcode_list: {str(e)}")
                barcode_list = None
                
        # Look for feature/gene list
        gene_list = None
        for key in ['feature_list', 'gene_list', 'features', 'genes']:
            if key in f:
                try:
                    gene_list = f[key][:]
                    if isinstance(gene_list[0], bytes):
                        gene_list = [g.decode('utf-8') for g in gene_list]
                    else:
                        gene_list = [str(g) for g in gene_list]
                    print(f"Loaded {len(gene_list)} genes using key '{key}'")
                    break
                except Exception as e:
                    print(f"Error loading {key}: {str(e)}")
        
        # For molecule level data, we need to convert to a matrix
        molecule_barcodes = None
        if 'barcodes' in f:
            try:
                molecule_barcodes = f['barcodes'][:]
                print(f"Loaded {len(molecule_barcodes)} molecule barcode indices")
            except Exception as e:
                print(f"Error loading molecule barcodes: {str(e)}")
                
        molecule_features = None
        for key in ['features', 'genes']:
            if key in f:
                try:
                    molecule_features = f[key][:]
                    print(f"Loaded {len(molecule_features)} molecule feature indices using key '{key}'")
                    break
                except Exception as e:
                    print(f"Error loading molecule features from {key}: {str(e)}")
        
        # Get molecule counts (often implied as 1 per entry in molecule-level data)
        molecule_counts = None
        if 'counts' in f:
            try:
                molecule_counts = f['counts'][:]
                print(f"Loaded {len(molecule_counts)} molecule counts")
            except Exception as e:
                print(f"Error loading molecule counts: {str(e)}")
                # If no counts, assume each entry is 1 count
                if molecule_barcodes is not None:
                    print("Assuming 1 count per molecule")
                    molecule_counts = np.ones(len(molecule_barcodes), dtype=np.int32)
        
        # If we have all required molecule-level information, construct matrix
        if molecule_barcodes is not None and molecule_features is not None and molecule_counts is not None:
            print("Constructing sparse count matrix from molecule-level data...")
            
            # Determine unique cell and gene indices if not provided
            if barcode_list is None:
                unique_barcodes = np.unique(molecule_barcodes)
                n_cells = len(unique_barcodes)
                cell_idx_map = {bc: i for i, bc in enumerate(unique_barcodes)}
                print(f"Found {n_cells} unique cells from molecule data")
            else:
                n_cells = len(barcode_list)
                # Assume molecule_barcodes contains indices into barcode_list
                cell_idx_map = None  # Not needed if already indices
                print(f"Using {n_cells} cells from barcode list")
            
            if gene_list is None:
                unique_features = np.unique(molecule_features)
                n_genes = len(unique_features)
                gene_idx_map = {g: i for i, g in enumerate(unique_features)}
                print(f"Found {n_genes} unique genes from molecule data")
            else:
                n_genes = len(gene_list)
                # Assume molecule_features contains indices into gene_list
                gene_idx_map = None  # Not needed if already indices
                print(f"Using {n_genes} genes from gene list")
            
            # Create row and column indices for sparse matrix
            if cell_idx_map is not None:
                row_indices = np.array([cell_idx_map[bc] for bc in molecule_barcodes])
            else:
                row_indices = molecule_barcodes
                
            if gene_idx_map is not None:
                col_indices = np.array([gene_idx_map[g] for g in molecule_features])
            else:
                col_indices = molecule_features
            
            # Check indices are within bounds
            if np.max(row_indices) >= n_cells or np.min(row_indices) < 0:
                print(f"Warning: Cell indices out of bounds. Max: {np.max(row_indices)}, Min: {np.min(row_indices)}")
                n_cells = np.max(row_indices) + 1
                
            if np.max(col_indices) >= n_genes or np.min(col_indices) < 0:
                print(f"Warning: Gene indices out of bounds. Max: {np.max(col_indices)}, Min: {np.min(col_indices)}")
                n_genes = np.max(col_indices) + 1
            
            # Create sparse matrix in COO format, then convert to CSR
            print(f"Creating sparse matrix of shape ({n_cells}, {n_genes})...")
            matrix = sp.coo_matrix((molecule_counts, (row_indices, col_indices)), shape=(n_cells, n_genes))
            matrix = matrix.tocsr()
            
            # Free memory
            del molecule_barcodes, molecule_features, molecule_counts
            gc.collect()
            
            # Create cell and gene names if needed
            if barcode_list is None:
                barcode_list = [f"cell_{i}" for i in range(n_cells)]
                
            if gene_list is None:
                gene_list = [f"gene_{i}" for i in range(n_genes)]
                
            # Ensure we have the right number of names
            if len(barcode_list) < n_cells:
                print(f"Warning: Not enough cell names ({len(barcode_list)}) for matrix size ({n_cells})")
                barcode_list.extend([f"cell_{i+len(barcode_list)}" for i in range(n_cells - len(barcode_list))])
                
            if len(gene_list) < n_genes:
                print(f"Warning: Not enough gene names ({len(gene_list)}) for matrix size ({n_genes})")
                gene_list.extend([f"gene_{i+len(gene_list)}" for i in range(n_genes - len(gene_list))])
            
            # Create AnnData object
            print("Creating AnnData object...")
            adata = anndata.AnnData(X=matrix, 
                                  obs=pd.DataFrame(index=barcode_list[:n_cells]),
                                  var=pd.DataFrame(index=gene_list[:n_genes]))
            
            # Calculate basic QC metrics
            print("Calculating QC metrics...")
            adata.var['n_cells'] = np.array((adata.X > 0).sum(axis=0)).flatten()
            adata.obs['n_genes'] = np.array((adata.X > 0).sum(axis=1)).flatten()
            adata.obs['total_counts'] = np.array(adata.X.sum(axis=1)).flatten()
            
            # Filter data
            print("Filtering data...")
            print(f"Before filtering: {adata.shape[0]} cells, {adata.shape[1]} genes")
            # Filter genes
            sc.pp.filter_genes(adata, min_cells=min_cells)
            # Filter cells
            sc.pp.filter_cells(adata, min_genes=min_genes)
            
            print(f"After filtering: {adata.shape[0]} cells, {adata.shape[1]} genes")
            
            # Write to file with compression options
            print(f"Writing to {output_file}...")
            adata.write_h5ad(output_file, compression=compression, compression_opts=compression_opts)
            print("Conversion complete!")
            return adata
        else:
            print("Could not find all required data for matrix construction.")
            missing = []
            if molecule_barcodes is None:
                missing.append("cell indices")
            if molecule_features is None:
                missing.append("gene indices")
            if molecule_counts is None:
                missing.append("UMI counts")
            print(f"Missing information: {', '.join(missing)}")
            return None

def main():
    args = parse_args()
    convert_pipseq_to_h5ad(
        args.input, 
        args.output, 
        args.min_cells, 
        args.min_genes,
        args.compression,
        args.compression_opts,
        args.verbose
    )

if __name__ == "__main__":
    main()