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
import os

def parse_args():
    parser = argparse.ArgumentParser(description='Convert 10X .h5 to h5ad format')
    parser.add_argument('--input', '-i', required=True, help='Input .h5 file from cellranger')
    parser.add_argument('--output', '-o', required=True, help='Output .h5ad file')
    parser.add_argument('--sample', '-s', required=True, help='Sample name')   
    return parser.parse_args()

def convert_10x_to_h5ad(input, sample, output):
    # Load the Cell Ranger h5 file
    adata = sc.read_10x_h5(input)

    # Add some metadata
    adata.obs['sample'] = sample

    # Save as h5ad
    adata.write(output)
    print(f'Converted {os.path.basename(input)} to h5ad format')


def main():
    args = parse_args()
    convert_10x_to_h5ad(
        args.input,
        args.sample, 
        args.output ) 

if __name__ == "__main__":
    main()