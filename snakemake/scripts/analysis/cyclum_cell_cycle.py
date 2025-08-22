'''
This script annotates cell cycle stages using Cyclum. 
    Input: Filtered h5ad file
    Output: Cell cycle plots and annotated h5ad file
'''
import cyclum 
import cyclum.models
import cyclum.tuning 
import cyclum.illustration
import scanpy as sc 
import argparse
import os
import sklearn
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Cyclum Cell Cycle Analysis")
parser.add_argument("--input", type=str, required=False, help="Input h5ad file", default='/private/groups/russelllab/jodie/scRNAseq/scripts/snakemake_pipeline/results_kallisto_bustools/filtered_h5ad/kallisto_JW18DOX-Ctrl-1_P.h5ad')
parser.add_argument("--output", type=str, required=False, help="Output h5ad file", default='/private/groups/russelllab/jodie/scRNAseq/scripts/snakemake_pipeline/results_kallisto_bustools/filtered_h5ad/cyclum_kallisto_JW18DOX-Ctrl-1_P.h5ad')

args = parser.parse_args()

adata = sc.read_h5ad(args.input)
output = args.output
mtx=adata.X

# Train model 
model = cyclum.tuning.CyclumAutoTune(mtx)

# Train with more epochs for better cell detection
model.train(mtx, epochs=800, verbose=100, rate=2e-4)

# Extract the circular pseudotime (this represents cell cycle phase)
pseudotime = model.predict_pseudotime(mtx)

# Flatten pseudotime and create labels
pseudotime_flat = pseudotime.flatten()

# Check the pseudotime shape and range
print(f"Pseudotime shape: {pseudotime.shape}")
print(f"Pseudotime range: {pseudotime.min():.3f} to {pseudotime.max():.3f}")

def assign_cell_cycle_stage(pseudotime_flat):
    """Assign cell cycle stages based on pseudotime distribution"""
    # Sort to find balanced boundaries
    sorted_pt = np.sort(pseudotime_flat)
    n = len(sorted_pt)
    
    # Use approximate thirds for G1, S, G2/M phases
    boundary1 = sorted_pt[n//3]
    boundary2 = sorted_pt[2*n//3]
    
    stages = []
    for pt in pseudotime_flat:
        if pt <= boundary1:
            stages.append('g0/g1')
        elif pt <= boundary2:
            stages.append('s')
        else:
            stages.append('g2/m')
    
    return stages

# Flatten pseudotime and create labels
pseudotime_flat = pseudotime.flatten()
stages = assign_cell_cycle_stage(pseudotime_flat)

# Create a label dictionary like in the tutorial
label = {'stage': np.array(stages)}

# Check the distribution
unique_stages, counts = np.unique(stages, return_counts=True)
print("Cell cycle stage distribution:")
for stage, count in zip(unique_stages, counts):
    print(f"{stage}: {count} cells")

adata.obs['cyclum_stage'] = stages
adata.obs['cyclum_pseudotime'] = pseudotime_flat

print("\nAdded to adata.obs:")
print(f"cyclum_stage: {len(adata.obs['cyclum_stage'])} cells")
print(f"cyclum_pseudotime: {len(adata.obs['cyclum_pseudotime'])} cells")

# Define color map (exactly like tutorial)
color_map = {'stage': {"g0/g1": "red", "s": "green", "g2/m": "blue"},
             'subcluster': {"intact": "cyan", "perturbed": "violet"}}

# Create the circular cycle cycle plot
fig = cyclum.illustration.plot_round_distr_color(pseudotime_flat, label['stage'], color_map['stage'])
plt.savefig(output.replace('.h5ad', '_cyclum_cell_cycle.pdf'), dpi=300)

# Show elbow plot
elbow_fig = model.show_elbow()
plt.savefig(output.replace('.h5ad', '_cyclum_elbow.pdf'), dpi=300)

# Show bar plot
bar_fig = model.show_bar()
plt.savefig(output.replace('.h5ad', '_cyclum_bar.pdf'), dpi=300)

# Save the adata object:
adata.write_h5ad(output)