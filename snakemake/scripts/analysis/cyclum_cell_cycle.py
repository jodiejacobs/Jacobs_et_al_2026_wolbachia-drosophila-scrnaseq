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
from sklearn.mixture import GaussianMixture
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Cyclum Cell Cycle Analysis")
parser.add_argument("--input", type=str, required=False, help="Input h5ad file", default='/private/groups/russelllab/jodie/scRNAseq/scripts/snakemake_pipeline/results_kallisto_bustools/filtered_h5ad/kallisto_JW18DOX-Ctrl-1_P.h5ad')
parser.add_argument("--output", type=str, required=False, help="Output h5ad file", default='/private/groups/russelllab/jodie/scRNAseq/scripts/snakemake_pipeline/results_kallisto_bustools/filtered_h5ad/cyclum_kallisto_JW18DOX-Ctrl-1_P.h5ad')

args = parser.parse_args()

adata = sc.read_h5ad(args.input)
output = args.output
mtx = adata.X

# Train model 
model = cyclum.tuning.CyclumAutoTune(mtx)

# Train with more epochs for better cell detection
model.train(mtx, epochs=800, verbose=100, rate=2e-4)

# Extract the circular pseudotime (this represents cell cycle phase)
pseudotime = model.predict_pseudotime(mtx)
pseudotime_flat = pseudotime.flatten()

# Check the pseudotime shape and range
print(f"Pseudotime shape: {pseudotime.shape}")
print(f"Pseudotime range: {pseudotime.min():.3f} to {pseudotime.max():.3f}")

def assign_cell_cycle_stage_advanced(pseudotime_flat):
    """
    Assign cell cycle stages using GMM followed by validation and refinement
    """
    print("Step 1: Initial assignment using Gaussian Mixture Model...")
    
    # Convert pseudotime to 2D circular coordinates for proper clustering
    angles = pseudotime_flat * 2 * np.pi if pseudotime_flat.max() <= 1 else pseudotime_flat
    x = np.cos(angles)
    y = np.sin(angles)
    circular_coords = np.column_stack([x, y])
    
    # Initial GMM clustering with 3 components
    gmm = GaussianMixture(n_components=3, random_state=42, max_iter=200)
    initial_clusters = gmm.fit_predict(circular_coords)
    cluster_probs = gmm.predict_proba(circular_coords)
    
    print(f"GMM converged: {gmm.converged_}")
    
    # Map clusters to phases based on mean pseudotime ordering
    cluster_means = []
    for i in range(3):
        cluster_mask = initial_clusters == i
        mean_pseudotime = np.mean(pseudotime_flat[cluster_mask])
        cluster_means.append((i, mean_pseudotime))
    
    # Sort by mean pseudotime to establish biological order
    cluster_means.sort(key=lambda x: x[1])
    
    # Create mapping: lowest mean = G1, middle = S, highest = G2/M
    phase_mapping = {
        cluster_means[0][0]: 'g0/g1',  # Lowest pseudotime
        cluster_means[1][0]: 's',       # Middle pseudotime  
        cluster_means[2][0]: 'g2/m'     # Highest pseudotime
    }
    
    initial_phases = [phase_mapping[cluster] for cluster in initial_clusters]
    
    print("Step 2: Calculating confidence scores...")
    
    # Calculate confidence scores based on cluster probabilities
    max_probs = np.max(cluster_probs, axis=1)
    confidence_threshold = np.percentile(max_probs, 25)  # Bottom 25% are low confidence
    
    print("Step 3: Refining assignments for low-confidence cells...")
    
    # Refine low-confidence assignments using local neighborhood
    refined_phases = initial_phases.copy()
    low_conf_indices = np.where(max_probs < confidence_threshold)[0]
    
    print(f"Refining {len(low_conf_indices)} low-confidence assignments...")
    
    for idx in low_conf_indices:
        # Find nearest neighbors in pseudotime space
        distances = np.abs(pseudotime_flat - pseudotime_flat[idx])
        nearest_indices = np.argsort(distances)[1:11]  # 10 nearest neighbors (excluding self)
        
        # Get high-confidence neighbors
        high_conf_neighbors = [i for i in nearest_indices if max_probs[i] >= confidence_threshold]
        
        if len(high_conf_neighbors) >= 3:  # Need at least 3 confident neighbors
            neighbor_phases = [initial_phases[i] for i in high_conf_neighbors[:5]]  # Use top 5
            # Assign most common phase among confident neighbors
            most_common = max(set(neighbor_phases), key=neighbor_phases.count)
            refined_phases[idx] = most_common
    
    print("Step 4: Final validation...")
    
    # Final validation: ensure reasonable phase distribution
    phase_counts = pd.Series(refined_phases).value_counts()
    total_cells = len(refined_phases)
    
    print("Phase distribution after refinement:")
    for phase in ['g0/g1', 's', 'g2/m']:
        count = phase_counts.get(phase, 0)
        percentage = (count / total_cells) * 100
        print(f"  {phase}: {count} cells ({percentage:.1f}%)")
    
    # Check if any phase is severely under-represented (< 15%)
    min_expected = 0.15 * total_cells
    for phase in ['g0/g1', 's', 'g2/m']:
        if phase_counts.get(phase, 0) < min_expected:
            print(f"Warning: {phase} phase has fewer cells than expected (< 15%)")
    
    return refined_phases, max_probs

# Assign cell cycle stages using advanced method
stages, confidence_scores = assign_cell_cycle_stage_advanced(pseudotime_flat)

# Create a label dictionary like in the tutorial
label = {'stage': np.array(stages)}

# Check the distribution
unique_stages, counts = np.unique(stages, return_counts=True)
print("\nFinal cell cycle stage distribution:")
for stage, count in zip(unique_stages, counts):
    percentage = (count / len(stages)) * 100
    print(f"{stage}: {count} cells ({percentage:.1f}%)")

# Add to adata
adata.obs['cyclum_stage'] = stages
adata.obs['cyclum_pseudotime'] = pseudotime_flat
adata.obs['cyclum_confidence'] = confidence_scores

print("\nAdded to adata.obs:")
print(f"cyclum_stage: {len(adata.obs['cyclum_stage'])} cells")
print(f"cyclum_pseudotime: {len(adata.obs['cyclum_pseudotime'])} cells")
print(f"cyclum_confidence: {len(adata.obs['cyclum_confidence'])} cells")

# Define color map (exactly like tutorial)
color_map = {'stage': {"g0/g1": "red", "s": "green", "g2/m": "blue"}}

# Create the circular cell cycle plot
fig = cyclum.illustration.plot_round_distr_color(pseudotime_flat, label['stage'], color_map['stage'])
plt.savefig(output.replace('.h5ad', '_cyclum_cell_cycle.pdf'), dpi=300, bbox_inches='tight')
plt.close()

# Show elbow plot
elbow_fig = model.show_elbow()
plt.savefig(output.replace('.h5ad', '_cyclum_elbow.pdf'), dpi=300, bbox_inches='tight')
plt.close()

# Show bar plot
bar_fig = model.show_bar()
plt.savefig(output.replace('.h5ad', '_cyclum_bar.pdf'), dpi=300, bbox_inches='tight')
plt.close()

# Additional plot: confidence scores
plt.figure(figsize=(10, 6))
plt.subplot(1, 2, 1)
plt.hist(confidence_scores, bins=50, alpha=0.7, edgecolor='black')
plt.xlabel('Confidence Score')
plt.ylabel('Number of Cells')
plt.title('Assignment Confidence Distribution')

plt.subplot(1, 2, 2)
stage_conf = pd.DataFrame({'stage': stages, 'confidence': confidence_scores})
for stage in ['g0/g1', 's', 'g2/m']:
    stage_data = stage_conf[stage_conf['stage'] == stage]['confidence']
    plt.hist(stage_data, alpha=0.7, label=stage, bins=30)
plt.xlabel('Confidence Score')
plt.ylabel('Number of Cells')
plt.title('Confidence by Cell Cycle Stage')
plt.legend()

plt.tight_layout()
plt.savefig(output.replace('.h5ad', '_cyclum_confidence.pdf'), dpi=300, bbox_inches='tight')
plt.close()

# Save the adata object
adata.write_h5ad(output)

print(f"\nAnalysis complete! Results saved to: {output}")
print(f"Plots saved as: {output.replace('.h5ad', '_cyclum_*.pdf')}")