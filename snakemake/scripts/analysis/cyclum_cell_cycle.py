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
from sklearn.neighbors import NearestNeighbors
from scipy import stats
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

def assign_cell_cycle_stage_robust(pseudotime_flat):
    """
    Robust cell cycle assignment using multiple methods with gap-filling
    """
    print("Step 1: Quantile-based initial assignment...")
    
    # Normalize pseudotime to [0, 2π] for circular analysis
    if pseudotime_flat.max() <= 1:
        angles = pseudotime_flat * 2 * np.pi
    else:
        # Normalize to [0, 2π] range
        angles = ((pseudotime_flat - pseudotime_flat.min()) / 
                 (pseudotime_flat.max() - pseudotime_flat.min())) * 2 * np.pi
    
    # Method 1: Quantile-based assignment (ensures all cells get assigned)
    sorted_indices = np.argsort(angles)
    n_cells = len(angles)
    
    # Divide into three roughly equal groups
    g1_end = n_cells // 3
    s_end = 2 * n_cells // 3
    
    quantile_phases = [''] * n_cells
    for i, idx in enumerate(sorted_indices):
        if i < g1_end:
            quantile_phases[idx] = 'g0/g1'
        elif i < s_end:
            quantile_phases[idx] = 's'
        else:
            quantile_phases[idx] = 'g2/m'
    
    print("Step 2: Density-based boundary refinement...")
    
    # Method 2: Use kernel density estimation to find natural boundaries
    try:
        # Convert to circular coordinates for density estimation
        x_coords = np.cos(angles)
        y_coords = np.sin(angles)
        
        # Find density peaks using KDE
        kde = stats.gaussian_kde(angles.T)
        angle_range = np.linspace(0, 2*np.pi, 1000)
        density = kde(angle_range)
        
        # Find local minima (boundaries between phases)
        from scipy.signal import find_peaks
        # Find peaks in density
        peaks, _ = find_peaks(density, height=np.percentile(density, 50))
        # Find valleys (minima) between peaks
        valleys, _ = find_peaks(-density)
        
        if len(valleys) >= 2:
            # Sort valleys by position
            valley_angles = angle_range[valleys]
            valley_angles = np.sort(valley_angles)
            
            # Use the two most prominent valleys as boundaries
            if len(valley_angles) >= 2:
                boundary1 = valley_angles[0]
                boundary2 = valley_angles[1] if valley_angles[1] > valley_angles[0] else valley_angles[-1]
                
                # Ensure boundaries are properly spaced
                if boundary2 - boundary1 < np.pi:
                    boundary2 = valley_angles[-1] if len(valley_angles) > 2 else boundary1 + 2*np.pi/3
                
                # Assign phases based on density boundaries
                density_phases = []
                for angle in angles:
                    if angle < boundary1:
                        density_phases.append('g0/g1')
                    elif angle < boundary2:
                        density_phases.append('s')
                    else:
                        density_phases.append('g2/m')
                        
                print(f"Found density boundaries at {boundary1:.2f} and {boundary2:.2f}")
            else:
                density_phases = quantile_phases
                print("Using quantile-based assignment (insufficient density peaks)")
        else:
            density_phases = quantile_phases
            print("Using quantile-based assignment (no clear density boundaries)")
            
    except Exception as e:
        print(f"Density estimation failed: {e}")
        density_phases = quantile_phases
    
    print("Step 3: Consistency check and smoothing...")
    
    # Method 3: Local consistency smoothing
    final_phases = density_phases.copy()
    
    # Build nearest neighbors for smoothing
    nn = NearestNeighbors(n_neighbors=min(20, n_cells//10), metric='euclidean')
    circular_coords = np.column_stack([np.cos(angles), np.sin(angles)])
    nn.fit(circular_coords)
    
    # Smooth assignments based on local neighborhoods
    changes_made = 0
    for i in range(n_cells):
        # Find neighbors
        distances, indices = nn.kneighbors([circular_coords[i]])
        neighbor_indices = indices[0][1:]  # Exclude self
        
        # Get neighbor phases
        neighbor_phases = [final_phases[j] for j in neighbor_indices]
        
        # If current assignment disagrees with >60% of neighbors, change it
        current_phase = final_phases[i]
        phase_counts = {phase: neighbor_phases.count(phase) for phase in ['g0/g1', 's', 'g2/m']}
        most_common_phase = max(phase_counts, key=phase_counts.get)
        
        if (phase_counts[most_common_phase] > len(neighbor_phases) * 0.6 and 
            current_phase != most_common_phase):
            final_phases[i] = most_common_phase
            changes_made += 1
    
    print(f"Smoothing changed {changes_made} assignments")
    
    print("Step 4: Final validation and confidence scoring...")
    
    # Calculate confidence scores based on local consistency
    confidence_scores = np.zeros(n_cells)
    for i in range(n_cells):
        distances, indices = nn.kneighbors([circular_coords[i]])
        neighbor_indices = indices[0][1:]
        neighbor_phases = [final_phases[j] for j in neighbor_indices]
        
        # Confidence = fraction of neighbors with same phase
        same_phase_count = neighbor_phases.count(final_phases[i])
        confidence_scores[i] = same_phase_count / len(neighbor_phases)
    
    # Final validation
    phase_counts = pd.Series(final_phases).value_counts()
    total_cells = len(final_phases)
    
    print("Final phase distribution:")
    for phase in ['g0/g1', 's', 'g2/m']:
        count = phase_counts.get(phase, 0)
        percentage = (count / total_cells) * 100
        print(f"  {phase}: {count} cells ({percentage:.1f}%)")
    
    # Check for gaps (should be none with this method)
    unassigned = sum(1 for phase in final_phases if phase not in ['g0/g1', 's', 'g2/m'])
    if unassigned > 0:
        print(f"Warning: {unassigned} cells left unassigned")
    else:
        print("All cells successfully assigned!")
    
    return final_phases, confidence_scores

# Assign cell cycle stages using robust method
stages, confidence_scores = assign_cell_cycle_stage_robust(pseudotime_flat)

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

# Additional plot: confidence scores and phase boundaries
plt.figure(figsize=(15, 5))

# Plot 1: Confidence distribution
plt.subplot(1, 3, 1)
plt.hist(confidence_scores, bins=50, alpha=0.7, edgecolor='black')
plt.xlabel('Confidence Score')
plt.ylabel('Number of Cells')
plt.title('Assignment Confidence Distribution')

# Plot 2: Confidence by phase
plt.subplot(1, 3, 2)
stage_conf = pd.DataFrame({'stage': stages, 'confidence': confidence_scores})
for stage in ['g0/g1', 's', 'g2/m']:
    stage_data = stage_conf[stage_conf['stage'] == stage]['confidence']
    plt.hist(stage_data, alpha=0.7, label=stage, bins=30)
plt.xlabel('Confidence Score')
plt.ylabel('Number of Cells')
plt.title('Confidence by Cell Cycle Stage')
plt.legend()

# Plot 3: Pseudotime vs confidence
plt.subplot(1, 3, 3)
colors = [color_map['stage'][stage] for stage in stages]
plt.scatter(pseudotime_flat, confidence_scores, c=colors, alpha=0.6, s=10)
plt.xlabel('Pseudotime')
plt.ylabel('Confidence Score')
plt.title('Confidence vs Pseudotime')

plt.tight_layout()
plt.savefig(output.replace('.h5ad', '_cyclum_confidence.pdf'), dpi=300, bbox_inches='tight')
plt.close()

# Save the adata object
adata.write_h5ad(output)

print(f"\nAnalysis complete! Results saved to: {output}")
print(f"Plots saved as: {output.replace('.h5ad', '_cyclum_*.pdf')}")