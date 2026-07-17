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
from sklearn.neighbors import NearestNeighbors
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

def assign_cell_cycle_stage_simple(pseudotime_flat):
    """
    Simple, gap-free cell cycle assignment
    """
    print("Step 1: Converting to circular coordinates...")
    
    # Normalize pseudotime to [0, 2π] 
    if pseudotime_flat.max() <= 1:
        angles = pseudotime_flat * 2 * np.pi
    else:
        # Normalize to [0, 2π] range
        angles = ((pseudotime_flat - pseudotime_flat.min()) / 
                 (pseudotime_flat.max() - pseudotime_flat.min())) * 2 * np.pi
    
    print("Step 2: Direct angle-based assignment...")
    
    # Simple direct assignment based on angle ranges
    # This guarantees complete coverage with no gaps
    phases = []
    for angle in angles:
        # Normalize angle to [0, 2π] range
        normalized_angle = angle % (2 * np.pi)
        
        if normalized_angle < (2 * np.pi / 3):  # 0 to 2π/3
            phases.append('g0/g1')
        elif normalized_angle < (4 * np.pi / 3):  # 2π/3 to 4π/3
            phases.append('s')
        else:  # 4π/3 to 2π
            phases.append('g2/m')
    
    print("Step 3: Light smoothing at boundaries...")
    
    # Apply light smoothing only near the boundaries to avoid hard transitions
    boundary1 = 2 * np.pi / 3
    boundary2 = 4 * np.pi / 3
    boundary_width = np.pi / 12  # Small smoothing window
    
    # Build simple nearest neighbors for boundary smoothing
    n_cells = len(angles)
    if n_cells > 10:  # Only smooth if we have enough cells
        nn = NearestNeighbors(n_neighbors=min(10, n_cells//10))
        circular_coords = np.column_stack([np.cos(angles), np.sin(angles)])
        nn.fit(circular_coords)
        
        smoothed_phases = phases.copy()
        changes_made = 0
        
        for i, angle in enumerate(angles):
            normalized_angle = angle % (2 * np.pi)
            
            # Check if near boundaries
            near_boundary = (abs(normalized_angle - boundary1) < boundary_width or 
                           abs(normalized_angle - boundary2) < boundary_width or
                           abs(normalized_angle - 0) < boundary_width or
                           abs(normalized_angle - 2*np.pi) < boundary_width)
            
            if near_boundary:
                # Get neighbors and their phases
                distances, indices = nn.kneighbors([circular_coords[i]])
                neighbor_indices = indices[0][1:]  # Exclude self
                neighbor_phases = [phases[j] for j in neighbor_indices]
                
                # If most neighbors agree on a different phase, consider changing
                current_phase = phases[i]
                phase_counts = {}
                for phase in neighbor_phases:
                    phase_counts[phase] = phase_counts.get(phase, 0) + 1
                
                if phase_counts:
                    most_common = max(phase_counts, key=phase_counts.get)
                    # Only change if >70% of neighbors agree and it's different
                    if (phase_counts[most_common] > len(neighbor_phases) * 0.7 and 
                        most_common != current_phase):
                        smoothed_phases[i] = most_common
                        changes_made += 1
        
        phases = smoothed_phases
        print(f"Boundary smoothing changed {changes_made} assignments")
    
    print("Step 4: Calculate confidence scores...")
    
    # Simple confidence based on distance from phase boundaries
    confidence_scores = np.ones(len(angles))  # Start with high confidence
    
    for i, angle in enumerate(angles):
        normalized_angle = angle % (2 * np.pi)
        
        # Distance from nearest boundary
        dist_to_b1 = min(abs(normalized_angle - boundary1), 2*np.pi - abs(normalized_angle - boundary1))
        dist_to_b2 = min(abs(normalized_angle - boundary2), 2*np.pi - abs(normalized_angle - boundary2))
        dist_to_start = min(normalized_angle, 2*np.pi - normalized_angle)
        
        min_dist_to_boundary = min(dist_to_b1, dist_to_b2, dist_to_start)
        
        # Confidence decreases near boundaries
        max_dist = np.pi / 3  # Maximum distance from boundary in a phase
        confidence_scores[i] = min(1.0, min_dist_to_boundary / (boundary_width * 2))
    
    # Final validation
    phase_counts = pd.Series(phases).value_counts()
    total_cells = len(phases)
    
    print("Final phase distribution:")
    for phase in ['g0/g1', 's', 'g2/m']:
        count = phase_counts.get(phase, 0)
        percentage = (count / total_cells) * 100
        print(f"  {phase}: {count} cells ({percentage:.1f}%)")
    
    # Verify no gaps
    unassigned = sum(1 for phase in phases if phase not in ['g0/g1', 's', 'g2/m'])
    print(f"Unassigned cells: {unassigned} (should be 0)")
    
    return phases, confidence_scores

# Assign cell cycle stages using simple method
stages, confidence_scores = assign_cell_cycle_stage_simple(pseudotime_flat)

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
plt.figure(figsize=(12, 4))

plt.subplot(1, 3, 1)
plt.hist(confidence_scores, bins=30, alpha=0.7, edgecolor='black')
plt.xlabel('Confidence Score')
plt.ylabel('Number of Cells')
plt.title('Assignment Confidence')

plt.subplot(1, 3, 2)
colors = [color_map['stage'][stage] for stage in stages]
plt.scatter(pseudotime_flat, confidence_scores, c=colors, alpha=0.6, s=10)
plt.xlabel('Pseudotime')
plt.ylabel('Confidence Score')
plt.title('Confidence vs Pseudotime')

plt.subplot(1, 3, 3)
# Show phase distribution around the circle
angles = pseudotime_flat * 2 * np.pi if pseudotime_flat.max() <= 1 else ((pseudotime_flat - pseudotime_flat.min()) / (pseudotime_flat.max() - pseudotime_flat.min())) * 2 * np.pi
plt.hist(angles, bins=60, alpha=0.7, edgecolor='black')
plt.xlabel('Angle (radians)')
plt.ylabel('Number of Cells')
plt.title('Cell Distribution Around Circle')
plt.axvline(2*np.pi/3, color='red', linestyle='--', alpha=0.7, label='G1/S boundary')
plt.axvline(4*np.pi/3, color='green', linestyle='--', alpha=0.7, label='S/G2M boundary')
plt.legend()

plt.tight_layout()
plt.savefig(output.replace('.h5ad', '_cyclum_confidence.pdf'), dpi=300, bbox_inches='tight')
plt.close()

# Save the adata object
adata.write_h5ad(output)

print(f"\nAnalysis complete! Results saved to: {output}")
print(f"Plots saved as: {output.replace('.h5ad', '_cyclum_*.pdf')}")