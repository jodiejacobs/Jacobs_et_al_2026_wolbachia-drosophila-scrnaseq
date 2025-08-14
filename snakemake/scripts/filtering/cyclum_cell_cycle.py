import cyclum 
from cyclum import cyclum.models
from cyclum import cyclum.tuning 
import scanpy as sc 
import argparse
import os
import sklearn


adata = sc.read_h5ad(file_path)

# Train model 
model = cyclum.tuning.CyclumAutoTune(adata.raw.X.toarray())
pseudotime = model.predict_pseudotime(adata.raw.X.toarray())

# Extract Pseudotime
pseudotime = model.predict_pseudotime(adata.raw.X.toarray())

# Train with more epochs for better cell detection 
history = model.train(adata.raw.X.toarray(), epochs=800, verbose=100, rate=2e-4)

# Extract the circular pseudotime (this represents cell cycle phase)
pseudotime = model.predict_pseudotime(adata.raw.X.toarray())

# Flatten pseudotime and create labels
pseudotime_flat = pseudotime.flatten()

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

label = {'stage': np.array(stages)}

unique_stages, counts = np.unique(stages, return_counts=True)

adata.obs['cyclum_stage'] = stages
adata.obs['cyclum_pseudotime'] = pseudotime_flat

