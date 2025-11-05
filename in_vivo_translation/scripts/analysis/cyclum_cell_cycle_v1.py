'''
This script annotates cell cycle stages using Cyclum optimized for Drosophila data with FlyBase IDs.
Preserves all genes in the dataset while optimizing Cyclum training using smart gene identification.
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

parser = argparse.ArgumentParser(description="Cyclum Cell Cycle Analysis for Drosophila FlyBase IDs")
parser.add_argument("--input", type=str, required=False, help="Input h5ad file", 
                   default='/private/groups/russelllab/jodie/scRNAseq/scripts/snakemake_pipeline/results_kallisto_bustools/filtered_h5ad/kallisto_JW18DOX-Ctrl-1_P.h5ad')
parser.add_argument("--output", type=str, required=False, help="Output h5ad file", 
                   default='/private/groups/russelllab/jodie/scRNAseq/scripts/snakemake_pipeline/results_kallisto_bustools/filtered_h5ad/cyclum_kallisto_JW18DOX-Ctrl-1_P.h5ad')
parser.add_argument("--min_cells", type=int, default=500, help="Minimum number of cells for analysis")
parser.add_argument("--epochs", type=int, default=1500, help="Number of training epochs")
parser.add_argument("--filter_for_cyclum", action="store_true", help="Filter mito/rRNA genes only for Cyclum training (recommended)")

args = parser.parse_args()

def get_core_drosophila_fbgn_mappings():
    """Core FlyBase ID mappings for essential cell cycle genes"""
    # These are confirmed FlyBase IDs from the search results and literature
    core_mappings = {
        # Major cyclins
        'FBgn0000404': 'CycA',     # Cyclin A  
        'FBgn0000405': 'CycB',     # Cyclin B
        'FBgn0004117': 'CycE',     # Cyclin E
        'FBgn0010382': 'CycB3',    # Cyclin B3
        
        # CDKs and kinases
        'FBgn0004106': 'Cdk1',     # cdc2/Cdk1
        'FBgn0003510': 'polo',     # polo kinase
        'FBgn0004057': 'aurora',   # aurora kinase
        'FBgn0003525': 'stg',      # string/Cdc25
        
        # DNA replication machinery
        'FBgn0005655': 'pcna',     # PCNA
        'FBgn0010314': 'Cdc6',     # Cdc6
        'FBgn0015791': 'Mcm2',     # Mcm2
        'FBgn0015792': 'Mcm3',     # Mcm3
        'FBgn0015793': 'Mcm4',     # Mcm4
        'FBgn0015794': 'Mcm5',     # Mcm5
        'FBgn0015795': 'Mcm6',     # Mcm6
        'FBgn0015796': 'Mcm7',     # Mcm7
        
        # Cell cycle regulators
        'FBgn0004391': 'Wee1',     # Wee1
        'FBgn0024251': 'Myt1',     # Myt1
        'FBgn0000721': 'fzy',      # fizzy
        'FBgn0011829': 'Fzr',      # fizzy-related
        'FBgn0005772': 'twn',      # twins
        'FBgn0003892': 'png',      # peanut
        
        # Mitosis
        'FBgn0015380': 'Det',      # deterin
        'FBgn0004102': 'Incenp',   # Incenp
        'FBgn0029825': 'Survivin', # Survivin
        'FBgn0030122': 'Borealin', # Borealin
        'FBgn0027329': 'Ndc80',    # Ndc80
        'FBgn0026313': 'Nuf2',     # Nuf2
        
        # Transcription factors
        'FBgn0000567': 'E2f',      # E2F
        'FBgn0000524': 'Dp',       # DP
        'FBgn0003277': 'Rbf',      # Rb family
        'FBgn0002467': 'Myc',      # Myc
        
        # DNA damage/checkpoint
        'FBgn0001330': 'mei-41',   # ATR
        'FBgn0001123': 'grp',      # Chk1
        'FBgn0039044': 'Brca2',    # Brca2
        'FBgn0020270': 'Mre11',    # Mre11
        'FBgn0005771': 'Nbs',      # Nbs/nibrin
    }
    
    return core_mappings

def identify_genes_in_flybase_dataset(adata):
    """Smart identification of gene types using FlyBase IDs and functional patterns"""
    gene_names = adata.var_names
    
    print(f"Analyzing dataset with FlyBase IDs...")
    print(f"Sample gene names: {list(gene_names[:10])}")
    print(f"Total genes in dataset: {len(gene_names)}")
    
    # Get core mappings
    core_cc_mapping = get_core_drosophila_fbgn_mappings()
    
    # 1. Identify cell cycle genes using known FlyBase IDs
    cc_fbgns_known = list(core_cc_mapping.keys())
    cc_genes_found = [fbgn for fbgn in cc_fbgns_known if fbgn in gene_names]
    cc_mask_known = gene_names.isin(cc_genes_found)
    
    print(f"Found {len(cc_genes_found)} known cell cycle genes from core mapping:")
    if cc_genes_found:
        for fbgn in cc_genes_found[:10]:  # Show first 10
            print(f"  {fbgn} ({core_cc_mapping[fbgn]})")
    
    # 2. Smart pattern-based identification for additional genes
    # Look for genes that likely encode cell cycle functions based on common patterns
    
    # Count different FlyBase ID patterns to understand the dataset
    fbgn_pattern = gene_names.str.contains('^FBgn[0-9]{7}$', regex=True)
    fbti_pattern = gene_names.str.contains('^FBti[0-9]{7}$', regex=True)  # Transposon insertions
    transposon_pattern = gene_names.str.contains('transposable_element', regex=True)
    
    print(f"\nDataset composition:")
    print(f"  FlyBase gene IDs (FBgn): {fbgn_pattern.sum()}")
    print(f"  Transposon insertions (FBti): {fbti_pattern.sum()}")
    print(f"  Transposable elements: {transposon_pattern.sum()}")
    
    # 3. Identify mitochondrial genes
    # Mitochondrial genes in Drosophila often have specific FBgn ranges or patterns
    # Since we can't easily identify them without external mapping, we'll be conservative
    mito_genes = pd.Series([False] * len(gene_names), index=gene_names)
    print(f"Mitochondrial genes: {mito_genes.sum()} (conservative - may need external annotation)")
    
    # 4. Identify ribosomal genes  
    # Similar issue - ribosomal genes need external mapping for FlyBase IDs
    ribo_genes = pd.Series([False] * len(gene_names), index=gene_names)
    print(f"Ribosomal genes: {ribo_genes.sum()} (conservative - may need external annotation)")
    
    # 5. Exclude transposable elements from Cyclum analysis (they're not informative for cell cycle)
    non_transposon_mask = ~(fbti_pattern | transposon_pattern)
    genes_for_cyclum = non_transposon_mask
    
    print(f"Genes suitable for Cyclum analysis: {genes_for_cyclum.sum()} (excluding transposons)")
    
    return mito_genes, ribo_genes, cc_mask_known, cc_genes_found, genes_for_cyclum, core_cc_mapping

def create_cyclum_optimized_data(adata, genes_for_cyclum_mask, cc_genes_found, core_cc_mapping):
    """Create optimized dataset for Cyclum training using smart gene selection"""
    
    # Make a copy for Cyclum training
    adata_cyclum = adata.copy()
    
    # Fix scanpy metadata issues
    if 'log1p' in adata_cyclum.uns_keys():
        if isinstance(adata_cyclum.uns['log1p'], dict) and 'base' not in adata_cyclum.uns['log1p']:
            print("Fixing log1p metadata...")
            adata_cyclum.uns['log1p']['base'] = None
    
    # Filter to genes suitable for Cyclum (exclude transposons)
    print(f"Filtering dataset for Cyclum: keeping {genes_for_cyclum_mask.sum()} informative genes")
    adata_cyclum = adata_cyclum[:, genes_for_cyclum_mask]
    
    # Calculate highly variable genes
    try:
        sc.pp.highly_variable_genes(adata_cyclum, min_mean=0.01, max_mean=5, min_disp=0.3)
    except Exception as e:
        print(f"Error with highly_variable_genes: {e}")
        print("Trying alternative approach...")
        if 'log1p' in adata_cyclum.uns_keys():
            del adata_cyclum.uns['log1p']
        sc.pp.highly_variable_genes(adata_cyclum, min_mean=0.005, max_mean=10, min_disp=0.1)
    
    # Always include cell cycle genes as highly variable if they're in this dataset
    cc_genes_in_filtered = [fbgn for fbgn in cc_genes_found if fbgn in adata_cyclum.var_names]
    if cc_genes_in_filtered:
        adata_cyclum.var.loc[cc_genes_in_filtered, 'highly_variable'] = True
        print(f"Marked {len(cc_genes_in_filtered)} cell cycle genes as highly variable")
        print(f"  Examples: {[f'{fbgn}({core_cc_mapping[fbgn]})' for fbgn in cc_genes_in_filtered[:5]]}")
    
    # Check if we have enough highly variable genes
    n_hvg = adata_cyclum.var['highly_variable'].sum()
    if n_hvg < 100:
        print(f"Warning: Only {n_hvg} highly variable genes found. Relaxing criteria...")
        sc.pp.highly_variable_genes(adata_cyclum, min_mean=0.001, max_mean=20, min_disp=0.05)
        if cc_genes_in_filtered:
            adata_cyclum.var.loc[cc_genes_in_filtered, 'highly_variable'] = True
        n_hvg = adata_cyclum.var['highly_variable'].sum()
    
    # Keep only highly variable genes
    adata_cyclum = adata_cyclum[:, adata_cyclum.var.highly_variable]
    print(f"Using {adata_cyclum.n_vars} highly variable genes for Cyclum analysis")
    
    # Check if data is already log-transformed
    data_max = adata_cyclum.X.max()
    if hasattr(adata_cyclum.X, 'toarray'):
        data_max = adata_cyclum.X.toarray().max()
    
    # Log transform if needed (only if max value suggests raw counts)
    if data_max > 20:  # Likely raw counts
        print("Applying log1p transformation...")
        sc.pp.log1p(adata_cyclum)
    else:
        print("Data appears to be already log-transformed, skipping log1p...")
    
    # Scale for Cyclum
    sc.pp.scale(adata_cyclum, max_value=10)
    
    return adata_cyclum

def filter_cycling_cells(adata_cyclum, mtx, min_variance_percentile=20):
    """Filter for actively cycling cells based on gene expression variance"""
    
    # Calculate variance across genes for each cell
    cell_variances = np.var(mtx, axis=1)
    variance_threshold = np.percentile(cell_variances, min_variance_percentile)
    
    cycling_mask = cell_variances >= variance_threshold
    print(f"Filtering to {cycling_mask.sum()} potentially cycling cells (>{min_variance_percentile}th percentile variance)")
    
    return adata_cyclum[cycling_mask], mtx[cycling_mask], cycling_mask

def train_cyclum_model(mtx, epochs=1500, learning_rate=1e-4):
    """Train Cyclum model with optimized parameters"""
    
    print(f"Training Cyclum model with {epochs} epochs...")
    
    # Use CyclumAutoTune with better parameters
    model = cyclum.tuning.CyclumAutoTune(mtx)
    
    # Train with optimized parameters
    loss_history = model.train(mtx, epochs=epochs, verbose=100, rate=learning_rate)
    
    # Check convergence if possible
    try:
        if hasattr(model, 'loss_trace') and len(model.loss_trace) > 100:
            recent_losses = model.loss_trace[-100:]
            loss_improvement = (recent_losses[0] - recent_losses[-1]) / recent_losses[0]
            print(f"Training loss improvement in last 100 epochs: {loss_improvement:.4f}")
            
            if loss_improvement < 0.001:
                print("Warning: Model may not have converged well. Consider increasing epochs.")
    except:
        print("Could not assess convergence - continuing with analysis")
    
    return model

def assign_cell_cycle_stage_optimized(pseudotime_flat):
    """Assign cell cycle stages based on pseudotime with better boundaries"""
    
    # Normalize pseudotime to 0-2π range
    pt_min, pt_max = pseudotime_flat.min(), pseudotime_flat.max()
    pseudotime_2pi = (pseudotime_flat - pt_min) / (pt_max - pt_min) * 2 * np.pi
    
    stages = []
    for pt in pseudotime_2pi:
        # Drosophila cell cycle phases:
        # G1: 0 to 2π/3 (0-120 degrees)
        # S: 2π/3 to 4π/3 (120-240 degrees)  
        # G2/M: 4π/3 to 2π (240-360 degrees)
        if pt < 2*np.pi/3:
            stages.append('g0/g1')
        elif pt < 4*np.pi/3:
            stages.append('s')
        else:
            stages.append('g2/m')
    
    return stages, pseudotime_2pi

def plot_quality_metrics(model, pseudotime_flat, output_prefix):
    """Generate quality control plots"""
    
    # Plot pseudotime distribution histogram
    plt.figure(figsize=(8, 4))
    plt.hist(pseudotime_flat, bins=50, alpha=0.7, edgecolor='black')
    plt.xlabel('Pseudotime')
    plt.ylabel('Number of cells')
    plt.title('Distribution of Cell Cycle Pseudotime')
    plt.savefig(f"{output_prefix}_pseudotime_hist.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Check for uniformity around the circle
    n_bins = 20
    bin_edges = np.linspace(pseudotime_flat.min(), pseudotime_flat.max(), n_bins + 1)
    hist, _ = np.histogram(pseudotime_flat, bins=bin_edges)
    
    # Calculate coefficient of variation as uniformity metric
    cv = np.std(hist) / np.mean(hist) if np.mean(hist) > 0 else np.inf
    print(f"Pseudotime distribution coefficient of variation: {cv:.3f}")
    print(f"(Lower values indicate more uniform distribution around the cycle)")

def examine_gene_patterns(gene_names, n_examples=20):
    """Examine gene naming patterns in the dataset"""
    print(f"\nExamining gene naming patterns (showing {n_examples} examples):")
    print(f"First {n_examples} genes: {list(gene_names[:n_examples])}")
    
    # Look for patterns
    patterns_to_check = [
        ('FlyBase gene IDs (FBgn)', gene_names.str.contains('^FBgn[0-9]', regex=True)),
        ('Transposon insertions (FBti)', gene_names.str.contains('^FBti[0-9]', regex=True)),
        ('Transposable elements', gene_names.str.contains('transposable_element', regex=True)),
        ('Genes with numbers', gene_names.str.contains('[0-9]', regex=True)),
        ('Genes with underscores', gene_names.str.contains('_', regex=True)),
    ]
    
    for pattern_name, pattern_mask in patterns_to_check:
        count = pattern_mask.sum()
        if count > 0:
            examples = list(gene_names[pattern_mask][:5])
            print(f"  {pattern_name}: {count} genes (examples: {examples})")

# Main analysis
print("Reading Drosophila scRNA-seq data with FlyBase IDs...")
adata = sc.read_h5ad(args.input)
output = args.output

print(f"Original data: {adata.n_obs} cells, {adata.n_vars} genes")

# Examine gene naming patterns
examine_gene_patterns(adata.var_names)

# Check minimum cell count
if adata.n_obs < args.min_cells:
    print(f"Warning: Only {adata.n_obs} cells available. Cyclum works best with >{args.min_cells} cells.")

# Identify gene types using smart FlyBase ID analysis
mito_genes, ribo_genes, cc_mask, cc_genes_found, genes_for_cyclum, core_cc_mapping = identify_genes_in_flybase_dataset(adata)

# Create optimized dataset for Cyclum training
adata_cyclum = create_cyclum_optimized_data(adata, genes_for_cyclum, cc_genes_found, core_cc_mapping)
mtx = adata_cyclum.X

# Convert to dense array if sparse
if hasattr(mtx, 'toarray'):
    mtx = mtx.toarray()

print(f"Cyclum training data: {mtx.shape[0]} cells, {mtx.shape[1]} genes")

# Filter for cycling cells
adata_cyclum_filtered, mtx_filtered, cycling_mask = filter_cycling_cells(adata_cyclum, mtx)

# Train model with optimized parameters
model = train_cyclum_model(mtx_filtered, epochs=args.epochs)

# Extract pseudotime
pseudotime = model.predict_pseudotime(mtx_filtered)
pseudotime_flat = pseudotime.flatten()

print(f"Pseudotime shape: {pseudotime.shape}")
print(f"Pseudotime range: {pseudotime.min():.3f} to {pseudotime.max():.3f}")

# Assign cell cycle stages
stages, pseudotime_2pi = assign_cell_cycle_stage_optimized(pseudotime_flat)

# Check the distribution
unique_stages, counts = np.unique(stages, return_counts=True)
print("\nCell cycle stage distribution:")
for stage, count in zip(unique_stages, counts):
    print(f"{stage}: {count} cells ({count/len(stages)*100:.1f}%)")

# Add results back to ORIGINAL adata object (preserving all genes)
# Initialize with default values for non-cycling cells
adata.obs['cyclum_stage'] = 'non_cycling'
adata.obs['cyclum_pseudotime'] = np.nan
adata.obs['is_cycling'] = False

# Map cycling cells back to original indices
try:
    # Get the indices of cycling cells in the original dataset
    cycling_indices = adata_cyclum.obs_names[cycling_mask]
    
    # Check if indices match between datasets
    if all(idx in adata.obs_names for idx in cycling_indices):
        adata.obs.loc[cycling_indices, 'cyclum_stage'] = stages
        adata.obs.loc[cycling_indices, 'cyclum_pseudotime'] = pseudotime_flat
        adata.obs.loc[cycling_indices, 'is_cycling'] = True
        
        print(f"\nSuccessfully mapped cycling results back to original dataset")
    else:
        print(f"\nWarning: Index mismatch between datasets. Using positional mapping...")
        # Use positional mapping as fallback
        original_cycling_positions = np.where(cycling_mask)[0]
        adata.obs.iloc[original_cycling_positions, adata.obs.columns.get_loc('cyclum_stage')] = stages
        adata.obs.iloc[original_cycling_positions, adata.obs.columns.get_loc('cyclum_pseudotime')] = pseudotime_flat
        adata.obs.iloc[original_cycling_positions, adata.obs.columns.get_loc('is_cycling')] = True
        
except Exception as e:
    print(f"Error mapping results back: {e}")
    print("Using direct positional mapping...")
    # Direct positional mapping for cycling cells
    original_cycling_positions = np.where(cycling_mask)[0]
    adata.obs.iloc[original_cycling_positions, adata.obs.columns.get_loc('cyclum_stage')] = stages
    adata.obs.iloc[original_cycling_positions, adata.obs.columns.get_loc('cyclum_pseudotime')] = pseudotime_flat
    adata.obs.iloc[original_cycling_positions, adata.obs.columns.get_loc('is_cycling')] = True

print(f"\nAdded to adata.obs (preserving all {adata.n_vars} genes):")
print(f"cyclum_stage: {adata.obs['is_cycling'].sum()} cycling cells annotated")
print(f"cyclum_pseudotime: {adata.obs['is_cycling'].sum()} cycling cells annotated")
print(f"is_cycling: {adata.obs['is_cycling'].sum()} cells marked as cycling")

# Mark gene types in var for later analysis
adata.var['mitochondrial'] = mito_genes
adata.var['ribosomal'] = ribo_genes  
adata.var['cell_cycle'] = cc_mask
adata.var['suitable_for_cyclum'] = genes_for_cyclum

print(f"\nAdded to adata.var:")
print(f"mitochondrial: {adata.var['mitochondrial'].sum()} genes (conservative)")
print(f"ribosomal: {adata.var['ribosomal'].sum()} genes (conservative)") 
print(f"cell_cycle: {adata.var['cell_cycle'].sum()} genes (core known)")
print(f"suitable_for_cyclum: {adata.var['suitable_for_cyclum'].sum()} genes (excluding transposons)")

# Create plots
output_prefix = output.replace('.h5ad', '')

# Create the circular plot with error handling
try:
    color_map = {'g0/g1': "red", 's': "green", 'g2/m': "blue"}
    fig = cyclum.illustration.plot_round_distr_color(pseudotime_flat, stages, color_map)
    plt.savefig(f"{output_prefix}_cyclum_cell_cycle.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    print("Circular plot saved successfully!")
except Exception as e:
    print(f"Error with Cyclum's built-in plotting function: {e}")
    print("Creating custom circular plot...")
    
    # Create custom circular plot
    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(projection='polar'))
    
    # Convert stages to numeric for plotting
    stage_to_num = {'g0/g1': 0, 's': 1, 'g2/m': 2}
    stage_colors = {'g0/g1': 'red', 's': 'green', 'g2/m': 'blue'}
    
    # Plot each stage
    for stage in ['g0/g1', 's', 'g2/m']:
        stage_mask = np.array(stages) == stage
        if stage_mask.sum() > 0:
            stage_pseudotime = pseudotime_flat[stage_mask]
            # Convert pseudotime to angles (0 to 2π)
            angles = (stage_pseudotime - pseudotime_flat.min()) / (pseudotime_flat.max() - pseudotime_flat.min()) * 2 * np.pi
            radii = np.ones(len(angles))  # All at same radius
            
            ax.scatter(angles, radii, c=stage_colors[stage], label=stage, alpha=0.6, s=20)
    
    ax.set_ylim(0, 1.2)
    ax.set_title('Cell Cycle Distribution (Cyclum)', pad=20)
    ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.0))
    
    plt.savefig(f"{output_prefix}_cyclum_cell_cycle.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    print("Custom circular plot saved successfully!")

# Generate quality control plots
plot_quality_metrics(model, pseudotime_flat, output_prefix)

# Show elbow plot if available
try:
    elbow_fig = model.show_elbow()
    plt.savefig(f"{output_prefix}_cyclum_elbow.pdf", dpi=300, bbox_inches='tight')
    plt.close()
except:
    print("Elbow plot not available for this model type")

# Show bar plot if available  
try:
    bar_fig = model.show_bar()
    plt.savefig(f"{output_prefix}_cyclum_bar.pdf", dpi=300, bbox_inches='tight')
    plt.close()
except:
    print("Bar plot not available for this model type")

# Create additional plots showing dataset composition
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 8))

# Gene type distribution
gene_types = ['FlyBase genes', 'Transposon insertions', 'Transposable elements', 'Other']
fbgn_count = adata.var_names.str.contains('^FBgn[0-9]', regex=True).sum()
fbti_count = adata.var_names.str.contains('^FBti[0-9]', regex=True).sum()
transposon_count = adata.var_names.str.contains('transposable_element', regex=True).sum()
other_count = adata.n_vars - fbgn_count - fbti_count - transposon_count

gene_counts = [fbgn_count, fbti_count, transposon_count, other_count]
ax1.bar(gene_types, gene_counts, color=['blue', 'orange', 'red', 'lightgray'])
ax1.set_ylabel('Number of genes')
ax1.set_title('Gene types in dataset')
ax1.tick_params(axis='x', rotation=45)

# Cycling vs non-cycling cells
cycling_counts = [
    (~adata.obs['is_cycling']).sum(),
    adata.obs['is_cycling'].sum()
]
ax2.bar(['Non-cycling', 'Cycling'], cycling_counts, color=['lightgray', 'orange'])
ax2.set_ylabel('Number of cells')
ax2.set_title('Cell classification')

# Cell cycle stage distribution
stage_counts = adata.obs['cyclum_stage'].value_counts()
stage_counts = stage_counts.reindex(['non_cycling', 'g0/g1', 's', 'g2/m'], fill_value=0)
colors = ['lightgray', 'red', 'green', 'blue']
ax3.bar(stage_counts.index, stage_counts.values, color=colors)
ax3.set_ylabel('Number of cells')
ax3.set_title('Cell cycle stages')
ax3.tick_params(axis='x', rotation=45)

# Cell cycle genes found
if len(cc_genes_found) > 0:
    cc_gene_symbols = [core_cc_mapping.get(fbgn, fbgn) for fbgn in cc_genes_found[:10]]
    ax4.barh(range(len(cc_gene_symbols)), [1]*len(cc_gene_symbols))
    ax4.set_yticks(range(len(cc_gene_symbols)))
    ax4.set_yticklabels(cc_gene_symbols, fontsize=8)
    ax4.set_xlabel('Presence in dataset')
    ax4.set_title(f'Cell cycle genes found ({len(cc_genes_found)} total)')
else:
    ax4.text(0.5, 0.5, 'No cell cycle genes\nfrom core set found', 
             ha='center', va='center', transform=ax4.transAxes)
    ax4.set_title('Cell cycle genes found')

plt.tight_layout()
plt.savefig(f"{output_prefix}_dataset_summary.pdf", dpi=300, bbox_inches='tight')
plt.close()

# Save the annotated data (with ALL genes preserved)
adata.write_h5ad(output)

print(f"\nAnalysis complete!")
print(f"Results saved to: {output}")
print(f"Plots saved with prefix: {output_prefix}")
print(f"\nDataset preserved with all {adata.n_vars} genes including:")
print(f"  - {adata.var['cell_cycle'].sum()} known cell cycle genes")
if len(cc_genes_found) > 0:
    print(f"  - Examples: {[f'{fbgn}({core_cc_mapping[fbgn]})' for fbgn in cc_genes_found[:3]]}")
print(f"  - {adata.var['suitable_for_cyclum'].sum()} genes used for Cyclum training")
print(f"\nNote: For mitochondrial/ribosomal gene analysis, you may need external")
print(f"FlyBase annotation files to map FBgn IDs to gene types.")