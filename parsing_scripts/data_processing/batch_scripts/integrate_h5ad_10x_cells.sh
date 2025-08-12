#!/bin/bash
#SBATCH --job-name=10x_integrate_h5ad
#SBATCH --output=/private/groups/russelllab/jodie/scRNAseq/10x/cellranger_results_v4/v1_by_conditon/logs/out/10x_integrate_h5ad_%j.out
#SBATCH --error=/private/groups/russelllab/jodie/scRNAseq/10x/cellranger_results_v4/v1_by_conditon/logs/error/10x_integrate_h5ad_%j.err
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=64G
#SBATCH --partition=short
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jomojaco@ucsc.edu

# Print some job information
echo "Starting at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Job name: $SLURM_JOB_NAME"
echo "Node: $SLURMD_NODENAME"

# Set input and output paths (modify these as needed)
INPUT_DIR="/private/groups/russelllab/jodie/scRNAseq/10x/cellranger_results_v3/h5ad_output"
OUTPUT_DIR="/private/groups/russelllab/jodie/scRNAseq/10x/cellranger_results_v4/v1_by_conditon/integrated_h5ad"
PYTHON_SCRIPT="/private/groups/russelllab/jodie/scRNAseq/scRNAseq_methods_comparison/data/scripts/integrate_h5ad_by_condition_v2.py"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Run the integration script to process files by sample type
python $PYTHON_SCRIPT \
  --input_dir "$INPUT_DIR" \
  --output_dir "$OUTPUT_DIR" \
  --sample_type_pattern "^([^-]+)-([^-]+)" \
  --batch_key "batch" \
  --min_cells 3 \
  --min_genes 200 \
  --n_pcs 50 \
  --n_neighbors 20 \
  --method "both" \
  --calculate_titer

# Print completion message
echo "Job completed at: $(date)"
echo "Integrated files saved to: $OUTPUT_DIR"