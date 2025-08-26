#!/bin/bash
#SBATCH --job-name=ovary_pip_summary_statistics
#SBATCH --output=/private/groups/russelllab/jodie/scRNAseq/pipseq/ovary_pipseq/v1_by_condition/logs/out/summary_stats_%j.out
#SBATCH --error=/private/groups/russelllab/jodie/scRNAseq/pipseq/ovary_pipseq/v1_by_condition/logs/error/summary_stats_%j.err
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=short

# Print some job information
echo "Starting at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Job name: $SLURM_JOB_NAME"
echo "Node: $SLURMD_NODENAME"

# Set input and output paths (modify these as needed)
INTEGRATED_DIR="/private/groups/russelllab/jodie/scRNAseq/pipseq/ovary_pipseq/v1_by_condition/integrated_h5ad"
RESULTS_DIR="/private/groups/russelllab/jodie/scRNAseq/pipseq/ovary_pipseq/v1_by_condition/summary_statistics"
SCRIPT_PATH="/private/groups/russelllab/jodie/scRNAseq/scRNAseq_methods_comparison/analysis/summary_statistics.py"

# Create results directory if it doesn't exist
mkdir -p $RESULTS_DIR

# Find all integrated h5ad files
INTEGRATED_FILES=$(find $INTEGRATED_DIR -name "*_integrated_bbknn.h5ad")

# Process each integrated file
for INTEGRATED_FILE in $INTEGRATED_FILES; do
    # Extract the condition name from the filename
    FILENAME=$(basename "$INTEGRATED_FILE")
    CONDITION=${FILENAME%%_integrated_bbknn.h5ad}
    
    echo "Processing condition: $CONDITION"
    
    # Create output directory for this condition
    CONDITION_DIR="$RESULTS_DIR/$CONDITION"
    mkdir -p "$CONDITION_DIR"
    
    # Run the summary statistics script
    echo "Running summary statistics analysis for $CONDITION..."
    python $SCRIPT_PATH "$INTEGRATED_FILE" "$CONDITION_DIR"
    
    # Check if analysis was successful
    if [ $? -eq 0 ]; then
        echo "Successfully analyzed $CONDITION"
    else
        echo "Error analyzing $CONDITION"
        echo "Check error log for details"
    fi
    
    echo "Results for $CONDITION saved to $CONDITION_DIR"
    echo "----------------------------------------"
done

# Print completion message
echo "Job completed at: $(date)"
echo "All analyses saved to: $RESULTS_DIR"