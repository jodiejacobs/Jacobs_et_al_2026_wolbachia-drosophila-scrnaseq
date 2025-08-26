#!/bin/bash
# This script will run Cell Ranger count on multiple samples in parallel with improved error handling

# Set up paths
DATA_DIR="/private/groups/russelllab/jodie/sequencing_data/scRNAseq/10x_cells/Jacobs_10525_241105A9"
OUTPUT_DIR="/private/groups/russelllab/jodie/scRNAseq/10x/cellranger_results_v3"
REF_DIR="/private/groups/russelllab/jodie/reference_genomes/10X_reference_genome/Drosophila_melanogaster_wMel_ASM1658442v1"
SCRIPTS_DIR="$OUTPUT_DIR/scripts"

# Create directories
mkdir -p $OUTPUT_DIR
mkdir -p $SCRIPTS_DIR
mkdir -p $OUTPUT_DIR/logs

# Go to data directory
cd $DATA_DIR
echo "Working directory: $(pwd)"
echo "Hostname: $(hostname)"
echo "Date: $(date)"

# Find all sample prefixes WITHOUT S numbers
SAMPLES=$(ls *_L006_R1_001.fastq.gz | sed 's/_S[0-9]*_L006_R1_001.fastq.gz//' | sort | uniq)

echo "Detected samples (without S numbers):"
echo "$SAMPLES"

# Process each sample
for SAMPLE in $SAMPLES; do
    echo "Submitting job for sample: $SAMPLE"
    
    # Get the FULL filename pattern (with S number) for Cell Ranger
    FASTQ_PATTERN=$(ls ${SAMPLE}_S*_L006_R1_001.fastq.gz | head -n 1 | sed 's/_L006_R1_001.fastq.gz//')
    S_NUMBER=$(echo $FASTQ_PATTERN | grep -o "_S[0-9]*" | sed 's/_S//')
    
    echo "FASTQ pattern for Cell Ranger: $FASTQ_PATTERN (S number: $S_NUMBER)"
    
    # Create a submission script for this sample
    TEMP_SCRIPT="$SCRIPTS_DIR/${SAMPLE}_cellranger.sh"
    
    cat > $TEMP_SCRIPT << EOF
#!/bin/bash
#SBATCH --job-name=${SAMPLE}_cr
#SBATCH --partition=long
#SBATCH --mail-user=jomojaco@ucsc.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=$OUTPUT_DIR/logs/${SAMPLE}_cellranger_%j.out
#SBATCH --error=$OUTPUT_DIR/logs/${SAMPLE}_cellranger_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128gb
#SBATCH --time=72:00:00

set -e  # Exit on any error
set -o pipefail  # Exit if any command in a pipe fails

pwd; hostname; date

# Create output directory - use SLURM_JOB_ID to make unique directory
WORK_DIR="${OUTPUT_DIR}/\${SLURM_JOB_ID}_${SAMPLE}"
mkdir -p \$WORK_DIR
cd \$WORK_DIR

echo "Working in: \$(pwd)"

# Monitor disk space and memory usage throughout the job
(while true; do
    echo "\$(date): Disk usage: \$(df -h . | tail -n 1)"
    echo "\$(date): Memory usage: \$(free -h)"
    sleep 300  # Check every 5 minutes
done) > resource_monitor.log 2>&1 &
MONITOR_PID=\$!

# Ensure monitor process is killed when script exits
trap "kill \$MONITOR_PID 2>/dev/null || true" EXIT

# Run Cell Ranger count with explicit --sample flag
echo "Starting Cell Ranger for sample: ${SAMPLE} (S number: $S_NUMBER)"
cellranger count \\
    --id=${SAMPLE} \\
    --transcriptome=$REF_DIR \\
    --fastqs=$DATA_DIR \\
    --sample=${SAMPLE} \\
    --localcores=24 \\
    --localmem=96 \\
    --disable-ui \\
    --nosecondary \\
    --create-bam=true

# Verify the run completed successfully
if [ ! -f "${SAMPLE}/outs/filtered_feature_bc_matrix.h5" ]; then
    echo "ERROR: Cell Ranger failed to produce output matrix for ${SAMPLE}"
    exit 1
fi

# Copy results to final location instead of moving
echo "Copying results to final destination"
mkdir -p $OUTPUT_DIR/${SAMPLE}
cp -r ${SAMPLE}/outs $OUTPUT_DIR/${SAMPLE}/
echo "Results copied to $OUTPUT_DIR/${SAMPLE}/outs"

# Convert h5 to h5ad for Scanpy compatibility
if [ -f "${SAMPLE}/outs/filtered_feature_bc_matrix.h5" ]; then
    # Create and activate a Python environment with scanpy
    source \$(conda info --base)/etc/profile.d/conda.sh
    conda activate scanpy  # Your environment with scanpy installed
    
    # Run Python script to convert h5 to h5ad
    python -c "
import scanpy as sc
import os

# Load the Cell Ranger h5 file
adata = sc.read_10x_h5('${SAMPLE}/outs/filtered_feature_bc_matrix.h5')

# Add some metadata
adata.obs['sample'] = '${SAMPLE}'

# Save as h5ad
adata.write('${SAMPLE}/outs/${SAMPLE}.h5ad')
print(f'Converted {os.path.basename(\"${SAMPLE}\")} to h5ad format')

# Also save to final destination
adata.write('$OUTPUT_DIR/${SAMPLE}/outs/${SAMPLE}.h5ad')
print(f'Copied h5ad file to final destination')
"
    conda deactivate
fi

# List all files to confirm they exist
echo "Listing files in output directory:"
ls -la ${SAMPLE}/outs/
ls -la $OUTPUT_DIR/${SAMPLE}/outs/

echo "Cell Ranger count completed successfully for $SAMPLE"
date

# Cleanup the temporary working directory if copy was successful
if [ -f "$OUTPUT_DIR/${SAMPLE}/outs/${SAMPLE}.h5ad" ]; then
    echo "Cleaning up temporary work directory"
    cd $OUTPUT_DIR
    rm -rf \$WORK_DIR
else
    echo "Warning: Keeping work directory as final copy verification failed"
fi

exit 0
EOF

    # Make script executable
    chmod +x $TEMP_SCRIPT
    
    # Submit job
    sbatch $TEMP_SCRIPT
done

echo "All Cell Ranger count jobs submitted"
echo "Monitor job progress with: squeue -u \$USER"