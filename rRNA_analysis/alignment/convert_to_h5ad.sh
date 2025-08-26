#!/bin/bash
#SBATCH --job-name=pipseq_convert
#SBATCH --output=pipseq_convert_%A.log
#SBATCH --error=pipseq_convert_%A.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --partition=short

# Define base directories
RESULTS_DIR="/private/groups/russelllab/jodie/scRNAseq/pipseq/ovary_pipseq/results"
OUTPUT_DIR="/private/groups/russelllab/jodie/scRNAseq/pipseq/ovary_pipseq/h5ad_output"
SCRIPT_PATH="/private/groups/russelllab/jodie/scRNAseq/pipseq/20250325_convert_pipseq_to_h5ad.py"

# Create output directory
mkdir -p $OUTPUT_DIR

# Find all sample directories
cd $RESULTS_DIR
SAMPLES=$(find . -maxdepth 1 -mindepth 1 -type d -exec basename {} \;)

# Check if we found any samples
if [ -z "$SAMPLES" ]; then
    echo "ERROR: No sample directories found in $RESULTS_DIR"
    exit 1
fi

# Process each sample by submitting a separate job
for SAMPLE in $SAMPLES; do
    echo "Submitting job for sample: $SAMPLE"
    
    # Create a temporary submission script for this sample
    TEMP_SCRIPT="/tmp/convert_${SAMPLE}.sh"
    
    cat > $TEMP_SCRIPT << EOF
#!/bin/bash
#SBATCH --job-name=conv_${SAMPLE}
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --output=convert_${SAMPLE}_%j.log
#SBATCH --time=01:00:00

pwd; hostname; date

# Input and output file paths
INPUT_FILE="${RESULTS_DIR}/${SAMPLE}/molecule_info.h5"
OUTPUT_FILE="${OUTPUT_DIR}/${SAMPLE}.h5ad"

# Check if the input file exists
if [ ! -f "\$INPUT_FILE" ]; then
    echo "ERROR: Input file not found: \$INPUT_FILE"
    exit 1
fi

echo "Converting \$INPUT_FILE to \$OUTPUT_FILE"

# Run the conversion script
python ${SCRIPT_PATH} \\
    --input "\$INPUT_FILE" \\
    --output "\$OUTPUT_FILE" \\
    --min_cells 3 \\
    --min_genes 200 \\
    --compression gzip \\
    --compression_opts 4 \\
    --verbose

echo "Finished processing ${SAMPLE}"
date
EOF

    # Submit the job
    sbatch $TEMP_SCRIPT
    
    # Optional: remove the temporary script after submission
    # rm $TEMP_SCRIPT
done

echo "All jobs submitted"