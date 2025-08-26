#!/bin/bash
# This is the main submission script that will create individual jobs

# Define the base directories
FASTQ_DIR="/private/groups/russelllab/jodie/sequencing_data/scRNAseq/PIP-seq-ovary/Jacobs_10706_250124B9"
STAR_INDEX="/private/groups/russelllab/jodie/scRNAseq/pipseq/references/Drosophila_melanogaster_wMel_combined_buildmapref"
OUTPUT_DIR="/private/groups/russelllab/jodie/scRNAseq/pipseq/ovary_pipseq/results"
PIPSEEKER="/private/home/jomojaco/pipseeker/pipseeker-v3.3.0-linux/pipseeker"

# Find all sample prefixes
cd $FASTQ_DIR
SAMPLES=$(ls *_L007_R1_001.fastq.gz | sed 's/_S[0-9]*_L007_R1_001.fastq.gz//' | sort | uniq)

# Process each sample by submitting a separate job
for SAMPLE in $SAMPLES; do
    echo "Submitting job for sample: $SAMPLE"
    
    # Create a temporary submission script for this sample
    TEMP_SCRIPT="/tmp/pipseeker_${SAMPLE}.sh"
    
    cat > $TEMP_SCRIPT << EOF
#!/bin/bash
#SBATCH --job-name=pip_${SAMPLE}
#SBATCH --partition=medium
#SBATCH --mail-user=jomojaco@ucsc.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --output=pipseeker_${SAMPLE}_%j.log
#SBATCH --time=04:00:00

pwd; hostname; date

# Create sample-specific output directory
SAMPLE_OUTPUT="${OUTPUT_DIR}/${SAMPLE}"
mkdir -p \$SAMPLE_OUTPUT

# Run PIPseeker for this sample
${PIPSEEKER} full \\
    --chemistry V \\
    --fastq ${FASTQ_DIR}/${SAMPLE} \\
    --star-index-path ${STAR_INDEX} \\
    --output-path \$SAMPLE_OUTPUT \\
    --threads 16

echo "Finished processing ${SAMPLE}"
date
EOF

    # Submit the job
    sbatch $TEMP_SCRIPT
    
    # Optional: remove the temporary script after submission
    # rm $TEMP_SCRIPT
done

echo "All jobs submitted"
