#!/bin/bash
#SBATCH --job-name=bam_sort
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --output=bam_sort_%j.out
#SBATCH --error=bam_sort_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jomojaco@ucsc.edu

# sbatch scripts/plot_subsample.sh $INPUT

REF_fasta="/private/groups/russelllab/jodie/scRNAseq/wolbachia_titer/rRNA_analysis/genome/rRNA_reference_Dmel_wMel.fasta"
REF_bed="/private/groups/russelllab/jodie/scRNAseq/wolbachia_titer/rRNA_analysis/Dmel_wMel_rRNA.gene.bed"
SCRIPT="/private/groups/russelllab/jodie/scRNAseq/wolbachia_titer/rRNA_analysis/scripts/plot_coverage.py"

INPUTBAM=$1
OUTDIR=$(dirname "$INPUTBAM")
SUB="0.0001 0.001 0.01"# 0.05 0.1 0.25 0.5

# Sort and index the original BAM if needed
BAM="${INPUTBAM%.bam}_sorted.bam"
if [ ! -f "${BAM}.bai" ]; then
    echo "Sorting and indexing original BAM..."
    samtools sort -@ $SLURM_CPUS_PER_TASK -o $BAM $INPUTBAM
    samtools index $BAM
fi

for S in $SUB; do  
    echo "Subsampling $S"

    # Use $S (not $SUB) in filenames
    READS="${INPUTBAM%.bam}_reads_${S}.bam"
    ALIGN="${INPUTBAM%.bam}_aligned_${S}.bam"
    COV="${INPUTBAM%.bam}_coverage_${S}.tsv"

    # Extract and subsample reads overlapping rRNA regions
    samtools view -b -L $REF_bed $BAM | samtools view -s $S -b > $READS

    # Convert BAM to FASTQ for minimap2, then align and sort
    samtools fastq $READS | \
    minimap2 -ax map-ont $REF_fasta - | \
    samtools view -Sb | samtools sort -o $ALIGN
    samtools index $ALIGN

    # Calculate depth
    samtools depth $ALIGN > $COV

    # Report the number of reads
    READ_COUNT=$(samtools view -c $ALIGN)
    echo "Subsampled $READ_COUNT reads for subsample $S"
    
    echo "Completed subsampling $S"
done



