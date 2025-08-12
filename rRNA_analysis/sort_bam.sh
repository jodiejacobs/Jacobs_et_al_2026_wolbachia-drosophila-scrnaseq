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

# To run:
# mamba activate sra-tools
# for bam in /path/to/directory/*.bam; do
#     sbatch sort_bam.sh "$bam"
# done

REF_fasta="/private/groups/russelllab/jodie/scRNAseq/wolbachia_titer/rRNA_analysis/genome/rRNA_reference_Dmel_wMel.fasta"
REF_bed="/private/groups/russelllab/jodie/scRNAseq/wolbachia_titer/rRNA_analysis/Dmel_wMel_rRNA.gene.bed"

INPUTBAM=$1
BAM="${INPUTBAM%.bam}_sorted.bam"
READS="${INPUTBAM%.bam}_reads.bam"
ALIGN="${INPUTBAM%.bam}_aligned.bam"
COV="${INPUTBAM%.bam}_coverage.tsv"

samtools sort -@ $SLURM_CPUS_PER_TASK -o $BAM $INPUTBAM
samtools index $BAM

# Extract reads overlapping rRNA regions
samtools view -b -L $REF_bed $BAM > $READS

# Convert BAM to FASTQ for minimap2, then align and sort
samtools fastq $READS | \
minimap2 -ax map-ont $REF_fasta - | \
samtools view -Sb | samtools sort -o $ALIGN
samtools index $ALIGN

# Calculate depth
samtools depth $ALIGN > $COV