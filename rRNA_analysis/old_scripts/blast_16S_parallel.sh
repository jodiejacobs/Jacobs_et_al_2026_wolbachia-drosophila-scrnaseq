#!/bin/bash
#SBATCH --job-name=blast_16S
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=64
#SBATCH --output=logs/blast_16S_%j.out
#SBATCH --error=logs/blast_16S_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=long
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jomojaco@ucsc.edu

# To run:
# sbatch blast_16S.sh fasta

# Set variables
tmp_dir="/tmp/blast_${SLURM_JOB_ID}"
fasta="$1"  # Input FASTA file
sample_name=$(basename "$fasta" .fasta)
output_dir="$2"  # Output directory
blast_db="/private/groups/russelllab/jodie/blast_db/16S_ribosomal_RNA_plus_wolbachia"

# Create temporary directory
mkdir -p $tmp_dir

# Split FASTA into 100 parts
echo "Splitting FASTA file into parts..."
seqkit split -p 100 $fasta -O $tmp_dir

# Run BLAST on each part in parallel
echo "Running BLAST on split files..."
for split_file in ${tmp_dir}/*.fasta; do
    part_name=$(basename "$split_file" .fasta)
    blastn -db $blast_db \
           -query "$split_file" \
           -out "${tmp_dir}/${part_name}.blast" \
           -num_threads 1 \
           -outfmt 6 \
           -max_target_seqs 5 \
           -evalue 1e-5 &
done

# Wait for all BLAST jobs to complete
wait

# Merge all BLAST results
echo "Merging BLAST results..."
cat ${tmp_dir}/*.blast > "${output_dir}/${sample_name}.blast"

# Clean up temporary files
rm -rf $tmp_dir

echo "BLAST analysis complete. Results saved to ${output_dir}/${sample_name}.blast"