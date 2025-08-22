#!/bin/bash
#SBATCH --job-name=blast_16S
#SBATCH --time=24:00:00
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

fasta=$1

# Set BLAST database path
# export BLASTDB=/private/groups/russelllab/jodie/blast_db/

# Blast the reads aligned to the rRNA genes to look for bacterial contamination
echo "Processing $fasta"

# Extract sample name (everything before first underscore)
sample_name=$(basename "$fasta" | cut -d'_' -f1)
output_dir=$(dirname "$fasta")

blastn -db /private/groups/russelllab/jodie/blast_db/16S_ribosomal_RNA_plus_wolbachia \
        -query "$fasta" \
        -out "${output_dir}/${sample_name}.blast" \
        -num_threads $SLURM_CPUS_PER_TASK \
        -outfmt 6 \
        -max_target_seqs 5 \
        -evalue 1e-5

echo "Completed BLAST for $sample_name"

# Look up organism information for ALL hits with percent identity
echo "Top organisms found in $sample_name (with % identity):"

# Create temporary files
temp_organisms=$(mktemp)
temp_identities=$(mktemp)

# Extract organism names and percent identities for ALL blast hits
while read line; do
    subject_id=$(echo "$line" | cut -f2)
    percent_identity=$(echo "$line" | cut -f3)
    organism=$(blastdbcmd -db /private/groups/russelllab/jodie/blast_db/16S_ribosomal_RNA_plus_wolbachia -entry "$subject_id" -outfmt "%t" | head -1)
    echo "$organism ($percent_identity% identity)" >> "$temp_organisms"
done < "${output_dir}/${sample_name}.blast"

# Sort, count, and display
sort "$temp_organisms" | uniq -c | sort -nr | tee "${output_dir}/${sample_name}.blast.summary"

# Also create a detailed summary with average identity per organism
echo "Average percent identity by organism:" > "${output_dir}/${sample_name}.blast.detailed_summary"
awk '{print $2 "\t" $3}' "${output_dir}/${sample_name}.blast" | while read subject_id percent_identity; do
    organism=$(blastdbcmd -db /private/groups/russelllab/jodie/blast_db/16S_ribosomal_RNA_plus_wolbachia -entry "$subject_id" -outfmt "%t" | head -1)
    echo -e "$organism\t$percent_identity"
done | sort | awk '
{
    organism = $1; 
    for(i=2; i<NF; i++) organism = organism " " $i; 
    identity = $NF;
    sum[organism] += identity; 
    count[organism]++
} 
END {
    for(org in sum) 
        printf "%s\t%.1f%% avg identity\t%d reads\n", org, sum[org]/count[org], count[org]
}' | sort -k3 -nr >> "${output_dir}/${sample_name}.blast.detailed_summary"

# Clean up
rm "$temp_organisms" "$temp_identities"

echo ""
echo "BLAST analysis completed for $sample_name"