#!/bin/bash
#SBATCH --job-name=combine_blast
#SBATCH --time=06:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=16
#SBATCH --output=bam_sort_%j.out
#SBATCH --error=bam_sort_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jomojaco@ucsc.edu

# Script to combine BLAST results, summarize top hits, and extract FASTA sequences
# Usage: bash combine_blast_extract_fasta.sh blast_files.txt output_directory [top_n]

#ex. sbatch scripts/combine_blast_extract_fasta.sh scripts/blast_files.txt results/combined_analysis 100


# Check arguments
if [ $# -lt 2 ]; then
    echo "Usage: bash combine_blast_extract_fasta.sh <blast_files.txt> <output_directory> [top_n]"
    echo "Example: bash combine_blast_extract_fasta.sh blast_files.txt results/combined_analysis 100"
    exit 1
fi

# Set variables
blast_files_list="$1"
output_dir="$2"
blast_db="/private/groups/russelllab/jodie/blast_db/16S_ribosomal_RNA_plus_wolbachia"
top_n=${3:-100}  # Default to 100 if not specified

echo "=== Combining BLAST Results and Extracting Top Hits ==="
echo "BLAST files list: $blast_files_list"
echo "Output directory: $output_dir"
echo "Top N hits to extract: $top_n"

# Check if blast files list exists
if [ ! -f "$blast_files_list" ]; then
    echo "ERROR: BLAST files list not found: $blast_files_list"
    exit 1
fi

# Create output directory
mkdir -p "$output_dir"

# Read and validate blast files
echo "Reading BLAST files from: $blast_files_list"
blast_files=()
while IFS= read -r line || [ -n "$line" ]; do
    # Skip empty lines and comments
    [[ -z "$line" || "$line" =~ ^[[:space:]]*# ]] && continue
    
    # Check if file exists
    if [ -f "$line" ]; then
        blast_files+=("$line")
    else
        echo "WARNING: File not found: $line"
    fi
done < "$blast_files_list"

file_count=${#blast_files[@]}

if [ $file_count -eq 0 ]; then
    echo "ERROR: No valid .blast files found in $blast_files_list"
    exit 1
fi

echo "Found $file_count valid BLAST result files"

# Combine all BLAST results
echo "Combining all BLAST results..."
combined_blast="${output_dir}/all_samples_combined.blast"
cat "${blast_files[@]}" > "$combined_blast"

total_hits=$(wc -l < "$combined_blast")
echo "Total BLAST hits across all samples: $total_hits"

# Extract unique subject IDs and their hit counts
echo "Analyzing subject ID frequencies..."
subject_counts="${output_dir}/subject_hit_counts.txt"
cut -f2 "$combined_blast" | sort | uniq -c | sort -nr > "$subject_counts"

unique_subjects=$(wc -l < "$subject_counts")
echo "Unique subject sequences: $unique_subjects"

# Get top N subject IDs
echo "Extracting top $top_n subject IDs..."
top_subjects="${output_dir}/top_${top_n}_subjects.txt"
head -n "$top_n" "$subject_counts" | awk '{print $2}' > "$top_subjects"

echo "Top 10 most frequent hits:"
head -10 "$subject_counts"

# Create detailed summary with organism names for top hits
echo "Looking up organism information for top hits..."
top_summary="${output_dir}/top_${top_n}_summary.txt"

echo "# Top $top_n BLAST hits across all samples" > "$top_summary"
echo "# Format: Hit_Count Subject_ID Organism_Description Avg_Identity" >> "$top_summary"

# Process each top subject ID
counter=0
while read subject_id; do
    counter=$((counter + 1))
    echo "Processing $counter/$top_n: $subject_id"
    
    # Get hit count for this subject
    hit_count=$(grep -w "$subject_id" "$subject_counts" | awk '{print $1}')
    
    # Get organism description
    organism=$(blastdbcmd -db "$blast_db" -entry "$subject_id" -outfmt "%t" 2>/dev/null | head -1)
    if [ -z "$organism" ]; then
        organism="Unknown organism"
    fi
    
    # Calculate average percent identity for this subject
    avg_identity=$(awk -v sid="$subject_id" '$2 == sid {sum += $3; count++} END {if(count > 0) printf "%.2f", sum/count; else print "0"}' "$combined_blast")
    
    # Write to summary
    echo -e "${hit_count}\t${subject_id}\t${organism}\t${avg_identity}%" >> "$top_summary"
    
done < "$top_subjects"

echo "Summary created: $top_summary"

# Extract FASTA sequences for top hits
echo "Extracting FASTA sequences for top $top_n hits..."
top_fasta="${output_dir}/top_${top_n}_sequences.fasta"

# Use blastdbcmd to extract sequences
blastdbcmd -db "$blast_db" -entry_batch "$top_subjects" -outfmt "%f" > "$top_fasta" 2>/dev/null

# Check if FASTA extraction worked
seq_count=$(grep -c "^>" "$top_fasta" 2>/dev/null || echo "0")
echo "Extracted $seq_count FASTA sequences"

if [ $seq_count -eq 0 ]; then
    echo "WARNING: No sequences extracted. Trying alternative method..."
    
    # Alternative method: extract one by one
    > "$top_fasta"  # Clear file
    while read subject_id; do
        blastdbcmd -db "$blast_db" -entry "$subject_id" -outfmt "%f" 2>/dev/null >> "$top_fasta"
    done < "$top_subjects"
    
    seq_count=$(grep -c "^>" "$top_fasta" 2>/dev/null || echo "0")
    echo "Alternative extraction: $seq_count sequences"
fi

# Create additional analyses
echo "Creating organism family summary..."
family_summary="${output_dir}/organism_families.txt"
echo "=== Organism Family Analysis ===" > "$family_summary"

# Extract and count organism families/genera
awk -F'\t' '{if(NF >= 3 && $3 != "Unknown organism") print $3}' "$top_summary" | \
sed 's/strain.*//' | sed 's/16S.*//' | \
awk '{print $1, $2}' | sort | uniq -c | sort -nr >> "$family_summary"

echo ""
echo "=== Analysis Complete ==="
echo "Files created:"
echo "  Combined BLAST results: $combined_blast"
echo "  Subject hit counts: $subject_counts"
echo "  Top $top_n summary: $top_summary"
echo "  Top $top_n FASTA sequences: $top_fasta"
echo "  Organism families: $family_summary"
echo ""
echo "Summary statistics:"
echo "  Total BLAST hits: $total_hits"
echo "  Unique subjects: $unique_subjects"
echo "  FASTA sequences extracted: $seq_count"

# Show quick preview of results
echo ""
echo "=== Top 5 Most Frequent Hits ==="
head -5 "$top_summary"
