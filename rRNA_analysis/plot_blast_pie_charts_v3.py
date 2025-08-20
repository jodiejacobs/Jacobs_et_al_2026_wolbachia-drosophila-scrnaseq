#!/bin/bash

# BLAST output summarizer
# Usage: ./summarize_blast.sh <blast_file> <blast_db_path>

if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <blast_file> <blast_db_path>"
    echo "Example: $0 results.blast /private/groups/russelllab/jodie/blast_db/16S_ribosomal_RNA_plus_wolbachia"
    exit 1
fi

blast_file="$1"
blast_db="$2"
sample_name=$(basename "$blast_file" .blast)
output_dir=$(dirname "$blast_file")

if [[ ! -f "$blast_file" ]]; then
    echo "Error: BLAST file $blast_file not found"
    exit 1
fi

echo "Summarizing BLAST results for $sample_name..."

# Create temporary file for organism data
temp_organisms=$(mktemp)

# Extract organism names and percent identities for ALL blast hits
echo "Looking up organism information..."
while read line; do
    subject_id=$(echo "$line" | cut -f2)
    percent_identity=$(echo "$line" | cut -f3)
    full_title=$(blastdbcmd -db "$blast_db" -entry "$subject_id" -outfmt "%t" 2>/dev/null | head -1)
    
    # Extract organism name from the full title
    if [[ -n "$full_title" ]]; then
        # Try to extract organism name - common patterns:
        # "Wolbachia pipientis strain wMel 16S ribosomal RNA, complete sequence"
        # "Magnetospirillum gryphiswaldense MSR-1 16S ribosomal RNA, partial sequence"
        organism=$(echo "$full_title" | sed -E 's/^([A-Za-z]+ [a-z]+).*$/\1/' | sed 's/ strain.*//; s/ 16S.*//; s/ rRNA.*//')
        
        # If extraction failed or looks like an accession, use the subject_id
        if [[ -z "$organism" || "$organism" =~ ^[A-Z]+_?[0-9] ]]; then
            organism="$subject_id"
        fi
    else
        organism="$subject_id"
    fi
    
    echo "$organism ($percent_identity% identity)" >> "$temp_organisms"
done < "$blast_file"

# Create main summary: sort, count, and display
echo "Creating summary..."
sort "$temp_organisms" | uniq -c | sort -nr > "${output_dir}/${sample_name}.blast.summary"

echo "Top organisms found in $sample_name:"
head -20 "${output_dir}/${sample_name}.blast.summary"

# Create detailed summary with average identity per organism
echo "Creating detailed summary with average identities..."
awk '{print $2 "\t" $3}' "$blast_file" | while read subject_id percent_identity; do
    full_title=$(blastdbcmd -db "$blast_db" -entry "$subject_id" -outfmt "%t" 2>/dev/null | head -1)
    
    # Extract organism name from the full title
    if [[ -n "$full_title" ]]; then
        organism=$(echo "$full_title" | sed -E 's/^([A-Za-z]+ [a-z]+).*$/\1/' | sed 's/ strain.*//; s/ 16S.*//; s/ rRNA.*//')
        if [[ -z "$organism" || "$organism" =~ ^[A-Z]+_?[0-9] ]]; then
            organism="$subject_id"
        fi
    else
        organism="$subject_id"
    fi
    
    echo -e "$organism\t$percent_identity"
done | awk '
{
    # Join all fields except the last one as organism name
    organism = $1; 
    for(i=2; i<NF; i++) organism = organism " " $i; 
    identity = $NF;
    
    sum[organism] += identity; 
    count[organism]++
} 
END {
    for(org in sum) 
        printf "%6d %s (%.1f%% avg identity)\n", count[org], org, sum[org]/count[org]
}' | sort -k1,1nr > "${output_dir}/${sample_name}.blast.detailed_summary"

# Clean up
rm "$temp_organisms"

echo ""
echo "Files created:"
echo "  ${output_dir}/${sample_name}.blast.summary - Simple count summary"
echo "  ${output_dir}/${sample_name}.blast.detailed_summary - Detailed with average identities"
echo ""
echo "BLAST analysis completed for $sample_name"
