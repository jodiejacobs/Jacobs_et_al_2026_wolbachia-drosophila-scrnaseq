#!/bin/bash
#SBATCH --job-name=blast_16S
#SBATCH --time=48:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=256
#SBATCH --output=logs/blast_16S_%j.out
#SBATCH --error=logs/blast_16S_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=long
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jomojaco@ucsc.edu

# To run:
# sbatch blast_16S.sh input.fasta output_directory
# As of 7/2/25 This is the final working version of the script

# Check arguments
if [ $# -ne 2 ]; then
    echo "Usage: sbatch blast_16S.sh <input.fasta> <output_directory>"
    exit 1
fi

# Set variables
tmp_dir="/tmp/blast_${SLURM_JOB_ID}"
fasta="$1"  # Input FASTA file
output_dir="$2"  # Output directory
sample_name=$(basename "$fasta" .fasta)
blast_db="/private/groups/russelllab/jodie/blast_db/16S_ribosomal_RNA_plus_wolbachia"

echo "=== BLAST 16S Analysis ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Input FASTA: $fasta"
echo "Sample name: $sample_name"
echo "Output directory: $output_dir"
echo "Temporary directory: $tmp_dir"
echo "Available CPUs: $SLURM_CPUS_PER_TASK"

# Check inputs
if [ ! -f "$fasta" ]; then
    echo "ERROR: Input FASTA file not found: $fasta"
    exit 1
fi

# Create directories
mkdir -p "$output_dir"
mkdir -p "$tmp_dir"
mkdir -p "${tmp_dir}/lookup_chunks"

# Check BLAST database
if [ ! -f "${blast_db}.nin" ] && [ ! -f "${blast_db}.00.nin" ]; then
    echo "ERROR: BLAST database not found: $blast_db"
    echo "Looking for database files..."
    ls -la $(dirname $blast_db)/*$(basename $blast_db)* 2>/dev/null || echo "No database files found"
    exit 1
fi

# Count input sequences
seq_count=$(grep -c "^>" "$fasta")
echo "Input sequences: $seq_count"

# Split FASTA into parts (use fewer parts for better parallelization)
num_parts=256  # Match your CPU count
echo "Splitting FASTA file into $num_parts parts..."
seqkit split -p $num_parts "$fasta" -O "$tmp_dir"

# Check if split worked
split_files=(${tmp_dir}/*.fasta)
if [ ${#split_files[@]} -eq 0 ]; then
    echo "ERROR: No split files created"
    exit 1
fi

echo "Created ${#split_files[@]} split files"

# Function to run BLAST on a single file
run_blast() {
    local split_file="$1"
    local part_name=$(basename "$split_file" .fasta)
    local blast_output="${tmp_dir}/${part_name}.blast"
    
    if [ -f "$split_file" ] && [ -s "$split_file" ]; then
        blastn -db "$blast_db" \
               -query "$split_file" \
               -out "$blast_output" \
               -num_threads 1 \
               -outfmt 6 \
               -max_target_seqs 5 \
               -evalue 1e-5 \
               -word_size 11
        
        # Check if output was created
        if [ -f "$blast_output" ]; then
            local hits=$(wc -l < "$blast_output" 2>/dev/null || echo "0")
            echo "Completed $part_name: $hits hits"
        else
            echo "Warning: No output for $part_name"
            touch "$blast_output"  # Create empty file
        fi
    else
        echo "Warning: Empty or missing split file: $split_file"
        touch "${tmp_dir}/${part_name}.blast"
    fi
}

# Export function for parallel execution
export -f run_blast
export tmp_dir blast_db

# Run BLAST on all parts in parallel
echo "Running BLAST on split files..."
printf '%s\n' "${split_files[@]}" | xargs -n 1 -P $SLURM_CPUS_PER_TASK -I {} bash -c 'run_blast "$@"' _ {}

# Wait a moment for all processes to finish writing
sleep 5

# Merge all BLAST results
echo "Merging BLAST results..."
blast_results=(${tmp_dir}/*.blast)
total_hits=0

# Check each result file
for blast_file in "${blast_results[@]}"; do
    if [ -f "$blast_file" ] && [ -s "$blast_file" ]; then
        hits=$(wc -l < "$blast_file")
        total_hits=$((total_hits + hits))
        echo "$(basename $blast_file): $hits hits"
    fi
done

# Merge results
if [ $total_hits -gt 0 ]; then
    cat ${tmp_dir}/*.blast > "${output_dir}/${sample_name}.blast"
    echo "Successfully merged $total_hits total hits"
    
    # Show summary of results
    echo "=== BLAST Results Summary ==="
    echo "Total hits: $total_hits"
    
    # Count unique queries with hits
    unique_queries=$(cut -f1 "${output_dir}/${sample_name}.blast" | sort -u | wc -l)
    echo "Queries with hits: $unique_queries / $seq_count"
    
    # Check for Wolbachia hits
    wolbachia_hits=$(grep -i wolbachia "${output_dir}/${sample_name}.blast" | wc -l)
    if [ $wolbachia_hits -gt 0 ]; then
        echo "Wolbachia hits found: $wolbachia_hits"
        echo "Top Wolbachia hits:"
        grep -i wolbachia "${output_dir}/${sample_name}.blast" | head -5 | cut -f1-4
    fi
    
    # Look up organism information for ALL hits with percent identity - PARALLELIZED VERSION
    echo "Top organisms found in $sample_name (with % identity):"
    
    # Split BLAST results into chunks for parallel organism lookup
    chunk_size=1000  # Process 1000 lines per chunk
    split -l $chunk_size "${output_dir}/${sample_name}.blast" "${tmp_dir}/lookup_chunks/blast_chunk_"
    
    # Function to lookup organisms for a chunk of BLAST results
    lookup_chunk() {
        local chunk_file="$1"
        local chunk_name=$(basename "$chunk_file")
        local output_file="${tmp_dir}/${chunk_name}_organisms.txt"
        
        # Process each line in the chunk
        while read line; do
            subject_id=$(echo "$line" | cut -f2)
            percent_identity=$(echo "$line" | cut -f3)
            organism=$(blastdbcmd -db "$blast_db" -entry "$subject_id" -outfmt "%t" 2>/dev/null | head -1)
            if [ -n "$organism" ]; then
                echo "$organism ($percent_identity% identity)"
            fi
        done < "$chunk_file" > "$output_file"
    }
    
    # Export function for parallel execution
    export -f lookup_chunk
    export blast_db tmp_dir
    
    # Run organism lookups in parallel
    lookup_chunks=(${tmp_dir}/lookup_chunks/blast_chunk_*)
    if [ ${#lookup_chunks[@]} -gt 0 ]; then
        # Use moderate parallelization to avoid overwhelming blastdbcmd
        lookup_procs=$((SLURM_CPUS_PER_TASK / 8))
        if [ $lookup_procs -lt 1 ]; then
            lookup_procs=1
        elif [ $lookup_procs -gt 32 ]; then
            lookup_procs=32
        fi
        
        echo "Running organism lookup with $lookup_procs parallel processes on ${#lookup_chunks[@]} chunks..."
        
        printf '%s\n' "${lookup_chunks[@]}" | xargs -n 1 -P $lookup_procs -I {} bash -c 'lookup_chunk "$@"' _ {}
        
        # Wait for completion
        sleep 2
        
        # Merge organism results and create summary
        echo "Merging organism lookup results..."
        cat ${tmp_dir}/lookup_chunks/*_organisms.txt > "${tmp_dir}/all_organisms.txt"
        
        # Sort, count, and display (same format as original)
        sort "${tmp_dir}/all_organisms.txt" | uniq -c | sort -nr | tee "${output_dir}/${sample_name}.blast.summary"
        
        # Create detailed summary with average identity per organism - PARALLELIZED VERSION
        echo "Average percent identity by organism:" > "${output_dir}/${sample_name}.blast.detailed_summary"
        
        # Split the identity calculation into chunks too
        split -l $chunk_size "${output_dir}/${sample_name}.blast" "${tmp_dir}/lookup_chunks/identity_chunk_"
        
        # Function to calculate identities for a chunk
        calc_identity_chunk() {
            local chunk_file="$1"
            local chunk_name=$(basename "$chunk_file")
            local output_file="${tmp_dir}/${chunk_name}_identities.txt"
            
            awk '{print $2 "\t" $3}' "$chunk_file" | while read subject_id percent_identity; do
                organism=$(blastdbcmd -db "$blast_db" -entry "$subject_id" -outfmt "%t" 2>/dev/null | head -1)
                if [ -n "$organism" ]; then
                    echo -e "$organism\t$percent_identity"
                fi
            done > "$output_file"
        }
        
        # Export function
        export -f calc_identity_chunk
        
        # Run identity calculations in parallel
        identity_chunks=(${tmp_dir}/lookup_chunks/identity_chunk_*)
        printf '%s\n' "${identity_chunks[@]}" | xargs -n 1 -P $lookup_procs -I {} bash -c 'calc_identity_chunk "$@"' _ {}
        
        # Wait and merge identity results
        sleep 2
        cat ${tmp_dir}/lookup_chunks/*_identities.txt | sort | awk '
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
        
    else
        echo "Warning: No lookup chunks created"
    fi
    
    echo ""
    echo "BLAST analysis completed for $sample_name"
    
else
    echo "No BLAST hits found"
    touch "${output_dir}/${sample_name}.blast"
fi

# Clean up temporary files
echo "Cleaning up temporary directory..."
rm -rf "$tmp_dir"

echo "=== BLAST analysis complete ==="
echo "Results saved to: ${output_dir}/${sample_name}.blast"