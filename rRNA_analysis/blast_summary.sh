blast_db="/private/groups/russelllab/jodie/blast_db/16S_ribosomal_RNA_plus_wolbachia"

for summary_file in results/blast_16S/10X/*.blast.summary; do
    echo "Processing $(basename "$summary_file")..."
    
    # Create backup
    cp "$summary_file" "$summary_file.backup"
    
    # Process each line and replace accession with organism name
    sed -E 's/^( *[0-9]+) ([A-Z_0-9.]+) (\([0-9.]+% identity\))$/\1 \2 \3/' "$summary_file.backup" | \
    while read count subject_id rest; do
        if [[ -n "$subject_id" && "$subject_id" != "" ]]; then
            full_title=$(blastdbcmd -db "$blast_db" -entry "$subject_id" -outfmt "%t" 2>/dev/null | head -1)
            if [[ -n "$full_title" ]]; then
                organism=$(echo "$full_title" | sed -E 's/^([A-Za-z]+ [a-z]+).*$/\1/' | sed 's/ strain.*//; s/ 16S.*//; s/ rRNA.*//')
                if [[ -z "$organism" || "$organism" =~ ^[A-Z]+_?[0-9] ]]; then
                    organism="$subject_id"
                fi
            else
                organism="$subject_id"
            fi
            echo "$count $organism $rest"
        fi
    done > "$summary_file"
done
