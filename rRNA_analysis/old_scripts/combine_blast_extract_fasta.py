#!/usr/bin/env python3

import pandas as pd
import subprocess
import os
import glob
from collections import defaultdict, Counter
from pathlib import Path

def run_blastdbcmd(blast_db, subject_ids, output_file):
    """Extract FASTA sequences using blastdbcmd"""
    try:
        # Write subject IDs to temporary file
        temp_file = f"{output_file}.temp_ids"
        with open(temp_file, 'w') as f:
            for sid in subject_ids:
                f.write(f"{sid}\n")
        
        # Run blastdbcmd
        cmd = [
            'blastdbcmd', '-db', blast_db, 
            '-entry_batch', temp_file, 
            '-outfmt', '%f'
        ]
        
        with open(output_file, 'w') as out_f:
            result = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True)
        
        # Clean up temp file
        os.remove(temp_file)
        
        if result.returncode == 0:
            return True
        else:
            print(f"blastdbcmd error: {result.stderr}")
            return False
            
    except Exception as e:
        print(f"Error running blastdbcmd: {e}")
        return False

def get_organism_name(blast_db, subject_id):
    """Get organism name for a subject ID"""
    try:
        cmd = ['blastdbcmd', '-db', blast_db, '-entry', subject_id, '-outfmt', '%t']
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode == 0:
            return result.stdout.strip().split('\n')[0]
        else:
            return "Unknown organism"
    except:
        return "Unknown organism"

def analyze_blast_results(blast_files_list, output_dir, blast_db, top_n=100):
    """Main function to analyze BLAST results"""
    
    print(f"=== Combining BLAST Results and Extracting Top {top_n} Hits ===")
    print(f"BLAST files list: {blast_files_list}")
    print(f"Output directory: {output_dir}")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Read blast files list
    blast_files = []
    with open(blast_files_list, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip empty lines and comments
            if not line or line.startswith('#'):
                continue
            if os.path.exists(line):
                blast_files.append(line)
            else:
                print(f"WARNING: File not found: {line}")
    
    if not blast_files:
        print(f"ERROR: No valid .blast files found in {blast_files_list}")
        return
    
    print(f"Found {len(blast_files)} valid BLAST result files")
    
    # BLAST column names
    columns = [
        'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'
    ]
    
    # Read and combine all BLAST results
    print("Reading and combining BLAST results...")
    all_results = []
    sample_info = {}
    
    for blast_file in blast_files:
        sample_name = os.path.splitext(os.path.basename(blast_file))[0]
        try:
            df = pd.read_csv(blast_file, sep='\t', header=None, names=columns)
            df['sample'] = sample_name
            all_results.append(df)
            sample_info[sample_name] = len(df)
            print(f"  {sample_name}: {len(df)} hits")
        except Exception as e:
            print(f"  ERROR reading {blast_file}: {e}")
    
    if not all_results:
        print("ERROR: No valid BLAST results found")
        return
    
    # Combine all results
    combined_df = pd.concat(all_results, ignore_index=True)
    total_hits = len(combined_df)
    print(f"Total BLAST hits across all samples: {total_hits}")
    
    # Save combined results
    combined_file = os.path.join(output_dir, "all_samples_combined.blast")
    combined_df[columns].to_csv(combined_file, sep='\t', header=False, index=False)
    print(f"Combined results saved to: {combined_file}")
    
    # Analyze subject ID frequencies
    print("Analyzing subject ID frequencies...")
    subject_counts = combined_df['sseqid'].value_counts()
    unique_subjects = len(subject_counts)
    print(f"Unique subject sequences: {unique_subjects}")
    
    # Get top N subjects
    top_subjects = subject_counts.head(top_n)
    print(f"Top {len(top_subjects)} most frequent subjects identified")
    
    # Calculate average identity for each top subject
    print("Calculating average identities...")
    avg_identities = combined_df.groupby('sseqid')['pident'].mean()
    
    # Get organism names for top subjects (this might take a while)
    print(f"Looking up organism information for top {top_n} hits...")
    organism_info = {}
    
    for i, (subject_id, count) in enumerate(top_subjects.items(), 1):
        print(f"  Processing {i}/{len(top_subjects)}: {subject_id}")
        organism_info[subject_id] = get_organism_name(blast_db, subject_id)
    
    # Create summary
    summary_data = []
    for subject_id, count in top_subjects.items():
        avg_identity = avg_identities[subject_id]
        organism = organism_info.get(subject_id, "Unknown organism")
        summary_data.append({
            'hit_count': count,
            'subject_id': subject_id,
            'organism': organism,
            'avg_identity': avg_identity
        })
    
    # Save summary
    summary_df = pd.DataFrame(summary_data)
    summary_file = os.path.join(output_dir, f"top_{top_n}_summary.txt")
    
    with open(summary_file, 'w') as f:
        f.write(f"# Top {top_n} BLAST hits across all samples\n")
        f.write("# Format: Hit_Count\tSubject_ID\tOrganism_Description\tAvg_Identity\n")
        for _, row in summary_df.iterrows():
            f.write(f"{row['hit_count']}\t{row['subject_id']}\t{row['organism']}\t{row['avg_identity']:.2f}%\n")
    
    print(f"Summary saved to: {summary_file}")
    
    # Extract FASTA sequences
    print(f"Extracting FASTA sequences for top {top_n} hits...")
    fasta_file = os.path.join(output_dir, f"top_{top_n}_sequences.fasta")
    
    if run_blastdbcmd(blast_db, top_subjects.index.tolist(), fasta_file):
        # Count extracted sequences
        with open(fasta_file, 'r') as f:
            seq_count = sum(1 for line in f if line.startswith('>'))
        print(f"Extracted {seq_count} FASTA sequences")
    else:
        print("ERROR: Failed to extract FASTA sequences")
        seq_count = 0
    
    # Create organism family analysis
    print("Creating organism family analysis...")
    family_data = []
    for _, row in summary_df.iterrows():
        organism = row['organism']
        if organism != "Unknown organism":
            # Extract genus and species (first two words)
            parts = organism.split()
            if len(parts) >= 2:
                genus_species = f"{parts[0]} {parts[1]}"
                family_data.append(genus_species)
    
    family_counts = Counter(family_data)
    family_file = os.path.join(output_dir, "organism_families.txt")
    with open(family_file, 'w') as f:
        f.write("=== Organism Family Analysis ===\n")
        for genus_species, count in family_counts.most_common():
            f.write(f"{count:6d} {genus_species}\n")
    
    # Summary statistics
    print("\n=== Analysis Complete ===")
    print(f"Files created:")
    print(f"  Combined BLAST results: {os.path.join(output_dir, 'all_samples_combined.blast')}")
    print(f"  Top {top_n} summary: {summary_file}")
    print(f"  Top {top_n} FASTA sequences: {fasta_file}")
    print(f"  Organism families: {family_file}")
    print(f"")
    print(f"Summary statistics:")
    print(f"  Total BLAST hits: {total_hits}")
    print(f"  Unique subjects: {unique_subjects}")
    print(f"  FASTA sequences extracted: {seq_count}")
    
    # Show top 5 results
    print("\n=== Top 5 Most Frequent Hits ===")
    for i, (_, row) in enumerate(summary_df.head(5).iterrows()):
        print(f"{i+1}. {row['organism']} ({row['hit_count']} hits, {row['avg_identity']:.1f}% avg identity)")
    
    return summary_df

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 3:
        print("Usage: python combine_blast_extract_fasta.py <blast_files.txt> <output_directory> [top_n]")
        print("Example: python combine_blast_extract_fasta.py blast_files.txt results/combined_analysis 100")
        sys.exit(1)
    
    blast_files_list = sys.argv[1]
    output_dir = sys.argv[2]
    top_n = int(sys.argv[3]) if len(sys.argv) > 3 else 100
    blast_db = "/private/groups/russelllab/jodie/blast_db/16S_ribosomal_RNA_plus_wolbachia"
    
    analyze_blast_results(blast_files_list, output_dir, blast_db, top_n)

