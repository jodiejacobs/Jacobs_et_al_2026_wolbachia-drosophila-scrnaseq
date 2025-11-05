#!/usr/bin/env python3
import argparse
import gzip
import os

def load_t2g(t2g_file):
    """Load transcript-to-gene mapping"""
    transcript_to_gene = {}
    with open(t2g_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                transcript_to_gene[parts[0]] = parts[1]
    return transcript_to_gene

def load_transcripts(transcripts_file):
    """Load transcript list"""
    transcript_list = []
    # Handle both plain and gzipped
    try:
        with gzip.open(transcripts_file, 'rt') as f:
            transcript_list = [line.strip() for line in f]
    except:
        with open(transcripts_file, 'r') as f:
            transcript_list = [line.strip() for line in f]
    return transcript_list

def load_ec(ec_file):
    """Load equivalence classes"""
    ec_to_transcripts = {}
    with open(ec_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                ec_id = int(parts[0])
                transcript_indices = [int(x) for x in parts[1].split(',')]
                ec_to_transcripts[ec_id] = transcript_indices
    return ec_to_transcripts

def filter_bus(target_gene, t2g_file, transcripts_file, ec_file, bus_text_file, output_file):
    """Filter BUS file for target gene"""
    
    print(f"Filtering for gene: {target_gene}")
    print("=" * 60)
    
    # Load mappings
    print("[1/4] Loading transcript-to-gene mapping...")
    transcript_to_gene = load_t2g(t2g_file)
    print(f"  Loaded {len(transcript_to_gene):,} mappings")
    
    # Find target transcripts
    print(f"[2/4] Finding transcripts for {target_gene}...")
    target_transcripts = {tid for tid, gid in transcript_to_gene.items() if gid == target_gene}
    print(f"  Found {len(target_transcripts)} transcripts")
    if len(target_transcripts) == 0:
        print(f"ERROR: No transcripts found for gene {target_gene}")
        exit(1)
    
    # Load transcript list
    print("[3/4] Loading transcript list...")
    transcript_list = load_transcripts(transcripts_file)
    print(f"  Loaded {len(transcript_list):,} transcripts")
    
    # Get target indices
    target_indices = {i for i, tid in enumerate(transcript_list) if tid in target_transcripts}
    print(f"  Target indices: {sorted(target_indices)}")
    
    # Load ECs and find target ones
    ec_to_transcripts = load_ec(ec_file)
    print(f"  Loaded {len(ec_to_transcripts):,} equivalence classes")
    
    target_ecs = set()
    for ec_id, trans_indices in ec_to_transcripts.items():
        if any(idx in target_indices for idx in trans_indices):
            target_ecs.add(ec_id)
    
    print(f"  Found {len(target_ecs)} target ECs")
    if len(target_ecs) == 0:
        print(f"WARNING: No equivalence classes found for {target_gene}")
        print("This gene may not be expressed in your data.")
    
    # Filter BUS records
    print("[4/4] Filtering BUS records...")
    kept = 0
    total = 0
    
    with open(bus_text_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            total += 1
            parts = line.strip().split()
            if len(parts) >= 3:
                ec = int(parts[2])
                if ec in target_ecs:
                    outfile.write(line)
                    kept += 1
            
            if total % 5000000 == 0:
                print(f"  Processed {total:,} records, kept {kept:,}...")
    
    print(f"\n" + "=" * 60)
    print(f"RESULTS:")
    print(f"  Total records: {total:,}")
    print(f"  Kept records: {kept:,} ({100*kept/total:.3f}%)")
    print(f"  Saved to: {output_file}")
    print("=" * 60)
    
    return kept, total

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter BUS file by gene')
    parser.add_argument('--gene', required=True, help='Target gene ID')
    parser.add_argument('--t2g', required=True, help='Transcripts to genes file')
    parser.add_argument('--transcripts', required=True, help='Transcripts file')
    parser.add_argument('--ec', required=True, help='Equivalence class file')
    parser.add_argument('--bus_text', required=True, help='BUS text file')
    parser.add_argument('--output', required=True, help='Output filtered BUS text file')
    
    args = parser.parse_args()
    
    filter_bus(args.gene, args.t2g, args.transcripts, args.ec, args.bus_text, args.output)
