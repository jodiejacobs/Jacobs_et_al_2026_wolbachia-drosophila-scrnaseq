#!/usr/bin/env python3
"""
Filter BUS records by gene name.
Extracts all BUS records (barcode, UMI, EC) that map to a specific gene.
"""

import argparse
import sys
from collections import defaultdict

def parse_t2g(t2g_file, target_gene):
    """Find all transcripts that map to the target gene."""
    transcripts = []
    with open(t2g_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 2:
                continue
            transcript, gene = fields[0], fields[1]
            if gene == target_gene:
                transcripts.append(transcript)
    
    if not transcripts:
        print(f"ERROR: Gene '{target_gene}' not found in t2g file", file=sys.stderr)
        sys.exit(1)
    
    print(f"Found {len(transcripts)} transcript(s) for gene {target_gene}: {transcripts}", file=sys.stderr)
    return transcripts

def get_transcript_indices(transcripts_file, target_transcripts):
    """Get 0-based indices for target transcripts from transcripts.txt."""
    indices = []
    with open(transcripts_file, 'r') as f:
        for idx, line in enumerate(f):
            transcript = line.strip()
            if transcript in target_transcripts:
                indices.append(idx)
    
    if not indices:
        print(f"ERROR: No target transcripts found in transcripts.txt", file=sys.stderr)
        sys.exit(1)
    
    print(f"Transcript indices: {indices}", file=sys.stderr)
    return set(indices)

def find_target_ecs(ec_file, target_indices):
    """Find equivalence classes containing any of the target transcript indices."""
    target_ecs = set()
    
    with open(ec_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 2:
                continue
            
            ec_num = int(fields[0])
            transcript_list = [int(x) for x in fields[1].split(',')]
            
            # Check if any target transcript is in this EC
            if any(idx in target_indices for idx in transcript_list):
                target_ecs.add(ec_num)
    
    if not target_ecs:
        print(f"ERROR: No equivalence classes found for target transcripts", file=sys.stderr)
        sys.exit(1)
    
    print(f"Found {len(target_ecs)} equivalence classes containing target transcripts", file=sys.stderr)
    return target_ecs

def filter_bus_records(bus_text_file, target_ecs, output_file):
    """Filter BUS records for target equivalence classes."""
    records_found = 0
    
    with open(bus_text_file, 'r') as fin, open(output_file, 'w') as fout:
        for line in fin:
            fields = line.strip().split('\t')
            if len(fields) < 4:
                continue
            
            barcode, umi, ec, count = fields[0], fields[1], int(fields[2]), fields[3]
            
            if ec in target_ecs:
                fout.write(line)
                records_found += 1
    
    print(f"Extracted {records_found} BUS records for target gene", file=sys.stderr)
    return records_found

def main():
    parser = argparse.ArgumentParser(description='Filter BUS records by gene name')
    parser.add_argument('--gene', required=True, help='Target gene name')
    parser.add_argument('--t2g', required=True, help='Transcript-to-gene mapping file')
    parser.add_argument('--transcripts', required=True, help='transcripts.txt file')
    parser.add_argument('--ec', required=True, help='matrix.ec file')
    parser.add_argument('--bus_text', required=True, help='BUS text file (from bustools text)')
    parser.add_argument('--output', required=True, help='Output filtered BUS text file')
    
    args = parser.parse_args()
    
    print(f"Processing gene: {args.gene}", file=sys.stderr)
    
    # Step 1: Find transcripts for target gene
    target_transcripts = parse_t2g(args.t2g, args.gene)
    
    # Step 2: Get transcript indices
    target_indices = get_transcript_indices(args.transcripts, target_transcripts)
    
    # Step 3: Find ECs containing those transcripts
    target_ecs = find_target_ecs(args.ec, target_indices)
    
    # Step 4: Filter BUS records
    records_found = filter_bus_records(args.bus_text, target_ecs, args.output)
    
    if records_found == 0:
        print(f"WARNING: No reads found for gene {args.gene}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Successfully filtered {records_found} records for gene {args.gene}", file=sys.stderr)

if __name__ == '__main__':
    main()