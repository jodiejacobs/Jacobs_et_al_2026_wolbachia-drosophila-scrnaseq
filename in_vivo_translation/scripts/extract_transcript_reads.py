#!/usr/bin/env python3
"""
Extract FASTQ reads that mapped to a specific transcript from kallisto BUS output

Usage:
    # Option 1: Use transcript/gene name (searches t2g file)
    python extract_transcript_reads.py \
        --bus output.bus \
        --t2g transcripts_to_genes.txt \
        --gene GQX67_05945 \
        --r1 R1.fastq.gz \
        --r2 R2.fastq.gz \
        --output extracted_reads
    
    # Option 2: Use transcript index directly
    python extract_transcript_reads.py \
        --bus output.bus \
        --transcript-idx 1189 \
        --r1 R1.fastq.gz \
        --r2 R2.fastq.gz \
        --output extracted_reads
    
    # Option 3: Use already-converted BUS text file
    python extract_transcript_reads.py \
        --bus-text output.txt \
        --transcript-idx 1189 \
        --r1 R1.fastq.gz \
        --r2 R2.fastq.gz \
        --output extracted_reads
"""

import argparse
import gzip
import sys
import subprocess
from pathlib import Path


def get_transcript_index(t2g_file, transcript_name=None, gene_name=None):
    """Get the numeric index of a transcript from kallisto's t2g file"""
    print(f"Searching t2g file for transcript/gene...", file=sys.stderr)
    with open(t2g_file) as f:
        for idx, line in enumerate(f):
            fields = line.strip().split('\t')
            if len(fields) >= 2:
                t_id, g_id = fields[0], fields[1]
                if (transcript_name and t_id == transcript_name) or \
                   (gene_name and g_id == gene_name):
                    print(f"Found {t_id} -> {g_id} at index {idx}", file=sys.stderr)
                    return idx, t_id, g_id
    return None, None, None


def bus_to_text(bus_file, transcript_idx=None):
    """Convert BUS file to text, optionally filtering for a transcript"""
    import tempfile
    import os
    
    print("Converting BUS to text...", file=sys.stderr)
    
    # Create temporary file for bustools output
    temp_fd, temp_path = tempfile.mkstemp(suffix='.txt', prefix='bus_text_')
    os.close(temp_fd)
    
    try:
        # Run bustools text with -o flag
        cmd = ['bustools', 'text', '-o', temp_path, bus_file]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        # Parse the output file
        barcode_umi_set = set()
        total_lines = 0
        
        with open(temp_path, 'r') as f:
            for line in f:
                if not line.strip():
                    continue
                total_lines += 1
                fields = line.strip().split('\t')
                if len(fields) >= 4:
                    cb, umi, tid, count = fields[0], fields[1], int(fields[2]), int(fields[3])
                    if transcript_idx is None or tid == transcript_idx:
                        barcode_umi_set.add((cb, umi))
        
        print(f"Processed {total_lines:,} BUS records", file=sys.stderr)
        return barcode_umi_set
    
    except subprocess.CalledProcessError as e:
        print(f"Error running bustools: {e}", file=sys.stderr)
        print(f"stderr: {e.stderr}", file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError:
        print("Error: bustools not found. Make sure it's installed and in your PATH", file=sys.stderr)
        sys.exit(1)
    finally:
        # Clean up temporary file
        if os.path.exists(temp_path):
            os.remove(temp_path)


def parse_bus_text_for_transcript(bus_text_file, transcript_idx):
    """Extract barcode-UMI pairs for a specific transcript from text BUS file"""
    print(f"Parsing BUS text file for transcript index {transcript_idx}...", file=sys.stderr)
    barcode_umi_set = set()
    total_lines = 0
    
    with open(bus_text_file) as f:
        for line in f:
            total_lines += 1
            fields = line.strip().split()
            if len(fields) >= 4:
                cb, umi, tid, count = fields[0], fields[1], int(fields[2]), int(fields[3])
                if tid == transcript_idx:
                    barcode_umi_set.add((cb, umi))
    
    print(f"Processed {total_lines:,} BUS records", file=sys.stderr)
    return barcode_umi_set


def parse_fastq(filename):
    """Generator to parse FASTQ files (gzipped or not)"""
    open_func = gzip.open if filename.endswith('.gz') else open
    with open_func(filename, 'rt') as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline()
            plus = f.readline()
            qual = f.readline()
            yield header, seq, qual


def extract_reads(r1_file, r2_file, barcode_umi_set, output_prefix, barcode_len=16, umi_len=12, debug=False):
    """Extract reads matching the barcode-UMI set"""
    
    print(f"\nExtracting reads from FASTQs...", file=sys.stderr)
    print(f"  Barcode length: {barcode_len}", file=sys.stderr)
    print(f"  UMI length: {umi_len}", file=sys.stderr)
    print(f"  Target barcode-UMI pairs: {len(barcode_umi_set):,}", file=sys.stderr)
    
    # Show first few target pairs for debugging
    if debug and barcode_umi_set:
        print(f"\n  First few target barcode-UMI pairs:", file=sys.stderr)
        for i, (bc, umi) in enumerate(list(barcode_umi_set)[:5]):
            print(f"    {bc}\t{umi}", file=sys.stderr)
    
    out_r1 = gzip.open(f"{output_prefix}_R1.fastq.gz", 'wt')
    out_r2 = gzip.open(f"{output_prefix}_R2.fastq.gz", 'wt')
    
    r1_gen = parse_fastq(r1_file)
    r2_gen = parse_fastq(r2_file)
    
    extracted_count = 0
    total_count = 0
    
    umi_start = barcode_len
    umi_end = barcode_len + umi_len
    
    for (h1, s1, q1), (h2, s2, q2) in zip(r1_gen, r2_gen):
        total_count += 1
        
        # Extract barcode and UMI from R1 (remove newlines, don't strip whitespace from sequence)
        barcode = s1[:barcode_len].rstrip('\n')
        umi = s1[umi_start:umi_end].rstrip('\n')
        
        # Debug: show first few reads
        if debug and total_count <= 5:
            print(f"\n  Read {total_count}:", file=sys.stderr)
            print(f"    Barcode: {barcode}", file=sys.stderr)
            print(f"    UMI: {umi}", file=sys.stderr)
            print(f"    Match: {(barcode, umi) in barcode_umi_set}", file=sys.stderr)
        
        if (barcode, umi) in barcode_umi_set:
            out_r1.write(h1)
            out_r1.write(s1)
            out_r1.write("+\n")
            out_r1.write(q1)
            
            out_r2.write(h2)
            out_r2.write(s2)
            out_r2.write("+\n")
            out_r2.write(q2)
            
            extracted_count += 1
        
        if total_count % 1000000 == 0:
            print(f"  Processed {total_count:,} reads, extracted {extracted_count:,}", file=sys.stderr)
    
    out_r1.close()
    out_r2.close()
    
    print(f"\nExtracted {extracted_count:,} / {total_count:,} read pairs ({100*extracted_count/total_count:.2f}%)", file=sys.stderr)
    return extracted_count


def main():
    parser = argparse.ArgumentParser(
        description='Extract FASTQ reads that mapped to a specific transcript from kallisto BUS output',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    # Input files
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--bus', help='Input BUS file (binary)')
    input_group.add_argument('--bus-text', help='Input BUS file (already converted to text)')
    
    # Transcript identification
    transcript_group = parser.add_mutually_exclusive_group(required=True)
    transcript_group.add_argument('--transcript-idx', type=int, help='Transcript index (0-based)')
    transcript_group.add_argument('--transcript', help='Transcript name (requires --t2g)')
    transcript_group.add_argument('--gene', help='Gene name (requires --t2g)')
    
    parser.add_argument('--t2g', help='Transcripts to genes file (required if using --transcript or --gene)')
    
    # FASTQ files
    parser.add_argument('--r1', required=True, help='Input R1 FASTQ file (with barcodes/UMIs)')
    parser.add_argument('--r2', required=True, help='Input R2 FASTQ file (with reads)')
    
    # Output
    parser.add_argument('-o', '--output', required=True, help='Output prefix for extracted FASTQs')
    
    # Technology parameters
    parser.add_argument('--barcode-length', type=int, default=16, 
                       help='Barcode length in bp (default: 16 for 10xv3/PIPseq)')
    parser.add_argument('--umi-length', type=int, default=12,
                       help='UMI length in bp (default: 12 for 10xv3/PIPseq)')
    parser.add_argument('--debug', action='store_true',
                       help='Print debug information about first few reads')
    
    args = parser.parse_args()
    
    # Validate arguments
    if (args.transcript or args.gene) and not args.t2g:
        parser.error("--t2g is required when using --transcript or --gene")
    
    # Determine transcript index
    transcript_idx = args.transcript_idx
    
    if args.transcript or args.gene:
        transcript_idx, t_id, g_id = get_transcript_index(
            args.t2g, 
            transcript_name=args.transcript,
            gene_name=args.gene
        )
        if transcript_idx is None:
            print(f"Error: Could not find transcript/gene in t2g file", file=sys.stderr)
            sys.exit(1)
        print(f"Using transcript index: {transcript_idx} ({t_id} -> {g_id})", file=sys.stderr)
    else:
        print(f"Using transcript index: {transcript_idx}", file=sys.stderr)
    
    # Parse BUS file for target transcript
    if args.bus:
        barcode_umi_set = bus_to_text(args.bus, transcript_idx)
    else:
        barcode_umi_set = parse_bus_text_for_transcript(args.bus_text, transcript_idx)
    
    print(f"Found {len(barcode_umi_set):,} unique barcode-UMI pairs for transcript {transcript_idx}", file=sys.stderr)
    
    if len(barcode_umi_set) == 0:
        print("Warning: No reads found for this transcript!", file=sys.stderr)
        sys.exit(0)
    
    # Extract reads from FASTQs
    extracted = extract_reads(
        args.r1, 
        args.r2, 
        barcode_umi_set, 
        args.output,
        barcode_len=args.barcode_length,
        umi_len=args.umi_length,
        debug=args.debug
    )
    
    print(f"\nDone! Output files:", file=sys.stderr)
    print(f"  {args.output}_R1.fastq.gz", file=sys.stderr)
    print(f"  {args.output}_R2.fastq.gz", file=sys.stderr)


if __name__ == "__main__":
    main()