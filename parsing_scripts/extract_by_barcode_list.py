#!/usr/bin/env python3
"""
Extract FASTQ reads by cell barcode list

Usage:
    python extract_by_barcode_list.py \
        --barcodes barcodes.txt \
        --r1 R1.fastq.gz \
        --r2 R2.fastq.gz \
        --output extracted_reads \
        --barcode-length 16
"""

import argparse
import gzip
import sys


def load_barcodes(barcode_file):
    """Load barcodes from file"""
    barcodes = set()
    with open(barcode_file) as f:
        for line in f:
            bc = line.strip()
            if bc:
                # Remove any suffix like -1
                bc = bc.split('-')[0]
                barcodes.add(bc)
    return barcodes


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


def extract_reads(r1_file, r2_file, barcodes, output_prefix, barcode_len=16):
    """Extract reads matching the barcode set"""
    
    print(f"Extracting reads from FASTQs...", file=sys.stderr)
    print(f"  Barcode length: {barcode_len}", file=sys.stderr)
    print(f"  Target barcodes: {len(barcodes):,}", file=sys.stderr)
    
    out_r1 = gzip.open(f"{output_prefix}_R1.fastq.gz", 'wt')
    out_r2 = gzip.open(f"{output_prefix}_R2.fastq.gz", 'wt')
    
    r1_gen = parse_fastq(r1_file)
    r2_gen = parse_fastq(r2_file)
    
    extracted_count = 0
    total_count = 0
    
    for (h1, s1, q1), (h2, s2, q2) in zip(r1_gen, r2_gen):
        total_count += 1
        
        # Extract barcode from R1
        barcode = s1[:barcode_len].rstrip('\n')
        
        if barcode in barcodes:
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
        description='Extract FASTQ reads by cell barcode list',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('--barcodes', required=True, 
                       help='File with cell barcodes (one per line)')
    parser.add_argument('--r1', required=True, 
                       help='Input R1 FASTQ file (with barcodes)')
    parser.add_argument('--r2', required=True, 
                       help='Input R2 FASTQ file (with reads)')
    parser.add_argument('-o', '--output', required=True, 
                       help='Output prefix for extracted FASTQs')
    parser.add_argument('--barcode-length', type=int, default=16,
                       help='Barcode length in bp (default: 16)')
    
    args = parser.parse_args()
    
    # Load barcodes
    print(f"Loading barcodes from {args.barcodes}...", file=sys.stderr)
    barcodes = load_barcodes(args.barcodes)
    print(f"Loaded {len(barcodes):,} unique barcodes", file=sys.stderr)
    
    if len(barcodes) == 0:
        print("Error: No barcodes loaded!", file=sys.stderr)
        sys.exit(1)
    
    # Extract reads
    extracted = extract_reads(
        args.r1, 
        args.r2, 
        barcodes, 
        args.output,
        barcode_len=args.barcode_length
    )
    
    print(f"\nDone! Output files:", file=sys.stderr)
    print(f"  {args.output}_R1.fastq.gz", file=sys.stderr)
    print(f"  {args.output}_R2.fastq.gz", file=sys.stderr)


if __name__ == "__main__":
    main()
