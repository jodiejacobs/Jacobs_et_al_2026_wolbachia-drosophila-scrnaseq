#!/usr/bin/env python3

import sys
import os
import subprocess
import math
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq

def run_mafft(input_fasta, output_fasta):
    """Run MAFFT alignment and save to output."""
    with open(output_fasta, "w") as out_fh:
        subprocess.run(["mafft", "--auto", input_fasta], stdout=out_fh, check=True)

def calculate_info_content(alignment):
    """Calculate per-position information content (in bits)."""
    bases = ['A', 'C', 'G', 'T']
    results = []
    for i in range(alignment.get_alignment_length()):
        column = alignment[:, i].upper()
        total = sum(column.count(b) for b in bases)
        if total == 0:
            results.append(0.0)
            continue
        entropy = 0.0
        for b in bases:
            p = column.count(b) / total
            if p > 0:
                entropy -= p * math.log2(p)
        ic = 2.0 - entropy
        results.append(ic)
    return results

def sliding_a_content(sequence, window=50, step=1):
    """Calculate A-content in 50bp sliding windows."""
    seq = str(sequence).upper()
    results = []
    for i in range(0, len(seq) - window + 1, step):
        window_seq = seq[i:i+window]
        a_content = window_seq.count('A') / window
        results.append((i + 1, i + window, a_content))
    return results

def main(original_fasta, consensus_fasta):
    base = os.path.splitext(os.path.basename(original_fasta))[0]

    # Combine original + consensus
    combined_fasta = f"{base}_with_consensus.fasta"
    with open(combined_fasta, "w") as out_fh:
        for rec in SeqIO.parse(original_fasta, "fasta"):
            SeqIO.write(rec, out_fh, "fasta")
        for rec in SeqIO.parse(consensus_fasta, "fasta"):
            SeqIO.write(rec, out_fh, "fasta")

    # Align with MAFFT
    aligned_fasta = f"{base}_aligned.fasta"
    run_mafft(combined_fasta, aligned_fasta)

    # Info content from alignment
    aln = AlignIO.read(aligned_fasta, "fasta")
    info_content = calculate_info_content(aln)
    with open(f"{base}_info_content.txt", "w") as f:
        for i, val in enumerate(info_content, 1):
            f.write(f"{i}\t{val:.4f}\n")

    # A content from original consensus (not aligned)
    consensus_seq = SeqIO.read(consensus_fasta, "fasta").seq
    a_content_windows = sliding_a_content(consensus_seq)
    with open(f"{base}_Acontent_windows.txt", "w") as f:
        for start, end, val in a_content_windows:
            f.write(f"{start}\t{end}\t{val:.4f}\n")

    print(f"✅ Finished processing: {base}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python analyze_rRNA.py <original.fasta> <consensus.fasta>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])