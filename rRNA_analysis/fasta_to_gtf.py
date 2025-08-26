#!/usr/bin/env python3

import sys
from Bio import SeqIO
import re

def fasta_to_gtf(fasta_file, output_gtf):
    """
    Convert FASTA file to GTF format for 16S rRNA sequences
    """
    
    # First pass: collect all sequence IDs and lengths from FASTA
    seq_info = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        # Use the exact sequence ID as it appears in the FASTA
        seq_id = record.id
        seq_length = len(record.seq)
        seq_info[seq_id] = {
            'length': seq_length,
            'description': record.description
        }
        print(f"Found sequence: '{seq_id}' with length {seq_length}")
    
    with open(output_gtf, 'w') as gtf_out:
        gene_counter = 1
        transcript_counter = 1
        
        for seq_id, info in seq_info.items():
            seq_length = info['length']
            description = info['description']
            
            # Extract organism name from description
            # Remove the sequence ID from description to get clean organism name
            organism_clean = description.replace(seq_id, "").strip()
            
            # Create clean locus tag from organism name
            # Extract genus and species for locus tag
            organism_parts = organism_clean.split()
            if len(organism_parts) >= 2:
                genus = organism_parts[0]
                species = organism_parts[1] if len(organism_parts) > 1 else "sp"
                locus_tag = f"{genus}_{species}_{gene_counter:05d}"
            else:
                locus_tag = f"rRNA_{gene_counter:05d}"
            
            # Clean locus tag (remove special characters)
            locus_tag = re.sub(r'[^a-zA-Z0-9_]', '_', locus_tag)
            
            # Gene ID and transcript ID
            gene_id = f"gene_{gene_counter:05d}"
            transcript_id = f"transcript_{transcript_counter:05d}"
            
            # Common attributes
            common_attrs = (
                f'gene_id "{gene_id}"; '
                f'locus_tag "{locus_tag}"; '
                f'product "16S ribosomal RNA"; '
                f'gene_biotype "rRNA"; '
                f'transcript_biotype "rRNA"'
            )
            
            transcript_attrs = (
                f'gene_id "{gene_id}"; '
                f'transcript_id "{transcript_id}"; '
                f'locus_tag "{locus_tag}"; '
                f'product "16S ribosomal RNA"; '
                f'gene_biotype "rRNA"; '
                f'transcript_biotype "rRNA"'
            )
            
            exon_attrs = (
                f'gene_id "{gene_id}"; '
                f'transcript_id "{transcript_id}"; '
                f'locus_tag "{locus_tag}"; '
                f'product "16S ribosomal RNA"; '
                f'gene_biotype "rRNA"; '
                f'transcript_biotype "rRNA"; '
                f'exon_number "1"'
            )
            
            # Write gene entry - use exact seq_id as chromosome name
            gtf_out.write(f"{seq_id}\tFASTA_derived\tgene\t1\t{seq_length}\t.\t+\t.\t{common_attrs};\n")
            
            # Write transcript entry
            gtf_out.write(f"{seq_id}\tFASTA_derived\ttranscript\t1\t{seq_length}\t.\t+\t.\t{transcript_attrs};\n")
            
            # Write exon entry
            gtf_out.write(f"{seq_id}\tFASTA_derived\texon\t1\t{seq_length}\t.\t+\t.\t{exon_attrs};\n")
            
            gene_counter += 1
            transcript_counter += 1

def validate_fasta_gtf_consistency(fasta_file, gtf_file):
    """
    Validate that GTF coordinates match FASTA sequence lengths
    """
    print("\nValidating FASTA-GTF consistency...")
    
    # Get FASTA lengths
    fasta_lengths = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        fasta_lengths[record.id] = len(record.seq)
    
    # Check GTF coordinates
    issues = []
    with open(gtf_file, 'r') as gtf:
        for line_num, line in enumerate(gtf, 1):
            if line.startswith('#') or line.strip() == '':
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 5:
                continue
                
            chrom = fields[0]
            end_pos = int(fields[4])
            
            if chrom in fasta_lengths:
                if end_pos > fasta_lengths[chrom]:
                    issues.append(f"Line {line_num}: {chrom} end={end_pos} > FASTA_length={fasta_lengths[chrom]}")
            else:
                issues.append(f"Line {line_num}: Chromosome '{chrom}' not found in FASTA")
    
    if issues:
        print("ISSUES FOUND:")
        for issue in issues[:10]:  # Show first 10 issues
            print(f"  {issue}")
        if len(issues) > 10:
            print(f"  ... and {len(issues) - 10} more issues")
        return False
    else:
        print("✓ All GTF coordinates are consistent with FASTA sequences")
        return True

def main():
    if len(sys.argv) not in [3, 4]:
        print("Usage: python fasta_to_gtf.py <input.fasta> <output.gtf> [--validate]")
        sys.exit(1)
    
    fasta_file = sys.argv[1]
    output_gtf = sys.argv[2]
    validate = len(sys.argv) == 4 and sys.argv[3] == "--validate"
    
    print(f"Converting {fasta_file} to GTF format...")
    fasta_to_gtf(fasta_file, output_gtf)
    print(f"GTF file created: {output_gtf}")
    
    if validate:
        validate_fasta_gtf_consistency(fasta_file, output_gtf)

if __name__ == "__main__":
    main()