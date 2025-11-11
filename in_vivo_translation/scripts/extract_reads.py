# extract_reads.py
import gzip
import sys

# Load target barcode+UMI combinations
targets = set()
with open('GQX67_05945_records.txt') as f:
    for line in f:
        barcode, umi, ec, count = line.strip().split('\t')
        targets.add((barcode, umi))

print(f"Found {len(targets)} unique barcode+UMI combinations", file=sys.stderr)

# Get your original FASTQ files - adjust paths as needed
r1_path = '/private/groups/russelllab/jodie/sequencing_data/scRNAseq/PIPseq-primary_cells/ubkhc-wMel-2_*_R1_001.fastq.gz'
r2_path = '/private/groups/russelllab/jodie/sequencing_data/scRNAseq/PIPseq-primary_cells/ubkhc-wMel-2_*_R2_001.fastq.gz'

# Use glob to find files
import glob
r1_files = glob.glob(r1_path)
r2_files = glob.glob(r2_path)

if not r1_files or not r2_files:
    print(f"ERROR: Could not find FASTQ files", file=sys.stderr)
    sys.exit(1)

print(f"Reading from:\n  {r1_files[0]}\n  {r2_files[0]}", file=sys.stderr)

extracted = 0
with gzip.open(r1_files[0], 'rt') as r1, \
     gzip.open(r2_files[0], 'rt') as r2, \
     open('GQX67_05945_R1.fastq', 'w') as out1, \
     open('GQX67_05945_R2.fastq', 'w') as out2:
    
    while True:
        # Read 4 lines from each file (one FASTQ record)
        r1_lines = [r1.readline() for _ in range(4)]
        r2_lines = [r2.readline() for _ in range(4)]
        
        if not r1_lines[0]:  # EOF
            break
        
        # Extract barcode (0-16) and UMI (16-28) from R1 sequence
        seq = r1_lines[1].strip()
        barcode = seq[0:16]
        umi = seq[16:28]
        
        if (barcode, umi) in targets:
            out1.writelines(r1_lines)
            out2.writelines(r2_lines)
            extracted += 1

print(f"Extracted {extracted} read pairs", file=sys.stderr)
