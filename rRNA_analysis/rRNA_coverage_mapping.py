#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from Bio import AlignIO, SeqIO
import pysam
from collections import defaultdict
import subprocess
import sys
import os
import re

# python scripts/rRNA_coverage_mapping.py \
# /private/groups/russelllab/jodie/scRNAseq/pipseq/ovary_pipseq/results/mei_P26_uninf-1/star_out.bam \
# /private/groups/russelllab/jodie/scRNAseq/wolbachia_titer/rRNA_analysis/genome/Dmelanogaster_rRNA.fasta


def identify_rRNA_elements(bam_file):
    """Identify and classify rRNA elements in BAM file"""
    
    # Get all reference names
    result = subprocess.run(['samtools', 'idxstats', bam_file], 
                          capture_output=True, text=True)
    all_refs = [line.split('\t')[0] for line in result.stdout.strip().split('\n') 
                if line and not line.startswith('*')]
    
    # Classify rRNA elements
    rRNA_groups = {
        '18S': [],
        '28S': [],
        '5.8S': [],
        '5S': [],
        '16S_mt': [],
        '12S_mt': [],
        '16S_wMel': [],
        '23S_wMel': [],
        '5S_wMel': [],
        'other_rRNA': []
    }
    
    # Classification patterns
    patterns = {
        '18S': [r'18s', r'ssu', r'small.*ribosom', r'rRNA.*18'],
        '28S': [r'28s', r'lsu', r'large.*ribosom', r'rRNA.*28'],
        '5.8S': [r'5\.8s', r'5_8s'],
        '5S': [r'(?<!\.8)5s(?!.*wmel|.*wolbachia)', r'rRNA.*5(?!.*8)'],
        '16S_mt': [r'16s.*mt', r'16s.*mito', r'mitochondr.*16s'],
        '12S_mt': [r'12s.*mt', r'12s.*mito', r'mitochondr.*12s'],
        '16S_wMel': [r'16s.*wmel', r'16s.*wolbachia', r'wmel.*16s'],
        '23S_wMel': [r'23s.*wmel', r'23s.*wolbachia', r'wmel.*23s'],
        '5S_wMel': [r'5s.*wmel', r'5s.*wolbachia', r'wmel.*5s']
    }
    
    # Classify each reference
    for ref in all_refs:
        ref_lower = ref.lower()
        classified = False
        
        for rRNA_type, pattern_list in patterns.items():
            for pattern in pattern_list:
                if re.search(pattern, ref_lower):
                    rRNA_groups[rRNA_type].append(ref)
                    classified = True
                    break
            if classified:
                break
        
        # Check for generic rRNA patterns if not classified
        if not classified:
            if any(keyword in ref_lower for keyword in ['rrna', 'rdna', 'ribosom']):
                rRNA_groups['other_rRNA'].append(ref)
            # Include numbered rRNA gene references (like your 211000022278279 examples)
            elif ref.startswith('211000022') or 'rDNA' in ref:
                rRNA_groups['other_rRNA'].append(ref)
    
    # Remove empty groups
    rRNA_groups = {k: v for k, v in rRNA_groups.items() if v}
    
    return rRNA_groups

def get_coverage_for_references(bam_file, ref_list):
    """Get coverage data for a list of references"""
    
    all_coverage = []
    ref_lengths = []
    
    bamfile = pysam.AlignmentFile(bam_file, "rb")
    
    for ref in ref_list:
        try:
            # Get coverage for this reference
            coverage = list(bamfile.count_coverage(ref))
            total_coverage = np.sum(coverage, axis=0)
            all_coverage.append(total_coverage)
            ref_lengths.append(len(total_coverage))
        except Exception as e:
            print(f"Warning: Could not get coverage for {ref}: {e}")
            continue
    
    bamfile.close()
    return all_coverage, ref_lengths

def align_to_consensus(coverage_list, ref_lengths, consensus_length):
    """Align multiple coverage arrays to consensus length using interpolation"""
    
    aligned_coverage = np.zeros(consensus_length)
    total_weight = np.zeros(consensus_length)
    
    for coverage, ref_len in zip(coverage_list, ref_lengths):
        if len(coverage) == 0:
            continue
            
        # Create position mapping from reference to consensus
        if ref_len != consensus_length:
            # Interpolate coverage to match consensus length
            old_positions = np.linspace(0, consensus_length-1, ref_len)
            new_positions = np.arange(consensus_length)
            interpolated_coverage = np.interp(new_positions, old_positions, coverage)
        else:
            interpolated_coverage = coverage
        
        # Add to aligned coverage with equal weighting
        aligned_coverage += interpolated_coverage
        total_weight += 1
    
    # Average coverage across all references
    aligned_coverage = np.where(total_weight > 0, aligned_coverage / total_weight, 0)
    
    return aligned_coverage

def load_consensus_sequences(consensus_file):
    """Load consensus sequences from FASTA file"""
    
    consensus_seqs = {}
    
    try:
        for record in SeqIO.parse(consensus_file, "fasta"):
            # Try to match consensus sequence names to rRNA types
            seq_name = record.id.lower()
            
            if '18s' in seq_name:
                consensus_seqs['18S'] = (record.id, len(record.seq))
            elif '28s' in seq_name:
                consensus_seqs['28S'] = (record.id, len(record.seq))
            elif '5.8s' in seq_name or '5_8s' in seq_name:
                consensus_seqs['5.8S'] = (record.id, len(record.seq))
            elif '5s' in seq_name and 'wmel' not in seq_name:
                consensus_seqs['5S'] = (record.id, len(record.seq))
            elif '16s' in seq_name and 'mt' in seq_name:
                consensus_seqs['16S_mt'] = (record.id, len(record.seq))
            elif '12s' in seq_name and 'mt' in seq_name:
                consensus_seqs['12S_mt'] = (record.id, len(record.seq))
            elif '16s' in seq_name and 'wmel' in seq_name:
                consensus_seqs['16S_wMel'] = (record.id, len(record.seq))
            elif '23s' in seq_name and 'wmel' in seq_name:
                consensus_seqs['23S_wMel'] = (record.id, len(record.seq))
            elif '5s' in seq_name and 'wmel' in seq_name:
                consensus_seqs['5S_wMel'] = (record.id, len(record.seq))
            else:
                # Use the sequence name as-is for other sequences
                consensus_seqs[record.id] = (record.id, len(record.seq))
                
    except Exception as e:
        print(f"Error loading consensus file: {e}")
        return None
    
    return consensus_seqs

def plot_rRNA_coverage(bam_file, consensus_file=None, output_dir="rRNA_coverage_plots"):
    """Plot coverage for each rRNA type mapped to consensus sequences"""
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    print("Identifying rRNA elements in BAM file...")
    rRNA_groups = identify_rRNA_elements(bam_file)
    
    if not rRNA_groups:
        print("No rRNA elements found in BAM file!")
        return
    
    print(f"Found rRNA groups: {list(rRNA_groups.keys())}")
    for group, refs in rRNA_groups.items():
        print(f"  {group}: {len(refs)} references")
    
    # Load consensus sequences if provided
    consensus_seqs = None
    if consensus_file and os.path.exists(consensus_file):
        print(f"Loading consensus sequences from {consensus_file}")
        consensus_seqs = load_consensus_sequences(consensus_file)
    
    # Process each rRNA group
    all_results = {}
    
    for rRNA_type, ref_list in rRNA_groups.items():
        print(f"\nProcessing {rRNA_type} ({len(ref_list)} references)...")
        
        # Get coverage for all references in this group
        coverage_list, ref_lengths = get_coverage_for_references(bam_file, ref_list)
        
        if not coverage_list:
            print(f"No coverage data for {rRNA_type}")
            continue
        
        # Determine consensus length
        if consensus_seqs and rRNA_type in consensus_seqs:
            consensus_length = consensus_seqs[rRNA_type][1]
            consensus_name = consensus_seqs[rRNA_type][0]
        else:
            # Use median length of references
            consensus_length = int(np.median(ref_lengths))
            consensus_name = f"{rRNA_type}_consensus"
        
        # Align coverage to consensus
        aligned_coverage = align_to_consensus(coverage_list, ref_lengths, consensus_length)
        
        # Store results
        all_results[rRNA_type] = {
            'coverage': aligned_coverage,
            'consensus_length': consensus_length,
            'num_references': len(ref_list),
            'reference_names': ref_list,
            'consensus_name': consensus_name
        }
        
        # Create individual plot
        plt.figure(figsize=(12, 6))
        
        positions = np.arange(consensus_length)
        plt.plot(positions, aligned_coverage, linewidth=1, color='blue')
        plt.fill_between(positions, aligned_coverage, alpha=0.3, color='blue')
        
        plt.xlabel('Position in Consensus Sequence')
        plt.ylabel('Average Coverage Depth')
        plt.title(f'{rRNA_type} Coverage Profile\n'
                 f'({len(ref_list)} genomic copies, consensus length: {consensus_length:,} bp)')
        plt.grid(True, alpha=0.3)
        
        # Add statistics
        mean_cov = np.mean(aligned_coverage)
        max_cov = np.max(aligned_coverage)
        plt.text(0.02, 0.98, f'Mean: {mean_cov:.1f}x\nMax: {max_cov:.1f}x', 
                transform=plt.gca().transAxes, verticalalignment='top',
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
        
        plt.tight_layout()
        
        # Save individual plot
        plot_file = os.path.join(output_dir, f'{rRNA_type}_coverage.png')
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.show()
        
        print(f"Plot saved: {plot_file}")
    
    # Create summary comparison plot
    if len(all_results) > 1:
        create_summary_coverage_plot(all_results, output_dir)
    
    # Save coverage data
    save_coverage_data(all_results, output_dir)
    
    print(f"\nAnalysis complete! Results saved in: {output_dir}")
    return all_results

def create_summary_coverage_plot(all_results, output_dir):
    """Create a summary plot comparing all rRNA types"""
    
    fig, axes = plt.subplots(len(all_results), 1, figsize=(14, 3*len(all_results)))
    if len(all_results) == 1:
        axes = [axes]
    
    colors = plt.cm.Set3(np.linspace(0, 1, len(all_results)))
    
    for i, (rRNA_type, data) in enumerate(all_results.items()):
        coverage = data['coverage']
        consensus_length = data['consensus_length']
        num_refs = data['num_references']
        
        positions = np.arange(consensus_length)
        axes[i].plot(positions, coverage, linewidth=1, color=colors[i])
        axes[i].fill_between(positions, coverage, alpha=0.3, color=colors[i])
        
        axes[i].set_ylabel('Coverage')
        axes[i].set_title(f'{rRNA_type} ({num_refs} copies, {consensus_length:,} bp)')
        axes[i].grid(True, alpha=0.3)
        
        # Add mean coverage
        mean_cov = np.mean(coverage)
        axes[i].axhline(y=mean_cov, color='red', linestyle='--', alpha=0.7, 
                       label=f'Mean: {mean_cov:.1f}x')
        axes[i].legend()
    
    axes[-1].set_xlabel('Position in Consensus Sequence')
    plt.suptitle('rRNA Coverage Profiles - All Types', fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    summary_file = os.path.join(output_dir, 'rRNA_coverage_summary.png')
    plt.savefig(summary_file, dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"Summary plot saved: {summary_file}")

def save_coverage_data(all_results, output_dir):
    """Save coverage data to text files"""
    
    for rRNA_type, data in all_results.items():
        coverage = data['coverage']
        
        # Save as tab-delimited file
        output_file = os.path.join(output_dir, f'{rRNA_type}_coverage_data.txt')
        with open(output_file, 'w') as f:
            f.write("Position\tCoverage\n")
            for i, cov in enumerate(coverage):
                f.write(f"{i+1}\t{cov:.2f}\n")
        
        print(f"Coverage data saved: {output_file}")

def main():
    """Main function for command line execution"""
    
    if len(sys.argv) < 2:
        print("Usage: python rRNA_coverage_mapping.py <bam_file> [consensus_file] [output_dir]")
        print("\nExample:")
        print("  python rRNA_coverage_mapping.py alignment.bam")
        print("  python rRNA_coverage_mapping.py alignment.bam consensus_sequence.fasta")
        print("  python rRNA_coverage_mapping.py alignment.bam consensus.fasta rRNA_plots/")
        print("\nArguments:")
        print("  bam_file       : BAM alignment file (required)")
        print("  consensus_file : FASTA file with consensus sequences (optional)")
        print("  output_dir     : Output directory (default: rRNA_coverage_plots)")
        sys.exit(1)
    
    # Get arguments
    bam_file = sys.argv[1]
    consensus_file = sys.argv[2] if len(sys.argv) > 2 else None
    output_dir = sys.argv[3] if len(sys.argv) > 3 else "rRNA_coverage_plots"
    
    # Check if BAM file exists
    if not os.path.exists(bam_file):
        print(f"Error: BAM file '{bam_file}' not found!")
        sys.exit(1)
    
    # Check if BAM is indexed
    if not os.path.exists(bam_file + ".bai") and not os.path.exists(bam_file.replace(".bam", ".bai")):
        print(f"Warning: BAM index not found. Creating index...")
        try:
            subprocess.run(['samtools', 'index', bam_file], check=True)
        except subprocess.CalledProcessError:
            print("Failed to index BAM file. Please index manually with 'samtools index'")
            sys.exit(1)
    
    # Check consensus file if provided
    if consensus_file and not os.path.exists(consensus_file):
        print(f"Warning: Consensus file '{consensus_file}' not found. Using median lengths.")
        consensus_file = None
    
    print(f"Starting rRNA coverage analysis...")
    print(f"BAM file: {bam_file}")
    print(f"Consensus file: {consensus_file or 'None (using median lengths)'}")
    print(f"Output directory: {output_dir}")
    print("-" * 50)
    
    try:
        results = plot_rRNA_coverage(bam_file, consensus_file, output_dir)
        
        print("\n" + "=" * 50)
        print("COVERAGE ANALYSIS COMPLETE!")
        print(f"Analyzed {len(results)} rRNA types")
        print(f"Files saved in: {output_dir}/")
        print("=" * 50)
        
    except Exception as e:
        print(f"Error during analysis: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()