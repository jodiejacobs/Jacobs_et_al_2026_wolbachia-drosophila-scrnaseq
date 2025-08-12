#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from Bio import AlignIO
import pysam
from collections import defaultdict
import subprocess
import os
import sys

def get_nucleotide_counts_from_bam(bam_file, reference_name):
    """Extract nucleotide counts at each position from BAM file"""
    bamfile = pysam.AlignmentFile(bam_file, "rb")
    
    # Get reference length
    ref_length = bamfile.get_reference_length(reference_name)
    
    # Initialize count matrix: [position][nucleotide]
    counts = defaultdict(lambda: defaultdict(int))
    
    # Iterate through alignments
    for read in bamfile.fetch(reference_name):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
            
        # Check if query_sequence exists
        if read.query_sequence is None:
            continue
            
        # Get aligned pairs (query_pos, ref_pos) - without requiring MD tag
        aligned_pairs = read.get_aligned_pairs(matches_only=False, with_seq=False)
        
        for query_pos, ref_pos in aligned_pairs:
            if ref_pos is None:  # Insertion in query
                continue
            if query_pos is None:  # Deletion in query
                counts[ref_pos]['D'] += 1  # Deletion
                continue
                
            # Get the actual base from the read
            if query_pos < len(read.query_sequence):
                read_base = read.query_sequence[query_pos].upper()
                counts[ref_pos][read_base] += 1
    
    bamfile.close()
    return counts, ref_length

def calculate_information_content_from_bam(bam_file, reference_name):
    """Calculate information content using BAM file nucleotide counts"""
    counts, ref_length = get_nucleotide_counts_from_bam(bam_file, reference_name)
    
    ic_scores = []
    base_counts = []
    
    for pos in range(ref_length):
        pos_counts = counts[pos]
        
        # Get total coverage (excluding deletions for IC calculation)
        total = sum(count for base, count in pos_counts.items() if base != 'D')
        
        if total == 0:
            ic_scores.append(0)
            base_counts.append({'A': 0, 'T': 0, 'G': 0, 'C': 0, 'total': 0})
            continue
        
        # Calculate frequencies for ATGC only
        frequencies = []
        pos_base_counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
        
        for base in ['A', 'T', 'G', 'C']:
            count = pos_counts.get(base, 0)
            pos_base_counts[base] = count
            if count > 0:
                freq = count / total
                frequencies.append(freq)
        
        pos_base_counts['total'] = total
        base_counts.append(pos_base_counts)
        
        # Calculate information content: IC = 2 + sum(p_i * log2(p_i))
        if len(frequencies) > 0:
            ic = 2 + sum(f * np.log2(f + 1e-10) for f in frequencies)
            ic_scores.append(max(0, ic))
        else:
            ic_scores.append(0)
    
    return ic_scores, base_counts

def plot_a_content_windows_from_bam(bam_file, reference_name, window_size=50):
    """Calculate A content in windows using BAM counts"""
    ic_scores, base_counts = calculate_information_content_from_bam(bam_file, reference_name)
    
    a_content = []
    positions = []
    
    for i in range(0, len(base_counts) - window_size + 1, window_size//2):
        window_end = min(i + window_size, len(base_counts))
        
        total_a = sum(base_counts[j]['A'] for j in range(i, window_end))
        total_bases = sum(base_counts[j]['total'] for j in range(i, window_end))
        
        if total_bases > 0:
            a_fraction = total_a / total_bases
            a_content.append(a_fraction * 100)
            positions.append(i + window_size//2)
    
    return positions, a_content

def create_individual_rRNA_report(bam_file, reference_name, window_size=50, output_dir="rRNA_reports"):
    """Create comprehensive report for a single rRNA element"""
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"Analyzing {reference_name}...")
    
    try:
        # Calculate information content and base counts
        ic_scores, base_counts = calculate_information_content_from_bam(bam_file, reference_name)
        
        # Get coverage
        bamfile = pysam.AlignmentFile(bam_file, "rb")
        coverage = list(bamfile.count_coverage(reference_name))
        total_coverage = np.sum(coverage, axis=0)
        bamfile.close()
        
        # Calculate A content in windows
        a_positions, a_content = plot_a_content_windows_from_bam(bam_file, reference_name, window_size)
        
        # Get nucleotide composition data
        positions = list(range(len(base_counts)))
        a_counts = [base_counts[i]['A'] for i in positions]
        t_counts = [base_counts[i]['T'] for i in positions]
        g_counts = [base_counts[i]['G'] for i in positions]
        c_counts = [base_counts[i]['C'] for i in positions]
        total_counts = [base_counts[i]['total'] for i in positions]
        
        # Calculate percentages
        a_pct = [a_counts[i]/total_counts[i]*100 if total_counts[i] > 0 else 0 for i in range(len(positions))]
        t_pct = [t_counts[i]/total_counts[i]*100 if total_counts[i] > 0 else 0 for i in range(len(positions))]
        g_pct = [g_counts[i]/total_counts[i]*100 if total_counts[i] > 0 else 0 for i in range(len(positions))]
        c_pct = [c_counts[i]/total_counts[i]*100 if total_counts[i] > 0 else 0 for i in range(len(positions))]
        
        # Create comprehensive figure
        fig = plt.figure(figsize=(16, 12))
        
        # Main analysis (3 subplots)
        gs = fig.add_gridspec(4, 2, height_ratios=[1, 1, 1, 1], hspace=0.3, wspace=0.3)
        
        # 1. Information content
        ax1 = fig.add_subplot(gs[0, :])
        ax1.plot(range(len(ic_scores)), ic_scores, linewidth=1, color='green')
        ax1.set_ylabel('Information Content (bits)')
        ax1.set_title(f'{reference_name}: Information Content from BAM reads')
        ax1.grid(True, alpha=0.3)
        ax1.set_ylim(0, 2.1)
        
        # 2. Coverage
        ax2 = fig.add_subplot(gs[1, :])
        ax2.plot(range(len(total_coverage)), total_coverage, linewidth=1, color='blue')
        ax2.fill_between(range(len(total_coverage)), total_coverage, alpha=0.3, color='blue')
        ax2.set_ylabel('Coverage Depth')
        ax2.set_title(f'{reference_name}: Coverage Depth')
        ax2.grid(True, alpha=0.3)
        
        # 3. A content
        ax3 = fig.add_subplot(gs[2, :])
        ax3.plot(a_positions, a_content, linewidth=2, marker='o', markersize=2, color='red')
        ax3.axhline(y=25, color='gray', linestyle='--', alpha=0.5, label='Expected (25%)')
        ax3.set_ylabel('A Content (%)')
        ax3.set_title(f'{reference_name}: A Content per {window_size}bp Window')
        ax3.grid(True, alpha=0.3)
        ax3.legend()
        
        # 4. Nucleotide composition (left)
        ax4 = fig.add_subplot(gs[3, 0])
        ax4.fill_between(positions, 0, a_pct, alpha=0.7, color='red', label='A')
        ax4.fill_between(positions, a_pct, [a_pct[i] + t_pct[i] for i in range(len(positions))], 
                         alpha=0.7, color='blue', label='T')
        ax4.fill_between(positions, [a_pct[i] + t_pct[i] for i in range(len(positions))], 
                         [a_pct[i] + t_pct[i] + g_pct[i] for i in range(len(positions))], 
                         alpha=0.7, color='green', label='G')
        ax4.fill_between(positions, [a_pct[i] + t_pct[i] + g_pct[i] for i in range(len(positions))], 
                         [a_pct[i] + t_pct[i] + g_pct[i] + c_pct[i] for i in range(len(positions))], 
                         alpha=0.7, color='orange', label='C')
        ax4.set_xlabel('Position')
        ax4.set_ylabel('Nucleotide Composition (%)')
        ax4.set_title(f'{reference_name}: Base Composition')
        ax4.legend()
        ax4.set_ylim(0, 100)
        
        # 5. Summary statistics (right)
        ax5 = fig.add_subplot(gs[3, 1])
        ax5.axis('off')
        
        # Calculate summary stats
        mean_coverage = np.mean(total_coverage)
        median_coverage = np.median(total_coverage)
        mean_ic = np.mean(ic_scores)
        mean_a_content = np.mean(a_content) if a_content else 0
        total_length = len(positions)
        high_ic_positions = sum(1 for ic in ic_scores if ic > 1.5)
        
        stats_text = f"""
        SUMMARY STATISTICS
        
        Length: {total_length:,} bp
        
        Coverage:
        • Mean: {mean_coverage:.1f}x
        • Median: {median_coverage:.1f}x
        • Max: {max(total_coverage):,}x
        
        Information Content:
        • Mean IC: {mean_ic:.3f} bits
        • Highly conserved positions (IC>1.5): {high_ic_positions}
        • Conservation rate: {high_ic_positions/total_length*100:.1f}%
        
        Composition:
        • Mean A content: {mean_a_content:.1f}%
        • Overall A%: {sum(a_counts)/sum(total_counts)*100:.1f}%
        • Overall T%: {sum(t_counts)/sum(total_counts)*100:.1f}%
        • Overall G%: {sum(g_counts)/sum(total_counts)*100:.1f}%
        • Overall C%: {sum(c_counts)/sum(total_counts)*100:.1f}%
        """
        
        ax5.text(0.05, 0.95, stats_text, transform=ax5.transAxes, fontsize=10,
                verticalalignment='top', fontfamily='monospace',
                bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray", alpha=0.8))
        
        plt.suptitle(f'{reference_name} - Comprehensive rRNA Analysis Report', fontsize=16, fontweight='bold')
        
        # Save figure
        output_file = os.path.join(output_dir, f'{reference_name}_report.png')
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close(fig)  # Close figure to save memory
        
        print(f"Report saved: {output_file}")
        
        # Return data for further analysis if needed
        return {
            'reference_name': reference_name,
            'ic_scores': ic_scores,
            'coverage': total_coverage,
            'a_content': a_content,
            'base_counts': base_counts,
            'summary_stats': {
                'length': total_length,
                'mean_coverage': mean_coverage,
                'median_coverage': median_coverage,
                'mean_ic': mean_ic,
                'mean_a_content': mean_a_content,
                'conservation_rate': high_ic_positions/total_length*100
            }
        }
        
    except Exception as e:
        print(f"Error analyzing {reference_name}: {e}")
        return None

def generate_all_rRNA_reports(bam_file, window_size=25, output_dir="rRNA_reports"):
    """Generate individual reports for all rRNA elements in BAM file"""
    
    # Get all reference names
    result = subprocess.run(['samtools', 'idxstats', bam_file], 
                          capture_output=True, text=True)
    all_refs = [line.split('\t')[0] for line in result.stdout.strip().split('\n') 
                if line and not line.startswith('*')]
    
    # Filter for likely rRNA elements (exclude major chromosomes)
    major_chromosomes = {'2L', '2R', '3L', '3R', '4', 'X', 'Y', 'MT', 'mitochondrion_genome'}
    scaffold_keywords = ['Scaffold', 'scaffold', 'contig', 'random', 'chrUn', 'hap']
    
    rRNA_refs = []
    for ref in all_refs:
        # Skip major chromosomes
        if ref in major_chromosomes:
            continue
        # Skip scaffolds unless they might be rRNA
        if any(keyword in ref for keyword in scaffold_keywords):
            continue
        # Keep rRNA-related references
        if any(keyword in ref.lower() for keyword in ['rrna', 'rdna', 'ribosom', '18s', '28s', '5.8s', '5s', '16s', '12s', '23s']):
            rRNA_refs.append(ref)
        # Keep Wolbachia references
        elif ref.startswith('NZ_') or 'wmel' in ref.lower() or 'wolbachia' in ref.lower():
            rRNA_refs.append(ref)
        # Keep references that look like rRNA gene IDs
        elif ref.startswith('211000022') or 'rDNA' in ref:
            rRNA_refs.append(ref)
    
    print(f"Found {len(all_refs)} total references")
    print(f"Filtering to {len(rRNA_refs)} potential rRNA elements: {rRNA_refs}")
    
    if not rRNA_refs:
        print("No rRNA elements found! You may need to specify them manually.")
        print("All references found:")
        for ref in all_refs:
            print(f"  {ref}")
        return {}
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Generate report for each element
    all_results = {}
    
    for ref in rRNA_refs:
        result = create_individual_rRNA_report(bam_file, ref, window_size, output_dir)
        if result:
            all_results[ref] = result
    
    # Create summary comparison if we have results
    if all_results:
        create_summary_comparison(all_results, output_dir)
    
    print(f"\nAll reports completed! Check the '{output_dir}' directory")
    print(f"Individual reports: {list(all_results.keys())}")
    
    return all_results

def create_summary_comparison(all_results, output_dir):
    """Create a summary comparison across all rRNA elements"""
    
    if not all_results:
        return
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 10))
    
    # Prepare data
    names = list(all_results.keys())
    lengths = [all_results[name]['summary_stats']['length'] for name in names]
    mean_coverages = [all_results[name]['summary_stats']['mean_coverage'] for name in names]
    mean_ics = [all_results[name]['summary_stats']['mean_ic'] for name in names]
    conservation_rates = [all_results[name]['summary_stats']['conservation_rate'] for name in names]
    
    colors = plt.cm.Set3(np.linspace(0, 1, len(names)))
    
    # 1. Sequence lengths
    bars1 = ax1.bar(range(len(names)), lengths, color=colors)
    ax1.set_xlabel('rRNA Element')
    ax1.set_ylabel('Length (bp)')
    ax1.set_title('rRNA Element Lengths')
    ax1.set_xticks(range(len(names)))
    ax1.set_xticklabels(names, rotation=45, ha='right')
    
    # Add value labels on bars
    for i, (bar, length) in enumerate(zip(bars1, lengths)):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(lengths)*0.01,
                f'{length:,}', ha='center', va='bottom', fontsize=9)
    
    # 2. Mean coverage
    bars2 = ax2.bar(range(len(names)), mean_coverages, color=colors)
    ax2.set_xlabel('rRNA Element')
    ax2.set_ylabel('Mean Coverage')
    ax2.set_title('Mean Coverage Depth')
    ax2.set_xticks(range(len(names)))
    ax2.set_xticklabels(names, rotation=45, ha='right')
    
    # 3. Mean information content
    bars3 = ax3.bar(range(len(names)), mean_ics, color=colors)
    ax3.set_xlabel('rRNA Element')
    ax3.set_ylabel('Mean Information Content (bits)')
    ax3.set_title('Mean Information Content')
    ax3.set_xticks(range(len(names)))
    ax3.set_xticklabels(names, rotation=45, ha='right')
    ax3.set_ylim(0, 2)
    
    # 4. Conservation rates
    bars4 = ax4.bar(range(len(names)), conservation_rates, color=colors)
    ax4.set_xlabel('rRNA Element')
    ax4.set_ylabel('Conservation Rate (%)')
    ax4.set_title('Highly Conserved Positions (IC > 1.5)')
    ax4.set_xticks(range(len(names)))
    ax4.set_xticklabels(names, rotation=45, ha='right')
    ax4.set_ylim(0, 100)
    
    plt.suptitle('rRNA Elements Comparison Summary', fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    # Save summary
    summary_file = os.path.join(output_dir, 'rRNA_summary_comparison.png')
    plt.savefig(summary_file, dpi=300, bbox_inches='tight')
    plt.close(fig)  # Close figure to save memory
    
    print(f"Summary comparison saved: {summary_file}")

def main():
    """Main function for command line execution"""
    
    # Parse command line arguments
    if len(sys.argv) < 2:
        print("Usage: python individual_rRNA_reports.py <bam_file> [window_size] [output_dir] [specific_ref]")
        print("\nExample:")
        print("  python individual_rRNA_reports.py alignment.bam 25")
        print("  python individual_rRNA_reports.py alignment.bam 50 my_output/")
        print("  python individual_rRNA_reports.py alignment.bam 25 output/ rDNA")
        print("\nArguments:")
        print("  bam_file     : Path to BAM alignment file (required)")
        print("  window_size  : Window size for A content analysis (default: 25)")
        print("  output_dir   : Output directory for reports (default: rRNA_reports)")
        print("  specific_ref : Analyze only this reference (optional)")
        sys.exit(1)
    
    # Get arguments
    bam_file = sys.argv[1]
    window_size = int(sys.argv[2]) if len(sys.argv) > 2 else 25
    output_dir = sys.argv[3] if len(sys.argv) > 3 else "rRNA_reports"
    specific_ref = sys.argv[4] if len(sys.argv) > 4 else None
    
    # Check if BAM file exists
    if not os.path.exists(bam_file):
        print(f"Error: BAM file '{bam_file}' not found!")
        sys.exit(1)
    
    # Check if BAM index exists
    if not os.path.exists(bam_file + ".bai") and not os.path.exists(bam_file.replace(".bam", ".bai")):
        print(f"Warning: BAM index not found. Attempting to create index...")
        try:
            # Try to index the BAM file
            subprocess.run(['samtools', 'index', bam_file], check=True)
        except subprocess.CalledProcessError:
            print(f"Error: Cannot index BAM file directly. This often happens with STAR output.")
            print(f"Attempting to fix and re-index...")
            
            # Create a fixed BAM file
            fixed_bam = bam_file.replace('.bam', '_fixed.bam')
            try:
                # Sort and fix the BAM file
                print("Sorting and fixing BAM file...")
                subprocess.run(['samtools', 'sort', '-o', fixed_bam, bam_file], check=True)
                subprocess.run(['samtools', 'index', fixed_bam], check=True)
                bam_file = fixed_bam  # Use the fixed BAM file
                print(f"Successfully created fixed BAM: {fixed_bam}")
            except subprocess.CalledProcessError as e:
                print(f"Error: Could not fix BAM file: {e}")
                print("Please try:")
                print(f"  samtools sort -o {bam_file.replace('.bam', '_sorted.bam')} {bam_file}")
                print(f"  samtools index {bam_file.replace('.bam', '_sorted.bam')}")
                sys.exit(1)
    
    print(f"Starting rRNA analysis...")
    print(f"BAM file: {bam_file}")
    print(f"Window size: {window_size}")
    print(f"Output directory: {output_dir}")
    if specific_ref:
        print(f"Analyzing specific reference: {specific_ref}")
    print("-" * 50)
    
    # Run the analysis
    try:
        if specific_ref:
            # Analyze just one reference
            os.makedirs(output_dir, exist_ok=True)
            result = create_individual_rRNA_report(bam_file, specific_ref, window_size, output_dir)
            if result:
                print(f"\nAnalysis complete for {specific_ref}!")
                print(f"Report saved in: {output_dir}/")
            else:
                print(f"Failed to analyze {specific_ref}")
        else:
            # Analyze all rRNA elements
            results = generate_all_rRNA_reports(bam_file, window_size, output_dir)
            
            print("\n" + "=" * 50)
            print("ANALYSIS COMPLETE!")
            print(f"Reports generated for {len(results)} rRNA elements")
            print(f"Output files saved in: {output_dir}/")
            print("=" * 50)
        
    except Exception as e:
        print(f"Error during analysis: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
    