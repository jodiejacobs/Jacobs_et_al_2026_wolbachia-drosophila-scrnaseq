#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from Bio import AlignIO
from Bio.Align import AlignInfo
from collections import Counter
import sys
import os

def calculate_information_content_from_alignment(alignment_file):
    """Calculate information content per position from multiple sequence alignment"""
    
    # Read the alignment
    alignment = AlignIO.read(alignment_file, "fasta")
    align_length = alignment.get_alignment_length()
    num_sequences = len(alignment)
    
    print(f"Loaded alignment with {num_sequences} sequences and {align_length} positions")
    
    # Calculate information content per position
    ic_scores = []
    base_frequencies = []
    coverage_per_position = []
    
    for pos in range(align_length):
        # Get column at this position
        column = alignment[:, pos]
        
        # Count bases (excluding gaps)
        bases = [base.upper() for base in column if base != '-' and base != 'N']
        coverage_per_position.append(len(bases))
        
        if len(bases) == 0:
            ic_scores.append(0)
            base_frequencies.append({'A': 0, 'T': 0, 'G': 0, 'C': 0, 'total': 0})
            continue
        
        # Count each nucleotide
        base_counts = Counter(bases)
        total_bases = len(bases)
        
        # Calculate frequencies for A, T, G, C
        frequencies = []
        freq_dict = {'A': 0, 'T': 0, 'G': 0, 'C': 0, 'total': total_bases}
        
        for base in ['A', 'T', 'G', 'C']:
            count = base_counts.get(base, 0)
            freq = count / total_bases if total_bases > 0 else 0
            freq_dict[base] = freq
            if freq > 0:
                frequencies.append(freq)
        
        base_frequencies.append(freq_dict)
        
        # Calculate information content: IC = 2 + sum(p_i * log2(p_i))
        # Maximum IC is 2 bits (when only one nucleotide is present)
        if len(frequencies) > 0:
            ic = 2 + sum(f * np.log2(f + 1e-10) for f in frequencies)
            ic_scores.append(max(0, ic))
        else:
            ic_scores.append(0)
    
    return ic_scores, base_frequencies, coverage_per_position

def generate_consensus_sequence(alignment_file, threshold=0.5):
    """Generate consensus sequence from alignment"""
    
    alignment = AlignIO.read(alignment_file, "fasta")
    align_length = alignment.get_alignment_length()
    
    consensus = []
    
    for pos in range(align_length):
        column = alignment[:, pos]
        
        # Count bases (excluding gaps)
        bases = [base.upper() for base in column if base != '-' and base != 'N']
        
        if len(bases) == 0:
            consensus.append('N')
            continue
        
        # Get most common base
        base_counts = Counter(bases)
        most_common_base, count = base_counts.most_common(1)[0]
        
        # Use most common base if it's above threshold, otherwise N
        if count / len(bases) >= threshold:
            consensus.append(most_common_base)
        else:
            consensus.append('N')
    
    return ''.join(consensus)

def calculate_a_content_windows(base_frequencies, window_size=50):
    """Calculate A content in sliding windows"""
    
    a_content = []
    positions = []
    
    for i in range(0, len(base_frequencies) - window_size + 1, window_size//2):
        window_end = min(i + window_size, len(base_frequencies))
        
        # Sum A content and total bases in window
        total_a = sum(base_frequencies[j]['A'] * base_frequencies[j]['total'] 
                     for j in range(i, window_end))
        total_bases = sum(base_frequencies[j]['total'] for j in range(i, window_end))
        
        if total_bases > 0:
            a_fraction = total_a / total_bases
            a_content.append(a_fraction * 100)
            positions.append(i + window_size//2)
    
    return positions, a_content

def plot_alignment_analysis(alignment_file, window_size=50, output_dir="alignment_analysis"):
    """Create comprehensive analysis plots from alignment"""
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"Analyzing alignment: {alignment_file}")
    
    # Calculate information content
    ic_scores, base_frequencies, coverage = calculate_information_content_from_alignment(alignment_file)
    
    # Generate consensus
    consensus = generate_consensus_sequence(alignment_file)
    
    # Calculate A content windows
    a_positions, a_content = calculate_a_content_windows(base_frequencies, window_size)
    
    # Prepare nucleotide composition data
    positions = list(range(len(base_frequencies)))
    a_freq = [base_frequencies[i]['A'] * 100 for i in positions]
    t_freq = [base_frequencies[i]['T'] * 100 for i in positions]
    g_freq = [base_frequencies[i]['G'] * 100 for i in positions]
    c_freq = [base_frequencies[i]['C'] * 100 for i in positions]
    
    # Create comprehensive figure
    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(4, 2, height_ratios=[1, 1, 1, 1], hspace=0.3, wspace=0.3)
    
    # 1. Information content
    ax1 = fig.add_subplot(gs[0, :])
    ax1.plot(positions, ic_scores, linewidth=1, color='green')
    ax1.set_ylabel('Information Content (bits)')
    ax1.set_title('Information Content per Position (from Alignment)')
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(0, 2.1)
    
    # 2. Coverage (number of sequences contributing to each position)
    ax2 = fig.add_subplot(gs[1, :])
    ax2.plot(positions, coverage, linewidth=1, color='blue')
    ax2.fill_between(positions, coverage, alpha=0.3, color='blue')
    ax2.set_ylabel('Number of Sequences')
    ax2.set_title('Sequence Coverage per Position')
    ax2.grid(True, alpha=0.3)
    
    # 3. A content windows
    ax3 = fig.add_subplot(gs[2, :])
    ax3.plot(a_positions, a_content, linewidth=2, marker='o', markersize=2, color='red')
    ax3.axhline(y=25, color='gray', linestyle='--', alpha=0.5, label='Expected (25%)')
    ax3.set_ylabel('A Content (%)')
    ax3.set_title(f'A Content per {window_size}bp Window')
    ax3.grid(True, alpha=0.3)
    ax3.legend()
    
    # 4. Nucleotide composition
    ax4 = fig.add_subplot(gs[3, 0])
    ax4.fill_between(positions, 0, a_freq, alpha=0.7, color='red', label='A')
    ax4.fill_between(positions, a_freq, [a_freq[i] + t_freq[i] for i in range(len(positions))], 
                     alpha=0.7, color='blue', label='T')
    ax4.fill_between(positions, [a_freq[i] + t_freq[i] for i in range(len(positions))], 
                     [a_freq[i] + t_freq[i] + g_freq[i] for i in range(len(positions))], 
                     alpha=0.7, color='green', label='G')
    ax4.fill_between(positions, [a_freq[i] + t_freq[i] + g_freq[i] for i in range(len(positions))], 
                     [a_freq[i] + t_freq[i] + g_freq[i] + c_freq[i] for i in range(len(positions))], 
                     alpha=0.7, color='orange', label='C')
    ax4.set_xlabel('Alignment Position')
    ax4.set_ylabel('Nucleotide Frequency (%)')
    ax4.set_title('Base Composition per Position')
    ax4.legend()
    ax4.set_ylim(0, 100)
    
    # 5. Summary statistics
    ax5 = fig.add_subplot(gs[3, 1])
    ax5.axis('off')
    
    # Calculate summary stats
    alignment = AlignIO.read(alignment_file, "fasta")
    num_sequences = len(alignment)
    align_length = len(ic_scores)
    mean_ic = np.mean(ic_scores)
    mean_coverage = np.mean(coverage)
    high_ic_positions = sum(1 for ic in ic_scores if ic > 1.5)
    mean_a_content = np.mean(a_content) if a_content else 0
    
    # Overall composition
    total_a = sum(base_frequencies[i]['A'] * base_frequencies[i]['total'] for i in range(len(base_frequencies)))
    total_t = sum(base_frequencies[i]['T'] * base_frequencies[i]['total'] for i in range(len(base_frequencies)))
    total_g = sum(base_frequencies[i]['G'] * base_frequencies[i]['total'] for i in range(len(base_frequencies)))
    total_c = sum(base_frequencies[i]['C'] * base_frequencies[i]['total'] for i in range(len(base_frequencies)))
    total_bases = total_a + total_t + total_g + total_c
    
    stats_text = f"""
    ALIGNMENT SUMMARY
    
    Sequences: {num_sequences}
    Alignment length: {align_length:,} bp
    
    Conservation:
    • Mean IC: {mean_ic:.3f} bits
    • Highly conserved (IC>1.5): {high_ic_positions}
    • Conservation rate: {high_ic_positions/align_length*100:.1f}%
    
    Coverage:
    • Mean sequences/position: {mean_coverage:.1f}
    • Min coverage: {min(coverage)}
    • Max coverage: {max(coverage)}
    
    Overall Composition:
    • A: {total_a/total_bases*100:.1f}%
    • T: {total_t/total_bases*100:.1f}%
    • G: {total_g/total_bases*100:.1f}%
    • C: {total_c/total_bases*100:.1f}%
    
    Windows:
    • Mean A content: {mean_a_content:.1f}%
    """
    
    ax5.text(0.05, 0.95, stats_text, transform=ax5.transAxes, fontsize=10,
            verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray", alpha=0.8))
    
    plt.suptitle(f'rRNA Alignment Analysis - {os.path.basename(alignment_file)}', 
                 fontsize=16, fontweight='bold')
    
    # Save main figure
    output_file = os.path.join(output_dir, 'rRNA_alignment_analysis.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()
    
    # Save consensus sequence
    consensus_file = os.path.join(output_dir, 'consensus_sequence.fasta')
    with open(consensus_file, 'w') as f:
        f.write(f">Consensus_rRNA\n{consensus}\n")
    
    print(f"Analysis complete!")
    print(f"Main plot saved: {output_file}")
    print(f"Consensus sequence saved: {consensus_file}")
    print(f"Consensus length: {len(consensus)} bp")
    
    return {
        'ic_scores': ic_scores,
        'base_frequencies': base_frequencies,
        'coverage': coverage,
        'consensus': consensus,
        'a_content': a_content,
        'summary_stats': {
            'num_sequences': num_sequences,
            'align_length': align_length,
            'mean_ic': mean_ic,
            'conservation_rate': high_ic_positions/align_length*100,
            'mean_coverage': mean_coverage
        }
    }

def plot_conservation_heatmap(alignment_file, output_dir="alignment_analysis"):
    """Create a heatmap showing conservation across sequences and positions"""
    
    alignment = AlignIO.read(alignment_file, "fasta")
    align_length = alignment.get_alignment_length()
    num_sequences = len(alignment)
    
    # Create matrix: sequences x positions
    seq_matrix = np.zeros((num_sequences, align_length))
    
    # Encode bases: A=1, T=2, G=3, C=4, gap/N=0
    base_encoding = {'A': 1, 'T': 2, 'G': 3, 'C': 4, '-': 0, 'N': 0}
    
    for i, seq_record in enumerate(alignment):
        for j, base in enumerate(str(seq_record.seq)):
            seq_matrix[i, j] = base_encoding.get(base.upper(), 0)
    
    # Create heatmap
    plt.figure(figsize=(20, max(8, num_sequences * 0.3)))
    
    # Custom colormap
    colors = ['white', 'red', 'blue', 'green', 'orange']  # gap, A, T, G, C
    from matplotlib.colors import ListedColormap
    cmap = ListedColormap(colors)
    
    plt.imshow(seq_matrix, cmap=cmap, aspect='auto', interpolation='nearest')
    plt.colorbar(label='Base (A=1, T=2, G=3, C=4, gap=0)', ticks=[0, 1, 2, 3, 4])
    
    plt.xlabel('Alignment Position')
    plt.ylabel('Sequence')
    plt.title('rRNA Alignment Heatmap')
    
    # Add sequence names if not too many
    if num_sequences <= 50:
        seq_names = [seq.id[:20] for seq in alignment]  # Truncate long names
        plt.yticks(range(num_sequences), seq_names)
    
    plt.tight_layout()
    
    heatmap_file = os.path.join(output_dir, 'alignment_heatmap.png')
    plt.savefig(heatmap_file, dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"Heatmap saved: {heatmap_file}")

def main():
    """Main function for command line execution"""
    
    if len(sys.argv) < 2:
        print("Usage: python alignment_information_content.py <alignment_file> [window_size] [output_dir]")
        print("\nExample:")
        print("  python alignment_information_content.py rRNA_aligned.fasta 50")
        print("  python alignment_information_content.py rRNA_aligned.fasta 25 my_output/")
        print("\nArguments:")
        print("  alignment_file : MAFFT alignment file in FASTA format (required)")
        print("  window_size    : Window size for A content analysis (default: 50)")
        print("  output_dir     : Output directory (default: alignment_analysis)")
        sys.exit(1)
    
    # Get arguments
    alignment_file = sys.argv[1]
    window_size = int(sys.argv[2]) if len(sys.argv) > 2 else 50
    output_dir = sys.argv[3] if len(sys.argv) > 3 else "alignment_analysis"
    
    # Check if alignment file exists
    if not os.path.exists(alignment_file):
        print(f"Error: Alignment file '{alignment_file}' not found!")
        sys.exit(1)
    
    print(f"Starting alignment analysis...")
    print(f"Alignment file: {alignment_file}")
    print(f"Window size: {window_size}")
    print(f"Output directory: {output_dir}")
    print("-" * 50)
    
    try:
        # Run main analysis
        results = plot_alignment_analysis(alignment_file, window_size, output_dir)
        
        # Create heatmap
        plot_conservation_heatmap(alignment_file, output_dir)
        
        print("\n" + "=" * 50)
        print("ANALYSIS COMPLETE!")
        print(f"Files saved in: {output_dir}/")
        print("=" * 50)
        
    except Exception as e:
        print(f"Error during analysis: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
