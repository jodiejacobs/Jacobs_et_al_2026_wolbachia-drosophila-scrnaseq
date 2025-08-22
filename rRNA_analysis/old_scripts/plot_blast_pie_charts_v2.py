#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import argparse
from pathlib import Path
import sys
import random
import seaborn as sns
from concurrent.futures import ThreadPoolExecutor
import multiprocessing as mp
from threading import Lock

# Set matplotlib for better compatibility and performance
matplotlib.use('Agg')  # Use non-interactive backend for speed
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['font.family'] = 'sans-serif'

# Global variables for threading
colors_file = '/private/groups/russelllab/jodie/scRNAseq/wolbachia_titer/rRNA_analysis/scripts/genus_colors.csv'
genus_colors = {}
palette_hex = []
color_lock = Lock()  # Thread-safe color assignment

def initialize_colors():
    """Initialize color palette once"""
    global genus_colors, palette_hex
    try:
        genus_colors = pd.read_csv(colors_file, index_col=0)['color'].to_dict()
    except FileNotFoundError:
        genus_colors = {}
    
    palette = sns.color_palette("husl", 80)
    palette_hex = [f"#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}" for r, g, b in palette]
    random.shuffle(palette_hex)

def get_color_for_genus(genus):
    """Thread-safe color assignment"""
    global genus_colors, palette_hex
    
    if genus.startswith('Other ('):
        return '#cccccc'
    
    with color_lock:
        if genus in genus_colors:
            return genus_colors[genus]
        else:
            # Find unused color
            used_colors = set(genus_colors.values())
            for color in palette_hex:
                if color not in used_colors:
                    genus_colors[genus] = color
                    return color
            # Fallback if all colors used
            color = palette_hex[len(genus_colors) % len(palette_hex)]
            genus_colors[genus] = color
            return color

def load_summary_file(filepath):
    """Optimized loading using pandas read_csv with regex"""
    try:
        # Read file as text first for preprocessing
        with open(filepath, 'r') as f:
            lines = [line.strip() for line in f if line.strip() and not line.startswith('#')]
        
        if not lines:
            return None
            
        # Parse lines more efficiently
        data = []
        for line in lines:
            parts = line.split(None, 1)
            if len(parts) >= 2:
                try:
                    count = int(parts[0])
                    rest = parts[1]
                    
                    # Extract percent identity more efficiently
                    if '(' in rest and '% identity)' in rest:
                        last_paren = rest.rfind('(')
                        organism = rest[:last_paren].strip()
                        identity_part = rest[last_paren+1:rest.rfind('% identity)')]
                        identity = float(identity_part)
                    else:
                        organism = rest
                        identity = None
                    
                    data.append((count, organism, identity))
                except (ValueError, IndexError):
                    continue
        
        if data:
            # Create DataFrame directly from list of tuples (faster)
            df = pd.DataFrame(data, columns=['count', 'organism', 'identity'])
            return df
        else:
            return None
            
    except Exception as e:
        print(f"Error loading {filepath}: {e}")
        return None

def extract_genus_vectorized(names):
    """Vectorized genus extraction using pandas string operations"""
    return names.str.split().str[0]

def process_single_file(args_tuple):
    """Process a single file - designed for threading"""
    summary_file, output_dir = args_tuple
    
    sample_name = Path(summary_file).stem.replace('.blast.summary', '').replace('.blast', '')
    print(f"Processing {sample_name}...")
    
    df = load_summary_file(summary_file)
    if df is None:
        print(f"  Skipping {summary_file} - no valid data")
        return None
    
    # Filter high-confidence hits efficiently
    if 'identity' in df.columns and df['identity'].notna().any():
        high_conf_mask = df['identity'] >= 95.0
        high_conf_df = df[high_conf_mask].copy()
        if len(high_conf_df) == 0:
            print(f"Warning: No hits ≥95% identity for {sample_name}, using all data")
            high_conf_df = df.copy()
        else:
            print(f"Filtered {sample_name}: {len(high_conf_df)}/{len(df)} hits ≥95% identity")
    else:
        high_conf_df = df.copy()
    
    # Vectorized genus extraction
    high_conf_df['genus'] = extract_genus_vectorized(high_conf_df['organism'])
    
    # Efficient groupby operations
    if 'identity' in high_conf_df.columns and high_conf_df['identity'].notna().any():
        genus_stats = high_conf_df.groupby('genus').agg({
            'count': 'sum',
            'identity': ['mean', 'min', 'max']
        }).round(1)
        genus_stats.columns = ['count', 'avg_identity', 'min_identity', 'max_identity']
        genus_stats = genus_stats.reset_index()
    else:
        genus_stats = high_conf_df.groupby('genus')['count'].sum().reset_index()
        genus_stats[['avg_identity', 'min_identity', 'max_identity']] = None
    
    genus_stats = genus_stats.sort_values('count', ascending=False)
    
    # Calculate percentages
    total_reads = genus_stats['count'].sum()
    genus_stats['percentage'] = (genus_stats['count'] / total_reads) * 100
    
    # Group small genera efficiently
    threshold_pct = 2.0
    threshold_count = total_reads * (threshold_pct / 100)
    major_mask = genus_stats['count'] >= threshold_count
    df_major = genus_stats[major_mask].copy()
    df_minor = genus_stats[~major_mask]
    
    if len(df_minor) > 0:
        other_count = df_minor['count'].sum()
        other_pct = (other_count / total_reads) * 100
        
        if 'avg_identity' in df_minor.columns and df_minor['avg_identity'].notna().any():
            other_stats = {
                'genus': f'Other ({len(df_minor)} genera)',
                'count': other_count,
                'percentage': other_pct,
                'avg_identity': df_minor['avg_identity'].mean(),
                'min_identity': df_minor['min_identity'].min(),
                'max_identity': df_minor['max_identity'].max()
            }
        else:
            other_stats = {
                'genus': f'Other ({len(df_minor)} genera)',
                'count': other_count,
                'percentage': other_pct,
                'avg_identity': None,
                'min_identity': None,
                'max_identity': None
            }
        
        df_major = pd.concat([df_major, pd.DataFrame([other_stats])], ignore_index=True)
    
    df_major = df_major.sort_values('count', ascending=False)
    
    # Generate colors efficiently using thread-safe function
    colors = [get_color_for_genus(genus) for genus in df_major['genus']]
    
    # Create plots with optimized settings
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8))
    
    # Pie chart
    wedges, texts, autotexts = ax1.pie(df_major['count'], 
                                      labels=None,
                                      colors=colors,
                                      autopct='%1.1f%%',
                                      startangle=90)
    
    ax1.set_title(f'{sample_name}\n16S rRNA Genus Distribution (≥95% identity)\n({total_reads:,} total reads, ≥2% shown)', 
                  fontsize=14, fontweight='bold')
    
    # Optimize text rendering
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontweight('bold')
        autotext.set_fontsize(10)
    
    # Create legend labels efficiently
    legend_labels = []
    for _, row in df_major.iterrows():
        if pd.notna(row.get('avg_identity')):
            legend_labels.append(f"{row['genus']} ({row['count']:,}, {row['percentage']:.1f}%, {row['avg_identity']:.1f}% ID)")
        else:
            legend_labels.append(f"{row['genus']} ({row['count']:,}, {row['percentage']:.1f}%)")
    
    ax1.legend(wedges, legend_labels, 
              title="Genus (reads, %, avg identity)",
              loc="center left", 
              bbox_to_anchor=(1, 0, 0.5, 1),
              fontsize=9)
    
    # Bar chart
    y_pos = np.arange(len(df_major))
    bars = ax2.barh(y_pos, df_major['count'], color=colors)
    ax2.set_yticks(y_pos)
    
    # Y-axis labels
    y_labels = []
    for _, row in df_major.iterrows():
        if pd.notna(row.get('avg_identity')):
            y_labels.append(f"{row['genus']}\n({row['avg_identity']:.1f}% ID)")
        else:
            y_labels.append(row['genus'])
    ax2.set_yticklabels(y_labels, fontsize=9)
    
    ax2.set_xlabel('Number of Reads')
    ax2.set_title(f'{sample_name}\nRead Counts by Genus (≥95% identity, ≥2%)\n(with average % identity)', 
                  fontsize=12, fontweight='bold')
    ax2.grid(axis='x', alpha=0.3)
    
    # Add count labels efficiently
    for bar, count in zip(bars, df_major['count']):
        ax2.text(bar.get_width() + total_reads*0.01, bar.get_y() + bar.get_height()/2, 
                f'{count:,}', ha='left', va='center', fontsize=9)
    
    plt.tight_layout()
    
    # Save only SVG (faster than PNG)
    safe_sample_name = sample_name.replace(':', '_').replace('-', '_')
    plt.savefig(f"{output_dir}/{safe_sample_name}_16S_genus_pie_chart_95pct.svg", bbox_inches='tight')
    plt.close()
    
    return {
        'sample': sample_name,
        'total_reads': total_reads,
        'num_genera': len(genus_stats),
        'top_genus': genus_stats.iloc[0]['genus'],
        'top_genus_pct': genus_stats.iloc[0]['percentage'],
        'top_genus_identity': genus_stats.iloc[0].get('avg_identity')
    }

def main():
    parser = argparse.ArgumentParser(description='Create pie charts from BLAST summary files')
    parser.add_argument('summary_files', nargs='+', help='BLAST summary files (*.blast.summary)')
    parser.add_argument('--output-dir', default='blast_pie_charts', help='Output directory')
    parser.add_argument('--threads', type=int, default=mp.cpu_count()//2, 
                       help='Number of parallel threads (default: half of CPU cores)')
    
    args = parser.parse_args()
    
    # Initialize colors once
    initialize_colors()
    
    # Create output directory
    Path(args.output_dir).mkdir(exist_ok=True)
    
    # Prepare arguments for multiprocessing
    process_args = [(summary_file, args.output_dir) for summary_file in args.summary_files]
    
    # Process files in parallel using threads (shares memory)
    print(f"Processing {len(args.summary_files)} files using {args.threads} threads...")
    
    all_stats = []
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        results = executor.map(process_single_file, process_args)
        all_stats = [result for result in results if result is not None]
    
    # Save updated colors
    if genus_colors:
        df_colors = pd.DataFrame(list(genus_colors.items()), columns=['genus', 'color'])
        df_colors.to_csv(colors_file, index=False)
    
    # Create summary table
    if all_stats:
        summary_df = pd.DataFrame(all_stats)
        summary_df.to_csv(f"{args.output_dir}/blast_summary_stats.csv", index=False)
        print(f"\nSummary statistics saved to {args.output_dir}/blast_summary_stats.csv")
        print("\nOverview:")
        for _, row in summary_df.iterrows():
            if pd.notna(row.get('top_genus_identity')):
                print(f"  {row['sample']}: {row['total_reads']:,} reads, {row['num_genera']} genera, "
                      f"top: {row['top_genus']} ({row['top_genus_pct']:.1f}%, {row['top_genus_identity']:.1f}% ID)")
            else:
                print(f"  {row['sample']}: {row['total_reads']:,} reads, {row['num_genera']} genera, "
                      f"top: {row['top_genus']} ({row['top_genus_pct']:.1f}%)")
    
    print(f"\nAll pie charts saved to {args.output_dir}/")

if __name__ == "__main__":
    main()
