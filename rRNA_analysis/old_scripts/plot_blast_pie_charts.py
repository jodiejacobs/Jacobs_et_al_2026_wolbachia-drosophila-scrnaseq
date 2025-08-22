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


# for summary in data/pipseq/ovary/*.summary; do      echo "Processing $(basename $summary)...";     python scripts/plot_blast_pie_charts.py "$summary" --output-dir blast_plots; done

# Set matplotlib for better compatibility
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['font.family'] = 'sans-serif'


# Set up plot colors:
colors_file = '/private/groups/russelllab/jodie/scRNAseq/wolbachia_titer/rRNA_analysis/scripts/genus_colors.csv'
genus_colors = pd.read_csv(colors_file, index_col=0)['color'].to_dict()

# For generating new colors if needed
palette = sns.color_palette("husl", 80) 
palette_hex = [f"#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}" for r, g, b in palette]
random.shuffle(palette_hex)

def generate_seaborn_color():
    used_colors = set(genus_colors.values())
    
    while True:
        color = palette_hex[generate_seaborn_color.counter % len(palette_hex)]
        generate_seaborn_color.counter += 1
        
        if color not in used_colors:
            return color

if not hasattr(generate_seaborn_color, 'counter'):
    generate_seaborn_color.counter = len(genus_colors)  # Start from the next index

def load_summary_file(filepath):
    """Load BLAST summary file and return DataFrame"""
    try:
        # Read the summary file - format: "count organism_name (percent% identity)"
        data = []
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    # Parse: "569 Wolbachia pipientis strain wMel 16S ribosomal RNA, complete sequence (100.000% identity)"
                    parts = line.split(None, 1)  # Split on first whitespace only
                    if len(parts) >= 2:
                        count = int(parts[0])
                        rest = parts[1]
                        
                        # Extract percent identity from parentheses
                        if '(' in rest and '% identity)' in rest:
                            # Find the last occurrence of ( to handle organism names with parentheses
                            last_paren = rest.rfind('(')
                            organism = rest[:last_paren].strip()
                            identity_str = rest[last_paren+1:].replace('% identity)', '').strip()
                            try:
                                identity = float(identity_str)
                            except ValueError:
                                print(f"Warning: Could not parse identity '{identity_str}' for organism '{organism}'")
                                identity = None
                        else:
                            # Fallback if no identity info
                            organism = rest
                            identity = None
                        
                        data.append({
                            'count': count, 
                            'organism': organism,
                            'identity': identity
                        })
        
        if data:
            df = pd.DataFrame(data)
            return df
        else:
            print(f"No valid data found in {filepath}")
            return None
            
    except Exception as e:
        print(f"Error loading {filepath}: {e}")
        return None

def extract_genus(name):
    """Extract genus (first word) from organism name"""
    # Split the name and take the first word (genus)
    genus = name.split()[0]
    return genus

def create_pie_chart(df, sample_name, output_dir):
    """Create pie chart for organism distribution"""
    
    # Filter to only high-confidence hits (>=95% identity)
    if 'identity' in df.columns and df['identity'].notna().any():
        original_count = len(df)
        high_conf_df = df[df['identity'] >= 95.0].copy()
        if len(high_conf_df) == 0:
            print(f"Warning: No hits ≥95% identity for {sample_name}, using all data")
            high_conf_df = df.copy()
        else:
            print(f"Filtered {sample_name}: {len(high_conf_df)}/{original_count} hits ≥95% identity")
    else:
        print(f"No identity data for {sample_name}, using all hits")
        high_conf_df = df.copy()
    
    # Extract genus from organism names
    high_conf_df['genus'] = high_conf_df['organism'].apply(extract_genus)
    
    # Group by genus and calculate statistics
    if 'identity' in high_conf_df.columns and high_conf_df['identity'].notna().any():
        genus_stats = high_conf_df.groupby('genus').agg({
            'count': 'sum',
            'identity': ['mean', 'min', 'max']
        }).round(1)
        # Flatten column names
        genus_stats.columns = ['count', 'avg_identity', 'min_identity', 'max_identity']
        genus_stats = genus_stats.reset_index()
    else:
        genus_stats = high_conf_df.groupby('genus')['count'].sum().reset_index()
        genus_stats['avg_identity'] = None
        genus_stats['min_identity'] = None
        genus_stats['max_identity'] = None
    
    genus_stats = genus_stats.sort_values('count', ascending=False)
    
    # Calculate percentages
    total_reads = genus_stats['count'].sum()
    genus_stats['percentage'] = (genus_stats['count'] / total_reads) * 100
    
    # Group genera with <1% of total reads into "Other"
    threshold_pct = 2.0  # 1% threshold
    threshold_count = total_reads * (threshold_pct / 100)
    df_major = genus_stats[genus_stats['count'] >= threshold_count].copy()
    df_minor = genus_stats[genus_stats['count'] < threshold_count]
    
    if len(df_minor) > 0:
        other_count = df_minor['count'].sum()
        other_pct = (other_count / total_reads) * 100
        
        if 'avg_identity' in df_minor.columns and df_minor['avg_identity'].notna().any():
            other_avg_identity = df_minor['avg_identity'].mean()
            other_min_identity = df_minor['min_identity'].min()
            other_max_identity = df_minor['max_identity'].max()
        else:
            other_avg_identity = None
            other_min_identity = None
            other_max_identity = None
            
        other_row = pd.DataFrame({
            'genus': [f'Other ({len(df_minor)} genera)'],
            'count': [other_count],
            'percentage': [other_pct],
            'avg_identity': [other_avg_identity],
            'min_identity': [other_min_identity],
            'max_identity': [other_max_identity]
        })
        
        df_major = pd.concat([df_major, other_row], ignore_index=True)
    
    # Sort by count for consistent ordering
    df_major = df_major.sort_values('count', ascending=False)
    
    # Create figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8))
    
    # Define genus-to-color dictionary
    colors = []
    for genus in df_major['genus']:
        if genus.startswith('Other ('):
            # Special color for the "Other" category
            colors.append('#cccccc')  # light gray
        elif genus in genus_colors:
            colors.append(genus_colors[genus])
        else:
            # Generate random color and add to dictionary
            random_color = generate_seaborn_color()
            genus_colors[genus] = random_color
            colors.append(random_color)

    # Pie chart
    wedges, texts, autotexts = ax1.pie(df_major['count'], 
                                      labels=None,
                                      colors=colors,
                                      autopct='%1.1f%%',
                                      startangle=90)
    
    ax1.set_title(f'{sample_name}\n16S rRNA Genus Distribution (≥95% identity)\n({total_reads:,} total reads, ≥1% shown)', 
                  fontsize=14, fontweight='bold')
    
    # Make percentage text more readable
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontweight('bold')
        autotext.set_fontsize(10)
    
    # Legend with counts and percentages
    legend_labels = []
    for idx, row in df_major.iterrows():
        if pd.notna(row.get('avg_identity')):
            label = f"{row['genus']} ({row['count']:,}, {row['percentage']:.1f}%, {row['avg_identity']:.1f}% ID)"
        else:
            label = f"{row['genus']} ({row['count']:,}, {row['percentage']:.1f}%)"
        legend_labels.append(label)
    
    ax1.legend(wedges, legend_labels, 
              title="Genus (reads, %, avg identity)",
              loc="center left", 
              bbox_to_anchor=(1, 0, 0.5, 1),
              fontsize=9)
    
    # Bar chart
    bars = ax2.barh(range(len(df_major)), df_major['count'], color=colors)
    ax2.set_yticks(range(len(df_major)))
    
    # Y-axis labels with identity info
    y_labels = []
    for idx, row in df_major.iterrows():
        if pd.notna(row.get('avg_identity')):
            y_labels.append(f"{row['genus']}\n({row['avg_identity']:.1f}% ID)")
        else:
            y_labels.append(row['genus'])
    ax2.set_yticklabels(y_labels, fontsize=9)
    
    ax2.set_xlabel('Number of Reads')
    ax2.set_title(f'{sample_name}\nRead Counts by Genus (≥95% identity, ≥1%)\n(with average % identity)', 
                  fontsize=12, fontweight='bold')
    ax2.grid(axis='x', alpha=0.3)
    
    # Add count labels on bars
    for i, (bar, count) in enumerate(zip(bars, df_major['count'])):
        ax2.text(bar.get_width() + total_reads*0.01, bar.get_y() + bar.get_height()/2, 
                f'{count:,}', ha='left', va='center', fontsize=9)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save plots
    safe_sample_name = sample_name.replace(':', '_').replace('-', '_')
    # plt.savefig(f"{output_dir}/{safe_sample_name}_16S_genus_pie_chart_95pct.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_dir}/{safe_sample_name}_16S_genus_pie_chart_95pct.svg", bbox_inches='tight')
    plt.close()
    
    # Create summary statistics
    summary_stats = {
        'sample': sample_name,
        'total_reads': total_reads,
        'num_genera': len(genus_stats),
        'top_genus': genus_stats.iloc[0]['genus'],
        'top_genus_pct': genus_stats.iloc[0]['percentage'],
        'top_genus_identity': genus_stats.iloc[0].get('avg_identity')
    }
    
    # Convert dictionary back to DataFrame and save
    df_colors = pd.DataFrame(list(genus_colors.items()), columns=['genus', 'color'])
    df_colors.to_csv(colors_file, index=False)

    return summary_stats

def main():
    parser = argparse.ArgumentParser(description='Create pie charts from BLAST summary files')
    parser.add_argument('summary_files', nargs='+', help='BLAST summary files (*.blast.summary)')
    parser.add_argument('--output-dir', default='blast_pie_charts', help='Output directory')
    
    args = parser.parse_args()
    
    # Create output directory
    Path(args.output_dir).mkdir(exist_ok=True)
    
    all_stats = []
    
    for summary_file in args.summary_files:
        sample_name = Path(summary_file).stem.replace('.blast.summary', '').replace('.blast', '')
        print(f"Processing {sample_name}...")
        
        df = load_summary_file(summary_file)
        if df is not None:
            stats = create_pie_chart(df, sample_name, args.output_dir)
            all_stats.append(stats)
            print(f"  Created pie chart: {sample_name}_16S_genus_pie_chart_95pct.png/.svg")
        else:
            print(f"  Skipping {summary_file} - no valid data")
    
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