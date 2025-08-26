#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
from pathlib import Path
import sys

def load_coverage_file(filepath):
    """Load coverage file and return DataFrame"""
    try:
        df = pd.read_csv(filepath, sep='\t', names=['region', 'position', 'coverage'])
        return df
    except Exception as e:
        print(f"Error loading {filepath}: {e}")
        return None

def get_sample_color(sample_name, color_dict):
    """Get color for sample based on pattern matching"""
    for pattern, color in color_dict.items():
        if pattern in sample_name:
            return color
    return 'gray'  # default color

def plot_coverage_by_region(coverage_files, output_dir="plots"):
    """Create line plots for each rRNA region"""
    
    # Define color dictionary for cell types
    color_dict = {
        'JW18DOX-SV': '#F3CA40', # Saffron
        'JW18DOX-Ctrl': '#3D2B56', # Giants orange
        'JW18wMel-SV': '#35ABAF', # Verdris
        'JW18wMel-Ctrl': '#2E294E', # Space cadet
        'mei-P26-uninf': '#E28413', # Light Orange
        'mei-P26-wMel': '#C2470A', # Deep Orange  
        'OreR-uninf': '#2C497F', # Medium Blue
        'OreR-wMel': '#71A1C9', # Dusty Blue      
    }
    
    # Create output directory
    Path(output_dir).mkdir(exist_ok=True)
    
    # Collect all data
    all_data = []
    
    for file_path in coverage_files:
        sample_name = Path(file_path).stem.replace('_coverage', '')
        df = load_coverage_file(file_path)
        
        if df is not None:
            df['sample'] = sample_name
            all_data.append(df)
        else:
            print(f"Skipping {file_path}")
    
    if not all_data:
        print("No valid coverage files found!")
        return
    
    # Combine all data
    combined_df = pd.concat(all_data, ignore_index=True)
    
    # Get unique regions
    regions = combined_df['region'].unique()
    print(f"Found {len(regions)} rRNA regions")
    # Create plots for each region
    for region in regions:
        region_data = combined_df[combined_df['region'] == region]
        
        # Create figure
        plt.figure(figsize=(12, 6))
        
        # Plot each sample as a line
        samples = region_data['sample'].unique()
        
        for sample in samples:
            sample_data = region_data[region_data['sample'] == sample]
            color = get_sample_color(sample, color_dict)
            plt.plot(sample_data['position'], sample_data['coverage'], 
                    label=sample, linewidth=1.5, color=color)
        
        plt.xlabel('Position')
        plt.ylabel('Coverage')
        plt.title(f'{region}')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.grid(True, alpha=0.3)
        
        # Save linear plot with vector text
        safe_region_name = region.replace(':', '_').replace('-', '_')
        plt.tight_layout()
        plt.savefig(f"{output_dir}/{safe_region_name}.png", dpi=300, bbox_inches='tight')
        plt.savefig(f"{output_dir}/{safe_region_name}.svg", bbox_inches='tight')
        plt.close()
        
        # Create log plot
        plt.figure(figsize=(12, 6))
        
        for sample in samples:
            sample_data = region_data[region_data['sample'] == sample]
            color = get_sample_color(sample, color_dict)
            # Add 1 to avoid log(0) issues
            plt.plot(sample_data['position'], sample_data['coverage'] + 1, 
                    label=sample, linewidth=1.5, color=color)
        
        plt.xlabel('Position')
        plt.ylabel('Coverage + 1 (log scale)')
        plt.yscale('log')
        plt.title(f'{region} (Log Scale)')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        # plt.savefig(f"{output_dir}/{safe_region_name}_log.png", dpi=300, bbox_inches='tight')
        plt.savefig(f"{output_dir}/{safe_region_name}_log.svg", bbox_inches='tight')
        plt.close()
        
        print(f"Saved: {safe_region_name}.png/.svg and {safe_region_name}_log.png/.svg")

def main():
    parser = argparse.ArgumentParser(description='Plot rRNA coverage lines for each region')
    parser.add_argument('file_list', help='Text file with list of coverage.tsv files')
    parser.add_argument('--output-dir', default='rRNA_plots', help='Output directory')
    
    args = parser.parse_args()
    
    # Read file list
    try:
        with open(args.file_list, 'r') as f:
            coverage_files = [line.strip() for line in f if line.strip()]
    except FileNotFoundError:
        print(f"Error: File list '{args.file_list}' not found")
        sys.exit(1)
    
    print(f"Processing {len(coverage_files)} coverage files")
    plot_coverage_by_region(coverage_files, args.output_dir)
    print("Done!")

if __name__ == "__main__":
    main()
