#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
from pathlib import Path
import sys

# Configure matplotlib for Adobe Illustrator compatibility
plt.rcParams['svg.fonttype'] = 'none'  # Keep fonts as text, not paths
plt.rcParams['font.size'] = 6  # 6pt font size
plt.rcParams['axes.titlesize'] = 6
plt.rcParams['axes.labelsize'] = 6
plt.rcParams['xtick.labelsize'] = 6
plt.rcParams['ytick.labelsize'] = 6
plt.rcParams['legend.fontsize'] = 6
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']

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
        'JW18DOX-Ctrl': '#F46036', # Giants orange
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
        
        # Create figure with enough space for legend (plot area will be 3x2)
        fig = plt.figure(figsize=(5, 2))
        ax = fig.add_axes([0.1, 0.15, 0.6, 0.75])  # [left, bottom, width, height] in figure coordinates
        
        # Plot each sample as a line
        samples = region_data['sample'].unique()
        
        for sample in samples:
            sample_data = region_data[region_data['sample'] == sample]
            color = get_sample_color(sample, color_dict)
            ax.plot(sample_data['position'], sample_data['coverage'], 
                   label=sample, linewidth=0.75, color=color)
        
        ax.set_xlabel('Position')
        ax.set_ylabel('Coverage')
        ax.set_title(f'{region}')
        
        # Position legend outside plot area
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)
        ax.grid(True, alpha=0.3, linewidth=0.5)
        
        # Reduce tick density for cleaner look
        ax.locator_params(axis='x', nbins=4)
        ax.locator_params(axis='y', nbins=4)
        
        # Save linear plot
        safe_region_name = region.replace(':', '_').replace('-', '_')
        plt.savefig(f"{output_dir}/{safe_region_name}.svg", 
                   bbox_inches='tight', pad_inches=0.05)
        plt.close()
        
        # Create log plot with same setup (plot area will be 3x2)
        fig = plt.figure(figsize=(5, 2))
        ax = fig.add_axes([0.1, 0.15, 0.6, 0.75])  # [left, bottom, width, height] in figure coordinates
        
        for sample in samples:
            sample_data = region_data[region_data['sample'] == sample]
            color = get_sample_color(sample, color_dict)
            # Add 1 to avoid log(0) issues
            ax.plot(sample_data['position'], sample_data['coverage'] + 1, 
                   label=sample, linewidth=0.75, color=color)
        
        ax.set_xlabel('Position')
        ax.set_ylabel('Coverage + 1 (log scale)')
        ax.set_yscale('log')
        ax.set_title(f'{region} (Log Scale)')
        
        # Position legend outside plot area
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)
        ax.grid(True, alpha=0.3, linewidth=0.5)
        
        # Reduce tick density for cleaner look
        ax.locator_params(axis='x', nbins=4)
        
        plt.savefig(f"{output_dir}/{safe_region_name}_log.svg", 
                   bbox_inches='tight', pad_inches=0.05)
        plt.close()
        
        print(f"Saved: {safe_region_name}.svg and {safe_region_name}_log.svg")

def main():
    parser = argparse.ArgumentParser(description='Plot rRNA coverage lines for each region (Illustrator-ready)')
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
    main()#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
from pathlib import Path
import sys

# Configure matplotlib for Adobe Illustrator compatibility
plt.rcParams['svg.fonttype'] = 'none'  # Keep fonts as text, not paths
plt.rcParams['font.size'] = 6  # 6pt font size
plt.rcParams['axes.titlesize'] = 6
plt.rcParams['axes.labelsize'] = 6
plt.rcParams['xtick.labelsize'] = 6
plt.rcParams['ytick.labelsize'] = 6
plt.rcParams['legend.fontsize'] = 6
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']

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
        'JW18DOX-Ctrl': '#F46036', # Giants orange
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
        
        # Create figure with enough space for legend (plot area will be 3x2)
        fig = plt.figure(figsize=(5, 2))
        ax = fig.add_axes([0.1, 0.15, 0.6, 0.75])  # [left, bottom, width, height] in figure coordinates
        
        # Plot each sample as a line
        samples = region_data['sample'].unique()
        
        for sample in samples:
            sample_data = region_data[region_data['sample'] == sample]
            color = get_sample_color(sample, color_dict)
            ax.plot(sample_data['position'], sample_data['coverage'], 
                   label=sample, linewidth=0.75, color=color)
        
        ax.set_xlabel('Position')
        ax.set_ylabel('Coverage')
        ax.set_title(f'{region}')
        
        # Position legend outside plot area
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)
        ax.grid(True, alpha=0.3, linewidth=0.5)
        
        # Reduce tick density for cleaner look
        ax.locator_params(axis='x', nbins=4)
        ax.locator_params(axis='y', nbins=4)
        
        # Save linear plot
        safe_region_name = region.replace(':', '_').replace('-', '_')
        plt.savefig(f"{output_dir}/{safe_region_name}.svg", 
                   bbox_inches='tight', pad_inches=0.05)
        plt.close()
        
        # Create log plot with same setup (plot area will be 3x2)
        fig = plt.figure(figsize=(5, 2))
        ax = fig.add_axes([0.1, 0.15, 0.6, 0.75])  # [left, bottom, width, height] in figure coordinates
        
        for sample in samples:
            sample_data = region_data[region_data['sample'] == sample]
            color = get_sample_color(sample, color_dict)
            # Add 1 to avoid log(0) issues
            ax.plot(sample_data['position'], sample_data['coverage'] + 1, 
                   label=sample, linewidth=0.75, color=color)
        
        ax.set_xlabel('Position')
        ax.set_ylabel('Coverage + 1 (log scale)')
        ax.set_yscale('log')
        ax.set_title(f'{region} (Log Scale)')
        
        # Position legend outside plot area
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)
        ax.grid(True, alpha=0.3, linewidth=0.5)
        
        # Reduce tick density for cleaner look
        ax.locator_params(axis='x', nbins=4)
        
        plt.savefig(f"{output_dir}/{safe_region_name}_log.svg", 
                   bbox_inches='tight', pad_inches=0.05)
        plt.close()
        
        print(f"Saved: {safe_region_name}.svg and {safe_region_name}_log.svg")

def main():
    parser = argparse.ArgumentParser(description='Plot rRNA coverage lines for each region (Illustrator-ready)')
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
