#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import os
from pathlib import Path

# Set style for publication-quality plots
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

def plot_rrna_info_content(file_path, output_dir="plots"):
    """
    Plot information content for a single rRNA file
    """
    # Extract rRNA type from filename
    filename = Path(file_path).name
    rrna_type = filename.split('.')[1].replace('_info_content', '')
    
    # Read the data
    df = pd.read_csv(file_path, sep='\t', header=None, names=['Position', 'Information_Content'])
    
    # Create the plot
    fig, ax = plt.subplots(figsize=(24, 6))
    
    # Plot the line
    ax.plot(df['Position'], df['Information_Content'], 
            linewidth=1.5, alpha=0.8)
    
    # Customize the plot
    ax.set_xlabel('Position', fontsize=12, fontweight='bold')
    ax.set_ylabel('Information Content (bits)', fontsize=12, fontweight='bold')
    ax.set_title(f'Information Content Along {rrna_type} rRNA Sequence', 
                fontsize=14, fontweight='bold', pad=20)
    
    # Add grid for better readability
    ax.grid(True, alpha=0.3)
    
    # Set y-axis limits if needed (information content is typically 0-2 bits)
    ax.set_ylim(bottom=0)
    
    # Adjust layout
    plt.tight_layout()
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Save the plot
    output_file = os.path.join(output_dir, f"Dmel_rRNA_{rrna_type}_info_content.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.savefig(output_file.replace('.png', '.pdf'), bbox_inches='tight')
    
    print(f"Plot saved: {output_file}")
    plt.show()
    
    return df

def plot_all_rrna_comparison(file_pattern, output_dir="plots"):
    """
    Create a comparison plot with all rRNA types on the same axes
    """
    files = glob.glob(file_pattern)
    
    if not files:
        print(f"No files found matching pattern: {file_pattern}")
        return
    
    fig, ax = plt.subplots(figsize=(24, 6))
    
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    
    for i, file_path in enumerate(sorted(files)):
        # Extract rRNA type from filename
        filename = Path(file_path).name
        rrna_type = filename.split('.')[1].replace('_info_content', '')
        
        # Read the data
        df = pd.read_csv(file_path, sep='\t', header=None, names=['Position', 'Information_Content'])
        
        # Plot the line
        ax.plot(df['Position'], df['Information_Content'], 
                linewidth=1.5, 
                label=f'{rrna_type} rRNA', color=colors[i % len(colors)], alpha=0.8)
    
    # Customize the plot
    ax.set_xlabel('Position', fontsize=12, fontweight='bold')
    ax.set_ylabel('Information Content (bits)', fontsize=12, fontweight='bold')
    ax.set_title('Information Content Comparison Across D. melanogaster rRNA Subunits', 
                fontsize=14, fontweight='bold', pad=20)
    
    # Add legend
    ax.legend(frameon=True, fancybox=True, shadow=True)
    
    # Add grid for better readability
    ax.grid(True, alpha=0.3)
    
    # Set y-axis limits
    ax.set_ylim(bottom=0)
    
    # Adjust layout
    plt.tight_layout()
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Save the comparison plot
    output_file = os.path.join(output_dir, "Dmel_rRNA_all_comparison_info_content.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.savefig(output_file.replace('.png', '.pdf'), bbox_inches='tight')
    
    print(f"Comparison plot saved: {output_file}")
    plt.show()

def main():
    # Define the file pattern for your rRNA files
    base_path = "/private/groups/russelllab/jodie/scRNAseq/wolbachia_titer/rRNA_analysis/"
    file_pattern = os.path.join(base_path, "Dmel_rRNA.*_info_content.txt")
    
    print("Creating individual plots for each rRNA subunit...")
    
    # Get all matching files
    files = glob.glob(file_pattern)
    
    if not files:
        print(f"No files found matching pattern: {file_pattern}")
        print("Please check the file path and pattern.")
        return
    
    # Create individual plots
    for file_path in sorted(files):
        print(f"Processing: {file_path}")
        df = plot_rrna_info_content(file_path)
        print(f"Data shape: {df.shape}")
        print(f"Info content range: {df['Information_Content'].min():.4f} - {df['Information_Content'].max():.4f}")
        print("-" * 50)
    
    print("\nCreating comparison plot...")
    plot_all_rrna_comparison(file_pattern)
    
    print("\nAll plots completed!")

if __name__ == "__main__":
    main()