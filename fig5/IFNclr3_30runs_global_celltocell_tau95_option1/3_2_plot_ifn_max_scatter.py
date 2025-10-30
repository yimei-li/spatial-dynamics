#!/usr/bin/env python3
"""
Script to create a scatter plot of IFN maximum values vs DIPBst values.
Each point represents a DIPBst value with its corresponding maximum IFN level.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob
import re

def extract_dipbst_from_filename(filename):
    """Extract DIPBst value from filename."""
    match = re.search(r'DIPBst(\d+)', filename)
    if match:
        return int(match.group(1))
    return None

def plot_ifn_max_scatter():
    """Create scatter plot of IFN maximum values vs DIPBst values."""
    
    # Set style for better looking plots
    plt.style.use('default')
    
    # Get all summary files
    summary_files = glob.glob('summary_DIPBst*.csv')
    summary_files.sort()
    
    print(f"Found {len(summary_files)} summary files")
    
    # Create the plot
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Store data for analysis
    dipbst_values = []
    max_ifn_values = []
    
    # Collect data from all files
    for file in summary_files:
        try:
            df = pd.read_csv(file)
            dipbst_value = extract_dipbst_from_filename(file)
            ifn_level = df['Global IFN Concentration Per Cell']
            max_ifn = np.max(ifn_level)
            
            dipbst_values.append(dipbst_value)
            max_ifn_values.append(max_ifn)
            
            print(f"DIPBst = {dipbst_value}, Max IFN = {max_ifn:.6f}")
            
        except Exception as e:
            print(f"Error processing {file}: {e}")
    
    # Sort by DIPBst value for proper ordering
    sorted_data = sorted(zip(dipbst_values, max_ifn_values))
    dipbst_values = [x[0] for x in sorted_data]
    max_ifn_values = [x[1] for x in sorted_data]
    
    # Find the maximum IFN value
    max_ifn_overall = max(max_ifn_values)
    max_ifn_index = max_ifn_values.index(max_ifn_overall)
    max_ifn_dipbst = dipbst_values[max_ifn_index]
    
    print(f"\nMaximum IFN value: {max_ifn_overall:.6f} at DIPBst = {max_ifn_dipbst}")
    
    # Plot each point with vertical line
    for i, (dipbst, max_ifn) in enumerate(zip(dipbst_values, max_ifn_values)):
        if max_ifn == max_ifn_overall:
            # Red for maximum IFN
            color = 'red'
            marker_size = 100
            linewidth = 3
            alpha = 1.0
            label = f'DIPBst = {dipbst} (Max IFN)'
        else:
            # Black for others
            color = 'black'
            marker_size = 60
            linewidth = 1.5
            alpha = 0.7
            label = None
        
        # Plot vertical line from 0 to max IFN
        ax.plot([dipbst, dipbst], [0, max_ifn], 
                color=color, linewidth=linewidth, alpha=alpha)
        
        # Plot point at max IFN
        ax.scatter(dipbst, max_ifn, 
                  color=color, s=marker_size, alpha=alpha, 
                  zorder=5, label=label)
    
    # Customize the plot
    ax.set_xlabel('DIPBst Value', fontsize=16, fontweight='bold')
    ax.set_ylabel('Maximum IFN Concentration Per Cell', fontsize=16, fontweight='bold')
    ax.set_title('Maximum IFN Levels vs DIPBst Values\n(30 runs averaged for each DIPBst value)', 
                fontsize=18, fontweight='bold', pad=20)
    
    # Add grid
    ax.grid(True, alpha=0.3, linestyle='--')
    
    # Set x-axis to show all DIPBst values
    ax.set_xticks(dipbst_values)
    ax.set_xticklabels(dipbst_values, rotation=45)
    
    # Create legend for the maximum point
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        legend = ax.legend(loc='upper left', fontsize=12,
                          frameon=True, fancybox=True, shadow=True)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save the plot
    plt.savefig('ifn_max_vs_dipbst_scatter.png', 
                dpi=300, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    print("\nScatter plot saved as 'ifn_max_vs_dipbst_scatter.png'")
    
    # Also save as PDF
    plt.savefig('ifn_max_vs_dipbst_scatter.pdf', 
                bbox_inches='tight',
                facecolor='white', edgecolor='none')
    print("Scatter plot also saved as 'ifn_max_vs_dipbst_scatter.pdf'")
    
    # Show the plot
    plt.show()

if __name__ == "__main__":
    plot_ifn_max_scatter() 