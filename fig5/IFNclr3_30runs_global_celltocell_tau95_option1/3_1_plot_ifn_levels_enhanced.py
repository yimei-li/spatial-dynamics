#!/usr/bin/env python3
"""
Enhanced script to plot IFN levels over time from summary CSV files.
Each line represents a different DIPBst value with better visualization.
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

def plot_ifn_levels_enhanced():
    """Plot IFN levels over time for all summary files with enhanced styling."""
    
    # Set style for better looking plots
    plt.style.use('default')
    
    # Get all summary files
    summary_files = glob.glob('summary_DIPBst*.csv')
    summary_files.sort()
    
    print(f"Found {len(summary_files)} summary files")
    
    # Create the plot with larger size
    fig, ax = plt.subplots(figsize=(14, 10))
    
    # Store data for analysis
    all_data = []
    max_ifn_data = []
    
    # First pass: collect all data to find the maximum IFN line
    for file in summary_files:
        try:
            df = pd.read_csv(file)
            dipbst_value = extract_dipbst_from_filename(file)
            ifn_level = df['Global IFN Concentration Per Cell']
            
            all_data.append((dipbst_value, ifn_level))
            max_ifn_data.append((dipbst_value, np.max(ifn_level)))
            
        except Exception as e:
            print(f"Error processing {file}: {e}")
    
    # Find the DIPBst value with the maximum IFN
    max_ifn_dipbst = max(max_ifn_data, key=lambda x: x[1])[0]
    print(f"DIPBst with maximum IFN: {max_ifn_dipbst}")
    
    # Sort by DIPBst value for color gradient
    all_data.sort(key=lambda x: x[0])
    
    # Create blue color gradient (lightest for lowest DIPBst)
    dipbst_values = [data[0] for data in all_data]
    min_dipbst = min(dipbst_values)
    max_dipbst = max(dipbst_values)
    
    # Plot each summary file
    for i, (dipbst_value, ifn_level) in enumerate(all_data):
        try:
            # Read the CSV file to get time
            file = f"summary_DIPBst{dipbst_value}.csv"
            df = pd.read_csv(file)
            time = df['Time']
            
            # Determine color: red for max IFN, blue gradient for others
            if dipbst_value == max_ifn_dipbst:
                color = 'red'
                linewidth = 3.0
                alpha = 1.0
                label = f'DIPBst = {dipbst_value} (Max IFN)'
            else:
                # Blue gradient: lightest for lowest DIPBst
                normalized_value = (dipbst_value - min_dipbst) / (max_dipbst - min_dipbst)
                color = plt.cm.Blues(0.3 + 0.6 * normalized_value)  # Start from 0.3 to avoid too light
                linewidth = 2.0
                alpha = 0.7
                label = f'DIPBst = {dipbst_value}'
            
            # Plot the line
            line, = ax.plot(time, ifn_level, 
                           label=label, 
                           color=color, 
                           linewidth=linewidth,
                           alpha=alpha)
            
            print(f"Plotted DIPBst = {dipbst_value} (color: {'red' if dipbst_value == max_ifn_dipbst else 'blue'})")
                
        except Exception as e:
            print(f"Error processing DIPBst {dipbst_value}: {e}")
    
    # Customize the plot
    ax.set_xlabel('Time', fontsize=16, fontweight='bold')
    ax.set_ylabel('Global IFN Concentration Per Cell', fontsize=16, fontweight='bold')
    ax.set_title('IFN Levels Over Time for Different DIPBst Values\n(30 runs averaged for each DIPBst value)', 
                fontsize=18, fontweight='bold', pad=20)
    
    # Add grid
    ax.grid(True, alpha=0.3, linestyle='--')
    
    # Create legend with better positioning
    legend = ax.legend(loc='center left', 
                      bbox_to_anchor=(1.02, 0.5),
                      fontsize=11,
                      frameon=True,
                      fancybox=True,
                      shadow=True)
    
    # Adjust layout to prevent label cutoff
    plt.tight_layout()
    
    # Save the plot with high resolution
    plt.savefig('ifn_levels_comparison_enhanced.png', 
                dpi=300, 
                bbox_inches='tight',
                facecolor='white',
                edgecolor='none')
    print("Enhanced plot saved as 'ifn_levels_comparison_enhanced.png'")
    
    # Also save as PDF for vector graphics
    plt.savefig('ifn_levels_comparison_enhanced.pdf', 
                bbox_inches='tight',
                facecolor='white',
                edgecolor='none')
    print("Enhanced plot also saved as 'ifn_levels_comparison_enhanced.pdf'")
    
    # Show the plot
    plt.show()

if __name__ == "__main__":
    plot_ifn_levels_enhanced() 