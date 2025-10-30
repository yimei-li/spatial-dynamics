#!/usr/bin/env python3
"""
Script to plot IFN levels over time from summary CSV files.
Each line represents a different DIPBst value.
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

def plot_ifn_levels():
    """Plot IFN levels over time for all summary files."""
    
    # Get all summary files
    summary_files = glob.glob('summary_DIPBst*.csv')
    summary_files.sort()
    
    print(f"Found {len(summary_files)} summary files")
    
    # Create the plot
    plt.figure(figsize=(12, 8))
    
    # Colors for different lines
    colors = plt.cm.viridis(np.linspace(0, 1, len(summary_files)))
    
    # Plot each summary file
    for i, file in enumerate(summary_files):
        try:
            # Read the CSV file
            df = pd.read_csv(file)
            
            # Extract DIPBst value from filename
            dipbst_value = extract_dipbst_from_filename(file)
            
            # Get time and IFN level data
            time = df['Time']
            
            # Look for IFN-related columns
            ifn_columns = [col for col in df.columns if 'IFN' in col.upper()]
            print(f"File {file}: IFN columns found: {ifn_columns}")
            
            # Use the first IFN column found (usually 'Global IFN Concentration Per Cell')
            if ifn_columns:
                ifn_level = df[ifn_columns[0]]
                
                # Plot the line
                plt.plot(time, ifn_level, 
                        label=f'DIPBst = {dipbst_value}', 
                        color=colors[i], 
                        linewidth=2)
                
                print(f"Plotted DIPBst = {dipbst_value}")
            else:
                print(f"Warning: No IFN columns found in {file}")
                
        except Exception as e:
            print(f"Error processing {file}: {e}")
    
    # Customize the plot
    plt.xlabel('Time', fontsize=14)
    plt.ylabel('Global IFN Concentration Per Cell', fontsize=14)
    plt.title('IFN Levels Over Time for Different DIPBst Values', fontsize=16)
    plt.grid(True, alpha=0.3)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
    
    # Adjust layout to prevent label cutoff
    plt.tight_layout()
    
    # Save the plot
    plt.savefig('ifn_levels_comparison.png', dpi=300, bbox_inches='tight')
    print("Plot saved as 'ifn_levels_comparison.png'")
    
    # Show the plot
    plt.show()

if __name__ == "__main__":
    plot_ifn_levels()