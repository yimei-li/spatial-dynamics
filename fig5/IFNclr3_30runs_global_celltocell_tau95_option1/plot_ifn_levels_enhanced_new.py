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
import os

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
    
    # Get all summary files from the current directory
    summary_dir = '.'
    summary_files = [os.path.join(summary_dir, f) for f in os.listdir(summary_dir) if f.startswith('summary_DIPBst') and f.endswith('.csv') and os.path.isfile(os.path.join(summary_dir, f))]
    summary_files.sort()
    
    print(f"Found {len(summary_files)} summary files")
    
    # Create the plot
    plt.figure(figsize=(12, 8))
    
    # Colors for different lines
    colors = plt.cm.viridis(np.linspace(0, 1, len(summary_files)))
    
    # Track maximum IFN values to find the highest
    max_ifn_values = []
    dipbst_values = []
    
    # Plot each summary file
    for i, file in enumerate(summary_files):
        try:
            df = pd.read_csv(file)
            dipbst_value = extract_dipbst_from_filename(file)
            
            if 'Global IFN Concentration Per Cell' in df.columns:
                # Find the maximum IFN value for this DIPBst
                max_ifn = df['Global IFN Concentration Per Cell'].max()
                max_ifn_values.append(max_ifn)
                dipbst_values.append(dipbst_value)
                
                # Plot the line
                plt.plot(df['Time'], df['Global IFN Concentration Per Cell'], 
                        color=colors[i], 
                        linewidth=2, 
                        label=f'DIPBst = {dipbst_value}')
                
                print(f"DIPBst = {dipbst_value}, Max IFN = {max_ifn:.6f}")
            else:
                print(f"Warning: No 'Global IFN Concentration Per Cell' column found in {file}")
                
        except Exception as e:
            print(f"Error processing {file}: {e}")
    
    # Find the DIPBst with maximum IFN
    if max_ifn_values:
        max_idx = np.argmax(max_ifn_values)
        max_dipbst = dipbst_values[max_idx]
        max_ifn = max_ifn_values[max_idx]
        
        # Replot the maximum line in red
        for i, file in enumerate(summary_files):
            df = pd.read_csv(file)
            dipbst_value = extract_dipbst_from_filename(file)
            
            if dipbst_value == max_dipbst and 'Global IFN Concentration Per Cell' in df.columns:
                plt.plot(df['Time'], df['Global IFN Concentration Per Cell'], 
                        color='red', 
                        linewidth=3, 
                        label=f'DIPBst = {dipbst_value} (Max IFN)')
                break
        
        print(f"\nMaximum IFN: DIPBst = {max_dipbst}, IFN = {max_ifn:.6f}")
    
    # Customize the plot
    plt.xlabel('Time', fontsize=14)
    plt.ylabel('IFN Level', fontsize=14)
    plt.title('IFN Levels Over Time for Different DIPBst Values', fontsize=16, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    
    # Save the plot
    plt.savefig('ifn_levels_comparison_enhanced_new.png', dpi=300, bbox_inches='tight')
    plt.savefig('ifn_levels_comparison_enhanced_new.pdf', bbox_inches='tight')
    plt.show()
    
    print("\nPlots saved as:")
    print("- ifn_levels_comparison_enhanced_new.png")
    print("- ifn_levels_comparison_enhanced_new.pdf")

if __name__ == '__main__':
    plot_ifn_levels_enhanced() 