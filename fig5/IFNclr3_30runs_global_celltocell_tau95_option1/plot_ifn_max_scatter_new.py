#!/usr/bin/env python3
"""
Script to create a stem plot of IFN maximum values vs DIPBst values.
Each point represents a DIPBst value with its corresponding maximum IFN level.
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

def plot_ifn_max_stem():
    """Create stem plot of IFN maximum values vs DIPBst values."""
    
    # Set style for better looking plots
    plt.style.use('default')
    
    # Get all summary files from the current directory
    summary_dir = '.'
    summary_files = [os.path.join(summary_dir, f) for f in os.listdir(summary_dir) if f.startswith('summary_DIPBst') and f.endswith('.csv') and os.path.isfile(os.path.join(summary_dir, f))]
    summary_files.sort()
    
    print(f"Found {len(summary_files)} summary files")
    
    # Extract DIPBst values and maximum IFN values
    dipbst_values = []
    max_ifn_values = []
    
    for file in summary_files:
        try:
            df = pd.read_csv(file)
            dipbst_value = extract_dipbst_from_filename(file)
            
            if 'Global IFN Concentration Per Cell' in df.columns:
                max_ifn = df['Global IFN Concentration Per Cell'].max()
                dipbst_values.append(dipbst_value)
                max_ifn_values.append(max_ifn)
                print(f"DIPBst = {dipbst_value}, Max IFN = {max_ifn:.6f}")
            else:
                print(f"Warning: No 'Global IFN Concentration Per Cell' column found in {file}")
                
        except Exception as e:
            print(f"Error processing {file}: {e}")
    
    if not dipbst_values:
        print("No valid data found!")
        return
    
    # Find the maximum IFN value
    max_idx = np.argmax(max_ifn_values)
    max_dipbst = dipbst_values[max_idx]
    max_ifn = max_ifn_values[max_idx]
    
    print(f"\nMaximum IFN: DIPBst = {max_dipbst}, IFN = {max_ifn:.6f}")
    
    # Create the stem plot
    plt.figure(figsize=(12, 8))
    
    # Plot stems for all points (black)
    for i, (dipbst, max_ifn) in enumerate(zip(dipbst_values, max_ifn_values)):
        if dipbst == max_dipbst:
            # Red for maximum
            plt.plot([dipbst, dipbst], [0, max_ifn], 'r-', linewidth=3, alpha=0.7)
            plt.plot(dipbst, max_ifn, 'ro', markersize=10, label=f'DIPBst = {dipbst} (Max IFN)' if i == max_idx else "")
        else:
            # Black for others
            plt.plot([dipbst, dipbst], [0, max_ifn], 'k-', linewidth=2, alpha=0.5)
            plt.plot(dipbst, max_ifn, 'ko', markersize=8)
    
    # Customize the plot
    plt.xlabel('DIPBst Value', fontsize=14)
    plt.ylabel('IFN Maximum Value', fontsize=14)
    plt.title('IFN Maximum Values vs DIPBst Values', fontsize=16, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    
    # Save the plot
    plt.savefig('ifn_max_vs_dipbst_stem_new.png', dpi=300, bbox_inches='tight')
    plt.savefig('ifn_max_vs_dipbst_stem_new.pdf', bbox_inches='tight')
    plt.show()
    
    print("\nPlots saved as:")
    print("- ifn_max_vs_dipbst_stem_new.png")
    print("- ifn_max_vs_dipbst_stem_new.pdf")

if __name__ == '__main__':
    plot_ifn_max_stem() 