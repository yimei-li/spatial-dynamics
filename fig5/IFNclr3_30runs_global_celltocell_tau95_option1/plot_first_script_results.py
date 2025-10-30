#!/usr/bin/env python3
"""
Plot results from the first script (DIPBst 100-1600)
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
import matplotlib.cm as cm

# Set font and style
plt.rcParams['font.size'] = 12
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['xtick.major.width'] = 1.5
plt.rcParams['ytick.major.width'] = 1.5

def main():
    # Get all summary files from the first script
    summary_files = [f for f in os.listdir('.') if f.startswith('summary_DIPBst') and f.endswith('.csv') and os.path.isfile(f)]
    summary_files.sort()
    
    print(f"Found {len(summary_files)} summary files")
    
    # Read and process data
    dipbst_values = []
    max_ifn_values = []
    all_data = []
    
    for file in summary_files:
        try:
            df = pd.read_csv(file)
            dipbst = int(file.split('_')[1].replace('DIPBst', '').replace('.csv', ''))
            
            if 'Global IFN Concentration Per Cell' in df.columns:
                max_ifn = df['Global IFN Concentration Per Cell'].max()
                dipbst_values.append(dipbst)
                max_ifn_values.append(max_ifn)
                all_data.append(df)
                print(f"DIPBst {dipbst}: Max IFN = {max_ifn:.6f}")
        except Exception as e:
            print(f"Error processing {file}: {e}")
    
    if not dipbst_values:
        print("No valid data found!")
        return
    
    # Sort by DIPBst values
    sorted_data = sorted(zip(dipbst_values, max_ifn_values))
    dipbst_values = [x[0] for x in sorted_data]
    max_ifn_values = [x[1] for x in sorted_data]
    
    # Find maximum IFN value
    max_ifn_overall = max(max_ifn_values)
    max_dipbst = dipbst_values[max_ifn_values.index(max_ifn_overall)]
    
    print(f"\nMaximum IFN: {max_ifn_overall:.6f} at DIPBst = {max_dipbst}")
    
    # Create figure with subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    # Plot 1: IFN levels over time for all DIPBst values
    print("\nGenerating IFN time series plot...")
    
    # Combine all data for time series
    combined_df = pd.concat(all_data, ignore_index=True)
    time_points = combined_df['Time'].unique()
    time_points.sort()
    
    # Create color map using matplotlib's blues
    min_dipbst = min(dipbst_values)
    max_dipbst_range = max(dipbst_values) - min_dipbst
    blues_cmap = cm.get_cmap('Blues')
    
    # Plot time series for each DIPBst
    for i, dipbst in enumerate(dipbst_values):
        # Get data for this DIPBst
        dipbst_data = []
        for df in all_data:
            if int(df['BURST_SIZE_D'].iloc[0]) == dipbst:
                dipbst_data.append(df)
        
        if dipbst_data:
            combined_dipbst = pd.concat(dipbst_data, ignore_index=True)
            ifn_values = []
            for time in time_points:
                time_data = combined_dipbst[combined_dipbst['Time'] == time]
                if not time_data.empty:
                    ifn_values.append(time_data['Global IFN Concentration Per Cell'].mean())
                else:
                    ifn_values.append(0)
            
            # Calculate color based on DIPBst value using blues colormap
            if max_dipbst_range > 0:
                color_intensity = (dipbst - min_dipbst) / max_dipbst_range
            else:
                color_intensity = 0.5
            
            # Get color from blues colormap (0.3 to 0.9 for better contrast)
            color = blues_cmap(0.3 + 0.6 * color_intensity)
            
            # Plot with different colors
            if dipbst == max_dipbst:
                ax1.plot(time_points, ifn_values, color='red', linewidth=3, label=f'DIPBst {dipbst} (Max)', zorder=10)
            else:
                ax1.plot(time_points, ifn_values, color=color, alpha=0.8, linewidth=1.5)
    
    ax1.set_xlabel('Time', fontsize=14)
    ax1.set_ylabel('Global IFN Concentration Per Cell', fontsize=14)
    ax1.set_title('IFN Levels Over Time - First Script Results (DIPBst 100-1600)\nColor: Light Blue (Low DIPBst) → Dark Blue (High DIPBst)', fontsize=16, fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: DIPBst vs Max IFN scatter plot
    print("Generating DIPBst vs Max IFN scatter plot...")
    
    # Plot each stem individually with different colors
    for i, (dipbst, ifn_val) in enumerate(zip(dipbst_values, max_ifn_values)):
        if dipbst == max_dipbst:
            # Highlight maximum point in red
            ax2.plot([dipbst], [ifn_val], 'ro', markersize=10, 
                     label=f'Max IFN: {max_ifn_overall:.6f} at DIPBst {max_dipbst}')
            ax2.vlines(dipbst, 0, ifn_val, colors='red', linewidth=2)
        else:
            # Use black for all other points
            ax2.plot([dipbst], [ifn_val], 'ko', markersize=6)
            ax2.vlines(dipbst, 0, ifn_val, colors='black', linewidth=1, alpha=0.7)
    
    # Set x-axis ticks to show all DIPBst values at 45 degree angle
    ax2.set_xticks(dipbst_values)
    ax2.set_xticklabels(dipbst_values, rotation=45, ha='right')
    
    ax2.set_xlabel('DIPBst Value', fontsize=14)
    ax2.set_ylabel('Maximum IFN Concentration', fontsize=14)
    ax2.set_title('DIPBst vs Maximum IFN - First Script Results', fontsize=16, fontweight='bold')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Adjust layout and save
    plt.tight_layout()
    
    # Save plots
    output_png = 'first_script_results.png'
    output_pdf = 'first_script_results.pdf'
    
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    plt.savefig(output_pdf, bbox_inches='tight')
    
    print(f"\nPlots saved as:")
    print(f"  {output_png}")
    print(f"  {output_pdf}")
    
    # Print summary statistics
    print(f"\n=== SUMMARY STATISTICS ===")
    print(f"Total DIPBst values: {len(dipbst_values)}")
    print(f"DIPBst range: {min(dipbst_values)} - {max(dipbst_values)}")
    print(f"IFN range: {min(max_ifn_values):.6f} - {max(max_ifn_values):.6f}")
    print(f"Optimal DIPBst: {max_dipbst} (IFN: {max_ifn_overall:.6f})")
    
    plt.show()

if __name__ == "__main__":
    main() 