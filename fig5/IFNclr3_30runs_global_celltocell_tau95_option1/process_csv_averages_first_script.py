#!/usr/bin/env python3
"""
Process CSV files from the first script directory and generate summary files
"""

import os
import pandas as pd
import numpy as np
from collections import defaultdict
import re

def extract_dipbst_from_folder(folder_name):
    """Extract DIPBst value from folder name"""
    match = re.search(r'DIPBst(\d+)', folder_name)
    if match:
        return int(match.group(1))
    return None

def main():
    # Directory containing the simulation results
    base_dir = 'IFNclr3_30runs_global_celltocell_tau95_option1'
    
    # Group folders by DIPBst value
    dipbst_groups = defaultdict(list)
    
    # Find all folders that contain CSV files
    for item in os.listdir(base_dir):
        if os.path.isdir(os.path.join(base_dir, item)) and item.startswith(('4', '5', '6', '7', '8', '9')):
            # Check if folder contains CSV files
            folder_path = os.path.join(base_dir, item)
            csv_files = [f for f in os.listdir(folder_path) if f.endswith('.csv')]
            if csv_files:
                dipbst = extract_dipbst_from_folder(item)
                if dipbst:
                    dipbst_groups[dipbst].append(item)
    
    print(f"Found {len(dipbst_groups)} different DIPBst values")
    
    # Process each DIPBst group
    for dipbst in sorted(dipbst_groups.keys()):
        folders = dipbst_groups[dipbst]
        print(f"Processing DIPBst {dipbst} with {len(folders)} folders...")
        
        # Collect all CSV files for this DIPBst
        all_data = []
        
        for folder in folders:
            folder_path = os.path.join(base_dir, folder)
            csv_files = [f for f in os.listdir(folder_path) if f.endswith('.csv')]
            for csv_file in csv_files:
                file_path = os.path.join(folder_path, csv_file)
                try:
                    df = pd.read_csv(file_path)
                    all_data.append(df)
                except Exception as e:
                    print(f"Error reading {file_path}: {e}")
        
        if all_data:
            # Combine all data
            combined_df = pd.concat(all_data, ignore_index=True)
            
            # Calculate averages for numeric columns
            summary_data = {}
            
            for column in combined_df.columns:
                if combined_df[column].dtype in ['int64', 'float64']:
                    summary_data[column] = combined_df[column].mean()
                else:
                    # For non-numeric columns, take the first value
                    summary_data[column] = combined_df[column].iloc[0]
            
            # Create summary DataFrame
            summary_df = pd.DataFrame([summary_data])
            
            # Save summary file
            output_filename = f'summary_DIPBst{dipbst}.csv'
            summary_df.to_csv(output_filename, index=False)
            print(f"Created {output_filename}")
            
            # Print max IFN value
            if 'Global IFN Concentration Per Cell' in summary_df.columns:
                max_ifn = summary_df['Global IFN Concentration Per Cell'].iloc[0]
                print(f"  Max IFN for DIPBst {dipbst}: {max_ifn}")
        else:
            print(f"No valid CSV data found for DIPBst {dipbst}")
    
    print("\nSummary files generation completed!")

if __name__ == "__main__":
    main() 