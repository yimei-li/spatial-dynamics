#!/usr/bin/env python3
"""
Script to process CSV files from 420 folders, group them by DIPBst values,
and calculate averages for each group of 30 folders to generate 14 summary CSV files.
"""

import os
import pandas as pd
import re
import numpy as np
from collections import defaultdict
import glob

def extract_dipbst_from_folder_name(folder_name):
    """
    Extract DIPBst value from folder name.
    Example: "1_Dinit0_DIPBst100_noJ_Vinit1_VBst50_Global_mdbk_times502_tau95_ifnBothFold1.00_grid50_VStimulateIFNtrue"
    should return 100
    """
    match = re.search(r'DIPBst(\d+)', folder_name)
    if match:
        return int(match.group(1))
    return None

def get_all_folders():
    """
    Get all folders in the current directory that match the pattern.
    """
    folders = []
    for item in os.listdir('.'):
        if os.path.isdir(item) and 'DIPBst' in item:
            folders.append(item)
    return sorted(folders)

def group_folders_by_dipbst(folders):
    """
    Group folders by their DIPBst values.
    """
    groups = defaultdict(list)
    for folder in folders:
        dipbst = extract_dipbst_from_folder_name(folder)
        if dipbst is not None:
            groups[dipbst].append(folder)
    return groups

def process_csv_file(csv_path):
    """
    Read and process a single CSV file.
    """
    try:
        df = pd.read_csv(csv_path)
        return df
    except Exception as e:
        print(f"Error reading {csv_path}: {e}")
        return None

def calculate_averages_for_group(folders, group_name):
    """
    Calculate averages for a group of folders (should be 30 folders).
    """
    print(f"Processing group {group_name} with {len(folders)} folders...")
    
    all_dataframes = []
    
    for folder in folders:
        csv_path = os.path.join(folder, 'simulation_output.csv')
        if os.path.exists(csv_path):
            df = process_csv_file(csv_path)
            if df is not None:
                all_dataframes.append(df)
        else:
            print(f"Warning: {csv_path} not found")
    
    if not all_dataframes:
        print(f"Error: No valid CSV files found for group {group_name}")
        return None
    
    print(f"Successfully loaded {len(all_dataframes)} CSV files for group {group_name}")
    
    # Calculate averages across all dataframes
    # For each column, calculate the mean across all dataframes
    averaged_data = {}
    
    # Get the first dataframe to get column names
    first_df = all_dataframes[0]
    
    for column in first_df.columns:
        if column == 'Time':
            # Keep time column as is (should be the same across all files)
            averaged_data[column] = first_df[column]
        else:
            # Only average numeric columns, keep non-numeric as the first value
            try:
                # Try converting all values to float
                column_values = []
                for df in all_dataframes:
                    if column in df.columns:
                        column_values.append(pd.to_numeric(df[column], errors='coerce').values)
                # If at least one value is numeric, do the mean
                if any([np.issubdtype(np.array(vals).dtype, np.number) for vals in column_values]):
                    stacked_values = np.column_stack(column_values)
                    averaged_data[column] = np.nanmean(stacked_values, axis=1)
                else:
                    # Non-numeric column, just use the first file's value
                    averaged_data[column] = first_df[column]
            except Exception as e:
                # If any error, fallback to first file's value
                averaged_data[column] = first_df[column]
    
    # Create the averaged dataframe
    averaged_df = pd.DataFrame(averaged_data)
    
    return averaged_df

def main():
    """
    Main function to process all folders and generate summary CSV files.
    """
    print("Starting CSV processing...")
    
    # Get all folders
    folders = get_all_folders()
    print(f"Found {len(folders)} folders")
    
    # Group folders by DIPBst values
    groups = group_folders_by_dipbst(folders)
    print(f"Found {len(groups)} different DIPBst values: {sorted(groups.keys())}")
    
    # Process each group
    for dipbst_value, group_folders in sorted(groups.items()):
        print(f"\nProcessing DIPBst = {dipbst_value}")
        print(f"Number of folders in this group: {len(group_folders)}")
        
        if len(group_folders) != 30:
            print(f"Warning: Expected 30 folders, found {len(group_folders)}")
        
        # Calculate averages for this group
        averaged_df = calculate_averages_for_group(group_folders, f"DIPBst{dipbst_value}")
        
        if averaged_df is not None:
            # Save the averaged data
            output_filename = f"summary_DIPBst{dipbst_value}.csv"
            averaged_df.to_csv(output_filename, index=False)
            print(f"Saved {output_filename}")
        else:
            print(f"Failed to process group DIPBst{dipbst_value}")
    
    print("\nProcessing complete!")

if __name__ == "__main__":
    main() 