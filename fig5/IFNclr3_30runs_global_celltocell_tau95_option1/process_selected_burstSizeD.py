#!/usr/bin/env python3
"""
Script to process only selected burstSizeD (DIPBst) values, average each group of 30 csv files,
and output summary csv files to an output folder.
"""
import os
import pandas as pd
import re
import numpy as np
from collections import defaultdict

# 只处理这11个DIPBst值
burstSizeD_list = [650, 680, 685, 690, 695, 700, 705, 710, 715, 720, 750]

# 输出文件夹
output_dir = 'output'
os.makedirs(output_dir, exist_ok=True)

def extract_dipbst_from_folder_name(folder_name):
    match = re.search(r'DIPBst(\d+)', folder_name)
    if match:
        return int(match.group(1))
    return None

def get_all_folders():
    folders = []
    for item in os.listdir('.'):
        if os.path.isdir(item) and 'DIPBst' in item:
            folders.append(item)
    return sorted(folders)

def group_folders_by_dipbst(folders, selected_list):
    groups = defaultdict(list)
    for folder in folders:
        dipbst = extract_dipbst_from_folder_name(folder)
        if dipbst is not None and dipbst in selected_list:
            groups[dipbst].append(folder)
    return groups

def process_csv_file(csv_path):
    try:
        df = pd.read_csv(csv_path)
        return df
    except Exception as e:
        print(f"Error reading {csv_path}: {e}")
        return None

def calculate_averages_for_group(folders, group_name):
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
    averaged_data = {}
    first_df = all_dataframes[0]
    for column in first_df.columns:
        if column == 'Time':
            averaged_data[column] = first_df[column]
        else:
            try:
                column_values = []
                for df in all_dataframes:
                    if column in df.columns:
                        column_values.append(pd.to_numeric(df[column], errors='coerce').values)
                if any([np.issubdtype(np.array(vals).dtype, np.number) for vals in column_values]):
                    stacked_values = np.column_stack(column_values)
                    averaged_data[column] = np.nanmean(stacked_values, axis=1)
                else:
                    averaged_data[column] = first_df[column]
            except Exception as e:
                averaged_data[column] = first_df[column]
    averaged_df = pd.DataFrame(averaged_data)
    return averaged_df

def main():
    print("Starting selected burstSizeD processing...")
    folders = get_all_folders()
    print(f"Found {len(folders)} folders")
    groups = group_folders_by_dipbst(folders, burstSizeD_list)
    print(f"Found {len(groups)} selected DIPBst values: {sorted(groups.keys())}")
    for dipbst_value in burstSizeD_list:
        group_folders = groups.get(dipbst_value, [])
        print(f"\nProcessing DIPBst = {dipbst_value}")
        print(f"Number of folders in this group: {len(group_folders)}")
        if len(group_folders) != 30:
            print(f"Warning: Expected 30 folders, found {len(group_folders)}")
        averaged_df = calculate_averages_for_group(group_folders, f"DIPBst{dipbst_value}")
        if averaged_df is not None:
            output_filename = os.path.join(output_dir, f"summary_DIPBst{dipbst_value}.csv")
            averaged_df.to_csv(output_filename, index=False)
            print(f"Saved {output_filename}")
        else:
            print(f"Failed to process group DIPBst{dipbst_value}")
    print("\nProcessing complete!")

if __name__ == "__main__":
    main() 