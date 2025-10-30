#!/usr/bin/env python3
"""
Script to process simulation_output.csv files from IFNclr3_30runs_global_celltocell_tau95_option1_2,
average each group of 30 for each DIPBst value, and output summary csv files to output folder.
"""
import os
import pandas as pd
import re
import numpy as np
from collections import defaultdict

# 只处理这11个DIPBst值
burstSizeD_list = [650, 680, 685, 690, 695, 700, 705, 710, 715, 720, 750]

# 输入和输出文件夹
input_dir = 'IFNclr3_30runs_global_celltocell_tau95_option1_2'
output_dir = os.path.join(input_dir, 'output')
os.makedirs(output_dir, exist_ok=True)

def extract_dipbst_from_folder_name(folder_name):
    match = re.search(r'DIPBst(\d+)', folder_name)
    if match:
        return int(match.group(1))
    return None

def group_folders_by_dipbst(folders, selected_list):
    grouped = defaultdict(list)
    for folder in folders:
        dipbst = extract_dipbst_from_folder_name(folder)
        if dipbst in selected_list:
            grouped[dipbst].append(folder)
    return grouped

def main():
    # 获取所有子文件夹
    all_folders = [f for f in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, f))]
    grouped = group_folders_by_dipbst(all_folders, burstSizeD_list)
    print(f"Found DIPBst groups: { {k: len(v) for k,v in grouped.items()} }")

    for dipbst in burstSizeD_list:
        folders = grouped.get(dipbst, [])
        if len(folders) == 0:
            print(f"Warning: No folders found for DIPBst={dipbst}")
            continue
        if len(folders) != 30:
            print(f"Warning: DIPBst={dipbst} has {len(folders)} folders (expected 30)")
        all_dataframes = []
        for folder in folders:
            csv_path = os.path.join(input_dir, folder, 'simulation_output.csv')
            if os.path.exists(csv_path):
                try:
                    df = pd.read_csv(csv_path)
                    all_dataframes.append(df)
                except Exception as e:
                    print(f"Error reading {csv_path}: {e}")
            else:
                print(f"Missing: {csv_path}")
        if not all_dataframes:
            print(f"No data for DIPBst={dipbst}")
            continue
        # 对所有dataframe按列做平均，只对数值型列
        first_df = all_dataframes[0]
        averaged_data = {}
        for column in first_df.columns:
            if column == 'Time':
                averaged_data[column] = first_df[column]
            else:
                try:
                    column_values = [pd.to_numeric(df[column], errors='coerce').values for df in all_dataframes if column in df.columns]
                    if any([np.issubdtype(np.array(vals).dtype, np.number) for vals in column_values]):
                        stacked_values = np.column_stack(column_values)
                        averaged_data[column] = np.nanmean(stacked_values, axis=1)
                    else:
                        averaged_data[column] = first_df[column]
                except Exception as e:
                    averaged_data[column] = first_df[column]
        summary_df = pd.DataFrame(averaged_data)
        out_path = os.path.join(output_dir, f'summary_DIPBst{dipbst}.csv')
        summary_df.to_csv(out_path, index=False)
        print(f"Saved summary for DIPBst={dipbst} to {out_path}")

if __name__ == '__main__':
    main() 