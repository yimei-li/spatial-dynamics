import os
import re
import pandas as pd
import glob
import sys

def find_replicate_dirs(root_path):
    """Finds all replicate directories."""
    pattern = os.path.join(root_path, 'IFNclr3_30runs_global_celltocell_tau95_option1*')
    all_paths = glob.glob(pattern)
    rep_dirs = [p for p in all_paths if os.path.isdir(p)]
    
    if not rep_dirs:
        print(f"Error: No replicate directories found matching pattern: {pattern}", file=sys.stderr)
        sys.exit(1)
        
    return rep_dirs

def extract_replicate_id(dir_name):
    """Extracts replicate ID from a directory name."""
    match = re.search(r'_(\d+)$', dir_name)
    if match:
        return int(match.group(1))
    # Check if it's the first replicate which might not have a suffix
    if 'option1' in os.path.basename(dir_name) and not re.search(r'_\d+$', dir_name):
        return 1
    return None

def process_simulation_dir(sim_dir):
    """Processes a single simulation directory to extract peak IFN data."""
    base_name = os.path.basename(sim_dir)
    
    # Extract burst size
    bst_match = re.search(r'DIPBst(\d+)', base_name)
    if not bst_match:
        return None
    burst_size_DIP = int(bst_match.group(1))

    # Extract seed
    seed = -1
    seed_file = os.path.join(sim_dir, 'seed.txt')
    try:
        with open(seed_file, 'r') as f:
            seed = int(f.read().strip())
    except (IOError, ValueError):
        print(f"Warning: Could not read seed from {seed_file}", file=sys.stderr)

    # Process simulation output data
    data_file = os.path.join(sim_dir, 'simulation_output.csv')
    if not os.path.exists(data_file):
        print(f"Warning: Data file not found: {data_file}", file=sys.stderr)
        return None
        
    try:
        df = pd.read_csv(data_file)
        if df.empty:
            print(f"Warning: Data file is empty: {data_file}", file=sys.stderr)
            return None
        
        # The CSV contains a time series of global IFN
        time_col = 'Time'
        ifn_col = 'Global IFN Concentration Per Cell'

        if time_col not in df.columns or ifn_col not in df.columns:
            print(f"Warning: '{time_col}' or '{ifn_col}' column not found in {data_file}", file=sys.stderr)
            return None

        peak_idx = df[ifn_col].idxmax()
        peak_IFN = df.loc[peak_idx, ifn_col]
        t_peak_IFN = df.loc[peak_idx, time_col]
        
        return {
            'burst_size_DIP': burst_size_DIP,
            'seed': seed,
            'peak_IFN': peak_IFN,
            't_peak_IFN': t_peak_IFN,
            'scenario': 'baseline'
        }
    except Exception as e:
        print(f"Error processing {data_file}: {e}", file=sys.stderr)
        return None

def main():
    """Main function to orchestrate the processing."""
    results = []
    
    # Look for replicate directories inside the main experiment folder and also as siblings
    base_dir = 'IFNclr3_30runs_global_celltocell_tau95_option1'
    
    # Case 1: Replicate folders are inside the base folder
    # e.g., IFNclr3_30runs_global_celltocell_tau95_option1/IFNclr3_30runs_global_celltocell_tau95_option1_2/
    nested_replicate_dirs = find_replicate_dirs(base_dir)

    # Case 2: Replicate folders are siblings to the base folder (as seen in some logs)
    sibling_replicate_dirs = find_replicate_dirs('.')
    
    replicate_dirs = sorted(list(set(nested_replicate_dirs + sibling_replicate_dirs)))

    if not replicate_dirs:
        print("Error: No replicate directories found in known locations.", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(replicate_dirs)} replicate directories to process.")

    for rep_dir in replicate_dirs:
        replicate_id = extract_replicate_id(rep_dir)
        if replicate_id is None:
            print(f"Warning: Could not determine replicate ID for {rep_dir}", file=sys.stderr)
            continue
            
        sim_dirs = glob.glob(os.path.join(rep_dir, '[0-9]*_DIPBst*'))
        
        for sim_dir in sim_dirs:
            data = process_simulation_dir(sim_dir)
            if data:
                data['replicate_id'] = replicate_id
                results.append(data)

    if not results:
        print("Error: No simulation results could be processed.", file=sys.stderr)
        sys.exit(1)
        
    # Create and save the final DataFrame
    final_df = pd.DataFrame(results)

    # Sort by burst size
    final_df = final_df.sort_values(by='burst_size_DIP').reset_index(drop=True)

    # Ensure all required columns are present
    required_cols = ['burst_size_DIP', 'replicate_id', 'seed', 'peak_IFN', 't_peak_IFN', 'scenario']
    for col in required_cols:
        if col not in final_df.columns:
            final_df[col] = None # Add missing columns if any, though it should not happen with this logic.
            
    final_df = final_df[required_cols] # Order columns correctly
    
    output_filename = 'ifn_peak_vs_dipburst_baseline.csv'
    final_df.to_csv(output_filename, index=False, encoding='utf-8')
    
    print(f"Successfully processed {len(results)} simulations.")
    print(f"Results saved to {output_filename}")
    
    # Verification checks
    print("\nVerification:")
    print(f"Total rows: {len(final_df)}")
    print("Number of replicates per burst size:")
    print(final_df.groupby('burst_size_DIP')['replicate_id'].nunique())


if __name__ == '__main__':
    main()
