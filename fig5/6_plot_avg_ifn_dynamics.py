import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import os
import re
from tqdm import tqdm

def find_all_simulation_dirs(root_dir='.'):
    """Find all simulation directories across all replicate folders."""
    # Correcting the path based on the user-provided directory structure
    base_path = os.path.join(root_dir, 'IFNclr3_30runs_global_celltocell_tau95_option1', 'IFNclr3_30runs_global_celltocell_tau95_option1')
    
    if not os.path.isdir(base_path):
        print(f"Error: The target directory does not exist: {base_path}")
        return []

    sim_pattern = os.path.join(base_path, '[0-9]*_DIPBst*')
    all_sim_dirs = glob.glob(sim_pattern)
    return all_sim_dirs

def group_simulations_by_burst_size(sim_dirs, target_burst_sizes):
    """Groups simulation file paths by their DIP burst size."""
    grouped_files = {size: [] for size in target_burst_sizes}
    
    for sim_dir in sim_dirs:
        base_name = os.path.basename(sim_dir)
        bst_match = re.search(r'DIPBst(\d+)', base_name)
        if bst_match:
            burst_size = int(bst_match.group(1))
            if burst_size in grouped_files:
                csv_path = os.path.join(sim_dir, 'simulation_output.csv')
                if os.path.exists(csv_path):
                    grouped_files[burst_size].append(csv_path)
    
    return grouped_files

def process_data(grouped_files):
    """Calculate average dynamics and average max IFN for each burst size."""
    avg_dynamics = {}
    avg_max_ifn = {}

    for burst_size, paths in tqdm(grouped_files.items(), desc="Processing burst sizes"):
        if not paths:
            continue

        # Process for dynamics (left plot)
        all_dfs = [pd.read_csv(p) for p in paths]
        combined_df = pd.concat(all_dfs)
        
        # Ensure time points are consistent for averaging
        # This assumes a common time grid or we can interpolate if needed
        # For simplicity, we group by time and average.
        time_series_avg = combined_df.groupby('Time')['Global IFN Concentration Per Cell'].mean()
        avg_dynamics[burst_size] = time_series_avg

        # Process for max IFN (right plot)
        max_ifns = [df['Global IFN Concentration Per Cell'].max() for df in all_dfs]
        avg_max_ifn[burst_size] = np.mean(max_ifns)
        
    return avg_dynamics, avg_max_ifn

def plot_results(avg_dynamics, avg_max_ifn):
    """Generate a plot exactly reproducing the provided style."""
    sns.set_style("white")
    fig, ax1 = plt.subplots(1, 1, figsize=(10, 10))

    # Find the burst size with the highest average max IFN to highlight it
    burst_size_with_max_ifn = max(avg_max_ifn, key=avg_max_ifn.get)
    highlight_yield = burst_size_with_max_ifn / 50
    print(f"Highlighting Relative DIP Yield {highlight_yield:g}, which has the highest average peak IFN.")

    # --- Plot: IFN Dynamics ---
    cmap = plt.get_cmap('Blues')
    burst_sizes = sorted(avg_dynamics.keys())
    colors = cmap(np.linspace(0.3, 1.0, len(burst_sizes)))

    for i, burst_size in enumerate(burst_sizes):
        ts = avg_dynamics[burst_size]
        color = 'red' if burst_size == burst_size_with_max_ifn else colors[i]
        linewidth = 6.0 if burst_size == burst_size_with_max_ifn else 2.0
        alpha = 1.0 if burst_size == burst_size_with_max_ifn else 0.8
        
        relative_yield = burst_size / 50
        # Use only the number for the label, and a title for the legend
        ax1.plot(ts.index, ts.values, color=color, linewidth=linewidth, alpha=alpha, label=f'{relative_yield:g}')

    ax1.set_xlabel('Time', fontsize=24)
    ax1.set_ylabel('IFN Concentration', fontsize=24)
    ax1.tick_params(axis='both', which='major', labelsize=18)
    ax1.set_xlim(0, 500)
    ax1.legend(title='Relative DIP Yield', loc='upper right', frameon=True, title_fontsize='20', fontsize='18')

    # Set y-axis limits
    ax1.set_ylim(0, 0.22)
    
    plt.tight_layout()
    
    output_filename = '6_reproduced_avg_ifn_plot.png'
    plt.savefig(output_filename, dpi=300)
    print(f"\nPlot saved as {output_filename}")


def main():
    target_burst_sizes = list(range(100, 1601, 100))
    
    print("Finding all simulation directories...")
    all_sim_dirs = find_all_simulation_dirs()
    
    if not all_sim_dirs:
        print("Error: No simulation directories found.")
        return
        
    print(f"Found {len(all_sim_dirs)} total simulation runs.")
    
    grouped_files = group_simulations_by_burst_size(all_sim_dirs, target_burst_sizes)
    
    avg_dynamics, avg_max_ifn = process_data(grouped_files)
    
    if not avg_dynamics:
        print("Error: No data processed. Check if CSV files and correct burst sizes are present.")
        return
        
    plot_results(avg_dynamics, avg_max_ifn)

if __name__ == '__main__':
    main()
