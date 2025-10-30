import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import os
import re
from tqdm import tqdm
from pygam import LinearGAM, s, intercept
from statsmodels.nonparametric.smoothers_lowess import lowess

# --- Plotting Configuration ---
sns.set_style("white")
TITLE_FONTSIZE = 36
LABEL_FONTSIZE = 48
TICK_FONTSIZE = 36
LEGEND_TITLE_FONTSIZE = 30
LEGEND_FONTSIZE = 27

# --- Data Processing Functions (from script 6) ---
def find_all_simulation_dirs(root_dir='.'):
    """Find all simulation directories based on the known structure."""
    base_path = os.path.join(root_dir, 'IFNclr3_30runs_global_celltocell_tau95_option1', 'IFNclr3_30runs_global_celltocell_tau95_option1')
    if not os.path.isdir(base_path): return []
    return glob.glob(os.path.join(base_path, '[0-9]*_DIPBst*'))

def group_simulations_by_burst_size(sim_dirs, target_burst_sizes):
    grouped_files = {size: [] for size in target_burst_sizes}
    for sim_dir in sim_dirs:
        bst_match = re.search(r'DIPBst(\d+)', os.path.basename(sim_dir))
        if bst_match and int(bst_match.group(1)) in grouped_files:
            csv_path = os.path.join(sim_dir, 'simulation_output.csv')
            if os.path.exists(csv_path):
                grouped_files[int(bst_match.group(1))].append(csv_path)
    return grouped_files

def process_avg_data(grouped_files):
    avg_dynamics, avg_max_ifn = {}, {}
    for burst_size, paths in tqdm(grouped_files.items(), desc="Averaging Dynamics"):
        if not paths: continue
        all_dfs = [pd.read_csv(p) for p in paths]
        combined_df = pd.concat(all_dfs)
        avg_dynamics[burst_size] = combined_df.groupby('Time')['Global IFN Concentration Per Cell'].mean()
        avg_max_ifn[burst_size] = np.mean([df['Global IFN Concentration Per Cell'].max() for df in all_dfs])
    return avg_dynamics, avg_max_ifn

def process_peak_data(grouped_files):
    """Processes individual runs to get peak IFN data for the right panel."""
    peak_data = []
    for burst_size, paths in tqdm(grouped_files.items(), desc="Extracting Peak IFN"):
        if not paths: continue
        for csv_path in paths:
            try:
                df = pd.read_csv(csv_path)
                if not df.empty:
                    peak_idx = df['Global IFN Concentration Per Cell'].idxmax()
                    peak_ifn = df.loc[peak_idx, 'Global IFN Concentration Per Cell']
                    peak_data.append({'burst_size_DIP': burst_size, 'peak_IFN': peak_ifn})
            except Exception as e:
                print(f"Warning: Could not process {csv_path}: {e}")
    return pd.DataFrame(peak_data)


# --- Analysis Functions (from script 5) ---
def find_optimum(model, X_data):
    grid = np.linspace(X_data.min(), X_data.max(), 500)
    preds = model.predict(grid)
    return grid[np.argmax(preds)]

def main():
    # --- Part 1: IFN Dynamics Plot (Left) ---
    target_burst_sizes = list(range(100, 1601, 100))
    all_sim_dirs = find_all_simulation_dirs()
    grouped_files = group_simulations_by_burst_size(all_sim_dirs, target_burst_sizes)
    avg_dynamics, avg_max_ifn = process_avg_data(grouped_files)

    # --- Part 2: Peak IFN Analysis Plot (Right) ---
    # Generate peak data on the fly from the same filtered set of runs
    df_peak = process_peak_data(grouped_files)
    if df_peak.empty:
        print("Error: No peak IFN data could be generated for the right panel.")
        return
        
    VIRION_BURST_SIZE = 50
    df_peak['burst_size_ratio'] = df_peak['burst_size_DIP'] / VIRION_BURST_SIZE
    X_peak = df_peak[['burst_size_ratio']].values
    y_peak = df_peak['peak_IFN'].values
    gam = LinearGAM(s(0, basis='ps', lam=0.6)).fit(X_peak, y_peak)
    loess_fit = lowess(y_peak, X_peak.ravel(), frac=0.5, it=1)
    b_star_observed = find_optimum(gam, X_peak)
    
    # --- Combined Plotting ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(30, 15))

    # Plot 1: IFN Dynamics
    burst_size_with_max_ifn = max(avg_max_ifn, key=avg_max_ifn.get)
    cmap = plt.get_cmap('Blues')
    burst_sizes = sorted(avg_dynamics.keys())
    colors = cmap(np.linspace(0.3, 1.0, len(burst_sizes)))

    for i, burst_size in enumerate(burst_sizes):
        ts = avg_dynamics[burst_size]
        color = 'red' if burst_size == burst_size_with_max_ifn else colors[i]
        linewidth = 6.0 if burst_size == burst_size_with_max_ifn else 2.0
        relative_yield = burst_size / 50
        ax1.plot(ts.index, ts.values, color=color, linewidth=linewidth, alpha=0.8, label=f'{relative_yield:g}')

    ax1.set_xlabel('Time', fontsize=LABEL_FONTSIZE)
    ax1.set_ylabel('IFN Concentration', fontsize=LABEL_FONTSIZE)
    ax1.tick_params(axis='both', which='major', labelsize=TICK_FONTSIZE)
    ax1.set_xlim(0, 500)
    ax1.set_ylim(0, 0.22)
    ax1.legend(title='Relative\nDIP Yield', loc='upper right', frameon=True, title_fontsize=LEGEND_TITLE_FONTSIZE, fontsize=LEGEND_FONTSIZE)

    # Plot 2: Peak IFN Analysis
    ax2.scatter(df_peak['burst_size_ratio'], df_peak['peak_IFN'], alpha=0.3, label='Stochastic Replicates', color='steelblue')
    grid = gam.generate_X_grid(term=0, n=500)
    ax2.plot(grid, gam.predict(grid), color='red', linewidth=3, label='GAM Fit')
    ax2.fill_between(grid[:, 0], gam.confidence_intervals(grid)[:, 0], gam.confidence_intervals(grid)[:, 1], color='red', alpha=0.2, label='95% Confidence Band')
    ax2.plot(loess_fit[:, 0], loess_fit[:, 1], color='green', linestyle='--', linewidth=3, label='LOESS Fit')
    ax2.axvline(b_star_observed, color='black', linestyle=':', linewidth=2.5, label=f'Optimum Yield = {b_star_observed:g}')
    
    ax2.set_xlabel('Relative DIP Yield', fontsize=LABEL_FONTSIZE)
    ax2.set_ylabel('Peak IFN', fontsize=LABEL_FONTSIZE)
    ax2.tick_params(axis='both', which='major', labelsize=TICK_FONTSIZE)
    ax2.legend(loc='upper right', fontsize=LEGEND_FONTSIZE)

    plt.tight_layout()
    plt.savefig('7_combined_ifn_analysis.png', dpi=300)
    print("\nCombined plot saved as '7_combined_ifn_analysis.png'")

if __name__ == '__main__':
    main()
