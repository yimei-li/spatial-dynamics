import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

def create_ifn10_integrated_linear():
    kjumpr_0_file = "kjumpr_0_local_dip_averaged_simulation_output.csv"
    base_dir = "kjumpr_0.01_vary_ifn_range_local_dip"
    ifn10_folder = "10_Dinit0_DIPBst50_JRand_Vinit1_VBst100_IFN10_mdbk_times500_tau12_ifnBothFold1.00_grid76_VStimulateIFNtrue"
    ifn10_file = os.path.join(base_dir, ifn10_folder, "simulation_output.csv")
    
    fig, ax = plt.subplots(figsize=(18, 12))
    
    try:
        df_ifn10 = pd.read_csv(ifn10_file)
        ax.plot(df_ifn10['Time'], df_ifn10['Percentage Dead Cells'], 
               color='darkred', linewidth=10, 
               label='IFN10 Dead (with $\mathbf{1\%\ jump}$)', linestyle='-', alpha=0.9)
        ax.plot(df_ifn10['Time'], df_ifn10['Percentage Antiviral Cells'], 
               color='darkblue', linewidth=10, 
               label='IFN10 Antiviral (with $\mathbf{1\%\ jump}$)', linestyle='-', alpha=0.9)
    except Exception as e:
        print(f"Error reading IFN10 file: {e}")
    
    try:
        df_baseline = pd.read_csv(kjumpr_0_file)
        ax.plot(df_baseline['Time'], df_baseline['Percentage Dead Cells'], 
               color='darkorange', linewidth=12, 
               label='IFN10Dead ($\mathbf{no\ jump}$)', linestyle='--', alpha=0.95)
        ax.plot(df_baseline['Time'], df_baseline['Percentage Antiviral Cells'], 
               color='darkturquoise', linewidth=12, 
               label='IFN10Antiviral ($\mathbf{no\ jump}$)', linestyle='--', alpha=0.95)
    except Exception as e:
        print(f"Error reading baseline file: {e}")
    
    ax.text(0.02, 0.95, 'Red/Orange: Dead Cells', transform=ax.transAxes,
            fontsize=26, color='darkred', weight='bold', 
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax.text(0.02, 0.88, 'Blue/Turquoise: Antiviral Cells', transform=ax.transAxes,
            fontsize=26, color='darkblue', weight='bold',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    ax.set_xlabel('Time (hours)', fontsize=34, fontweight='normal')
    ax.set_ylabel('Percentage of Cells (%)', fontsize=34, fontweight='normal')
    ax.set_title('IFN10 Only: Dead and Antiviral Cells Dynamics (Linear Scale)\nComparing 1% Jump vs No Jump', 
                fontsize=38, fontweight='bold', pad=30)
    ax.tick_params(axis='both', which='major', labelsize=32)
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=22, 
             frameon=True, fancybox=True, shadow=True, ncol=1)
    ax.grid(True, linestyle='--', alpha=0.3)
    ax.set_facecolor('#f8f8f8')
    ax.set_ylim([0, 100])
    
    plt.tight_layout()
    output_filename = 'ifn10_only_integrated_linear_plot.png'
    plt.savefig(output_filename, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {output_filename}")

def create_ifn10_integrated_log10():
    kjumpr_0_file = "kjumpr_0_local_dip_averaged_simulation_output.csv"
    base_dir = "kjumpr_0.01_vary_ifn_range_local_dip"
    ifn10_folder = "10_Dinit0_DIPBst50_JRand_Vinit1_VBst100_IFN10_mdbk_times500_tau12_ifnBothFold1.00_grid76_VStimulateIFNtrue"
    ifn10_file = os.path.join(base_dir, ifn10_folder, "simulation_output.csv")
    
    fig, ax = plt.subplots(figsize=(18, 12))
    
    try:
        df_ifn10 = pd.read_csv(ifn10_file)
        dead_values = df_ifn10['Percentage Dead Cells'].values
        antiviral_values = df_ifn10['Percentage Antiviral Cells'].values
        dead_values[dead_values == 0] = 1e-6
        antiviral_values[antiviral_values == 0] = 1e-6
        
        ax.plot(df_ifn10['Time'], dead_values, 
               color='darkred', linewidth=10, 
               label='IFN10 Dead (with $\mathbf{1\%\ jump}$)', linestyle='-', alpha=0.9)
        ax.plot(df_ifn10['Time'], antiviral_values, 
               color='darkblue', linewidth=10, 
               label='IFN10 Antiviral (with $\mathbf{1\%\ jump}$)', linestyle='-', alpha=0.9)
    except Exception as e:
        print(f"Error reading IFN10 file: {e}")
    
    try:
        df_baseline = pd.read_csv(kjumpr_0_file)
        dead_values = df_baseline['Percentage Dead Cells'].values
        antiviral_values = df_baseline['Percentage Antiviral Cells'].values
        dead_values[dead_values == 0] = 1e-6
        antiviral_values[antiviral_values == 0] = 1e-6
        
        ax.plot(df_baseline['Time'], dead_values, 
               color='darkorange', linewidth=12, 
               label='IFN10Dead ($\mathbf{no\ jump}$)', linestyle='--', alpha=0.95)
        ax.plot(df_baseline['Time'], antiviral_values, 
               color='darkturquoise', linewidth=12, 
               label='IFN10Antiviral ($\mathbf{no\ jump}$)', linestyle='--', alpha=0.95)
    except Exception as e:
        print(f"Error reading baseline file: {e}")
    
    ax.set_yscale('log')
    
    ax.text(0.02, 0.95, 'Red/Orange: Dead Cells', transform=ax.transAxes,
            fontsize=26, color='darkred', weight='bold', 
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax.text(0.02, 0.88, 'Blue/Turquoise: Antiviral Cells', transform=ax.transAxes,
            fontsize=26, color='darkblue', weight='bold',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    ax.set_xlabel('Time (hours)', fontsize=34, fontweight='normal')
    ax.set_ylabel('Percentage of Cells (%) - Log10 Scale', fontsize=34, fontweight='normal')
    ax.set_title('IFN10 Only: Dead and Antiviral Cells Dynamics (Log10 Scale)\nComparing 1% Jump vs No Jump', 
                fontsize=38, fontweight='bold', pad=30)
    ax.tick_params(axis='both', which='major', labelsize=32)
    ax.set_ylim([1e-6, 100])
    ax.set_yticks([1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100])
    ax.set_yticklabels(['1e-6', '1e-5', '1e-4', '1e-3', '0.01', '0.1', '1', '10', '100'])
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=22, 
             frameon=True, fancybox=True, shadow=True, ncol=1)
    ax.grid(True, which='both', linestyle='--', alpha=0.3)
    ax.grid(True, which='minor', linestyle=':', alpha=0.2)
    ax.set_facecolor('#f8f8f8')
    
    plt.tight_layout()
    output_filename = 'ifn10_only_integrated_log10_plot.png'
    plt.savefig(output_filename, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {output_filename}")

def create_ifn10_integrated_dualaxis_linear():
    kjumpr_0_file = "kjumpr_0_local_dip_averaged_simulation_output.csv"
    base_dir = "kjumpr_0.01_vary_ifn_range_local_dip"
    ifn10_folder = "10_Dinit0_DIPBst50_JRand_Vinit1_VBst100_IFN10_mdbk_times500_tau12_ifnBothFold1.00_grid76_VStimulateIFNtrue"
    ifn10_file = os.path.join(base_dir, ifn10_folder, "simulation_output.csv")
    
    fig, ax_left = plt.subplots(figsize=(20, 12))
    ax_right = ax_left.twinx()
    
    try:
        df_ifn10 = pd.read_csv(ifn10_file)
        ax_left.plot(df_ifn10['Time'], df_ifn10['Percentage Dead Cells'],
            color='darkred', linewidth=10, linestyle='-', alpha=0.9,
            label='IFN10 Dead (with $\mathbf{1\%\ jump}$)')
        ax_right.plot(df_ifn10['Time'], df_ifn10['Percentage Antiviral Cells'],
            color='darkblue', linewidth=10, linestyle='-', alpha=0.9,
            label='IFN10 Antiviral (with $\mathbf{1\%\ jump}$)')
    except Exception as e:
        print(f"Error reading IFN10 file: {e}")
    
    try:
        df_baseline = pd.read_csv(kjumpr_0_file)
        ax_left.plot(df_baseline['Time'], df_baseline['Percentage Dead Cells'],
            color='darkorange', linewidth=12, linestyle='--', alpha=0.95,
            label='IFN10Dead ($\mathbf{no\ jump}$)')
        ax_right.plot(df_baseline['Time'], df_baseline['Percentage Antiviral Cells'],
            color='darkturquoise', linewidth=12, linestyle='--', alpha=0.95,
            label='IFN10Antiviral ($\mathbf{no\ jump}$)')
    except Exception as e:
        print(f"Error reading baseline file: {e}")
    
    ax_left.set_xlabel('Time (hours)', fontsize=34, fontweight='normal')
    ax_left.set_ylabel('Percentage Dead Cells (%)', fontsize=34, color='darkred')
    ax_right.set_ylabel('Percentage Antiviral Cells (%)', fontsize=34, color='darkblue', rotation=270, labelpad=20)
    
    ax_left.set_ylim(1, 20)
    ax_right.set_ylim(1, 100)
    
    ax_left.tick_params(axis='both', which='major', labelsize=32, colors='black')
    ax_right.tick_params(axis='y', labelsize=32, colors='black')
    
    ax_left.set_title('IFN10 Only: Dead (1–20%) and Antiviral (1–100%)\nDual Y-axes, Comparing 1% Jump vs No Jump',
        fontsize=38, fontweight='bold', pad=30)
    ax_left.grid(False)
    
    lines_left, labels_left = ax_left.get_legend_handles_labels()
    lines_right, labels_right = ax_right.get_legend_handles_labels()
    ax_left.legend(lines_left + lines_right, labels_left + labels_right,
        bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=22,
        frameon=True, fancybox=True, shadow=True, ncol=1)
    
    plt.tight_layout()
    output_filename = 'ifn10_only_integrated_dualaxis_linear_plot.png'
    plt.savefig(output_filename, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {output_filename}")

def create_ifn1_to_ifn10_dualaxis_linear():
    base_dir = "kjumpr_0.01_vary_ifn_range_local_dip"
    kjumpr_0_file = "kjumpr_0_local_dip_averaged_simulation_output.csv"

    ifn_files = []
    ifn_labels = []
    for i in range(1, 11):
        folder = f"{i}_Dinit0_DIPBst50_JRand_Vinit1_VBst100_IFN{i}_mdbk_times500_tau12_ifnBothFold1.00_grid76_VStimulateIFNtrue"
        csv_path = os.path.join(base_dir, folder, "simulation_output.csv")
        if os.path.exists(csv_path):
            ifn_files.append(csv_path)
            ifn_labels.append(f"IFN{i}")
        else:
            print(f"Warning: {csv_path} not found")

    fig, ax_left = plt.subplots(figsize=(20, 12))
    ax_right = ax_left.twinx()

    dead_colors = []
    anti_colors = []
    for i in range(10):
        t = i / 9.0
        # mistyrose -> red
        r = 1.0
        g = 0.894 * (1 - t)
        b = 0.882 * (1 - t)
        dead_colors.append((r, g, b))
        # alice blue -> blue
        r = 0.941 * (1 - t)
        g = 0.973 * (1 - t)
        b = 1.0
        anti_colors.append((r, g, b))

    for idx, (csv_path, label) in enumerate(zip(ifn_files, ifn_labels)):
        try:
            df = pd.read_csv(csv_path)
            ax_left.plot(df['Time'], df['Percentage Dead Cells'],
                color=dead_colors[idx], linewidth=3, linestyle='-', alpha=0.9)
            ax_right.plot(df['Time'], df['Percentage Antiviral Cells'],
                color=anti_colors[idx], linewidth=3, linestyle='-', alpha=0.9)
        except Exception as e:
            print(f"Error reading {csv_path}: {e}")

    baseline_dead_line = None
    baseline_anti_line = None
    try:
        df_baseline = pd.read_csv(kjumpr_0_file)
        baseline_dead_line = ax_left.plot(df_baseline['Time'], df_baseline['Percentage Dead Cells'],
            color='darkorange', linewidth=10, linestyle='-', alpha=0.95,
            label='IFN10Dead ($\mathbf{no\ jump}$)')[0]
        baseline_anti_line = ax_right.plot(df_baseline['Time'], df_baseline['Percentage Antiviral Cells'],
            color='darkturquoise', linewidth=10, linestyle='-', alpha=0.95,
            label='IFN10Antiviral ($\mathbf{no\ jump}$)')[0]
    except Exception as e:
        print(f"Error reading baseline file: {e}")

    ax_left.set_xlabel('Time (hours)', fontsize=34, fontweight='normal')
    ax_left.set_ylabel('Percentage Dead Cells (%)', fontsize=34, color='darkred')
    ax_right.set_ylabel('Percentage Antiviral Cells (%)', fontsize=34, color='darkblue', rotation=270, labelpad=20)

    ax_left.set_ylim(1, 30)
    ax_right.set_ylim(1, 100)

    ax_left.tick_params(axis='both', which='major', labelsize=30, colors='black')
    ax_right.tick_params(axis='y', labelsize=30, colors='black')

    ax_left.set_title('IFN1–IFN10: Dead (1–30%) and Antiviral (1–100%)\nDual Y-axes (Linear) with Baseline',
        fontsize=38, fontweight='bold', pad=26)
    ax_left.grid(False)

    from matplotlib.lines import Line2D
    legend_elements = []
    dead_gradient = Line2D([0], [0], color='red', linewidth=3, linestyle='-',
                          label='IFN range 1→10\nDead ($\mathbf{1\%\ jump}$)\n(light to dark)')
    legend_elements.append(dead_gradient)
    anti_gradient = Line2D([0], [0], color='blue', linewidth=3, linestyle='-',
                          label='IFN range 1→10\nAntiviral ($\mathbf{1\%\ jump}$)\n(light to dark)')
    legend_elements.append(anti_gradient)
    if baseline_dead_line:
        legend_elements.append(baseline_dead_line)
    if baseline_anti_line:
        legend_elements.append(baseline_anti_line)
    
    ax_left.legend(handles=legend_elements, loc='upper left', fontsize=16,
        frameon=True, fancybox=True, shadow=True, ncol=1,
        facecolor='white', edgecolor='black')

    plt.tight_layout()
    output_filename = 'ifn1_to_ifn10_integrated_dualaxis_linear_plot.png'
    plt.savefig(output_filename, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {output_filename}")

if __name__ == "__main__":
    print("Generating plots...")
    create_ifn10_integrated_linear()
    create_ifn10_integrated_log10()
    create_ifn10_integrated_dualaxis_linear()
    create_ifn1_to_ifn10_dualaxis_linear()
    print("Done")
