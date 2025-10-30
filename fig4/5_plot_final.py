
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
import matplotlib.lines as mlines

# Load the CSV file
file_path = './simulation_output.csv'
df = pd.read_csv(file_path)

# Extract values for the title
max_global_IFN = df['max_global_IFN'].iloc[0]
v_pfu_initial = df['v_pfu_initial'].iloc[0]
d_pfu_initial = df['d_pfu_initial'].iloc[0]
RHO= df['RHO'].iloc[0]
BURST_SIZE = df['BURST_SIZE'].iloc[0]
DIP_BURST_PCT = df['DIP_BURST_PCT'].iloc[0]

# Determine the prefix based on max_global_IFN
prefix = "Vero" if max_global_IFN == -1 else "MDBK"

# Set up the plotting grid: 1 row, 1 column
fig, ax = plt.subplots(1, 1, figsize=(8, 6))

# Data from experiment
experiment_vero_50 = [0,0, 0.704717105910202, 1.795833997483578, 3.5070185136826266, 4.288143985511917, 3.5190749602357716]

# 1. Plaque Percentage over Time
# Define time points for markers based on the given data
time_points_for_markers = [0, 24, 48, 72, 96, 120, 144]

# Filter plaque percentage values at specific time points
df_filtered = df[df['Time'].isin(time_points_for_markers)]
plaque_percentage_values_at_markers = df_filtered['Plaque Percentage'].values

# Plot the line for Simulation Plaque Percentage (including red triangle markers)
ax.plot(time_points_for_markers, plaque_percentage_values_at_markers, color='red', linewidth=5, alpha=0.6)

# Add the triangle markers at the specific time points for the Simulation Plaque Percentage
ax.scatter(time_points_for_markers, plaque_percentage_values_at_markers, marker='^', color='red', s=150)

# Add the new experiment data with 'v' marker (inverted triangle)
ax.scatter(time_points_for_markers, experiment_vero_50, color='black', marker='v', s=150)

# Set the title, y-label, and x-ticks
ax.set_ylabel('Dead Cell, %', fontsize=20)
ax.set_xticks(time_points_for_markers)  # Set x-axis ticks
ax.set_xlabel('Time (hours)', fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=16)


# Create custom legend items with both lines and markers
simulation_line = mlines.Line2D([], [], color='red', marker='^', markersize=15, label='Simulation Result', linewidth=4)
experiment_line = mlines.Line2D([], [], color='black', marker='v', markersize=15, label='Experiment', linestyle='None')


# Customize the legend to show both lines and markers
ax.legend(handles=[simulation_line, experiment_line], loc='upper left', fontsize=18)

# Create a filename using the same content, but replace \n with _ for a valid file name
filename = f"5_{prefix}_RHO={RHO}_VInt={v_pfu_initial}_DInt={d_pfu_initial}_VBt={BURST_SIZE}_DBt={DIP_BURST_PCT}.png"

# Determine output path
output_dir = None
if len(sys.argv) >= 2:
    # sys.argv[1] is expected to be an absolute or relative output directory path
    output_dir = sys.argv[1]

if output_dir is not None:
    try:
        os.makedirs(output_dir, exist_ok=True)
    except Exception:
        pass
    # Save standardized filename for Shiny integration
    output_file = os.path.join(output_dir, 'comparison_plot.png')
else:
    # Fallback: original behavior in current working directory
    output_folder = './'
    output_file = os.path.join(output_folder, filename)

plt.savefig(output_file, dpi=300, bbox_inches='tight')

# Display the plots
# Optional display (commented out for headless execution)
# plt.show()

print(f"Plot saved to {output_file}")

