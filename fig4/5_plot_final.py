import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
import matplotlib.lines as mlines

data = pd.read_csv('./simulation_output.csv')

max_ifn = data['max_global_IFN'].iloc[0]
v_init = data['v_pfu_initial'].iloc[0]
d_init = data['d_pfu_initial'].iloc[0]
rho_val = data['RHO'].iloc[0]
burst_size = data['BURST_SIZE'].iloc[0]
dip_burst = data['DIP_BURST_PCT'].iloc[0]

if max_ifn == -1:
    prefix = "Vero"
else:
    prefix = "MDBK"

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(1, 1, 1)

vero_exp = [0, 0, 0.704717105910202, 1.795833997483578, 3.5070185136826266, 4.288143985511917, 3.5190749602357716]
t_pts = [0, 24, 48, 72, 96, 120, 144]

selected = data[data['Time'].isin(t_pts)]
sim_data = selected['Plaque Percentage'].values

ax.plot(t_pts, sim_data, color='red', linewidth=5, alpha=0.6)
ax.scatter(t_pts, sim_data, marker='^', color='red', s=150)
ax.scatter(t_pts, vero_exp, color='black', marker='v', s=150)

ax.set_ylabel('Dead Cell, %', fontsize=20)
ax.set_xticks(t_pts)
ax.set_xlabel('Time (hours)', fontsize=20)
ax.tick_params(labelsize=16)

sim_line = mlines.Line2D([], [], color='red', marker='^', markersize=15, label='Simulation Result', linewidth=4)
exp_line = mlines.Line2D([], [], color='black', marker='v', markersize=15, label='Experiment', linestyle='None')
ax.legend(handles=[sim_line, exp_line], loc='upper left', fontsize=18)

fname = "5_%s_RHO=%s_VInt=%s_DInt=%s_VBt=%s_DBt=%s.png" % (prefix, rho_val, v_init, d_init, burst_size, dip_burst)

output_dir = None
if len(sys.argv) >= 2:
    output_dir = sys.argv[1]

if output_dir is not None:
    if not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir)
        except:
            pass
    output_file = os.path.join(output_dir, 'comparison_plot.png')
else:
    output_file = fname

plt.savefig(output_file, dpi=300, bbox_inches='tight')

print("Plot saved to " + output_file)
