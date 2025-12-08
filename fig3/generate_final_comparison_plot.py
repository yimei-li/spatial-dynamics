import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

matplotlib.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Helvetica']
matplotlib.rcParams['axes.unicode_minus'] = False

N_RUNS = 30
files = []
for i in range(1, 31):
    files.append('result%d_simulation_output.csv' % i)
exp_file = '../infection_counts_by_time.csv'
times = [7, 13, 19, 25]

colors = {}
colors['virion'] = '#d62728'
colors['dip'] = '#2ca02c'
colors['both'] = '#ff7f0e'
colors['sus'] = '#000000'

def closest_row(df, target_time, col_time='Time'):
    diff = (df[col_time] - target_time).abs()
    idx = diff.idxmin()
    return df.loc[[idx]]

if os.path.exists(exp_file):
    df_exp = pd.read_csv(exp_file)
    exp = {}
    exp['time'] = times
    exp['virion'] = []
    exp['dip'] = []
    exp['both'] = []
    exp['sus'] = []
    for t in times:
        r = closest_row(df_exp, t, 'time')
        exp['virion'].append(int(r['virion_counts'].iloc[0]))
        exp['dip'].append(int(r['dip_counts'].iloc[0]))
        exp['both'].append(int(r['both_infected_counts'].iloc[0]))
        exp['sus'].append(int(r['susceptible_counts'].iloc[0]))
else:
    exp = {}
    exp['time'] = times
    exp['virion'] = [1, 183, 574, 983]
    exp['dip'] = [0, 0, 10, 35]
    exp['both'] = [0, 1, 68, 713]
    exp['sus'] = [2877, 2694, 2326, 1147]

exp_init_sus = exp['sus'][0]
print("exp initial sus at t=7: %d" % exp_init_sus)

per_file = {}
per_file['virion'] = []
per_file['dip'] = []
per_file['both'] = []
per_file['sus'] = []

for path in files:
    try:
        df = pd.read_csv(path)
    except:
        continue
    
    vals_v = []
    vals_d = []
    vals_b = []
    vals_s = []
    
    r_init = closest_row(df, times[0], 'Time')
    if len(r_init) == 0:
        continue
    gs_init = int(r_init['GRID_SIZE'].iloc[0])
    tot_init = gs_init * gs_init
    sim_init_sus = float(r_init['Percentage Susceptible Cells'].iloc[0]) * tot_init / 100.0
    
    offset = exp_init_sus - sim_init_sus
    
    for t in times:
        r = closest_row(df, t, 'Time')
        if len(r) == 0:
            vals_v.append(np.nan)
            vals_d.append(np.nan)
            vals_b.append(np.nan)
            vals_s.append(np.nan)
            continue
        gs = int(r['GRID_SIZE'].iloc[0])
        tot = gs * gs
        vir = float(r['virionOnlyInfected'].iloc[0])
        dip = float(r['dipOnlyInfected'].iloc[0])
        both_val = float(r['bothInfected'].iloc[0])
        sus = float(r['Percentage Susceptible Cells'].iloc[0]) * tot / 100.0
        
        sus_norm = sus + offset
        
        vals_v.append(vir)
        vals_d.append(dip)
        vals_b.append(both_val)
        vals_s.append(sus_norm)
    
    per_file['virion'].append(vals_v)
    per_file['dip'].append(vals_d)
    per_file['both'].append(vals_b)
    per_file['sus'].append(vals_s)

for k in per_file:
    if len(per_file[k]) > 0:
        per_file[k] = np.array(per_file[k], dtype=float)
    else:
        per_file[k] = np.empty((0, len(times)))

n_runs = per_file['virion'].shape[0]
print("loaded %d runs" % n_runs)

means = {}
for k in per_file:
    arr = per_file[k]
    if arr.size > 0:
        m = np.nanmean(arr, axis=0)
        means[k] = list(m)
    else:
        means[k] = []
        for _ in times:
            means[k].append(np.nan)

ci_lower = {}
ci_upper = {}
for k in per_file:
    arr = per_file[k]
    if arr.size == 0:
        ci_lower[k] = np.full(len(times), np.nan)
        ci_upper[k] = np.full(len(times), np.nan)
    else:
        lower = []
        upper = []
        for j in range(len(times)):
            col_data = arr[:, j]
            col_data = col_data[np.isfinite(col_data)]
            if col_data.size == 0:
                lower.append(np.nan)
                upper.append(np.nan)
            else:
                lower.append(np.percentile(col_data, 2.5))
                upper.append(np.percentile(col_data, 97.5))
        ci_lower[k] = np.array(lower)
        ci_upper[k] = np.array(upper)

fig = plt.figure(figsize=(12, 10))
fig.patch.set_facecolor('white')

ax1 = plt.subplot(2, 2, 1)
ax2 = plt.subplot(2, 2, 2)
ax3 = plt.subplot(2, 2, 3)
ax4 = plt.subplot(2, 2, 4)

ax1.set_facecolor('white')
ax2.set_facecolor('white')
ax3.set_facecolor('white')
ax4.set_facecolor('white')

virion_color = colors['virion']
virion_exp = exp['virion']
virion_sim = means['virion']
virion_lo = ci_lower['virion']
virion_hi = ci_upper['virion']
has_ci = False
for x in virion_lo:
    if not np.isnan(x):
        has_ci = True
        break
if not has_ci:
    for x in virion_hi:
        if not np.isnan(x):
            has_ci = True
            break
if has_ci:
    ax1.fill_between(times, virion_lo, virion_hi, color='#cccccc', alpha=0.7, label='95% Interval')
ax1.plot(times, virion_exp, color=virion_color, linewidth=2.5, linestyle='-', marker='o', markersize=8, label='Experimental Data')
ax1.plot(times, virion_sim, color=virion_color, linewidth=2.5, linestyle='--', marker='s', markersize=8, label='Average Results')
ax1.set_title('Virion-Infected Cells', fontsize=20, fontweight='bold')
ax1.set_ylabel('Virion-Infected Cells Count', fontsize=18)
ax1.legend(fontsize=14, loc='upper left')
ax1.tick_params(axis='x', labelsize=16)
ax1.tick_params(axis='y', labelsize=16)
ax1.grid(False)
ax1.set_xticks(times)
ax1.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))

dip_color = colors['dip']
dip_exp = exp['dip']
dip_sim = means['dip']
dip_lower = ci_lower['dip']
dip_upper = ci_upper['dip']
has_valid_ci = False
for i in range(len(dip_lower)):
    if not np.isnan(dip_lower[i]) or not np.isnan(dip_upper[i]):
        has_valid_ci = True
        break
if has_valid_ci:
    ax2.fill_between(times, dip_lower, dip_upper, color='#cccccc', alpha=0.7, label='95% Interval')
ax2.plot(times, dip_exp, c=dip_color, lw=2.5, ls='-', marker='o', ms=8, label='Experimental Data')
ax2.plot(times, dip_sim, c=dip_color, lw=2.5, ls='--', marker='s', ms=8, label='Average Results')
ax2.set_title('DIP-Infected Cells', fontsize=20, fontweight='bold')
ax2.set_ylabel('DIP-Infected Cells Count', fontsize=18)
ax2.legend(fontsize=14, loc='upper left')
ax2.tick_params(axis='x', labelsize=16)
ax2.tick_params(axis='y', labelsize=16)
ax2.grid(False)
ax2.set_xticks(times)
ax2.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))

both_clr = colors['both']
both_exp_data = exp['both']
both_sim_data = means['both']
both_ci_low = ci_lower['both']
both_ci_high = ci_upper['both']
skip_fill = True
if len(both_ci_low) > 0:
    for x in both_ci_low:
        if not np.isnan(x):
            skip_fill = False
            break
    if not skip_fill:
        for x in both_ci_high:
            if not np.isnan(x):
                skip_fill = False
                break
if not skip_fill:
    ax3.fill_between(times, both_ci_low, both_ci_high, color='#cccccc', alpha=0.7, label='95% Interval')
ax3.plot(times, both_exp_data, color=both_clr, linewidth=2.5, linestyle='-', marker='o', markersize=8, label='Experimental Data')
ax3.plot(times, both_sim_data, color=both_clr, linewidth=2.5, linestyle='--', marker='s', markersize=8, label='Average Results')
ax3.set_title('Dual-Infected Cells', fontsize=20, fontweight='bold')
ax3.set_xlabel('Time (hours)', fontsize=18)
ax3.set_ylabel('Dual-Infected Cells Count', fontsize=18)
ax3.legend(fontsize=14, loc='upper left')
ax3.tick_params(labelsize=16)
ax3.grid(False)
ax3.set_xticks(times)
loc3 = MaxNLocator(nbins=5, integer=True)
ax3.yaxis.set_major_locator(loc3)

sus_clr = colors['sus']
sus_exp_vals = exp['sus']
sus_sim_vals = means['sus']
sus_lo_bound = ci_lower['sus']
sus_hi_bound = ci_upper['sus']
can_plot_ci = False
if len(sus_lo_bound) == len(sus_hi_bound):
    for idx in range(len(sus_lo_bound)):
        if not (np.isnan(sus_lo_bound[idx]) and np.isnan(sus_hi_bound[idx])):
            can_plot_ci = True
            break
if can_plot_ci:
    ax4.fill_between(times, sus_lo_bound, sus_hi_bound, color='#cccccc', alpha=0.7, label='95% Interval')
ax4.plot(times, sus_exp_vals, color=sus_clr, linewidth=2.5, linestyle='-', marker='o', markersize=8, label='Experimental Data')
ax4.plot(times, sus_sim_vals, color=sus_clr, linewidth=2.5, linestyle='--', marker='s', markersize=8, label='Average Results')
ax4.set_title('Susceptible Cells', fontsize=20, fontweight='bold')
ax4.set_xlabel('Time (hours)', fontsize=18)
ax4.set_ylabel('Susceptible Cells Count', fontsize=18)
ax4.legend(fontsize=14, loc='lower left')
ax4.tick_params(axis='both', which='major', labelsize=16)
ax4.grid(False)
ax4.set_xticks(times)
ax4.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))

plt.tight_layout()

png_name = 'final_comparison_with_ci.png'
pdf_name = 'final_comparison_with_ci.pdf'
fig.savefig(png_name, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
fig.savefig(pdf_name, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
print("Saved: " + png_name)
print("Saved: " + pdf_name)
plt.close()
