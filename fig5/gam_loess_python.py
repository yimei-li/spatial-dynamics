import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
# package for gam
from pygam import LinearGAM, s
from statsmodels.nonparametric.smoothers_lowess import lowess


# import data
df = pd.read_csv('./ifn_peak_vs_dipburst_baseline.csv')
##################
# the raw data has all different values of burst_size_DIP, you can either just use the raw data, which in my paper what I did initially, or you can use the multiples of 50 for burst_size_DIP, which is what i did in this script. I think this makes more sense, as we care about ratio, and we want our input data to have same steps.
df = df[df['burst_size_DIP'] % 50 == 0].copy()
df['relative_DIP'] = df['burst_size_DIP'] / 50
X = df[['relative_DIP']].values
y = df['peak_IFN'].values


# use pygam, there are other packaged I tried, like mgcv, good too
# I did not use scipy.UnivariateSpline
# n_splines is k in R, the larger, the more flexible (model), but the more overfitting
# lam is sp in R, the larger, the more smooth (model), but the more underfitting
# R default: k=10, sp=auto (run get_r_params.R to get exact sp value？？？)
# To match R exactly: use n_splines=10 and lam from R's auto-selected value
# Current: n_splines=10, lam=270.127 (from R run)


########## Run this to match R exactly ################################
# gam = LinearGAM(s(0, basis='ps', n_splines=10, lam=270.127)).fit(X, y)
########## Run this to use default parameters in python ###############
gam = LinearGAM(s(0, basis='ps')).fit(X, y)




grid = gam.generate_X_grid(term=0, n=500) # 500 or more points, the smoother the line
preds = gam.predict(grid)


# find maximum of the line
optimum_idx = np.argmax(preds)
optimum_x_axis = round(grid[optimum_idx, 0], 1)
optimum_legend = f"Optimum DIP Yield"


######################################
# gam's 95% CI
ci = gam.confidence_intervals(grid, width=0.95) 
ci_lower = ci[:, 0]
ci_upper = ci[:, 1]

######################################
# loess fit
# loess_result = lowess(y, X.ravel(), frac=0.6667, return_sorted=True) # similar to R's default loess
loess_result = lowess(y, X.ravel(), return_sorted=True)
loess_x = loess_result[:, 0]
loess_y = loess_result[:, 1]

# fonts size
plt.rcParams['font.size'] = 12
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
fig, ax = plt.subplots(figsize=(6, 6), dpi=300)

ax.fill_between(grid[:, 0], ci_lower, ci_upper, color='red', alpha=0.4, label='_nolegend_')
ax.scatter(X.ravel(), y, color='#3780b4', alpha=0.2, s=10) # light blue, order above lines 


ax.plot(grid[:, 0], preds, color='red', linewidth=1.5)
ax.plot(loess_x, loess_y, color='green', linestyle='--', linewidth=1.5)
ax.axvline(x=optimum_x_axis, color='black', linestyle=':', linewidth=1)

ax.set_ylim(0, 0.5)
ax.set_xlim(X.min(), X.max())
ax.set_xlabel('Relative DIP Burst Size')
ax.set_ylabel('Peak IFN')
ax.set_facecolor('white')
ax.grid(False)

# same as theme_bw() in R
for spine in ax.spines.values():
    spine.set_color('black')
    spine.set_linewidth(0.5)

legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#3780b4', 
           markersize=8, label='Stochastic Replicates'),
    Line2D([0], [0], color='red', linewidth=1.5, label='GAM Fit'),
    Line2D([0], [0], color='green', linestyle='--', linewidth=1.5, label='LOESS Fit'),
    Line2D([0], [0], color='black', linestyle=':', linewidth=1, label=optimum_legend),
    Patch(facecolor='red', alpha=0.6, edgecolor='red', label='95% Confidence Band')
]
ax.legend(handles=legend_elements, loc='lower right', bbox_to_anchor=(0.98, 0.02),
          frameon=False, fontsize=12, handlelength=2)
plt.tight_layout()
plt.savefig('./fig5_1125_python.png', dpi=300, bbox_inches='tight', 
            facecolor='white', edgecolor='none')
