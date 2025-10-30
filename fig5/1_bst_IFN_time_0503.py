import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import re
import numpy as np

# === 设置 ===
folder_name = "IFNclr5_loop_burstSizeD_global_celltocell_tau95_option1"  # ⚠️ 只改这里
select_all = False
selected_bursts = [10] + [50] + list(range(100, 4201, 100))  # ✅ 只画 1, 100, ..., 2000

    
# === 自动信息提取 ===
tau_match = re.search(r'tau(\d+)', folder_name)
option_match = re.search(r'option(\d+)', folder_name)
tau_val = tau_match.group(1) if tau_match else "NA"
option_val = option_match.group(1) if option_match else "NA"
base_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), folder_name)

# 判断 loop 类型
if 'burstSizeV' in folder_name:
    burst_label = "BurstV"
    pattern = r'VBst(\d+)'
    cmap = cm.get_cmap('Reds')
    norm = mcolors.Normalize(vmin=10, vmax=200)
elif 'burstSizeD' in folder_name:
    burst_label = "BurstD"
    pattern = r'DIPBst(\d+)'
    cmap = cm.get_cmap('Blues')  # ✅ Matplotlib 自带，专为这个设计


    burst_values = selected_bursts

    norm = mcolors.Normalize(vmin=min(burst_values), vmax=max(burst_values))  # ✅ 实际数据范围

else:
    raise ValueError("Cannot determine loop type from folder name.")

# 找到所有子文件夹
subfolders = [f for f in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, f)) and re.match(r'^\d+_', f)]

# 提取 burst 值并排序
burst_folders = []
for folder in subfolders:
    match = re.search(pattern, folder)
    if match:
        burst_value = int(match.group(1))
        if select_all or burst_value in selected_bursts:
            burst_folders.append((burst_value, folder))

burst_folders.sort()

# 记录 IFN 数据
ifn_results = []
max_ifn_lookup = {}
for burst_value, folder in burst_folders:
    file_path = os.path.join(base_dir, folder, 'simulation_output.csv')
    if os.path.exists(file_path):
        df = pd.read_csv(file_path)
        max_ifn_lookup[burst_value] = df['Global IFN Concentration Per Cell'].max()

max_ifn_burst_value = max(max_ifn_lookup, key=max_ifn_lookup.get)

# 创建左右图布局
# 创建左右图布局
from matplotlib.gridspec import GridSpec

fig = plt.figure(figsize=(14, 6))
gs = GridSpec(1, 2, width_ratios=[1.2, 1], wspace=0.4)  # ✅ 增加 wspace 空隙
main_ax = fig.add_subplot(gs[0])
bar_ax = fig.add_subplot(gs[1])

# 左图：IFN dynamics 曲线
for burst_value, folder in burst_folders:
    file_path = os.path.join(base_dir, folder, 'simulation_output.csv')
    if os.path.exists(file_path):
        df = pd.read_csv(file_path)
        time = df['Time']
        ifn_concentration = df['Global IFN Concentration Per Cell']
        is_max = (burst_value == max_ifn_burst_value)
        color = 'red' if is_max else cmap(norm(burst_value))
        linewidth = 3.5 if is_max else (2.5 if is_max else 1.5)
        z = 15 if is_max else (10 if is_max else 1)
        alpha = 1.0

        main_ax.plot(time, ifn_concentration, label=f'{burst_label} {burst_value}',
                     color=color, linewidth=linewidth, zorder=z, alpha=alpha)

        # ✅ 动态高亮最大 IFN 的线
        # if is_max:
        #     mid_index = len(time) // 2
        #     main_ax.text(
        #         time.iloc[mid_index],
        #         ifn_concentration.iloc[mid_index],
        #         f'{burst_value} ← Max IFN',
        #         fontsize=9,
        #         va='center',
        #         ha='right',
        #         color='red',
        #         fontweight='bold',
        #         bbox=dict(facecolor='white', edgecolor='navy', boxstyle='round,pad=0.3')
        #     )

        

        ifn_results.append({
            'BurstSize': burst_value,
            'End_IFN_Concentration': ifn_concentration.iloc[-1],
            'Max_IFN_Concentration': max_ifn_lookup[burst_value]
        })
    else:
        print(f'⚠️ Warning: {file_path} not found.')

main_ax.set_xlabel('Time')
main_ax.set_ylabel('Global IFN Concentration Per Cell')
main_ax.set_title(f'IFN Dynamics (Time ≤ 500) Across {burst_label}s\n(τ={tau_val}, option={option_val})')

# ✅ 修改 legend 样式：小字体，多列，图外显示
main_ax.legend(
    fontsize=6,
    loc='upper left',
    bbox_to_anchor=(1.02, 1),
    frameon=False
)

main_ax.grid(True)

# 右图：高线图
# 右图：高线图
results_df = pd.DataFrame(ifn_results)
df_bar = results_df.copy()
max_ifn = df_bar['Max_IFN_Concentration'].max()

for x, y in zip(df_bar['BurstSize'], df_bar['Max_IFN_Concentration']):
    color = 'red' if y == max_ifn else 'black'
    bar_ax.vlines(x=x, ymin=0, ymax=y, color=color, linewidth=1.5)
    bar_ax.plot(x, y, 'o', color=color)
    # bar_ax.text(x, y + 0.05, f'{y:.2f}', ha='center', va='bottom', fontsize=9, color=color)

bar_ax.set_xticks(sorted(df_bar['BurstSize']))
bar_ax.tick_params(axis='x', labelrotation=45, labelsize=7)  # ✅ 旋转字体 + 缩小字号

bar_ax.set_ylabel('Max IFN Concentration')
bar_ax.set_title('Max IFN vs. ' + "Burst Size of DIPs", pad=17)
bar_ax.set_xticks(sorted(df_bar['BurstSize']))
bar_ax.grid(True, alpha=0.3)

# 保存合图
fig.tight_layout(rect=[0, 0, 0.9, 1])  # ← 更新这一行
combined_path = os.path.join(base_dir, f"1_maxIFN_vs_{burst_label}_{folder_name}.png")
plt.savefig(combined_path, dpi=400)
print(f"✅ Saved combined figure to: {combined_path}")
plt.show()

# 保存 CSV
results_df = results_df[['BurstSize', 'End_IFN_Concentration', 'Max_IFN_Concentration']]
results_df['Order'] = results_df['Max_IFN_Concentration'].rank(method='min', ascending=False).astype(int)
row_max = results_df.loc[results_df['Max_IFN_Concentration'].idxmax()].copy()
row_max['Order'] = "Max IFN BurstSize"
results_df = pd.concat([results_df, pd.DataFrame([row_max])], ignore_index=True)

csv_path = os.path.join(base_dir, f"maxIFN_burst_{folder_name}.csv")
results_df.to_csv(csv_path, index=False)
print(f"\n✅ Saved IFN summary to: {csv_path}")
print("\nIFN Summary per Burst Size:")
print(results_df.to_string(index=False))
