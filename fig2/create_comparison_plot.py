import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import sys
import os
matplotlib.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'SimHei', 'DejaVu Sans']
matplotlib.rcParams['axes.unicode_minus'] = False

# 获取simulation文件夹路径（如果作为参数提供）
simulation_folder = sys.argv[1] if len(sys.argv) > 1 else '.'
simulation_csv_path = os.path.join(simulation_folder, 'simulation_output.csv')

# 读取simulation数据 (从指定文件夹的simulation_output.csv文件)
try:
    df_sim_output = pd.read_csv(simulation_csv_path)
    print("✅ Successfully loaded simulation data from simulation_output.csv")
    
    # 从simulation_output.csv数据创建simulation数据
    data_simulation = {
        'time': [7, 13, 19, 25],  # 固定的时间点
        'virion_counts': [0, 0, 0, 0],  # 初始化
        'dip_counts': [0, 0, 0, 0],  # 初始化
        'both_infected_counts': [0, 0, 0, 0],  # 初始化
        'susceptible_counts': [0, 0, 0, 0]  # 初始化
    }
    
    # 从simulation_output.csv数据中提取对应时间点的数据
    for i, target_time in enumerate([7, 13, 19, 25]):
        # 找到最接近的时间点
        closest_row = df_sim_output.iloc[(df_sim_output['Time'] - target_time).abs().argsort()[:1]]
        if not closest_row.empty:
            # 直接读取实际细胞数 (这些列已经是实际细胞数，不是百分比)
            virion_count = closest_row['virionOnlyInfected'].iloc[0]
            dip_count = closest_row['dipOnlyInfected'].iloc[0]
            both_count = closest_row['bothInfected'].iloc[0]
            susceptible_percent = closest_row['Percentage Susceptible Cells'].iloc[0]
            
            # 计算susceptible细胞数 (从百分比转换)
            grid_size = closest_row['GRID_SIZE'].iloc[0]  # 从CSV文件读取GRID_SIZE
            total_cells = grid_size * grid_size
            susceptible_count = int(susceptible_percent * total_cells / 100)
            
            data_simulation['virion_counts'][i] = int(virion_count)
            data_simulation['dip_counts'][i] = int(dip_count)
            data_simulation['both_infected_counts'][i] = int(both_count)
            data_simulation['susceptible_counts'][i] = susceptible_count
    
    print("📊 Simulation data loaded:")
    print(f"   Time points: {data_simulation['time']}")
    print(f"   Virion counts: {data_simulation['virion_counts']}")
    print(f"   DIP counts: {data_simulation['dip_counts']}")
    print(f"   Both infected counts: {data_simulation['both_infected_counts']}")
    print(f"   Susceptible counts: {data_simulation['susceptible_counts']}")
    
except FileNotFoundError:
    print("❌ Error: simulation_output.csv not found.")
    print("Please run the simulation first to generate simulation_output.csv file.")
    print("Or copy simulation_output.csv from a specific simulation folder.")
    exit(1)
except Exception as e:
    print(f"❌ Error reading simulation_output.csv: {e}")
    print("Please check if simulation_output.csv exists and is readable.")
    exit(1)

# 读取实验数据 (从CSV文件)
try:
    df_csv = pd.read_csv('infection_counts_by_time.csv')
    print("✅ Successfully loaded experimental data from infection_counts_by_time.csv")
    
    # 从CSV数据创建实验数据
    data_experimental = {
        'time': [7, 13, 19, 25],  # 固定的时间点
        'virion_counts': [0, 0, 0, 0],  # 初始化
        'dip_counts': [0, 0, 0, 0],  # 初始化
        'both_infected_counts': [0, 0, 0, 0],  # 初始化
        'susceptible_counts': [0, 0, 0, 0]  # 初始化
    }
    
    # 从CSV数据中提取对应时间点的数据
    for i, target_time in enumerate([7, 13, 19, 25]):
        # 找到最接近的时间点
        closest_row = df_csv.iloc[(df_csv['time'] - target_time).abs().argsort()[:1]]
        if not closest_row.empty:
            data_experimental['virion_counts'][i] = closest_row['virion_counts'].iloc[0]
            data_experimental['dip_counts'][i] = closest_row['dip_counts'].iloc[0]
            data_experimental['both_infected_counts'][i] = closest_row['both_infected_counts'].iloc[0]
            # 实验数据的susceptible counts使用simulation的GRID_SIZE计算
            # 先读取simulation数据来获取GRID_SIZE
            sim_closest_row = df_sim_output.iloc[(df_sim_output['Time'] - target_time).abs().argsort()[:1]]
            if not sim_closest_row.empty:
                grid_size = sim_closest_row['GRID_SIZE'].iloc[0]
                total_cells = grid_size * grid_size
                # 使用公式: GRID_SIZE*GRID_SIZE - (total_cells - susceptible_counts)
                # 其中total_cells是实验数据中的total_cells，susceptible_counts是实验数据中的susceptible_counts
                exp_total_cells = closest_row['total_cells'].iloc[0]
                exp_susceptible = closest_row['susceptible_counts'].iloc[0]
                data_experimental['susceptible_counts'][i] = total_cells - (exp_total_cells - exp_susceptible)
            else:
                # 如果找不到对应的simulation数据，使用原始值
                data_experimental['susceptible_counts'][i] = closest_row['susceptible_counts'].iloc[0]
    
    print("📊 Experimental data loaded:")
    print(f"   Time points: {data_experimental['time']}")
    print(f"   Virion counts: {data_experimental['virion_counts']}")
    print(f"   DIP counts: {data_experimental['dip_counts']}")
    print(f"   Both infected counts: {data_experimental['both_infected_counts']}")
    print(f"   Susceptible counts: {data_experimental['susceptible_counts']}")
    
except FileNotFoundError:
    print("❌ Error: infection_counts_by_time.csv not found. Using default hardcoded experimental data.")
    data_experimental = {
        'time': [7, 13, 19, 25],
        'virion_counts': [0, 3, 253, 2033],
        'dip_counts': [0, 0, 0, 21],
        'both_infected_counts': [0, 0, 0, 55],
        'susceptible_counts': [7020, 7017, 6767, 4911]
    }
except Exception as e:
    print(f"❌ Error reading infection_counts_by_time.csv: {e}. Using default hardcoded experimental data.")
    data_experimental = {
        'time': [7, 13, 19, 25],
        'virion_counts': [0, 3, 253, 2033],
        'dip_counts': [0, 0, 0, 21],
        'both_infected_counts': [0, 0, 0, 55],
        'susceptible_counts': [7020, 7017, 6767, 4911]
    }

# 创建DataFrame
df_experimental = pd.DataFrame(data_experimental)
df_simulation = pd.DataFrame(data_simulation)

# 计算log值 (避免log(0)，使用log(1)代替)
def safe_log(x):
    return np.log(np.maximum(x, 1))

# 应用log转换
for col in ['virion_counts', 'dip_counts', 'both_infected_counts', 'susceptible_counts']:
    df_experimental[f'{col}_log'] = safe_log(df_experimental[col])
    df_simulation[f'{col}_log'] = safe_log(df_simulation[col])

# 创建对比图 - Log版本
fig_log, axes_log = plt.subplots(2, 2, figsize=(15, 12))
fig_log.suptitle('Experimental Data vs Simulation Results - Log Scale Comparison', fontsize=16, fontweight='bold')

# 创建对比图 - 非Log版本
fig_linear, axes_linear = plt.subplots(2, 2, figsize=(15, 12))
fig_linear.suptitle('Experimental Data vs Simulation Results - Linear Scale Comparison', fontsize=16, fontweight='bold')

# 颜色设置 - 按照用户要求
colors = ['#d62728', '#ff7f0e', '#2ca02c', '#000000']  # 红色, 深黄色, 绿色, 黑色
line_styles = ['-', '--']
labels = ['Experimental Data', 'Simulation Results']

# 确保虚线代表simulation
solid_line = '-'
dashed_line = '--'

# 1. Virion感染对比 - Log版本
ax1_log = axes_log[0, 0]
ax1_log.plot(df_experimental['time'], df_experimental['virion_counts_log'], 
         color=colors[0], linewidth=3, 
         label=labels[0], linestyle=solid_line)
ax1_log.plot(df_simulation['time'], df_simulation['virion_counts_log'], 
         color=colors[0], linewidth=3, 
         label=labels[1], linestyle=dashed_line)
ax1_log.set_title('Virion-Infected Cells (Log)', fontsize=24, fontweight='bold')
ax1_log.set_xlabel('Time (hours)', fontsize=22)
ax1_log.set_ylabel('Log(Cell Count)', fontsize=22)
ax1_log.legend(fontsize=20)
ax1_log.tick_params(axis='both', which='major', labelsize=20)
ax1_log.grid(True, alpha=0.3)

# 1. Virion感染对比 - 线性版本
ax1_linear = axes_linear[0, 0]
line1 = ax1_linear.plot(df_experimental['time'], df_experimental['virion_counts'], 
         color=colors[0], linewidth=3, 
         label=labels[0], linestyle=solid_line, marker='o', markersize=8)
line2 = ax1_linear.plot(df_simulation['time'], df_simulation['virion_counts'], 
         color=colors[0], linewidth=3, 
         label=labels[1], linestyle=dashed_line, marker='s', markersize=8)

# 添加数值标注
for i, (x, y1, y2) in enumerate(zip(df_experimental['time'], df_experimental['virion_counts'], df_simulation['virion_counts'])):
    ax1_linear.annotate(f'{int(y1)}', (x, y1), textcoords="offset points", xytext=(0,10), ha='center', fontsize=14, fontweight='bold')
    ax1_linear.annotate(f'{int(y2)}', (x, y2), textcoords="offset points", xytext=(0,-15), ha='center', fontsize=14, fontweight='bold')

ax1_linear.set_title('Virion-Infected Cells (Linear)', fontsize=24, fontweight='bold')
ax1_linear.set_xlabel('Time (hours)', fontsize=22)
ax1_linear.set_ylabel('Cell Count', fontsize=22)
ax1_linear.legend(fontsize=20)
ax1_linear.tick_params(axis='both', which='major', labelsize=20)
ax1_linear.grid(True, alpha=0.3)

# 2. DIP感染对比 - Log版本
ax2_log = axes_log[0, 1]
ax2_log.plot(df_experimental['time'], df_experimental['dip_counts_log'], 
         color=colors[2], linewidth=3, 
         label=labels[0], linestyle=solid_line)
ax2_log.plot(df_simulation['time'], df_simulation['dip_counts_log'], 
         color=colors[2], linewidth=3, 
         label=labels[1], linestyle=dashed_line)
ax2_log.set_title('DIP-Infected Cells (Log)', fontsize=24, fontweight='bold')
ax2_log.set_xlabel('Time (hours)', fontsize=22)
ax2_log.set_ylabel('Log(Cell Count)', fontsize=22)
ax2_log.legend(fontsize=20)
ax2_log.tick_params(axis='both', which='major', labelsize=20)
ax2_log.grid(True, alpha=0.3)

# 2. DIP感染对比 - 线性版本
ax2_linear = axes_linear[0, 1]
line1 = ax2_linear.plot(df_experimental['time'], df_experimental['dip_counts'], 
         color=colors[2], linewidth=3, 
         label=labels[0], linestyle=solid_line, marker='o', markersize=8)
line2 = ax2_linear.plot(df_simulation['time'], df_simulation['dip_counts'], 
         color=colors[2], linewidth=3, 
         label=labels[1], linestyle=dashed_line, marker='s', markersize=8)

# 添加数值标注
for i, (x, y1, y2) in enumerate(zip(df_experimental['time'], df_experimental['dip_counts'], df_simulation['dip_counts'])):
    ax2_linear.annotate(f'{int(y1)}', (x, y1), textcoords="offset points", xytext=(0,10), ha='center', fontsize=14, fontweight='bold')
    ax2_linear.annotate(f'{int(y2)}', (x, y2), textcoords="offset points", xytext=(0,-15), ha='center', fontsize=14, fontweight='bold')

ax2_linear.set_title('DIP-Infected Cells (Linear)', fontsize=24, fontweight='bold')
ax2_linear.set_xlabel('Time (hours)', fontsize=22)
ax2_linear.set_ylabel('Cell Count', fontsize=22)
ax2_linear.legend(fontsize=20)
ax2_linear.tick_params(axis='both', which='major', labelsize=20)
ax2_linear.grid(True, alpha=0.3)

# 3. 双重感染对比 - Log版本
ax3_log = axes_log[1, 0]
ax3_log.plot(df_experimental['time'], df_experimental['both_infected_counts_log'], 
         color=colors[1], linewidth=3, 
         label=labels[0], linestyle=solid_line)
ax3_log.plot(df_simulation['time'], df_simulation['both_infected_counts_log'], 
         color=colors[1], linewidth=3, 
         label=labels[1], linestyle=dashed_line)
ax3_log.set_title('Dual-Infected Cells (Log)', fontsize=24, fontweight='bold')
ax3_log.set_xlabel('Time (hours)', fontsize=22)
ax3_log.set_ylabel('Log(Cell Count)', fontsize=22)
ax3_log.legend(fontsize=20)
ax3_log.tick_params(axis='both', which='major', labelsize=20)
ax3_log.grid(True, alpha=0.3)

# 3. 双重感染对比 - 线性版本
ax3_linear = axes_linear[1, 0]
line1 = ax3_linear.plot(df_experimental['time'], df_experimental['both_infected_counts'], 
         color=colors[1], linewidth=3, 
         label=labels[0], linestyle=solid_line, marker='o', markersize=8)
line2 = ax3_linear.plot(df_simulation['time'], df_simulation['both_infected_counts'], 
         color=colors[1], linewidth=3, 
         label=labels[1], linestyle=dashed_line, marker='s', markersize=8)

# 添加数值标注
for i, (x, y1, y2) in enumerate(zip(df_experimental['time'], df_experimental['both_infected_counts'], df_simulation['both_infected_counts'])):
    ax3_linear.annotate(f'{int(y1)}', (x, y1), textcoords="offset points", xytext=(0,10), ha='center', fontsize=14, fontweight='bold')
    ax3_linear.annotate(f'{int(y2)}', (x, y2), textcoords="offset points", xytext=(0,-15), ha='center', fontsize=14, fontweight='bold')

ax3_linear.set_title('Dual-Infected Cells (Linear)', fontsize=24, fontweight='bold')
ax3_linear.set_xlabel('Time (hours)', fontsize=22)
ax3_linear.set_ylabel('Cell Count', fontsize=22)
ax3_linear.legend(fontsize=20)
ax3_linear.tick_params(axis='both', which='major', labelsize=20)
ax3_linear.grid(True, alpha=0.3)

# 4. 易感细胞对比 - Log版本
ax4_log = axes_log[1, 1]
ax4_log.plot(df_experimental['time'], df_experimental['susceptible_counts_log'], 
         color=colors[3], linewidth=3, 
         label=labels[0], linestyle=solid_line)
ax4_log.plot(df_simulation['time'], df_simulation['susceptible_counts_log'], 
         color=colors[3], linewidth=3, 
         label=labels[1], linestyle=dashed_line)
ax4_log.set_title('Susceptible Cells (Log)', fontsize=24, fontweight='bold')
ax4_log.set_xlabel('Time (hours)', fontsize=22)
ax4_log.set_ylabel('Log(Cell Count)', fontsize=22)
ax4_log.legend(fontsize=20)
ax4_log.tick_params(axis='both', which='major', labelsize=20)
ax4_log.grid(True, alpha=0.3)

# 4. 易感细胞对比 - 线性版本
ax4_linear = axes_linear[1, 1]
line1 = ax4_linear.plot(df_experimental['time'], df_experimental['susceptible_counts'], 
         color=colors[3], linewidth=3, 
         label=labels[0], linestyle=solid_line, marker='o', markersize=8)
line2 = ax4_linear.plot(df_simulation['time'], df_simulation['susceptible_counts'], 
         color=colors[3], linewidth=3, 
         label=labels[1], linestyle=dashed_line, marker='s', markersize=8)

# 添加数值标注
for i, (x, y1, y2) in enumerate(zip(df_experimental['time'], df_experimental['susceptible_counts'], df_simulation['susceptible_counts'])):
    ax4_linear.annotate(f'{int(y1)}', (x, y1), textcoords="offset points", xytext=(0,10), ha='center', fontsize=14, fontweight='bold')
    ax4_linear.annotate(f'{int(y2)}', (x, y2), textcoords="offset points", xytext=(0,-15), ha='center', fontsize=14, fontweight='bold')

ax4_linear.set_title('Susceptible Cells (Linear)', fontsize=24, fontweight='bold')
ax4_linear.set_xlabel('Time (hours)', fontsize=22)
ax4_linear.set_ylabel('Cell Count', fontsize=22)
ax4_linear.legend(fontsize=20)
ax4_linear.tick_params(axis='both', which='major', labelsize=20)
ax4_linear.grid(True, alpha=0.3)

# 设置x轴刻度
x_ticks = [7, 13, 19, 25]
for ax in axes_log.flat:
    ax.set_xticks(x_ticks)
    ax.set_xticklabels([str(x) for x in x_ticks])
for ax in axes_linear.flat:
    ax.set_xticks(x_ticks)
    ax.set_xticklabels([str(x) for x in x_ticks])

# 调整布局
fig_log.tight_layout()
fig_linear.tight_layout()

# 保存图片到指定的输出文件夹
output_log_png = os.path.join(simulation_folder, 'comparison_plot_log.png')
output_log_pdf = os.path.join(simulation_folder, 'comparison_plot_log.pdf')
output_linear_png = os.path.join(simulation_folder, 'comparison_plot_linear.png')
output_linear_pdf = os.path.join(simulation_folder, 'comparison_plot_linear.pdf')

fig_log.savefig(output_log_png, dpi=300, bbox_inches='tight')
fig_log.savefig(output_log_pdf, bbox_inches='tight')
fig_linear.savefig(output_linear_png, dpi=300, bbox_inches='tight')
fig_linear.savefig(output_linear_pdf, bbox_inches='tight')

# 打印数值对比
print("=== Numerical Comparison (Log Scale) ===")
print("\nVirion-Infected Cells (Log):")
for i, t in enumerate(df_experimental['time']):
    print(f"Time {t}h: Experimental={df_experimental['virion_counts_log'].iloc[i]:.2f}, Simulation={df_simulation['virion_counts_log'].iloc[i]:.2f}")

print("\nDIP-Infected Cells (Log):")
for i, t in enumerate(df_experimental['time']):
    print(f"Time {t}h: Experimental={df_experimental['dip_counts_log'].iloc[i]:.2f}, Simulation={df_simulation['dip_counts_log'].iloc[i]:.2f}")

print("\nDual-Infected Cells (Log):")
for i, t in enumerate(df_experimental['time']):
    print(f"Time {t}h: Experimental={df_experimental['both_infected_counts_log'].iloc[i]:.2f}, Simulation={df_simulation['both_infected_counts_log'].iloc[i]:.2f}")

print("\nSusceptible Cells (Log):")
for i, t in enumerate(df_experimental['time']):
    print(f"Time {t}h: Experimental={df_experimental['susceptible_counts_log'].iloc[i]:.2f}, Simulation={df_simulation['susceptible_counts_log'].iloc[i]:.2f}")

print("\n✅ Comparison plot saved as 'comparison_plot_log.png' and 'comparison_plot_log.pdf'") 