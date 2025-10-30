#!/bin/bash

# Author: Yimei Li  
# Affiliation: Princeton University, Grenfell Lab / teVelthuis Lab / Levin Lab  
# Year: 2024  
# Copyright: © 2024 Yimei Li. All rights reserved.  
# License: [Specify if applicable, e.g., MIT, GPL, or "Proprietary"]  

echo "Hello, World!"
SOURCE_FILE="mdbk_small_vero_0502.go"
BINARY_FILE="mdbk_small_vero_0502"

# 如果已有可执行文件，则删除
if [ -f "${BINARY_FILE}" ]; then
    echo "Executable ${BINARY_FILE} found. Removing..."
    rm -f "${BINARY_FILE}"
fi

# 编译 Go 程序
echo "Building ${BINARY_FILE}..."
go build -o "${BINARY_FILE}" "${SOURCE_FILE}"
if [ $? -ne 0 ]; then
    echo "Build failed! Exiting."
    exit 1
fi
echo "Build successful: ${BINARY_FILE}"

# 固定参数 - IFNclr3 情况
meanLysisTime=12
burstSizeV=50
tau=95 
kJumpR=NaN
ifnBothFold=1.0
rho=0.026
virion_half_life=0.0
dip_half_life=0.0
ifn_half_life=3.0  # IFNclr3 使用 3.0
option=1
videotype="states"  # 视频类型

# 模拟选项
particleSpreadOption="celltocell"
ifnSpreadOption="global"
dipOption=true

# burstSizeD 列表 - 从100到1400，步长100
burstSizeD_list=($(seq 100 100 1400))

# 每个burstSizeD运行的次数
RUNS_PER_BURST_SIZE=30

# 当前脚本文件名
this_script_name="$(basename "$0")"

# ==== 一开始，记录现有的文件/文件夹 ====
initial_items=($(ls))

# 创建结果汇总文件
summary_file="IFNclr3_30runs_summary.csv"
echo "BurstSize,End_IFN_Concentration,Max_IFN_Concentration,Order" > "${summary_file}"

# 循环每个 burstSizeD
for burstSizeD in "${burstSizeD_list[@]}"; do
    echo "-------------------------------------------"
    echo "Running ${RUNS_PER_BURST_SIZE} simulations for burstSizeD=${burstSizeD}"
    echo "  TAU                   = ${tau}"
    echo "  burstSizeV            = ${burstSizeV}"
    echo "  burstSizeD            = ${burstSizeD}"
    echo "  meanLysisTime         = ${meanLysisTime}"
    echo "  kJumpR                = ${kJumpR}"
    echo "  ifnBothFold           = ${ifnBothFold}"
    echo "  rho                   = ${rho}"
    echo "  virion_half_life      = ${virion_half_life}"
    echo "  dip_half_life         = ${dip_half_life}"
    echo "  ifn_half_life         = ${ifn_half_life}"
    echo "  option                = ${option}"
    echo "  particleSpreadOption  = ${particleSpreadOption}"
    echo "  ifnSpreadOption       = ${ifnSpreadOption}"
    echo "  dipOption             = ${dipOption}"
    echo "  videotype             = ${videotype}"
    echo "-------------------------------------------"

    # 存储这个burstSizeD的所有结果
    end_ifn_values=()
    max_ifn_values=()

    # 运行30次仿真
    for run in $(seq 1 ${RUNS_PER_BURST_SIZE}); do
        echo "  Run ${run}/${RUNS_PER_BURST_SIZE} for burstSizeD=${burstSizeD}"

        # 执行仿真
        ./${BINARY_FILE} \
            --tau="${tau}" \
            --burstSizeV="${burstSizeV}" \
            --burstSizeD="${burstSizeD}" \
            --meanLysisTime="${meanLysisTime}" \
            --kJumpR="${kJumpR}" \
            --ifnBothFold="${ifnBothFold}" \
            --rho="${rho}" \
            --virion_half_life="${virion_half_life}" \
            --dip_half_life="${dip_half_life}" \
            --ifn_half_life="${ifn_half_life}" \
            --particleSpreadOption="${particleSpreadOption}" \
            --ifnSpreadOption="${ifnSpreadOption}" \
            --dipOption="${dipOption}" \
            --option="${option}" \
            --videotype="${videotype}" 

        # 找到最新生成的输出文件夹
        latest_folder=$(ls -td output_* 2>/dev/null | head -1)
        
        if [ -n "${latest_folder}" ] && [ -d "${latest_folder}" ]; then
            # 读取simulation_output.csv的最后一行来获取最终IFN浓度
            if [ -f "${latest_folder}/simulation_output.csv" ]; then
                # 获取最后一行的Global IFN Concentration Per Cell和max_global_IFN
                last_line=$(tail -1 "${latest_folder}/simulation_output.csv")
                end_ifn=$(echo "${last_line}" | cut -d',' -f5)
                max_ifn=$(echo "${last_line}" | cut -d',' -f23)
                
                # 添加到数组
                end_ifn_values+=("${end_ifn}")
                max_ifn_values+=("${max_ifn}")
                
                echo "    Run ${run}: End IFN = ${end_ifn}, Max IFN = ${max_ifn}"
            fi
        fi
    done

    # 计算平均值
    if [ ${#end_ifn_values[@]} -gt 0 ]; then
        # 计算End IFN平均值
        end_ifn_sum=0
        for val in "${end_ifn_values[@]}"; do
            end_ifn_sum=$(echo "${end_ifn_sum} + ${val}" | bc -l)
        done
        end_ifn_avg=$(echo "scale=6; ${end_ifn_sum} / ${#end_ifn_values[@]}" | bc -l)

        # 计算Max IFN平均值
        max_ifn_sum=0
        for val in "${max_ifn_values[@]}"; do
            max_ifn_sum=$(echo "${max_ifn_sum} + ${val}" | bc -l)
        done
        max_ifn_avg=$(echo "scale=6; ${max_ifn_sum} / ${#max_ifn_values[@]}" | bc -l)

        # 写入汇总文件
        echo "${burstSizeD}.0,${end_ifn_avg},${max_ifn_avg},0" >> "${summary_file}"
        
        echo "✅ Average for burstSizeD=${burstSizeD}: End IFN = ${end_ifn_avg}, Max IFN = ${max_ifn_avg}"
    else
        echo "❌ No valid results for burstSizeD=${burstSizeD}"
    fi
done

echo "🎯 All simulations completed."

# ==== 结束后，找出新生成的文件夹 ====
final_items=($(ls))

# 比较 initial_items 和 final_items
generated_folders=()
for item in "${final_items[@]}"; do
    if [[ ! " ${initial_items[@]} " =~ " ${item} " ]]; then
        if [ -d "$item" ]; then
            generated_folders+=("$item")
        fi
    fi
done

# ==== 新建总文件夹并移动 ====

# 自动找 IFNclr3_30runs 名字
base_folder="IFNclr3_30runs_${ifnSpreadOption}_${particleSpreadOption}_tau${tau}_option${option}"

final_folder="${base_folder}"
count=1
while [ -d "${final_folder}" ]; do
    count=$((count + 1))
    final_folder="${base_folder}_${count}"
done

# 创建最终总文件夹
mkdir "${final_folder}"
echo "📂 Created folder: ${final_folder}"

# 移动所有仿真输出
for folder in "${generated_folders[@]}"; do
    mv "${folder}" "${final_folder}/"
    echo "Moved ${folder} -> ${final_folder}/"
done

# 把这份 Bash 脚本复制进去
cp "$0" "${final_folder}/${this_script_name}"
echo "📄 Copied this script (${this_script_name}) into ${final_folder}/"

# 把 SOURCE_FILE（Go 源文件）也复制进去
cp "${SOURCE_FILE}" "${final_folder}/${SOURCE_FILE}"
echo "📄 Copied source file (${SOURCE_FILE}) into ${final_folder}/"

# 把汇总文件也复制进去
cp "${summary_file}" "${final_folder}/${summary_file}"
echo "📄 Copied summary file (${summary_file}) into ${final_folder}/"

echo "🚀 All output folders + script + source file + summary have been moved and saved into ${final_folder}."
echo "📊 Summary file: ${final_folder}/${summary_file}" 