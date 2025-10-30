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

# 固定参数
meanLysisTime=12
burstSizeV=50
tau=95 
kJumpR=NaN
ifnBothFold=1.0
rho=0.026
virion_half_life=0.0
dip_half_life=0.0
ifn_half_life=0.0
option=1
videotype="states"  # 视频类型

# 模拟选项
particleSpreadOption="celltocell"
ifnSpreadOption="global"
dipOption=true

# 复制当前脚本与 Go 源码文件
cp "${SCRIPT_PATH}" "${BASE_OUTPUT_FOLDER}/$(basename "$SCRIPT_PATH")"
cp "${SOURCE_FILE}" "${BASE_OUTPUT_FOLDER}/"


# burstSizeD 列表
burstSizeD_list=($(seq 0 100 4200))


# 记录生成的文件夹
generated_folders=()

# 当前脚本文件名
this_script_name="$(basename "$0")"

# ==== 一开始，记录现有的文件/文件夹 ====
initial_items=($(ls))

# 循环每个 burstSizeD
for burstSizeD in "${burstSizeD_list[@]}"; do
    echo "-------------------------------------------"
    echo "Running test simulation with the following parameters:"
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

    echo "✅ Simulation for burstSizeD=${burstSizeD} finished."

    # 记录仿真输出的文件夹名
    folder_name="output_burstSizeD${burstSizeD}"
    if [ -d "${folder_name}" ]; then
        generated_folders+=("${folder_name}")
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

# 自动找 loop_burstsizeD_optionX 名字
base_folder="loop_burstSizeD_${ifnSpreadOption}_${particleSpreadOption}_tau${tau}_option${option}"

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

echo "🚀 All output folders + script + source file have been moved and saved into ${final_folder}."

