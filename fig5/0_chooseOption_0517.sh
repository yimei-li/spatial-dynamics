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
tau=95
burstSizeV=50
burstSizeD=1600
meanLysisTime=12.0 # 12.0
kJumpR=NaN
ifnBothFold=1.0
rho=0.026
virion_half_life=0.0 
dip_half_life=0.0 
ifn_half_life=0.0 
option=1  # ✅ 显式传入的初始化方式
videotype="IFNconcentration"  # 视频类型
dipSynthesisAdvantage=4

# 模拟选项
particleSpreadOption="celltocell"
ifnSpreadOption="global"
dipOption="true"


# 复制当前脚本与 Go 源码文件
cp "${SCRIPT_PATH}" "${BASE_OUTPUT_FOLDER}/$(basename "$SCRIPT_PATH")"
cp "${SOURCE_FILE}" "${BASE_OUTPUT_FOLDER}/"


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
    --videotype="${videotype}" \
    -dipSynthesisAdvantage $dipSynthesisAdvantage 

echo "✅ Test simulation finished."
