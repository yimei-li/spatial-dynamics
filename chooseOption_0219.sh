#!/bin/bash
# 脚本说明：
# 1. 如果已编译好的可执行文件存在，则删除之，并重新编译 Go 程序。
# 2. 根据选项组合生成一个子文件夹（例如 "jumprandomly_local_diptrue"），
#    然后在该子文件夹内调用上级目录的可执行文件运行仿真，
#    仿真结果会保存在这个子文件夹内。

# 设置 Go 源代码文件名和生成的可执行文件名
SOURCE_FILE="mdbk_small_vero_0219.go"
BINARY_FILE="mdbk_small_vero_0219"

# 如果已有可执行文件，则先删除
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

# 固定参数设置
tau=12
burstSizeV=50
burstSizeD=100
meanLysisTime=12.0
kJumpR=0.5
ifnBothFold=1.0
rho=0.026
percentRandJumpR=1.0
virion_half_life=0.0 #3.2
dip_half_life=0.0 # 3.2
ifn_half_life=10.0
# 选项参数：这里测试某一特定组合
# 可能的 particleSpreadOption: "celltocell"、"jumprandomly"、"jumpradius"、"partition"
# 可能的 ifnSpreadOption: "global"、"local"、"noIFN"
# 可能的 dipOption: "true" 或 "false"
particleSpreadOption="jumprandomly"
ifnSpreadOption="local"
dipOption="true"


# 根据选项组合构造基础子文件夹名称
BASE_FOLDER="${particleSpreadOption}_${ifnSpreadOption}_dip${dipOption}_ifnclr${ifn_half_life}_vclr${virion_half_life}_dclr${dip_half_life}"

# 如果文件夹已存在，则增加序号
FOLDER="${BASE_FOLDER}"
COUNT=1
while [ -d "${FOLDER}" ]; do
    FOLDER="${BASE_FOLDER}_${COUNT}"
    ((COUNT++))
done


mkdir -p "${FOLDER}"
echo "Output will be saved in subfolder: ${FOLDER}"

echo "-------------------------------------------"
echo "Running test simulation with the following parameters:"
echo "  TAU                   = ${tau}"
echo "  burstSizeV            = ${burstSizeV}"
echo "  burstSizeD            = ${burstSizeD}"
echo "  meanLysisTime         = ${meanLysisTime}"
echo "  kJumpR                = ${kJumpR}"
echo "  ifnBothFold           = ${ifnBothFold}"
echo "  rho                   = ${rho}"
echo "  percentRandJumpR      = ${percentRandJumpR}"
echo "  particleSpreadOption  = ${particleSpreadOption}"
echo "  ifnSpreadOption       = ${ifnSpreadOption}"
echo "  dipOption             = ${dipOption}"
echo "-------------------------------------------"

# 切换到子文件夹中运行仿真（注意可执行文件在上级目录）
pushd "${FOLDER}" > /dev/null
../${BINARY_FILE} \
    --tau="${tau}" \
    --burstSizeV="${burstSizeV}" \
    --burstSizeD="${burstSizeD}" \
    --meanLysisTime="${meanLysisTime}" \
    --kJumpR="${kJumpR}" \
    --ifnBothFold="${ifnBothFold}" \
    --rho="${rho}" \
    --particleSpreadOption="${particleSpreadOption}" \
    --percentRandJumpR="${percentRandJumpR}" \
    --ifnSpreadOption="${ifnSpreadOption}" \
    --dipOption="${dipOption}" \
    --virion_half_life="${virion_half_life}" \
    --dip_half_life="${dip_half_life}" \
    --ifn_half_life="${ifn_half_life}"
popd > /dev/null




echo "Test simulation finished. Results are saved in subfolder: ${FOLDER}"
