#!/bin/bash
# 脚本说明：如果已经编译好了可执行文件，则直接运行，否则先编译。
# 然后循环运行仿真，传入不同的 TAU 参数以及其他选项。

# 设置 Go 源代码文件名和生成的可执行文件名
SOURCE_FILE="mdbk_small_0218.go"
BINARY_FILE="mdbk_small_0218"

# 测试是否存在可执行文件
if [ -f "${BINARY_FILE}" ]; then
    echo "Executable ${BINARY_FILE} found. Removing..."
    rm -f "${BINARY_FILE}"
fi

# 重新构建
echo "Building ${BINARY_FILE}..."
go build -o "${BINARY_FILE}" "${SOURCE_FILE}"
if [ $? -ne 0 ]; then
    echo "Build failed! Exiting."
    exit 1
fi

echo "Build successful: ${BINARY_FILE}"

# 定义需要循环的 TAU 值数组
taus=(12) # default 95

# 固定参数，可以根据需要修改
burstSizeV=50
burstSizeD=100
meanLysisTime=12.0
kJumpR=0.5
ifnBothFold=1.0
rho=0.026

# 新增选项参数：
# particleSpreadOption 可选值："celltocell"、"jumprandomly"、"jumpradius"、"partition"
particleSpreadOption="jumprandomly"
# ifnSpreadOption 可选值："global"、"local"、"noIFN"
ifnSpreadOption="noIFN"
# dipOption 可选值：true 或 false
dipOption=true
# 如果选择随机跳跃，可设置百分比
percentRandJumpR=1.0

# 循环运行仿真，每次传入不同的 TAU 值
for tau in "${taus[@]}"; do
    echo "-------------------------------------------"
    echo "Running simulation with TAU = ${tau}"
    ./${BINARY_FILE} \
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
        --dipOption="${dipOption}"
    echo "Executing: ./${BINARY_FILE} --ifnSpreadOption=${ifnSpreadOption}"

    echo "Simulation with TAU = ${tau} finished."
done

echo "All simulations completed."
