#!/bin/bash
# 脚本说明：此脚本会循环 tau=12,24,...120，为每个 tau 生成一个独立的输出文件夹，并运行 18 种仿真组合。

# 设置 Go 源代码文件名和生成的可执行文件名
SOURCE_FILE="mdbk_small_vero_0219.go"
BINARY_FILE="mdbk_small_vero_0219"

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

# 需要循环的 tau 值
tau_values=(12 24 36 48 60 72 84 96 108 120)

# 固定参数
burstSizeV=50
burstSizeD=100
meanLysisTime=12.0
kJumpR=0.5
ifnBothFold=1.0
rho=0.026
percentRandJumpR=1.0
virion_half_life=0.0  # 3.2
dip_half_life=0.0     # 3.2
ifn_half_life=0.0     # 3.0

# 需要遍历的选项
particleSpreadOptions=("celltocell" "jumprandomly" "jumpradius")
ifnSpreadOptions=("global" "local" "noIFN")
dipOptions=("true" "false")

# 循环不同的 tau
for tau in "${tau_values[@]}"; do
    BASE_OUTPUT_FOLDER="combine18_tau${tau}"

    # 如果该文件夹已存在，则加上后缀 _1, _2, ... 直到找到一个不存在的名称
    if [ -d "${BASE_OUTPUT_FOLDER}" ]; then
        counter=1
        while [ -d "${BASE_OUTPUT_FOLDER}_${counter}" ]; do
            counter=$((counter+1))
        done
        BASE_OUTPUT_FOLDER="${BASE_OUTPUT_FOLDER}_${counter}"
    fi
    mkdir -p "${BASE_OUTPUT_FOLDER}"
    echo "Output will be saved in folder: ${BASE_OUTPUT_FOLDER}"

    # 进入输出文件夹，这样生成的所有子文件夹都会在 BASE_OUTPUT_FOLDER 内
    pushd "${BASE_OUTPUT_FOLDER}" > /dev/null

    runCount=1

    # 遍历所有组合，共 18 种可能
    for particleSpread in "${particleSpreadOptions[@]}"; do
        for ifnSpread in "${ifnSpreadOptions[@]}"; do
            for dip in "${dipOptions[@]}"; do
                # 生成子文件夹
                SUBFOLDER="${particleSpread}_${ifnSpread}_dip${dip}"
                mkdir -p "${SUBFOLDER}"

                echo "-------------------------------------------"
                echo "Run ${runCount}:"
                echo "  TAU = ${tau}"
                echo "  particleSpreadOption = ${particleSpread}"
                echo "  ifnSpreadOption = ${ifnSpread}"
                echo "  dipOption = ${dip}"
                echo "  Output subfolder: ${SUBFOLDER}"

                # 进入子文件夹并运行仿真
                pushd "${SUBFOLDER}" > /dev/null
                ../../${BINARY_FILE} \
                    --tau="${tau}" \
                    --burstSizeV="${burstSizeV}" \
                    --burstSizeD="${burstSizeD}" \
                    --meanLysisTime="${meanLysisTime}" \
                    --kJumpR="${kJumpR}" \
                    --ifnBothFold="${ifnBothFold}" \
                    --rho="${rho}" \
                    --particleSpreadOption="${particleSpread}" \
                    --percentRandJumpR="${percentRandJumpR}" \
                    --ifnSpreadOption="${ifnSpread}" \
                    --dipOption="${dip}" \
                    --virion_half_life="${virion_half_life}" \
                    --dip_half_life="${dip_half_life}" \
                    --ifn_half_life="${ifn_half_life}"

                echo "Executed: ../../${BINARY_FILE} --tau=${tau} --particleSpreadOption=${particleSpread} --ifnSpreadOption=${ifnSpread} --dipOption=${dip}"
                echo "Run ${runCount} finished."
                popd > /dev/null

                ((runCount++))
            done
        done
    done

    popd > /dev/null
done

echo "All simulations completed."
