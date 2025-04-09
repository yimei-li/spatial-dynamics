#!/bin/bash
# 脚本说明：如果已经编译好了可执行文件，则直接运行，否则先编译。
# 然后循环运行仿真，传入不同的 particleSpreadOption、ifnSpreadOption 和 dipOption 组合（共 18 种组合）。
# 每一种组合将在输出文件夹下生成一个子文件夹，子文件夹命名方式为：
#   <particleSpreadOption>_<ifnSpreadOption>_dip<dipOption>
# 如果输出文件夹（默认 combine18）已存在，则在名称后面加上 _数字 后缀。
# 在对应的子文件夹中运行仿真（调用上级目录的可执行文件），仿真结果会保存在对应子文件夹内。

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

# 设置基本输出文件夹名称


# 固定参数（可根据需要修改）
tau=12
burstSizeV=50
burstSizeD=100
meanLysisTime=12.0
kJumpR=0.5
ifnBothFold=1.0
rho=0.026
percentRandJumpR=1.0
virion_half_life=0.0 # 3.2
dip_half_life=0.0 # 3.0
ifn_half_life=50.0 # 3.0


BASE_OUTPUT_FOLDER="combine18_ifnclr${ifn_half_life}_vclr${virion_half_life}_dclr${dip_half_life}"



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


# 定义需要循环的三个选项数组
particleSpreadOptions=("celltocell" "jumprandomly" "jumpradius")
ifnSpreadOptions=("global" "local" "noIFN")
dipOptions=("true" "false")

runCount=1

# 进入输出文件夹，这样生成的所有子文件夹都会在 BASE_OUTPUT_FOLDER 内
pushd "${BASE_OUTPUT_FOLDER}" > /dev/null

# 循环遍历所有组合，共 18 种可能
for particleSpread in "${particleSpreadOptions[@]}"; do
  for ifnSpread in "${ifnSpreadOptions[@]}"; do
    for dip in "${dipOptions[@]}"; do
      # 根据选项组合生成子文件夹名称
      SUBFOLDER="${particleSpread}_${ifnSpread}_dip${dip}"
      mkdir -p "${SUBFOLDER}"
      
      echo "-------------------------------------------"
      echo "Run ${runCount}:"
      echo "  TAU = ${tau}"
      echo "  particleSpreadOption = ${particleSpread}"
      echo "  ifnSpreadOption = ${ifnSpread}"
      echo "  dipOption = ${dip}"
      echo "  Output subfolder: ${SUBFOLDER}"
      
      # 进入子文件夹，并运行仿真（调用上级目录中的可执行文件）
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

echo "All simulations completed. Outputs are saved in the '${BASE_OUTPUT_FOLDER}' folder."
