#!/bin/bash

# Author: Yimei Li  
# Affiliation: Princeton University, Grenfell Lab / teVelthuis Lab / Levin Lab  
# Year: 2024  
# Copyright: © 2024 Yimei Li. All rights reserved.  
# License: [Specify if applicable, e.g., MIT, GPL, or "Proprietary"]  

#==============================================================================
# VIRAL INFECTION SIMULATION WITH DIPs - PARAMETER CONFIGURATION
#==============================================================================
# Based on Howat et al. (2006) spatial model extended with DIPs and finite-rate IFN diffusion
# Model: Hexagonal lattice, 7 cell states, hourly updates, spatially explicit infection dynamics
# 
# CORE PARAMETERS AND THEIR BIOLOGICAL MEANINGS:
#
# INFECTION DYNAMICS:
# - rho: Base infection probability for MDBK cells (ρ = 0.026, from Howat et al.)
#   * For Vero cells: P_infection = ρ (no IFN effects)
#   * For MDBK cells: P_infection = ρ · exp(-αI) where α = 1.5, I = IFN concentration
# - meanLysisTime: Mean duration before virion/both infected cells die (12 hours, std dev 3 hours)
# - dvgRecoveryTime: Mean duration before DVG-only infected cells return to susceptible (3 hours, std dev 1 hour)
#   * INFECTED_VIRION/INFECTED_BOTH: lyse and die (release particles)
#   * INFECTED_DIP: return to SUSCEPTIBLE state (no particle release)
# - burstSizeV: Virion burst size per lysed cell (normally distributed around 50 pfu)
# - burstSizeD: DIP burst size per lysed cell (normally distributed around 100 pfu)
#
# IFN IMMUNE RESPONSE:
# - tau: Mean delay time for IFN-triggered antiviral state transition (12 hours, std dev 3 hours)
#   * Antiviral cells form local barriers preventing viral spread
#   * IFN decays with 3-hour half-life
# - ifnBothFold: IFN production multiplier for co-infected cells (should be 10.0 based on paper)
#   * DIP-only cells: 5-fold more IFN than virion-only cells  
#   * Both-infected cells: 10-fold more IFN than virion-only cells
# - ifnSpreadOption: IFN diffusion mechanism
#   * "local": Finite-rate diffusion within specified radius (biologically realistic)
#   * "global": Uniform IFN concentration across entire grid
#   * "noIFN": No interferon response (Vero cell behavior)
#
# PARTICLE TRANSMISSION:
# - particleSpreadOption: How viral particles spread between cells
#   * "celltocell": Only to neighboring cells (distance-weighted on hexagonal lattice)
#   * "jumprandomly": All particles can jump to any location (non-local spread)
#   * "partition": Mixed mode - kJumpR fraction jumps randomly, rest cell-to-cell
#   * "jumpradius": Particles jump within specified radius
# - kJumpR: Fraction of particles that spread non-locally (0-1, used with "partition" mode)
#
# DIP EFFECTS:
# - dipOption: Enable/disable Defective Interfering Particles
#   * "true": DIPs compete with virions, enhance IFN response, provide interference
#   * "false": Virion-only simulation (reproduces Howat et al. baseline)
#
# CELL DYNAMICS:
# - Dead cells regrow when healthy neighbors present (mean 24 hours, std dev 6 hours)
# - Regrown cells become susceptible again (reinfectable)
# - INFECTED_VIRION/INFECTED_BOTH cells: lyse and die after ~12±3 hours (release particles)
# - INFECTED_DIP cells: return to SUSCEPTIBLE state after ~3±1 hours (DVG recovery, no particle release)
# - Model captures infection → lysis/recovery → IFN → antiviral → death/regrowth cycle
#
# EXPERIMENTAL FEATURES:
# - enableParticleRemoval: Remove viral particles outside IFN protection range
#   * "false": Normal simulation (recommended)
#   * "true": Experimental particle removal at specified timepoint
# - removalTimepoint: When to remove particles (hours, only if enableParticleRemoval=true)
# - ifnThreshold: IFN concentration threshold for particle removal
# - removeVirionAndDIP: Remove both virions and DIPs (true) or only virions (false)
#
# EXPERIMENTAL FEATURES - VIRAL PARTICLE REMOVAL EXPERIMENT:
# This experiment tests the hypothesis that DIPs provide protective effects by removing
# viral particles outside the IFN protection range at specific timepoints.
# - enableParticleRemoval: Enable/disable the experiment ("true"/"false")
# - removalTimepoint: When to remove particles (hours, typically 72h)
# - ifnThreshold: IFN concentration threshold defining "protection range"
# - removeVirionAndDIP: Remove both particle types ("true") or only virions ("false")
# Expected result: Without DIPs, particle removal should reduce infection significantly.
# With DIPs, infection may persist despite removal, indicating DIP protective effects.
#
# SIMULATION CONTROL:
# - useFixedSeed: Random number generation control
#   * "no"/"false": Different results each run (for statistical analysis)
#   * "yes"/"true": Fixed seed=42 (for reproducibility)
#
# EXPOSURE / INOCULUM HETEROGENEITY (Exposure Mask)
# - Purpose: represent unexposed or low-coverage regions of a well that remain effectively non-infectable.
# - Parameter: unexposedAreaFraction in [0.0, 1.0], the fraction of area treated as non-exposed.
# - Default: 0.10 (≈10% area), reflecting typical edge effects, bubbles, plating gaps.
# - Activation rule: Effective only when videotype = "baltes". For all other videotypes it is auto-set to 0.0 (disabled).
# - Visualization intent: These regions would appear black and never become infected in baltes-style outputs.
# - Note: This script-level flag documents the assumption and can be consumed by downstream tooling.
#   The current Go binary does not read this flag; future integration may apply a spatially correlated mask
#   (clustered/edge-like rather than uniform random) with near-zero susceptibility inside masked areas.
#
# GRID AND TIMING:
# - Grid size: 50×50 cells (small scale) or 360×360 (full scale, matches Howat et al.)
# - Simulation time: ~500 hours (fixed in code)
# - Time step: 1 hour per frame
# - Initial inoculum: V_PFU_INITIAL = 1 virion, D_PFU_INITIAL varies by condition
#==============================================================================

echo "🦠 Viral Infection Simulation with DIPs"
SOURCE_FILE="mdbk_small_vero_0818.go"
BINARY_FILE="mdbk_small_vero_0818"

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



# CASE 4 SPECIFIC PARAMETERS
case4_spread_mode="burst"   # Case 4 spread mode: "burst" (traditional lysis) or "continuous" (persistent production) - CHANGE THIS LINE TO SWITCH MODES
# if continuous mode, set the following parameters
continuousProductionRateV=125    # Virion production rate per timestep for Case 4 continuous mode
continuousProductionRateD=500    # DIP production rate per timestep for Case 4 continuous mode (V:D = 1:4 base ratio)  
continuousIncubationPeriod=7    # Hours before cells start producing (Case 4)
continuousLysisTime=18.0        # Lysis time for continuous production cells (Case 4)

# CORE BIOLOGICAL PARAMETERS (based on Materials and Methods)
tau=0                           # IFN response delay parameter (paper: τ=12h, but code implementation may differ)
burstSizeV=100 # 790                  # Virion burst size - USER SUGGESTED: increased to 200 for exponential growth
burstSizeD=100  # 100                  # DIP burst size (increased to make distribution difference more visible)
v_pfu_initial=1                # Initial virion count - USER REQUESTED: direct value setting
d_pfu_initial=30                # Initial DIP count (0 for no DIP, 1+ for DIP experiments) - Case 4 hardcoded to 30  
# For case 4: Use a single dipInitRange M; we draw from [M-M/4, M+M/4] with sd=M/8.

meanLysisTime=6             # Mean infection duration - USER SUGGESTED: increased for stability (was 15h)
dvgRecoveryTime=0           # Mean recovery time for DVG-only infected cells - FROM SUCCESSFUL CONFIG 16
kJumpR=NaN                       # Fraction of particles spreading non-locally - SMALL amount for partition mode
burstRadius=11                   # Burst radius (neighbor circles) - INCREASED for wider spread
# DIP RADIUS CONFIGURATION
dipRadius=27                  # Absolute DIP spread radius for bursts (cells)
ifnBothFold=1.0                 # IFN multiplier for co-infected cells (paper: should be 10.0)
rho=0.8 #0.026                         # Base infection probability - REASONABLE for DIP infection testing

# PARTICLE AND IFN DECAY RATES  
virion_half_life=4.0            # Virion decay rate - USER SUGGESTED: NO DECAY (infinite survival)
dip_half_life=4.0               # DIP decay rate - INCREASED for longer survival
ifn_half_life=0.0               # IFN decay rate (paper: 3-hour half-life, but 0 = no decay here)

# SIMULATION CONTROL
option=4                       # Cell initialization method - options: 1, 2, 3, 4 (4=Case 4 with choice of burst/continuous)

# Case 4 mode is set directly in parameters above (line 112: case4_spread_mode)
# No interactive selection needed - change case4_spread_mode variable to "burst" or "continuous"
if [ "$option" -eq 4 ]; then
    echo "🎯 Case 4 detected: Using ${case4_spread_mode} mode"
    if [ "$case4_spread_mode" = "continuous" ]; then
        echo "   Parameters: continuousProductionRateV=${continuousProductionRateV}, continuousProductionRateD=${continuousProductionRateD}"
        echo "   Incubation: ${continuousIncubationPeriod}h, Lysis time: ${continuousLysisTime}h"
    else
        echo "   Parameters: burstSizeV=${burstSizeV}, burstSizeD=${burstSizeD}"
    fi
    echo ""
fi

videotype="baltes"    # Output visualization type: can be "states", "IFNconcentration", "IFNonlyLargerThanZero", "antiviralState", or "particles"

# BIOLOGICAL PROCESS OPTIONS
particleSpreadOption="celltocell"  # Viral spread mechanism: can be "celltocell", "jumprandomly", "partition", or "jumpradius"
ifnSpreadOption="noIFN"           # IFN diffusion: can be "local", "global", or "noIFN"
dipOption="true"                 # Enable DIPs: can be "true" or "false"
# Virion-only burst mode for VIRION-infected cells:
#   - "both": VIRION-only infected cells release both virions and DIPs (original behavior)
#   - "virionOnly": VIRION-only infected cells release only virions (no DIPs)
virionBurstMode="virionOnly" # "both" or "virionOnly"


# DIP hotspot mode for Case 4 initialization
# - dipHotspotMode: "random" (choose DIP hotspot within burstRadius around center) or "fixed"
# - If "fixed", set dipHotspotX, dipHotspotY (0..GRID_SIZE-1). Cells in UNEXPOSED will be shifted to nearest unmasked.
dipHotspotMode="fixed"
dipHotspotX=50
dipHotspotY=50
dipInitRange=5
# WELL EXPOSURE MASK (baltes-only)
# - Fraction [0.0–1.0] of lattice considered non-exposed (non-infectable) due to edge effects/bubbles/coverage gaps
# - Only effective when videotype = "baltes"; otherwise forced to 0.0 (disabled)
# - Seed control: uses the same RNG as the simulation (useFixedSeed/randomSeed). When useFixedSeed=yes, mask is reproducible; when no, mask is randomized each run.
# Visualization-only overlay fraction (baltes-only). Does not affect simulation logic
unexposedAreaFraction=0.0
unexposedSetAreaFraction=0.20

# Blocking Cell-Free Virions VIRAL PARTICLE REMOVAL EXPERIMENT (for testing DIP protective effects)
# This experiment tests the hypothesis that DIPs provide protective effects against infection spread
# by removing viral particles outside the IFN protection range at specific timepoints
enableParticleRemoval="false"    # Enable/disable the experiment: can be "true" or "false"
removalTimepoint=72               # Timepoint for particle removal (hours, only used if enableParticleRemoval=true)
removeVirionAndDIP="true"         # Particle removal scope: "true" (both virions and DIPs) or "false" (only virions)
ifnThreshold=0.1                  # IFN concentration threshold for protection (0.0-1.0+, only used if enableParticleRemoval=true)

# REPRODUCIBILITY CONTROL
useFixedSeed="yes"                # Random seed control: can be "yes"/"true" (reproducible) or "no"/"false" (different results)

# Set randomSeed value based on useFixedSeed setting
if [[ "$useFixedSeed" == "yes" || "$useFixedSeed" == "true" ]]; then
    randomSeed=42
    echo "🎲 Fixed random seed enabled: randomSeed = ${randomSeed} (reproducible results)"
else
    randomSeed=-1
    echo "🎲 Random seed disabled: randomSeed = ${randomSeed} (different results each time)"
fi

# Enforce exposure mask applicability and range before validation
# - Clamp unexposedAreaFraction to [0.0, 1.0]
# - Disable automatically unless videotype = "baltes"
if command -v bc > /dev/null 2>&1; then
    lt0=$(echo "${unexposedAreaFraction} < 0" | bc)
    gt1=$(echo "${unexposedAreaFraction} > 1" | bc)
    if [[ "$lt0" == "1" ]]; then unexposedAreaFraction=0.0; fi
    if [[ "$gt1" == "1" ]]; then unexposedAreaFraction=1.0; fi
fi
if [[ "$videotype" != "baltes" ]]; then
    unexposedAreaFraction=0.0
fi

# Parameter validation function
validate_parameters() {
    echo "-------------------------------------------"
    echo "Validating parameter combinations..."
    
    # Check if kJumpR is a numerical value (not NaN)
    if [[ "$kJumpR" != "NaN" ]]; then
        # Use bc for floating point comparison (if available, otherwise use awk)
        if command -v bc > /dev/null 2>&1; then
            kJumpR_gt_0=$(echo "$kJumpR > 0" | bc)
            kJumpR_eq_1=$(echo "$kJumpR == 1" | bc)
        else
            # Use awk as fallback
            kJumpR_gt_0=$(awk "BEGIN {print ($kJumpR > 0)}")
            kJumpR_eq_1=$(awk "BEGIN {print ($kJumpR == 1)}")
        fi
        
        # If kJumpR > 0 and != 1, must use partition mode
        if [[ "$kJumpR_gt_0" == "1" && "$kJumpR_eq_1" == "0" ]]; then
            if [[ "$particleSpreadOption" != "partition" ]]; then
                echo "❌ ERROR: Invalid parameter combination!"
                echo "   When kJumpR > 0 and kJumpR != 1 (kJumpR = $kJumpR),"
                echo "   particleSpreadOption MUST be 'partition' for mixed transmission."
                echo "   Current particleSpreadOption: '$particleSpreadOption'"
                echo "   "
                echo "   Valid combinations:"
                echo "   - kJumpR = NaN + particleSpreadOption = 'celltocell' (cell-to-cell only)"
                echo "   - kJumpR = 1.0 + particleSpreadOption = 'jumprandomly' (free jump only)"
                echo "   - kJumpR ∈ (0,1) + particleSpreadOption = 'partition' (mixed mode)"
                echo "-------------------------------------------"
                exit 1
            fi
        # If kJumpR = 1, should use jumprandomly mode
        elif [[ "$kJumpR_eq_1" == "1" ]]; then
            if [[ "$particleSpreadOption" != "jumprandomly" ]]; then
                echo "⚠️  WARNING: When kJumpR = 1.0, particleSpreadOption should be 'jumprandomly'"
                echo "   Current: kJumpR = $kJumpR, particleSpreadOption = '$particleSpreadOption'"
                echo "   Continuing anyway, but consider using 'jumprandomly' for clarity."
            fi
        fi
    else
        # kJumpR = NaN, should use celltocell mode
        if [[ "$particleSpreadOption" != "celltocell" ]]; then
            echo "⚠️  WARNING: When kJumpR = NaN, particleSpreadOption should be 'celltocell'"
            echo "   Current: kJumpR = $kJumpR, particleSpreadOption = '$particleSpreadOption'"
            echo "   Continuing anyway, but consider using 'celltocell' for clarity."
        fi
    fi
    
    # Check if viral particle removal experiment parameters are set when feature is disabled
    if [[ "$enableParticleRemoval" == "false" ]]; then
        echo "Checking viral particle removal experiment parameters..."
        
        # Default values for particle removal experiment
        DEFAULT_REMOVAL_TIMEPOINT=72
        DEFAULT_REMOVE_VIRION_AND_DIP="true"
        DEFAULT_IFN_THRESHOLD=0.1
        
        WARNINGS_FOUND=false
        
        # Check if removalTimepoint is non-default
        if [[ "$removalTimepoint" != "$DEFAULT_REMOVAL_TIMEPOINT" ]]; then
            echo "⚠️  WARNING: removalTimepoint is set to $removalTimepoint but enableParticleRemoval=false"
            echo "   This parameter will be ignored. Set enableParticleRemoval=true to use particle removal."
            WARNINGS_FOUND=true
        fi
        
        # Check if removeVirionAndDIP is non-default
        if [[ "$removeVirionAndDIP" != "$DEFAULT_REMOVE_VIRION_AND_DIP" ]]; then
            echo "⚠️  WARNING: removeVirionAndDIP is set to $removeVirionAndDIP but enableParticleRemoval=false"
            echo "   This parameter will be ignored. Set enableParticleRemoval=true to use particle removal."
            WARNINGS_FOUND=true
        fi
        
        # Check if ifnThreshold is non-default
        if [[ "$ifnThreshold" != "$DEFAULT_IFN_THRESHOLD" ]]; then
            echo "⚠️  WARNING: ifnThreshold is set to $ifnThreshold but enableParticleRemoval=false"
            echo "   This parameter will be ignored. Set enableParticleRemoval=true to use particle removal."
            WARNINGS_FOUND=true
        fi
        
        if [[ "$WARNINGS_FOUND" == "true" ]]; then
            echo "   "
            echo "   NOTE: You have configured parameters for the 'Viral Particle Removal Experiment'"
            echo "   but the experiment is disabled (enableParticleRemoval=false)."
            echo "   To run the particle removal experiment, set enableParticleRemoval=true"
            echo "   "
        fi
    else
        echo "✅ Viral particle removal experiment enabled - all related parameters will be used."
    fi

    # Validate DIP hotspot fixed mode inputs
    if [[ "$dipHotspotMode" == "fixed" ]]; then
        if ! [[ "$dipHotspotX" =~ ^[0-9]+$ ]] || ! [[ "$dipHotspotY" =~ ^[0-9]+$ ]]; then
            echo "⚠️  WARNING: dipHotspotMode=fixed but dipHotspotX/dipHotspotY are not non-negative integers"
            echo "   Current: dipHotspotX=$dipHotspotX, dipHotspotY=$dipHotspotY"
            echo "   The Go program will fallback to nearest valid/unmasked cell if out of bounds."
        fi
    fi
    
    echo "✅ Parameter validation passed!"
    echo "-------------------------------------------"
}



# Script path variable
SCRIPT_PATH="$0"

# Validate parameter combinations
validate_parameters

echo "-------------------------------------------"
echo "Running viral infection simulation with the following configuration:"
echo "  TAU                     = ${tau}"
echo "  burstSizeV              = ${burstSizeV}"
echo "  burstSizeD              = ${burstSizeD}"
echo "  v_pfu_initial           = ${v_pfu_initial}"
echo "  d_pfu_initial           = ${d_pfu_initial}"
echo "  meanLysisTime           = ${meanLysisTime}"
echo "  dvgRecoveryTime         = ${dvgRecoveryTime}"
echo "  kJumpR                  = ${kJumpR}"
echo "  burstRadius             = ${burstRadius}"
echo "  ifnBothFold             = ${ifnBothFold}"
echo "  rho                     = ${rho}"
echo "  virion_half_life        = ${virion_half_life}"
echo "  dip_half_life           = ${dip_half_life}"
echo "  ifn_half_life           = ${ifn_half_life}"
echo "  option                  = ${option}"

# Case 4 specific parameters display
if [ "$option" -eq 4 ]; then
    echo "  🎯 CASE 4 PARAMETERS:"
    echo "    case4_spread_mode     = ${case4_spread_mode}"
    if [ "$case4_spread_mode" = "continuous" ]; then
            echo "    continuousProductionRateV       = ${continuousProductionRateV}"
    echo "    continuousProductionRateD       = ${continuousProductionRateD}"
    echo "    continuousIncubationPeriod      = ${continuousIncubationPeriod}"
        echo "    continuousLysisTime   = ${continuousLysisTime}"
    fi
fi

echo "  particleSpreadOption    = ${particleSpreadOption}"
echo "  ifnSpreadOption         = ${ifnSpreadOption}"
echo "  dipOption               = ${dipOption}"
echo "  videotype               = ${videotype}"
echo "  unexposedAreaFraction   = ${unexposedAreaFraction}  # Exposure mask (baltes-only)"
echo "  dipHotspotMode          = ${dipHotspotMode}  # Case 4 DIP hotspot init"
echo "  dipHotspotX,Y           = ${dipHotspotX}, ${dipHotspotY}  # Used when mode=fixed"
echo "  enableParticleRemoval   = ${enableParticleRemoval}"
echo "  removalTimepoint        = ${removalTimepoint}"
echo "  removeVirionAndDIP      = ${removeVirionAndDIP}"
echo "  ifnThreshold            = ${ifnThreshold}"
echo "  useFixedSeed            = ${useFixedSeed} (randomSeed = ${randomSeed})"
echo "-------------------------------------------"

# 执行仿真
./${BINARY_FILE} \
    -tau="${tau}" \
    -burstSizeV="${burstSizeV}" \
    -burstSizeD="${burstSizeD}" \
    -v_pfu_initial="${v_pfu_initial}" \
    -d_pfu_initial="${d_pfu_initial}" \
    -meanLysisTime="${meanLysisTime}" \
    -dvgRecoveryTime="${dvgRecoveryTime}" \
    -kJumpR="${kJumpR}" \
    -burstRadius="${burstRadius}" \
    -ifnBothFold="${ifnBothFold}" \
    -rho="${rho}" \
    -virion_half_life="${virion_half_life}" \
    -dip_half_life="${dip_half_life}" \
    -ifn_half_life="${ifn_half_life}" \
    -particleSpreadOption="${particleSpreadOption}" \
    -ifnSpreadOption="${ifnSpreadOption}" \
    -dipOption="${dipOption}" \
    -option="${option}" \
    -videotype="${videotype}" \
    -unexposedAreaFraction="${unexposedAreaFraction}" \
    -unexposedSetAreaFraction="${unexposedSetAreaFraction}" \
    -enableParticleRemoval="${enableParticleRemoval}" \
    -removalTimepoint="${removalTimepoint}" \
    -removeVirionAndDIP="${removeVirionAndDIP}" \
    -ifnThreshold="${ifnThreshold}" \
    -randomSeed="${randomSeed}" \
    -virionBurstMode="${virionBurstMode}" \
    -dipRadius="${dipRadius}" \
    -dipRestrictToHotspot="false" \
    -dipInitRange="${dipInitRange}" \
    -dipHotspotMode="${dipHotspotMode}" \
    -dipHotspotX="${dipHotspotX}" \
    -dipHotspotY="${dipHotspotY}" \
    -continuousMode=$([ "$case4_spread_mode" = "continuous" ] && echo "true" || echo "false") \
    -continuousProductionRateV="${continuousProductionRateV}" \
-continuousProductionRateD="${continuousProductionRateD}" \
-continuousIncubationPeriod="${continuousIncubationPeriod}" \
    -continuousLysisTime="${continuousLysisTime}" \
	-lambdaDip 0.3

echo "✅ Simulation finished."

# Copy configuration script to the generated simulation folder
SIMULATION_FOLDER=$(ls -dt *_D* 2>/dev/null | head -1)
if [ -n "${SIMULATION_FOLDER}" ] && [ -d "${SIMULATION_FOLDER}" ]; then
    cp "${SCRIPT_PATH}" "${SIMULATION_FOLDER}/"
    echo "📋 Configuration script copied to: ${SIMULATION_FOLDER}/$(basename "$SCRIPT_PATH")"
    echo "🖼️  Check ${SIMULATION_FOLDER}/simulation_frames_combined.png for simulation results"
    
    # Copy experimental result images to simulation folder
    echo "📂 Copying experimental result images..."
    if [ -f "7_hours.png" ]; then
        cp "7_hours.png" "${SIMULATION_FOLDER}/"
        echo "   ✅ Copied 7_hours.png"
    else
        echo "   ⚠️  7_hours.png not found in main directory"
    fi
    
    if [ -f "13_hours.png" ]; then
        cp "13_hours.png" "${SIMULATION_FOLDER}/"
        echo "   ✅ Copied 13_hours.png"
    else
        echo "   ⚠️  13_hours.png not found in main directory"
    fi
    
    if [ -f "19_hours.png" ]; then
        cp "19_hours.png" "${SIMULATION_FOLDER}/"
        echo "   ✅ Copied 19_hours.png"
    else
        echo "   ⚠️  19_hours.png not found in main directory"
    fi
    
    if [ -f "25_hours_comprehensive.png" ]; then
        cp "25_hours_comprehensive.png" "${SIMULATION_FOLDER}/25_hours.png"
        echo "   ✅ Copied 25_hours_comprehensive.png as 25_hours.png"
    else
        echo "   ⚠️  25_hours_comprehensive.png not found in main directory"
    fi
    
    # Create 4x2 composite image with existing and selected frame images
    echo "-------------------------------------------"
    echo "🖼️  Creating 4x2 composite image..."
    
    # Define the image files to combine
    LEFT_IMAGES=(
        "${SIMULATION_FOLDER}/7_hours.png"
        "${SIMULATION_FOLDER}/13_hours.png"
        "${SIMULATION_FOLDER}/19_hours.png"
        "${SIMULATION_FOLDER}/25_hours.png"
    )
    
    RIGHT_IMAGES=(
        "${SIMULATION_FOLDER}/simulation_7_hours.png"
        "${SIMULATION_FOLDER}/simulation_13_hours.png"
        "${SIMULATION_FOLDER}/simulation_19_hours.png"
        "${SIMULATION_FOLDER}/simulation_25_hours.png"
    )
    
    # Check if ImageMagick is available
    if command -v magick > /dev/null 2>&1; then
        MAGICK_CMD="magick"
    elif command -v convert > /dev/null 2>&1; then
        MAGICK_CMD="convert"
    else
        echo "⚠️  ImageMagick not found. Cannot create composite image."
        echo "   Please install ImageMagick to generate the 4x2 composite image."
        echo "   Install with: brew install imagemagick (macOS) or apt install imagemagick (Ubuntu)"
        echo "📁 Results saved in the simulation-generated folder"  
        echo "🖼️  Check */selected_frames_combined.png for visualization results"
        exit 0
    fi
    
    # Check if all required images exist and show detailed status
    echo "🔍 Checking for required images..."
    ALL_IMAGES_EXIST=true
    MISSING_IMAGES=()
    FOUND_IMAGES=()
    
    for img in "${LEFT_IMAGES[@]}" "${RIGHT_IMAGES[@]}"; do
        if [ ! -f "$img" ]; then
            ALL_IMAGES_EXIST=false
            MISSING_IMAGES+=("$(basename "$img")")
            echo "   ❌ Missing: $(basename "$img")"
        else
            FOUND_IMAGES+=("$(basename "$img")")
            echo "   ✅ Found: $(basename "$img")"
        fi
    done
    
    echo "📊 Status: Found ${#FOUND_IMAGES[@]}/8 required images"
    
    if [ "$ALL_IMAGES_EXIST" = true ]; then
        # Create the 4x2 composite image
        OUTPUT_COMPOSITE="${SIMULATION_FOLDER}/composite_4x2_comparison.png"
        
        # Using ImageMagick to create a 4x2 grid (4 rows, 2 columns)
        # Resize all images to 800x800 to ensure consistent dimensions and 1:1 aspect ratio
        # Left column: original hour images, Right column: selected frame images
        if [ "$MAGICK_CMD" = "magick" ]; then
            # ImageMagick 7 syntax - resize all images to 800x800 before combining
            magick \
                \( "${LEFT_IMAGES[0]}" -resize 800x800! "${RIGHT_IMAGES[0]}" -resize 800x800! +append \) \
                \( "${LEFT_IMAGES[1]}" -resize 800x800! "${RIGHT_IMAGES[1]}" -resize 800x800! +append \) \
                \( "${LEFT_IMAGES[2]}" -resize 800x800! "${RIGHT_IMAGES[2]}" -resize 800x800! +append \) \
                \( "${LEFT_IMAGES[3]}" -resize 800x800! "${RIGHT_IMAGES[3]}" -resize 800x800! +append \) \
                -append "${OUTPUT_COMPOSITE}"
        else
            # ImageMagick 6 syntax - resize all images to 800x800 before combining
            convert \
                \( "${LEFT_IMAGES[0]}" -resize 800x800! "${RIGHT_IMAGES[0]}" -resize 800x800! +append \) \
                \( "${LEFT_IMAGES[1]}" -resize 800x800! "${RIGHT_IMAGES[1]}" -resize 800x800! +append \) \
                \( "${LEFT_IMAGES[2]}" -resize 800x800! "${RIGHT_IMAGES[2]}" -resize 800x800! +append \) \
                \( "${LEFT_IMAGES[3]}" -resize 800x800! "${RIGHT_IMAGES[3]}" -resize 800x800! +append \) \
                -append "${OUTPUT_COMPOSITE}"
        fi
        
        if [ $? -eq 0 ]; then
            echo "✅ Successfully created 4x2 composite image: $(basename "$OUTPUT_COMPOSITE")"
            echo "   Left column: Experimental results (7h, 13h, 19h, 25h) - resized to 800x800"
            echo "   Right column: Simulation results (7h, 13h, 19h, 25h) - resized to 800x800"
            echo "   All images now have consistent 1:1 aspect ratio and equal dimensions"
        else
            echo "❌ Failed to create composite image"
        fi
    else
        echo "⚠️  Cannot create complete 4x2 composite. Missing ${#MISSING_IMAGES[@]} files:"
        printf '   - %s\n' "${MISSING_IMAGES[@]}"
        echo ""
        echo "💡 To fix this issue:"
        echo "   1. Make sure experimental result images (7_hours.png, 13_hours.png, etc.) exist in main directory"
        echo "   2. Check that simulation runs to completion and generates simulation_*_hours.png files"
        echo "   3. Re-run the simulation after ensuring all images are available"
        echo ""
        echo "📋 Available files for manual composite creation:"
        for img in "${FOUND_IMAGES[@]}"; do
            echo "   ✅ $img"
        done
    fi
    
    echo "-------------------------------------------"
else
    echo "📁 Results saved in the simulation-generated folder"  
    echo "🖼️  Check */simulation_frames_combined.png for simulation results"
fi
