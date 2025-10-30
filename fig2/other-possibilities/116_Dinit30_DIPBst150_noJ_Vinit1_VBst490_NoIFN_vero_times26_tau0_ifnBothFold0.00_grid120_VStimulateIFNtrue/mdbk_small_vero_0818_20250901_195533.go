/*
 * Author: Yimei Li
 * Affiliation: Princeton University, Grenfell Lab / teVelthuis Lab / Levin Lab
 * Year: 2024
 * Copyright: © 2024 Yimei Li. All rights reserved.
 * License: Proprietary. All rights reserved.
 *
 * Usage: Used to generate IFN dynamics plots and summary statistics for simulations in my PhD thesis.
 *
 * NOTE: Burst size is fixed. The burst size for virions (burst size V) is a constant value in this simulation.
 */

//: A simpler script with a higher-resolution image for easier inspection. I changed the GRID_SIZE to a smaller value, like 10, and set HOWAT_V_PFU_INITIA to 1, a small initial number of virions. I also increased CELL_SIZE to 10 for a clearer image. Additionally, I changed the "option" from 3 to 2 so you can modify the initial virion location. In option 2, I set the initial location at [4][5] ( "g.localVirions[4][5]++" ) which you can change. RHO is 1, and the virus spreads to neighboring cells weighted by distance.

package main

import (
	"bytes"
	"encoding/csv"
	"flag"
	"fmt"
	"image"
	"image/color"
	"image/draw"
	"image/jpeg"
	"image/png"
	"io/ioutil"
	"log"
	"math"
	"math/rand"
	"os" // Used for file operations
	"os/exec"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"
	"time"

	"github.com/icza/mjpeg"
	"github.com/wcharczuk/go-chart/v2" // Used for plotting the graph
	"github.com/wcharczuk/go-chart/v2/drawing"
	"golang.org/x/image/font"
	"golang.org/x/image/font/basicfont"
	"golang.org/x/image/math/fixed"
)

// Constant definitions
const (
	TIME_STEPS = 26  // Number of time steps (0-25, so 26 total)
	GRID_SIZE  = 120 // Size of the grid

	FRAME_RATE   = 1          // Frame rate for the video
	OUTPUT_VIDEO = "0421.mp4" // Output video file name
	CELL_SIZE    = 4          // The size of each hexagonal cell
	TIMESTEP     = 1          // Time step size

)

// Define flag variables (note they are all pointer types)
var (

	// Option parameters
	// Particle spread option: can be "celltocell", "jumprandomly", "partition" or "jumpradius"
	flag_particleSpreadOption = flag.String("particleSpreadOption", "celltocell", "Particle spread option: celltocell, jumprandomly, jumpradius, or partition")
	// If jumprandomly is selected, this parameter represents the random jump ratio (0~1)

	// IFN spread option: can be "global", "local" or "noIFN"
	flag_ifnSpreadOption = flag.String("ifnSpreadOption", "local", "IFN spread option: global, local, or noIFN")
	// DIP option: if true then enable DIP, if false then disable DIP
	flag_dipOption = flag.Bool("dipOption", true, "DIP option: if true then enable DIP, if false then disable DIP")

	flag_burstSizeV       = flag.Int("burstSizeV", 50, "Number of virions released when a cell lyses")
	flag_burstSizeD       = flag.Int("burstSizeD", 100, "Number of DIPs released when a cell lyses")
	flag_meanLysisTime    = flag.Float64("meanLysisTime", 12.0, "Mean lysis time for virion/both infected cells")
	flag_dvgRecoveryTime  = flag.Float64("dvgRecoveryTime", 3.0, "Mean recovery time for DVG-only infected cells")
	flag_kJumpR           = flag.Float64("kJumpR", 0.5, "Parameter for cell-to-cell jump randomness")
	flag_tau              = flag.Int("tau", 12, "TAU value (e.g., lysis time)")
	flag_ifnBothFold      = flag.Float64("ifnBothFold", 1.0, "Fold effect for IFN stimulation")
	flag_rho              = flag.Float64("rho", 0.026, "Infection rate constant")
	flag_virion_half_life = flag.Float64("virion_half_life", 3.2, "Virion clearance rate (e.g., 3.2 d^-1)")
	flag_dip_half_life    = flag.Float64("dip_half_life", 3.2, "DIP clearance rate (e.g., 3.2 d^-1)")
	flag_ifn_half_life    = flag.Float64("ifn_half_life", 4.0, "IFN clearance rate (e.g., 3.0 d^-1)")
	flag_option           = flag.Int("option", 2, "Option for infection initialization (e.g., 1, 2, 3)")
	flag_burstRadius      = flag.Int("burstRadius", 3, "Burst radius (number of neighbor circles) - Controls how far virions and DIPs spread from infected cells")

	// Case 4 continuous production mode parameters
	flag_continuousMode             = flag.Bool("continuousMode", false, "Enable continuous production mode for case 4")
	flag_continuousProductionRateV  = flag.Int("continuousProductionRateV", 50, "Virion production rate per timestep for case 4 continuous mode")
	flag_continuousProductionRateD  = flag.Int("continuousProductionRateD", 25, "DIP production rate per timestep for case 4 continuous mode")
	flag_continuousIncubationPeriod = flag.Int("continuousIncubationPeriod", 6, "Hours before cells start producing (case 4 continuous mode)")
	flag_continuousLysisTime        = flag.Float64("continuousLysisTime", 20.0, "Lysis time for continuous production cells")

	// DIP infection probability parameter
	flag_lambdaDip = flag.Float64("lambdaDip", 30.0, "Poisson distribution lambda parameter for DIP infection probability")

	flag_v_pfu_initial = flag.Float64("v_pfu_initial", 1.0, "Initial PFU count for virions")
	flag_d_pfu_initial = flag.Float64("d_pfu_initial", 0.0, "Initial PFU count for DIPs")
	flag_videotype     = flag.String("videotype", "states", "Video type: states, IFNconcentration, IFNonlyLargerThanZero, antiviralState, particles, baltes")
	// Exposure mask (baltes-only): fraction of area treated as non-exposed (uniformly sampled)
	flag_unexposedAreaFraction = flag.Float64("unexposedAreaFraction", 0.0, "Fraction [0-1] of area treated as non-exposed/uninfectable (baltes-only; uniform)")
	// Visualization-only overlay (baltes-only): fraction of cells drawn as black, without affecting simulation state
	flag_unexposedSetAreaFraction = flag.Float64("unexposedSetAreaFraction", 0.0, "Fraction [0-1] of cells visually overlaid as black in baltes rendering (display-only)")

	// New experimental parameters for viral particle removal
	flag_enableParticleRemoval = flag.Bool("enableParticleRemoval", false, "Enable removal of viral particles outside IFN range")
	flag_removalTimepoint      = flag.Int("removalTimepoint", 72, "Timepoint (in hours) to remove viral particles outside IFN range")
	flag_ifnThreshold          = flag.Float64("ifnThreshold", 0.1, "IFN concentration threshold to determine cells outside IFN range")
	flag_removeVirionAndDIP    = flag.Bool("removeVirionAndDIP", true, "Remove both virions and DIPs (true) or only virions (false)")

	// VIRION-only burst mode: control whether virion-only infected cells release DIPs
	flag_virionBurstMode = flag.String("virionBurstMode", "virionOnly", "For VIRION-only infected cells: 'both' to release virions and DIPs (original), 'virionOnly' to release only virions (no DIPs)")

	// Random seed parameter
	flag_randomSeed = flag.Int64("randomSeed", -1, "Random seed for reproducible results (-1 for random seed based on time)")

	// DIP radius parameter
	flag_dipRadius = flag.Int("dipRadius", 10, "Absolute DIP spread radius for bursts (cells)")

	// Restrict DIPs to the initial fixed hotspot only (no DIP production anywhere)
	flag_dipRestrictToHotspot = flag.Bool("dipRestrictToHotspot", false, "If true, disable all DIP release during bursts; DIPs only exist at the initial fixed hotspot")

	// DIP hotspot initialization (case 4) controls
	flag_dipHotspotMode = flag.String("dipHotspotMode", "random", "DIP hotspot mode for case 4: 'random' (default) or 'fixed'")
	flag_dipHotspotX    = flag.Int("dipHotspotX", -1, "Fixed DIP hotspot X index (0..GRID_SIZE-1), used when dipHotspotMode='fixed'")
	flag_dipHotspotY    = flag.Int("dipHotspotY", -1, "Fixed DIP hotspot Y index (0..GRID_SIZE-1), used when dipHotspotMode='fixed'")

	// Initial DIP seeding range (case 4). If >=0, sample from mean=M with range [M - M/4, M + M/4], sd=M/8; if <0, falls back per rules
	flag_dipInitRange = flag.Int("dipInitRange", -1, "Target initial DIPs at hotspot (case 4). Draw from [M-M/4, M+M/4] with sd=M/8; set to -1 to disable")
)

// Particle spread related
var (
	particleSpreadOption  string  // "celltocell", "jumprandomly", "jumpradius"
	jumpRadiusV           int     // e.g., when "jumpradius" is selected, set to 5
	jumpRadiusD           int     // same as above
	jumpRandomly          bool    // whether to use random jump (true when "jumprandomly" is selected)
	k_JumpR               float64 // random jump ratio
	par_celltocell_random bool
)

// IFN spread related
var (
	ifnSpreadOption string // "global",
	//  "local", "noIFN"
	IFN_wave_radius int  // if ifnSpreadOption=="local", e.g., set to 10; "global" or "noIFN" set to 0
	ifnWave         bool // whether to enable IFN wave
)

// DIP related
var (
	dipOption bool // true to enable DIP, false to disable DIP

	// When DIP is enabled, default DIP-related ratios remain default; when disabled, set to 0
	D_only_IFN_stimulate_ratio float64 = 5.0 * ifnBothFold
	BOTH_IFN_stimulate_ratio   float64 = 10.0 * ifnBothFold
)

// Viral particle removal experiment related
var (
	enableParticleRemoval bool    // whether to enable viral particle removal experiment
	removalTimepoint      int     // timepoint at which to remove viral particles
	ifnThreshold          float64 // IFN concentration threshold for determining cells outside IFN range
	removeVirionAndDIP    bool    // whether to remove both virions and DIPs or only virions
)

// VIRION-only burst mode
var (
	virionBurstMode string // "both" or "virionOnly"
)

// Random seed related
var (
	randomSeed int64 // random seed for reproducible results (-1 for time-based seed)
)

// Global variables
var (
	// particleSpreadOption  = "jumpradius" // options: "celltocell", "jumprandomly", "jumpradius"
	// PartionParticleSpreadOption = false // options: "true" or "false"

	//ifnSpreadOption = "global" // options: "global", "local" or "noIFN"
	//dipOption =

	BURST_SIZE_V  int    // CHANGE 50 Number of virions released when a cell lyses
	BURST_SIZE_D  int    // CHANGE 100 // Number of DIPs released when a cell lyses
	VStimulateIFN = true // CHANGE if false then usually only DIP stimulate IFN in this situlation, not virion
	//jumpRandomly          = true // CHANGE
	//jumpRadiusV           = 0    // CHANGE Virion jump radius
	//jumpRadiusD           = 0    // CHANGE DIP jump radius
	//IFN_wave_radius       = 10   // CHANGE 10
	// this is true only when jumpRandomly is true

	TAU         int // 95
	ifnBothFold = 1.0
	//D_only_IFN_stimulate_ratio = 5.0 * ifnBothFold  // D/V *R *D_only_IFN_stimulate_ratio
	//BOTH_IFN_stimulate_ratio = 10.0 * ifnBothFold // D/V *R *D_only_IFN_stimulate_ratio

	// option    = 2        // Option for infection initialization
	//videotype = "states" // "states" // color in "states" or "IFNconcentration" or "IFNonlyLargerThanZero" or "antiviralState" or "particles"
	RHO    float64 //0.026    //0.02468  // 0.09 Infection rate constant
	option int
	// radius 10 of grid has 331 cells
	R int
	// radius 10 of grid has 331 cells, originally infected cell increases R IFN,
	ALPHA     = 1.0  // Parameter for infection probability (set to 1.5)
	lambdaDip = 30.0 // DIP infection probability parameter

	REGROWTH_MEAN              = 24.0  // Mean time for regrowth
	REGROWTH_STD               = 6.0   // Standard deviation for regrowth time
	MEAN_LYSIS_TIME            float64 // Mean lysis time for virion/both infected cells
	STANDARD_LYSIS_TIME        float64 // Standard deviation for lysis time for virion/both infected cells
	MEAN_DVG_RECOVERY_TIME     float64 // Mean recovery time for DVG-only infected cells
	STANDARD_DVG_RECOVERY_TIME float64 // Standard deviation for recovery time for DVG-only infected cells
	maxGlobalIFN               = -1.0  // used to track maximum IFN value
	globalIFN                  = -1.0  // global IFN concentration
	globalIFNperCell           = 0.0
	IFN_DELAY                  = 5
	STD_IFN_DELAY              = 1

	// allowVirionJump = jumpRadiusV > 0 || jumpRandomly // Allow virions to jump to other cells
	// allowDIPJump = jumpRadiusD > 0 || jumpRandomly // Allow DIPs to jump to other cells
	allowVirionJump bool
	allowDIPJump    bool
	//ifnWave = IFN_wave_radius > 0

	yMax          float64
	xMax          = float64(TIME_STEPS)
	ticksInterval float64 // Interval for X-axis ticks

	adjusted_DIP_IFN_stimulate   float64
	perParticleInfectionChance_V float64
	totalDeadFromBoth            int
	totalDeadFromV               int
	virionDiffusionRate          int
	dipDiffusionRate             int

	virion_half_life float64 //= 0.0 // 3.2 // ~4 d^-1 => half-life ~4.2 hours
	dip_half_life    float64 //= 0.0 // 3.2 // ~4 d^-1 => half-life ~4.2 hours
	ifn_half_life    float64 //= 0.0 // 3.0 // ~3 d^-1 => half-life ~5.5 hours
	videotype        string
	dipAdvantage     float64 // DIP advantage = burstSizeD / burstSizeV
)

// Cell state definitions
const (
	SUSCEPTIBLE     = 0 // Susceptible state
	INFECTED_VIRION = 1 // Infected by virion
	INFECTED_DIP    = 5 // Infected by DIP
	INFECTED_BOTH   = 6 // Infected by both virion and DIP
	DEAD            = 2 // Dead state
	ANTIVIRAL       = 3 // Antiviral state
	REGROWTH        = 4 // Regrowth state
	// Case 4 continuous production states
	INFECTED_VIRION_CONTINUOUS = 7 // Mature virion-infected cell, continuously producing
	INFECTED_DIP_CONTINUOUS    = 8 // Mature DIP-infected cell, continuously producing
	INFECTED_BOTH_CONTINUOUS   = 9 // Mature co-infected cell, continuously producing
	// UNEXPOSED: permanently non-exposed cell (never changes state)
	UNEXPOSED = 10
)

// Grid structure for storing the simulation state
type Grid struct {
	state                  [GRID_SIZE][GRID_SIZE]int        // State of the cells in the grid
	localVirions           [GRID_SIZE][GRID_SIZE]int        // Number of virions in each cell
	localDips              [GRID_SIZE][GRID_SIZE]int        // Number of DIPs in each cell
	IFNConcentration       [GRID_SIZE][GRID_SIZE]float64    // IFN concentration in each cell
	timeSinceInfectVorBoth [GRID_SIZE][GRID_SIZE]int        // Time since infection for each cell
	timeSinceInfectDIP     [GRID_SIZE][GRID_SIZE]int        // Time since infection for each cell
	timeSinceDead          [GRID_SIZE][GRID_SIZE]int        // Time since death for each cell
	timeSinceRegrowth      [GRID_SIZE][GRID_SIZE]int        // Time since regrowth for each cell
	timeSinceSusceptible   [GRID_SIZE][GRID_SIZE]int        // Time since cell became susceptible
	neighbors1             [GRID_SIZE][GRID_SIZE][6][2]int  // Neighbors at distance 1 (6 neighbors)
	neighbors2             [GRID_SIZE][GRID_SIZE][12][2]int // Neighbors at distance 2 (12 neighbors)
	neighbors3             [GRID_SIZE][GRID_SIZE][18][2]int // Neighbors at distance 3 (18 neighbors)
	neighbors4             [GRID_SIZE][GRID_SIZE][24][2]int // Neighbors at distance 4 (24 neighbors)
	neighbors5             [GRID_SIZE][GRID_SIZE][30][2]int // Neighbors at distance 5 (30 neighbors)
	neighbors6             [GRID_SIZE][GRID_SIZE][36][2]int // Neighbors at distance 6 (36 neighbors)
	neighbors7             [GRID_SIZE][GRID_SIZE][42][2]int // Neighbors at distance 7 (42 neighbors)
	neighbors8             [GRID_SIZE][GRID_SIZE][48][2]int // Neighbors at distance 8 (48 neighbors)
	neighbors9             [GRID_SIZE][GRID_SIZE][54][2]int // Neighbors at distance 9 (54 neighbors)
	neighbors10            [GRID_SIZE][GRID_SIZE][60][2]int // Neighbors at distance 10 (60 neighbors)
	neighborsBurstArea     [GRID_SIZE][GRID_SIZE][][2]int   // Neighbors within burst radius (configurable)
	neighborsIFNArea       [GRID_SIZE][GRID_SIZE][][2]int   // Neighbors within IFN wave radius
	stateChanged           [GRID_SIZE][GRID_SIZE]bool       // Flag to indicate if the state of a cell has changed
	antiviralDuration      [GRID_SIZE][GRID_SIZE]int        // Duration of antiviral state
	previousStates         [GRID_SIZE][GRID_SIZE]int        // Previous state of the cell
	antiviralFlag          [GRID_SIZE][GRID_SIZE]bool       // Flag to indicate if the cell is in the antiviral state
	timeSinceAntiviral     [GRID_SIZE][GRID_SIZE]int        // Time since the cell entered the antiviral state
	antiviralCellCount     int                              // Number of cells in the antiviral state
	totalAntiviralTime     int
	intraWT                [GRID_SIZE][GRID_SIZE]int // IntraWT
	intraDVG               [GRID_SIZE][GRID_SIZE]int // IntraDVG
	// Exposure mask: true marks cells as non-exposed/uninfectable (baltes-only)
	unexposedMask          [GRID_SIZE][GRID_SIZE]bool
	allowJumpRandomly      [][]bool
	totalRandomJumpVirions int                       // record total number of randomly jumping Virions
	totalRandomJumpDIPs    int                       // record total number of randomly jumping DIPs
	lysisThreshold         [GRID_SIZE][GRID_SIZE]int // fixed lysis time for each cell (virion/both infected)
	dipLysisThreshold      [GRID_SIZE][GRID_SIZE]int // fixed lysis time for each DIP-infected cell
	dipClearanceThreshold  [GRID_SIZE][GRID_SIZE]int // time steps until DIP-only infected cells become susceptible
	burstRadius            int                       // configurable burst radius for virus and DIP spread

	// Case 4 continuous production mode fields
	continuousMode             bool                       // whether continuous production mode is enabled
	continuousProductionRateV  int                        // virion production rate per timestep for continuous mode
	continuousProductionRateD  int                        // DIP production rate per timestep for continuous mode
	continuousIncubationPeriod int                        // hours before cells start producing in continuous mode
	continuousLysisTime        float64                    // lysis time for continuous production cells
	infectionTime              [GRID_SIZE][GRID_SIZE]int  // timestep when cell was infected (for incubation)
	isProducing                [GRID_SIZE][GRID_SIZE]bool // whether cell is actively producing
	initOption                 int                        // case number (1,2,3,4)

	// Per-cell DIP half-life (hours), sampled at initialization from N(mean=*flag_dip_half_life, std=2)
	dipHalfLife [GRID_SIZE][GRID_SIZE]float64
}

// Initialize the infection state
func (g *Grid) initializeInfection(option int) {
	// Set random seed - use provided seed or current time for randomness
	if randomSeed >= 0 {
		rand.Seed(randomSeed)
		fmt.Printf("Using fixed random seed: %d\n", randomSeed)
	} else {
		seed := time.Now().UnixNano()
		rand.Seed(seed)
		fmt.Printf("Using time-based random seed: %d\n", seed)
	}

	vInit := int(math.Round(*flag_v_pfu_initial))
	dInit := int(math.Round(*flag_d_pfu_initial))

	switch option {
	case 1:
		if vInit > 0 {
			g.localVirions[25][25] = vInit
		} else {
			fmt.Printf("v_pfu_initial < 0: %.2f\n", *flag_v_pfu_initial)
		}
		if dInit > 0 {
			g.localDips[25][25] = dInit
		} else {
			fmt.Printf("d_pfu_initial < 0: %.2f\n", *flag_d_pfu_initial)
		}
	case 2:
		if vInit > 0 && dInit > 0 {
			g.state[25][25] = INFECTED_BOTH
		} else if vInit > 0 {
			g.state[25][25] = INFECTED_VIRION
		} else if dInit > 0 {
			g.state[25][25] = INFECTED_DIP
		}
		g.localVirions[25][25] = vInit
		g.localDips[25][25] = dInit

	case 3:
		for k := 0; k < vInit; k++ {
			i := rand.Intn(GRID_SIZE)
			j := rand.Intn(GRID_SIZE)
			g.localVirions[i][j]++
		}
		for k := 0; k < dInit; k++ {
			i := rand.Intn(GRID_SIZE)
			j := rand.Intn(GRID_SIZE)
			g.localDips[i][j]++
		}
	case 4:
		// Place virions at configurable number of positions clustered around the center
		centerX := GRID_SIZE / 2
		centerY := GRID_SIZE / 2

		// Set state based on continuous mode
		if g.continuousMode {
			g.state[centerX][centerY] = INFECTED_VIRION_CONTINUOUS
			fmt.Printf("🌱 Initial cell set to INFECTED_VIRION_CONTINUOUS at (%d,%d)\n", centerX, centerY)
		} else {
			g.state[centerX][centerY] = INFECTED_VIRION
			fmt.Printf("🌱 Initial cell set to INFECTED_VIRION at (%d,%d)\n", centerX, centerY)
		}

		g.localVirions[centerX][centerY] = vInit // Add actual virion particles

		// 在中心的 burstRadius 半径内随机选择一个 DIP 热点 (hx, hy)
		r := g.burstRadius
		if r < 1 {
			r = 1
		}
		var hx, hy int
		if *flag_dipHotspotMode == "fixed" && *flag_dipHotspotX >= 0 && *flag_dipHotspotX < GRID_SIZE && *flag_dipHotspotY >= 0 && *flag_dipHotspotY < GRID_SIZE {
			// Use fixed hotspot
			hx, hy = *flag_dipHotspotX, *flag_dipHotspotY
		} else {
			// Build ring cells around center within radius r and choose randomly
			var burstArea [][2]int
			for rad := 1; rad <= r; rad++ {
				ring := generateHexRing(centerX, centerY, rad)
				for _, nb := range ring {
					nx, ny := nb[0], nb[1]
					if nx >= 0 && nx < GRID_SIZE && ny >= 0 && ny < GRID_SIZE {
						burstArea = append(burstArea, [2]int{nx, ny})
					}
				}
			}
			if len(burstArea) == 0 {
				burstArea = append(burstArea, [2]int{centerX, centerY})
			}
			idxHot := rand.Intn(len(burstArea))
			hx, hy = burstArea[idxHot][0], burstArea[idxHot][1]
		}

		// Move DIP hotspot if masked
		hx, hy = g.findNearestUnmasked(hx, hy)

		// 2) 初始 DIPs 数量仅由 d_pfu_initial 控制
		centerDIPs := 0
		if *flag_d_pfu_initial >= 0 {
			centerDIPs = int(math.Round(*flag_d_pfu_initial))
		} else {
			centerDIPs = 0
		}

		if centerDIPs == 0 {
			fmt.Printf("🚫 d_pfu_initial==0: no DIPs seeded for case 4 (hotspot at %d,%d)\n", hx, hy)
		} else if *flag_dipHotspotMode == "fixed" {
			// 固定热点：dipInitRange 控制铺开半径；<=0 表示单点
			initR := *flag_dipInitRange
			if initR <= 0 {
				g.localDips[hx][hy] += centerDIPs
				fmt.Printf("🎯 Hotspot at (%d,%d): placed %d DIPs at single point (initRange=%d)\n", hx, hy, centerDIPs, initR)
			} else {
				// 在热点为中心、半径 initR 内按距离加权分布
				hotArea := make([][2]int, 0, 1+6*initR*(initR+1)/2)
				hotArea = append(hotArea, [2]int{hx, hy})
				for rad := 1; rad <= initR; rad++ {
					ring := generateHexRing(hx, hy, rad)
					for _, nb := range ring {
						nx, ny := nb[0], nb[1]
						if nx >= 0 && nx < GRID_SIZE && ny >= 0 && ny < GRID_SIZE {
							hotArea = append(hotArea, [2]int{nx, ny})
						}
					}
				}
				weights := make([]float64, len(hotArea))
				totalW := 0.0
				for idx, cell := range hotArea {
					d := getHexDistanceBetweenPoints(hx, hy, cell[0], cell[1])
					w := 1.0 / (float64(d) + 0.1)
					weights[idx] = w
					totalW += w
				}
				distributed := 0
				for idx, cell := range hotArea {
					share := int(math.Floor(float64(centerDIPs) * (weights[idx] / totalW)))
					if share > 0 {
						g.localDips[cell[0]][cell[1]] += share
						distributed += share
					}
				}
				left := centerDIPs - distributed
				if left > 0 && len(hotArea) > 0 {
					indices := make([]int, len(hotArea))
					for i2 := range indices {
						indices[i2] = i2
					}
					rand.Shuffle(len(indices), func(a, b int) { indices[a], indices[b] = indices[b], indices[a] })
					k := 0
					for left > 0 {
						idx := indices[k%len(indices)]
						cell := hotArea[idx]
						g.localDips[cell[0]][cell[1]]++
						left--
						k++
					}
				}
				fmt.Printf("🎯 Hotspot at (%d,%d): distributed %d DIPs within initRange=%d (distance-weighted)\n", hx, hy, centerDIPs, initR)
			}
		} else {
			// 随机热点：dipInitRange>=0 则用其作为铺开半径，否则使用 burstRadius r
			initR := r
			if *flag_dipInitRange >= 0 {
				initR = *flag_dipInitRange
			}
			// 在热点为中心、半径 initR 内按距离加权分布
			hotArea := make([][2]int, 0)
			hotArea = append(hotArea, [2]int{hx, hy})
			for rad := 1; rad <= initR; rad++ {
				ring := generateHexRing(hx, hy, rad)
				for _, nb := range ring {
					nx, ny := nb[0], nb[1]
					if nx >= 0 && nx < GRID_SIZE && ny >= 0 && ny < GRID_SIZE {
						hotArea = append(hotArea, [2]int{nx, ny})
					}
				}
			}

			weights := make([]float64, len(hotArea))
			totalW := 0.0
			for idx, cell := range hotArea {
				d := getHexDistanceBetweenPoints(hx, hy, cell[0], cell[1])
				w := 1.0 / (float64(d) + 0.1)
				weights[idx] = w
				totalW += w
			}

			distributed := 0
			for idx, cell := range hotArea {
				share := int(math.Floor(float64(centerDIPs) * (weights[idx] / totalW)))
				if share > 0 {
					g.localDips[cell[0]][cell[1]] += share
					distributed += share
				}
			}
			left := centerDIPs - distributed
			if left > 0 && len(hotArea) > 0 {
				indices := make([]int, len(hotArea))
				for i2 := range indices {
					indices[i2] = i2
				}
				rand.Shuffle(len(indices), func(a, b int) { indices[a], indices[b] = indices[b], indices[a] })
				k := 0
				for left > 0 {
					idx := indices[k%len(indices)]
					cell := hotArea[idx]
					g.localDips[cell[0]][cell[1]]++
					left--
					k++
				}
			}
		}

		// 不做额外随机邻居撒点，仅热点单点放置 DIPs

		// Record intracellular virus counts for continuous mode
		if g.continuousMode {
			g.intraWT[centerX][centerY] = 1       // Initial intracellular wild-type virus count
			g.intraDVG[centerX][centerY] = 0      // No DVG initially
			g.infectionTime[centerX][centerY] = 0 // Record infection time
		}

	}

}

// Function to generate ticks dynamically
func generateTicks(xMax float64, interval float64) []chart.Tick {
	var ticks []chart.Tick
	for value := 0.0; value <= xMax; value += interval {
		label := fmt.Sprintf("%.0f", value) // Format the label as an integer
		ticks = append(ticks, chart.Tick{
			Value: value,
			Label: label,
		})
	}
	return ticks
}

// Initialize the grid, setting all cells to SUSCEPTIBLE
func (g *Grid) initialize() {
	for i := 0; i < GRID_SIZE; i++ {
		for j := 0; j < GRID_SIZE; j++ {
			g.state[i][j] = SUSCEPTIBLE
			g.stateChanged[i][j] = false // Initialize as unchanged
			g.timeSinceInfectVorBoth[i][j] = -1
			g.timeSinceDead[i][j] = -1
			g.timeSinceRegrowth[i][j] = -1
			g.IFNConcentration[i][j] = 0
			g.antiviralDuration[i][j] = -1
			g.timeSinceSusceptible[i][j] = 0
			g.previousStates[i][j] = -1
			g.antiviralFlag[i][j] = false
			g.timeSinceAntiviral[i][j] = -1
			g.intraWT[i][j] = 0
			g.intraDVG[i][j] = 0
			g.lysisThreshold[i][j] = -1
			g.dipLysisThreshold[i][j] = -1
			g.dipClearanceThreshold[i][j] = -1

			// Initialize per-cell DIP half-life from Normal(mean=*flag_dip_half_life, std=2)
			// Clamp to a small positive minimum to avoid division by zero or negative values
			val := *flag_dip_half_life + 2.0*rand.NormFloat64()
			// Round to integer hours
			val = math.Round(val)
			if val < 1.0 {
				val = 1.0
			}
			g.dipHalfLife[i][j] = val
		}
	}

	// Initialize exposure mask (uniform sampling) if enabled for baltes
	maskFraction := *flag_unexposedAreaFraction
	if videotype == "baltes" && maskFraction > 0.0 {
		if maskFraction < 0.0 {
			maskFraction = 0.0
		}
		if maskFraction > 1.0 {
			maskFraction = 1.0
		}
		totalCells := GRID_SIZE * GRID_SIZE
		target := int(math.Round(maskFraction * float64(totalCells)))
		if target > 0 {
			// Build a list of all indices and sample without replacement uniformly
			indices := make([]int, totalCells)
			for idx := 0; idx < totalCells; idx++ {
				indices[idx] = idx
			}
			rand.Shuffle(totalCells, func(a, b int) { indices[a], indices[b] = indices[b], indices[a] })
			for k := 0; k < target; k++ {
				idx := indices[k]
				i := idx / GRID_SIZE
				j := idx % GRID_SIZE
				g.unexposedMask[i][j] = true
				g.state[i][j] = UNEXPOSED
			}
			fmt.Printf("Exposure mask initialized (uniform): fraction=%.3f, cells=%d\n", maskFraction, target)
		}
	}

	fmt.Println("Grid initialized")

}

// Ensure the entire canvas is initialized with uniform background color
func fillBackground(img *image.RGBA, bgColor color.Color) {
	for y := 0; y < img.Bounds().Dy(); y++ {
		for x := 0; x < img.Bounds().Dx(); x++ {
			img.Set(x, y, bgColor)
		}
	}
}

func (g *Grid) calculateDiffusionRates() (float64, float64) {
	totalVirions := 0
	totalVirionDiffusion := 0
	totalDIPs := 0
	totalDIPDiffusion := 0

	for i := 0; i < GRID_SIZE; i++ {
		for j := 0; j < GRID_SIZE; j++ {
			virionMoves := g.localVirions[i][j]
			dipMoves := g.localDips[i][j]

			// Sum particles moved out for diffusion
			totalVirionDiffusion += virionMoves
			totalDIPDiffusion += dipMoves

			// Add to total particles in the grid
			totalVirions += g.localVirions[i][j]
			totalDIPs += g.localDips[i][j]
		}
	}

	virionDiffusionRate := float64(totalVirionDiffusion) / float64(totalVirions)
	dipDiffusionRate := float64(totalDIPDiffusion) / float64(totalDIPs)
	return virionDiffusionRate, dipDiffusionRate
}

// Function to get the nth figure number in the folder
func getNextFigureNumber(outputFolder string) int {
	files, err := os.ReadDir(outputFolder)
	if err != nil {
		log.Fatalf("Failed to read output folder: %v", err)
	}
	count := 0
	for _, file := range files {
		if strings.HasSuffix(file.Name(), ".png") {
			count++
		}
	}
	return count + 1 // Return the next number
}

// Logic to determine IFN spreading type
func getIFNType() string {
	if IFN_wave_radius == 0 && globalIFN > 0 {
		return "Global"
	} else if IFN_wave_radius == 10 {
		return "IFNJ10"
	} else {
		return "NoIFN"
	}
}

func generateFolderName(
	no int,
	jumpRandomly bool,
	jumpRadiusD int,
	jumpRadiusV int,
	burstSizeD int,
	burstSizeV int,
	ifnWaveRadius int,
	TAU int,
	timeSteps int,
) string {
	// Determine Dinit naming part (keep at most 2 decimal places)
	dInit := fmt.Sprintf("Dinit%s", strconv.FormatFloat(*flag_d_pfu_initial, 'f', -1, 64))

	// Determine D naming part
	dName := ""
	if jumpRandomly {
		dName = fmt.Sprintf("DIPBst%d_JRand", burstSizeD)
	} else if jumpRadiusD > 0 {
		dName = fmt.Sprintf("DIPBst%d_J%d", burstSizeD, jumpRadiusD)
	} else if jumpRadiusD == 0 {
		dName = fmt.Sprintf("DIPBst%d_noJ", burstSizeD)
	} else {
		if burstSizeD == 0 && *flag_d_pfu_initial == 0 && D_only_IFN_stimulate_ratio == 0 && jumpRadiusD == 0 {
			dName = "NoDIP"
		} else {
			dName = fmt.Sprintf("DIPBst%d", burstSizeD)
		}
	}

	// Determine Vinit naming part (keep at most 2 decimal places)
	vInit := ""
	if *flag_v_pfu_initial > 0 {
		vInit = fmt.Sprintf("Vinit%s", strconv.FormatFloat(*flag_v_pfu_initial, 'f', -1, 64))
	} else if jumpRandomly {
		vInit = "JRand"
	} else if jumpRadiusV > 0 {
		vInit = fmt.Sprintf("J%d", jumpRadiusV)
	} else {
		vInit = "noJ"
	}

	vName := fmt.Sprintf("VBst%d", burstSizeV)

	// Determine IFN naming part
	ifnName := ""
	if TAU == 0 {
		ifnName = "NoIFN"
	} else if ifnWaveRadius == 0 {
		ifnName = "Global"
	} else {
		ifnName = fmt.Sprintf("IFN%d", ifnWaveRadius)
	}

	cellType := ""
	if TAU > 0 {
		cellType = "mdbk"
	} else {
		cellType = "vero"
	}

	folderName := fmt.Sprintf("%d_%s_%s_%s_%s_%s_%s_times%d_tau%d_ifnBothFold%.2f_grid%d_VStimulateIFN%t",
		no, dInit, dName, vInit, vName, ifnName, cellType, timeSteps, TAU, ifnBothFold, GRID_SIZE, VStimulateIFN)

	return folderName
}

// Combine images into one row
func combineImagesHorizontally(images []*image.RGBA) *image.RGBA {
	if len(images) == 0 {
		return nil
	}

	// Calculate the width and height of the combined image
	totalWidth := 0
	maxHeight := 0
	for _, img := range images {
		totalWidth += img.Bounds().Dx() // accumulate width
		if img.Bounds().Dy() > maxHeight {
			maxHeight = img.Bounds().Dy() // calculate maximum height
		}
	}

	// Create the combined image
	combinedImg := image.NewRGBA(image.Rect(0, 0, totalWidth, maxHeight))
	offsetX := 0
	for _, img := range images {
		rect := img.Bounds()
		draw.Draw(combinedImg, image.Rect(offsetX, 0, offsetX+rect.Dx(), rect.Dy()), img, rect.Min, draw.Src)
		offsetX += rect.Dx()
	}

	return combinedImg
}

// Save PNG image
func savePNGImage(img *image.RGBA, filename string) {
	file, err := os.Create(filename)
	if err != nil {
		log.Fatalf("Failed to create file %s: %v", filename, err)
	}
	defer file.Close()

	err = png.Encode(file, img)
	if err != nil {
		log.Fatalf("Failed to encode PNG: %v", err)
	}
}
func contains(arr []int, val int) bool {
	for _, v := range arr {
		if v == val {
			return true
		}
	}
	return false
}

// Function to calculate the maximum value from multiple datasets
func calculateMax(data ...[]float64) float64 {
	max := 0.0
	for _, dataset := range data {
		for _, value := range dataset {
			if value > max {
				max = value
			}
		}
	}
	return max
}
func LogValueFormatter(v interface{}) string {
	if value, ok := v.(float64); ok && value > 0 {
		return fmt.Sprintf("%.2f", math.Log10(value))
	}
	return "0"
}
func clampValues(data []float64, min, max float64) []float64 {
	clamped := make([]float64, len(data))
	for i, v := range data {
		if v < min {
			clamped[i] = min
		} else if v > max {
			clamped[i] = max
		} else {
			clamped[i] = v
		}
	}
	return clamped
}

// Modified function definition
func createInfectionGraph(frameNum int, virionOnly, dipOnly, both []float64, showLegend bool) *image.RGBA {
	graphWidth := GRID_SIZE * CELL_SIZE * 2
	graphHeight := 200

	if frameNum < 1 {
		log.Fatalf("Not enough data to render the graph: frameNum = %d", frameNum)
	}

	virionOnly = clampValues(virionOnly, 0.00, yMax)
	dipOnly = clampValues(dipOnly, 0.00, yMax)
	both = clampValues(both, 0.00, yMax)

	// Dynamically set legend name
	var series []chart.Series

	series = []chart.Series{
		chart.ContinuousSeries{
			Name:    "Infected by Virion Only",
			XValues: createTimeSeries(frameNum),
			YValues: virionOnly,
			Style:   chart.Style{StrokeColor: chart.ColorRed, StrokeWidth: 6.0},
		},
		chart.ContinuousSeries{
			Name:    "Infected by DIP Only",
			XValues: createTimeSeries(frameNum),
			YValues: dipOnly,
			Style:   chart.Style{StrokeColor: chart.ColorGreen, StrokeWidth: 6.0},
		},
		chart.ContinuousSeries{
			Name:    "Infected by Both",
			XValues: createTimeSeries(frameNum),
			YValues: both,
			Style:   chart.Style{StrokeColor: drawing.Color{R: 255, G: 165, B: 0, A: 255}, StrokeWidth: 8.0},
		},
	}

	graph := chart.Chart{
		Width:  459, // int(float64(GRID_SIZE*CELL_SIZE) * 1.51)
		Height: 100,
		XAxis: chart.XAxis{
			Style: chart.Style{FontSize: 10.0},
			ValueFormatter: func(v interface{}) string {
				return fmt.Sprintf("%d", int(v.(float64)))
			},
			Ticks: generateTicks(xMax, ticksInterval),
		},
		YAxis: chart.YAxis{
			Style: chart.Style{FontSize: 10.0},
		},
		Series: series,
	}

	buffer := bytes.NewBuffer([]byte{})
	err := graph.Render(chart.PNG, buffer)
	if err != nil {
		log.Printf("Failed to render graph: %v", err)
		// Return a simple colored rectangle instead of crashing
		return image.NewRGBA(image.Rect(0, 0, 459, 100))
	}

	graphImg, _, err := image.Decode(buffer)
	if err != nil {
		log.Fatalf("Failed to decode graph image: %v", err)
	}

	rgbaImg := image.NewRGBA(image.Rect(0, 0, graphWidth, graphHeight))
	draw.Draw(rgbaImg, rgbaImg.Bounds(), graphImg, image.Point{}, draw.Src)

	return rgbaImg
}

// saveCurrentGoFile saves the current Go source file into the specified output folder.
// saveCurrentGoFile saves the current Go source file with its original name and a timestamp.
func saveCurrentGoFile(outputFolder string) {
	_, currentFile, _, ok := runtime.Caller(0)
	if !ok {
		log.Println("Unable to get current Go file path")
		return
	}

	// Get the base name of the current Go file (without path)
	originalFileName := filepath.Base(currentFile)

	// Generate timestamp
	timestamp := time.Now().Format("20060102_150405") // Format: YYYYMMDD_HHMMSS

	// Target filename: original filename_timestamp.go
	newFileName := fmt.Sprintf("%s_%s.go", originalFileName[:len(originalFileName)-3], timestamp)
	outputFilePath := filepath.Join(outputFolder, newFileName)
	// Read Go file content
	content, err := ioutil.ReadFile(currentFile)
	if err != nil {
		log.Printf("cant read file %s: %v\n", currentFile, err)
		return
	}

	// Ensure target folder exists
	if err := os.MkdirAll(outputFolder, os.ModePerm); err != nil {
		log.Printf("cant make outputfolder %s: %v\n", outputFolder, err)
		return
	}

	// Write file
	err = ioutil.WriteFile(outputFilePath, content, 0644)
	if err != nil {
		log.Printf("cant save file to %s: %v\n", outputFilePath, err)
		return
	}

	log.Printf("file successfully saved in %s\n", outputFilePath)
}

func getNextFolderNumber(basePath string) int {
	files, err := os.ReadDir(basePath)
	if err != nil {
		log.Fatalf("Failed to read directory %s: %v", basePath, err)
	}

	maxNumber := 0
	for _, file := range files {
		if file.IsDir() {
			// Try to parse number from folder name
			var folderNumber int
			_, err := fmt.Sscanf(file.Name(), "%d", &folderNumber)
			if err == nil && folderNumber > maxNumber {
				maxNumber = folderNumber
			}
		}
	}
	return maxNumber + 1 // Return next available number
}

func transformToLogScale(data []float64) []float64 {
	transformed := make([]float64, len(data))
	for i, value := range data {
		if value > 0 {
			transformed[i] = math.Log10(value)
		} else {
			transformed[i] = math.Log10(0.0001) // handle log(0) case, use a very small value instead
		}
	}
	return transformed
}

func createTimeSeries(frameNum int) []float64 {
	if frameNum < 1 {
		return []float64{0, 1} // Return a default time series if not enough data
	}

	timeSeries := make([]float64, frameNum+1)
	for i := 0; i <= frameNum; i++ {
		timeSeries[i] = float64(i)
	}
	return timeSeries
}

// Function to calculate total virions in the grid
func (g *Grid) totalVirions() int {
	totalVirions := 0
	for i := 0; i < GRID_SIZE; i++ {
		for j := 0; j < GRID_SIZE; j++ {
			totalVirions += g.localVirions[i][j]
		}
	}
	return totalVirions
}

// Function to calculate total DIPs in the grid
func (g *Grid) totalDIPs() int {
	totalDIPs := 0
	for i := 0; i < GRID_SIZE; i++ {
		for j := 0; j < GRID_SIZE; j++ {
			totalDIPs += g.localDips[i][j]
		}
	}
	return totalDIPs
}

// Function to calculate the total number of regrowth cells in the grid
func (g *Grid) calculateRegrowthCount() int {
	regrowthCells := 0
	for i := 0; i < GRID_SIZE; i++ {
		for j := 0; j < GRID_SIZE; j++ {
			if g.state[i][j] == REGROWTH {
				regrowthCells++
				g.timeSinceRegrowth[i][j] += TIMESTEP
			}
		}
	}
	return regrowthCells
}

// Function to calculate the percentage of susceptible cells in the grid
func (g *Grid) calculateSusceptiblePercentage() float64 {
	totalCells := GRID_SIZE * GRID_SIZE
	susceptibleCells := 0
	for i := 0; i < GRID_SIZE; i++ {
		for j := 0; j < GRID_SIZE; j++ {
			if g.state[i][j] == SUSCEPTIBLE {
				susceptibleCells++
				g.timeSinceSusceptible[i][j] += TIMESTEP
			}
		}
	}
	return (float64(susceptibleCells) / float64(totalCells)) * 100
}

// Function to calculate the percentage of regrowthed or antiviral cells
func (g *Grid) calculateRegrowthedOrAntiviralPercentage() float64 {
	totalCells := GRID_SIZE * GRID_SIZE
	regrowthedOrAntiviralCells := 0
	for i := 0; i < GRID_SIZE; i++ {
		for j := 0; j < GRID_SIZE; j++ {
			if g.state[i][j] == REGROWTH || g.state[i][j] == ANTIVIRAL {
				regrowthedOrAntiviralCells++
			}
		}
	}
	return (float64(regrowthedOrAntiviralCells) / float64(totalCells)) * 100
}

// Function to calculate the percentage of infected cells (both virion and DIP infections)
func (g *Grid) calculateInfectedPercentage() float64 {
	totalCells := GRID_SIZE * GRID_SIZE
	infectedCells := 0
	for i := 0; i < GRID_SIZE; i++ {
		for j := 0; j < GRID_SIZE; j++ {
			if g.state[i][j] == INFECTED_VIRION || g.state[i][j] == INFECTED_DIP || g.state[i][j] == INFECTED_DIP_CONTINUOUS || g.state[i][j] == INFECTED_BOTH ||
				g.state[i][j] == INFECTED_VIRION_CONTINUOUS || g.state[i][j] == INFECTED_DIP || g.state[i][j] == INFECTED_DIP_CONTINUOUS || g.state[i][j] == INFECTED_BOTH_CONTINUOUS {
				infectedCells++
			}
		}
	}
	return (float64(infectedCells) / float64(totalCells)) * 100
}

// Function to calculate the percentage of DIP-only infected cells
func (g *Grid) calculateInfectedDIPOnlyPercentage() float64 {
	totalCells := GRID_SIZE * GRID_SIZE
	infectedDIPOnlyCells := 0
	for i := 0; i < GRID_SIZE; i++ {
		for j := 0; j < GRID_SIZE; j++ {
			if g.state[i][j] == INFECTED_DIP || g.state[i][j] == INFECTED_DIP_CONTINUOUS {
				infectedDIPOnlyCells++
			}
		}
	}
	return (float64(infectedDIPOnlyCells) / float64(totalCells)) * 100
}

// Function to calculate the percentage of cells infected by both virions and DIPs
func (g *Grid) calculateInfectedBothPercentage() float64 {
	totalCells := GRID_SIZE * GRID_SIZE
	infectedBothCells := 0
	for i := 0; i < GRID_SIZE; i++ {
		for j := 0; j < GRID_SIZE; j++ {
			if g.state[i][j] == INFECTED_BOTH || g.state[i][j] == INFECTED_BOTH_CONTINUOUS {
				infectedBothCells++
			}
		}
	}
	return (float64(infectedBothCells) / float64(totalCells)) * 100
}

// Function to calculate the percentage of antiviral cells (if antiviral state is modeled)
func (g *Grid) calculateAntiviralPercentage() float64 {
	totalCells := GRID_SIZE * GRID_SIZE
	antiviralCells := 0
	for i := 0; i < GRID_SIZE; i++ {
		for j := 0; j < GRID_SIZE; j++ {
			if g.state[i][j] == ANTIVIRAL {
				antiviralCells++
			}
		}
	}
	return (float64(antiviralCells) / float64(totalCells)) * 100
}

// Function to calculate the percentage of uninfected cells (susceptible and regrowth cells)
func (g *Grid) calculateUninfectedPercentage() float64 {
	totalCells := GRID_SIZE * GRID_SIZE
	uninfectedCells := 0
	for i := 0; i < GRID_SIZE; i++ {
		for j := 0; j < GRID_SIZE; j++ {
			if g.state[i][j] == SUSCEPTIBLE || g.state[i][j] == REGROWTH {
				uninfectedCells++
			}
		}
	}
	return (float64(uninfectedCells) / float64(totalCells)) * 100
}

// Function to calculate plaque percentage (for simplicity, counting dead cells as plaques)
func (g *Grid) calculatePlaquePercentage() float64 {
	totalCells := GRID_SIZE * GRID_SIZE
	plaqueCells := 0
	for i := 0; i < GRID_SIZE; i++ {
		for j := 0; j < GRID_SIZE; j++ {
			if g.state[i][j] == DEAD {
				plaqueCells++
			}
		}
	}
	return (float64(plaqueCells) / float64(totalCells)) * 100
}

// Function to calculate the percentage of dead cells
func calculateDeadCellPercentage(grid [GRID_SIZE][GRID_SIZE]int) float64 {
	totalCells := GRID_SIZE * GRID_SIZE
	deadCells := 0
	for i := 0; i < GRID_SIZE; i++ {
		for j := 0; j < GRID_SIZE; j++ {
			if grid[i][j] == DEAD {
				deadCells++
			}
		}
	}
	return (float64(deadCells) / float64(totalCells)) * 100
}

// Function to calculate the number of cells infected by virion only
func (g *Grid) calculateVirionOnlyInfected() int {
	virionOnlyInfected := 0
	for i := 0; i < GRID_SIZE; i++ {
		for j := 0; j < GRID_SIZE; j++ {
			if g.state[i][j] == INFECTED_VIRION || g.state[i][j] == INFECTED_VIRION_CONTINUOUS {
				virionOnlyInfected++
			}
		}
	}
	return virionOnlyInfected
}

// Function to calculate the number of cells infected by DIP only
func (g *Grid) calculateDipOnlyInfected() int {
	dipOnlyInfected := 0
	for i := 0; i < GRID_SIZE; i++ {
		for j := 0; j < GRID_SIZE; j++ {
			if g.state[i][j] == INFECTED_DIP || g.state[i][j] == INFECTED_DIP_CONTINUOUS {
				dipOnlyInfected++
			}
		}
	}
	return dipOnlyInfected
}

// Function to calculate the number of cells infected by both virion and DIP
func (g *Grid) calculateBothInfected() int {
	bothInfected := 0
	for i := 0; i < GRID_SIZE; i++ {
		for j := 0; j < GRID_SIZE; j++ {
			if g.state[i][j] == INFECTED_BOTH || g.state[i][j] == INFECTED_BOTH_CONTINUOUS {
				bothInfected++
			}
		}
	}
	return bothInfected
}

var precomputedRing [][2]int

func precomputeRing(radius int) [][2]int {
	var offsets [][2]int
	for dx := -radius; dx <= radius; dx++ {
		for dy := -radius; dy <= radius; dy++ {
			if dx*dx+dy*dy <= radius*radius {
				offsets = append(offsets, [2]int{dx, dy})
			}
		}
	}
	rand.Shuffle(len(offsets), func(i, j int) { offsets[i], offsets[j] = offsets[j], offsets[i] })
	return offsets
}

func precomputeIFNArea(radius int) [][2]int {
	var area [][2]int
	for di := -radius; di <= radius; di++ {
		for dj := -radius; dj <= radius; dj++ {
			distance := math.Sqrt(float64(di*di + dj*dj))
			// Include only cells within the radius
			if distance <= float64(radius) {
				area = append(area, [2]int{di, dj})
			}
		}
	}
	return area
}

// Add this new function, based on the competition mechanism from the paper

// Calculate neighbor relationships
func (g *Grid) initializeNeighbors() {
	// Initialize neighbors for all cells
	for i := 0; i < GRID_SIZE; i++ {
		for j := 0; j < GRID_SIZE; j++ {
			// Initialize fixed neighbor distances (1-10) using hexagonal neighbor calculation
			for radius := 1; radius <= 10; radius++ {
				neighbors := generateHexRing(i, j, radius)

				// Assign to appropriate neighbor array based on radius
				switch radius {
				case 1:
					copy(g.neighbors1[i][j][:], neighbors[:min(len(neighbors), 6)])
				case 2:
					copy(g.neighbors2[i][j][:], neighbors[:min(len(neighbors), 12)])
				case 3:
					copy(g.neighbors3[i][j][:], neighbors[:min(len(neighbors), 18)])
				case 4:
					copy(g.neighbors4[i][j][:], neighbors[:min(len(neighbors), 24)])
				case 5:
					copy(g.neighbors5[i][j][:], neighbors[:min(len(neighbors), 30)])
				case 6:
					copy(g.neighbors6[i][j][:], neighbors[:min(len(neighbors), 36)])
				case 7:
					copy(g.neighbors7[i][j][:], neighbors[:min(len(neighbors), 42)])
				case 8:
					copy(g.neighbors8[i][j][:], neighbors[:min(len(neighbors), 48)])
				case 9:
					copy(g.neighbors9[i][j][:], neighbors[:min(len(neighbors), 54)])
				case 10:
					copy(g.neighbors10[i][j][:], neighbors[:min(len(neighbors), 60)])
				}
			}

			// Initialize burst area neighbors based on configurable radius (virus and DIP use same radius)
			var burstAreaNeighbors [][2]int
			for radius := 1; radius <= g.burstRadius; radius++ {
				neighbors := generateHexRing(i, j, radius)
				for _, neighbor := range neighbors {
					if neighbor[0] >= 0 && neighbor[0] < GRID_SIZE && neighbor[1] >= 0 && neighbor[1] < GRID_SIZE {
						burstAreaNeighbors = append(burstAreaNeighbors, neighbor)
					}
				}
			}
			g.neighborsBurstArea[i][j] = burstAreaNeighbors

			// Initialize IFN area neighbors if enabled
			if ifnWave == true {
				precomputedIFNArea := precomputeIFNArea(IFN_wave_radius)
				var ifnAreaNeighbors [][2]int
				for _, offset := range precomputedIFNArea {
					newI, newJ := i+offset[0], j+offset[1]
					if newI >= 0 && newI < GRID_SIZE && newJ >= 0 && newJ < GRID_SIZE {
						ifnAreaNeighbors = append(ifnAreaNeighbors, [2]int{newI, newJ})
					}
				}
				g.neighborsIFNArea[i][j] = ifnAreaNeighbors
			}
		}
	}

	fmt.Println("Neighbors initialized")
}

// Generate neighbors in a hexagonal ring at specified radius
func generateHexRing(i, j, radius int) [][2]int {
	var neighbors [][2]int
	for dx := -radius; dx <= radius; dx++ {
		for dy := -radius; dy <= radius; dy++ {
			// Skip the center point
			if dx == 0 && dy == 0 {
				continue
			}
			// Calculate hexagonal distance
			hexDist := getHexDistance(dx, dy)
			// Only include neighbors at exactly the specified radius
			if hexDist == radius {
				neighbors = append(neighbors, [2]int{i + dx, j + dy})
			}
		}
	}
	return neighbors
}

// Calculate hexagonal distance between two points
func getHexDistance(dx, dy int) int {
	// Convert to cubic coordinate system to calculate true hexagonal distance
	q := dx
	r := dy - (dx+(dx&1))/2
	s := -q - r

	// Hexagonal distance is the maximum of the absolute values of the three cubic coordinates
	return max(abs(q), abs(r), abs(s))
}

// Calculate hexagonal distance between two points
func getHexDistanceBetweenPoints(x1, y1, x2, y2 int) int {
	dx := x2 - x1
	dy := y2 - y1
	return getHexDistance(dx, dy)
}

// Helper function for maximum of three integers
func max(a, b, c int) int {
	if a > b {
		if a > c {
			return a
		}
		return c
	}
	if b > c {
		return b
	}
	return c
}

// Helper function for absolute value
func abs(x int) int {
	if x < 0 {
		return -x
	}
	return x
}

// Helper function for minimum of two integers
func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

// Distance-weighted burst logic - "近多远少" (closer neighbors get more particles)
func (g *Grid) handleBurstWithConfigurableRadius(i, j int) {
	// Use neighborsBurstArea which contains all neighbors within burstRadius
	availableNeighbors := g.neighborsBurstArea[i][j]

	if len(availableNeighbors) == 0 {
		fmt.Printf("No available neighbors for burst at (%d,%d)\n", i, j)
		return
	}

	// Calculate distance-based weights for each neighbor
	// Weight formula: weight = 1 / (distance^1.5) - closer neighbors get exponentially more particles
	neighborWeights := make([]float64, len(availableNeighbors))
	totalWeight := 0.0

	for idx, neighbor := range availableNeighbors {
		ni, nj := neighbor[0], neighbor[1]
		distance := getHexDistanceBetweenPoints(i, j, ni, nj)
		if distance == 0 {
			distance = 1 // Avoid division by zero
		}
		// Distance-based weight: closer = higher weight
		weight := 1.0 / math.Pow(float64(distance), 1.5)
		neighborWeights[idx] = weight
		totalWeight += weight
	}

	// Determine burst sizes based on cell type and intracellular DIP content (following 0716 logic)
	burstSizeV := BURST_SIZE_V
	burstSizeD := 0 // Default: no DIPs

	// Calculate adjusted burst size for DIPs based on local ratio (like 0716 version)
	totalVirionsAtCell := g.localVirions[i][j]
	totalDIPsAtCell := g.localDips[i][j]
	adjustedBurstSizeD := 0

	if totalVirionsAtCell > 0 {
		dipVirionRatio := float64(totalDIPsAtCell) / float64(totalVirionsAtCell)
		adjustedBurstSizeD = BURST_SIZE_D + int(float64(BURST_SIZE_D)*dipVirionRatio)
	}

	// Adjust based on cell type (following 0716 logic - always use adjustedBurstSizeD)
	if g.state[i][j] == INFECTED_VIRION {
		// Virion-only infected cells: release virions only, no DVG/DIP release
		burstSizeV = BURST_SIZE_V
		if virionBurstMode == "both" {
			burstSizeD = adjustedBurstSizeD
			fmt.Printf("🔍 Cell (%d,%d) is INFECTED_VIRION, mode=both, will release %d DIPs\n", i, j, burstSizeD)
		} else {
			burstSizeD = 0
			fmt.Printf("🔍 Cell (%d,%d) is INFECTED_VIRION, mode=virionOnly, will release 0 DIPs\n", i, j)
		}
	} else if g.state[i][j] == INFECTED_BOTH {
		// Co-infected cells: produce both virions and DIPs
		burstSizeV = BURST_SIZE_V
		burstSizeD = adjustedBurstSizeD
	} else if g.state[i][j] == INFECTED_DIP || g.state[i][j] == INFECTED_DIP_CONTINUOUS {
		// DIP-only infected cells: no production (will be cleared)
		burstSizeV = 0
		burstSizeD = 0
	}

	fmt.Printf("🦠 Distance-weighted burst at (%d,%d): %d virions, %d DIPs to %d neighbors (radius=%d)\n",
		i, j, burstSizeV, burstSizeD, len(availableNeighbors), g.burstRadius)

	// Distribute particles based on distance weights
	virionsDistributed := 0
	dipsDistributed := 0

	for idx, neighbor := range availableNeighbors {
		ni, nj := neighbor[0], neighbor[1]
		distance := getHexDistanceBetweenPoints(i, j, ni, nj)

		// Calculate weighted distribution
		weight := neighborWeights[idx]
		proportion := weight / totalWeight

		// Distribute virions based on weight
		virionsToAdd := int(math.Round(float64(BURST_SIZE_V) * proportion))
		if virionsToAdd > 0 {
			g.localVirions[ni][nj] += virionsToAdd
			virionsDistributed += virionsToAdd
		}

		// Distribute DIPs based on weight
		dipsToAdd := int(math.Round(float64(burstSizeD) * proportion))
		if dipsToAdd > 0 {
			g.localDips[ni][nj] += dipsToAdd
			dipsDistributed += dipsToAdd
		}

		// Debug: Show distribution for close neighbors
		if distance <= 3 {
			fmt.Printf("  📍 Distance=%d → %d virions, %d DIPs (weight=%.3f)\n",
				distance, virionsToAdd, dipsToAdd, weight)
		}
	}

	// Handle any remaining particles due to rounding
	remainingVirions := burstSizeV - virionsDistributed
	remainingDips := burstSizeD - dipsDistributed

	// Distribute remaining particles to the closest neighbors (distance=1)
	for _, neighbor := range availableNeighbors {
		if remainingVirions <= 0 && remainingDips <= 0 {
			break
		}
		ni, nj := neighbor[0], neighbor[1]
		distance := getHexDistanceBetweenPoints(i, j, ni, nj)

		if distance == 1 { // Give remaining particles to immediate neighbors
			if remainingVirions > 0 {
				g.localVirions[ni][nj]++
				remainingVirions--
			}
			if remainingDips > 0 {
				g.localDips[ni][nj]++
				remainingDips--
			}
		}
	}

	fmt.Printf("  💊 Distributed: %d/%d virions, %d/%d DIPs\n",
		virionsDistributed+(burstSizeV-remainingVirions), burstSizeV,
		dipsDistributed+(burstSizeD-remainingDips), burstSizeD)
}

// Continuous production logic for Case 4
func (g *Grid) handleContinuousProduction(i, j, frameNum int) {
	// Only handle continuous mode states, skip burst mode states
	if g.state[i][j] != INFECTED_VIRION_CONTINUOUS &&
		g.state[i][j] != INFECTED_DIP_CONTINUOUS &&
		g.state[i][j] != INFECTED_BOTH_CONTINUOUS {
		return // Skip burst mode states
	}

	fmt.Printf("🔍 handleContinuousProduction called for cell (%d,%d) with state %d at frame %d\n", i, j, g.state[i][j], frameNum)

	// Check if cell is mature enough to start producing
	if !g.isProducing[i][j] {
		if frameNum-g.infectionTime[i][j] >= g.continuousIncubationPeriod {
			g.isProducing[i][j] = true
			fmt.Printf("🌱 Cell (%d,%d) matured and started continuous production at frame %d\n", i, j, frameNum)
		} else {
			return // Not yet mature
		}
	}

	// Check if it's time for lysis (if lysis time is set)
	if g.continuousLysisTime > 0 {
		if float64(frameNum-g.infectionTime[i][j]) >= g.continuousLysisTime {
			// Time for lysis - transition to DEAD state
			g.state[i][j] = DEAD
			g.stateChanged[i][j] = true
			g.isProducing[i][j] = false
			fmt.Printf("💀 Continuous production cell (%d,%d) lysed after %.1f hours\n", i, j, g.continuousLysisTime)
			return
		}
	}

	// Continuous production: release particles every timestep
	virionsToRelease := g.continuousProductionRateV
	dipsToRelease := g.continuousProductionRateD

	// Adjust production based on cell type and intracellular virus counts
	switch g.state[i][j] {
	case INFECTED_VIRION_CONTINUOUS:
		// Virion-only infected cells: continuous production respects virionBurstMode
		if g.intraWT[i][j] > 0 {
			virionsToRelease = int(float64(virionsToRelease) * float64(g.intraWT[i][j]))
		} else {
			virionsToRelease = 0
		}
		if virionBurstMode == "both" {
			if g.intraDVG[i][j] > 0 {
				dipsToRelease = int(float64(dipsToRelease) * float64(g.intraDVG[i][j]))
			} else {
				dipsToRelease = g.continuousProductionRateD
			}
		} else {
			dipsToRelease = 0
		}
	case INFECTED_DIP_CONTINUOUS:
		// DIP-only infected cells: NO PRODUCTION - they will be cleared back to susceptible
		// This follows the same mechanism as burst mode (handleDipOnlyClearance)
		virionsToRelease = 0
		dipsToRelease = 0
		return // Exit early - no production for DIP-only cells
	case INFECTED_BOTH_CONTINUOUS:
		// Co-infected cells: produce both virions and DIPs using base production rates
		// Scale production by intracellular virus counts
		if g.intraWT[i][j] > 0 {
			virionsToRelease = int(float64(virionsToRelease) * float64(g.intraWT[i][j]))
		} else {
			virionsToRelease = 0 // No virion production if no intracellular virions
		}
		if g.intraDVG[i][j] > 0 {
			dipsToRelease = int(float64(dipsToRelease) * float64(g.intraDVG[i][j]))
		} else {
			dipsToRelease = 0 // No DIP production if no intracellular DIPs
		}
	}

	fmt.Printf("🔄 Continuous production at (%d,%d): %d virions, %d DIPs (intraWT=%d, intraDVG=%d, state=%d, frame %d)\n",
		i, j, virionsToRelease, dipsToRelease, g.intraWT[i][j], g.intraDVG[i][j], g.state[i][j], frameNum)

	// Use the same distance-weighted distribution as burst mode
	g.distributeContinuousParticles(i, j, virionsToRelease, dipsToRelease)

	// Update intracellular virus counts based on production (simplified replication model)
	if g.continuousMode {
		switch g.state[i][j] {
		case INFECTED_VIRION_CONTINUOUS:
			// Virion replication: increase intracellular wild-type virus count
			if g.intraWT[i][j] > 0 {
				g.intraWT[i][j] += 1 // Simple replication model
			}
		case INFECTED_BOTH_CONTINUOUS:
			// Co-infection: both virus types replicate
			if g.intraWT[i][j] > 0 {
				g.intraWT[i][j] += 1 // Wild-type virus replication
			}
			if g.intraDVG[i][j] > 0 {
				g.intraDVG[i][j] += 1 // DVG replication
			}
		}
	}
}

// Distribute particles using continuous production with distance weights
func (g *Grid) distributeContinuousParticles(i, j, virions, dips int) {
	availableNeighbors := g.neighborsBurstArea[i][j]

	if len(availableNeighbors) == 0 {
		fmt.Printf("No available neighbors for continuous production at (%d,%d)\n", i, j)
		return
	}

	// Calculate distance-based weights (same as burst mode)
	neighborWeights := make([]float64, len(availableNeighbors))
	totalWeight := 0.0

	for idx, neighbor := range availableNeighbors {
		ni, nj := neighbor[0], neighbor[1]
		distance := getHexDistanceBetweenPoints(i, j, ni, nj)
		if distance == 0 {
			distance = 1
		}
		weight := 1.0 / math.Pow(float64(distance), 1.5)
		neighborWeights[idx] = weight
		totalWeight += weight
	}

	// Distribute particles based on distance weights
	for idx, neighbor := range availableNeighbors {
		ni, nj := neighbor[0], neighbor[1]
		weight := neighborWeights[idx]
		proportion := weight / totalWeight

		// Distribute virions
		virionsToAdd := int(math.Round(float64(virions) * proportion))
		if virionsToAdd > 0 {
			g.localVirions[ni][nj] += virionsToAdd
		}

		// Distribute DIPs
		dipsToAdd := int(math.Round(float64(dips) * proportion))
		if dipsToAdd > 0 {
			g.localDips[ni][nj] += dipsToAdd
		}
	}
}

// Handle burst or continuous production based on Case 4 mode
func (g *Grid) handleViralProduction(i, j, frameNum int) {
	// Check if this is Case 4 and continuous mode is enabled
	if g.initOption == 4 && g.continuousMode {
		// Use continuous production mode
		fmt.Printf("🔧 handleViralProduction: Case 4 continuous mode enabled, calling handleContinuousProduction\n")
		g.handleContinuousProduction(i, j, frameNum)
	} else {
		// Use traditional burst mode (all cases including Case 4 burst mode)
		fmt.Printf("🔧 handleViralProduction: Using burst mode (initOption=%d, continuousMode=%t), calling handleCase4Burst for cell (%d,%d)\n", g.initOption, g.continuousMode, i, j)
		g.handleCase4Burst(i, j, BURST_SIZE_V, BURST_SIZE_D, k_JumpR)
	}
}

// Handle Case 4 burst with 0716 logic (transplanted from 0716 version)
func (g *Grid) handleCase4Burst(i, j, burstSizeV, burstSizeD int, kJumpR float64) {
	// Calculate adjusted burst size for DIPs based on local ratio (like 0716 version)
	totalVirionsAtCell := g.localVirions[i][j]
	totalDIPsAtCell := g.localDips[i][j]
	adjustedBurstSizeD := 0

	if totalVirionsAtCell > 0 {
		dipVirionRatio := float64(totalDIPsAtCell) / float64(totalVirionsAtCell)
		adjustedBurstSizeD = burstSizeD + int(float64(burstSizeD)*dipVirionRatio)
	}

	// Under hotspot restriction, prevent DIPs from virion-only bursts when no DIPs present
	if *flag_dipRestrictToHotspot && totalDIPsAtCell == 0 && g.state[i][j] == INFECTED_VIRION {
		adjustedBurstSizeD = 0
	}

	// Guard: DIP-only cells should not release any virions (and we also set DIPs to 0 here; DIP-only clears elsewhere)
	if g.state[i][j] == INFECTED_DIP || g.state[i][j] == INFECTED_DIP_CONTINUOUS {
		burstSizeV = 0
		adjustedBurstSizeD = 0
	}

	// If restricted, forbid any DIP release during bursts; else apply virion-only rule
	if *flag_dipRestrictToHotspot {
		// Under restriction, only virion-only bursts are blocked from releasing DIPs
		if g.state[i][j] == INFECTED_VIRION {
			adjustedBurstSizeD = 0
		}
	} else {
		// Virion-only infected cells: DIP release depends on virionBurstMode
		if g.state[i][j] == INFECTED_VIRION {
			if virionBurstMode == "virionOnly" {
				adjustedBurstSizeD = 0
			}
		}
	}

	// Enforce VIRION-only pre-lysis cells release 0 DIPs even after state changes to DEAD
	// Use previousStates captured right before death to infer the infection type at lysis
	preBurstState := g.previousStates[i][j]
	if preBurstState == INFECTED_VIRION {
		adjustedBurstSizeD = 0
	}
	// NOTE: If preBurstState == INFECTED_BOTH, keep normal DIP release (burstSizeD, adjustedBurstSizeD)
	if preBurstState == INFECTED_BOTH {
		if burstSizeD > 0 {
			if adjustedBurstSizeD <= 0 {
				fmt.Printf("VIOLATION_PRE_BURST_BOTH: at (%d,%d) burstSizeD=%d adjustedBurstSizeD=%d (expected >0)\n",
					i, j, burstSizeD, adjustedBurstSizeD)
				// Fallback: ensure minimum DIP release equals BURST_SIZE_D when BOTH
				adjustedBurstSizeD = burstSizeD
			} else {
				fmt.Printf("ASSERT_PRE_BURST_BOTH: at (%d,%d) burstSizeD=%d adjustedBurstSizeD=%d OK\n",
					i, j, burstSizeD, adjustedBurstSizeD)
			}
		}
	}

	// DEBUG: log state and adjustedBurstSizeD at burst time (case 4)
	fmt.Printf("DEBUG Burst state=%d at (%d,%d): burstSizeV=%d, adjustedBurstSizeD=%d, virionBurstMode=%s\n",
		g.state[i][j], i, j, burstSizeV, adjustedBurstSizeD, virionBurstMode)

	// Get ALL neighbors from radius 1 to burstRadius (supports >10 via dynamic ring generation)
	var neighbors [][2]int
	var neighborsForDIP [][2]int
	radius := g.burstRadius
	radiusForDIP := *flag_dipRadius // DIP uses its own absolute radius
	// Enforce and assert: when pre-lysis state is BOTH, use dipRadius for DIP burst radius
	if preBurstState == INFECTED_BOTH {
		if radiusForDIP != *flag_dipRadius {
			radiusForDIP = *flag_dipRadius
		}
		fmt.Printf("ASSERT_BOTH_DIP_RADIUS: at (%d,%d) radiusD=%d dipRadiusFlag=%d\n", i, j, radiusForDIP, *flag_dipRadius)
	}
	if radius < 1 {
		radius = 1
	}
	// Allow up to 30; >10 rings will be generated dynamically below
	if radius > 30 {
		radius = 30
	}
	if radiusForDIP < 1 {
		radiusForDIP = 1
	}
	// Removed the upper limit restriction for DIP radius

	// Collect ALL neighbors from radius 1 to the specified radius for virions (using precomputed arrays)
	for r := 1; r <= radius; r++ {
		switch r {
		case 1:
			for n := 0; n < 6; n++ {
				neighbor := g.neighbors1[i][j][n]
				if neighbor[0] >= 0 && neighbor[0] < GRID_SIZE && neighbor[1] >= 0 && neighbor[1] < GRID_SIZE {
					neighbors = append(neighbors, neighbor)
				}
			}
		case 2:
			for n := 0; n < 12; n++ {
				neighbor := g.neighbors2[i][j][n]
				if neighbor[0] >= 0 && neighbor[0] < GRID_SIZE && neighbor[1] >= 0 && neighbor[1] < GRID_SIZE {
					neighbors = append(neighbors, neighbor)
				}
			}
		case 3:
			for n := 0; n < 18; n++ {
				neighbor := g.neighbors3[i][j][n]
				if neighbor[0] >= 0 && neighbor[0] < GRID_SIZE && neighbor[1] >= 0 && neighbor[1] < GRID_SIZE {
					neighbors = append(neighbors, neighbor)
				}
			}
		case 4:
			for n := 0; n < 24; n++ {
				neighbor := g.neighbors4[i][j][n]
				if neighbor[0] >= 0 && neighbor[0] < GRID_SIZE && neighbor[1] >= 0 && neighbor[1] < GRID_SIZE {
					neighbors = append(neighbors, neighbor)
				}
			}
		case 5:
			for n := 0; n < 30; n++ {
				neighbor := g.neighbors5[i][j][n]
				if neighbor[0] >= 0 && neighbor[0] < GRID_SIZE && neighbor[1] >= 0 && neighbor[1] < GRID_SIZE {
					neighbors = append(neighbors, neighbor)
				}
			}
		case 6:
			for n := 0; n < 36; n++ {
				neighbor := g.neighbors6[i][j][n]
				if neighbor[0] >= 0 && neighbor[0] < GRID_SIZE && neighbor[1] >= 0 && neighbor[1] < GRID_SIZE {
					neighbors = append(neighbors, neighbor)
				}
			}
		case 7:
			for n := 0; n < 42; n++ {
				neighbor := g.neighbors7[i][j][n]
				if neighbor[0] >= 0 && neighbor[0] < GRID_SIZE && neighbor[1] >= 0 && neighbor[1] < GRID_SIZE {
					neighbors = append(neighbors, neighbor)
				}
			}
		case 8:
			for n := 0; n < 48; n++ {
				neighbor := g.neighbors8[i][j][n]
				if neighbor[0] >= 0 && neighbor[0] < GRID_SIZE && neighbor[1] >= 0 && neighbor[1] < GRID_SIZE {
					neighbors = append(neighbors, neighbor)
				}
			}
		case 9:
			for n := 0; n < 54; n++ {
				neighbor := g.neighbors9[i][j][n]
				if neighbor[0] >= 0 && neighbor[0] < GRID_SIZE && neighbor[1] >= 0 && neighbor[1] < GRID_SIZE {
					neighbors = append(neighbors, neighbor)
				}
			}
		case 10:
			for n := 0; n < 60; n++ {
				neighbor := g.neighbors10[i][j][n]
				if neighbor[0] >= 0 && neighbor[0] < GRID_SIZE && neighbor[1] >= 0 && neighbor[1] < GRID_SIZE {
					neighbors = append(neighbors, neighbor)
				}
			}
		default:
			// r > 10: dynamically generate exact ring r
			for dx := -r; dx <= r; dx++ {
				for dy := -r; dy <= r; dy++ {
					if dx == 0 && dy == 0 {
						continue
					}
					ni, nj := i+dx, j+dy
					if ni >= 0 && ni < GRID_SIZE && nj >= 0 && nj < GRID_SIZE {
						if getHexDistance(dx, dy) == r {
							neighbors = append(neighbors, [2]int{ni, nj})
						}
					}
				}
			}
		}
	}

	// Collect DIP neighbors within a circular radius (no sector/annular ring)
	for r := 1; r <= radiusForDIP; r++ {
		switch r {
		case 1:
			for n := 0; n < 6; n++ {
				neighbor := g.neighbors1[i][j][n]
				if neighbor[0] >= 0 && neighbor[0] < GRID_SIZE && neighbor[1] >= 0 && neighbor[1] < GRID_SIZE {
					neighborsForDIP = append(neighborsForDIP, neighbor)
				}
			}
		case 2:
			for n := 0; n < 12; n++ {
				neighbor := g.neighbors2[i][j][n]
				if neighbor[0] >= 0 && neighbor[0] < GRID_SIZE && neighbor[1] >= 0 && neighbor[1] < GRID_SIZE {
					neighborsForDIP = append(neighborsForDIP, neighbor)
				}
			}
		case 3:
			for n := 0; n < 18; n++ {
				neighbor := g.neighbors3[i][j][n]
				if neighbor[0] >= 0 && neighbor[0] < GRID_SIZE && neighbor[1] >= 0 && neighbor[1] < GRID_SIZE {
					neighborsForDIP = append(neighborsForDIP, neighbor)
				}
			}
		case 4:
			for n := 0; n < 24; n++ {
				neighbor := g.neighbors4[i][j][n]
				if neighbor[0] >= 0 && neighbor[0] < GRID_SIZE && neighbor[1] >= 0 && neighbor[1] < GRID_SIZE {
					neighborsForDIP = append(neighborsForDIP, neighbor)
				}
			}
		case 5:
			for n := 0; n < 30; n++ {
				neighbor := g.neighbors5[i][j][n]
				if neighbor[0] >= 0 && neighbor[0] < GRID_SIZE && neighbor[1] >= 0 && neighbor[1] < GRID_SIZE {
					neighborsForDIP = append(neighborsForDIP, neighbor)
				}
			}
		case 6:
			for n := 0; n < 36; n++ {
				neighbor := g.neighbors6[i][j][n]
				if neighbor[0] >= 0 && neighbor[0] < GRID_SIZE && neighbor[1] >= 0 && neighbor[1] < GRID_SIZE {
					neighborsForDIP = append(neighborsForDIP, neighbor)
				}
			}
		case 7:
			for n := 0; n < 42; n++ {
				neighbor := g.neighbors7[i][j][n]
				if neighbor[0] >= 0 && neighbor[0] < GRID_SIZE && neighbor[1] >= 0 && neighbor[1] < GRID_SIZE {
					neighborsForDIP = append(neighborsForDIP, neighbor)
				}
			}
		case 8:
			for n := 0; n < 48; n++ {
				neighbor := g.neighbors8[i][j][n]
				if neighbor[0] >= 0 && neighbor[0] < GRID_SIZE && neighbor[1] >= 0 && neighbor[1] < GRID_SIZE {
					neighborsForDIP = append(neighborsForDIP, neighbor)
				}
			}
		case 9:
			for n := 0; n < 54; n++ {
				neighbor := g.neighbors9[i][j][n]
				if neighbor[0] >= 0 && neighbor[0] < GRID_SIZE && neighbor[1] >= 0 && neighbor[1] < GRID_SIZE {
					neighborsForDIP = append(neighborsForDIP, neighbor)
				}
			}
		case 10:
			for n := 0; n < 60; n++ {
				neighbor := g.neighbors10[i][j][n]
				if neighbor[0] >= 0 && neighbor[0] < GRID_SIZE && neighbor[1] >= 0 && neighbor[1] < GRID_SIZE {
					neighborsForDIP = append(neighborsForDIP, neighbor)
				}
			}
		default:
			for dx := -r; dx <= r; dx++ {
				for dy := -r; dy <= r; dy++ {
					if dx == 0 && dy == 0 {
						continue
					}
					ni, nj := i+dx, j+dy
					if ni >= 0 && ni < GRID_SIZE && nj >= 0 && nj < GRID_SIZE {
						if getHexDistance(dx, dy) == r {
							neighborsForDIP = append(neighborsForDIP, [2]int{ni, nj})
						}
					}
				}
			}
		}
	}

	fmt.Printf("Case 4 burst at [%d][%d] with radiusV=%d, radiusD=%d, using %d virion neighbors, %d DIP neighbors, burstSizeV=%d, adjustedBurstSizeD=%d\n",
		i, j, radius, radiusForDIP, len(neighbors), len(neighborsForDIP), burstSizeV, adjustedBurstSizeD)

	// Distribute virions using original radius
	if len(neighbors) > 0 {
		// Group neighbors by distance for weighted distribution
		neighborsByDistance := make(map[float64][][2]int)
		for _, neighbor := range neighbors {
			ni, nj := neighbor[0], neighbor[1]
			if ni >= 0 && ni < GRID_SIZE && nj >= 0 && nj < GRID_SIZE {
				dx := float64(ni - i)
				dy := float64(nj - j)
				distance := math.Sqrt(dx*dx + dy*dy)
				neighborsByDistance[distance] = append(neighborsByDistance[distance], neighbor)
			}
		}

		// Calculate total weight (inverse distance weighting)
		totalWeight := 0.0
		for distance := range neighborsByDistance {
			weight := 1.0 / (distance + 0.1) // Closer neighbors get higher weight
			totalWeight += weight * float64(len(neighborsByDistance[distance]))
		}

		// Distribute virions by distance
		virionsByDistance := make(map[float64]int)
		for distance := range neighborsByDistance {
			weight := 1.0 / (distance + 0.1)
			virionsForDistance := int(math.Floor(float64(burstSizeV) *
				(weight * float64(len(neighborsByDistance[distance]))) / totalWeight))
			virionsByDistance[distance] = virionsForDistance
		}

		// Distribute virions to neighbors
		for distance, neighborsAtDistance := range neighborsByDistance {
			virionsForThisDistance := virionsByDistance[distance]

			if len(neighborsAtDistance) > 0 {
				// Shuffle order within this ring to avoid directional bias
				shuffled := make([][2]int, len(neighborsAtDistance))
				copy(shuffled, neighborsAtDistance)
				rand.Shuffle(len(shuffled), func(a, b int) { shuffled[a], shuffled[b] = shuffled[b], shuffled[a] })

				virionsPerNeighbor := virionsForThisDistance / len(shuffled)
				remainingVirions := virionsForThisDistance % len(shuffled)

				for idx, neighbor := range shuffled {
					ni, nj := neighbor[0], neighbor[1]
					virionsToAdd := virionsPerNeighbor

					if idx < remainingVirions {
						virionsToAdd++
					}

					g.localVirions[ni][nj] += virionsToAdd
				}
			}
		}
	}

	// Distribute DIPs to neighbors with the SAME distance-weighted strategy as virions
	if len(neighborsForDIP) > 0 {
		// Group DIP neighbors by distance
		dipNeighborsByDistance := make(map[float64][][2]int)
		for _, neighbor := range neighborsForDIP {
			ni, nj := neighbor[0], neighbor[1]
			if ni >= 0 && ni < GRID_SIZE && nj >= 0 && nj < GRID_SIZE {
				dx := float64(ni - i)
				dy := float64(nj - j)
				distance := math.Sqrt(dx*dx + dy*dy)
				dipNeighborsByDistance[distance] = append(dipNeighborsByDistance[distance], neighbor)
			}
		}

		// Calculate total weight (inverse distance weighting)
		dipTotalWeight := 0.0
		for distance := range dipNeighborsByDistance {
			weight := 1.0 / (distance + 0.1)
			dipTotalWeight += weight * float64(len(dipNeighborsByDistance[distance]))
		}

		// Determine DIPs per distance bucket
		dipsByDistance := make(map[float64]int)
		for distance := range dipNeighborsByDistance {
			weight := 1.0 / (distance + 0.1)
			dipsForDistance := int(math.Floor(float64(adjustedBurstSizeD) *
				(weight * float64(len(dipNeighborsByDistance[distance]))) / dipTotalWeight))
			dipsByDistance[distance] = dipsForDistance
		}

		// Distribute DIPs to neighbors within each distance bucket
		distributedDIPs := 0
		for distance, neighborsAtDistance := range dipNeighborsByDistance {
			dipsForThisDistance := dipsByDistance[distance]
			if len(neighborsAtDistance) > 0 {
				// Shuffle order within this ring to avoid directional bias
				shuffled := make([][2]int, len(neighborsAtDistance))
				copy(shuffled, neighborsAtDistance)
				rand.Shuffle(len(shuffled), func(a, b int) { shuffled[a], shuffled[b] = shuffled[b], shuffled[a] })

				dipsPerNeighbor := dipsForThisDistance / len(shuffled)
				remaining := dipsForThisDistance % len(shuffled)
				for idx, neighbor := range shuffled {
					ni, nj := neighbor[0], neighbor[1]
					dipsToAdd := dipsPerNeighbor
					if idx < remaining {
						dipsToAdd++
					}
					if dipsToAdd > 0 {
						g.localDips[ni][nj] += dipsToAdd
						distributedDIPs += dipsToAdd
					}
				}
			}
		}

		// Handle rounding leftovers: if any remain, give to closest ring (distance minimum)
		remainingDIPs := adjustedBurstSizeD - distributedDIPs
		if remainingDIPs > 0 {
			// find minimal distance bucket
			minDist := math.MaxFloat64
			for d := range dipNeighborsByDistance {
				if d < minDist {
					minDist = d
				}
			}
			if neighborsAtMin, ok := dipNeighborsByDistance[minDist]; ok {
				// Randomize starting index to avoid fixed-direction bias
				start := 0
				if len(neighborsAtMin) > 1 {
					start = rand.Intn(len(neighborsAtMin))
				}
				idx := 0
				for remainingDIPs > 0 && len(neighborsAtMin) > 0 {
					spot := neighborsAtMin[(start+idx)%len(neighborsAtMin)]
					ni, nj := spot[0], spot[1]
					g.localDips[ni][nj]++
					remainingDIPs--
					idx++
				}
			}
		}
	}
	fmt.Printf("Case 4 burst completed - distributed virions to %d neighbors, DIPs to %d neighbors\n", len(neighbors), len(neighborsForDIP))
}

// Helper function to clear viral particles from dead cell locations
func (g *Grid) clearParticlesFromDeadCells() {
	for i := 0; i < GRID_SIZE; i++ {
		for j := 0; j < GRID_SIZE; j++ {
			if g.state[i][j] == DEAD {
				// Clear all extracellular virions and DIPs from dead cell locations
				g.localVirions[i][j] = 0
				g.localDips[i][j] = 0
			}
		}
	}
}

// Generate DIP clearance time using normal distribution (mean=2, std=1)
func (g *Grid) generateDipClearanceTime() int {
	// Generate time using normal distribution with mean=2, std=1
	clearanceTime := int(rand.NormFloat64()*1.0 + 2.0)
	// Ensure minimum clearance time of 1 hour
	if clearanceTime < 1 {
		clearanceTime = 1
	}
	return clearanceTime
}

// Handle DIP-only infected cells clearance (become susceptible after mean=2±1 hours if still DIP-only)
func (g *Grid) handleDipOnlyClearance(frameNum int) {
	dipOnlyClearedCount := 0

	for i := 0; i < GRID_SIZE; i++ {
		for j := 0; j < GRID_SIZE; j++ {
			// Handle both burst mode (INFECTED_DIP) and continuous mode (INFECTED_DIP_CONTINUOUS)
			if g.state[i][j] == INFECTED_DIP || g.state[i][j] == INFECTED_DIP_CONTINUOUS {
				// Set clearance threshold if not already set
				if g.dipClearanceThreshold[i][j] == -1 {
					g.dipClearanceThreshold[i][j] = g.generateDipClearanceTime()
				}

				// Check if DIP clearance time has been reached
				if g.timeSinceInfectDIP[i][j] >= g.dipClearanceThreshold[i][j] {
					// Clear DIP-only infected cell back to susceptible
					g.state[i][j] = SUSCEPTIBLE
					g.timeSinceInfectDIP[i][j] = -1
					g.dipClearanceThreshold[i][j] = -1
					g.timeSinceSusceptible[i][j] = 0
					g.isProducing[i][j] = false // Reset continuous production flag
					dipOnlyClearedCount++
				}
			}
		}
	}

	if dipOnlyClearedCount > 0 {
		fmt.Printf("🔄 Frame %d: %d DIP-only infected cells cleared and became susceptible\n", frameNum, dipOnlyClearedCount)
	}
}

// Test function to verify that dead cells have no viral particles
func (g *Grid) testDeadCellParticleClearance(frameNum int) {
	deadCellsWithParticles := 0
	totalDeadCells := 0

	for i := 0; i < GRID_SIZE; i++ {
		for j := 0; j < GRID_SIZE; j++ {
			if g.state[i][j] == DEAD {
				totalDeadCells++
				if g.localVirions[i][j] > 0 || g.localDips[i][j] > 0 {
					deadCellsWithParticles++
					fmt.Printf("⚠️  Frame %d: Dead cell at (%d,%d) has %d virions and %d DIPs!\n",
						frameNum, i, j, g.localVirions[i][j], g.localDips[i][j])
				}
			}
		}
	}

	if deadCellsWithParticles == 0 && totalDeadCells > 0 {
		fmt.Printf("✅ Frame %d: All %d dead cells have 0 viral particles (test passed)\n", frameNum, totalDeadCells)
	} else if totalDeadCells == 0 {
		// No dead cells to test - this is normal in early frames
	} else {
		fmt.Printf("❌ Frame %d: %d out of %d dead cells still have viral particles (test failed)\n",
			frameNum, deadCellsWithParticles, totalDeadCells)
	}
}

// Update the state of the grid at each time step
func (g *Grid) update(frameNum int) {
	newGrid := g.state

	if ifnWave == true {
		for i := 0; i < GRID_SIZE; i++ {
			for j := 0; j < GRID_SIZE; j++ {
				g.stateChanged[i][j] = false

			}
		}

		// Step 3: Update max global IFN if needed
		if globalIFN < 0 {
			globalIFN = -1.0
		}
		if globalIFN > maxGlobalIFN {

			maxGlobalIFN = globalIFN

		}
		fmt.Printf("Global IFN concentration: %.2f\n", globalIFN)

		// Traverse the grid
		for i := 0; i < GRID_SIZE; i++ {
			for j := 0; j < GRID_SIZE; j++ {
				// Skip UNEXPOSED cells entirely (never change)
				if g.state[i][j] == UNEXPOSED {
					continue
				}
				// Only consider cells that are in the SUSCEPTIBLE or REGROWTH state

				var regional_sumIFN float64
				neighborsCount := len(g.neighborsIFNArea[i][j])

				if ifn_half_life != 0 {
					for i := 0; i < GRID_SIZE; i++ {
						for j := 0; j < GRID_SIZE; j++ {
							// Update IFN amount using half-life formula
							factorIFN := math.Pow(0.5, float64(TIMESTEP)/ifn_half_life)
							g.IFNConcentration[i][j] *= factorIFN
							// Remove IFN if concentration is below threshold
							if g.IFNConcentration[i][j] < (1.0 / (float64(GRID_SIZE) * float64(GRID_SIZE))) {
								g.IFNConcentration[i][j] = 0
							}
						}
					}
				}

				// Sum the IFN concentration within the IFN area
				for _, neighbor := range g.neighborsIFNArea[i][j] {
					ni, nj := neighbor[0], neighbor[1]

					regional_sumIFN += g.IFNConcentration[ni][nj]
				}

				// Calculate the average IFN concentration if there are neighbors within the radius
				var regionalAverageIFN float64
				if neighborsCount > 0 {
					regionalAverageIFN = regional_sumIFN / float64(neighborsCount)
				} else {
					regionalAverageIFN = 0 // Default to 0 if no neighbors, though this should rarely occur
				}

				if g.state[i][j] == SUSCEPTIBLE || g.state[i][j] == REGROWTH || g.state[i][j] == INFECTED_DIP || g.state[i][j] == INFECTED_DIP_CONTINUOUS {
					if g.IFNConcentration[i][j] > 0 && TAU > 0 {

						if g.antiviralDuration[i][j] <= -1 {
							g.antiviralDuration[i][j] = int(rand.NormFloat64()*float64(TAU)/4 + float64(TAU))
							g.timeSinceAntiviral[i][j] = 0
						} else if g.timeSinceAntiviral[i][j] <= int(g.antiviralDuration[i][j]) {
							g.timeSinceAntiviral[i][j] += TIMESTEP
						} else {

							g.previousStates[i][j] = g.state[i][j]
							newGrid[i][j] = ANTIVIRAL

							g.timeSinceAntiviral[i][j] = -2
							g.totalAntiviralTime += g.antiviralDuration[i][j]
							if g.state[i][j] == ANTIVIRAL && !g.antiviralFlag[i][j] {
								g.antiviralFlag[i][j] = true
								g.antiviralCellCount++
							}

						}
					}

					if g.state[i][j] == SUSCEPTIBLE || g.state[i][j] == REGROWTH {
						// Check if the cell is infected by virions or DIPs
						if g.localVirions[i][j] > 0 || g.localDips[i][j] > 0 {
							// Calculate the infection probabilities
							if R == 0 || TAU == 0 {
								perParticleInfectionChance_V = RHO
							} else if VStimulateIFN == true && R > 0 { // R=1
								perParticleInfectionChance_V = RHO * math.Exp(-ALPHA*(regionalAverageIFN/float64(R)))
							} else if !VStimulateIFN { // usually only DIP stimulate IFN in this situlation
								perParticleInfectionChance_V = RHO * math.Exp(-ALPHA*(regionalAverageIFN))
							}
							var probabilityVInfection, probabilityDInfection float64

							// Virion infection probability
							probabilityVInfection = 1 - math.Pow(1-perParticleInfectionChance_V, float64(g.localVirions[i][j]))
							infectedByVirion := rand.Float64() <= probabilityVInfection

							// DIP infection probability
							probabilityDInfection = 1 - math.Pow(1-(RHO*math.Exp(-ALPHA*(regionalAverageIFN))), float64(g.localDips[i][j]))
							infectedByDip := rand.Float64() <= probabilityDInfection

							// Determine the infection state based on virion and DIP infection
							if infectedByVirion && infectedByDip {
								if g.continuousMode {
									newGrid[i][j] = INFECTED_BOTH_CONTINUOUS
								} else {
									newGrid[i][j] = INFECTED_BOTH
								}
								g.timeSinceSusceptible[i][j] = -1
								g.timeSinceRegrowth[i][j] = -1
								// Record intracellular virus counts for continuous mode
								if g.continuousMode {
									// Calculate actual number of infecting virions and DIPs
									infectingVirions := int(math.Round(float64(g.localVirions[i][j]) * probabilityVInfection))
									infectingDIPs := int(math.Round(float64(g.localDips[i][j]) * probabilityDInfection))
									if infectingVirions > 0 {
										g.intraWT[i][j] += infectingVirions
									}
									if infectingDIPs > 0 {
										g.intraDVG[i][j] += infectingDIPs
									}
									g.infectionTime[i][j] = frameNum
								}
							} else if infectedByVirion {
								if g.continuousMode {
									newGrid[i][j] = INFECTED_VIRION_CONTINUOUS
								} else {
									newGrid[i][j] = INFECTED_VIRION
								}
								g.timeSinceSusceptible[i][j] = -1
								g.timeSinceRegrowth[i][j] = -1
								// Record intracellular virus count for continuous mode
								if g.continuousMode {
									// Calculate actual number of infecting virions
									infectingVirions := int(math.Round(float64(g.localVirions[i][j]) * probabilityVInfection))
									if infectingVirions > 0 {
										g.intraWT[i][j] += infectingVirions
									}
									g.infectionTime[i][j] = frameNum
								}
							} else if infectedByDip {
								if g.continuousMode {
									newGrid[i][j] = INFECTED_DIP_CONTINUOUS
								} else {
									newGrid[i][j] = INFECTED_DIP
								}
								g.timeSinceSusceptible[i][j] = -1
								g.timeSinceRegrowth[i][j] = -1
								// Record intracellular DVG count for continuous mode
								if g.continuousMode {
									// Calculate actual number of infecting DIPs
									infectingDIPs := int(math.Round(float64(g.localDips[i][j]) * probabilityDInfection))
									if infectingDIPs > 0 {
										g.intraDVG[i][j] += infectingDIPs
									}
									g.infectionTime[i][j] = frameNum
								}
							}
						}

						// Mark the state as changed if the cell is infected
						if newGrid[i][j] != g.state[i][j] {
							g.stateChanged[i][j] = true
						}
					}

				}

			}
		}

		// Process infected cells
		for i := 0; i < GRID_SIZE; i++ {
			for j := 0; j < GRID_SIZE; j++ {

				var regional_sumIFN float64

				// Sum the IFN concentration within the IFN area
				for _, neighbor := range g.neighborsIFNArea[i][j] {
					ni, nj := neighbor[0], neighbor[1]
					regional_sumIFN += g.IFNConcentration[ni][nj]
				}

				// Note: regionalAverageIFN is not used in this section

				if g.state[i][j] == INFECTED_VIRION || g.state[i][j] == INFECTED_DIP || g.state[i][j] == INFECTED_DIP_CONTINUOUS || g.state[i][j] == INFECTED_BOTH ||
					g.state[i][j] == INFECTED_VIRION_CONTINUOUS || g.state[i][j] == INFECTED_DIP || g.state[i][j] == INFECTED_DIP_CONTINUOUS || g.state[i][j] == INFECTED_BOTH_CONTINUOUS {
					fmt.Printf("🔍 DEBUG: Processing infected cell at (%d,%d) with state %d at frame %d\n", i, j, g.state[i][j], frameNum)

					// Handle burst mode cells (lysis logic)
					if g.state[i][j] == INFECTED_VIRION || g.state[i][j] == INFECTED_BOTH {
						if g.lysisThreshold[i][j] == -1 {
							g.lysisThreshold[i][j] = int(rand.NormFloat64()*STANDARD_LYSIS_TIME + MEAN_LYSIS_TIME)
						}
						g.timeSinceInfectVorBoth[i][j] += TIMESTEP
						g.timeSinceInfectDIP[i][j] = -1

						// Check if the cell should lyse and release virions and DIPs
						if g.lysisThreshold[i][j] > 0 && g.timeSinceInfectVorBoth[i][j] >= g.lysisThreshold[i][j] {

							// After lysis, the cell becomes DEAD and virions and DIPs are spread to neighbors
							if g.state[i][j] == INFECTED_VIRION {
								totalDeadFromV++ // Increase INFECTED_VIRION death count
							} else if g.state[i][j] == INFECTED_BOTH {
								totalDeadFromBoth++ // Increase INFECTED_BOTH death count
							}

							prevState := g.state[i][j]
							g.previousStates[i][j] = prevState // Save the previous state before death
							newGrid[i][j] = DEAD
							g.state[i][j] = DEAD
							g.timeSinceDead[i][j] = 0
							g.timeSinceInfectVorBoth[i][j] = -1
							g.timeSinceInfectDIP[i][j] = -1
							g.lysisThreshold[i][j] = -1

							///////////// for k_jumpR percent cells that jump reandomly
							if par_celltocell_random == true {
								// Calculate adjusted burst size for DIPs based on local ratio
								totalVirionsAtCell := g.localVirions[i][j]
								totalDIPsAtCell := g.localDips[i][j]
								adjustedBurstSizeD := BURST_SIZE_D
								if totalVirionsAtCell > 0 {
									dipVirionRatio := float64(totalDIPsAtCell) / float64(totalVirionsAtCell)
									adjustedBurstSizeD += int(float64(BURST_SIZE_D) * dipVirionRatio)
								}
								if prevState == INFECTED_VIRION {
									adjustedBurstSizeD = 0
								}
								//  ---------------------------------------
								// Partition mode: split particles between random jump and cell-to-cell
								randomVirions := int(math.Floor(float64(BURST_SIZE_V) * k_JumpR))
								virionsForLocalDiffusion := BURST_SIZE_V - randomVirions

								randomDIPs := int(math.Floor(float64(adjustedBurstSizeD) * k_JumpR))
								dipsForLocalDiffusion := adjustedBurstSizeD - randomDIPs

								// Handle random jumps
								for v := 0; v < randomVirions; v++ {
									ni, nj := rand.Intn(GRID_SIZE), rand.Intn(GRID_SIZE)
									g.localVirions[ni][nj]++
									g.totalRandomJumpVirions++
								}
								for d := 0; d < randomDIPs; d++ {
									ni, nj := rand.Intn(GRID_SIZE), rand.Intn(GRID_SIZE)
									g.localDips[ni][nj]++
									g.totalRandomJumpDIPs++
								}

								// Handle local diffusion
								// Handle local diffusion with localVirions & localDIPs (keep original logic unchanged)
								if virionsForLocalDiffusion > 0 || dipsForLocalDiffusion > 0 {
									// Calculate the total number of valid neighbors
									totalNeighbors := 0

									// Count valid neighbors from neighbors1
									for _, dir := range g.neighbors1[i][j] {
										if dir != [2]int{-1, -1} {
											totalNeighbors++
										}
									}
									// Count valid neighbors from neighbors2
									for _, dir := range g.neighbors2[i][j] {
										if dir != [2]int{-1, -1} {
											totalNeighbors++
										}
									}
									// Count valid neighbors from neighbors3
									for _, dir := range g.neighbors3[i][j] {
										if dir != [2]int{-1, -1} {
											totalNeighbors++
										}
									}

									if totalNeighbors == 0 {
										return
									}

									// Calculate the distribution based on the ratio √3 : 2√3 : 3
									sqrt3 := math.Sqrt(3)
									ratio1 := 1.0               // sqrt3     // Weight for neighbors1
									ratio2 := 1.0 / 2           // 2 * sqrt3 // Weight for neighbors2
									ratio3 := 1.0 / (3 / sqrt3) // 3.0       // Weight for neighbors3
									totalRatio := ratio1*float64(len(g.neighbors1[i][j])) +
										ratio2*float64(len(g.neighbors2[i][j])) +
										ratio3*float64(len(g.neighbors3[i][j]))

									// Calculate virions for each neighbor group
									virionsForNeighbors1 := int(math.Floor(float64(virionsForLocalDiffusion) * (ratio1 * float64(len(g.neighbors1[i][j]))) / totalRatio))
									virionsForNeighbors2 := int(math.Floor(float64(virionsForLocalDiffusion) * (ratio2 * float64(len(g.neighbors2[i][j]))) / totalRatio))
									virionsForNeighbors3 := int(math.Floor(float64(virionsForLocalDiffusion) * (ratio3 * float64(len(g.neighbors3[i][j]))) / totalRatio))

									// Calculate remaining virions
									remainingVirions := virionsForLocalDiffusion - (virionsForNeighbors1 + virionsForNeighbors2 + virionsForNeighbors3)

									// Distribute remaining virions based on ratio
									for remainingVirions > 0 {
										randVal := rand.Float64() * totalRatio
										if randVal < ratio1 && len(g.neighbors1[i][j]) > 0 {
											virionsForNeighbors1++
										} else if randVal < (ratio1+ratio2) && len(g.neighbors2[i][j]) > 0 {
											virionsForNeighbors2++
										} else if len(g.neighbors3[i][j]) > 0 {
											virionsForNeighbors3++
										}
										remainingVirions--
									}

									// Calculate DIPs for each neighbor group (same logic)
									dipsForNeighbors1 := int(math.Floor(float64(dipsForLocalDiffusion) * (ratio1 * float64(len(g.neighbors1[i][j]))) / totalRatio))
									dipsForNeighbors2 := int(math.Floor(float64(dipsForLocalDiffusion) * (ratio2 * float64(len(g.neighbors2[i][j]))) / totalRatio))
									dipsForNeighbors3 := int(math.Floor(float64(dipsForLocalDiffusion) * (ratio3 * float64(len(g.neighbors3[i][j]))) / totalRatio))

									remainingDIPs := dipsForLocalDiffusion - (dipsForNeighbors1 + dipsForNeighbors2 + dipsForNeighbors3)

									for remainingDIPs > 0 {
										randVal := rand.Float64() * totalRatio
										if randVal < ratio1 && len(g.neighbors1[i][j]) > 0 {
											dipsForNeighbors1++
										} else if randVal < (ratio1+ratio2) && len(g.neighbors2[i][j]) > 0 {
											dipsForNeighbors2++
										} else if len(g.neighbors3[i][j]) > 0 {
											dipsForNeighbors3++
										}
										remainingDIPs--
									}

									// Distribute virions to neighbors1
									for _, dir := range g.neighbors1[i][j] {
										ni, nj := dir[0], dir[1]
										if dir != [2]int{-1, -1} && ni >= 0 && ni < GRID_SIZE && nj >= 0 && nj < GRID_SIZE {
											if g.state[ni][nj] == SUSCEPTIBLE {
												g.localVirions[ni][nj] += virionsForNeighbors1 / len(g.neighbors1[i][j])
												g.localDips[ni][nj] += dipsForNeighbors1 / len(g.neighbors1[i][j])
											}
										}
									}

									// Distribute virions to neighbors2
									for _, dir := range g.neighbors2[i][j] {
										ni, nj := dir[0], dir[1]
										if dir != [2]int{-1, -1} && ni >= 0 && ni < GRID_SIZE && nj >= 0 && nj < GRID_SIZE {
											g.localVirions[ni][nj] += virionsForNeighbors2 / len(g.neighbors2[i][j])
											g.localDips[ni][nj] += dipsForNeighbors2 / len(g.neighbors2[i][j])
										}
									}

									// Distribute virions to neighbors3
									for _, dir := range g.neighbors3[i][j] {
										ni, nj := dir[0], dir[1]
										if dir != [2]int{-1, -1} && ni >= 0 && ni < GRID_SIZE && nj >= 0 && nj < GRID_SIZE {
											g.localVirions[ni][nj] += virionsForNeighbors3 / len(g.neighbors3[i][j])
											g.localDips[ni][nj] += dipsForNeighbors3 / len(g.neighbors3[i][j])
										}
									}
								}

							} else if par_celltocell_random == false {
								//////////////////////////////
								fmt.Println("parition particles jump celltocell and randomly is false")

								if !allowVirionJump && !allowDIPJump {
									fmt.Println("Virion and DIP jump are both disabled, using viral production logic")
									// Use the new viral production function (burst or continuous based on case 4 mode)
									g.handleViralProduction(i, j, frameNum)
									// Skip the old complex logic by jumping to the end of this condition
									goto skipOldLogic

									// Calculate the distribution of virions and DIPs to each neighbor based on the ratio √3 : 2√3 : 3
									sqrt3 := math.Sqrt(3)
									ratio1 := 1.0               // sqrt3     // Weight for neighbors1
									ratio2 := 1.0 / 2           // 2 * sqrt3 // Weight for neighbors2
									ratio3 := 1.0 / (3 / sqrt3) // 3.0 // Weight for neighbors3
									totalRatio := ratio1*float64(len(g.neighbors1[i][j])) + ratio2*float64(len(g.neighbors2[i][j])) + ratio3*float64(len(g.neighbors3[i][j]))
									// if infected by virion or infected by both:
									// Calculate the number of virions and DIPs assigned to each type of neighbor
									virionsForNeighbors1 := int(math.Floor(float64(BURST_SIZE_V) * (ratio1 * float64(len(g.neighbors1[i][j]))) / totalRatio))
									virionsForNeighbors2 := int(math.Floor(float64(BURST_SIZE_V) * (ratio2 * float64(len(g.neighbors2[i][j]))) / totalRatio))
									virionsForNeighbors3 := int(math.Floor(float64(BURST_SIZE_V) * (ratio3 * float64(len(g.neighbors3[i][j]))) / totalRatio))

									// Calculate the remaining virions and DIPs
									remainingVirions := BURST_SIZE_V - (virionsForNeighbors1 + virionsForNeighbors2 + virionsForNeighbors3)

									// // Randomly distribute the remaining virions based on the ratio
									for remainingVirions > 0 {
										randVal := rand.Float64() * totalRatio
										if randVal < ratio1 && len(g.neighbors1[i][j]) > 0 {
											virionsForNeighbors1++
										} else if randVal < (ratio1+ratio2) && len(g.neighbors2[i][j]) > 0 {
											virionsForNeighbors2++
										} else if len(g.neighbors3[i][j]) > 0 {
											virionsForNeighbors3++
										}
										remainingVirions--
									}
									// if infected by vrion only or both:

									totalVirionsAtCell := g.localVirions[i][j]
									totalDIPsAtCell := g.localDips[i][j]

									// Ensure we avoid division by zero
									adjustedBurstSizeD := 0
									if totalVirionsAtCell > 0 {
										// Adjust BURST_SIZE_D based on the DIP-to-virion ratio at this cell
										dipVirionRatio := float64(totalDIPsAtCell) / float64(totalVirionsAtCell)
										adjustedBurstSizeD = BURST_SIZE_D + int(math.Floor(float64(BURST_SIZE_D)*dipVirionRatio))
									}

									// Distribute DIPs to neighbors based on the adjusted BURST_SIZE_D
									dipsForNeighbors1 := int(math.Floor(float64(adjustedBurstSizeD) * (ratio1 * float64(len(g.neighbors1[i][j]))) / totalRatio))
									dipsForNeighbors2 := int(math.Floor(float64(adjustedBurstSizeD) * (ratio2 * float64(len(g.neighbors2[i][j]))) / totalRatio))
									dipsForNeighbors3 := int(math.Floor(float64(adjustedBurstSizeD) * (ratio3 * float64(len(g.neighbors3[i][j]))) / totalRatio))
									remainingDips := adjustedBurstSizeD - (dipsForNeighbors1 + dipsForNeighbors2 + dipsForNeighbors3)

									// Randomly distribute the remaining DIPs based on the ratio
									for remainingDips > 0 {
										randVal := rand.Float64() * totalRatio
										if randVal < ratio1 && len(g.neighbors1[i][j]) > 0 {
											dipsForNeighbors1++
										} else if randVal < (ratio1+ratio2) && len(g.neighbors2[i][j]) > 0 {
											dipsForNeighbors2++
										} else if len(g.neighbors3[i][j]) > 0 {
											dipsForNeighbors3++
										}
										remainingDips--
									}
									// Distribute virions and DIPs to neighbors1
									for _, dir := range g.neighbors1[i][j] {
										ni, nj := dir[0], dir[1]
										if dir != [2]int{-1, -1} && ni >= 0 && ni < GRID_SIZE && nj >= 0 && nj < GRID_SIZE {
											if g.state[ni][nj] == SUSCEPTIBLE {
												g.localVirions[ni][nj] += virionsForNeighbors1 / len(g.neighbors1[i][j])
												g.localDips[ni][nj] += dipsForNeighbors1 / len(g.neighbors1[i][j])
											}
										}
									}

									// Distribute virions and DIPs to neighbors2
									for _, dir := range g.neighbors2[i][j] {
										ni, nj := dir[0], dir[1]
										if dir != [2]int{-1, -1} && ni >= 0 && ni < GRID_SIZE && nj >= 0 && nj < GRID_SIZE {
											g.localVirions[ni][nj] += virionsForNeighbors2 / len(g.neighbors2[i][j])
											g.localDips[ni][nj] += dipsForNeighbors2 / len(g.neighbors2[i][j])
										}
									}

									// Distribute virions and DIPs to neighbors3
									for _, dir := range g.neighbors3[i][j] {
										ni, nj := dir[0], dir[1]
										if dir != [2]int{-1, -1} && ni >= 0 && ni < GRID_SIZE && nj >= 0 && nj < GRID_SIZE {
											g.localVirions[ni][nj] += virionsForNeighbors3 / len(g.neighbors3[i][j])
											g.localDips[ni][nj] += dipsForNeighbors3 / len(g.neighbors3[i][j])
										}
									}

								} else { // "Jump" case for either virions, DIPs, or both
									fmt.Println("Virion and DIP jump are allowed to JUMP")
									if allowVirionJump {
										totalVirionsAtCell := g.localVirions[i][j]
										totalDIPsAtCell := g.localDips[i][j]
										adjustedBurstSizeD := 0
										if totalVirionsAtCell > 0 {
											dipVirionRatio := float64(totalDIPsAtCell) / float64(totalVirionsAtCell)
											adjustedBurstSizeD = BURST_SIZE_D + int(math.Floor(float64(BURST_SIZE_D)*dipVirionRatio))
										}
										if jumpRandomly {
											for v := 0; v < BURST_SIZE_V; v++ {
												ni := rand.Intn(GRID_SIZE) // Randomly select a row
												nj := rand.Intn(GRID_SIZE) // Randomly select a column

												// Apply the virion jump
												g.localVirions[ni][nj]++
											}

											// DIP jump randomly to any location
											for d := 0; d < adjustedBurstSizeD; d++ {
												ni := rand.Intn(GRID_SIZE) // Randomly select a row
												nj := rand.Intn(GRID_SIZE) // Randomly select a column

												// Apply the DIP jump
												g.localDips[ni][nj]++
											}
										} else {

											totalVirionsAtCell := g.localVirions[i][j]
											totalDIPsAtCell := g.localDips[i][j]
											adjustedBurstSizeD := 0
											if totalVirionsAtCell > 0 {
												dipVirionRatio := float64(totalDIPsAtCell) / float64(totalVirionsAtCell)
												adjustedBurstSizeD = BURST_SIZE_D + int(math.Floor(float64(BURST_SIZE_D)*dipVirionRatio))
											}

											// Virion jump logic
											virionTargets := make([]int, BURST_SIZE_V)
											for v := 0; v < BURST_SIZE_V; v++ {
												virionTargets[v] = rand.Intn(len(g.neighborsBurstArea[i][j]))
											}

											// Apply virion jumps
											for _, targetIndex := range virionTargets {
												spot := g.neighborsBurstArea[i][j][targetIndex]
												ni, nj := spot[0], spot[1]

												// Ensure the jump target is valid
												if ni < 0 || ni >= GRID_SIZE || nj < 0 || nj >= GRID_SIZE {
													// fmt.Printf("Skipping invalid jump target (%d, %d) from (%d, %d)\n", ni, nj, i, j)
													continue
												}

												// Apply the virion jump
												g.localVirions[ni][nj]++
											}

											// DIP jump logic
											dipTargets := make([]int, adjustedBurstSizeD)
											for d := 0; d < adjustedBurstSizeD; d++ {
												dipTargets[d] = rand.Intn(len(g.neighborsBurstArea[i][j]))
											}

											// Apply DIP jumps
											for _, targetIndex := range dipTargets {
												spot := g.neighborsBurstArea[i][j][targetIndex]
												ni, nj := spot[0], spot[1]

												// Ensure the jump target is valid
												if ni < 0 || ni >= GRID_SIZE || nj < 0 || nj >= GRID_SIZE {
													//fmt.Printf("Skipping invalid jump target (%d, %d) from (%d, %d)\n", ni, nj, i, j)
													continue
												}

												// Apply the DIP jump
												g.localDips[ni][nj]++
											}
										}
									}

									if allowDIPJump {
										totalVirionsAtCell := g.localVirions[i][j]
										totalDIPsAtCell := g.localDips[i][j]
										adjustedBurstSizeD := 0
										if totalVirionsAtCell > 0 {
											dipVirionRatio := float64(totalDIPsAtCell) / float64(totalVirionsAtCell)
											adjustedBurstSizeD = BURST_SIZE_D + int(math.Floor(float64(BURST_SIZE_D)*dipVirionRatio))
										}

										if jumpRandomly {
											go func() {
												for d := 0; d < adjustedBurstSizeD; d++ {
													ni := rand.Intn(GRID_SIZE)
													nj := rand.Intn(GRID_SIZE)
													g.localDips[ni][nj]++
												}
											}()
										} else {
											dipTargets := make([]int, adjustedBurstSizeD)
											for d := 0; d < adjustedBurstSizeD; d++ {
												dipTargets[d] = rand.Intn(len(g.neighborsBurstArea[i][j]))
											}
											go func() {
												for _, targetIndex := range dipTargets {
													spot := g.neighborsBurstArea[i][j][targetIndex]
													ni, nj := spot[0], spot[1]
													if ni >= 0 && ni < GRID_SIZE && nj >= 0 && nj < GRID_SIZE {
														g.localDips[ni][nj]++
													}
												}
											}()
										}
									}

								}
							}
						}
					}
				skipOldLogic:
					// Handle continuous mode cells (production logic)
					if g.state[i][j] == INFECTED_VIRION_CONTINUOUS || g.state[i][j] == INFECTED_DIP || g.state[i][j] == INFECTED_DIP_CONTINUOUS || g.state[i][j] == INFECTED_BOTH_CONTINUOUS {
						fmt.Printf("🚀 DEBUG: Found continuous state cell at (%d,%d) with state %d at frame %d\n", i, j, g.state[i][j], frameNum)
						// Use continuous production logic
						g.handleViralProduction(i, j, frameNum)
					}

					// update infected only by DIP or only by virions cells become "infected by both"
					if g.state[i][j] == INFECTED_VIRION || g.state[i][j] == INFECTED_DIP || g.state[i][j] == INFECTED_DIP_CONTINUOUS {

						if g.stateChanged[i][j] == false {
							// Check if the cell is infected by virions or DIPs

							if g.localVirions[i][j] > 0 || g.localDips[i][j] > 0 {
								// Calculate the infection probabilities

								if R == 0 || TAU == 0 {
									perParticleInfectionChance_V = RHO
								} else {
									if VStimulateIFN == true { // R=1
										perParticleInfectionChance_V = RHO * math.Exp(-ALPHA*(globalIFNperCell/float64(R)))
									} else if VStimulateIFN == false { // usually only DIP stimulate IFN in this situlation
										perParticleInfectionChance_V = RHO * math.Exp(-ALPHA*(globalIFNperCell))
									}
								}
								var probabilityVInfection, probabilityDInfection float64

								// Virion infection probability
								probabilityVInfection = 1 - math.Pow(1-perParticleInfectionChance_V, float64(g.localVirions[i][j]))
								infectedByVirion := rand.Float64() <= probabilityVInfection

								// DIP infection probability
								probabilityDInfection = 1 - math.Pow(1-(RHO*math.Exp(-ALPHA*(globalIFNperCell))), float64(g.localDips[i][j]))
								infectedByDip := rand.Float64() <= probabilityDInfection

								// Handle co-infection of already infected cells
								if g.state[i][j] == INFECTED_VIRION {
									if infectedByDip {
										fmt.Printf("COINFECT: frame %d cell (%d,%d) VIRION->BOTH by DIP; localVirions=%d localDIPs=%d pV=%.6f pD=%.6f\n",
											frameNum, i, j, g.localVirions[i][j], g.localDips[i][j], probabilityVInfection, probabilityDInfection)
										newGrid[i][j] = INFECTED_BOTH // Virion + DIP = Both
									}
									// Otherwise keep INFECTED_VIRION state
								} else if g.state[i][j] == INFECTED_DIP || g.state[i][j] == INFECTED_DIP_CONTINUOUS {
									if infectedByVirion {
										fmt.Printf("COINFECT: frame %d cell (%d,%d) DIP->BOTH by VIRION; localVirions=%d localDIPs=%d pV=%.6f pD=%.6f\n",
											frameNum, i, j, g.localVirions[i][j], g.localDips[i][j], probabilityVInfection, probabilityDInfection)
										newGrid[i][j] = INFECTED_BOTH // DIP + Virion = Both
									}
									// Otherwise keep INFECTED_DIP state
								}
							}

						}

						if g.state[i][j] == INFECTED_VIRION || g.state[i][j] == INFECTED_BOTH {

							if g.timeSinceInfectVorBoth[i][j] > IFN_DELAY+int(math.Floor(rand.NormFloat64()*float64(STD_IFN_DELAY))) && TAU > 0 {
								adjusted_DIP_IFN_stimulate := 1.0
								// if g.intraWT[i][j] > 0 {
								// 	dvgWtRatio := float64(g.intraDVG[i][j]) / float64(g.intraWT[i][j])
								// 	if g.state[i][j] == INFECTED_BOTH {
								// 		adjusted_DIP_IFN_stimulate *= dvgWtRatio * BOTH_IFN_stimulate_ratio
								// 	}
								// }
								adjusted_DIP_IFN_stimulate = BOTH_IFN_stimulate_ratio
								var totalIncreaseAmount float64
								if VStimulateIFN == true {

									if g.state[i][j] == INFECTED_VIRION {

										totalIncreaseAmount = float64(R) * float64(TIMESTEP) * ifnBothFold
									} else if g.state[i][j] == INFECTED_BOTH {
										totalIncreaseAmount = (float64(R) + adjusted_DIP_IFN_stimulate) * float64(TIMESTEP)
									}
								} else if VStimulateIFN == false {

									if g.state[i][j] == INFECTED_VIRION {

										totalIncreaseAmount = 0.0
									} else if g.state[i][j] == INFECTED_BOTH {
										totalIncreaseAmount = (adjusted_DIP_IFN_stimulate) * float64(TIMESTEP)
									}
									fmt.Println("totalIncreaseAmount", totalIncreaseAmount)
								}

								cellCount := len(g.neighborsIFNArea[i][j])

								if cellCount > 0 {
									averageIncreaseAmount := totalIncreaseAmount / float64(cellCount)

									for _, offset := range g.neighborsIFNArea[i][j] {
										ni, nj := offset[0], offset[1]

										if ni >= 0 && ni < GRID_SIZE && nj >= 0 && nj < GRID_SIZE {
											g.IFNConcentration[ni][nj] += averageIncreaseAmount
											globalIFN += averageIncreaseAmount
										}
									}
								}
							}

						}

						if g.state[i][j] == INFECTED_DIP || g.state[i][j] == INFECTED_DIP_CONTINUOUS {
							// Set DVG recovery threshold if not already set
							if g.dipLysisThreshold[i][j] == -1 {
								g.dipLysisThreshold[i][j] = int(rand.NormFloat64()*STANDARD_DVG_RECOVERY_TIME + MEAN_DVG_RECOVERY_TIME)
							}

							g.timeSinceInfectDIP[i][j] += TIMESTEP

							// Check if DVG-infected cell should recover (return to susceptible)
							if g.dipLysisThreshold[i][j] > 0 && g.timeSinceInfectDIP[i][j] >= g.dipLysisThreshold[i][j] {
								// DVG-infected cell returns to susceptible state (no particle release)
								newGrid[i][j] = SUSCEPTIBLE
								g.timeSinceInfectDIP[i][j] = -1
								g.dipLysisThreshold[i][j] = -1
								g.timeSinceSusceptible[i][j] = 0
							} else if g.timeSinceInfectDIP[i][j] > IFN_DELAY+int(math.Floor(rand.NormFloat64()*float64(STD_IFN_DELAY))) && TAU > 0 {
								// Continue producing IFN while infected
								// adjusted_DIP_IFN_stimulate := float64(g.intraDVG[i][j]) * D_only_IFN_stimulate_ratio
								adjusted_DIP_IFN_stimulate := D_only_IFN_stimulate_ratio
								totalIncreaseAmount := adjusted_DIP_IFN_stimulate * float64(TIMESTEP)

								cellCount := len(g.neighborsIFNArea[i][j])
								if cellCount > 0 {
									averageIncreaseAmount := totalIncreaseAmount / float64(cellCount)
									for _, offset := range g.neighborsIFNArea[i][j] {
										ni, nj := offset[0], offset[1]

										if ni >= 0 && ni < GRID_SIZE && nj >= 0 && nj < GRID_SIZE {
											g.IFNConcentration[ni][nj] += averageIncreaseAmount
											globalIFN += averageIncreaseAmount
										}
									}
								}
							}
						}

					}

				}
			}
		}
		// Handle potentially regrowing dead cells
		for i := 0; i < GRID_SIZE; i++ {
			for j := 0; j < GRID_SIZE; j++ {
				if g.state[i][j] == DEAD {
					g.timeSinceDead[i][j] += TIMESTEP

					// Check if any neighboring cells are susceptible, allowing for regrowth
					canRegrow := false
					neighbors := g.neighbors1[i][j]

					// Iterate over the neighbors and check if any are SUSCEPTIBLE
					for _, neighbor := range neighbors {
						ni, nj := neighbor[0], neighbor[1]

						// Ensure the neighbor indices are valid (within grid bounds)
						if ni >= 0 && ni < GRID_SIZE && nj >= 0 && nj < GRID_SIZE {

							//if g.timeSinceSusceptible[ni][nj]+g.timeSinceAntiviral[ni][nj] > int(math.Floor(rand.NormFloat64()*REGROWTH_STD+REGROWTH_MEAN)) || g.timeSinceRegrowth[ni][nj]+g.timeSinceAntiviral[ni][nj] > int(math.Floor(rand.NormFloat64()*REGROWTH_STD+REGROWTH_MEAN)) {
							//	canRegrow = true
							//  break
							//}
							if g.state[ni][nj] == SUSCEPTIBLE || g.state[ni][nj] == ANTIVIRAL {
								canRegrow = true
								break

							}
						}
					}

					// If the conditions are met, the cell regrows
					if canRegrow && g.timeSinceDead[i][j] >= int(rand.NormFloat64()*REGROWTH_STD+REGROWTH_MEAN) {
						newGrid[i][j] = REGROWTH
						g.timeSinceRegrowth[i][j] = 0
						g.timeSinceDead[i][j] = -1

					}

				}
			}
		}

		// Apply the updated grid state
		g.state = newGrid

		// Calculate and log the total virions and DIPs for each time step
		totalVirions, totalDIPs := g.totalVirions(), g.totalDIPs()
		fmt.Printf("Time step %d: Total Virions = %d, Total DIPs = %d\n", frameNum, totalVirions, totalDIPs)

		// Additional calculations based on simulation parameters for tracking purposes
		regrowthCount := g.calculateRegrowthCount()
		susceptiblePercentage := g.calculateSusceptiblePercentage()

		regrowthedOrAntiviralPercentage := g.calculateRegrowthedOrAntiviralPercentage()
		infectedPercentage := g.calculateInfectedPercentage()
		infectedDIPOnlyPercentage := g.calculateInfectedDIPOnlyPercentage()
		infectedBothPercentage := g.calculateInfectedBothPercentage()
		antiviralPercentage := g.calculateAntiviralPercentage()
		deadCellPercentage := calculateDeadCellPercentage(g.state)
		uninfectedPercentage := g.calculateUninfectedPercentage()
		plaquePercentage := g.calculatePlaquePercentage()

		// Log additional data as necessary
		fmt.Printf("Regrowth Count: %d, Susceptible: %.2f%%\n", regrowthCount, susceptiblePercentage)
		fmt.Printf("Regrowthed or Antiviral: %.2f%%, Infected: %.2f%%, DIP Only: %.2f%%, Both Infected: %.2f%%, Antiviral: %.2f%%\n",
			regrowthedOrAntiviralPercentage, infectedPercentage, infectedDIPOnlyPercentage, infectedBothPercentage, antiviralPercentage)
		fmt.Printf("Dead: %.2f%%, Uninfected: %.2f%%, Plaque: %.2f%%\n", deadCellPercentage, uninfectedPercentage, plaquePercentage)

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	} else if ifnWave == false { // ifnWave == false

		for i := 0; i < GRID_SIZE; i++ {
			for j := 0; j < GRID_SIZE; j++ {
				g.stateChanged[i][j] = false
				g.IFNConcentration[i][j] = globalIFN / float64(GRID_SIZE*GRID_SIZE)
			}
		}
		if globalIFN < 0 {
			globalIFN = -1.0
		}
		// Step 3: Update max global IFN if needed
		if globalIFN > maxGlobalIFN {

			maxGlobalIFN = globalIFN

		}

		// Traverse the grid
		for i := 0; i < GRID_SIZE; i++ {
			for j := 0; j < GRID_SIZE; j++ {
				// Skip UNEXPOSED cells entirely (never change)
				if g.state[i][j] == UNEXPOSED {
					continue
				}
				// Only consider cells that are in the SUSCEPTIBLE or REGROWTH state

				if g.state[i][j] == SUSCEPTIBLE || g.state[i][j] == REGROWTH || g.state[i][j] == INFECTED_DIP || g.state[i][j] == INFECTED_DIP_CONTINUOUS {
					if g.IFNConcentration[i][j] > 0 && TAU > 0 {

						if g.antiviralDuration[i][j] == -1 {
							g.antiviralDuration[i][j] = int(math.Floor(rand.NormFloat64()*float64(TAU)/4 + float64(TAU)))
							g.timeSinceAntiviral[i][j] = 0
						} else if g.timeSinceAntiviral[i][j] <= int(g.antiviralDuration[i][j]) {
							g.timeSinceAntiviral[i][j] += TIMESTEP
						} else {

							g.previousStates[i][j] = g.state[i][j]
							newGrid[i][j] = ANTIVIRAL
							g.timeSinceAntiviral[i][j] = -2
							g.totalAntiviralTime += g.antiviralDuration[i][j]
							if g.state[i][j] == ANTIVIRAL && !g.antiviralFlag[i][j] {
								g.antiviralFlag[i][j] = true
								g.antiviralCellCount++
							}

						}

					}

					if g.state[i][j] == SUSCEPTIBLE || g.state[i][j] == REGROWTH {
						// Check if the cell is infected by virions or DIPs
						if g.localVirions[i][j] > 0 || g.localDips[i][j] > 0 {
							// Calculate the infection probabilities

							if R == 0 || TAU == 0 {
								perParticleInfectionChance_V = RHO

							} else {
								if VStimulateIFN == true { // R=1
									perParticleInfectionChance_V = RHO * math.Exp(-ALPHA*(globalIFNperCell/float64(R)))
								} else if VStimulateIFN == false { // usually only DIP stimulate IFN in this situlation
									perParticleInfectionChance_V = RHO * math.Exp(-ALPHA*(globalIFNperCell))
								}
							}

							var probabilityVInfection, probabilityDInfection float64

							// Virion infection probability
							probabilityVInfection = 1 - math.Pow(1-perParticleInfectionChance_V, float64(g.localVirions[i][j]))
							infectedByVirion := rand.Float64() <= probabilityVInfection

							// DIP infection probability - use same logic as virion
							perParticleInfectionChance_D := perParticleInfectionChance_V
							probabilityDInfection = 1 - math.Pow(1-perParticleInfectionChance_D, float64(g.localDips[i][j]))
							infectedByDip := rand.Float64() <= probabilityDInfection

							// Determine the infection state based on virion and DIP infection
							if infectedByVirion && infectedByDip {
								if g.continuousMode {
									newGrid[i][j] = INFECTED_BOTH_CONTINUOUS
								} else {
									newGrid[i][j] = INFECTED_BOTH
								}
								g.timeSinceSusceptible[i][j] = -1
								g.timeSinceRegrowth[i][j] = -1
								// Record intracellular virus counts for continuous mode
								if g.continuousMode {
									// Calculate actual number of infecting virions and DIPs
									infectingVirions := int(math.Round(float64(g.localVirions[i][j]) * probabilityVInfection))
									infectingDIPs := int(math.Round(float64(g.localDips[i][j]) * probabilityDInfection))
									if infectingVirions > 0 {
										g.intraWT[i][j] += infectingVirions
									}
									if infectingDIPs > 0 {
										g.intraDVG[i][j] += infectingDIPs
									}
									g.infectionTime[i][j] = frameNum
								}
							} else if infectedByVirion {
								if g.continuousMode {
									newGrid[i][j] = INFECTED_VIRION_CONTINUOUS
								} else {
									newGrid[i][j] = INFECTED_VIRION
								}
								g.timeSinceSusceptible[i][j] = -1
								g.timeSinceRegrowth[i][j] = -1
								// Record intracellular virus count for continuous mode
								if g.continuousMode {
									// Calculate actual number of infecting virions
									infectingVirions := int(math.Round(float64(g.localVirions[i][j]) * probabilityVInfection))
									if infectingVirions > 0 {
										g.intraWT[i][j] += infectingVirions
									}
									g.infectionTime[i][j] = frameNum
								}
							} else if infectedByDip {
								if g.continuousMode {
									newGrid[i][j] = INFECTED_DIP_CONTINUOUS
								} else {
									newGrid[i][j] = INFECTED_DIP
								}
								g.timeSinceSusceptible[i][j] = -1
								g.timeSinceRegrowth[i][j] = -1
								// Record intracellular DVG count for continuous mode
								if g.continuousMode {
									// Calculate actual number of infecting DIPs
									infectingDIPs := int(math.Round(float64(g.localDips[i][j]) * probabilityDInfection))
									if infectingDIPs > 0 {
										g.intraDVG[i][j] += infectingDIPs
									}
									g.infectionTime[i][j] = frameNum
								}
							}
						}

						// Mark the state as changed if the cell is infected
						if newGrid[i][j] != g.state[i][j] {
							g.stateChanged[i][j] = true
						}
					}

				}

			}
		}

		// Process infected cells, no ifn wave, globally constant ifn
		for i := 0; i < GRID_SIZE; i++ {
			for j := 0; j < GRID_SIZE; j++ {
				if par_celltocell_random == true {

					allowRandomly := make([][]bool, GRID_SIZE)
					for i := range allowRandomly {
						allowRandomly[i] = make([]bool, GRID_SIZE)
					}

					// Calculate total number of cells allowed for random jumping based on k_JumpR
					totalCells := GRID_SIZE * GRID_SIZE

					randomJumpCells := int(math.Floor(float64(totalCells) * k_JumpR))

					// Randomly select randomJumpCells cells and mark them as allowRandomly
					selectedCells := make(map[[2]int]bool)
					for len(selectedCells) < randomJumpCells {
						ni := rand.Intn(GRID_SIZE)
						nj := rand.Intn(GRID_SIZE)
						selectedCells[[2]int{ni, nj}] = true
					}
					for pos := range selectedCells {
						allowRandomly[pos[0]][pos[1]] = true
					}

				}

				if g.state[i][j] == INFECTED_VIRION || g.state[i][j] == INFECTED_DIP || g.state[i][j] == INFECTED_DIP_CONTINUOUS || g.state[i][j] == INFECTED_BOTH ||
					g.state[i][j] == INFECTED_VIRION_CONTINUOUS || g.state[i][j] == INFECTED_BOTH_CONTINUOUS {
					fmt.Printf("🔍 DEBUG ifnWave=false: Processing infected cell at (%d,%d) with state %d at frame %d\n", i, j, g.state[i][j], frameNum)

					// update infected by V or BOTH cells become dead
					if g.state[i][j] == INFECTED_VIRION || g.state[i][j] == INFECTED_BOTH {

						if g.lysisThreshold[i][j] == -1 {
							g.lysisThreshold[i][j] = int(rand.NormFloat64()*STANDARD_LYSIS_TIME + MEAN_LYSIS_TIME)
						}
						g.timeSinceInfectVorBoth[i][j] += TIMESTEP
						g.timeSinceInfectDIP[i][j] = -1

						// Check if the cell should lyse and release virions and DIPs
						if g.timeSinceInfectVorBoth[i][j] > g.lysisThreshold[i][j] {
							if g.state[i][j] == INFECTED_VIRION {
								totalDeadFromV++ // Increase INFECTED_VIRION death count
							} else if g.state[i][j] == INFECTED_BOTH {
								totalDeadFromBoth++ // Increase INFECTED_BOTH death count
							}

							// After lysis, the cell becomes DEAD and virions and DIPs are spread to neighbors
							g.previousStates[i][j] = g.state[i][j] // Save the previous state before death
							newGrid[i][j] = DEAD
							g.state[i][j] = DEAD
							g.timeSinceDead[i][j] = 0
							g.timeSinceInfectVorBoth[i][j] = -1
							g.timeSinceInfectDIP[i][j] = -1
							g.lysisThreshold[i][j] = -1

							if par_celltocell_random == true {
								// Calculate adjusted burst size for DIPs based on local ratio
								totalVirionsAtCell := g.localVirions[i][j]
								totalDIPsAtCell := g.localDips[i][j]
								adjustedBurstSizeD := 0
								if totalVirionsAtCell > 0 {
									dipVirionRatio := float64(totalDIPsAtCell) / float64(totalVirionsAtCell)
									adjustedBurstSizeD = BURST_SIZE_D + int(float64(BURST_SIZE_D)*dipVirionRatio)
								}
								//  ---------------------------------------
								// Partition mode: split particles between random jump and cell-to-cell

								randomVirions := int(math.Floor(float64(BURST_SIZE_V) * k_JumpR))
								virionsForLocalDiffusion := BURST_SIZE_V - randomVirions

								randomDIPs := int(math.Floor(float64(adjustedBurstSizeD) * k_JumpR))
								dipsForLocalDiffusion := adjustedBurstSizeD - randomDIPs

								// Handle random jumps
								for v := 0; v < randomVirions; v++ {
									ni, nj := rand.Intn(GRID_SIZE), rand.Intn(GRID_SIZE)
									g.localVirions[ni][nj]++
									g.totalRandomJumpVirions++
								}
								for d := 0; d < randomDIPs; d++ {
									ni, nj := rand.Intn(GRID_SIZE), rand.Intn(GRID_SIZE)
									g.localDips[ni][nj]++
									g.totalRandomJumpDIPs++
								}

								// Handle local diffusion
								// Handle local diffusion with localVirions & localDIPs (keep original logic unchanged)
								if virionsForLocalDiffusion > 0 || dipsForLocalDiffusion > 0 {
									// Calculate the total number of valid neighbors
									totalNeighbors := 0

									// Count valid neighbors from neighbors1
									for _, dir := range g.neighbors1[i][j] {
										if dir != [2]int{-1, -1} {
											totalNeighbors++
										}
									}
									// Count valid neighbors from neighbors2
									for _, dir := range g.neighbors2[i][j] {
										if dir != [2]int{-1, -1} {
											totalNeighbors++
										}
									}
									// Count valid neighbors from neighbors3
									for _, dir := range g.neighbors3[i][j] {
										if dir != [2]int{-1, -1} {
											totalNeighbors++
										}
									}

									if totalNeighbors == 0 {
										return
									}

									// Calculate the distribution based on the ratio √3 : 2√3 : 3
									sqrt3 := math.Sqrt(3)
									ratio1 := 1.0               // sqrt3     // Weight for neighbors1
									ratio2 := 1.0 / 2           // 2 * sqrt3 // Weight for neighbors2
									ratio3 := 1.0 / (3 / sqrt3) // 3.0       // Weight for neighbors3
									totalRatio := ratio1*float64(len(g.neighbors1[i][j])) +
										ratio2*float64(len(g.neighbors2[i][j])) +
										ratio3*float64(len(g.neighbors3[i][j]))

									// Calculate virions for each neighbor group
									virionsForNeighbors1 := int(math.Floor(float64(virionsForLocalDiffusion) * (ratio1 * float64(len(g.neighbors1[i][j]))) / totalRatio))
									virionsForNeighbors2 := int(math.Floor(float64(virionsForLocalDiffusion) * (ratio2 * float64(len(g.neighbors2[i][j]))) / totalRatio))
									virionsForNeighbors3 := int(math.Floor(float64(virionsForLocalDiffusion) * (ratio3 * float64(len(g.neighbors3[i][j]))) / totalRatio))

									// Calculate remaining virions
									remainingVirions := virionsForLocalDiffusion - (virionsForNeighbors1 + virionsForNeighbors2 + virionsForNeighbors3)

									// Distribute remaining virions based on ratio
									for remainingVirions > 0 {
										randVal := rand.Float64() * totalRatio
										if randVal < ratio1 && len(g.neighbors1[i][j]) > 0 {
											virionsForNeighbors1++
										} else if randVal < (ratio1+ratio2) && len(g.neighbors2[i][j]) > 0 {
											virionsForNeighbors2++
										} else if len(g.neighbors3[i][j]) > 0 {
											virionsForNeighbors3++
										}
										remainingVirions--
									}

									// Calculate DIPs for each neighbor group (same logic)
									dipsForNeighbors1 := int(math.Floor(float64(dipsForLocalDiffusion) * (ratio1 * float64(len(g.neighbors1[i][j]))) / totalRatio))
									dipsForNeighbors2 := int(math.Floor(float64(dipsForLocalDiffusion) * (ratio2 * float64(len(g.neighbors2[i][j]))) / totalRatio))
									dipsForNeighbors3 := int(math.Floor(float64(dipsForLocalDiffusion) * (ratio3 * float64(len(g.neighbors3[i][j]))) / totalRatio))

									remainingDIPs := dipsForLocalDiffusion - (dipsForNeighbors1 + dipsForNeighbors2 + dipsForNeighbors3)

									for remainingDIPs > 0 {
										randVal := rand.Float64() * totalRatio
										if randVal < ratio1 && len(g.neighbors1[i][j]) > 0 {
											dipsForNeighbors1++
										} else if randVal < (ratio1+ratio2) && len(g.neighbors2[i][j]) > 0 {
											dipsForNeighbors2++
										} else if len(g.neighbors3[i][j]) > 0 {
											dipsForNeighbors3++
										}
										remainingDIPs--
									}
									// Distribute virions to neighbors1
									for _, dir := range g.neighbors1[i][j] {
										ni, nj := dir[0], dir[1]
										if dir != [2]int{-1, -1} && ni >= 0 && ni < GRID_SIZE && nj >= 0 && nj < GRID_SIZE {
											if g.state[ni][nj] == SUSCEPTIBLE {
												g.localVirions[ni][nj] += virionsForNeighbors1 / len(g.neighbors1[i][j])
												g.localDips[ni][nj] += dipsForNeighbors1 / len(g.neighbors1[i][j])
											}
										}
									}

									// Distribute virions to neighbors2
									for _, dir := range g.neighbors2[i][j] {
										ni, nj := dir[0], dir[1]
										if dir != [2]int{-1, -1} && ni >= 0 && ni < GRID_SIZE && nj >= 0 && nj < GRID_SIZE {
											g.localVirions[ni][nj] += virionsForNeighbors2 / len(g.neighbors2[i][j])
											g.localDips[ni][nj] += dipsForNeighbors2 / len(g.neighbors2[i][j])
										}
									}

									// Distribute virions to neighbors3
									for _, dir := range g.neighbors3[i][j] {
										ni, nj := dir[0], dir[1]
										if dir != [2]int{-1, -1} && ni >= 0 && ni < GRID_SIZE && nj >= 0 && nj < GRID_SIZE {
											g.localVirions[ni][nj] += virionsForNeighbors3 / len(g.neighbors3[i][j])
											g.localDips[ni][nj] += dipsForNeighbors3 / len(g.neighbors3[i][j])
										}
									}
								}

							} else if par_celltocell_random == false {
								if !allowVirionJump && !allowDIPJump {
									fmt.Println("Virion and DIP jump are both disabled, using viral production logic (2nd location)")
									// Use the new viral production function (burst or continuous based on case 4 mode)
									g.handleViralProduction(i, j, frameNum)
									// Skip the old complex logic by jumping to the end of this condition
									goto skipOldLogic2

									// Calculate the distribution of virions and DIPs to each neighbor based on the ratio √3 : 2√3 : 3
									sqrt3 := math.Sqrt(3)
									ratio1 := 1.0               // sqrt3     // Weight for neighbors1
									ratio2 := 1.0 / 2           // 2 * sqrt3 // Weight for neighbors2
									ratio3 := 1.0 / (3 / sqrt3) // 3.0 // Weight for neighbors3
									totalRatio := ratio1*float64(len(g.neighbors1[i][j])) + ratio2*float64(len(g.neighbors2[i][j])) + ratio3*float64(len(g.neighbors3[i][j]))

									// if infected by virion or infected by both:
									// Calculate the number of virions and DIPs assigned to each type of neighbor
									virionsForNeighbors1 := int(math.Floor(float64(BURST_SIZE_V) * (ratio1 * float64(len(g.neighbors1[i][j]))) / totalRatio))
									virionsForNeighbors2 := int(math.Floor(float64(BURST_SIZE_V) * (ratio2 * float64(len(g.neighbors2[i][j]))) / totalRatio))
									virionsForNeighbors3 := int(math.Floor(float64(BURST_SIZE_V) * (ratio3 * float64(len(g.neighbors3[i][j]))) / totalRatio))

									// Calculate the remaining virions and DIPs
									remainingVirions := BURST_SIZE_V - (virionsForNeighbors1 + virionsForNeighbors2 + virionsForNeighbors3)

									// // Randomly distribute the remaining virions based on the ratio
									for remainingVirions > 0 {
										randVal := rand.Float64() * totalRatio
										if randVal < ratio1 && len(g.neighbors1[i][j]) > 0 {
											virionsForNeighbors1++
										} else if randVal < (ratio1+ratio2) && len(g.neighbors2[i][j]) > 0 {
											virionsForNeighbors2++
										} else if len(g.neighbors3[i][j]) > 0 {
											virionsForNeighbors3++
										}
										remainingVirions--
									}
									// if infected by vrion only or both:

									totalVirionsAtCell := g.localVirions[i][j]
									totalDIPsAtCell := g.localDips[i][j]

									// Ensure we avoid division by zero
									adjustedBurstSizeD := 0
									if totalVirionsAtCell > 0 {
										// Adjust BURST_SIZE_D based on the DIP-to-virion ratio at this cell
										dipVirionRatio := float64(totalDIPsAtCell) / float64(totalVirionsAtCell)
										adjustedBurstSizeD = BURST_SIZE_D + int(math.Floor(float64(BURST_SIZE_D)*dipVirionRatio))
									}

									// Distribute DIPs to neighbors based on the adjusted BURST_SIZE_D
									dipsForNeighbors1 := int(math.Floor(float64(adjustedBurstSizeD) * (ratio1 * float64(len(g.neighbors1[i][j]))) / totalRatio))
									dipsForNeighbors2 := int(math.Floor(float64(adjustedBurstSizeD) * (ratio2 * float64(len(g.neighbors2[i][j]))) / totalRatio))
									dipsForNeighbors3 := int(math.Floor(float64(adjustedBurstSizeD) * (ratio3 * float64(len(g.neighbors3[i][j]))) / totalRatio))
									remainingDips := adjustedBurstSizeD - (dipsForNeighbors1 + dipsForNeighbors2 + dipsForNeighbors3)

									// Randomly distribute the remaining DIPs based on the ratio
									for remainingDips > 0 {
										randVal := rand.Float64() * totalRatio
										if randVal < ratio1 && len(g.neighbors1[i][j]) > 0 {
											dipsForNeighbors1++
										} else if randVal < (ratio1+ratio2) && len(g.neighbors2[i][j]) > 0 {
											dipsForNeighbors2++
										} else if len(g.neighbors3[i][j]) > 0 {
											dipsForNeighbors3++
										}
										remainingDips--
									}
									// Distribute virions and DIPs to neighbors1
									for _, dir := range g.neighbors1[i][j] {
										ni, nj := dir[0], dir[1]
										if dir != [2]int{-1, -1} && ni >= 0 && ni < GRID_SIZE && nj >= 0 && nj < GRID_SIZE {
											if g.state[ni][nj] == SUSCEPTIBLE {
												g.localVirions[ni][nj] += virionsForNeighbors1 / len(g.neighbors1[i][j])
												g.localDips[ni][nj] += dipsForNeighbors1 / len(g.neighbors1[i][j])
											}
										}
									}

									// Distribute virions and DIPs to neighbors2
									for _, dir := range g.neighbors2[i][j] {
										ni, nj := dir[0], dir[1]
										if dir != [2]int{-1, -1} && ni >= 0 && ni < GRID_SIZE && nj >= 0 && nj < GRID_SIZE {
											g.localVirions[ni][nj] += virionsForNeighbors2 / len(g.neighbors2[i][j])
											g.localDips[ni][nj] += dipsForNeighbors2 / len(g.neighbors2[i][j])
										}
									}

									// Distribute virions and DIPs to neighbors3
									for _, dir := range g.neighbors3[i][j] {
										ni, nj := dir[0], dir[1]
										if dir != [2]int{-1, -1} && ni >= 0 && ni < GRID_SIZE && nj >= 0 && nj < GRID_SIZE {
											g.localVirions[ni][nj] += virionsForNeighbors3 / len(g.neighbors3[i][j])
											g.localDips[ni][nj] += dipsForNeighbors3 / len(g.neighbors3[i][j])
										}
									}

								} else { // "Jump" case for either virions, DIPs, or both

									if allowVirionJump {
										totalVirionsAtCell := g.localVirions[i][j]
										totalDIPsAtCell := g.localDips[i][j]
										adjustedBurstSizeD := 0
										if totalVirionsAtCell > 0 {
											dipVirionRatio := float64(totalDIPsAtCell) / float64(totalVirionsAtCell)
											adjustedBurstSizeD = BURST_SIZE_D + int(math.Floor(float64(BURST_SIZE_D)*dipVirionRatio))
										}

										if jumpRandomly {
											for v := 0; v < BURST_SIZE_V; v++ {
												ni := rand.Intn(GRID_SIZE) // Randomly select a row
												nj := rand.Intn(GRID_SIZE) // Randomly select a column

												// Apply the virion jump
												g.localVirions[ni][nj]++
											}

											// DIP jump randomly to any location
											for d := 0; d < adjustedBurstSizeD; d++ {
												ni := rand.Intn(GRID_SIZE) // Randomly select a row
												nj := rand.Intn(GRID_SIZE) // Randomly select a column

												// Apply the DIP jump
												g.localDips[ni][nj]++
											}
										} else {
											// Virion jump logic
											virionTargets := make([]int, BURST_SIZE_V)
											for v := 0; v < BURST_SIZE_V; v++ {
												virionTargets[v] = rand.Intn(len(g.neighborsBurstArea[i][j]))
											}

											// Apply virion jumps
											for _, targetIndex := range virionTargets {
												spot := g.neighborsBurstArea[i][j][targetIndex]
												ni, nj := spot[0], spot[1]

												// Ensure the jump target is valid
												if ni < 0 || ni >= GRID_SIZE || nj < 0 || nj >= GRID_SIZE {
													// fmt.Printf("Skipping invalid jump target (%d, %d) from (%d, %d)\n", ni, nj, i, j)
													continue
												}

												// Apply the virion jump
												g.localVirions[ni][nj]++
											}

											// DIP jump logic
											dipTargets := make([]int, adjustedBurstSizeD)
											for d := 0; d < adjustedBurstSizeD; d++ {
												dipTargets[d] = rand.Intn(len(g.neighborsBurstArea[i][j]))
											}

											// Apply DIP jumps
											for _, targetIndex := range dipTargets {
												spot := g.neighborsBurstArea[i][j][targetIndex]
												ni, nj := spot[0], spot[1]

												// Ensure the jump target is valid
												if ni < 0 || ni >= GRID_SIZE || nj < 0 || nj >= GRID_SIZE {
													//fmt.Printf("Skipping invalid jump target (%d, %d) from (%d, %d)\n", ni, nj, i, j)
													continue
												}

												// Apply the DIP jump
												g.localDips[ni][nj]++
											}
										}
									}

									if allowDIPJump {
										totalVirionsAtCell := g.localVirions[i][j]
										totalDIPsAtCell := g.localDips[i][j]
										adjustedBurstSizeD := 0
										if totalVirionsAtCell > 0 {
											dipVirionRatio := float64(totalDIPsAtCell) / float64(totalVirionsAtCell)
											adjustedBurstSizeD = BURST_SIZE_D + int(math.Floor(float64(BURST_SIZE_D)*dipVirionRatio))
										}

										if jumpRandomly {
											go func() {
												for d := 0; d < adjustedBurstSizeD; d++ {
													ni := rand.Intn(GRID_SIZE)
													nj := rand.Intn(GRID_SIZE)
													g.localDips[ni][nj]++
												}
											}()
										} else {
											dipTargets := make([]int, adjustedBurstSizeD)
											for d := 0; d < adjustedBurstSizeD; d++ {
												dipTargets[d] = rand.Intn(len(g.neighborsBurstArea[i][j]))
											}
											go func() {
												for _, targetIndex := range dipTargets {
													spot := g.neighborsBurstArea[i][j][targetIndex]
													ni, nj := spot[0], spot[1]
													if ni >= 0 && ni < GRID_SIZE && nj >= 0 && nj < GRID_SIZE {
														g.localDips[ni][nj]++
													}
												}
											}()
										}
									}
								}
							}
						}
					}
				skipOldLogic2:
					// Handle continuous mode cells (production logic) - ifnWave = false branch
					if g.state[i][j] == INFECTED_VIRION_CONTINUOUS || g.state[i][j] == INFECTED_DIP || g.state[i][j] == INFECTED_DIP_CONTINUOUS || g.state[i][j] == INFECTED_BOTH_CONTINUOUS {
						fmt.Printf("🚀 DEBUG: Found continuous state cell at (%d,%d) with state %d at frame %d (ifnWave=false branch)\n", i, j, g.state[i][j], frameNum)
						// Use continuous production logic
						g.handleViralProduction(i, j, frameNum)
					}

					// update infected only by DIP or only by virions cells become infected by both
					if g.state[i][j] == INFECTED_VIRION || g.state[i][j] == INFECTED_DIP || g.state[i][j] == INFECTED_DIP_CONTINUOUS {

						if g.stateChanged[i][j] == false {
							// Check if the cell is infected by virions or DIPs

							if g.localVirions[i][j] > 0 || g.localDips[i][j] > 0 {
								// Calculate the infection probabilities

								if R == 0 || TAU == 0 {
									perParticleInfectionChance_V = RHO
								} else {
									if VStimulateIFN == true { // R=1
										perParticleInfectionChance_V = RHO * math.Exp(-ALPHA*(globalIFNperCell/float64(R)))
									} else if VStimulateIFN == false { // usually only DIP stimulate IFN in this situlation
										perParticleInfectionChance_V = RHO * math.Exp(-ALPHA*(globalIFNperCell))
									}
								}
								var probabilityVInfection, probabilityDInfection float64

								// Virion infection probability
								probabilityVInfection = 1 - math.Pow(1-perParticleInfectionChance_V, float64(g.localVirions[i][j]))
								infectedByVirion := rand.Float64() <= probabilityVInfection

								// DIP infection probability
								probabilityDInfection = 1 - math.Pow(1-(RHO*math.Exp(-ALPHA*(globalIFNperCell))), float64(g.localDips[i][j]))
								infectedByDip := rand.Float64() <= probabilityDInfection

								// Handle co-infection of already infected cells
								if g.state[i][j] == INFECTED_VIRION {
									if infectedByDip {
										fmt.Printf("COINFECT: frame %d cell (%d,%d) VIRION->BOTH by DIP; localVirions=%d localDIPs=%d pV=%.6f pD=%.6f\n",
											frameNum, i, j, g.localVirions[i][j], g.localDips[i][j], probabilityVInfection, probabilityDInfection)
										newGrid[i][j] = INFECTED_BOTH // Virion + DIP = Both
									}
									// Otherwise keep INFECTED_VIRION state
								} else if g.state[i][j] == INFECTED_DIP || g.state[i][j] == INFECTED_DIP_CONTINUOUS {
									if infectedByVirion {
										fmt.Printf("COINFECT: frame %d cell (%d,%d) DIP->BOTH by VIRION; localVirions=%d localDIPs=%d pV=%.6f pD=%.6f\n",
											frameNum, i, j, g.localVirions[i][j], g.localDips[i][j], probabilityVInfection, probabilityDInfection)
										newGrid[i][j] = INFECTED_BOTH // DIP + Virion = Both
									}
									// Otherwise keep INFECTED_DIP state
								}
							}

						}
						if g.state[i][j] == INFECTED_VIRION || g.state[i][j] == INFECTED_BOTH && TAU > 0 {

							if VStimulateIFN == true {
								if g.state[i][j] == INFECTED_VIRION {
									g.IFNConcentration[i][j] += float64(R) * float64(TIMESTEP) * ifnBothFold
								} else if g.state[i][j] == INFECTED_BOTH {

									// if g.intraWT[i][j] > 0 {
									// 	dvgWtRatio := float64(g.intraDVG[i][j]) / float64(g.intraWT[i][j])
									// 	if g.state[i][j] == INFECTED_BOTH {
									// 		adjusted_DIP_IFN_stimulate *= dvgWtRatio * BOTH_IFN_stimulate_ratio
									// 	}
									// }
									adjusted_DIP_IFN_stimulate = BOTH_IFN_stimulate_ratio
									g.IFNConcentration[i][j] += (float64(R) + adjusted_DIP_IFN_stimulate) * float64(TIMESTEP)
								}
							} else if VStimulateIFN == false {
								if g.state[i][j] == INFECTED_VIRION {
									// do nothing since virions do not stimulate IFN
								} else if g.state[i][j] == INFECTED_BOTH {

									// if g.intraWT[i][j] > 0 {
									// 	dvgWtRatio := float64(g.intraDVG[i][j]) / float64(g.intraWT[i][j])
									// 	if g.state[i][j] == INFECTED_BOTH {
									// 		adjusted_DIP_IFN_stimulate *= dvgWtRatio * BOTH_IFN_stimulate_ratio
									// 	}
									// }

									adjusted_DIP_IFN_stimulate = BOTH_IFN_stimulate_ratio

								}
								g.IFNConcentration[i][j] += (float64(R) + adjusted_DIP_IFN_stimulate) * float64(TIMESTEP)
							}

							globalIFN += g.IFNConcentration[i][j]

						}

						if g.state[i][j] == INFECTED_DIP || g.state[i][j] == INFECTED_DIP_CONTINUOUS {
							// Set DVG recovery threshold if not already set
							if g.dipLysisThreshold[i][j] == -1 {
								g.dipLysisThreshold[i][j] = int(rand.NormFloat64()*STANDARD_DVG_RECOVERY_TIME + MEAN_DVG_RECOVERY_TIME)
							}

							g.timeSinceInfectDIP[i][j] += TIMESTEP

							// Check if DVG-infected cell should recover (return to susceptible)
							if g.dipLysisThreshold[i][j] > 0 && g.timeSinceInfectDIP[i][j] >= g.dipLysisThreshold[i][j] {
								// DVG-infected cell returns to susceptible state (no particle release)
								newGrid[i][j] = SUSCEPTIBLE
								g.timeSinceInfectDIP[i][j] = -1
								g.dipLysisThreshold[i][j] = -1
								g.timeSinceSusceptible[i][j] = 0
							} else if g.timeSinceInfectDIP[i][j] > IFN_DELAY+int(math.Floor(rand.NormFloat64()*float64(STD_IFN_DELAY))) && TAU > 0 {
								// Continue producing IFN while infected
								//adjusted_DIP_IFN_stimulate := float64(g.intraDVG[i][j]) * D_only_IFN_stimulate_ratio
								adjusted_DIP_IFN_stimulate := D_only_IFN_stimulate_ratio
								g.IFNConcentration[i][j] += (float64(R) + adjusted_DIP_IFN_stimulate) * float64(TIMESTEP)
								globalIFN += g.IFNConcentration[i][j]
							}

						}

					}

				}
			}
		}
		// Handle potentially regrowing dead cells
		for i := 0; i < GRID_SIZE; i++ {
			for j := 0; j < GRID_SIZE; j++ {
				if g.state[i][j] == DEAD {
					g.timeSinceDead[i][j] += TIMESTEP

					// Check if any neighboring cells are susceptible, allowing for regrowth
					canRegrow := false
					neighbors := g.neighbors1[i][j]

					// Iterate over the neighbors and check if any are SUSCEPTIBLE
					for _, neighbor := range neighbors {
						ni, nj := neighbor[0], neighbor[1]

						// Ensure the neighbor indices are valid (within grid bounds)
						if ni >= 0 && ni < GRID_SIZE && nj >= 0 && nj < GRID_SIZE {

							if g.state[ni][nj] == SUSCEPTIBLE || g.state[ni][nj] == ANTIVIRAL {
								canRegrow = true
								break

							}

						}
					}

					// If the conditions are met, the cell regrows
					if canRegrow && g.timeSinceDead[i][j] >= int(rand.NormFloat64()*REGROWTH_STD+REGROWTH_MEAN) {
						newGrid[i][j] = REGROWTH
						g.timeSinceRegrowth[i][j] = 0
						g.timeSinceDead[i][j] = -1

					}

				}
			}
		}
		// IFN exponential decay

		if ifn_half_life != 0 {
			globalIFN = globalIFN * math.Pow(0.5, float64(TIMESTEP)/ifn_half_life)
			if globalIFN < (1.0 / (float64(GRID_SIZE) * float64(GRID_SIZE))) {
				globalIFN = 0
			}
		}

		globalIFNperCell = globalIFN / float64(GRID_SIZE*GRID_SIZE)
		// Apply the updated grid state
		g.state = newGrid

		// Calculate and log the total virions and DIPs for each time step
		totalVirions, totalDIPs := g.totalVirions(), g.totalDIPs()
		fmt.Printf("Time step %d: Total Virions = %d, Total DIPs = %d\n", frameNum, totalVirions, totalDIPs)

		// Additional calculations based on simulation parameters for tracking purposes
		regrowthCount := g.calculateRegrowthCount()
		susceptiblePercentage := g.calculateSusceptiblePercentage()

		regrowthedOrAntiviralPercentage := g.calculateRegrowthedOrAntiviralPercentage()
		infectedPercentage := g.calculateInfectedPercentage()
		infectedDIPOnlyPercentage := g.calculateInfectedDIPOnlyPercentage()
		infectedBothPercentage := g.calculateInfectedBothPercentage()
		antiviralPercentage := g.calculateAntiviralPercentage()
		deadCellPercentage := calculateDeadCellPercentage(g.state)
		uninfectedPercentage := g.calculateUninfectedPercentage()
		plaquePercentage := g.calculatePlaquePercentage()
		//virionDiffusionRate, dipDiffusionRate := g.calculateDiffusionRates()

		// Log additional data as necessary
		fmt.Printf("Regrowth Count: %d, Susceptible: %.2f%%", regrowthCount, susceptiblePercentage)
		fmt.Printf("Regrowthed or Antiviral: %.2f%%, Infected: %.2f%%, DIP Only: %.2f%%, Both Infected: %.2f%%, Antiviral: %.2f%%\n",
			regrowthedOrAntiviralPercentage, infectedPercentage, infectedDIPOnlyPercentage, infectedBothPercentage, antiviralPercentage)
		fmt.Printf("Dead: %.2f%%, Uninfected: %.2f%%, Plaque: %.2f%%\n", deadCellPercentage, uninfectedPercentage, plaquePercentage)
		//fmt.Printf("Virion Diffusion Rate: %d, DIP Diffusion Rate: %d\n", virionDiffusionRate, dipDiffusionRate)

	}

	// TIMESTEP = 1 hour. If 1 hour/step, use dt = 1.0

	if virion_half_life != 0 {
		for i := 0; i < GRID_SIZE; i++ {
			for j := 0; j < GRID_SIZE; j++ {
				// Update virus count using half-life formula
				factorV := math.Pow(0.5, float64(TIMESTEP)/virion_half_life)
				g.localVirions[i][j] = int(math.Floor(float64(g.localVirions[i][j])*factorV + 0.5))

				// Use per-cell DIP half-life
				hl := g.dipHalfLife[i][j]
				if hl > 0 {
					factorD := math.Pow(0.5, float64(TIMESTEP)/hl)
					g.localDips[i][j] = int(math.Floor(float64(g.localDips[i][j])*factorD + 0.5))
				}
			}
		}
	}

	// Clear any viral particles that may have accumulated on dead cell locations
	g.clearParticlesFromDeadCells()

	// Handle DIP-only infected cells clearance (become susceptible after mean=2±1 hours if still DIP-only)
	g.handleDipOnlyClearance(frameNum)

	// Test to verify dead cells have no particles (only run test every 6 hours to reduce output)
	if frameNum%6 == 0 {
		g.testDeadCellParticleClearance(frameNum)
	}

}

// Function to record simulation data into CSV at each timestep
// Function to record simulation data into CSV at each timestep
func (g *Grid) recordSimulationData(writer *csv.Writer, frameNum int) {
	totalVirions := g.totalVirions()
	totalDIPs := g.totalDIPs()
	deadCellPercentage := strconv.FormatFloat(calculateDeadCellPercentage(g.state), 'f', 6, 64)
	susceptiblePercentage := strconv.FormatFloat(g.calculateSusceptiblePercentage(), 'f', 6, 64)
	infectedPercentage := strconv.FormatFloat(g.calculateInfectedPercentage(), 'f', 6, 64)
	infectedDIPOnlyPercentage := strconv.FormatFloat(g.calculateInfectedDIPOnlyPercentage(), 'f', 6, 64)
	infectedBothPercentage := strconv.FormatFloat(g.calculateInfectedBothPercentage(), 'f', 6, 64)
	antiviralPercentage := strconv.FormatFloat(g.calculateAntiviralPercentage(), 'f', 6, 64)
	virionOnlyInfected := g.calculateVirionOnlyInfected()
	dipOnlyInfected := g.calculateDipOnlyInfected()
	bothInfected := g.calculateBothInfected()

	// Calculate DIP advantage = burstSizeD / burstSizeV
	dipAdvantage = float64(BURST_SIZE_D) / float64(BURST_SIZE_V)

	row := []string{
		strconv.Itoa(frameNum),
		strconv.FormatFloat(virion_half_life, 'f', 6, 64), // Add virion clearance rate
		strconv.FormatFloat(dip_half_life, 'f', 6, 64),    // Add DIP clearance rate
		strconv.FormatFloat(ifn_half_life, 'f', 6, 64),    // Add IFN clearance rate
		strconv.FormatFloat(globalIFN/float64(GRID_SIZE*GRID_SIZE), 'f', 6, 64),
		strconv.Itoa(totalVirions),
		strconv.Itoa(totalDIPs),
		deadCellPercentage,
		susceptiblePercentage,
		infectedPercentage,
		infectedDIPOnlyPercentage,
		infectedBothPercentage,
		antiviralPercentage,
		strconv.Itoa(g.calculateRegrowthCount()),
		strconv.FormatFloat(g.calculateSusceptiblePercentage(), 'f', 6, 64),
		strconv.FormatFloat(g.calculateRegrowthedOrAntiviralPercentage(), 'f', 6, 64),
		"variate, depending on radius 10 of IFN",
		"variate, depending on radius 10 of IFN",
		strconv.FormatFloat(RHO, 'f', 6, 64),
		strconv.Itoa(totalVirions + totalDIPs),
		strconv.FormatFloat(g.calculatePlaquePercentage(), 'f', 6, 64),
		strconv.FormatFloat(float64(maxGlobalIFN), 'f', 6, 64),
		"-1.0",
		strconv.FormatFloat(g.calculateUninfectedPercentage(), 'f', 6, 64),
		"0",
		strconv.Itoa(GRID_SIZE),
		strconv.Itoa(TIMESTEP),
		strconv.Itoa(IFN_DELAY),
		strconv.Itoa(STD_IFN_DELAY),
		strconv.FormatFloat(ALPHA, 'f', 6, 64),
		strconv.FormatFloat(RHO, 'f', 6, 64),
		strconv.FormatFloat(float64(TAU), 'f', 6, 64),
		strconv.Itoa(BURST_SIZE_V),
		strconv.FormatFloat(REGROWTH_MEAN, 'f', 6, 64),
		strconv.FormatFloat(REGROWTH_STD, 'f', 6, 64),
		strconv.Itoa(TIME_STEPS),
		strconv.FormatFloat(MEAN_LYSIS_TIME, 'f', 6, 64),
		strconv.FormatFloat(STANDARD_LYSIS_TIME, 'f', 6, 64),
		strconv.FormatFloat(float64(*flag_v_pfu_initial)/float64(GRID_SIZE*GRID_SIZE), 'f', 6, 64),
		strconv.FormatFloat(float64(*flag_d_pfu_initial)/float64(GRID_SIZE*GRID_SIZE), 'f', 6, 64),
		"-1.0",
		"-1.0",
		strconv.FormatFloat(float64(R), 'f', 6, 64),
		strconv.Itoa(BURST_SIZE_D),
		"-1.0",

		strconv.Itoa(option),
		strconv.FormatFloat(*flag_v_pfu_initial, 'f', -1, 64),
		strconv.FormatFloat(*flag_d_pfu_initial, 'f', -1, 64),
		strconv.Itoa(virionOnlyInfected),
		strconv.Itoa(dipOnlyInfected),
		strconv.Itoa(bothInfected),
		strconv.Itoa(totalDeadFromV),
		strconv.Itoa(totalDeadFromBoth),
		strconv.Itoa(virionDiffusionRate),
		strconv.Itoa(dipDiffusionRate),
		strconv.FormatFloat(k_JumpR, 'f', 6, 64),
		strconv.Itoa(jumpRadiusV),
		strconv.Itoa(jumpRadiusD),
		strconv.FormatBool(jumpRandomly),
		strconv.FormatBool(par_celltocell_random),
		strconv.FormatBool(allowVirionJump),
		strconv.FormatBool(allowDIPJump),
		strconv.Itoa(IFN_wave_radius),
		strconv.FormatBool(ifnWave),
		strconv.FormatFloat(ifnBothFold, 'f', 6, 64),
		strconv.FormatFloat(D_only_IFN_stimulate_ratio, 'f', 6, 64),
		strconv.FormatFloat(BOTH_IFN_stimulate_ratio, 'f', 6, 64),
		strconv.Itoa(g.totalRandomJumpVirions),        // New: total number of randomly jumping Virions
		strconv.Itoa(g.totalRandomJumpDIPs),           // New: total number of randomly jumping DIPs
		strconv.FormatFloat(dipAdvantage, 'f', 6, 64), // DIP advantage = burstSizeD / burstSizeV
	}

	writer.Write(row)
	writer.Flush()
}

// Convert the grid state into an image
func (g *Grid) gridToImage(videotype string) *image.RGBA {

	imgWidth := GRID_SIZE * CELL_SIZE * 2                       // Calculate the image width
	imgHeight := GRID_SIZE * CELL_SIZE * 2                      // Calculate the image height
	img := image.NewRGBA(image.Rect(0, 0, imgWidth, imgHeight)) // Create a new image
	if videotype == "states" {
		// Define colors for different states
		colors := map[int]color.Color{
			SUSCEPTIBLE:     color.RGBA{0, 0, 0, 255},       // Susceptible state: black
			INFECTED_VIRION: color.RGBA{255, 0, 0, 255},     // Infected by virion: red
			INFECTED_DIP:    color.RGBA{0, 255, 0, 255},     // Infected by DIP: green
			INFECTED_BOTH:   color.RGBA{255, 255, 0, 255},   // Infected by both: yellow
			DEAD:            color.RGBA{169, 169, 169, 255}, // Dead state: gray
			ANTIVIRAL:       color.RGBA{0, 0, 255, 255},     // Antiviral state: blue
			REGROWTH:        color.RGBA{128, 0, 128, 255},   // Regrowth state: purple
			UNEXPOSED:       color.RGBA{0, 0, 0, 255},       // UNEXPOSED: black (same as susceptible but frozen)
			// Continuous mode states (use same colors as burst mode for now)
			INFECTED_VIRION_CONTINUOUS: color.RGBA{255, 0, 0, 255},   // Infected by virion continuous: red
			INFECTED_DIP_CONTINUOUS:    color.RGBA{0, 255, 0, 255},   // Infected by DIP continuous: green
			INFECTED_BOTH_CONTINUOUS:   color.RGBA{255, 255, 0, 255}, // Infected by both continuous: yellow
		}
		fillBackground(img, color.RGBA{0, 0, 0, 255})
		for i := 0; i < GRID_SIZE; i++ {
			for j := 0; j < GRID_SIZE; j++ {
				x, y := calculateHexCenter(i, j)              // Calculate the center of each hexagon
				drawHexagon(img, x, y, colors[g.state[i][j]]) // Draw the hexagon based on the cell state
			}
		}
		// Return the image
	} else if videotype == "IFNconcentration" { // IFN concentration visualization
		black := color.RGBA{0, 0, 0, 255} // Default color (black)

		for i := 0; i < GRID_SIZE; i++ {
			for j := 0; j < GRID_SIZE; j++ {
				x, y := calculateHexCenter(i, j) // Calculate hexagon center coordinates
				ifnValue := g.IFNConcentration[i][j]

				var cellColor color.RGBA
				if ifnValue <= 0 {
					cellColor = black // IFN ≤ 0, black
				} else if ifnValue > 0 && ifnValue <= 1 {
					cellColor = color.RGBA{0, 0, 255, 255} // Blue
				} else if ifnValue > 1 && ifnValue <= 2 {
					cellColor = color.RGBA{0, 255, 0, 255} // Green
				} else if ifnValue > 2 && ifnValue <= 5 {
					cellColor = color.RGBA{255, 255, 0, 255} // Yellow
				} else if ifnValue > 5 && ifnValue <= 10 {
					cellColor = color.RGBA{255, 165, 0, 255} // Orange
				} else {
					cellColor = color.RGBA{255, 0, 0, 255} // Red
				}

				drawHexagon(img, x, y, cellColor)
			}
		}
	} else if videotype == "IFNonlyLargerThanZero" { // IFN concentration visualization
		red := color.RGBA{255, 0, 0, 255} // Cells with interferon > 0
		blue := color.RGBA{0, 0, 255, 255}
		black := color.RGBA{0, 0, 0, 255} // Default color for all other cells
		yellow := color.RGBA{255, 255, 0, 255}
		green := color.RGBA{0, 255, 0, 255}
		organge := color.RGBA{255, 165, 0, 255}
		for i := 0; i < GRID_SIZE; i++ {
			for j := 0; j < GRID_SIZE; j++ {
				x, y := calculateHexCenter(i, j) // Calculate the center of each hexagon

				// Apply color based on the specified conditions
				if g.timeSinceAntiviral[i][j] > g.antiviralDuration[i][j] {
					drawHexagon(img, x, y, blue) // blue for cells in antiviral state exceeding duration

				} else if g.timeSinceAntiviral[i][j] > 110 {
					drawHexagon(img, x, y, red) //

				} else if g.timeSinceAntiviral[i][j] > 90 {
					drawHexagon(img, x, y, organge) //

				} else if g.timeSinceAntiviral[i][j] > 70 {
					drawHexagon(img, x, y, green) //

				} else if g.timeSinceAntiviral[i][j] > 50 {
					drawHexagon(img, x, y, yellow) //

				} else {
					drawHexagon(img, x, y, black) // Black for all other cells
				}
			}
		}
	} else if videotype == "antiviralState" {

		blue := color.RGBA{0, 0, 255, 255}
		black := color.RGBA{0, 0, 0, 255} // Default color for all other cells

		for i := 0; i < GRID_SIZE; i++ {
			for j := 0; j < GRID_SIZE; j++ {
				x, y := calculateHexCenter(i, j) // Calculate the center of each hexagon

				// Apply color based on the specified conditions
				if g.timeSinceAntiviral[i][j] > g.antiviralDuration[i][j] {
					drawHexagon(img, x, y, blue) // blue for cells in antiviral state exceeding duration
				} else {
					drawHexagon(img, x, y, black) // Black for all other cells
				}
			}
		}
	} else if videotype == "particles" {

		fillBackground(img, color.RGBA{0, 0, 0, 255})
		for i := 0; i < GRID_SIZE; i++ {
			for j := 0; j < GRID_SIZE; j++ {
				x, y := calculateHexCenter(i, j)

				// Determine color based on particle presence
				hasVirion := g.localVirions[i][j] > 0
				hasDIP := g.localDips[i][j] > 0

				var particleColor color.Color
				switch {
				case hasVirion && hasDIP:
					particleColor = color.RGBA{255, 255, 0, 255} // Yellow (both present)
				case hasVirion:
					particleColor = color.RGBA{255, 0, 0, 255} // Red (Virion only)
				case hasDIP:
					particleColor = color.RGBA{0, 255, 0, 255} // Green (DIP only)
				default:
					particleColor = color.RGBA{0, 0, 0, 255} // Black (no particles)
				}

				drawHexagon(img, x, y, particleColor)

				// Optional: add particle count text at hexagon center
				//if hasVirion || hasDIP {
				//	label := fmt.Sprintf("V:%d\nD:%d", g.localVirions[i][j], g.localDips[i][j])
				//	addLabelCentered(img, x, y, label, color.White)
				//}

			}
		}

	} else if videotype == "baltes" {
		// Define colors for different states (same as "states" videotype)
		colors := map[int]color.Color{
			SUSCEPTIBLE:     color.RGBA{0, 0, 0, 255},       // Susceptible state: black
			INFECTED_VIRION: color.RGBA{255, 0, 0, 255},     // Infected by virion: red
			INFECTED_DIP:    color.RGBA{0, 255, 0, 255},     // Infected by DIP: green
			INFECTED_BOTH:   color.RGBA{255, 255, 0, 255},   // Infected by both: yellow
			DEAD:            color.RGBA{169, 169, 169, 255}, // Dead state: gray
			ANTIVIRAL:       color.RGBA{0, 0, 255, 255},     // Antiviral state: blue
			REGROWTH:        color.RGBA{128, 0, 128, 255},   // Regrowth state: purple
			UNEXPOSED:       color.RGBA{0, 0, 0, 255},       // UNEXPOSED: black (same as susceptible but frozen)
			// Continuous mode states (use same colors as burst mode for now)
			INFECTED_VIRION_CONTINUOUS: color.RGBA{255, 0, 0, 255},   // Infected by virion continuous: red
			INFECTED_DIP_CONTINUOUS:    color.RGBA{0, 255, 0, 255},   // Infected by DIP continuous: green
			INFECTED_BOTH_CONTINUOUS:   color.RGBA{255, 255, 0, 255}, // Infected by both continuous: yellow
		}

		fillBackground(img, color.RGBA{0, 0, 0, 255})

		// Precompute visual-only black overlay set based on flag_unexposedSetAreaFraction
		overlayFraction := *flag_unexposedSetAreaFraction
		useOverlay := overlayFraction > 0.0
		var overlayMask [][]bool
		if useOverlay {
			total := GRID_SIZE * GRID_SIZE
			target := int(math.Round(overlayFraction * float64(total)))
			if target < 0 {
				target = 0
			}
			if target > total {
				target = total
			}
			indices := make([]int, total)
			for idx := 0; idx < total; idx++ {
				indices[idx] = idx
			}
			rand.Shuffle(total, func(i, j int) { indices[i], indices[j] = indices[j], indices[i] })
			overlayMask = make([][]bool, GRID_SIZE)
			for i := 0; i < GRID_SIZE; i++ {
				overlayMask[i] = make([]bool, GRID_SIZE)
			}
			for k := 0; k < target; k++ {
				idx := indices[k]
				i := idx / GRID_SIZE
				j := idx % GRID_SIZE
				overlayMask[i][j] = true
			}
		}

		for i := 0; i < GRID_SIZE; i++ {
			for j := 0; j < GRID_SIZE; j++ {
				x, y := calculateHexCenter(i, j) // Calculate the center of each hexagon

				// Visualization-only overlay: if selected, draw black regardless of state
				if useOverlay && overlayMask[i][j] {
					drawHexagon(img, x, y, color.RGBA{0, 0, 0, 255})
					continue
				}

				var cellColor color.Color
				if g.state[i][j] == DEAD {
					if prevColor, exists := colors[g.previousStates[i][j]]; exists {
						cellColor = prevColor
					} else {
						cellColor = color.RGBA{169, 169, 169, 255}
					}
				} else {
					if currColor, exists := colors[g.state[i][j]]; exists {
						cellColor = currColor
					} else {
						cellColor = color.RGBA{0, 0, 0, 255}
					}
				}

				drawHexagon(img, x, y, cellColor)
			}
		}

	} else {
		fmt.Println("Error: Unknown videotype provided.")
	}

	return img // Return the image
}

func drawTextWithBackground(img *image.RGBA, x, y int, label string, textColor, borderColor, bgColor color.Color) {
	face := basicfont.Face7x13
	textWidth := len(label) * 7
	textHeight := 13

	// White background box
	bgRect := image.Rect(x-4, y-4, x+textWidth+4, y+textHeight+4)
	draw.Draw(img, bgRect, &image.Uniform{bgColor}, image.Point{}, draw.Src)

	// Text starting point
	point := fixed.Point26_6{
		X: fixed.I(x),
		Y: fixed.I(y + textHeight),
	}
	d := &font.Drawer{
		Dst:  img,
		Src:  image.NewUniform(textColor),
		Face: face,
		Dot:  point,
	}
	d.DrawString(label)
}

// addLabel draws a text label onto an image at the specified position.
func addLabel(img *image.RGBA, x, y int, label string, col color.Color) {
	point := fixed.Point26_6{
		X: fixed.I(x),
		Y: fixed.I(y),
	}
	d := &font.Drawer{
		Dst:  img,
		Src:  image.NewUniform(col),
		Face: basicfont.Face7x13, // Basic font for rendering
		Dot:  point,
	}
	d.DrawString(label)
}
func addStaticLegend(img *image.RGBA, startX, startY int) {
	// Keep original colors and label definitions unchanged
	legendItems := []string{
		"By both", "By DIP", "By Virion",
		"Antiviral", "Uninfected", "Plaque", "Regrowth",
	}
	legendColors := map[string]color.Color{
		"By both":    color.RGBA{255, 200, 0, 255},
		"By DIP":     color.RGBA{0, 255, 0, 255},
		"By Virion":  color.RGBA{255, 0, 0, 255},
		"Antiviral":  color.RGBA{0, 102, 255, 255},
		"Uninfected": color.RGBA{0, 0, 0, 255},
		"Plaque":     color.RGBA{84, 110, 122, 255},
		"Regrowth":   color.RGBA{128, 0, 128, 255},
	}

	// Calculate background box size (keep original logic)
	const (
		fontWidth   = 7
		lineSpacing = 17
		padding     = 4
	)

	maxLabelLen := 0
	for _, label := range legendItems {
		if len(label) > maxLabelLen {
			maxLabelLen = len(label)
		}
	}

	bgWidth := maxLabelLen*fontWidth + 10
	bgHeight := len(legendItems)*lineSpacing + 6

	// Draw background (keep white opaque)
	bgRect := image.Rect(
		startX-padding,
		startY-padding,
		startX+bgWidth,
		startY+bgHeight,
	)
	draw.Draw(img, bgRect, &image.Uniform{color.RGBA{255, 255, 255, 255}}, image.Point{}, draw.Src)

	// Draw legend items (keep original drawing logic)
	for i, label := range legendItems {
		yPos := startY + i*lineSpacing
		drawTextWithBackground(
			img,
			startX,
			yPos,
			label,
			legendColors[label],
			legendColors[label],
			color.RGBA{255, 255, 255, 255},
		)
	}
}

func (g *Grid) gridToImageWithGraph(frameNum int, virionOnly, dipOnly, both []float64, mode string, showLegend bool) *image.RGBA {
	const graphHeight = 100
	const spacing = 0

	gridImg := g.gridToImage(videotype)
	gridHeight := gridImg.Bounds().Dy()

	imgWidth := GRID_SIZE * CELL_SIZE * 2
	imgHeight := graphHeight + gridHeight + spacing
	canvas := image.NewRGBA(image.Rect(0, 0, imgWidth, imgHeight))

	graphImg := createInfectionGraph(frameNum, virionOnly, dipOnly, both, showLegend)
	draw.Draw(canvas, image.Rect(0, 0, imgWidth, graphHeight), graphImg, image.Point{}, draw.Src)
	draw.Draw(canvas, image.Rect(0, graphHeight+spacing, imgWidth, graphHeight+gridHeight+spacing), gridImg, image.Point{}, draw.Src)

	if showLegend {
		addStaticLegend(canvas, canvas.Bounds().Dx()-183, canvas.Bounds().Dy()-183)
	}

	return canvas
}

// Calculate the center of each hexagonal cell
func calculateHexCenter(i, j int) (int, int) {
	x := i * CELL_SIZE * 3 / 2                                                          // Calculate the x-coordinate
	y := int(float64(j)*CELL_SIZE*math.Sqrt(3) + float64(i%2)*CELL_SIZE*math.Sqrt(3)/2) // Calculate the y-coordinate
	return x, y                                                                         // Return the center coordinates
}

func drawHexagon(img *image.RGBA, x, y int, c color.Color) {
	var hex [6]image.Point
	for i := 0; i < 6; i++ {
		angle := math.Pi / 3 * float64(i) // Calculate the angle for each vertex of the hexagon
		hex[i] = image.Point{
			X: x + int(float64(CELL_SIZE)*math.Cos(angle)), // Calculate x-coordinate
			Y: y + int(float64(CELL_SIZE)*math.Sin(angle)), // Calculate y-coordinate
		}
	}
	fillHexagon(img, hex, c) // Fill the hexagon with the specified color
}

func fillHexagon(img *image.RGBA, hex [6]image.Point, c color.Color) {
	minX, minY, maxX, maxY := hex[0].X, hex[0].Y, hex[0].X, hex[0].Y // Initialize boundary values
	for _, p := range hex {
		if p.X < minX {
			minX = p.X // Update minimum x-coordinate
		}
		if p.Y < minY {
			minY = p.Y // Update minimum y-coordinate
		}
		if p.X > maxX {
			maxX = p.X // Update maximum x-coordinate
		}
		if p.Y > maxY {
			maxY = p.Y // Update maximum y-coordinate
		}
	}
	for x := minX; x <= maxX; x++ { // Iterate through x-coordinates
		for y := minY; y <= maxY; y++ { // Iterate through y-coordinates
			if isPointInHexagon(image.Point{x, y}, hex) { // Check if the point is inside the hexagon
				img.Set(x, y, c) // Set the color of the point
			}
		}
	}
}

func isPointInHexagon(p image.Point, hex [6]image.Point) bool {
	for i := 0; i < 6; i++ {
		j := (i + 1) % 6
		if (hex[j].X-hex[i].X)*(p.Y-hex[i].Y)-(hex[j].Y-hex[i].Y)*(p.X-hex[i].X) < 0 {
			return false // Return false if the point is outside the hexagon
		}
	}
	return true // Return true if the point is inside the hexagon
}

func main() {
	flag.Parse()
	fmt.Printf("Parsed ifnSpreadOption: %q\n", *flag_ifnSpreadOption)
	fmt.Printf("Parsed particleSpreadOption: %q\n", *flag_particleSpreadOption)

	// Assign parsed flag values to global variables (note dereferencing)
	BURST_SIZE_V = *flag_burstSizeV
	BURST_SIZE_D = *flag_burstSizeD
	MEAN_LYSIS_TIME = *flag_meanLysisTime
	STANDARD_LYSIS_TIME = MEAN_LYSIS_TIME / 4
	MEAN_DVG_RECOVERY_TIME = *flag_dvgRecoveryTime
	STANDARD_DVG_RECOVERY_TIME = MEAN_DVG_RECOVERY_TIME / 3 // 3±1 hours, so std = 1
	k_JumpR = *flag_kJumpR
	TAU = *flag_tau
	ifnBothFold = *flag_ifnBothFold
	RHO = *flag_rho
	lambdaDip = *flag_lambdaDip
	option = *flag_option

	// Special parameter overrides for case 4
	if option == 4 {

		STANDARD_LYSIS_TIME = MEAN_LYSIS_TIME / 4 // Recalculate standard deviation
	}

	virion_half_life = *flag_virion_half_life
	dip_half_life = *flag_dip_half_life
	ifn_half_life = *flag_ifn_half_life

	particleSpreadOption = *flag_particleSpreadOption
	ifnSpreadOption = *flag_ifnSpreadOption
	dipOption = *flag_dipOption
	// Recalculate dependent parameters (note that ifnBothFold is now float64, not *float64)
	D_only_IFN_stimulate_ratio = 5.0 * ifnBothFold
	BOTH_IFN_stimulate_ratio = 10.0 * ifnBothFold
	videotype = *flag_videotype

	// Exposure mask: enforce baltes-only activation
	if videotype != "baltes" {
		*flag_unexposedAreaFraction = 0.0
	}
	fmt.Printf("Exposure mask (uniform) fraction = %.3f (baltes-only)\n", *flag_unexposedAreaFraction)

	// Parse viral particle removal experiment parameters
	enableParticleRemoval = *flag_enableParticleRemoval
	ifnThreshold = *flag_ifnThreshold
	removalTimepoint = *flag_removalTimepoint
	removeVirionAndDIP = *flag_removeVirionAndDIP

	// VIRION-only burst mode
	virionBurstMode = *flag_virionBurstMode
	if virionBurstMode != "both" && virionBurstMode != "virionOnly" {
		log.Fatalf("Unknown virionBurstMode: %s (expected 'both' or 'virionOnly')", virionBurstMode)
	}

	// Parse random seed parameter
	randomSeed = *flag_randomSeed

	fmt.Printf("flag_videotype = %q\n", *flag_videotype)
	// Optional: print debug information
	fmt.Printf("Parameters:\n  burstSizeV = %d\n  burstSizeD = %d\n  MEAN_LYSIS_TIME = %.2f\n  kJumpR = %.2f\n  TAU = %d\n  ifnBothFold = %.2f\n  RHO = %.3f\n par_celltocell_random = %v\n",
		BURST_SIZE_V, BURST_SIZE_D, MEAN_LYSIS_TIME, k_JumpR, TAU, ifnBothFold, RHO, par_celltocell_random)

	// Print viral particle removal experiment parameters
	if enableParticleRemoval {
		fmt.Printf("Viral Particle Removal Experiment Enabled:\n  ifnThreshold = %.3f\n  removalTimepoint = %d hours\n  removeVirionAndDIP = %v\n", ifnThreshold, removalTimepoint, removeVirionAndDIP)
	}

	// --- Particle Diffusion Options ---
	particleSpreadOption = *flag_particleSpreadOption
	if particleSpreadOption == "celltocell" {
		jumpRadiusV = 0
		jumpRadiusD = 0
		jumpRandomly = false
		// k_JumpR = 0.0
		allowVirionJump = false
		allowDIPJump = false
		fmt.Println("flag main celltocell")
	} else if particleSpreadOption == "jumprandomly" {
		jumpRadiusV = 0
		jumpRadiusD = 0
		jumpRandomly = true
		// par_celltocell_random = false
		allowVirionJump = true
		allowDIPJump = true
		// k_JumpR = 1.0
		fmt.Println("flag main jump randomly")
	} else if particleSpreadOption == "jumpradius" {
		jumpRadiusV = 5
		jumpRadiusD = 5
		jumpRandomly = false
		allowVirionJump = true
		allowDIPJump = true
		// k_JumpR = 0.0
	} else if particleSpreadOption == "partition" {
		jumpRadiusV = 0
		jumpRadiusD = 0
		jumpRandomly = true
		par_celltocell_random = true
		allowVirionJump = true // Need to enable jumping
		allowDIPJump = true    // Need to enable jumping
		fmt.Println("DEBUG: par_celltocell_random set to", par_celltocell_random)

		k_JumpR = *flag_kJumpR
	} else {
		log.Fatalf("Unknown particleSpreadOption: %s", particleSpreadOption)
	}
	fmt.Println("\nParticle spread option settings:")
	fmt.Printf("  particleSpreadOption: %s\n", particleSpreadOption)
	fmt.Printf("  jumpRadiusV: %d, jumpRadiusD: %d, jumpRandomly: %v, k_JumpR: %.2f\n",
		jumpRadiusV, jumpRadiusD, jumpRandomly, k_JumpR)

	// --- IFN Propagation Options ---
	ifnSpreadOption = *flag_ifnSpreadOption

	switch ifnSpreadOption {

	case "global":
		IFN_wave_radius = 0
		ifnWave = false
		fmt.Printf("hello: ifnSpreadOption set to: %s, IFN_wave_radius: %d\n", ifnSpreadOption, IFN_wave_radius)

	case "local":
		IFN_wave_radius = 10
		ifnWave = true
		fmt.Printf("ummmm: ifnSpreadOption set to: %s, IFN_wave_radius: %d\n", ifnSpreadOption, IFN_wave_radius)

	case "noIFN":
		IFN_wave_radius = 0
		// Disable IFN: set IFN-related parameters to zero
		ifnBothFold = 0.0
		// Additionally in the model, R, ALPHA, IFN_DELAY, STD_IFN_DELAY, tau, etc. can be set to zero
		ifnWave = false
		ALPHA = 0.0
		IFN_DELAY = 0
		STD_IFN_DELAY = 0
		TAU = 0
		ifn_half_life = 0.0
	default:
		log.Fatalf("Unknown ifnSpreadOption: %s", ifnSpreadOption)
		fmt.Printf("ifnSpreadOption set to: %s, IFN_wave_radius: %d\n", ifnSpreadOption, IFN_wave_radius)

	}
	fmt.Println("\nIFN spread option settings:")
	fmt.Printf("  ifnSpreadOption: %s, IFN_wave_radius: %d, ifnBothFold: %.2f\n",
		ifnSpreadOption, IFN_wave_radius, ifnBothFold)
	fmt.Printf("flag_ifnSpreadOption = %q\n", *flag_ifnSpreadOption)
	// --- DIP Options ---
	dipOption = *flag_dipOption
	if dipOption {
		BURST_SIZE_D = *flag_burstSizeD
		// Keep D_only_IFN_stimulate_ratio default value
	} else {
		BURST_SIZE_D = 0
		D_only_IFN_stimulate_ratio = 0.0
	}
	fmt.Println("\nDIP option settings:")
	fmt.Printf("  dipOption: %v, BURST_SIZE_D: %d, D_only_IFN_stimulate_ratio: %.2f, BOTH_IFN_stimulate_ratio: %.2f\n",
		dipOption, BURST_SIZE_D, D_only_IFN_stimulate_ratio, BOTH_IFN_stimulate_ratio)

	// Simulation code can be integrated here later, this example only shows parameter setup
	fmt.Println("\nSimulation initialization complete.")
	var grid Grid

	// Set burst radius from flag
	grid.burstRadius = *flag_burstRadius

	// Set Case 4 continuous production parameters
	grid.continuousMode = *flag_continuousMode
	grid.continuousProductionRateV = *flag_continuousProductionRateV
	grid.continuousProductionRateD = *flag_continuousProductionRateD
	grid.continuousIncubationPeriod = *flag_continuousIncubationPeriod
	grid.continuousLysisTime = *flag_continuousLysisTime
	grid.initOption = *flag_option

	// Set random seed - use provided seed or current time for randomness
	if randomSeed >= 0 {
		rand.Seed(randomSeed)
		fmt.Printf("Main: Using fixed random seed: %d\n", randomSeed)
	} else {
		seed := time.Now().UnixNano()
		rand.Seed(seed)
		fmt.Printf("Main: Using time-based random seed: %d\n", seed)
	}
	// Dynamically set the value of R
	if VStimulateIFN {
		R = int(1 * ifnBothFold)
	} else {
		R = 0
	}
	grid.initialize()                // Initialize the grid
	grid.initializeNeighbors()       // Initialize the neighbors
	grid.initializeInfection(option) // Initialize the infection state

	switch {
	case TIME_STEPS > 1000:
		ticksInterval = 500.0
	case TIME_STEPS == 145:
		ticksInterval = 24.0
	case TIME_STEPS > 500:
		ticksInterval = 100.0
	case TIME_STEPS > 100:
		ticksInterval = 50.0
	case TIME_STEPS%24 == 0:
		ticksInterval = 24.0
	default:
		ticksInterval = 100.0
	}

	// Switch statement with conditional cases
	// Switch statement with conditional cases
	switch {
	case IFN_wave_radius == 0 && TAU == 12 && jumpRandomly == true:
		yMax = 0.2
	case IFN_wave_radius == 0 && TAU == 12 && jumpRandomly == true:
		yMax = 1.0
	case IFN_wave_radius == 0 && TAU == 12 && jumpRandomly == true:
		yMax = 0.03
	case IFN_wave_radius == 0 && TAU == 12 && jumpRandomly == true:
		yMax = 1.5
	case IFN_wave_radius == 10 && TAU == 12 && jumpRandomly == true:
		yMax = 20.0
	case IFN_wave_radius == 10 && TAU == 12 && jumpRandomly == true:
		yMax = 0.1
	case IFN_wave_radius == 0 && jumpRadiusV == 0 && jumpRadiusD == 0 && TAU == 12:
		yMax = 0.3
	case IFN_wave_radius == 0 && jumpRadiusV == 5 && jumpRadiusD == 0 && TAU == 12:
		yMax = 1.0
	case IFN_wave_radius == 0 && jumpRadiusV == 0 && jumpRadiusD == 5 && TAU == 12:
		yMax = 0.03
	case IFN_wave_radius == 0 && jumpRadiusV == 5 && jumpRadiusD == 5 && TAU == 12:
		yMax = 0.1
	case IFN_wave_radius == 10 && jumpRadiusV == 5 && jumpRadiusD == 5 && TAU == 12:
		yMax = 0.2
	case IFN_wave_radius == 0 && jumpRadiusV == 0 && jumpRadiusD == 0 && TAU == 24:
		yMax = 0.3
	case IFN_wave_radius == 0 && jumpRadiusV == 5 && jumpRadiusD == 5 && TAU == 24:
		yMax = 1.5
	case IFN_wave_radius == 10 && jumpRadiusV == 5 && jumpRadiusD == 5 && TAU == 24:
		yMax = 0.2

	case IFN_wave_radius == 0 && jumpRadiusV == 0 && jumpRadiusD == 0:
		yMax = 0.3
	case IFN_wave_radius == 0 && jumpRadiusV == 5 && jumpRadiusD == 5:
		yMax = 1.5
	case IFN_wave_radius == 10 && jumpRadiusV == 5 && jumpRadiusD == 5:
		yMax = 35.0

	default:
		yMax = -1.0 // Default value in case no conditions are met
	}

	folderNumber := getNextFolderNumber("./")

	// Call generateFolderName function to generate folder name
	outputFolder := generateFolderName(
		folderNumber, // Current folder number
		jumpRandomly, // DIP random jumping logic
		jumpRadiusD,  // DIP jump radius
		jumpRadiusV,  // Virion jump radius
		BURST_SIZE_D, // DIP burst size
		BURST_SIZE_V, // Virion burst size
		//V_PFU_INITIAL,   // Virion initial value
		//D_PFU_INITIAL,   // DIP initial value
		IFN_wave_radius, // IFN wave radius
		TAU,             // TAU value
		TIME_STEPS,      // Time steps
	)

	// Create folder
	os.Mkdir(outputFolder, os.ModePerm)

	err := os.MkdirAll(outputFolder, os.ModePerm)
	if err != nil {
		log.Fatalf("Failed to create folder: %v", err)
	}
	saveCurrentGoFile(outputFolder)
	csvFilePath := filepath.Join(outputFolder, "simulation_output.csv")
	videoFilePath := filepath.Join(outputFolder, "video.mp4")

	// Open a CSV file to record the infected states over time
	file, err := os.Create(csvFilePath)
	if err != nil {
		log.Fatalf("Failed to create CSV file: %v", err)
	}
	defer file.Close()

	writer := csv.NewWriter(file)
	defer writer.Flush()

	// Write the CSV headers
	headers := []string{
		"Time", "virion_half_life", "dip_half_life", "ifn_half_life", "Global IFN Concentration Per Cell", "Total Extracellular Virions",
		"Total Extracellular DIPs", "Percentage Dead Cells", "Percentage Susceptible Cells",
		"Percentage Infected Cells", "Percentage Infected DIP-only Cells",
		"Percentage Infected Both Cells", "Percentage Antiviral Cells",
		"Regrowth Count",
		"Percentage Susceptible and Antiviral (Real Susceptible cells without regrowthed ones) Cells",
		"Percentage Regrowthed or Regrowthed and Antiviral Cells",
		"Probability Virion Infection", "Probability DIP Infection",
		"Per Particle Infection Chance RHO", "Total Local Particles",
		"Plaque Percentage", "max_global_IFN", "time_all_cells_uninfected",
		"Percentage Uninfected Cells", "num_plaques", "GRID_SIZE", "TIMESTEP",
		"IFN_DELAY", "STD_IFN_DELAY", "ALPHA", "RHO", "TAU", "BURST_SIZE_V",
		"REGROWTH_MEAN", "REGROWTH_STD", "TIME_STEPS", "MEAN_LYSIS_TIME",
		"STANDARD_LYSIS_TIME", "init_v_pfu_per_cell", "init_d_pfu_per_cell",
		"MEAN_ANTI_TIME_Per_Cell", "STD_ANTI_TIME", "R", "BURST_SIZE_D", "H",
		"option", "d_pfu_initial", "v_pfu_initial", "virionOnlyInfected", "dipOnlyInfected",
		"bothInfected", "totalDeadFromV", "totalDeadFromBoth", "virionDiffusionRate", "dipDiffusionRate", "k_JumpR",
		"jumpRadiusV", "jumpRadiusD", "jumpRandomly", "par_celltocell_random",
		"allowVirionJump", "allowDIPJump", "IFN_wave_radius", "ifnWave",
		"ifnBothFold", "D_only_IFN_stimulate_ratio", "BOTH_IFN_stimulate_ratio",
		"totalRandomJumpVirions", "totalRandomJumpDIPs", "dipAdvantage",
	}

	err = writer.Write(headers)
	if err != nil {
		log.Fatalf("Failed to write CSV headers: %v", err)
	}

	// Create an MJPEG video writer
	videoWriter, err := mjpeg.New(videoFilePath, int32(GRID_SIZE*CELL_SIZE*2), int32(GRID_SIZE*CELL_SIZE*2), int32(FRAME_RATE))
	if err != nil {
		log.Fatalf("Failed to create MJPEG writer: %v", err) // Handle the error if the writer fails to create
	}
	defer videoWriter.Close() // Ensure the writer is closed when the program ends

	var buf bytes.Buffer                       // Buffer for JPEG encoding
	jpegOptions := &jpeg.Options{Quality: 100} // JPEG encoding options, quality set to 75

	var frameNumbers []int            // Slice to store frame numbers
	var deadCellPercentages []float64 // Slice to store dead cell percentages
	virionOnly := make([]float64, TIME_STEPS)
	dipOnly := make([]float64, TIME_STEPS)
	both := make([]float64, TIME_STEPS)
	// Ensure the first frame has valid values
	virionOnly[0] = 0.0
	dipOnly[0] = 0.0
	both[0] = 0.0
	// Output image save directory

	var extractedImages []*image.RGBA          // Store selected frame images
	selectedTimePoints := []int{7, 13, 19, 25} // Time points for saving simulation images

	for frameNum := 0; frameNum < TIME_STEPS; frameNum++ {

		grid.update(frameNum) // Update the grid state

		// Experimental viral particle removal (if enabled)
		grid.removeViralParticlesOutsideIFNRange(frameNum)

		// Call the function to record infected state counts at the specific frames
		grid.recordSimulationData(writer, frameNum)

		// Calculate and record the percentage of dead cells, excluding regrowth cells
		deadCellsPercentage := calculateDeadCellPercentage(grid.state)
		frameNumbers = append(frameNumbers, frameNum)                          // Record the current frame number
		deadCellPercentages = append(deadCellPercentages, deadCellsPercentage) // Record the percentage of dead cells

		// Calculate infection percentages
		virionOnly[frameNum] = float64(grid.calculateVirionOnlyInfected()) / float64(GRID_SIZE*GRID_SIZE) * 100
		dipOnly[frameNum] = float64(grid.calculateDipOnlyInfected()) / float64(GRID_SIZE*GRID_SIZE) * 100
		both[frameNum] = float64(grid.calculateBothInfected()) / float64(GRID_SIZE*GRID_SIZE) * 100

		// Check if current frame is at one of the selected time points
		for _, timePoint := range selectedTimePoints {
			if frameNum == timePoint {
				fmt.Printf("DEBUG: Saving simulation frame at frameNum=%d, timePoint=%d\n", frameNum, timePoint)
				// Create simulation result image
				img := grid.gridToImage(videotype)
				extractedImages = append(extractedImages, img)

				// Save individual frame image as simulation result
				individualFrameName := fmt.Sprintf("simulation_%d_hours.png", timePoint)
				savePNGImage(img, filepath.Join(outputFolder, individualFrameName))
				fmt.Printf("Saved simulation result frame: %s\n", individualFrameName)
			}
		}

		if frameNum > 1 {
			if frameNum%24 == 0 { // Save every 10 frames

				img := grid.gridToImageWithGraph(frameNum, virionOnly[:frameNum+1], dipOnly[:frameNum+1], both[:frameNum+1], videotype, false)

				extractedImages = append(extractedImages, img)
			}
		}

		// Log `y` values before feeding them to the graph
		log.Printf("Frame %d: Virion Only: %.2f%%, DIP Only: %.2f%%, Both: %.2f%%", frameNum, virionOnly[frameNum], dipOnly[frameNum], both[frameNum])
		// Generate the graph only if there are at least two frames of data
		var img *image.RGBA
		if frameNum > 0 {
			img = grid.gridToImageWithGraph(frameNum, virionOnly[:frameNum+1], dipOnly[:frameNum+1], both[:frameNum+1], videotype, true)
		} else {
			// For the first frame, only render the grid without the graph
			img = grid.gridToImage(videotype)
		}

		// Encode the image to JPEG format
		err = jpeg.Encode(&buf, img, jpegOptions)
		if err != nil {
			log.Fatalf("Failed to encode image: %v", err)
		}

		// Add the frame to the video
		err = videoWriter.AddFrame(buf.Bytes())
		if err != nil {
			log.Fatalf("Failed to add frame: %v", err)
		}
		buf.Reset() // Reset the buffer for the next frame

		if len(extractedImages) > 0 {
			combinedImage := combineImagesHorizontally(extractedImages)

			showLegend := false // ⬅️ Change here: control whether to add legend

			if showLegend {
				legendWidth := 180
				legendHeight := 80
				legendX := combinedImage.Bounds().Dx() - legendWidth - 20
				legendY := 20

				draw.Draw(combinedImage,
					image.Rect(legendX-5, legendY-5, legendX+legendWidth+5, legendY+legendHeight+5),
					&image.Uniform{color.RGBA{255, 255, 255, 200}},
					image.Point{}, draw.Over)

				addStaticLegend(combinedImage, legendX, legendY)
			}

			savePNGImage(combinedImage, filepath.Join(outputFolder, "selected_frames_combined.png"))
		}
	}
	log.Println("Video and graph saved successfully.") // Print a success message
	fmt.Println("ifnWave is ", ifnWave)

	// Generate comparison plots including composite_4x2_comparison.png
	generateComparisonPlots(outputFolder)
}

// Function to remove viral particles outside IFN range at specified timepoint (72 hours)
func (g *Grid) removeViralParticlesOutsideIFNRange(frameNum int) {
	// Check if this is the removal timepoint (72 hours)
	if !enableParticleRemoval || frameNum != removalTimepoint {
		return // IFN范围判定逻辑只在指定时间点（默认72小时）执行
	}

	// 只有在达到移除时间点时才打印这条消息和执行IFN范围判定
	fmt.Printf("=== Executing viral particle removal at Frame %d (%d hours) ===\n", frameNum, removalTimepoint)

	removedVirions := 0
	removedDIPs := 0
	totalCellsProcessed := 0

	for i := 0; i < GRID_SIZE; i++ {
		for j := 0; j < GRID_SIZE; j++ {
			// Check if there are any viral particles in this cell
			if g.localVirions[i][j] > 0 || g.localDips[i][j] > 0 {
				// Calculate local IFN concentration
				localIFN := 0.0
				if ifnWave && len(g.neighborsIFNArea[i][j]) > 0 {
					var regionalSumIFN float64
					neighborsCount := len(g.neighborsIFNArea[i][j])
					for _, neighbor := range g.neighborsIFNArea[i][j] {
						ni, nj := neighbor[0], neighbor[1]
						regionalSumIFN += g.IFNConcentration[ni][nj]
					}
					localIFN = regionalSumIFN / float64(neighborsCount)
				}

				// If IFN concentration is below threshold (outside IFN range), remove viral particles
				if localIFN < ifnThreshold {
					// Count particles to be removed
					removedVirions += g.localVirions[i][j]
					if removeVirionAndDIP {
						removedDIPs += g.localDips[i][j]
					}

					// Remove viral particles (but keep cell state unchanged)
					g.localVirions[i][j] = 0
					if removeVirionAndDIP {
						g.localDips[i][j] = 0
					}

					totalCellsProcessed++
				}
			}
		}
	}

	fmt.Printf("=== Frame %d (%dh): Viral Particle Removal Results ===\n", frameNum, removalTimepoint)
	fmt.Printf("    Removed Virions: %d\n", removedVirions)
	fmt.Printf("    Removed DIPs: %d\n", removedDIPs)
	fmt.Printf("    Cells processed (outside IFN range): %d\n", totalCellsProcessed)
	fmt.Printf("    IFN threshold used: %.3f\n", ifnThreshold)
	fmt.Printf("=== End of viral particle removal ===\n")
}

// Function to generate comparison plots (log and linear scale)
func generateComparisonPlots(outputFolder string) {
	// Check if simulation_output.csv exists
	simulationCSVPath := filepath.Join(outputFolder, "simulation_output.csv")
	if _, err := os.Stat(simulationCSVPath); os.IsNotExist(err) {
		fmt.Printf("⚠️  Cannot create comparison plots. Missing simulation_output.csv\n")
		return
	}

	// Check if infection_counts_by_time.csv exists in current directory
	experimentalCSVPath := "infection_counts_by_time.csv"
	if _, err := os.Stat(experimentalCSVPath); os.IsNotExist(err) {
		fmt.Printf("⚠️  Cannot create comparison plots. Missing infection_counts_by_time.csv\n")
		return
	}

	// Run the Python script to generate comparison plots
	// The script reads simulation_output.csv from the output folder and infection_counts_by_time.csv from current directory
	cmd := exec.Command("python3", "create_comparison_plot.py", outputFolder)
	cmd.Dir = "." // Run from current directory so it can find infection_counts_by_time.csv

	output, err := cmd.CombinedOutput()
	if err != nil {
		fmt.Printf("❌ Error running comparison plot script: %v\n", err)
		fmt.Printf("Output: %s\n", string(output))
		return
	}

	fmt.Printf("✅ Successfully generated comparison plots:\n")
	fmt.Printf("   - comparison_plot_log.png\n")
	fmt.Printf("   - comparison_plot_linear.png\n")
	fmt.Printf("   - comparison_plot_log.pdf\n")
	fmt.Printf("   - comparison_plot_linear.pdf\n")
	fmt.Printf("   - composite_4x2_comparison.png\n")
}

// Find the nearest unmasked cell to (i,j); returns the input if already unmasked
func (g *Grid) findNearestUnmasked(i, j int) (int, int) {
	if i >= 0 && i < GRID_SIZE && j >= 0 && j < GRID_SIZE {
		if !g.unexposedMask[i][j] {
			return i, j
		}
	}
	for rad := 1; rad < GRID_SIZE; rad++ {
		ring := generateHexRing(i, j, rad)
		for _, nb := range ring {
			nx, ny := nb[0], nb[1]
			if nx >= 0 && nx < GRID_SIZE && ny >= 0 && ny < GRID_SIZE {
				if !g.unexposedMask[nx][ny] {
					return nx, ny
				}
			}
		}
	}
	return i, j
}
