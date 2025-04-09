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
	TIME_STEPS    = 300        // Number of time steps
	GRID_SIZE     = 50         // Size of the grid
	V_PFU_INITIAL = 1          // Initial number of virions (viral particles)
	FRAME_RATE    = 1          // Frame rate for the video
	OUTPUT_VIDEO  = "0211.mp4" // Output video file name
	CELL_SIZE     = 4          // The size of each hexagonal cell
	TIMESTEP      = 1          // Time step size
	D_PFU_INITIAL = 0          // Initial number of DIPs (defective interfering particles)
)

// 定义 flag 变量（注意它们都是指针类型）
var (

	// 选项参数
	// 粒子扩散选项：可选 "celltocell"、"jumprandomly""partition" 或 "jumpradius"
	flag_particleSpreadOption = flag.String("particleSpreadOption", "jumprandomly", "Particle spread option: celltocell, jumprandomly, jumpradius, or partition")
	// 如果选择 jumprandomly，则该参数表示随机跳跃比例（0~1）
	flag_percentRandJumpR = flag.Float64("percentRandJumpR", 1.0, "Percentage of random jump for particles (0~1, used when jumprandomly is selected)")
	// IFN 传播选项：可选 "global"、"local" 或 "noIFN"
	flag_ifnSpreadOption = flag.String("ifnSpreadOption", "local", "IFN spread option: global, local, or noIFN")
	// DIP 选项：如果为 true，则启用 DIP；如果为 false，则禁用 DIP
	flag_dipOption = flag.Bool("dipOption", true, "DIP option: if true then enable DIP, if false then disable DIP")

	flag_burstSizeV       = flag.Int("burstSizeV", 50, "Number of virions released when a cell lyses")
	flag_burstSizeD       = flag.Int("burstSizeD", 100, "Number of DIPs released when a cell lyses")
	flag_meanLysisTime    = flag.Float64("meanLysisTime", 12.0, "Mean lysis time")
	flag_kJumpR           = flag.Float64("kJumpR", 0.5, "Parameter for cell-to-cell jump randomness")
	flag_tau              = flag.Int("tau", 12, "TAU value (e.g., lysis time)")
	flag_ifnBothFold      = flag.Float64("ifnBothFold", 1.0, "Fold effect for IFN stimulation")
	flag_rho              = flag.Float64("rho", 0.026, "Infection rate constant")
	flag_virion_half_life = flag.Float64("virion_half_life", 3.2, "Virion clearance rate (e.g., 3.2 d^-1)")
	flag_dip_half_life    = flag.Float64("dip_half_life", 3.2, "DIP clearance rate (e.g., 3.2 d^-1)")
	flag_ifn_half_life    = flag.Float64("ifn_half_life", 4.0, "IFN clearance rate (e.g., 3.0 d^-1)")
)

// 粒子扩散相关
var (
	particleSpreadOption  string  // "celltocell", "jumprandomly", "jumpradius"
	jumpRadiusV           int     // 例如：当选择 "jumpradius" 时，设为 5
	jumpRadiusD           int     // 同上
	jumpRandomly          bool    // 是否采用随机跳跃（当选择 "jumprandomly" 时为 true）
	k_JumpR               float64 // 随机跳跃比例
	par_celltocell_random bool
)

// IFN 传播相关
var (
	ifnSpreadOption string // "global", "local", "noIFN"
	IFN_wave_radius int    // 若 ifnSpreadOption=="local"，例如设为 10；"global" 或 "noIFN" 设为 0
	ifnWave         bool   // 是否启用 IFN 波

)

// DIP 相关
var (
	dipOption    bool // true 启用 DIP，false 禁用 DIP
	BURST_SIZE_D int
	// 当 DIP 启用时，默认 DIP 相关比例保持默认；禁用时置 0
	D_only_IFN_stimulate_ratio float64 = 5.0 * ifnBothFold
	BOTH_IFN_stimulate_ratio   float64 = 10.0 * ifnBothFold
)

// Global variables
var (
	// particleSpreadOption  = "jumpradius" // 可选："celltocell"、"jumprandomly"、"jumpradius"
	// PartionParticleSpreadOption = false // 可选："true" 或 "false"
	//percentRandJumpR = 0 // percentage of random jump particles.0~1.this number is larger than zero or not NA only when par_celltocell_random = true

	//ifnSpreadOption = "global" // 可选："global"、"local" or "noIFN"
	//dipOption =

	BURST_SIZE_V = 50 // CHANGE 50 Number of virions released when a cell lyses
	// BURST_SIZE_D          = 100  // CHANGE 100 // Number of DIPs released when a cell lyses
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

	option    = 2        // Option for infection initialization
	videotype = "states" // "states" // color in "states" or "IFN" or "IFNonlyLargerThanZero" or "antiviralState"
	RHO       = 0.026    //0.026    //0.02468  // 0.09 Infection rate constant
	// radius 10 of grid has 331 cells
	R int
	// radius 10 of grid has 331 cells,原来知识infected cell 增加R个IFN ，
	ALPHA = 1.0 // Parameter for infection probability (set to 1.5)

	REGROWTH_MEAN       = 24.0                // Mean time for regrowth
	REGROWTH_STD        = 6.0                 // Standard deviation for regrowth time
	MEAN_LYSIS_TIME     = 12.0                // Mean lysis time
	STANDARD_LYSIS_TIME = MEAN_LYSIS_TIME / 4 // Standard deviation for lysis time
	maxGlobalIFN        = -1.0                // 用于追踪最大IFN值
	globalIFN           = -1.0                // 全局IFN浓度
	globalIFNperCell    = 0.0
	IFN_DELAY           = 5
	STD_IFN_DELAY       = 1

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
	neighbors1             [GRID_SIZE][GRID_SIZE][6][2]int  // Neighbors at distance 1
	neighbors2             [GRID_SIZE][GRID_SIZE][6][2]int  // Neighbors at distance 2
	neighbors3             [GRID_SIZE][GRID_SIZE][6][2]int  // Neighbors at distance 3
	neighborsRingVirion    [GRID_SIZE][GRID_SIZE][60][2]int // Neighbors at distance 10 ring
	neighborsRingDIP       [GRID_SIZE][GRID_SIZE][60][2]int // Neighbors at distance 10 ring
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
	allowJumpRandomly      [][]bool
}

// Initialize the infection state
func (g *Grid) initializeInfection(option int) {
	rand.Seed(42)
	switch option {
	case 2:

		// Assign all virions to this cell
		if V_PFU_INITIAL > 0 {
			g.localVirions[25][25] = V_PFU_INITIAL
		} else {
			fmt.Printf("V_PFU_INITIAL<0, V_PFU_INITIAL: %d", V_PFU_INITIAL)
		}

		// Assign all DIPs to this cell
		if D_PFU_INITIAL > 0 {
			g.localDips[25][25] = D_PFU_INITIAL
		} else {
			fmt.Printf("D_PFU_INITIAL<0, D_PFU_INITIAL: %d", D_PFU_INITIAL)

		}

	case 3:
		// Random distribution
		for k := 0; k < V_PFU_INITIAL; k++ {
			i := rand.Intn(GRID_SIZE)
			j := rand.Intn(GRID_SIZE)
			g.localVirions[i][j]++
		}
		for k := 0; k < D_PFU_INITIAL; k++ {
			i := rand.Intn(GRID_SIZE)
			j := rand.Intn(GRID_SIZE)
			g.localDips[i][j]++
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

		}
	}
	fmt.Println("Grid initialized")
}

// 确保整个画布初始化为统一背景色
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
	vPFUInitial int,
	dPFUInitial int,
	ifnWaveRadius int,
	TAU int,
	timeSteps int,
) string {
	// 判断 Dinit 命名部分
	dInit := fmt.Sprintf("Dinit%d", dPFUInitial)

	// 判断 D 命名部分
	dName := ""
	if jumpRandomly {
		dName = fmt.Sprintf("DIPBst%d_JRand", burstSizeD)
	} else if jumpRadiusD > 0 {
		dName = fmt.Sprintf("DIPBst%d_J%d", burstSizeD, jumpRadiusD)
	} else if jumpRadiusD == 0 {
		dName = fmt.Sprintf("DIPBst%d_noJ", burstSizeD)
	} else {
		if burstSizeD == 0 && D_PFU_INITIAL == 0 && D_only_IFN_stimulate_ratio == 0 && jumpRadiusD == 0 {
			dName = "NoDIP"
		} else {
			dName = fmt.Sprintf("DIPBst%d", burstSizeD)
		}
	}

	// 判断 Vinit 命名部分
	vInit := ""
	if vPFUInitial > 0 {
		vInit = fmt.Sprintf("Vinit%d", vPFUInitial)
	} else if jumpRandomly {
		vInit = "JRand"
	} else if jumpRadiusV > 0 {
		vInit = fmt.Sprintf("J%d", jumpRadiusV)
	} else {
		vInit = "noJ"
	}

	vName := fmt.Sprintf("VBst%d", burstSizeV)

	// 判断 IFN 命名部分
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

func extractAndCombineFrames(grid *Grid, virionOnly, dipOnly, both []float64, selectedFrames []int, outputFolder string) {
	var extractedImages []*image.RGBA

	// 提取特定帧
	for _, timeStep := range selectedFrames {
		if timeStep >= TIME_STEPS {
			continue // 跳过超出时间步长的帧
		}

		// 创建特定帧的图像
		img := grid.gridToImageWithGraph(timeStep, virionOnly[:timeStep+1], dipOnly[:timeStep+1], both[:timeStep+1], videotype)

		// 将图像添加到切片
		extractedImages = append(extractedImages, img)

		// 保存单独的帧图像（可选）
		outputPath := filepath.Join(outputFolder, fmt.Sprintf("frame_%d.png", timeStep))
		savePNGImage(img, outputPath)
	}

	// 合并提取的图像
	combinedImage := combineImagesHorizontally(extractedImages)

	// 保存合并后的图像
	combinedImagePath := filepath.Join(outputFolder, "combined_selected_frames.png")
	savePNGImage(combinedImage, combinedImagePath)

	log.Printf("Combined selected frames saved at %s", combinedImagePath)
}

// 合并图像为一行
func combineImagesHorizontally(images []*image.RGBA) *image.RGBA {
	if len(images) == 0 {
		return nil
	}

	// 计算合并图像的宽高
	totalWidth := 0
	maxHeight := 0
	for _, img := range images {
		totalWidth += img.Bounds().Dx() // 累加宽度
		if img.Bounds().Dy() > maxHeight {
			maxHeight = img.Bounds().Dy() // 计算最大高度
		}
	}

	// 创建合并后的图像
	combinedImg := image.NewRGBA(image.Rect(0, 0, totalWidth, maxHeight))
	offsetX := 0
	for _, img := range images {
		rect := img.Bounds()
		draw.Draw(combinedImg, image.Rect(offsetX, 0, offsetX+rect.Dx(), rect.Dy()), img, rect.Min, draw.Src)
		offsetX += rect.Dx()
	}

	return combinedImg
}

// 保存 PNG 图像
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

func positionsAllowJumpRandomly(gridSize int, jumpRate float64) [][]bool {
	// 初始化全false的矩阵
	positions := make([][]bool, gridSize)
	for i := range positions {
		positions[i] = make([]bool, gridSize)
	}

	// 计算需要随机跳跃的细胞数量
	totalCells := gridSize * gridSize
	randomJumpCells := int(float64(totalCells) * jumpRate)

	// 生成唯一的随机坐标
	selectedCells := make(map[[2]int]bool)
	rand.Seed(time.Now().UnixNano()) // 设定随机种子
	for len(selectedCells) < randomJumpCells {
		ni := rand.Intn(gridSize)
		nj := rand.Intn(gridSize)
		selectedCells[[2]int{ni, nj}] = true
	}

	// 将选定的细胞标记为可跳跃
	for pos := range selectedCells {
		positions[pos[0]][pos[1]] = true
	}

	return positions
}

// 修改后的函数定义
func createInfectionGraph(frameNum int, virionOnly, dipOnly, both []float64) *image.RGBA {
	graphWidth := GRID_SIZE * CELL_SIZE * 2
	graphHeight := 200 // 图形区域高度

	// 确保有足够的点来绘制图形
	if frameNum < 1 {
		log.Fatalf("Not enough data to render the graph: frameNum = %d", frameNum)
	}

	//yMax := calculateMax(virionOnly, dipOnly, both)

	// Set a default value for yMax in case all values are zero

	// Apply clamping to each dataset
	virionOnly = clampValues(virionOnly, 0.00, yMax)
	dipOnly = clampValues(dipOnly, 0.00, yMax)
	both = clampValues(both, 0.00, yMax)

	// Plot the graph
	graph := chart.Chart{
		Width:  GRID_SIZE * CELL_SIZE * 1.51, // TODO
		Height: 100,
		XAxis: chart.XAxis{
			Name: "",
			Style: chart.Style{
				FontSize: 10.0,
			},
			Range: &chart.ContinuousRange{Min: 0, Max: xMax}, // Fixed X-axis range
			ValueFormatter: func(v interface{}) string {
				return fmt.Sprintf("%d", int(v.(float64))) // Format X-axis labels as integers
			},
			Ticks: generateTicks(xMax, ticksInterval),
		},
		YAxis: chart.YAxis{
			//Name: "Infected Cells (log10,%)",
			Name: "",
			Style: chart.Style{
				FontSize: 10.0,
			},
			Range: &chart.ContinuousRange{Min: 0.0, Max: yMax}, // Fixed Y-axis range
		},
		Series: []chart.Series{
			chart.ContinuousSeries{
				Name:    "Infected by Virion Only",
				XValues: createTimeSeries(frameNum),
				YValues: virionOnly,
				Style: chart.Style{
					StrokeColor: chart.ColorRed,
					StrokeWidth: 6.0,
				},
			},
			chart.ContinuousSeries{
				Name:    "Infected by DIP Only",
				XValues: createTimeSeries(frameNum),
				YValues: dipOnly,
				Style: chart.Style{
					StrokeColor: chart.ColorGreen,
					StrokeWidth: 6.0,
				},
			},
			chart.ContinuousSeries{
				Name:    "Infected by Both",
				XValues: createTimeSeries(frameNum),
				YValues: both,
				Style: chart.Style{
					StrokeColor: drawing.Color{R: 255, G: 165, B: 0, A: 255}, // Deep orange
					StrokeWidth: 8.0,
				},
			},
		},
	}

	// 渲染图形
	buffer := bytes.NewBuffer([]byte{})
	err := graph.Render(chart.PNG, buffer)
	if err != nil {
		log.Fatalf("Failed to render graph: %v", err)
	}

	// 解码图像
	graphImg, _, err := image.Decode(buffer)
	if err != nil {
		log.Fatalf("Failed to decode graph image: %v", err)
	}

	// 转换为 RGBA
	rgbaImg := image.NewRGBA(image.Rect(0, 0, graphWidth, graphHeight))
	draw.Draw(rgbaImg, rgbaImg.Bounds(), graphImg, image.Point{}, draw.Src)

	return rgbaImg
}

// saveCurrentGoFile saves the current Go source file into the specified output folder.
// saveCurrentGoFile saves the current Go source file with its original name and a timestamp.
func saveCurrentGoFile(outputFolder string) {
	_, currentFile, _, ok := runtime.Caller(0)
	if !ok {
		log.Println("无法获取当前 Go 文件路径")
		return
	}

	// 获取当前 Go 文件的基本名称（不带路径）
	originalFileName := filepath.Base(currentFile)

	// 生成时间戳
	timestamp := time.Now().Format("20060102_150405") // 格式：YYYYMMDD_HHMMSS

	// 目标文件名：原文件名_时间戳.go
	newFileName := fmt.Sprintf("%s_%s.go", originalFileName[:len(originalFileName)-3], timestamp)
	outputFilePath := filepath.Join(outputFolder, newFileName)

	// 读取 Go 文件内容
	content, err := ioutil.ReadFile(currentFile)
	if err != nil {
		log.Printf("无法读取文件 %s: %v\n", currentFile, err)
		return
	}

	// 确保目标文件夹存在
	if err := os.MkdirAll(outputFolder, os.ModePerm); err != nil {
		log.Printf("无法创建目录 %s: %v\n", outputFolder, err)
		return
	}

	// 写入文件
	err = ioutil.WriteFile(outputFilePath, content, 0644)
	if err != nil {
		log.Printf("无法保存文件到 %s: %v\n", outputFilePath, err)
		return
	}

	log.Printf("文件已成功保存到 %s\n", outputFilePath)
}

func getNextFolderNumber(basePath string) int {
	files, err := os.ReadDir(basePath)
	if err != nil {
		log.Fatalf("Failed to read directory %s: %v", basePath, err)
	}

	maxNumber := 0
	for _, file := range files {
		if file.IsDir() {
			// 尝试从文件夹名称解析出数字
			var folderNumber int
			_, err := fmt.Sscanf(file.Name(), "%d", &folderNumber)
			if err == nil && folderNumber > maxNumber {
				maxNumber = folderNumber
			}
		}
	}
	return maxNumber + 1 // 返回下一个可用编号
}

func transformToLogScale(data []float64) []float64 {
	transformed := make([]float64, len(data))
	for i, value := range data {
		if value > 0 {
			transformed[i] = math.Log10(value)
		} else {
			transformed[i] = math.Log10(0.0001) // 处理 log(0) 的情况，用很小的值代替
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
			if g.state[i][j] == INFECTED_VIRION || g.state[i][j] == INFECTED_DIP || g.state[i][j] == INFECTED_BOTH {
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
			if g.state[i][j] == INFECTED_DIP {
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
			if g.state[i][j] == INFECTED_BOTH {
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
			if g.state[i][j] == INFECTED_VIRION {
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
			if g.state[i][j] == INFECTED_DIP {
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
			if g.state[i][j] == INFECTED_BOTH {
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

// Calculate neighbor relationships
func (g *Grid) initializeNeighbors() {

	precomputedRingV := precomputeRing(jumpRadiusV)
	precomputedRingD := precomputeRing(jumpRadiusD)

	for i := 0; i < GRID_SIZE; i++ {
		for j := 0; j < GRID_SIZE; j++ {
			// Initialize the virion neighbors based on jumpRadiusV

			///////////////////////////////////////////
			indexV := 0
			for _, offset := range precomputedRingV {
				newI, newJ := i+offset[0], j+offset[1]
				// Ensure the new indices are within bounds
				if newI >= 0 && newI < GRID_SIZE && newJ >= 0 && newJ < GRID_SIZE {
					g.neighborsRingVirion[i][j][indexV] = [2]int{newI, newJ}
					indexV++
				}

				// Stop if we have filled all available spots
				if indexV >= len(g.neighborsRingVirion[i][j]) {
					break
				}
			}
			for ; indexV < len(g.neighborsRingVirion[i][j]); indexV++ {
				g.neighborsRingVirion[i][j][indexV] = [2]int{-1, -1}
			}

			// Initialize the DIP neighbors based on jumpRadiusD
			indexD := 0
			for _, offset := range precomputedRingD {
				newI, newJ := i+offset[0], j+offset[1]
				// Ensure the new indices are within bounds
				if newI >= 0 && newI < GRID_SIZE && newJ >= 0 && newJ < GRID_SIZE {
					g.neighborsRingDIP[i][j][indexD] = [2]int{newI, newJ}
					indexD++
				}
				// Stop if we have filled all available spots
				if indexD >= len(g.neighborsRingDIP[i][j]) {
					break
				}
			}

			for ; indexD < len(g.neighborsRingDIP[i][j]); indexD++ {
				g.neighborsRingDIP[i][j][indexD] = [2]int{-1, -1}
			}

			if ifnWave == true {

				precomputedIFNArea := precomputeIFNArea(IFN_wave_radius)
				// Initialize neighbors for IFN area
				var ifnAreaNeighbors [][2]int
				for _, offset := range precomputedIFNArea {
					newI, newJ := i+offset[0], j+offset[1]
					// Ensure the new indices are within grid bounds
					if newI >= 0 && newI < GRID_SIZE && newJ >= 0 && newJ < GRID_SIZE {
						ifnAreaNeighbors = append(ifnAreaNeighbors, [2]int{newI, newJ})
					}
				}
				g.neighborsIFNArea[i][j] = ifnAreaNeighbors

			}

		}

	}

	invalidNeighbor := [2]int{-1, -1} // Invalid neighbor coordinate

	for i := 0; i < GRID_SIZE; i++ {
		for j := 0; j < GRID_SIZE; j++ {

			if i%2 == 0 && j%2 == 0 {
				// Even centerX, even centerY
				// Neighbors at distance 1
				g.neighbors1[i][j] = [6][2]int{
					{i - 1, j},     // left up
					{i + 1, j},     // right up
					{i, j - 1},     // up
					{i, j + 1},     // down
					{i - 1, j + 1}, // left down
					{i + 1, j + 1}, // right down
				}
				// Neighbors at distance 2
				g.neighbors2[i][j] = [6][2]int{
					{i, j - 2},     // up
					{i, j + 2},     // down
					{i - 2, j - 1}, // left up
					{i + 2, j - 1}, // right up
					{i - 2, j + 1}, // left down
					{i + 2, j + 1}, // right down
				}
				// Neighbors at distance 3
				g.neighbors3[i][j] = [6][2]int{
					{i - 2, j},     // left
					{i + 2, j},     // right
					{i - 1, j - 1}, // up left
					{i + 1, j - 1}, // up right
					{i - 1, j - 2}, // up left
					{i + 1, j - 2}, // up right
				}
			} else if i%2 == 1 && j%2 == 0 {
				// Odd centerX, even centerY
				g.neighbors1[i][j] = [6][2]int{
					{i - 1, j},     // left up
					{i + 1, j},     // right up
					{i, j - 1},     // up
					{i, j + 1},     // down
					{i - 1, j + 1}, // left down
					{i + 1, j + 1}, // right down
				}
				g.neighbors2[i][j] = [6][2]int{
					{i, j - 2},     // up
					{i, j + 2},     // down
					{i - 2, j - 1}, // left up
					{i + 2, j - 1}, // right up
					{i - 2, j + 1}, // left down
					{i + 2, j + 1}, // right down
				}
				g.neighbors3[i][j] = [6][2]int{
					{i - 2, j},     // left
					{i + 2, j},     // right
					{i - 1, j - 1}, // up left
					{i + 1, j - 1}, // up right
					{i - 1, j + 2}, // down left
					{i + 1, j + 2}, // down right
				}
			} else if i%2 == 0 && j%2 == 1 {
				// Even centerX, odd centerY
				g.neighbors1[i][j] = [6][2]int{
					{i - 1, j},     // left up
					{i + 1, j},     // right up
					{i, j - 1},     // up
					{i, j + 1},     // down
					{i - 1, j + 1}, // left down
					{i + 1, j + 1}, // right down
				}
				g.neighbors2[i][j] = [6][2]int{
					{i, j - 2},     // up
					{i, j + 2},     // down
					{i - 2, j - 1}, // left up
					{i + 2, j - 1}, // right up
					{i - 2, j + 1}, // left down
					{i + 2, j + 1}, // right down
				}
				g.neighbors3[i][j] = [6][2]int{
					{i - 2, j},     // left
					{i + 2, j},     // right
					{i - 1, j - 1}, // up left
					{i + 1, j - 1}, // up right
					{i - 1, j - 2}, // down left
					{i + 1, j - 2}, // down right
				}
			} else if i%2 == 1 && j%2 == 1 {
				// Odd centerX, odd centerY
				g.neighbors1[i][j] = [6][2]int{
					{i - 1, j}, {i + 1, j}, {i, j - 1}, {i, j + 1}, {i - 1, j + 1}, {i + 1, j + 1},
				}
				g.neighbors2[i][j] = [6][2]int{
					{i, j - 2}, {i, j + 2}, {i - 2, j - 1}, {i + 2, j - 1}, {i - 2, j + 1}, {i + 2, j + 1},
				}
				g.neighbors3[i][j] = [6][2]int{
					{i - 2, j}, {i + 2, j}, {i - 1, j - 1}, {i + 1, j - 1}, {i - 1, j + 2}, {i + 1, j + 2},
				}
			}

			// Remove neighbors that are out of bounds by setting them to invalid values
			for n := 0; n < 6; n++ {
				if g.neighbors1[i][j][n][0] < 0 || g.neighbors1[i][j][n][0] >= GRID_SIZE || g.neighbors1[i][j][n][1] < 0 || g.neighbors1[i][j][n][1] >= GRID_SIZE {
					g.neighbors1[i][j][n] = invalidNeighbor
				}
				if g.neighbors2[i][j][n][0] < 0 || g.neighbors2[i][j][n][0] >= GRID_SIZE || g.neighbors2[i][j][n][1] < 0 || g.neighbors2[i][j][n][1] >= GRID_SIZE {
					g.neighbors2[i][j][n] = invalidNeighbor
				}
				if g.neighbors3[i][j][n][0] < 0 || g.neighbors3[i][j][n][0] >= GRID_SIZE || g.neighbors3[i][j][n][1] < 0 || g.neighbors3[i][j][n][1] >= GRID_SIZE {
					g.neighbors3[i][j][n] = invalidNeighbor
				}
			}
		}

	}

	fmt.Println("Neighbors initialized")

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

		//////////////// print the IFN concentration on the LEFT

		//for i := 0; i < GRID_SIZE; i++ {
		//for j := 0; j < 10; j++ {
		//fmt.Printf("Cell (%d, %d) IFN: %.5f\n", i, j, g.IFNConcentration[i][j])
		//}
		//}
		//////////////// print the IFN concentration on the LEFT

		// Step 3: Update max global IFN if needed
		if globalIFN > maxGlobalIFN {

			maxGlobalIFN = globalIFN

		}
		fmt.Printf("Global IFN concentration: %.2f\n", globalIFN)

		// Traverse the grid
		for i := 0; i < GRID_SIZE; i++ {
			for j := 0; j < GRID_SIZE; j++ {
				// Only consider cells that are in the SUSCEPTIBLE or REGROWTH state

				var regional_sumIFN float64
				neighborsCount := len(g.neighborsIFNArea[i][j])

				if ifn_half_life != 0 {
					for i := 0; i < GRID_SIZE; i++ {
						for j := 0; j < GRID_SIZE; j++ {
							// 使用半衰期公式更新 IFN 数量
							factorIFN := math.Pow(0.5, float64(TIMESTEP)/ifn_half_life)
							g.IFNConcentration[i][j] *= factorIFN
							// if g.IFNConcentration[i][j] < 1e-8 {
							// g.IFNConcentration[i][j] = 0
							//}
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

				if g.state[i][j] == SUSCEPTIBLE || g.state[i][j] == REGROWTH || g.state[i][j] == INFECTED_DIP {
					if g.IFNConcentration[i][j] > 0 && TAU > 0 {

						if g.antiviralDuration[i][j] <= -1 {
							g.antiviralDuration[i][j] = int(rand.NormFloat64()*float64(TAU)/4 + float64(TAU))
							g.timeSinceAntiviral[i][j] = 0
						} else if g.timeSinceAntiviral[i][j] <= int(g.antiviralDuration[i][j]) {
							g.timeSinceAntiviral[i][j] += TIMESTEP
						} else {

							g.previousStates[i][j] = g.state[i][j]
							newGrid[i][j] = ANTIVIRAL
							fmt.Printf("Cell (%d, %d) transitioned to ANTIVIRAL. IFN: %.5f, Previous State: %d\n",
								i, j, g.IFNConcentration[i][j], g.previousStates[i][j])

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
								newGrid[i][j] = INFECTED_BOTH
								g.timeSinceSusceptible[i][j] = -1
								g.timeSinceRegrowth[i][j] = -1
							} else if infectedByVirion {
								newGrid[i][j] = INFECTED_VIRION
								g.timeSinceSusceptible[i][j] = -1
								g.timeSinceRegrowth[i][j] = -1
							} else if infectedByDip {
								newGrid[i][j] = INFECTED_DIP
								g.timeSinceSusceptible[i][j] = -1
								g.timeSinceRegrowth[i][j] = -1
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
				neighborsCount := len(g.neighborsIFNArea[i][j])

				// Sum the IFN concentration within the IFN area
				for _, neighbor := range g.neighborsIFNArea[i][j] {
					ni, nj := neighbor[0], neighbor[1]
					regional_sumIFN += g.IFNConcentration[ni][nj]
				}

				// Calculate the average IFN concentration if there are neighbors within the radius
				var regionalAverageIFN float64
				if neighborsCount > 0 && TAU > 0 {
					regionalAverageIFN = regional_sumIFN / float64(neighborsCount)
					//fmt.Printf("Cell (%d, %d): Regional IFN Sum = %.2f, Average = %.2f\n", i, j, regional_sumIFN, regionalAverageIFN)

				} else {
					regionalAverageIFN = 0 // Default to 0 if no neighbors, though this should rarely occur
				}

				if g.state[i][j] == INFECTED_VIRION || g.state[i][j] == INFECTED_DIP || g.state[i][j] == INFECTED_BOTH {

					// update infected by V or BOTH cells become dead
					if g.state[i][j] == INFECTED_VIRION || g.state[i][j] == INFECTED_BOTH {
						g.timeSinceInfectVorBoth[i][j] += TIMESTEP
						g.timeSinceInfectDIP[i][j] = -1

						// Check if the cell should lyse and release virions and DIPs
						if g.timeSinceInfectVorBoth[i][j] > int(rand.NormFloat64()*STANDARD_LYSIS_TIME+MEAN_LYSIS_TIME) {

							// After lysis, the cell becomes DEAD and virions and DIPs are spread to neighbors
							if g.state[i][j] == INFECTED_VIRION {
								totalDeadFromV++ // 增加 INFECTED_VIRION 死亡计数
							} else if g.state[i][j] == INFECTED_BOTH {
								totalDeadFromBoth++ // 增加 INFECTED_BOTH 死亡计数
							}

							newGrid[i][j] = DEAD

							g.timeSinceDead[i][j] = 0
							g.timeSinceInfectVorBoth[i][j] = -1
							g.timeSinceInfectDIP[i][j] = -1
							////////// work on percentage of jump randomly and percentage of spread cell to cell without jump
							///////////// for k_jumpR percent cells that jump reandomly
							if par_celltocell_random == true {

								g.allowJumpRandomly = positionsAllowJumpRandomly(GRID_SIZE, k_JumpR)

								if g.allowJumpRandomly[i][j] {
									if allowVirionJump {
										if jumpRandomly {
											for v := 0; v < BURST_SIZE_V; v++ {
												ni := rand.Intn(GRID_SIZE) // Randomly select a row
												nj := rand.Intn(GRID_SIZE) // Randomly select a column

												// Apply the virion jump
												g.localVirions[ni][nj]++
											}

											// Additional DIP burst logic for infected by virion or both
											totalVirionsAtCell := g.localVirions[i][j]
											totalDIPsAtCell := g.localDips[i][j]

											// Ensure we avoid division by zero
											adjustedBurstSizeD := BURST_SIZE_D
											if totalVirionsAtCell > 0 {
												dipVirionRatio := float64(totalDIPsAtCell) / float64(totalVirionsAtCell)
												adjustedBurstSizeD = BURST_SIZE_D + int(float64(BURST_SIZE_D)*dipVirionRatio)
											}

											// DIP jump randomly to any location
											for d := 0; d < adjustedBurstSizeD; d++ {
												ni := rand.Intn(GRID_SIZE) // Randomly select a row
												nj := rand.Intn(GRID_SIZE) // Randomly select a column

												// Apply the DIP jump
												g.localDips[ni][nj]++
											}
										}
									}

									if allowDIPJump {
										if jumpRandomly {
											go func() {
												for d := 0; d < BURST_SIZE_D; d++ {
													ni := rand.Intn(GRID_SIZE)
													nj := rand.Intn(GRID_SIZE)
													g.localDips[ni][nj]++
												}
											}()
										}
									}
								} else {
									// if not allowRandomly, then jump to neighbors

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

									// If there are no valid neighbors, return early
									if totalNeighbors == 0 {
										return
									}

									// Calculate the distribution of virions and DIPs to each neighbor based on the ratio √3 : 2√3 : 3
									sqrt3 := math.Sqrt(3)
									ratio1 := 1.0               // sqrt3     // Weight for neighbors1
									ratio2 := 1.0 / 2           // 2 * sqrt3 // Weight for neighbors2
									ratio3 := 1.0 / (3 / sqrt3) // 3.0 // Weight for neighbors3
									totalRatio := ratio1*float64(len(g.neighbors1[i][j])) + ratio2*float64(len(g.neighbors2[i][j])) + ratio3*float64(len(g.neighbors3[i][j]))

									// if infected by virion or infected by both:
									// Calculate the number of virions and DIPs assigned to each type of neighbor
									virionsForNeighbors1 := int(float64(BURST_SIZE_V) * (ratio1 * float64(len(g.neighbors1[i][j]))) / totalRatio)
									virionsForNeighbors2 := int(float64(BURST_SIZE_V) * (ratio2 * float64(len(g.neighbors2[i][j]))) / totalRatio)
									virionsForNeighbors3 := int(float64(BURST_SIZE_V) * (ratio3 * float64(len(g.neighbors3[i][j]))) / totalRatio)

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
									adjustedBurstSizeD := BURST_SIZE_D

									// Adjust BURST_SIZE_D based on the DIP-to-virion ratio at this cell
									dipVirionRatio := float64(totalDIPsAtCell) / float64(totalVirionsAtCell)
									adjustedBurstSizeD = BURST_SIZE_D + int(float64(BURST_SIZE_D)*dipVirionRatio)

									// Distribute DIPs to neighbors based on the adjusted BURST_SIZE_D
									dipsForNeighbors1 := int(float64(adjustedBurstSizeD) * (ratio1 * float64(len(g.neighbors1[i][j]))) / totalRatio)
									dipsForNeighbors2 := int(float64(adjustedBurstSizeD) * (ratio2 * float64(len(g.neighbors2[i][j]))) / totalRatio)
									dipsForNeighbors3 := int(float64(adjustedBurstSizeD) * (ratio3 * float64(len(g.neighbors3[i][j]))) / totalRatio)
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

								}
							} else if par_celltocell_random == false {
								//////////////////////////////
								fmt.Println("parition particles jump celltocell and randomly is false")
								if !allowVirionJump && !allowDIPJump {
									fmt.Println("Virion and DIP jump are both disabled, NO JUMP")
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

									// If there are no valid neighbors, return early
									if totalNeighbors == 0 {
										return
									}

									// Calculate the distribution of virions and DIPs to each neighbor based on the ratio √3 : 2√3 : 3
									sqrt3 := math.Sqrt(3)
									ratio1 := 1.0               // sqrt3     // Weight for neighbors1
									ratio2 := 1.0 / 2           // 2 * sqrt3 // Weight for neighbors2
									ratio3 := 1.0 / (3 / sqrt3) // 3.0 // Weight for neighbors3
									totalRatio := ratio1*float64(len(g.neighbors1[i][j])) + ratio2*float64(len(g.neighbors2[i][j])) + ratio3*float64(len(g.neighbors3[i][j]))

									// if infected by virion or infected by both:
									// Calculate the number of virions and DIPs assigned to each type of neighbor
									virionsForNeighbors1 := int(float64(BURST_SIZE_V) * (ratio1 * float64(len(g.neighbors1[i][j]))) / totalRatio)
									virionsForNeighbors2 := int(float64(BURST_SIZE_V) * (ratio2 * float64(len(g.neighbors2[i][j]))) / totalRatio)
									virionsForNeighbors3 := int(float64(BURST_SIZE_V) * (ratio3 * float64(len(g.neighbors3[i][j]))) / totalRatio)

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
									adjustedBurstSizeD := BURST_SIZE_D

									// Adjust BURST_SIZE_D based on the DIP-to-virion ratio at this cell
									dipVirionRatio := float64(totalDIPsAtCell) / float64(totalVirionsAtCell)
									adjustedBurstSizeD = BURST_SIZE_D + int(float64(BURST_SIZE_D)*dipVirionRatio)

									// Distribute DIPs to neighbors based on the adjusted BURST_SIZE_D
									dipsForNeighbors1 := int(float64(adjustedBurstSizeD) * (ratio1 * float64(len(g.neighbors1[i][j]))) / totalRatio)
									dipsForNeighbors2 := int(float64(adjustedBurstSizeD) * (ratio2 * float64(len(g.neighbors2[i][j]))) / totalRatio)
									dipsForNeighbors3 := int(float64(adjustedBurstSizeD) * (ratio3 * float64(len(g.neighbors3[i][j]))) / totalRatio)
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
										if jumpRandomly {
											for v := 0; v < BURST_SIZE_V; v++ {
												ni := rand.Intn(GRID_SIZE) // Randomly select a row
												nj := rand.Intn(GRID_SIZE) // Randomly select a column

												// Apply the virion jump
												g.localVirions[ni][nj]++
											}

											// Additional DIP burst logic for infected by virion or both
											totalVirionsAtCell := g.localVirions[i][j]
											totalDIPsAtCell := g.localDips[i][j]

											// Ensure we avoid division by zero
											adjustedBurstSizeD := BURST_SIZE_D
											if totalVirionsAtCell > 0 {
												dipVirionRatio := float64(totalDIPsAtCell) / float64(totalVirionsAtCell)
												adjustedBurstSizeD = BURST_SIZE_D + int(float64(BURST_SIZE_D)*dipVirionRatio)
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
												virionTargets[v] = rand.Intn(len(g.neighborsRingVirion[i][j]))
											}

											// Apply virion jumps
											for _, targetIndex := range virionTargets {
												spot := g.neighborsRingVirion[i][j][targetIndex]
												ni, nj := spot[0], spot[1]

												// Ensure the jump target is valid
												if ni < 0 || ni >= GRID_SIZE || nj < 0 || nj >= GRID_SIZE {
													// fmt.Printf("Skipping invalid jump target (%d, %d) from (%d, %d)\n", ni, nj, i, j)
													continue
												}

												// Apply the virion jump
												g.localVirions[ni][nj]++
											}

											// Additional DIP burst logic for infected by virion or both
											totalVirionsAtCell := g.localVirions[i][j]
											totalDIPsAtCell := g.localDips[i][j]

											// Ensure we avoid division by zero
											adjustedBurstSizeD := BURST_SIZE_D
											if totalVirionsAtCell > 0 {
												dipVirionRatio := float64(totalDIPsAtCell) / float64(totalVirionsAtCell)
												adjustedBurstSizeD = BURST_SIZE_D + int(float64(BURST_SIZE_D)*dipVirionRatio)
											}

											// DIP jump logic
											dipTargets := make([]int, adjustedBurstSizeD)
											for d := 0; d < adjustedBurstSizeD; d++ {
												dipTargets[d] = rand.Intn(len(g.neighborsRingDIP[i][j]))
											}

											// Apply DIP jumps
											for _, targetIndex := range dipTargets {
												spot := g.neighborsRingDIP[i][j][targetIndex]
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
										if jumpRandomly {
											go func() {
												for d := 0; d < BURST_SIZE_D; d++ {
													ni := rand.Intn(GRID_SIZE)
													nj := rand.Intn(GRID_SIZE)
													g.localDips[ni][nj]++
												}
											}()
										} else {
											dipTargets := make([]int, BURST_SIZE_D)
											for d := 0; d < BURST_SIZE_D; d++ {
												dipTargets[d] = rand.Intn(len(g.neighborsRingDIP[i][j]))
											}
											go func() {
												for _, targetIndex := range dipTargets {
													spot := g.neighborsRingDIP[i][j][targetIndex]
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
					// update infected only by DIP or only by virions cells become "infected by both"
					if g.state[i][j] == INFECTED_VIRION || g.state[i][j] == INFECTED_DIP {

						if g.stateChanged[i][j] == false {
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
									newGrid[i][j] = INFECTED_BOTH
								} else if infectedByVirion {
									newGrid[i][j] = INFECTED_VIRION
								} else if infectedByDip {
									newGrid[i][j] = INFECTED_DIP
								}
							}

						}

						if g.state[i][j] == INFECTED_VIRION || g.state[i][j] == INFECTED_BOTH {

							if g.timeSinceInfectVorBoth[i][j] > IFN_DELAY+int(rand.NormFloat64()*float64(STD_IFN_DELAY)) && TAU > 0 {
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

											g.IFNConcentration[ni][nj] += averageIncreaseAmount

											globalIFN += averageIncreaseAmount

										}
									}
								}
							}

						}

						if g.state[i][j] == INFECTED_DIP {
							g.timeSinceInfectDIP[i][j] += TIMESTEP

							if g.timeSinceInfectDIP[i][j] > IFN_DELAY+int(rand.NormFloat64()*float64(STD_IFN_DELAY)) && TAU > 0 {
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

					neighbors := g.neighbors1[i][j]

					// Iterate over the neighbors and check if any are SUSCEPTIBLE
					for _, neighbor := range neighbors {
						ni, nj := neighbor[0], neighbor[1]

						// Ensure the neighbor indices are valid (within grid bounds)
						if ni >= 0 && ni < GRID_SIZE && nj >= 0 && nj < GRID_SIZE {

							if g.timeSinceSusceptible[ni][nj]+g.timeSinceAntiviral[ni][nj] > int(rand.NormFloat64()*REGROWTH_STD+REGROWTH_MEAN) || g.timeSinceRegrowth[ni][nj]+g.timeSinceAntiviral[ni][nj] > int(rand.NormFloat64()*REGROWTH_STD+REGROWTH_MEAN) {
								newGrid[i][j] = REGROWTH
								g.timeSinceRegrowth[i][j] = 0
								g.timeSinceDead[i][j] = -1
							}

						}
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
		fmt.Printf("Regrowth Count: %d, Susceptible: %.2f%%, Real Susceptible: %.2f%%\n", regrowthCount, susceptiblePercentage)
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

		// Step 3: Update max global IFN if needed
		if globalIFN > maxGlobalIFN {

			maxGlobalIFN = globalIFN

		}

		// Traverse the grid
		for i := 0; i < GRID_SIZE; i++ {
			for j := 0; j < GRID_SIZE; j++ {
				// Only consider cells that are in the SUSCEPTIBLE or REGROWTH state

				if g.state[i][j] == SUSCEPTIBLE || g.state[i][j] == REGROWTH || g.state[i][j] == INFECTED_DIP {
					if g.IFNConcentration[i][j] > 0 && TAU > 0 {

						if g.antiviralDuration[i][j] == -1 {
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

							// Determine the infection state based on virion and DIP infection
							if infectedByVirion && infectedByDip {
								newGrid[i][j] = INFECTED_BOTH
								g.timeSinceSusceptible[i][j] = -1
								g.timeSinceRegrowth[i][j] = -1
							} else if infectedByVirion {
								newGrid[i][j] = INFECTED_VIRION
								g.timeSinceSusceptible[i][j] = -1
								g.timeSinceRegrowth[i][j] = -1
							} else if infectedByDip {
								newGrid[i][j] = INFECTED_DIP
								g.timeSinceSusceptible[i][j] = -1
								g.timeSinceRegrowth[i][j] = -1
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

				if g.state[i][j] == INFECTED_VIRION || g.state[i][j] == INFECTED_DIP || g.state[i][j] == INFECTED_BOTH {

					// update infected by V or BOTH cells become dead
					if g.state[i][j] == INFECTED_VIRION || g.state[i][j] == INFECTED_BOTH {
						g.timeSinceInfectVorBoth[i][j] += TIMESTEP
						g.timeSinceInfectDIP[i][j] = -1

						// Check if the cell should lyse and release virions and DIPs
						if g.timeSinceInfectVorBoth[i][j] > int(rand.NormFloat64()*STANDARD_LYSIS_TIME+MEAN_LYSIS_TIME) {
							if g.state[i][j] == INFECTED_VIRION {
								totalDeadFromV++ // 增加 INFECTED_VIRION 死亡计数
							} else if g.state[i][j] == INFECTED_BOTH {
								totalDeadFromBoth++ // 增加 INFECTED_BOTH 死亡计数
							}
							// After lysis, the cell becomes DEAD and virions and DIPs are spread to neighbors
							newGrid[i][j] = DEAD

							g.timeSinceDead[i][j] = 0
							g.timeSinceInfectVorBoth[i][j] = -1
							g.timeSinceInfectDIP[i][j] = -1

							if !allowVirionJump && !allowDIPJump {
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

								// If there are no valid neighbors, return early
								if totalNeighbors == 0 {
									return
								}

								// Calculate the distribution of virions and DIPs to each neighbor based on the ratio √3 : 2√3 : 3
								sqrt3 := math.Sqrt(3)
								ratio1 := 1.0               // sqrt3     // Weight for neighbors1
								ratio2 := 1.0 / 2           // 2 * sqrt3 // Weight for neighbors2
								ratio3 := 1.0 / (3 / sqrt3) // 3.0 // Weight for neighbors3
								totalRatio := ratio1*float64(len(g.neighbors1[i][j])) + ratio2*float64(len(g.neighbors2[i][j])) + ratio3*float64(len(g.neighbors3[i][j]))

								// if infected by virion or infected by both:
								// Calculate the number of virions and DIPs assigned to each type of neighbor
								virionsForNeighbors1 := int(float64(BURST_SIZE_V) * (ratio1 * float64(len(g.neighbors1[i][j]))) / totalRatio)
								virionsForNeighbors2 := int(float64(BURST_SIZE_V) * (ratio2 * float64(len(g.neighbors2[i][j]))) / totalRatio)
								virionsForNeighbors3 := int(float64(BURST_SIZE_V) * (ratio3 * float64(len(g.neighbors3[i][j]))) / totalRatio)

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
								adjustedBurstSizeD := BURST_SIZE_D

								// Adjust BURST_SIZE_D based on the DIP-to-virion ratio at this cell
								dipVirionRatio := float64(totalDIPsAtCell) / float64(totalVirionsAtCell)
								adjustedBurstSizeD = BURST_SIZE_D + int(float64(BURST_SIZE_D)*dipVirionRatio)

								// Distribute DIPs to neighbors based on the adjusted BURST_SIZE_D
								dipsForNeighbors1 := int(float64(adjustedBurstSizeD) * (ratio1 * float64(len(g.neighbors1[i][j]))) / totalRatio)
								dipsForNeighbors2 := int(float64(adjustedBurstSizeD) * (ratio2 * float64(len(g.neighbors2[i][j]))) / totalRatio)
								dipsForNeighbors3 := int(float64(adjustedBurstSizeD) * (ratio3 * float64(len(g.neighbors3[i][j]))) / totalRatio)
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
									if jumpRandomly {
										for v := 0; v < BURST_SIZE_V; v++ {
											ni := rand.Intn(GRID_SIZE) // Randomly select a row
											nj := rand.Intn(GRID_SIZE) // Randomly select a column

											// Apply the virion jump
											g.localVirions[ni][nj]++
										}

										// Additional DIP burst logic for infected by virion or both
										totalVirionsAtCell := g.localVirions[i][j]
										totalDIPsAtCell := g.localDips[i][j]

										// Ensure we avoid division by zero
										adjustedBurstSizeD := BURST_SIZE_D
										if totalVirionsAtCell > 0 {
											dipVirionRatio := float64(totalDIPsAtCell) / float64(totalVirionsAtCell)
											adjustedBurstSizeD = BURST_SIZE_D + int(float64(BURST_SIZE_D)*dipVirionRatio)
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
											virionTargets[v] = rand.Intn(len(g.neighborsRingVirion[i][j]))
										}

										// Apply virion jumps
										for _, targetIndex := range virionTargets {
											spot := g.neighborsRingVirion[i][j][targetIndex]
											ni, nj := spot[0], spot[1]

											// Ensure the jump target is valid
											if ni < 0 || ni >= GRID_SIZE || nj < 0 || nj >= GRID_SIZE {
												// fmt.Printf("Skipping invalid jump target (%d, %d) from (%d, %d)\n", ni, nj, i, j)
												continue
											}

											// Apply the virion jump
											g.localVirions[ni][nj]++
										}

										// Additional DIP burst logic for infected by virion or both
										totalVirionsAtCell := g.localVirions[i][j]
										totalDIPsAtCell := g.localDips[i][j]

										// Ensure we avoid division by zero
										adjustedBurstSizeD := BURST_SIZE_D
										if totalVirionsAtCell > 0 {
											dipVirionRatio := float64(totalDIPsAtCell) / float64(totalVirionsAtCell)
											adjustedBurstSizeD = BURST_SIZE_D + int(float64(BURST_SIZE_D)*dipVirionRatio)
										}

										// DIP jump logic
										dipTargets := make([]int, adjustedBurstSizeD)
										for d := 0; d < adjustedBurstSizeD; d++ {
											dipTargets[d] = rand.Intn(len(g.neighborsRingDIP[i][j]))
										}

										// Apply DIP jumps
										for _, targetIndex := range dipTargets {
											spot := g.neighborsRingDIP[i][j][targetIndex]
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
									if jumpRandomly {
										go func() {
											for d := 0; d < BURST_SIZE_D; d++ {
												ni := rand.Intn(GRID_SIZE)
												nj := rand.Intn(GRID_SIZE)
												g.localDips[ni][nj]++
											}
										}()
									} else {
										dipTargets := make([]int, BURST_SIZE_D)
										for d := 0; d < BURST_SIZE_D; d++ {
											dipTargets[d] = rand.Intn(len(g.neighborsRingDIP[i][j]))
										}
										go func() {
											for _, targetIndex := range dipTargets {
												spot := g.neighborsRingDIP[i][j][targetIndex]
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
					// update infected only by DIP or only by virions cells become infected by both
					if g.state[i][j] == INFECTED_VIRION || g.state[i][j] == INFECTED_DIP {

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

								// Determine the infection state based on virion and DIP infection
								if infectedByVirion && infectedByDip {
									newGrid[i][j] = INFECTED_BOTH
								} else if infectedByVirion {
									newGrid[i][j] = INFECTED_VIRION
								} else if infectedByDip {
									newGrid[i][j] = INFECTED_DIP
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

						if g.state[i][j] == INFECTED_DIP {

							g.timeSinceInfectDIP[i][j] += TIMESTEP

							if g.timeSinceInfectDIP[i][j] > IFN_DELAY+int(rand.NormFloat64()*float64(STD_IFN_DELAY)) && TAU > 0 {

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

					neighbors := g.neighbors1[i][j]

					// Iterate over the neighbors and check if any are SUSCEPTIBLE
					for _, neighbor := range neighbors {
						ni, nj := neighbor[0], neighbor[1]

						// Ensure the neighbor indices are valid (within grid bounds)
						if ni >= 0 && ni < GRID_SIZE && nj >= 0 && nj < GRID_SIZE {

							if g.timeSinceSusceptible[ni][nj]+g.timeSinceAntiviral[ni][nj] > int(rand.NormFloat64()*REGROWTH_STD+REGROWTH_MEAN) || g.timeSinceRegrowth[ni][nj]+g.timeSinceAntiviral[ni][nj] > int(rand.NormFloat64()*REGROWTH_STD+REGROWTH_MEAN) {
								newGrid[i][j] = REGROWTH
								g.timeSinceRegrowth[i][j] = 0
								g.timeSinceDead[i][j] = -1
							}

						}
					}

				}
			}
		}
		// IFN exponential decay

		if ifn_half_life != 0 {
			globalIFN = globalIFN * math.Pow(0.5, float64(TIMESTEP)/ifn_half_life)
			if globalIFN < 1e-8 {
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
		fmt.Printf("Regrowth Count: %d, Susceptible: %.2f%%, Real Susceptible: %.2f%%\n", regrowthCount, susceptiblePercentage)
		fmt.Printf("Regrowthed or Antiviral: %.2f%%, Infected: %.2f%%, DIP Only: %.2f%%, Both Infected: %.2f%%, Antiviral: %.2f%%\n",
			regrowthedOrAntiviralPercentage, infectedPercentage, infectedDIPOnlyPercentage, infectedBothPercentage, antiviralPercentage)
		fmt.Printf("Dead: %.2f%%, Uninfected: %.2f%%, Plaque: %.2f%%\n", deadCellPercentage, uninfectedPercentage, plaquePercentage)
		//fmt.Printf("Virion Diffusion Rate: %d, DIP Diffusion Rate: %d\n", virionDiffusionRate, dipDiffusionRate)

	}

	// TIMESTEP = 1 hour. If 1 hour/step, use dt = 1.0

	if virion_half_life != 0 {
		for i := 0; i < GRID_SIZE; i++ {
			for j := 0; j < GRID_SIZE; j++ {
				// 使用半衰期公式更新病毒数量
				factorV := math.Pow(0.5, float64(TIMESTEP)/virion_half_life)
				g.localVirions[i][j] = int(float64(g.localVirions[i][j])*factorV + 0.5)

				if dip_half_life != 0 {
					factorD := math.Pow(0.5, float64(TIMESTEP)/dip_half_life)
					g.localDips[i][j] = int(float64(g.localDips[i][j])*factorD + 0.5)
				}
			}
		}
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
		strconv.FormatFloat(float64(V_PFU_INITIAL)/float64(GRID_SIZE*GRID_SIZE), 'f', 6, 64),
		strconv.FormatFloat(float64(D_PFU_INITIAL)/float64(GRID_SIZE*GRID_SIZE), 'f', 6, 64),
		"-1.0",
		"-1.0",
		strconv.FormatFloat(float64(R), 'f', 6, 64),
		strconv.Itoa(BURST_SIZE_D),
		"-1.0",

		strconv.Itoa(option),
		strconv.Itoa(D_PFU_INITIAL),
		strconv.Itoa(V_PFU_INITIAL),
		strconv.Itoa(virionOnlyInfected),
		strconv.Itoa(dipOnlyInfected),
		strconv.Itoa(bothInfected),
		strconv.Itoa(totalDeadFromV),
		strconv.Itoa(totalDeadFromBoth),
		strconv.Itoa(virionDiffusionRate),
		strconv.Itoa(dipDiffusionRate),
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
		}
		fillBackground(img, color.RGBA{0, 0, 0, 255})
		for i := 0; i < GRID_SIZE; i++ {
			for j := 0; j < GRID_SIZE; j++ {
				x, y := calculateHexCenter(i, j)              // Calculate the center of each hexagon
				drawHexagon(img, x, y, colors[g.state[i][j]]) // Draw the hexagon based on the cell state
			}
		}
		// Return the image
	} else if videotype == "IFN" { // IFN concentration visualization
		black := color.RGBA{0, 0, 0, 255} // 默认颜色（黑色）

		for i := 0; i < GRID_SIZE; i++ {
			for j := 0; j < GRID_SIZE; j++ {
				x, y := calculateHexCenter(i, j) // 计算六边形中心坐标
				ifnValue := g.IFNConcentration[i][j]

				var cellColor color.RGBA
				if ifnValue <= 0 {
					cellColor = black // IFN ≤ 0，黑色
				} else if ifnValue > 0 && ifnValue <= 1 {
					cellColor = color.RGBA{0, 0, 255, 255} // 蓝色
				} else if ifnValue > 1 && ifnValue <= 2 {
					cellColor = color.RGBA{0, 255, 0, 255} // 绿色
				} else if ifnValue > 2 && ifnValue <= 5 {
					cellColor = color.RGBA{255, 255, 0, 255} // 黄色
				} else if ifnValue > 5 && ifnValue <= 10 {
					cellColor = color.RGBA{255, 165, 0, 255} // 橙色
				} else {
					cellColor = color.RGBA{255, 0, 0, 255} // 红色
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
	} else {
		fmt.Println("Error: Unknown videotype provided.")
	}

	return img // Return the image
}

func drawTextWithBackground(img *image.RGBA, x, y int, text string, color color.Color, textColor color.Color, bgColor color.Color) {
	labelWidth := len(text) * 7 // 粗略估算文字宽度
	labelHeight := 16           // 文字高度

	// 绘制背景
	for i := x; i < x+labelWidth; i++ {
		for j := y - labelHeight; j < y; j++ {
			img.Set(i, j, bgColor) // 设置背景颜色
		}
	}

	// 绘制文字
	addLabel(img, x, y, text, textColor)
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
func addStaticLegend(img *image.RGBA, gridHeight int, graphHeight int) {
	// 定义图例项和背景颜色
	legendItems := []string{
		"Yellow: infected by both",
		"Green: infected by DIP",
		"Red: infected by Virion",
		"Blue: antiviral",
		"Black: uninfected cells",
		"Grey: dead cells",
		"Purple: regrowed cells",
	}
	legendColors := map[string]color.Color{
		"Yellow: infected by both": color.RGBA{183, 149, 11, 255}, // Yellow
		"Green: infected by DIP":   color.RGBA{0, 255, 0, 255},    // Green
		"Red: infected by Virion":  color.RGBA{255, 0, 0, 255},    // Red
		"Blue: antiviral":          color.RGBA{0, 102, 255, 255},  // Blue
		"Black: uninfected cells":  color.RGBA{0, 0, 0, 255},      // Black
		"Grey: dead cells":         color.RGBA{84, 110, 122, 255}, // Grey
		"Purple: regrowed cells":   color.RGBA{128, 0, 128, 255},  // Purple
	}
	bgColor := color.RGBA{255, 255, 255, 255} // 白色背景

	// 图例起始位置
	startX := gridHeight*2 - 900 // 靠近右上角，调整横坐标偏移
	startY := 15                 // 距离顶部一定距离

	columnWidth := 200 // 每列的宽度
	lineSpacing := 20  // 每行的间距

	// 绘制图例
	for i, label := range legendItems {
		// 定义布局规则
		var x, y int
		if i == 0 {
			// 第一列只有第一个项目
			x = startX
			y = startY
		} else {
			// 第二列和第三列各放三个项目
			column := (i - 1) / 3 // 从第二个项目开始计算列
			row := (i - 1) % 3    // 计算行

			x = startX + (column+1)*columnWidth // 列偏移
			y = startY + row*lineSpacing        // 行偏移
		}

		// 设置字体颜色为图例颜色
		textColor := legendColors[label]
		drawTextWithBackground(img, x, y, label, textColor, textColor, bgColor)
	}
}

func (g *Grid) gridToImageWithGraph(frameNum int, virionOnly, dipOnly, both []float64, videotype string) *image.RGBA {
	const graphHeight = 100 // 图表固定高度
	const spacing = 0       // 图表和网格之间的间距

	// 获取网格图像和实际高度
	gridImg := g.gridToImage(videotype)
	gridHeight := gridImg.Bounds().Dy()

	// 计算画布总高度
	imgWidth := GRID_SIZE * CELL_SIZE * 2
	imgHeight := graphHeight + gridHeight + spacing

	// 创建空白画布
	canvas := image.NewRGBA(image.Rect(0, 0, imgWidth, imgHeight))

	// 将图表绘制到画布上
	graphImg := createInfectionGraph(frameNum, virionOnly, dipOnly, both)
	draw.Draw(canvas, image.Rect(0, 0, imgWidth, graphHeight), graphImg, image.Point{}, draw.Src)

	// 将网格绘制到画布上，紧接图表下面
	draw.Draw(canvas, image.Rect(0, graphHeight+spacing, imgWidth, graphHeight+gridHeight+spacing), gridImg, image.Point{}, draw.Src)

	// 添加静态图例
	addStaticLegend(canvas, GRID_SIZE*CELL_SIZE, graphHeight)

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

	// 将解析后的 flag 值赋值给全局变量（注意解引用）
	BURST_SIZE_V = *flag_burstSizeV
	BURST_SIZE_D = *flag_burstSizeD
	MEAN_LYSIS_TIME = *flag_meanLysisTime
	k_JumpR = *flag_kJumpR
	TAU = *flag_tau
	ifnBothFold = *flag_ifnBothFold
	RHO = *flag_rho

	virion_half_life = *flag_virion_half_life
	dip_half_life = *flag_dip_half_life
	ifn_half_life = *flag_ifn_half_life

	// 重新计算依赖参数（注意此时 ifnBothFold 已经是 float64，而不是 *float64）
	D_only_IFN_stimulate_ratio = 5.0 * ifnBothFold
	BOTH_IFN_stimulate_ratio = 10.0 * ifnBothFold

	// 可选：打印调试信息
	fmt.Printf("Parameters:\n  burstSizeV = %d\n  burstSizeD = %d\n  MEAN_LYSIS_TIME = %.2f\n  kJumpR = %.2f\n  TAU = %d\n  ifnBothFold = %.2f\n  RHO = %.3f\n",
		BURST_SIZE_V, BURST_SIZE_D, MEAN_LYSIS_TIME, k_JumpR, TAU, ifnBothFold, RHO)

	// --- 粒子扩散选项 ---
	particleSpreadOption = *flag_particleSpreadOption
	if particleSpreadOption == "celltocell" {
		jumpRadiusV = 0
		jumpRadiusD = 0
		jumpRandomly = false
		k_JumpR = 0.0
		allowVirionJump = false
		allowDIPJump = false
		fmt.Println("flag main celltocell")
	} else if particleSpreadOption == "jumprandomly" {
		jumpRadiusV = 0
		jumpRadiusD = 0
		jumpRandomly = true
		par_celltocell_random = false
		allowVirionJump = true
		allowDIPJump = true
		k_JumpR = 1.0 //*flag_percentRandJumpR
		fmt.Println("flag main jump randomly")
	} else if particleSpreadOption == "jumpradius" {
		jumpRadiusV = 5
		jumpRadiusD = 5
		jumpRandomly = false
		allowVirionJump = true
		allowDIPJump = true
		k_JumpR = 0.0
	} else if particleSpreadOption == "partition" {
		jumpRadiusV = 0
		jumpRadiusD = 0
		jumpRandomly = true
		par_celltocell_random = true
		k_JumpR = *flag_percentRandJumpR
	} else {
		log.Fatalf("未知的 particleSpreadOption: %s", particleSpreadOption)
	}
	fmt.Println("\nParticle spread option settings:")
	fmt.Printf("  particleSpreadOption: %s\n", particleSpreadOption)
	fmt.Printf("  jumpRadiusV: %d, jumpRadiusD: %d, jumpRandomly: %v, k_JumpR: %.2f\n",
		jumpRadiusV, jumpRadiusD, jumpRandomly, k_JumpR)

	// --- IFN 传播选项 ---
	ifnSpreadOption = *flag_ifnSpreadOption

	fmt.Printf("BBBBBefore switch: ifnSpreadOption = %q\n", ifnSpreadOption)

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
		// 禁用 IFN：将 IFN 相关参数置零
		ifnBothFold = 0.0
		// 此外在模型中可将 R、ALPHA、IFN_DELAY、STD_IFN_DELAY、tau 等置零
		ifnWave = false
		ALPHA = 0.0
		IFN_DELAY = 0
		STD_IFN_DELAY = 0
		TAU = 0
		ifn_half_life = 0.0
	default:
		log.Fatalf("未知的 ifnSpreadOption: %s", ifnSpreadOption)
		fmt.Printf("ifnSpreadOption set to: %s, IFN_wave_radius: %d\n", ifnSpreadOption, IFN_wave_radius)

	}
	fmt.Println("\nIFN spread option settings:")
	fmt.Printf("  ifnSpreadOption: %s, IFN_wave_radius: %d, ifnBothFold: %.2f\n",
		ifnSpreadOption, IFN_wave_radius, ifnBothFold)
	fmt.Printf("flag_ifnSpreadOption = %q\n", *flag_ifnSpreadOption)
	// --- DIP 选项 ---
	dipOption = *flag_dipOption
	if dipOption {
		BURST_SIZE_D = 100
		// 保持 D_only_IFN_stimulate_ratio 默认值
	} else {
		BURST_SIZE_D = 0
		D_only_IFN_stimulate_ratio = 0.0
	}
	fmt.Println("\nDIP option settings:")
	fmt.Printf("  dipOption: %v, BURST_SIZE_D: %d, D_only_IFN_stimulate_ratio: %.2f, BOTH_IFN_stimulate_ratio: %.2f\n",
		dipOption, BURST_SIZE_D, D_only_IFN_stimulate_ratio, BOTH_IFN_stimulate_ratio)

	// 此处后续可以接入仿真代码，此示例仅展示参数设置
	fmt.Println("\nSimulation initialization complete.")
	var grid Grid

	rand.Seed(42) // Seed the random number generator
	// 动态设置 R 的值
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
	case TIME_STEPS > 100:
		ticksInterval = 100.0
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

	// 调用 generateFolderName 函数生成文件夹名称
	outputFolder := generateFolderName(
		folderNumber,    // 当前文件夹编号
		jumpRandomly,    // DIP 随机跳跃逻辑
		jumpRadiusD,     // DIP 跳跃半径
		jumpRadiusV,     // Virion 跳跃半径
		BURST_SIZE_D,    // DIP 爆发大小
		BURST_SIZE_V,    // Virion 爆发大小
		V_PFU_INITIAL,   // Virion 初始值
		D_PFU_INITIAL,   // DIP 初始值
		IFN_wave_radius, // IFN 波动半径
		TAU,             // TAU 值
		TIME_STEPS,      // 时间步长
	)

	// 创建文件夹
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
		"MEAN_ANTI_TIME_Per_Cell", "STD_ANTI_TIME", "R", "DIP_BURST_PCT", "H",
		"option", "d_pfu_initial", "v_pfu_initial", "virionOnlyInfected", "dipOnlyInfected",
		"bothInfected", "totalDeadFromV", "totalDeadFromBoth", "virionDiffusionRate", "dipDiffusionRate",
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
	// 输出图像保存目录

	var extractedImages []*image.RGBA // 存储选定帧的图像

	for frameNum := 0; frameNum < TIME_STEPS; frameNum++ {

		grid.update(frameNum) // Update the grid state

		// Call the function to record infected state counts at the specific frames
		grid.recordSimulationData(writer, frameNum)

		// Convert the grid state to an image
		//img := grid.gridToImage(videotype)

		// Calculate and record the percentage of dead cells, excluding regrowth cells
		deadCellsPercentage := calculateDeadCellPercentage(grid.state)
		frameNumbers = append(frameNumbers, frameNum)                          // Record the current frame number
		deadCellPercentages = append(deadCellPercentages, deadCellsPercentage) // Record the percentage of dead cells

		// Calculate infection percentages
		virionOnly[frameNum] = float64(grid.calculateVirionOnlyInfected()) / float64(GRID_SIZE*GRID_SIZE) * 100
		dipOnly[frameNum] = float64(grid.calculateDipOnlyInfected()) / float64(GRID_SIZE*GRID_SIZE) * 100
		both[frameNum] = float64(grid.calculateBothInfected()) / float64(GRID_SIZE*GRID_SIZE) * 100

		if frameNum > 1 {
			if frameNum%24 == 0 { // 每10帧保存一次
				img := grid.gridToImageWithGraph(frameNum, virionOnly[:frameNum+1], dipOnly[:frameNum+1], both[:frameNum+1], videotype)
				extractedImages = append(extractedImages, img)
			}
		}

		// Log `y` values before feeding them to the graph
		log.Printf("Frame %d: Virion Only: %.2f%%, DIP Only: %.2f%%, Both: %.2f%%", frameNum, virionOnly[frameNum], dipOnly[frameNum], both[frameNum])
		// Generate the graph only if there are at least two frames of data
		var img *image.RGBA
		if frameNum > 0 {
			img = grid.gridToImageWithGraph(frameNum, virionOnly[:frameNum+1], dipOnly[:frameNum+1], both[:frameNum+1], videotype)
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
			savePNGImage(combinedImage, filepath.Join(outputFolder, "selected_frames_combined.png"))

		}
	}
	log.Println("Video and graph saved successfully.") // Print a success message
	fmt.Println("ifnWave is ", ifnWave)
}
