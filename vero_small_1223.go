//: A simpler script with a higher-resolution image for easier inspection. I changed the GRID_SIZE to a smaller value, like 10, and set HOWAT_V_PFU_INITIA to 1, a small initial number of virions. I also increased CELL_SIZE to 10 for a clearer image. Additionally, I changed the "option" from 3 to 2 so you can modify the initial virion location. In option 2, I set the initial location at [4][5] ( "g.localVirions[4][5]++" ) which you can change. RHO is 1, and the virus spreads to neighboring cells weighted by distance.

package main

import (
	"encoding/csv"
	"fmt"
	"os" // Used for file operations
	"path/filepath"
	"strings"

	"bytes"
	"image"
	"image/color"
	"image/draw"
	"image/jpeg"
	"image/png"
	"log"
	"math"
	"math/rand"
	"strconv"

	"github.com/icza/mjpeg"
	"github.com/wcharczuk/go-chart/v2" // Used for plotting the graph
	"github.com/wcharczuk/go-chart/v2/drawing"
	"golang.org/x/image/font"
	"golang.org/x/image/font/basicfont"
	"golang.org/x/image/math/fixed"
)

// Constant definitions
const (
	GRID_SIZE     = 50  // Size of the grid
	TIME_STEPS    = 301 // Number of time steps
	V_PFU_INITIAL = 1   // Initial number of virions (viral particles)

	FRAME_RATE    = 1          // Frame rate for the video
	OUTPUT_VIDEO  = "1203.mp4" // Output video file name
	CELL_SIZE     = 4          // The size of each hexagonal cell
	TIMESTEP      = 1          // Time step size
	D_PFU_INITIAL = 0          // Initial number of DIPs (defective interfering particles)
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

// Global variables
var (
	option    = 2        // Option for infection initialization
	videotype = "states" // color in "states" or "IFN" or "IFNonlyLargerThanZero" or "antiviralState"
	RHO       = 0.026    //0.026    //0.02468  // 0.09 Infection rate constant
	// radius 10 of grid has 331 cells
	R int
	// radius 10 of grid has 331 cells,原来知识infected cell 增加R个IFN ，
	ALPHA                      = 0.0  // Parameter for infection probability (set to 1.5)
	BURST_SIZE_V               = 50   // CHANGE Number of virions released when a cell lyses
	BURST_SIZE_D               = 100  // CHANGE 100 // Number of DIPs released when a cell lyses
	REGROWTH_MEAN              = 24.0 // Mean time for regrowth
	REGROWTH_STD               = 6.0  // Standard deviation for regrowth time
	MEAN_LYSIS_TIME            = 12.0 // Mean lysis time
	STANDARD_LYSIS_TIME        = 3.0  // Standard deviation for lysis time
	maxGlobalIFN               = -1.0 // 用于追踪最大IFN值
	globalIFN                  = -1.0 // 全局IFN浓度
	globalIFNperCell           = 0.0
	IFN_DELAY                  = 0
	STD_IFN_DELAY              = 0
	TAU                        = 0 // 95
	ifnBothFold                = 1.0
	D_only_IFN_stimulate_ratio = 0.0 // D/V *R *D_only_IFN_stimulate_ratio
	BOTH_IFN_stimulate_ratio   = 0.0 // D/V *R *D_only_IFN_stimulate_ratio

	VStimulateIFN   = false
	jumpRandomly    = true
	IFN_wave_radius = 0
	jumpRadiusV     = 5                               //CHANGE                              // Virion jump radius
	jumpRadiusD     = 5                               //CHANGE                               // DIP jump radius
	allowVirionJump = jumpRadiusV > 0 || jumpRandomly // Allow virions to jump to other cells
	allowDIPJump    = jumpRadiusD > 0 || jumpRandomly // Allow DIPs to jump to other cells

	ifnWave = IFN_wave_radius > 0

	yMax          float64
	xMax          = float64(TIME_STEPS)
	ticksInterval float64 // Interval for X-axis ticks

	adjusted_DIP_IFN_stimulate   float64
	perParticleInfectionChance_V float64
	totalDeadFromBoth            int
	totalDeadFromV               int
	virionDiffusionRate          int
	dipDiffusionRate             int

	k_JumpR = 0.5
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

// Function to initialize selectedFrames dynamically based on TIME_STEPS

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

// Function to record infected states (INFECTED_VIRION, INFECTED_DIP, INFECTED_BOTH)
func (g *Grid) recordInfectedStates(frameNum int, file *os.File) {

	file, err := os.Create("simulation_output.csv")

	if err != nil {
		log.Fatalf("Failed to create CSV file: %v", err)
	}
	defer file.Close()

	writer := csv.NewWriter(file)
	defer writer.Flush()

	// Write the CSV headers
	headers := []string{
		"Time", "Global IFN Concentration Per Cell", "Total Extracellular Virions",
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

		"MEAN_ANTI_TIME_Per_Cell", "STD_ANTI_TIME", "R", "DIP_BURST_PCT", "H", "virion_clearance_rate",
		"dip_clearance_rate", "option",
		"d_pfu_initial", "v_pfu_initial", "virionOnlyInfected", "dipOnlyInfected", "bothInfected", "totalDeadfromV", "totalDeadfromBoth", "virionDiffusionRate", "dipDiffusionRate",
	}
	err = writer.Write(headers)
	if err != nil {
		log.Fatalf("Failed to write CSV headers: %v", err)
	}

	// Write simulation data for each timestep
	for timeStep := 0; timeStep < TIME_STEPS; timeStep++ {
		g.recordSimulationData(writer, timeStep)
	}
}

// Update the state of the grid at each time step
func (g *Grid) update(frameNum int) {
	newGrid := g.state
	if ifnWave == false { // ifnWave == false

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

				if g.state[i][j] == SUSCEPTIBLE || g.state[i][j] == REGROWTH {
					// Check if the cell is infected by virions or DIPs
					if g.localVirions[i][j] > 0 || g.localDips[i][j] > 0 {
						// Calculate the infection probabilities

						perParticleInfectionChance := RHO * math.Exp(-ALPHA*(float64(globalIFNperCell)))

						var probabilityVInfection, probabilityDInfection float64

						// Virion infection probability
						probabilityVInfection = 1 - math.Pow(1-perParticleInfectionChance, float64(g.localVirions[i][j]))
						infectedByVirion := rand.Float64() <= probabilityVInfection

						// DIP infection probability
						probabilityDInfection = 1 - math.Pow(1-perParticleInfectionChance, float64(g.localDips[i][j]))
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

				}

			}
		}

		allowRandomly := make([][]bool, GRID_SIZE)
		for i := range allowRandomly {
			allowRandomly[i] = make([]bool, GRID_SIZE)
		}

		// 根据 k_JumpR 计算允许随机跳跃的细胞总数
		totalCells := GRID_SIZE * GRID_SIZE
		randomJumpCells := int(float64(totalCells) * k_JumpR)

		// 随机选择 randomJumpCells 个细胞，标记为 allowRandomly
		selectedCells := make(map[[2]int]bool)
		for len(selectedCells) < randomJumpCells {
			ni := rand.Intn(GRID_SIZE)
			nj := rand.Intn(GRID_SIZE)
			selectedCells[[2]int{ni, nj}] = true
		}
		for pos := range selectedCells {
			allowRandomly[pos[0]][pos[1]] = true
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

							// After lysis, the cell becomes DEAD and virions and DIPs are spread to neighbors
							newGrid[i][j] = DEAD

							g.timeSinceDead[i][j] = 0
							g.timeSinceInfectVorBoth[i][j] = -1
							g.timeSinceInfectDIP[i][j] = -1

							if allowRandomly[i][j] {
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
						}
					}
					// update infected only by DIP or only by virions cells become infected by both
					if g.state[i][j] == INFECTED_VIRION || g.state[i][j] == INFECTED_DIP {

						if g.stateChanged[i][j] == false {
							// Check if the cell is infected by virions or DIPs
							if g.localVirions[i][j] > 0 || g.localDips[i][j] > 0 {
								// Calculate the infection probabilities
								perParticleInfectionChance := RHO * math.Exp(-ALPHA*(globalIFNperCell))

								var probabilityVInfection, probabilityDInfection float64

								// Virion infection probability
								probabilityVInfection = 1 - math.Pow(1-perParticleInfectionChance, float64(g.localVirions[i][j]))
								infectedByVirion := rand.Float64() <= probabilityVInfection

								// DIP infection probability
								probabilityDInfection = 1 - math.Pow(1-perParticleInfectionChance, float64(g.localDips[i][j]))
								infectedByDip := rand.Float64() <= probabilityDInfection

								// Determine the infection state based on virion and DIP infection
								if g.state[i][j] == INFECTED_VIRION && infectedByDip || g.state[i][j] == INFECTED_DIP && infectedByVirion {
									newGrid[i][j] = INFECTED_BOTH
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
	virionDiffusionRate, dipDiffusionRate := g.calculateDiffusionRates()

	// Log additional data as necessary
	fmt.Printf("Regrowth Count: %d, Susceptible: %.2f%%, Real Susceptible: %.2f%%\n", regrowthCount, susceptiblePercentage)
	fmt.Printf("Regrowthed or Antiviral: %.2f%%, Infected: %.2f%%, DIP Only: %.2f%%, Both Infected: %.2f%%, Antiviral: %.2f%%\n",
		regrowthedOrAntiviralPercentage, infectedPercentage, infectedDIPOnlyPercentage, infectedBothPercentage, antiviralPercentage)
	fmt.Printf("Dead: %.2f%%, Uninfected: %.2f%%, Plaque: %.2f%%\n", deadCellPercentage, uninfectedPercentage, plaquePercentage)
	fmt.Printf("Virion Diffusion Rate: %d, DIP Diffusion Rate: %d\n", virionDiffusionRate, dipDiffusionRate)

}

// Function to record simulation data into CSV at each timestep
func (g *Grid) recordSimulationData(writer *csv.Writer, frameNum int) {
	totalVirions := g.totalVirions()
	totalDIPs := g.totalDIPs()
	deadCellPercentage := strconv.FormatFloat(calculateDeadCellPercentage(g.state), 'f', 6, 64)
	susceptiblePercentage := strconv.FormatFloat(g.calculateSusceptiblePercentage(), 'f', 6, 64)
	infectedPercentage := strconv.FormatFloat(g.calculateInfectedPercentage(), 'f', 6, 64)
	infectedDIPOnlyPercentage := strconv.FormatFloat(g.calculateInfectedDIPOnlyPercentage(), 'f', 6, 64)
	infectedBothPercentage := strconv.FormatFloat(g.calculateInfectedBothPercentage(), 'f', 6, 64)
	antiviralPercentage := strconv.FormatFloat(g.calculateAntiviralPercentage(), 'f', 6, 64) // Placeholder; modify if antiviral state is modeled
	virionOnlyInfected := g.calculateVirionOnlyInfected()
	dipOnlyInfected := g.calculateDipOnlyInfected()
	bothInfected := g.calculateBothInfected()

	row := []string{
		strconv.Itoa(frameNum), //"Time"
		strconv.FormatFloat(globalIFN/float64(GRID_SIZE*GRID_SIZE), 'f', 6, 64), // "Global IFN Concentration Per Cell"
		strconv.Itoa(totalVirions),               // "Total Extracellular Virions"
		strconv.Itoa(totalDIPs),                  //"Total Extracellular DIPs"
		deadCellPercentage,                       // Percentage Dead Cells
		susceptiblePercentage,                    // Percentage Susceptible Cells
		infectedPercentage,                       // Percentage Infected Cells
		infectedDIPOnlyPercentage,                // Percentage Infected DIP-only Cells
		infectedBothPercentage,                   // Percentage Infected Both Cells
		antiviralPercentage,                      // Percentage Antiviral Cells
		strconv.Itoa(g.calculateRegrowthCount()), // Regrowth Count
		strconv.FormatFloat(g.calculateSusceptiblePercentage(), 'f', 6, 64),           //"Percentage Susceptible and Antiviral (Real Susceptible cells without regrowthed ones) Cells"
		strconv.FormatFloat(g.calculateRegrowthedOrAntiviralPercentage(), 'f', 6, 64), // "Percentage Regrowthed or Regrowthed and Antiviral Cells",
		"variate, depending on radius 10 of IFN",                                      // "Probability Virion Infection"
		"variate, depending on radius 10 of IFN",                                      // "Probability DIP Infection"
		strconv.FormatFloat(RHO, 'f', 6, 64),                                          // "Per Particle Infection Chance RHO"
		strconv.Itoa(totalVirions + totalDIPs),                                        // "Total Local Particles"
		strconv.FormatFloat(g.calculatePlaquePercentage(), 'f', 6, 64),                // "Plaque Percentage"
		strconv.FormatFloat(float64(maxGlobalIFN), 'f', 6, 64),                        // Max Global IFN (placeholder, not used in this simulation)
		"-1.0", // Time when all cells are uninfected (placeholder)
		strconv.FormatFloat(g.calculateUninfectedPercentage(), 'f', 6, 64), // Percentage Uninfected Cells
		"0",                                                  // Number of Plaques (not calculated, placeholder)
		strconv.Itoa(GRID_SIZE),                              // GRID_SIZE
		strconv.Itoa(TIMESTEP),                               // TIMESTEP
		strconv.Itoa(IFN_DELAY),                              // IFN_DELAY (placeholder)
		strconv.Itoa(STD_IFN_DELAY),                          // STD_IFN_DELAY (placeholder)
		strconv.FormatFloat(ALPHA, 'f', 6, 64),               // ALPHA
		strconv.FormatFloat(RHO, 'f', 6, 64),                 // RHO
		strconv.FormatFloat(float64(TAU), 'f', 6, 64),        // TAU (placeholder)
		strconv.Itoa(BURST_SIZE_V),                           // BURST_SIZE for Virions
		strconv.FormatFloat(REGROWTH_MEAN, 'f', 6, 64),       // REGROWTH_MEAN
		strconv.FormatFloat(REGROWTH_STD, 'f', 6, 64),        // REGROWTH_STD
		strconv.Itoa(TIME_STEPS),                             // TIME_STEPS
		strconv.FormatFloat(MEAN_LYSIS_TIME, 'f', 6, 64),     // MEAN_LYSIS_TIME
		strconv.FormatFloat(STANDARD_LYSIS_TIME, 'f', 6, 64), // STANDARD_LYSIS_TIME
		strconv.FormatFloat(float64(V_PFU_INITIAL)/float64(GRID_SIZE*GRID_SIZE), 'f', 6, 64), // Initial Virion PFU per cell (placeholder)
		strconv.FormatFloat(float64(D_PFU_INITIAL)/float64(GRID_SIZE*GRID_SIZE), 'f', 6, 64), // Initial DIP PFU per cell (placeholder)
		"-1.0", // Mean Antiviral Time (placeholder)
		"-1.0", // STD Antiviral Time (placeholder)
		strconv.FormatFloat(float64(R), 'f', 6, 64), // R (placeholder, optional parameter for specific models)
		strconv.Itoa(BURST_SIZE_D),                  // DIP Burst Percentage (placeholder)
		"-1.0",                                      // H (placeholder)
		"-1.0",                                      // Virion Clearance Rate (placeholder)
		"-1.0",                                      // DIP Clearance Rate (placeholder)
		strconv.Itoa(option),                        // Infection Initialization Option

		strconv.Itoa(D_PFU_INITIAL),      // Initial DIP PFU
		strconv.Itoa(V_PFU_INITIAL),      // Initial Virion PFU
		strconv.Itoa(virionOnlyInfected), // Number of virion-only infected cells
		strconv.Itoa(dipOnlyInfected),    // Number of DIP-only infected cells
		strconv.Itoa(bothInfected),       // Number of both-infected cells

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

		red := color.RGBA{255, 0, 0, 255}   // Cells with interferon > 0
		green := color.RGBA{0, 255, 0, 255} // Cells in antiviral state exceeding duration
		black := color.RGBA{0, 0, 0, 255}   // Default color for all other cells
		yellow := color.RGBA{255, 255, 0, 255}

		for i := 0; i < GRID_SIZE; i++ {
			for j := 0; j < GRID_SIZE; j++ {
				x, y := calculateHexCenter(i, j) // Calculate the center of each hexagon

				// Apply color based on the specified conditions
				if g.state[i][j] == INFECTED_VIRION || g.state[i][j] == INFECTED_DIP || g.state[i][j] == INFECTED_BOTH {
					drawHexagon(img, x, y, yellow) // Yellow for infected cells
				} else if g.timeSinceAntiviral[i][j] > g.antiviralDuration[i][j] {
					drawHexagon(img, x, y, red) // Green for cells in antiviral state exceeding duration
				} else if g.IFNConcentration[i][j] > 0 {
					drawHexagon(img, x, y, green) // Red for cells with interferon > 0
				} else {
					drawHexagon(img, x, y, black) // Black for all other cells
				}
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
		yMax = 35.0
	case IFN_wave_radius == 10 && TAU == 12 && jumpRandomly == true:
		yMax = 35.0
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
		yMax = 35.0
	case IFN_wave_radius == 10 && jumpRadiusV == 5 && jumpRadiusD == 5 && TAU == 24:
		yMax = 0.2

	case IFN_wave_radius == 0 && jumpRadiusV == 0 && jumpRadiusD == 0:
		yMax = 0.3
	case IFN_wave_radius == 0 && jumpRadiusV == 5 && jumpRadiusD == 5:
		yMax = 35.0
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
		"Time", "Global IFN Concentration Per Cell", "Total Extracellular Virions",
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
		"MEAN_ANTI_TIME_Per_Cell", "STD_ANTI_TIME", "R", "DIP_BURST_PCT", "H", "virion_clearance_rate",
		"dip_clearance_rate", "option",
		"d_pfu_initial", "v_pfu_initial", "virionOnlyInfected", "dipOnlyInfected", "bothInfected", "totalDeadFromV", "totalDeadFromBoth", "virionDiffusionRate", "dipDiffusionRate",
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

		if frameNum > 31 {
			if frameNum%50 == 0 { // 每10帧保存一次
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
			log.Printf("Selected frames combined and saved.")
		}
	}
	log.Println("Video and graph saved successfully.") // Print a success message
	fmt.Println("ifnWave is ", ifnWave)
}
