package main

import (
	"encoding/csv"
	"fmt"
	"os" // Used for file operations

	"bytes"
	"image"
	"image/color"
	"image/jpeg"
	"log"
	"math"
	"math/rand"
	"strconv"

	"github.com/icza/mjpeg"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/plotutil"
	"gonum.org/v1/plot/vg"
)

// Constant definitions
const (
	GRID_SIZE     = 360                // Size of the grid
	TIME_STEPS    = 145                // Number of time steps
	V_PFU_INITIAL = 500                // Initial number of virions (viral particles)
	D_PFU_INITIAL = 0                  // Initial number of DIPs (defective interfering particles)
	FRAME_RATE    = 5                  // Frame rate for the video
	OUTPUT_VIDEO  = "dip_mbdk_360.mp4" // Output video file name
	CELL_SIZE     = 1                  // The size of each hexagonal cell
	TIMESTEP      = 1                  // Time step size
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
	option              = 3     // Option for infection initialization
	RHO                 = 0.015 // 0.09 Infection rate constant
	R                   = 1
	ALPHA               = 1.5                         // Parameter for infection probability (set to 1.5)
	BURST_SIZE_V        = 100                         // Number of virions released when a cell lyses
	BURST_SIZE_D        = 0                           // Number of DIPs released when a cell lyses
	REGROWTH_MEAN       = 24.0                        // Mean time for regrowth
	REGROWTH_STD        = 6.0                         // Standard deviation for regrowth time
	MEAN_LYSIS_TIME     = 12.0                        // Mean lysis time
	STANDARD_LYSIS_TIME = 3.0                         // Standard deviation for lysis time
	maxGlobalIFN        = math.SmallestNonzeroFloat64 // 用于追踪最大IFN值
	globalIFN           = -1.0                        // 全局IFN浓度
	DIP_IFN_stimulate   = 0.0                         // DIPs对IFN的刺激
	IFN_DELAY           = 5
	STD_IFN_DELAY       = 1
	TAU                 = 120 // 95
)

// Grid structure for storing the simulation state
type Grid struct {
	state                [GRID_SIZE][GRID_SIZE]int       // State of the cells in the grid
	localVirions         [GRID_SIZE][GRID_SIZE]int       // Number of virions in each cell
	localDips            [GRID_SIZE][GRID_SIZE]int       // Number of DIPs in each cell
	IFNConcentration     [GRID_SIZE][GRID_SIZE]float64   // IFN concentration in each cell
	timeSinceInfection   [GRID_SIZE][GRID_SIZE]int       // Time since infection for each cell
	timeSinceDead        [GRID_SIZE][GRID_SIZE]int       // Time since death for each cell
	timeSinceRegrowth    [GRID_SIZE][GRID_SIZE]int       // Time since regrowth for each cell
	timeSinceSusceptible [GRID_SIZE][GRID_SIZE]int       // Time since cell became susceptible
	neighbors1           [GRID_SIZE][GRID_SIZE][6][2]int // Neighbors at distance 1
	neighbors2           [GRID_SIZE][GRID_SIZE][6][2]int // Neighbors at distance 2
	neighbors3           [GRID_SIZE][GRID_SIZE][6][2]int // Neighbors at distance 3
	stateChanged         [GRID_SIZE][GRID_SIZE]bool      // Flag to indicate if the state of a cell has changed
	antiviralDuration    [GRID_SIZE][GRID_SIZE]int       // Duration of antiviral state
	previousStates       [GRID_SIZE][GRID_SIZE]int       // Previous state of the cell
	antiviralFlag        [GRID_SIZE][GRID_SIZE]bool      // Flag to indicate if the cell is in the antiviral state
	timeSinceAntiviral   [GRID_SIZE][GRID_SIZE]int       // Time since the cell entered the antiviral state
	antiviralCellCount   int                             // Number of cells in the antiviral state
	totalAntiviralTime   int
}

// Initialize the grid, setting all cells to SUSCEPTIBLE
func (g *Grid) initialize() {
	for i := 0; i < GRID_SIZE; i++ {
		for j := 0; j < GRID_SIZE; j++ {
			g.state[i][j] = SUSCEPTIBLE
			g.stateChanged[i][j] = false // Initialize as unchanged
			g.timeSinceInfection[i][j] = -1
			g.timeSinceDead[i][j] = -1
			g.timeSinceRegrowth[i][j] = -1
			g.IFNConcentration[i][j] = 0
			g.antiviralDuration[i][j] = -1
			g.timeSinceSusceptible[i][j] = 0
			g.previousStates[i][j] = -1
			g.antiviralFlag[i][j] = false
			g.timeSinceAntiviral[i][j] = -2
		}
	}
	fmt.Println("Grid initialized")
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
			}
		}
	}
	return regrowthCells
}

func (g *Grid) calculatemeanAntiviralTimePerCell() float64 {
	if g.antiviralCellCount == 0 {
		return 0 // Avoid division by zero
	}
	return float64(g.totalAntiviralTime) / float64(g.antiviralCellCount)
}

// Function to calculate the percentage of susceptible cells in the grid
func (g *Grid) calculateSusceptiblePercentage() float64 {
	totalCells := GRID_SIZE * GRID_SIZE
	susceptibleCells := 0
	for i := 0; i < GRID_SIZE; i++ {
		for j := 0; j < GRID_SIZE; j++ {
			if g.state[i][j] == SUSCEPTIBLE {
				susceptibleCells++
			}
		}
	}
	return (float64(susceptibleCells) / float64(totalCells)) * 100
}

// calculateSusceptibleIncludingRegrowthPercentag
func (g *Grid) calculateRealSusceptiblePercentage() float64 {
	totalCells := GRID_SIZE * GRID_SIZE
	susceptibleCells := 0
	for i := 0; i < GRID_SIZE; i++ {
		for j := 0; j < GRID_SIZE; j++ {
			// Include all cells that are in the SUSCEPTIBLE state, regardless of regrowth
			if g.state[i][j] == SUSCEPTIBLE {
				susceptibleCells++
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

// Calculate neighbor relationships
func (g *Grid) initializeNeighbors() {
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
		"Percentage Infected Cells", "Percentage Infected DIP only Cells",
		"Percentage Infected Both Cells", "Percentage Antiviral Cells",
		"Regrowth Count",
		"Percentage Susceptible and Antiviral (Real Susceptible cells without regrowthed ones) Cells",
		"Percentage Regrowthed or Regrowthed and Antiviral Cells",
		"Probability Virion Infection", "Probability DIP Infection",
		"Per Particle Infection Chance", "Total Local Particles",
		"Plaque Percentage", "max_global_IFN", "time_all_cells_uninfected",
		"Percentage Uninfected Cells", "num_plaques", "GRID_SIZE", "TIMESTEP",
		"IFN_DELAY", "STD_IFN_DELAY", "ALPHA", "RHO", "TAU", "BURST_SIZE",
		"REGROWTH_MEAN", "REGROWTH_STD", "TIME_STEPS", "MEAN_LYSIS_TIME",
		"STANDARD_LYSIS_TIME", "init_v_pfu_per_cell", "init_d_pfu_per_cell",
		"MEAN_ANTI_TIME", "STD_ANTI_TIME", "R", "DIP_BURST_PCT", "H", "virion_clearance_rate",
		"dip_clearance_rate", "option",
		"d_pfu_initial", "v_pfu_initial", "virionOnlyInfected", "dipOnlyInfected", "bothInfected",
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

// Initialize the infection state
func (g *Grid) initializeInfection(option int) {
	switch option {
	case 2:
		// Manual initialization
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

// Update the state of the grid at each time step
func (g *Grid) update(frameNum int) {
	newGrid := g.state

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
					perParticleInfectionChance := RHO * math.Exp(-ALPHA*(float64(globalIFN)/float64(R)))

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
					} else if infectedByVirion {
						newGrid[i][j] = INFECTED_VIRION
					} else if infectedByDip {
						newGrid[i][j] = INFECTED_DIP
					}
				}

				// Mark the state as changed if the cell is infected
				if newGrid[i][j] != g.state[i][j] {
					g.stateChanged[i][j] = true
				}
			}

			// Antiviral logic
			if g.IFNConcentration[i][j] > 0 {

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
		}
	}

	// Process infected cells
	for i := 0; i < GRID_SIZE; i++ {
		for j := 0; j < GRID_SIZE; j++ {
			if g.state[i][j] == INFECTED_VIRION || g.state[i][j] == INFECTED_DIP || g.state[i][j] == INFECTED_BOTH {
				g.timeSinceInfection[i][j] += TIMESTEP

				// Check if the cell should lyse and release virions and DIPs
				if g.timeSinceInfection[i][j] > int(rand.NormFloat64()*STANDARD_LYSIS_TIME+MEAN_LYSIS_TIME) {

					// After lysis, the cell becomes DEAD and virions and DIPs are spread to neighbors
					newGrid[i][j] = DEAD
					g.timeSinceDead[i][j] = 0
					g.timeSinceInfection[i][j] = -1

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

					// Calculate the number of virions and DIPs assigned to each type of neighbor
					virionsForNeighbors1 := int(float64(BURST_SIZE_V) * (ratio1 * float64(len(g.neighbors1[i][j]))) / totalRatio)
					virionsForNeighbors2 := int(float64(BURST_SIZE_V) * (ratio2 * float64(len(g.neighbors2[i][j]))) / totalRatio)
					virionsForNeighbors3 := int(float64(BURST_SIZE_V) * (ratio3 * float64(len(g.neighbors3[i][j]))) / totalRatio)

					dipsForNeighbors1 := int(float64(BURST_SIZE_D) * (ratio1 * float64(len(g.neighbors1[i][j]))) / totalRatio)
					dipsForNeighbors2 := int(float64(BURST_SIZE_D) * (ratio2 * float64(len(g.neighbors2[i][j]))) / totalRatio)
					dipsForNeighbors3 := int(float64(BURST_SIZE_D) * (ratio3 * float64(len(g.neighbors3[i][j]))) / totalRatio)

					// Calculate the remaining virions and DIPs
					remainingVirions := BURST_SIZE_V - (virionsForNeighbors1 + virionsForNeighbors2 + virionsForNeighbors3)
					remainingDips := BURST_SIZE_D - (dipsForNeighbors1 + dipsForNeighbors2 + dipsForNeighbors3)

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

				if g.timeSinceInfection[i][j] > IFN_DELAY+int(rand.NormFloat64()*float64(STD_IFN_DELAY)) {

					g.IFNConcentration[i][j] += float64(R) * float64(TIMESTEP)
					globalIFN += g.IFNConcentration[i][j]

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
						if g.state[ni][nj] == SUSCEPTIBLE || g.state[ni][nj] == REGROWTH || g.state[ni][nj] == ANTIVIRAL {
							canRegrow = true
							break // Exit the loop as soon as a susceptible neighbor is found
						}
					}
				}

				// If the conditions are met, the cell regrows
				if canRegrow && g.timeSinceDead[i][j] >= int(rand.NormFloat64()*REGROWTH_STD+REGROWTH_MEAN) {
					newGrid[i][j] = REGROWTH
					g.timeSinceDead[i][j] = -1

				}
			}
		}
	}
	globalIFN = globalIFN / float64(GRID_SIZE*GRID_SIZE)

	// Apply the updated grid state
	g.state = newGrid

	// Calculate and log the total virions and DIPs for each time step
	totalVirions, totalDIPs := g.totalVirions(), g.totalDIPs()
	fmt.Printf("Time step %d: Total Virions = %d, Total DIPs = %d\n", frameNum, totalVirions, totalDIPs)

	// Additional calculations based on simulation parameters for tracking purposes
	regrowthCount := g.calculateRegrowthCount()
	susceptiblePercentage := g.calculateSusceptiblePercentage()
	realSusceptiblePercentage := g.calculateRealSusceptiblePercentage()
	regrowthedOrAntiviralPercentage := g.calculateRegrowthedOrAntiviralPercentage()
	infectedPercentage := g.calculateInfectedPercentage()
	infectedDIPOnlyPercentage := g.calculateInfectedDIPOnlyPercentage()
	infectedBothPercentage := g.calculateInfectedBothPercentage()
	antiviralPercentage := g.calculateAntiviralPercentage()
	deadCellPercentage := calculateDeadCellPercentage(g.state)
	uninfectedPercentage := g.calculateUninfectedPercentage()
	plaquePercentage := g.calculatePlaquePercentage()

	// Log additional data as necessary
	fmt.Printf("Regrowth Count: %d, Susceptible: %.2f%%, Real Susceptible: %.2f%%\n", regrowthCount, susceptiblePercentage, realSusceptiblePercentage)
	fmt.Printf("Regrowthed or Antiviral: %.2f%%, Infected: %.2f%%, DIP Only: %.2f%%, Both Infected: %.2f%%, Antiviral: %.2f%%\n",
		regrowthedOrAntiviralPercentage, infectedPercentage, infectedDIPOnlyPercentage, infectedBothPercentage, antiviralPercentage)
	fmt.Printf("Dead: %.2f%%, Uninfected: %.2f%%, Plaque: %.2f%%\n", deadCellPercentage, uninfectedPercentage, plaquePercentage)
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
	meanAntiviralTimePerCell := g.calculatemeanAntiviralTimePerCell()
	//probVirionInfection := RHO * math.Exp(-ALPHA*(globalIFN/float64(R)))
	//probDIPInfection := RHO * math.Exp(-ALPHA*(globalIFN/float64(R)))

	row := []string{
		strconv.Itoa(frameNum), // Time step
		strconv.FormatFloat(globalIFN/float64(GRID_SIZE*GRID_SIZE), 'f', 6, 64), // Global IFN Concentration Per Cell
		strconv.Itoa(totalVirions),               // Total Extracellular Virions
		strconv.Itoa(totalDIPs),                  // Total Extracellular DIPs
		deadCellPercentage,                       // Percentage Dead Cells
		susceptiblePercentage,                    // Percentage Susceptible Cells
		infectedPercentage,                       // Percentage Infected Cells
		infectedDIPOnlyPercentage,                // Percentage Infected DIP-only Cells
		infectedBothPercentage,                   // Percentage Infected Both Cells
		antiviralPercentage,                      // Percentage Antiviral Cells
		strconv.Itoa(g.calculateRegrowthCount()), // Regrowth Count
		strconv.FormatFloat(g.calculateRealSusceptiblePercentage(), 'f', 6, 64),       // Real Susceptible Cells without regrowth
		strconv.FormatFloat(g.calculateRegrowthedOrAntiviralPercentage(), 'f', 6, 64), // Regrowthed or Antiviral Cells
		"-1.0",                                 // Probability Virion Infection (placeholder)
		"-1.0",                                 // Probability DIP Infection (placeholder)
		strconv.FormatFloat(RHO, 'f', 6, 64),   // Per Particle Infection Chance
		strconv.Itoa(totalVirions + totalDIPs), // Total Local Particles (sum of virions and DIPs)
		strconv.FormatFloat(g.calculatePlaquePercentage(), 'f', 6, 64), // Plaque Percentage
		strconv.FormatFloat(float64(maxGlobalIFN), 'f', 6, 64),         // Max Global IFN (placeholder, not used in this simulation)
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
		strconv.FormatFloat(meanAntiviralTimePerCell, 'f', 6, 64),                            // Mean Antiviral Time (placeholder)
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
func (g *Grid) gridToImage() *image.RGBA {
	imgWidth := GRID_SIZE * CELL_SIZE * 2                       // Calculate the image width
	imgHeight := GRID_SIZE * CELL_SIZE * 2                      // Calculate the image height
	img := image.NewRGBA(image.Rect(0, 0, imgWidth, imgHeight)) // Create a new image

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

	for i := 0; i < GRID_SIZE; i++ {
		for j := 0; j < GRID_SIZE; j++ {
			x, y := calculateHexCenter(i, j)              // Calculate the center of each hexagon
			drawHexagon(img, x, y, colors[g.state[i][j]]) // Draw the hexagon based on the cell state
		}
	}
	return img // Return the image
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

	grid.initialize()                // Initialize the grid
	grid.initializeNeighbors()       // Initialize the neighbors
	grid.initializeInfection(option) // Initialize the infection state

	// Open a CSV file to record the infected states over time
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
		"IFN_DELAY", "STD_IFN_DELAY", "ALPHA", "RHO", "TAU", "BURST_SIZE",
		"REGROWTH_MEAN", "REGROWTH_STD", "TIME_STEPS", "MEAN_LYSIS_TIME",
		"STANDARD_LYSIS_TIME", "init_v_pfu_per_cell", "init_d_pfu_per_cell",
		"MEAN_ANTI_TIME_Per_Cell", "STD_ANTI_TIME", "R", "DIP_BURST_PCT", "H", "virion_clearance_rate",
		"dip_clearance_rate", "option",
		"d_pfu_initial", "v_pfu_initial", "virionOnlyInfected", "dipOnlyInfected", "bothInfected",
	}
	err = writer.Write(headers)
	if err != nil {
		log.Fatalf("Failed to write CSV headers: %v", err)
	}

	// Create a file for logging the total virions and DIPs counts
	virionDipFile, err := os.Create("virion_dip_counts.txt")
	if err != nil {
		log.Fatalf("Failed to create virion and DIP counts file: %v", err)
	}
	defer virionDipFile.Close()

	// Create an MJPEG video writer
	videoWriter, err := mjpeg.New(OUTPUT_VIDEO, int32(GRID_SIZE*CELL_SIZE*2), int32(GRID_SIZE*CELL_SIZE*2), int32(FRAME_RATE))
	if err != nil {
		log.Fatalf("Failed to create MJPEG writer: %v", err) // Handle the error if the writer fails to create
	}
	defer videoWriter.Close() // Ensure the writer is closed when the program ends

	var buf bytes.Buffer                      // Buffer for JPEG encoding
	jpegOptions := &jpeg.Options{Quality: 75} // JPEG encoding options, quality set to 75

	var frameNumbers []int            // Slice to store frame numbers
	var deadCellPercentages []float64 // Slice to store dead cell percentages

	// Create a file to store results, using RHO value in the file name
	fileName := fmt.Sprintf("DeadCell_RHO_%.2f.txt", RHO)
	deadCellFile, err := os.Create(fileName)
	if err != nil {
		log.Fatalf("Failed to create file: %v", err) // Handle error if file creation fails
	}
	defer deadCellFile.Close() // Ensure the file is closed when the program ends

	for frameNum := 0; frameNum < TIME_STEPS; frameNum++ {
		grid.update(frameNum) // Update the grid state

		// Call the function to record infected state counts at the specific frames
		grid.recordSimulationData(writer, frameNum)

		// Calculate and record the total virions and DIPs for each time step
		totalVirions, totalDIPs := grid.totalVirions(), grid.totalDIPs()
		_, err := virionDipFile.WriteString(fmt.Sprintf("Timestep %d: Total Virions = %d, Total DIPs = %d\n", frameNum, totalVirions, totalDIPs))
		if err != nil {
			log.Fatalf("Failed to write to virion/dip counts file: %v", err)
		}

		// Convert the grid state to an image
		img := grid.gridToImage()

		// Calculate and record the percentage of dead cells, excluding regrowth cells
		deadCellsPercentage := calculateDeadCellPercentage(grid.state)
		frameNumbers = append(frameNumbers, frameNum)                          // Record the current frame number
		deadCellPercentages = append(deadCellPercentages, deadCellsPercentage) // Record the percentage of dead cells

		// Output the percentage of dead cells at specific time steps
		if frameNum == 24 || frameNum == 48 || frameNum == 72 || frameNum == 96 || frameNum == 120 || frameNum == 144 {
			log.Printf("Timestep %d: Percentage of Dead Cells = %.2f%%\n", frameNum, deadCellsPercentage)
			_, err := deadCellFile.WriteString(fmt.Sprintf("Timestep %d: Percentage of Dead Cells = %.2f%%\n", frameNum, deadCellsPercentage))
			if err != nil {
				log.Fatalf("Failed to write to file: %v", err) // Handle error if writing to file fails
			}

		}

		// Encode the image to JPEG format
		err = jpeg.Encode(&buf, img, jpegOptions)
		if err != nil {
			log.Fatalf("Failed to encode image: %v", err) // Handle error if image encoding fails
		}

		// Add the frame to the video
		err = videoWriter.AddFrame(buf.Bytes())
		if err != nil {
			log.Fatalf("Failed to add frame: %v", err) // Handle error if adding frame to video fails
		}
		buf.Reset() // Reset the buffer for the next frame
	}

	// Create and save the graph showing the percentage of dead cells over time
	p := plot.New()

	p.Title.Text = "Percentage of Dead Cells Over Time" // Set the title of the graph
	p.X.Label.Text = "Timestep"                         // Set the label for the x-axis
	p.Y.Label.Text = "Percentage of Dead Cells"         // Set the label for the y-axis

	// Customize x-axis tick marks, marking every 24 time steps
	p.X.Tick.Marker = plot.TickerFunc(func(min, max float64) []plot.Tick {
		var ticks []plot.Tick
		for i := 0.0; i <= max; i += 24 {
			ticks = append(ticks, plot.Tick{Value: i, Label: strconv.Itoa(int(i))})
		}
		return ticks
	})

	points := make(plotter.XYs, len(frameNumbers)) // Prepare data points for the plot
	for i := range points {
		points[i].X = float64(frameNumbers[i]) // Set x-axis data to frame number
		points[i].Y = deadCellPercentages[i]   // Set y-axis data to percentage of dead cells
	}

	// Add the data points to the plot
	err = plotutil.AddLinePoints(p, "Dead Cells", points)
	if err != nil {
		log.Fatalf("Failed to add data points to plot: %v", err) // Handle error if adding data points fails
	}

	// Save the plot as an image file
	if err := p.Save(8*vg.Inch, 4*vg.Inch, "dead_cells_percentage_mbdk_360.png"); err != nil {
		log.Fatalf("Failed to save plot: %v", err) // Handle error if saving the plot fails
	}

	log.Println("Video and graph saved successfully.") // Print a success message
}
