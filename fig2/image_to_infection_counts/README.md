# Image to Infection Counts

I want to convert experimental fluorescence microscopy images into overlay imagesn ,then by color, convert into infection count CSV data.

References:
Liang Q, Yang J, Fan WTL and Lo WC. Patch Formation Driven by Stochastic Effects of Interaction
between Viruses and Defective Interfering Particles. PLoS Comput Biol. 2023; 19:e1011513. (doi:10.
1371/journal.pcbi.1011513



## Final Output
`13_infection_counts_by_time.csv`:
```
time,virion_counts,dip_counts,both_infected_counts,susceptible_counts,total_cells,virion_rate,dip_rate,both_infected_rate,susceptible_rate
```

## Source Images (Original Experimental Data)
- `7_hours.png` - 7 hours post infection
- `13_hours.png` - 13 hours post infection  
- `19_hours.png` - 19 hours post infection
- `25_hours_comprehensive.png` - 25 hours (all cell types)

comprehensive figure is combined (at t=25) from:
- `25_hours_DIP-infected.png` - 25 hours (DIP-infected cells only)
- `25_hours_Virion-infected.png` - 25 hours (Virion-infected cells only)

## Make into hexigon ; Processed Overlay Images
These show hexagonal grid analysis with color-coded cell classification:
- `2_7_hours_adjusted_overlay.png` 
- `2_13_hours_adjusted_overlay.png` 
- `2_19_hours_adjusted_overlay.png` 
- `2_25_hours_comprehensive_adjusted_overlay.png`
- `2_25_hours_DIP-infected_adjusted_overlay.png`
- `2_25_hours_Virion-infected_adjusted_overlay.png`


### `2_adjusted_calibrated_overlay.py`
Main image processing script that generates the overlay images:
- Loads PNG images
- Creates hexagonal grid overlay
- Detects cell colors (R=Virion, G=DIP, Y=Both, B=Susceptible)
- Generates overlay visualizations and CSV files

### `14_generate_latest_csv.py`
Generates the final adjusted `13_infection_counts_by_time.csv` with calibrated values.

Run: `python 14_generate_latest_csv.py`

### `1_final_calibrated_overlay.py`
Alternative/earlier version of overlay generator (produces different counts).

## Color Legend
- red: Virion-only infected cells
- green: DIP-only infected cells  
- yellow: Both virion and DIP infected cells
- black: Susceptible (uninfected) cells


