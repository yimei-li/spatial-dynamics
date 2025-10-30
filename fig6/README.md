# Figure 6: IFN1-10 Integrated Analysis

## Overview
This directory contains integrated analysis of IFN levels (IFN1 through IFN10) with dual-axis visualization.

## Key Files

### Main Results
- `ifn1_to_ifn10_integrated_dualaxis_linear_plot.png` - Main figure showing IFN1-10 integrated analysis
- `infection_counts_by_time.csv` - Time-series infection count data

### Analysis Scripts
- `create_ifn10_only_plots.py` - Script to generate IFN10-specific plots
- `mdbk_small_vero_0821.go` - Main simulation code (binary)

### Configuration
- `go.mod` - Go module dependencies
- `go.sum` - Go module checksums

## Complete Dataset

**Note:** Due to file size limitations, detailed simulation runs (~85GB total) are not included in this repository.

The complete dataset includes:
- Multiple kjump_r parameter variations
- Spatial dynamics videos (*.mp4)
- Full simulation outputs for each condition
- Per-timepoint snapshots (7h, 13h, 19h, 25h)

### Excluded Data Structure
```
fig6/detials/
├── kjumpr_0.001_noIFN/
│   ├── *_Dinit*_DIPBst*_JRand_Vinit*_VBst*/
│   │   ├── video.mp4
│   │   ├── simulation_*_hours.png
│   │   ├── simulation_output.csv
│   │   └── selected_frames_combined.png
│   └── ...
└── [other kjump_r conditions]/
```

## To Access Complete Data
The full dataset with all simulation runs, videos, and detailed outputs is available upon request.
Contact the authors or refer to the associated publication.

## File Structure (Repository)
```
fig6/
├── README.md                                              (this file)
├── ifn1_to_ifn10_integrated_dualaxis_linear_plot.png    (main figure)
├── infection_counts_by_time.csv                         (summary data)
├── create_ifn10_only_plots.py                           (analysis script)
├── go.mod, go.sum                                       (dependencies)
└── [detials/ directory excluded from repository]
```

