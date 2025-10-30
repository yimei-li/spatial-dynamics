# Figure 5: IFN Response Analysis

## Overview
This directory contains analysis of interferon (IFN) response in relation to DIP burst size variations.

## Key Files

### Summary Data
- `IFNclr3_30runs_summary.csv` - Main summary of 30 simulation runs
- `ifn_peak_vs_dipburst_baseline.csv` - Peak IFN levels vs DIP burst size baseline data
- `infection_counts_by_time.csv` - Time-series infection data

### Analysis Scripts
- `1_bst_IFN_time_0503.py` - IFN time-series analysis
- `4_process_peak_ifn.py` - Peak IFN processing
- `5_analyze_peak_ifn.py` - Peak IFN analysis
- `loop_burstSizeD_IFNclr3_30runs.sh` - Batch simulation script
- `loop_burstSizeD_0522.sh` - Parameter sweep script

### Key Results
- `5_peak_ifn_vs_dip_burst_size.png` / `.pdf` - Main figure: Peak IFN vs DIP burst size
- `5_peak_ifn_vs_dip_burst_size_ratio.pdf` - Ratio analysis
- `6_reproduced_avg_ifn_plot.png` - Average IFN reproduction
- `7_combined_ifn_analysis.png` - Combined analysis visualization

### Detailed Results
- `IFNclr3_30runs_global_celltocell_tau95_option1/` - Contains summary CSV files for different DIP burst sizes
  - `summary_DIPBst*.csv` - Summary data for specific burst sizes (650-1600)
  - `first_script_results.png` / `.pdf` - Initial analysis results
  - Analysis scripts: `plot_*.py`

## Complete Dataset

**Note:** Due to file size limitations, individual simulation runs (~585GB total) are not included in this repository.

The complete dataset includes:
- 30 runs × multiple DIP burst size conditions
- Video files (*.mp4) showing spatial dynamics
- Full simulation output files (*.go binary files)
- Detailed per-run CSV files

### To Access Complete Data
Contact the authors or refer to the associated publication for data availability.

## File Structure
```
fig5/
├── README.md                          (this file)
├── *_summary.csv                      (aggregated results)
├── *.py                               (analysis scripts)
├── *.sh                               (simulation batch scripts)
├── *.png, *.pdf                       (key figures)
└── IFNclr3_30runs_global_celltocell_tau95_option1/
    └── summary_DIPBst*.csv           (per-condition summaries)
```

## Excluded Files (Available Upon Request)
- Individual simulation runs: `IFNclr3_30runs_global_celltocell_tau95_option1_*/`
- Video animations: `*.mp4`
- Binary simulation files: `*.go`

