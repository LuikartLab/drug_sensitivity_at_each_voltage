# drug_sensitivity_at_each_voltage

Python tool to quantify and visualize voltage-gated channel **drug (blocker) sensitivity** from whole-cell current–voltage (I–V) recordings, comparing control (WT) and knockout (KO) cells.

## Overview

Given a directory of CSV files — one per drug/blocker — each containing I–V curve data from WT and KO cells recorded **with and without** the drug, the script computes, **at each voltage**, the percent of current blocked by the drug:

```
block% = 100 × (1 − mean_drug / mean_control)
```

for WT and KO separately (with the standard error propagated through the ratio), and renders a **heatmap** of drug sensitivity across the chosen voltage range (WT vs KO).

## Input format

- A folder containing one or more `.csv` files (one file per drug/blocker; the file name is used as the drug label).
- Each CSV is an **X–Y data table exported from GraphPad Prism**, where **X = voltage (mV)** and **Y = current amplitude** (e.g., pA) — i.e., one current–voltage (I–V) curve per cell/condition.
- The X column is read as `Voltage (mV)`. The Y (current) columns are labeled by genotype (`wt`, `dKO`) and by condition (control vs drug); the script matches them automatically by these labels.

## Outputs (written into the input folder)

- `Drug sensitivity at each voltage.csv` — per-voltage means, SEMs, and block% for WT and KO.
- A heatmap **TIFF** (requires `tifffile`) and a labeled **PNG** with colorbar (requires `matplotlib`), using a custom "fire/ice" colormap.

## Requirements

- Python 3
- `numpy`, `pandas` (required)
- `tifffile` (optional — TIFF heatmap), `matplotlib` (optional — labeled PNG)

```
pip install numpy pandas tifffile matplotlib
```

## Usage

```
python DrugSensitivityAtEachVoltage.py
```

You will be prompted for the directory of CSV files and the minimum/maximum voltage (mV) to analyze (press Enter for the full range).

## Notes

Development of this code was assisted by an AI tool under author direction and validation; see `Ai Assistance.md`.

## Citation / availability

Luikart Lab. *drug_sensitivity_at_each_voltage*. GitHub: https://github.com/LuikartLab/drug_sensitivity_at_each_voltage (archived on Zenodo upon release; DOI to be added).

## License

Released under the MIT License (see `LICENSE`).
