# KIC Spheroid Calcium Pipeline (Script 01/03)

This repository contains an R script for preprocessing calcium transient measurements generated using the **CyteSeer Cytometer (Vala Sciences)** and analyzed with the onboard **Cardiac_Two-Channel_TimeSeries_Extended** algorithm.

This script represents the first step (**01/03**) of a modular **KIC Spheroid Calcium analysis pipeline**.

The pipeline processes raw CyteSeer CSV exports from **3D cardiac spheroids**, reconstructs calcium waveforms, computes derived kinetic parameters, and exports structured datasets ready for downstream statistical analysis.

---

## What the pipeline does

Starting from raw **CyteSeerResults** export folders, the script:

- Recursively scans and loads multiple CyteSeerResults folders  
- Reads:
  - `Primary_DataTable_CStats.csv` (cell-level traces)
  - `_Cell Peak Measurements_DataTable.csv`
  - `_Whole Image Measurements_DataTable_CStats.csv`
  - `_Whole-Image Peak Measurements_DataTable.csv`
- Parses KIC waveform strings to reconstruct calcium time series
- Computes additional calcium kinetic parameters:
  - **Upstroke 10→25, 10→50, 10→75**
  - **Downstroke 90→75, 90→50, 90→25**
- Uses waveform-based estimation when available
- Applies fallback interpolation when waveform reconstruction is not possible
- Exports structured outputs for both **cell-level** and **whole-spheroid** analysis

All file paths are selected interactively via GUI dialogs (`tk_choose.dir`).

---

## Required inputs

### CyteSeer export folders

- One or multiple **CyteSeerResults*** directories  
- Files may be stored in nested subfolders (recursive search supported)

Each dataset must contain:

- `Primary_DataTable_CStats.csv`  
- `_Cell Peak Measurements_DataTable.csv`  
- `_Whole Image Measurements_DataTable_CStats.csv`  
- `_Whole-Image Peak Measurements_DataTable.csv`  

All files must be generated using the **Cardiac_Two-Channel_TimeSeries_Extended** algorithm within CyteSeer.

---

## Cleaned summary output

The script generates a structured output folder containing:

- `raw_traces`
- `raw_csv`
- `cell_peak_parameters`
- `whole_image_traces`
- `whole_image_csv`
- `whole_image_parameters`
- `summary`

The exported datasets contain:

- Original measurement values  
- Derived kinetic parameters  
- Per-peak data  
- Per-cell mean summaries  
- Whole-spheroid summary statistics (mean, SD, median, IQR, Num.Peaks)

The output is ready for:

- Statistical analysis  
- Spheroid-level aggregation  
- Condition-based comparisons  
- Downstream merging (**Script 02**)  

---

## Typical use cases

- Spheroid calcium-handling validation  
- hiPSC-derived disease modeling  
- Drug response studies (baseline vs treatment)  
- Patient-specific calcium kinetics comparison  
- Standardized preprocessing in collaborative cardiovascular projects  

---

## Position in the KIC Pipeline

This script is **Script 01** of a 3-step workflow:

- **Script 01 (this repository)** – Raw CyteSeer CSV processing and structured dataset generation  
- Script 02a – Single-batch MEANS merging and filtering
- Script 02b – Cross-batch merging
- Script 03 – Statistical analysis and visualization

---

## Methods Description

Calcium transient measurements from 3D cardiac spheroids were acquired using the **CyteSeer Cytometer (Vala Sciences)** and analyzed with the **Cardiac_Two-Channel_TimeSeries_Extended** algorithm. Raw CSV exports were processed using a custom R-based GUI pipeline. Waveform strings were parsed to reconstruct time series data, and additional upstroke and downstroke timing metrics were computed. The resulting structured dataset was exported for downstream statistical analysis and visualization.

---

## Authorship

This script was developed by **Michele Buono, Talitha Spanjersberg, Nikki Scheen, Nina van der Wilt** and can be used freely for research purposes, provided appropriate citation of the authors.

**Buono, M. F., Spanjersberg, T., Nikki, S., & van der Wilt, C. N. (2026). KIC Spheroid Calcium Pipeline (Script 01/03) – CyteSeer Preprocessing & Waveform Reconstruction (v.1.0). Zenodo. https://doi.org/10.5281/zenodo.18798891**

The overall workflow, structure, and clarity of the pipeline were iteratively refined with assistance from **OpenAI – ChatGPT 5.2**, which was used as a tool to improve code organization, documentation, and usability.
