# Wiang Nong Lom Wetland Hydrological Change Analysis

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19115866.svg)](https://doi.org/10.5281/zenodo.19115866)

Satellite-based monitoring of hydrological and vegetation changes at Wiang Nong Lom wetland, Thailand, following the Kaem Ling water storage project (2022).

## Overview

This repository contains the complete analysis workflow for characterising hydrological regime shifts using Sentinel-1 SAR and Sentinel-2 multispectral imagery (2019–2024). The study applies change-point detection and hydroperiod persistence mapping to evaluate post-intervention ecological trajectories in a data-scarce tropical wetland context.

**Key methods:**
- PELT change-point detection on SAR-derived inundation time series
- Piecewise linear threshold modelling (NDVI vs hydroperiod persistence)
- Landscape metrics (LPI, patch analysis)
- Pre/post intervention comparison using Mann-Whitney U tests

## Repository Structure

```
wiang-nong-lom_draft/
├── 01_gee_export.ipynb              # Google Earth Engine data export
├── 02_data_preprocessing.ipynb      # Build analysis cubes (inundation, NDVI)
├── 03_change_point_detection.ipynb  # PELT breakpoint detection
├── 04_inundation_transition.ipynb   # Hydroperiod persistence analysis
├── 05_vegetation_response.ipynb     # NDVI threshold modelling
├── 06_landscape_metrics.ipynb       # Patch metrics and LPI
├── requirements.txt                 # Python dependencies
└── README.md                        # This file
```

## Workflow

1. **Data Export (01):** Extract Sentinel-1 (VV backscatter) and Sentinel-2 (B4, B8) from GEE
2. **Preprocessing (02):** Build 3D data cubes, apply water threshold (-16 dB), compute NDVI
3. **Change Detection (03):** Identify structural breakpoints using PELT algorithm
4. **Transition Analysis (04):** Compute hydroperiod persistence, McNemar's test, transition matrices
5. **Vegetation Response (05):** Piecewise linear modelling, ecological risk zone mapping
6. **Landscape Metrics (06):** Connected component analysis, LPI calculation

## Data Requirements

**Input data (not included in repository):**
- Sentinel-1 GRD IW VV imagery (monthly median composites, 10m)
- Sentinel-2 L2A surface reflectance (B4, B8, cloud-masked, 10m)
- OSM wetland boundary polygon (`data/raw/wetland_mask.tif`)

**Generated outputs:**
- `data/processed/inundation_cube.nc` — Binary inundation masks (time, y, x)
- `data/processed/ndvi_cube.nc` — NDVI values (time, y, x)
- `data/processed/breakpoint.json` — Detected structural breakpoints
- `outputs/` — Figures and analysis results

## Data Availability

This repository includes:
- **Analysis code** (notebooks 01-06)
- **Processed data cubes** (`data/processed/*.nc`) — ready for analysis
- **Wetland mask** (`data/raw/wetland_mask.tif`)

Raw Sentinel-1/Sentinel-2 imagery is not included due to size. Users can:
1. Use the provided processed cubes directly (notebooks 02-06), OR
2. Run notebook 01 to download raw imagery from Google Earth Engine

## Installation

```bash
# Create conda environment
conda create -n wiang-nong-lom python=3.10
conda activate wiang-nong-lom

# Install dependencies
pip install -r requirements.txt

# Authenticate Google Earth Engine (for notebook 01)
earthengine authenticate
```

## Google Earth Engine Setup

Notebook 01 requires an active GEE project to access and export satellite data.

### 1. Create a GEE Project

1. Go to [Google Earth Engine](https://earthengine.google.com/)
2. Sign in with a Google account
3. Navigate to **Projects** → **New Project**
4. Note your **Project ID** (e.g., `wiang-nom-lom`)

### 2. Update Notebook 01

Open `01_gee_export.ipynb` and update the project initialization:

```python
# Line ~101 in notebook 01
ee.Initialize(project='YOUR_PROJECT_ID')  # Replace with your GEE project ID
```

### 3. Verify Authentication

```bash
python -c "import ee; ee.Initialize(); print('GEE authenticated successfully')"
```

**Note:** First-time users may need to complete the OAuth flow by visiting the provided URL and copying the authorization code.

## Usage

Run notebooks sequentially:

```bash
jupyter notebook 01_gee_export.ipynb
# ... continue through 06_landscape_metrics.ipynb
```

**Note:** Notebook 01 requires Google Earth Engine access and authentication. Subsequent notebooks use exported local files.

## Key Results

- **Primary breakpoint:** June 2021 (inundation extent series)
- **Vegetation threshold:** p* = 58% persistence (above which NDVI decline accelerates)
- **Critical zones:** 17.6% of wetland area classified as high vegetation stress
- **Landscape pattern:** LPI declined despite 34% wet season water area increase (infilling pattern)

## Citation

If using this code or data, please cite:

```
[dataset] Limsupaputtikun, K., Pailoplee, S. & Aruninta, A. (2026). 
Wiang Nong Lom wetland hydrological change analysis. Zenodo. 
https://doi.org/10.5281/zenodo.19115866
```

## License

MIT License — See LICENSE file for details.

## Acknowledgements

- Satellite data: European Space Agency Copernicus Programme (Sentinel-1, Sentinel-2)
- Processing platform: Google Earth Engine
- Change-point detection: Python `ruptures` library (Truong et al., 2020)
