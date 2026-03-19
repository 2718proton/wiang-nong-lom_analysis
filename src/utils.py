"""
Utility functions for Wiang Nong Lom hydro-ecological analysis.

This module provides shared functions used across multiple analysis notebooks.
"""

import numpy as np
import pandas as pd
import xarray as xr
import rasterio
from pathlib import Path
from datetime import datetime


# =============================================================================
# Configuration
# =============================================================================

# Study area bounding box [min_lon, min_lat, max_lon, max_lat]
BBOX = [99.98, 20.18, 100.05, 20.25]

# Thailand seasons
WET_MONTHS = [6, 7, 8, 9, 10]  # June - October
DRY_MONTHS = [12, 1, 2, 3, 4]  # December - April

# Default paths
DATA_DIR = Path('data')
RAW_DIR = DATA_DIR / 'raw'
PROCESSED_DIR = DATA_DIR / 'processed'
OUTPUT_DIR = Path('outputs')


# =============================================================================
# Data Loading
# =============================================================================

def load_raster(filepath):
    """
    Load a GeoTIFF file and return data with coordinates.

    Parameters
    ----------
    filepath : str or Path
        Path to GeoTIFF file

    Returns
    -------
    tuple : (data, x_coords, y_coords, crs, transform)
    """
    with rasterio.open(filepath) as src:
        data = src.read(1)
        transform = src.transform
        crs = src.crs
        height, width = data.shape

        # Get 1D coordinate arrays directly from transform
        x_coords = np.array([transform[2] + (col + 0.5) * transform[0] for col in range(width)])
        y_coords = np.array([transform[5] + (row + 0.5) * transform[4] for row in range(height)])

    return data, x_coords, y_coords, crs, transform


def load_inundation_cube():
    """Load the processed inundation cube."""
    return xr.open_dataarray(PROCESSED_DIR / 'inundation_cube.nc')


def load_ndvi_cube():
    """Load the processed NDVI cube."""
    return xr.open_dataarray(PROCESSED_DIR / 'ndvi_cube.nc')


def load_breakpoint_config():
    """Load breakpoint configuration."""
    config = pd.read_csv(PROCESSED_DIR / 'breakpoint_config.csv')
    return {
        'date': pd.Timestamp(config['breakpoint_date'].values[0]),
        'index': int(config['breakpoint_index'].values[0])
    }


# =============================================================================
# Date Parsing
# =============================================================================

def parse_date_from_filename(filename, prefix='S1_VV_'):
    """
    Extract date from filename like S1_VV_202006.tif or S2_202006.tif

    Parameters
    ----------
    filename : str or Path
        Filename to parse
    prefix : str
        Prefix before the date string

    Returns
    -------
    datetime
    """
    name = Path(filename).stem
    date_str = name.split('_')[-1]
    year = int(date_str[:4])
    month = int(date_str[4:6])
    return datetime(year, month, 1)


# =============================================================================
# Hydroperiod Metrics
# =============================================================================

def compute_persistence(inundation_cube):
    """
    Compute hydroperiod persistence (% of months inundated) per pixel.

    Parameters
    ----------
    inundation_cube : xr.DataArray
        3D cube with dims (time, y, x), values: 1=water, 0=non-water, -1=nodata

    Returns
    -------
    xr.DataArray : Persistence percentage per pixel
    """
    valid_mask = inundation_cube >= 0
    water_mask = inundation_cube == 1

    n_valid = valid_mask.sum(dim='time').astype(float)
    n_valid = n_valid.where(n_valid > 0)

    n_water = water_mask.where(valid_mask).sum(dim='time')
    persistence = (n_water / n_valid * 100)
    persistence.name = 'persistence'

    return persistence


def classify_hydroperiod(persistence, thresholds=(20, 60, 90)):
    """
    Classify pixels into hydroperiod categories based on persistence.

    Categories:
    - 0: Rarely inundated (0 to threshold[0])
    - 1: Seasonally inundated (threshold[0] to threshold[1])
    - 2: Frequently inundated (threshold[1] to threshold[2])
    - 3: Persistently inundated (threshold[2] to 100)
    - -1: NoData

    Parameters
    ----------
    persistence : xr.DataArray or np.ndarray
        Persistence values (0-100)
    thresholds : tuple
        Three thresholds defining category boundaries

    Returns
    -------
    np.ndarray : Classification array
    """
    if isinstance(persistence, xr.DataArray):
        p = persistence.values
    else:
        p = persistence

    classes = np.full_like(p, -1, dtype=np.int8)
    valid = ~np.isnan(p)

    classes[valid & (p < thresholds[0])] = 0
    classes[valid & (p >= thresholds[0]) & (p < thresholds[1])] = 1
    classes[valid & (p >= thresholds[1]) & (p < thresholds[2])] = 2
    classes[valid & (p >= thresholds[2])] = 3

    return classes


# =============================================================================
# Seasonal Analysis
# =============================================================================

def get_season_mask(times, months):
    """
    Create boolean mask for specified months.

    Parameters
    ----------
    times : array-like
        Array of datetime values
    months : list
        List of month numbers (1-12)

    Returns
    -------
    np.ndarray : Boolean mask
    """
    time_months = pd.to_datetime(times).month
    return np.isin(time_months, months)


def compute_seasonal_composite(cube, season_months):
    """
    Compute seasonal inundation probability map.

    Parameters
    ----------
    cube : xr.DataArray
        Inundation cube with time dimension
    season_months : list
        List of months defining the season

    Returns
    -------
    xr.DataArray : Probability of inundation (0-1) per pixel
    """
    season_mask = get_season_mask(cube.time.values, season_months)
    season_cube = cube.isel(time=season_mask)

    if len(season_cube.time) == 0:
        return None

    valid_mask = season_cube >= 0
    water_mask = season_cube == 1

    n_valid = valid_mask.sum(dim='time').astype(float)
    n_valid = n_valid.where(n_valid > 0)

    n_water = water_mask.where(valid_mask).sum(dim='time')
    probability = n_water / n_valid

    return probability


# =============================================================================
# Statistical Functions
# =============================================================================

def deseasonalize(series, period=12):
    """
    Remove seasonal component from time series using monthly averages.

    Parameters
    ----------
    series : pd.Series
        Time series with DatetimeIndex
    period : int
        Seasonality period (12 for monthly data)

    Returns
    -------
    pd.Series : Deseasonalized (anomaly) series
    """
    monthly_means = series.groupby(series.index.month).mean()
    seasonal = series.index.month.map(monthly_means)
    return series - seasonal


def winsorize(data, percentiles=(1, 99)):
    """
    Winsorize data by clipping to specified percentiles.

    Parameters
    ----------
    data : np.ndarray
        Data to winsorize
    percentiles : tuple
        Lower and upper percentiles for clipping

    Returns
    -------
    np.ndarray : Winsorized data
    """
    lower = np.nanpercentile(data, percentiles[0])
    upper = np.nanpercentile(data, percentiles[1])
    return np.clip(data, lower, upper)


def compute_binned_response(x, y, n_bins=20, ci=95, min_samples=10):
    """
    Compute binned response of y to x with confidence intervals.

    Parameters
    ----------
    x : np.ndarray
        Independent variable (e.g., hydroperiod persistence)
    y : np.ndarray
        Dependent variable (e.g., NDVI)
    n_bins : int
        Number of bins
    ci : float
        Confidence interval percentage
    min_samples : int
        Minimum samples per bin

    Returns
    -------
    pd.DataFrame : Binned statistics
    """
    from scipy import stats

    bins = np.linspace(np.nanmin(x), np.nanmax(x), n_bins + 1)
    bin_centers = (bins[:-1] + bins[1:]) / 2

    results = []
    z = stats.norm.ppf(1 - (100 - ci) / 200)

    for i in range(n_bins):
        mask = (x >= bins[i]) & (x < bins[i+1])
        bin_y = y[mask]
        bin_y = bin_y[~np.isnan(bin_y)]

        if len(bin_y) >= min_samples:
            mean_val = np.mean(bin_y)
            std_val = np.std(bin_y)
            n = len(bin_y)
            margin = z * std_val / np.sqrt(n)

            results.append({
                'bin_center': bin_centers[i],
                'mean': mean_val,
                'std': std_val,
                'ci_lower': mean_val - margin,
                'ci_upper': mean_val + margin,
                'n': n
            })

    return pd.DataFrame(results)
