"""
Utility functions for CoastCams analysis.

This module provides helper functions for smoothing, filtering,
and mathematical operations used throughout the CoastCams toolbox.
"""

import numpy as np
from scipy import ndimage, signal
from typing import Tuple, Optional


def smooth2d(data: np.ndarray, window_size: int = 3) -> np.ndarray:
    """
    Apply 2D smoothing to data using a moving average filter.

    Parameters
    ----------
    data : np.ndarray
        2D array to smooth
    window_size : int, optional
        Size of the smoothing window (default: 3)

    Returns
    -------
    np.ndarray
        Smoothed 2D array
    """
    if window_size < 1:
        return data

    # Create uniform kernel
    kernel = np.ones((window_size, window_size)) / (window_size ** 2)

    # Apply convolution with same mode to preserve size
    smoothed = signal.convolve2d(data, kernel, mode='same', boundary='symm')

    return smoothed


def local_maxima(data: np.ndarray, min_distance: int = 1) -> np.ndarray:
    """
    Find local maxima in 1D array.

    Parameters
    ----------
    data : np.ndarray
        1D array to find maxima in
    min_distance : int, optional
        Minimum distance between maxima (default: 1)

    Returns
    -------
    np.ndarray
        Indices of local maxima
    """
    from scipy.signal import find_peaks

    peaks, _ = find_peaks(data, distance=min_distance)
    return peaks


def filter_mean(data: np.ndarray, axis: int = 0) -> np.ndarray:
    """
    Remove mean along specified axis.

    Parameters
    ----------
    data : np.ndarray
        Input array
    axis : int, optional
        Axis along which to compute mean (default: 0)

    Returns
    -------
    np.ndarray
        Data with mean removed along specified axis
    """
    return data - np.mean(data, axis=axis, keepdims=True)


def linear_wave_celerity(period: float, depth: float, g: float = 9.81) -> float:
    """
    Calculate wave celerity using linear wave theory (dispersion relation).

    Parameters
    ----------
    period : float
        Wave period in seconds
    depth : float
        Water depth in meters
    g : float, optional
        Gravitational acceleration (default: 9.81 m/s^2)

    Returns
    -------
    float
        Wave celerity in m/s
    """
    if period <= 0 or depth <= 0:
        return 0.0

    omega = 2 * np.pi / period  # Angular frequency

    # Solve dispersion relation iteratively
    # omega^2 = g*k*tanh(k*h) where k is wavenumber
    k = omega / np.sqrt(g * depth)  # Initial guess (shallow water)

    # Newton-Raphson iteration
    for _ in range(10):
        k_new = omega**2 / (g * np.tanh(k * depth))
        if abs(k_new - k) < 1e-6:
            break
        k = k_new

    # Celerity C = omega / k
    celerity = omega / k

    return celerity


def calculate_depth_from_celerity(celerity: float, period: float,
                                  g: float = 9.81) -> float:
    """
    Calculate water depth from wave celerity and period using linear wave theory.

    Parameters
    ----------
    celerity : float
        Wave celerity in m/s
    period : float
        Wave period in seconds
    g : float, optional
        Gravitational acceleration (default: 9.81 m/s^2)

    Returns
    -------
    float
        Estimated water depth in meters
    """
    if celerity <= 0 or period <= 0:
        return 0.0

    omega = 2 * np.pi / period
    k = omega / celerity  # Wavenumber

    # Solve for depth from dispersion relation
    # omega^2 = g*k*tanh(k*h)
    # Iteratively solve for h
    h = celerity * period / 2  # Initial guess

    for _ in range(20):
        tanh_kh = (omega**2) / (g * k)
        if tanh_kh >= 1:
            tanh_kh = 0.9999
        kh = np.arctanh(tanh_kh)
        h_new = kh / k

        if abs(h_new - h) < 1e-6:
            break
        h = h_new

    return max(h, 0.0)


def moving_average(data: np.ndarray, window: int) -> np.ndarray:
    """
    Compute moving average of 1D array.

    Parameters
    ----------
    data : np.ndarray
        Input 1D array
    window : int
        Window size for moving average

    Returns
    -------
    np.ndarray
        Smoothed array
    """
    if window < 1:
        return data

    cumsum = np.cumsum(np.insert(data, 0, 0))
    return (cumsum[window:] - cumsum[:-window]) / window


def detrend_2d(data: np.ndarray, axis: int = -1) -> np.ndarray:
    """
    Remove linear trend from 2D data along specified axis.

    Parameters
    ----------
    data : np.ndarray
        2D input array
    axis : int, optional
        Axis along which to detrend (default: -1, last axis)

    Returns
    -------
    np.ndarray
        Detrended data
    """
    from scipy.signal import detrend
    return detrend(data, axis=axis, type='linear')


def butter_lowpass_filter(data: np.ndarray, cutoff: float, fs: float,
                          order: int = 5) -> np.ndarray:
    """
    Apply Butterworth low-pass filter.

    Parameters
    ----------
    data : np.ndarray
        Input signal
    cutoff : float
        Cutoff frequency in Hz
    fs : float
        Sampling frequency in Hz
    order : int, optional
        Filter order (default: 5)

    Returns
    -------
    np.ndarray
        Filtered signal
    """
    from scipy.signal import butter, filtfilt

    nyquist = 0.5 * fs
    normal_cutoff = cutoff / nyquist
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    filtered = filtfilt(b, a, data)

    return filtered


def butter_highpass_filter(data: np.ndarray, cutoff: float, fs: float,
                           order: int = 5) -> np.ndarray:
    """
    Apply Butterworth high-pass filter.

    Parameters
    ----------
    data : np.ndarray
        Input signal
    cutoff : float
        Cutoff frequency in Hz
    fs : float
        Sampling frequency in Hz
    order : int, optional
        Filter order (default: 5)

    Returns
    -------
    np.ndarray
        Filtered signal
    """
    from scipy.signal import butter, filtfilt

    nyquist = 0.5 * fs
    normal_cutoff = cutoff / nyquist
    b, a = butter(order, normal_cutoff, btype='high', analog=False)
    filtered = filtfilt(b, a, data)

    return filtered


def resample_timeseries(timestamps: np.ndarray, data: np.ndarray,
                       interval_minutes: int = 15) -> Tuple[np.ndarray, np.ndarray]:
    """
    Resample time series data to regular intervals.

    Parameters
    ----------
    timestamps : np.ndarray
        Array of datetime objects
    data : np.ndarray
        Data values corresponding to timestamps
    interval_minutes : int, optional
        Resampling interval in minutes (default: 15)

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        Resampled timestamps and data
    """
    import pandas as pd

    # Create pandas Series
    series = pd.Series(data, index=timestamps)

    # Resample to regular intervals and take mean
    resampled = series.resample(f'{interval_minutes}min').mean()

    return resampled.index.to_numpy(), resampled.values


def extract_timestamp_from_filename(filename: str,
                                    pattern: str = "S_%d_%Y%m%d%H%M") -> Optional[str]:
    """
    Extract timestamp from filename based on pattern.

    Parameters
    ----------
    filename : str
        Filename to parse
    pattern : str, optional
        Expected filename pattern (default: "S_%d_%Y%m%d%H%M")

    Returns
    -------
    Optional[str]
        Extracted timestamp string or None if parsing fails
    """
    import re
    from datetime import datetime

    # Extract date/time portion from filename
    # Expected format: S_1_YYYYMMDDHHSS.jpeg
    match = re.search(r'(\d{12})', filename)
    if match:
        return match.group(1)

    return None


def compute_wave_energy(wave_height: float, rho: float = 1025,
                       g: float = 9.81) -> float:
    """
    Compute wave energy per unit area.

    Parameters
    ----------
    wave_height : float
        Significant wave height in meters
    rho : float, optional
        Water density in kg/m^3 (default: 1025 for seawater)
    g : float, optional
        Gravitational acceleration (default: 9.81 m/s^2)

    Returns
    -------
    float
        Wave energy in J/m^2
    """
    # E = (1/8) * rho * g * H^2
    return (1/8) * rho * g * (wave_height ** 2)


def nan_helper(y: np.ndarray) -> Tuple[np.ndarray, callable]:
    """
    Helper function to handle NaN values in interpolation.

    Parameters
    ----------
    y : np.ndarray
        Array possibly containing NaN values

    Returns
    -------
    Tuple[np.ndarray, callable]
        Boolean array indicating NaNs and a function to get indices
    """
    return np.isnan(y), lambda z: z.nonzero()[0]


def interpolate_nans(data: np.ndarray) -> np.ndarray:
    """
    Interpolate NaN values in 1D array.

    Parameters
    ----------
    data : np.ndarray
        1D array with possible NaN values

    Returns
    -------
    np.ndarray
        Array with NaNs interpolated
    """
    nans, x = nan_helper(data)
    if not np.any(nans):
        return data

    data_copy = data.copy()
    data_copy[nans] = np.interp(x(nans), x(~nans), data[~nans])

    return data_copy
