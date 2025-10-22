"""
Wave parameter analysis for CoastCams.

This module extracts wave characteristics including height, period,
celerity, and energy from timestack images.
"""

import numpy as np
from scipy import signal, interpolate
from typing import Dict, Tuple, Optional, List
from .utils import (local_maxima, compute_wave_energy, linear_wave_celerity,
                   calculate_depth_from_celerity)


class WaveAnalyzer:
    """
    Analyze wave parameters from timestack data.

    Extracts wave height, period, celerity, energy, and other
    characteristics from preprocessed timestack images.
    """

    def __init__(self, config=None):
        """
        Initialize wave analyzer.

        Parameters
        ----------
        config : CoastCamsConfig, optional
            Configuration object with analysis parameters
        """
        if config is not None:
            self.dt = 1.0 / config.acquisition_frequency
            self.pixel_resolution = config.pixel_resolution
            self.min_period = config.wave_period_min
            self.max_period = config.wave_period_max
            self.gravity = config.gravity
        else:
            # Default parameters
            self.dt = 0.5  # Time step (seconds)
            self.pixel_resolution = 0.1  # meters/pixel
            self.min_period = 4.0  # seconds
            self.max_period = 25.0  # seconds
            self.gravity = 9.81  # m/s^2

    def analyze_timestack(self, timestack: np.ndarray,
                         cross_shore_positions: np.ndarray,
                         shorelines: list = None) -> Dict:
        """
        Perform complete wave analysis on timestack data.

        Parameters
        ----------
        timestack : np.ndarray
            2D timestack array (space x time)
        cross_shore_positions : np.ndarray
            Cross-shore pixel positions
        shorelines : list, optional
            List of shoreline arrays for each time step (for improved wave height estimation)

        Returns
        -------
        Dict
            Dictionary containing wave parameters
        """
        results = {}

        # Analyze timestack intensity to compute wave height
        # MATLAB approach: analyze intensity variations at breaking location
        # We'll use the median cross-shore position and analyze its time series

        # Find the location with maximum variability (breaking zone)
        if timestack.shape[0] > 0 and timestack.shape[1] > 1:
            # Compute standard deviation along time axis for each position
            std_profile = np.std(timestack, axis=1)

            # Find location of maximum variability (likely breaking zone)
            if len(std_profile) > 0 and not np.all(np.isnan(std_profile)):
                break_location = np.nanargmax(std_profile)

                # Extract time series at breaking location
                intensity_timeseries = timestack[break_location, :]

                # Detrend the signal
                x = np.arange(len(intensity_timeseries))
                p = np.polyfit(x, intensity_timeseries, 1)
                detrended = intensity_timeseries - (p[0] * x + p[1])
                detrended = detrended - np.mean(detrended)

                # Compute wave height using spectral method: Hs = 4 * σ
                # Scale by pixel resolution to convert intensity variations to physical units
                # Empirical scaling factor for intensity-based measurements
                sigma = np.std(detrended)

                # Convert normalized intensity (0-1 range) to meters
                # Typical wave amplitude corresponds to ~10-50 intensity units in normalized timestack
                # Scale factor: assume intensity range of 0.1 corresponds to ~1m wave height
                scale_factor = 10.0  # Empirical factor
                Hs = 4.0 * sigma * scale_factor * self.pixel_resolution

                # Clip to realistic range (0.1m to 10m)
                Hs = np.clip(Hs, 0.1, 10.0)

                # Mean wave height
                Hm = 0.63 * Hs

                # Debug output
                print(f"  Timestack analysis at position {break_location} ({cross_shore_positions[break_location]:.1f}m)")
                print(f"  Intensity std: {sigma:.4f}")
                print(f"  Computed Hs: {Hs:.3f} m, Hm: {Hm:.3f} m")

                results['mean_Hs'] = Hs
                results['mean_Hm'] = Hm
                results['intensity_std'] = sigma
                results['break_location'] = break_location
            else:
                results['mean_Hs'] = np.nan
                results['mean_Hm'] = np.nan
        else:
            # Fallback: analyze timestack intensity variations
            # (less accurate but works without shoreline data)
            wave_heights = []
            wave_periods = []
            wave_energies = []

            for i in range(timestack.shape[0]):
                timeseries = timestack[i, :]

                # Analyze this position
                params = self._analyze_timeseries(timeseries)

                # Scale intensity-based heights to reasonable values
                # Typical intensity range is 0-1, scale by pixel resolution and empirical factor
                scaled_height = params['significant_height'] * self.pixel_resolution * 10.0

                wave_heights.append(scaled_height)
                wave_periods.append(params['mean_period'])
                wave_energies.append(params['energy'])

            # Store arrays
            wave_heights_array = np.array(wave_heights)
            wave_periods_array = np.array(wave_periods)
            wave_energies_array = np.array(wave_energies)

            results['wave_heights'] = wave_heights_array
            results['wave_periods'] = wave_periods_array
            results['wave_energies'] = wave_energies_array

            # Compute aggregate statistics
            valid_heights = wave_heights_array[~np.isnan(wave_heights_array)]

            results['mean_Hs'] = np.mean(valid_heights) if len(valid_heights) > 0 else np.nan
            results['max_Hs'] = np.max(valid_heights) if len(valid_heights) > 0 else np.nan

        # Analyze wave periods from timestack
        if timestack.shape[0] > 0 and timestack.shape[1] > 1 and 'break_location' in results:
            # Use intensity time series at breaking location for period calculation
            break_loc = results['break_location']
            intensity_ts = timestack[break_loc, :]

            # Detrend
            x = np.arange(len(intensity_ts))
            p = np.polyfit(x, intensity_ts, 1)
            detrended_ts = intensity_ts - (p[0] * x + p[1])
            detrended_ts = detrended_ts - np.mean(detrended_ts)

            # Compute periods
            periods = self._compute_wave_periods_from_timeseries(detrended_ts)
            if len(periods) > 0:
                mean_period = np.mean(periods)
                results['mean_Tm'] = mean_period
                print(f"  Mean wave period: {mean_period:.2f} s ({len(periods)} waves detected)")

                # Create array of periods for bathymetry calculation
                # Use the mean period for all positions
                results['wave_periods'] = np.full(len(cross_shore_positions), mean_period)
            else:
                results['mean_Tm'] = np.nan
                results['wave_periods'] = np.full(len(cross_shore_positions), np.nan)
                print(f"  Mean wave period: No valid periods detected")
        else:
            # Fallback: analyze timestack for periods
            wave_periods = []
            for i in range(min(10, timestack.shape[0])):  # Sample first 10 positions
                timeseries = timestack[i, :]
                params = self._analyze_timeseries(timeseries)
                if not np.isnan(params['mean_period']):
                    wave_periods.append(params['mean_period'])

            if len(wave_periods) > 0:
                mean_period = np.mean(wave_periods)
                results['mean_Tm'] = mean_period
                results['wave_periods'] = np.full(len(cross_shore_positions), mean_period)
            else:
                results['mean_Tm'] = np.nan
                results['wave_periods'] = np.full(len(cross_shore_positions), np.nan)

        results['cross_shore_positions'] = cross_shore_positions

        return results

    def _analyze_timeseries(self, timeseries: np.ndarray) -> Dict:
        """
        Analyze a single time series for wave parameters.

        Parameters
        ----------
        timeseries : np.ndarray
            1D time series data

        Returns
        -------
        Dict
            Wave parameters for this time series
        """
        params = {}

        # Remove mean
        ts = timeseries - np.mean(timeseries)

        # Compute wave height using zero-crossing method
        heights = self._compute_wave_heights(ts)
        params['significant_height'] = np.mean(sorted(heights, reverse=True)[:max(1, len(heights)//3)])

        # Compute wave period
        periods = self._compute_wave_periods(ts)
        params['mean_period'] = np.mean(periods) if len(periods) > 0 else np.nan
        params['peak_period'] = self._compute_peak_period(ts)

        # Compute wave energy
        params['energy'] = compute_wave_energy(params['significant_height'])

        return params

    def _compute_wave_heights(self, timeseries: np.ndarray) -> List[float]:
        """
        Compute individual wave heights using zero-crossing method.

        Parameters
        ----------
        timeseries : np.ndarray
            1D time series (mean should be removed)

        Returns
        -------
        List[float]
            List of individual wave heights
        """
        # Find zero crossings
        zero_crossings = np.where(np.diff(np.sign(timeseries)))[0]

        heights = []

        # Measure height between successive crossings
        for i in range(len(zero_crossings) - 1):
            start_idx = zero_crossings[i]
            end_idx = zero_crossings[i + 1]

            segment = timeseries[start_idx:end_idx]

            if len(segment) > 0:
                height = np.max(segment) - np.min(segment)
                heights.append(height)

        return heights

    def _compute_wave_periods_from_timeseries(self, timeseries: np.ndarray) -> List[float]:
        """
        Compute wave periods from a detrended time series using zero-crossing method.
        This matches the MATLAB Wave_Char.m approach.

        Parameters
        ----------
        timeseries : np.ndarray
            Detrended time series

        Returns
        -------
        List[float]
            List of wave periods in seconds
        """
        # Find zero crossings (upward crossings)
        signs = np.sign(timeseries)
        zero_crossings = np.where((signs[:-1] <= 0) & (signs[1:] > 0))[0]

        periods = []

        # Measure time between successive crossings (full wave period)
        for i in range(len(zero_crossings) - 1):
            period = (zero_crossings[i + 1] - zero_crossings[i]) * self.dt

            # Filter by expected period range
            if self.min_period <= period <= self.max_period:
                periods.append(period)

        return periods

    def _compute_wave_periods(self, timeseries: np.ndarray) -> List[float]:
        """
        Compute individual wave periods using zero-crossing method.

        Parameters
        ----------
        timeseries : np.ndarray
            1D time series

        Returns
        -------
        List[float]
            List of wave periods in seconds
        """
        # Find zero crossings (upward crossings)
        signs = np.sign(timeseries)
        zero_crossings = np.where((signs[:-1] <= 0) & (signs[1:] > 0))[0]

        periods = []

        # Measure time between successive crossings
        for i in range(len(zero_crossings) - 1):
            period = (zero_crossings[i + 1] - zero_crossings[i]) * self.dt

            # Filter by expected period range
            if self.min_period <= period <= self.max_period:
                periods.append(period)

        return periods

    def _compute_peak_period(self, timeseries: np.ndarray) -> float:
        """
        Compute peak wave period using spectral analysis.

        Parameters
        ----------
        timeseries : np.ndarray
            1D time series

        Returns
        -------
        float
            Peak period in seconds
        """
        # Compute power spectral density
        frequencies, psd = signal.welch(timeseries, fs=1/self.dt,
                                       nperseg=min(256, len(timeseries)))

        # Filter frequency range corresponding to wave periods
        min_freq = 1.0 / self.max_period
        max_freq = 1.0 / self.min_period

        valid_idx = (frequencies >= min_freq) & (frequencies <= max_freq)
        valid_freqs = frequencies[valid_idx]
        valid_psd = psd[valid_idx]

        if len(valid_psd) == 0:
            return np.nan

        # Find peak frequency
        peak_idx = np.argmax(valid_psd)
        peak_freq = valid_freqs[peak_idx]

        # Convert to period
        peak_period = 1.0 / peak_freq if peak_freq > 0 else np.nan

        return peak_period

    def compute_wave_celerity(self, timestack: np.ndarray,
                             dx: float) -> np.ndarray:
        """
        Compute wave celerity (speed) from timestack using phase tracking.

        Parameters
        ----------
        timestack : np.ndarray
            2D timestack array (space x time)
        dx : float
            Spatial spacing (meters)

        Returns
        -------
        np.ndarray
            Wave celerity values (m/s) for each time point
        """
        # Use cross-correlation to find phase velocity
        celerities = []

        for t in range(1, timestack.shape[1]):
            # Compare successive time steps
            profile1 = timestack[:, t-1]
            profile2 = timestack[:, t]

            # Compute cross-correlation
            correlation = np.correlate(profile1, profile2, mode='full')
            lags = np.arange(-len(profile1) + 1, len(profile1))

            # Find peak correlation
            peak_idx = np.argmax(correlation)
            lag = lags[peak_idx]

            # Compute celerity
            if lag != 0:
                distance = lag * dx
                celerity = distance / self.dt
                celerities.append(celerity)
            else:
                celerities.append(np.nan)

        return np.array(celerities)

    def estimate_wavelength(self, period: float, celerity: float) -> float:
        """
        Estimate wavelength from period and celerity.

        Parameters
        ----------
        period : float
            Wave period (seconds)
        celerity : float
            Wave celerity (m/s)

        Returns
        -------
        float
            Wavelength (meters)
        """
        return period * celerity

    def compute_breaker_location(self, wave_heights: np.ndarray,
                                cross_shore_positions: np.ndarray) -> Tuple[float, int]:
        """
        Estimate wave breaking location.

        Parameters
        ----------
        wave_heights : np.ndarray
            Wave heights along cross-shore
        cross_shore_positions : np.ndarray
            Cross-shore positions (pixels or meters)

        Returns
        -------
        Tuple[float, int]
            Breaking position and index
        """
        # Find maximum wave height (breaking point)
        if len(wave_heights) == 0:
            return np.nan, -1

        # Smooth wave heights
        from scipy.ndimage import gaussian_filter1d
        smoothed_heights = gaussian_filter1d(wave_heights, sigma=2)

        # Find maximum
        break_idx = np.argmax(smoothed_heights)
        break_position = cross_shore_positions[break_idx]

        return break_position, break_idx

    def compute_roller_length(self, timestack: np.ndarray,
                             break_location: int,
                             threshold: float = 0.5) -> float:
        """
        Estimate roller length (foam region after breaking).

        Parameters
        ----------
        timestack : np.ndarray
            2D timestack array
        break_location : int
            Index of breaking location
        threshold : float, optional
            Intensity threshold for foam detection (default: 0.5)

        Returns
        -------
        float
            Roller length in pixels
        """
        if break_location < 0 or break_location >= timestack.shape[0]:
            return 0.0

        # Average intensity after breaking
        if break_location < timestack.shape[0] - 1:
            post_break = timestack[break_location:, :]
            mean_intensity = np.mean(post_break, axis=1)

            # Find where intensity drops below threshold
            below_threshold = mean_intensity < (threshold * np.max(mean_intensity))

            if np.any(below_threshold):
                roller_end = np.argmax(below_threshold)
                return float(roller_end)

        return 0.0

    def compute_spectral_moments(self, timeseries: np.ndarray) -> Dict:
        """
        Compute spectral moments for wave analysis.

        Parameters
        ----------
        timeseries : np.ndarray
            1D time series

        Returns
        -------
        Dict
            Spectral moments (m0, m1, m2)
        """
        # Compute power spectral density
        frequencies, psd = signal.welch(timeseries, fs=1/self.dt,
                                       nperseg=min(256, len(timeseries)))

        # Compute moments
        df = frequencies[1] - frequencies[0]

        m0 = np.sum(psd) * df  # Zeroth moment (variance)
        m1 = np.sum(frequencies * psd) * df  # First moment
        m2 = np.sum(frequencies**2 * psd) * df  # Second moment

        return {
            'm0': m0,
            'm1': m1,
            'm2': m2,
            'significant_height': 4 * np.sqrt(m0)  # Hs = 4 * sqrt(m0)
        }

    def analyze_wave_groups(self, timeseries: np.ndarray) -> Dict:
        """
        Analyze wave grouping characteristics.

        Parameters
        ----------
        timeseries : np.ndarray
            1D time series

        Returns
        -------
        Dict
            Wave grouping parameters
        """
        # Compute envelope using Hilbert transform
        analytic_signal = signal.hilbert(timeseries)
        envelope = np.abs(analytic_signal)

        # Find peaks in envelope (wave groups)
        min_distance = int(self.min_period / self.dt)
        peaks = local_maxima(envelope, min_distance=min_distance)

        # Compute group statistics
        if len(peaks) > 1:
            group_intervals = np.diff(peaks) * self.dt
            mean_group_period = np.mean(group_intervals)
        else:
            mean_group_period = np.nan

        return {
            'envelope': envelope,
            'num_groups': len(peaks),
            'mean_group_period': mean_group_period
        }
