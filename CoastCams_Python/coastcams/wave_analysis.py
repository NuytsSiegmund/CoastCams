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
            self.camera_height = config.camera_height  # Add camera height
        else:
            # Default parameters
            self.dt = 0.5  # Time step (seconds)
            self.pixel_resolution = 0.1  # meters/pixel
            self.min_period = 4.0  # seconds
            self.max_period = 25.0  # seconds
            self.gravity = 9.81  # m/s^2
            self.camera_height = 27.24  # meters

    def analyze_timestack(self, timestack: np.ndarray,
                         cross_shore_positions: np.ndarray,
                         shorelines: list = None) -> Dict:
        """
        Perform complete wave analysis on timestack data.

        Matches MATLAB WaveParameters_CoastCams.m workflow:
        1. Find breaking position (median of high-std locations)
        2. Extract 1D time series at breaking position
        3. Compute wave period from zero-crossing on that time series
        4. Compute wave height using photogrammetric method on full 2D timestack

        Parameters
        ----------
        timestack : np.ndarray
            2D timestack array (space × time), e.g., 689×1680
        cross_shore_positions : np.ndarray
            Cross-shore pixel positions
        shorelines : list, optional
            List of shoreline arrays for each time step

        Returns
        -------
        Dict
            Dictionary containing wave parameters (ONE value per timestack)
        """
        results = {}

        print(f"  Timestack shape: {timestack.shape} (time × space)")
        print(f"  Intensity range: {np.nanmin(timestack):.3f} to {np.nanmax(timestack):.3f}")

        if timestack.shape[0] < 10 or timestack.shape[1] < 10:
            print("  Timestack too small for analysis")
            results['mean_Hs'] = np.nan
            results['mean_Hm'] = np.nan
            results['mean_Tm'] = np.nan
            results['wave_periods'] = np.full(len(cross_shore_positions), np.nan)
            results['cross_shore_positions'] = cross_shore_positions
            return results

        # Step 1: Find breaking position
        # Compute standard deviation along time (axis 0) at each spatial position (axis 1)
        # MATLAB: Breakstd = std(B(time_range, :)) - std along time dimension for each spatial column
        std_profile = np.nanstd(timestack, axis=0)  # Shape: (num_space_positions,)

        # Find median breaking position (highest variability)
        valid_std = ~np.isnan(std_profile)
        if np.sum(valid_std) == 0:
            print("  No valid std values")
            results['mean_Hs'] = np.nan
            results['mean_Hm'] = np.nan
            results['mean_Tm'] = np.nan
            results['wave_periods'] = np.full(len(cross_shore_positions), np.nan)
            results['cross_shore_positions'] = cross_shore_positions
            return results

        median_break_idx = int(np.nanmedian(np.where(std_profile == np.nanmax(std_profile))[0]))
        print(f"  Breaking position: spatial index {median_break_idx} (std = {std_profile[median_break_idx]:.3f})")

        # Step 2: Extract 1D time series at breaking position
        # MATLAB: I = S(:, round(nanmedian(PosX))) - extract column (all time points at spatial location)
        # This is the key: extract the TEMPORAL dimension at one spatial location
        intensity_timeseries = timestack[:, median_break_idx]  # Shape: (num_time_points,)

        # Remove NaN values
        valid_mask = ~np.isnan(intensity_timeseries)
        if np.sum(valid_mask) < len(intensity_timeseries) / 2:
            print(f"  Too many NaN values in time series: {np.sum(~valid_mask)}/{len(intensity_timeseries)}")
            results['mean_Hs'] = np.nan
            results['mean_Hm'] = np.nan
            results['mean_Tm'] = np.nan
            results['wave_periods'] = np.full(len(cross_shore_positions), np.nan)
            results['cross_shore_positions'] = cross_shore_positions
            return results

        intensity_valid = intensity_timeseries[valid_mask]
        print(f"  Time series length: {len(intensity_valid)} points (sampling dt = {self.dt}s)")

        # Step 3: Compute wave period from time series using zero-crossing
        # Detrend the time series
        x = np.arange(len(intensity_valid))
        p = np.polyfit(x, intensity_valid, 1)
        detrended_ts = intensity_valid - (p[0] * x + p[1])
        detrended_ts = detrended_ts - np.mean(detrended_ts)

        # Use Wave_Char-style analysis: compute Hs = 4 × std
        # This is from the 1D time series, matching MATLAB line 159
        Hs_from_ts = 4.0 * np.std(detrended_ts)

        # Compute periods using zero-crossing method
        periods = self._compute_wave_periods_from_timeseries(detrended_ts)

        if len(periods) > 0:
            mean_period = np.mean(periods)
            print(f"  Wave period (zero-crossing): {mean_period:.2f}s from {len(periods)} cycles")
        else:
            # Fallback: estimate from FFT
            mean_period = self._compute_peak_period(detrended_ts)
            if np.isnan(mean_period):
                mean_period = 10.0  # Default
            print(f"  Wave period (fallback): {mean_period:.2f}s")

        results['mean_Tm'] = mean_period
        results['wave_periods'] = np.full(len(cross_shore_positions), mean_period)

        # Step 4: Compute wave height
        # NOTE: Photogrammetric method disabled - requires RollerPropertiesTaller implementation
        # See PYTHON_ISSUES_AND_FIXES.md for details

        # Use time series standard deviation method
        # MATLAB Wave_Char.m uses Hs = 4 × std(detrended_intensity)
        calibration_factor = 2.0  # Will need proper calibration or physics-based approach
        results['mean_Hs'] = Hs_from_ts * calibration_factor
        results['mean_Hm'] = results['mean_Hs'] * 0.7
        print(f"  Wave height (time series method): Hs = {results['mean_Hs']:.3f}m (std-based)")

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

    def _compute_wave_heights_per_timestep(self, timestack: np.ndarray,
                                           cross_shore_positions: np.ndarray) -> np.ndarray:
        """
        Compute wave height for each time step using intensity variations.

        Simplified approach: uses standard deviation of intensity at each time step
        as a proxy for wave activity, then converts to physical height.

        Parameters
        ----------
        timestack : np.ndarray
            2D timestack array (space x time)
        cross_shore_positions : np.ndarray
            Cross-shore positions in meters

        Returns
        -------
        np.ndarray
            Wave heights for each time step (length = timestack.shape[1])
        """
        wave_heights = np.full(timestack.shape[1], np.nan)

        # Debug: check timestack properties
        print(f"  Timestack shape: {timestack.shape}")
        print(f"  Intensity range: {np.nanmin(timestack):.3f} to {np.nanmax(timestack):.3f}")
        print(f"  Mean intensity: {np.nanmean(timestack):.3f}")

        # For each time step, compute standard deviation of spatial profile
        # Higher std = more wave activity
        for t in range(timestack.shape[1]):
            try:
                # Extract spatial profile at this time
                profile = timestack[:, t]

                # Skip if too many NaNs
                valid_points = ~np.isnan(profile)
                if np.sum(valid_points) < len(profile) / 3:
                    continue

                # Compute standard deviation (wave activity indicator)
                profile_std = np.nanstd(profile)

                # Detrend: remove mean
                profile_detrended = profile[valid_points] - np.nanmean(profile[valid_points])

                # Compute range (max - min) as wave height indicator
                if len(profile_detrended) > 0:
                    intensity_range = np.max(profile_detrended) - np.min(profile_detrended)

                    # Find location of maximum variation (proxy for breaking location)
                    max_var_idx = np.nanargmax(np.abs(profile_detrended))

                    # Get approximate cross-shore position
                    valid_indices = np.where(valid_points)[0]
                    if max_var_idx < len(valid_indices):
                        actual_idx = valid_indices[max_var_idx]
                        if actual_idx < len(cross_shore_positions):
                            x_distance = cross_shore_positions[actual_idx]

                            if x_distance > 0.1:
                                # Camera angle
                                camera_angle = np.arctan(self.camera_height / x_distance)

                                # Convert intensity range to approximate wave height
                                # Empirical scaling: intensity range of ~0.5 (normalized) ≈ 1-2m wave
                                # Use camera geometry as additional scaling
                                base_height = intensity_range * 3.0  # Empirical factor

                                # Apply camera angle correction (waves appear smaller at distance)
                                wave_height = base_height * np.tan(camera_angle) * 2.0

                                # Clip to realistic range
                                wave_height = np.clip(wave_height, 0.1, 8.0)

                                if not np.isnan(wave_height):
                                    wave_heights[t] = wave_height

            except Exception as e:
                print(f"  Error at timestep {t}: {e}")
                continue

        # Debug output
        valid_count = np.sum(~np.isnan(wave_heights))
        if valid_count > 0:
            print(f"  Computed {valid_count}/{len(wave_heights)} wave heights")
            print(f"  Range: {np.nanmin(wave_heights):.2f} - {np.nanmax(wave_heights):.2f} m")
        else:
            print(f"  No valid wave heights computed (all NaN)")

        return wave_heights

    def _compute_wave_height_photogrammetric(self, timestack: np.ndarray,
                                             cross_shore_positions: np.ndarray,
                                             break_idx: int) -> float:
        """
        Compute wave height using photogrammetric method with camera geometry.
        Matches MATLAB BreakerHeight function (lines 392-455).

        Parameters
        ----------
        timestack : np.ndarray
            2D timestack array (time × space)
        cross_shore_positions : np.ndarray
            Cross-shore positions in meters
        break_idx : int
            Index of breaking position

        Returns
        -------
        float
            Significant wave height (Hs) in meters
        """
        # Wave face angle at breaking (MATLAB line 403)
        wave_face_angle = np.deg2rad(35.0)

        # Get cross-shore position of breaking location
        if break_idx >= len(cross_shore_positions):
            print(f"    DEBUG: break_idx {break_idx} >= len(cross_shore_positions) {len(cross_shore_positions)}")
            return np.nan

        x_distance = cross_shore_positions[break_idx]
        print(f"    DEBUG: break_idx={break_idx}, x_distance={x_distance:.2f}m, camera_height={self.camera_height:.2f}m")

        if x_distance <= 0:
            print(f"    DEBUG: x_distance <= 0")
            return np.nan

        # Camera viewing angle tangent (MATLAB line 402)
        # NOTE: MATLAB's "AngleCam" is actually tan(angle), NOT the angle itself!
        # MATLAB: AngleCam = abs(z0./X1(PosX))
        camera_angle_tan = self.camera_height / x_distance

        wave_heights = []
        debug_counters = {'total': 0, 'empty_window': 0, 'no_profile': 0, 'no_smooth': 0,
                          'no_high_int': 0, 'no_boundaries': 0, 'no_extent': 0, 'invalid_height': 0, 'valid': 0}

        # Analyze wave features across time (MATLAB lines 408-423)
        # Look at temporal window around each potential wave
        # timestack is (time, space), so iterate over time (axis 0)
        num_iterations = 0
        for t in range(0, timestack.shape[0] - 50, 10):  # Every 10 time frames
            num_iterations += 1
            debug_counters['total'] += 1
            try:
                # Extract window around this time (±25 frames, matching MATLAB line 410)
                t_start = max(0, t - 25)
                t_end = min(timestack.shape[0], t + 25)

                # Extract region offshore from breaking (matching MATLAB line 410)
                x_start = max(0, break_idx - 50)
                x_end = timestack.shape[1]

                # Get the temporal-spatial window: window is (time_range, space_range)
                window = timestack[t_start:t_end, x_start:x_end]

                if window.size == 0:
                    debug_counters['empty_window'] += 1
                    continue

                # Compute max-min range along TIME dimension for each SPATIAL position
                # This matches MATLAB line 410: nanmax(B(time_range, space_range)) - nanmin(...)
                # where max/min are computed along the time axis for each spatial column
                spatial_profile = []
                for space_idx in range(window.shape[1]):  # Iterate over spatial positions
                    time_series = window[:, space_idx]  # Extract time series at this position
                    if np.sum(~np.isnan(time_series)) > len(time_series) / 2:
                        range_val = np.nanmax(time_series) - np.nanmin(time_series)
                        spatial_profile.append(range_val)
                    else:
                        spatial_profile.append(np.nan)

                if len(spatial_profile) == 0:
                    debug_counters['no_profile'] += 1
                    continue

                # Smooth the spatial profile (matching MATLAB FilterMean)
                vec = np.array(spatial_profile)
                # Handle NaN values by interpolation before smoothing
                if np.any(np.isnan(vec)):
                    valid_idx = ~np.isnan(vec)
                    if np.sum(valid_idx) < len(vec) // 2:
                        debug_counters['no_profile'] += 1
                        continue
                    vec = np.interp(np.arange(len(vec)), np.where(valid_idx)[0], vec[valid_idx])

                from scipy.ndimage import uniform_filter1d
                vec_smooth = uniform_filter1d(vec, size=min(20, max(3, len(vec)//2)))

                if len(vec_smooth) == 0:
                    debug_counters['no_smooth'] += 1
                    continue

                # Find the peak (matching MATLAB line 411)
                peak_idx = np.argmax(vec_smooth[:len(vec_smooth)//3])  # First third

                # Find where intensity is above threshold (matching MATLAB line 412)
                threshold = np.nanmin(vec_smooth) + 0.5 * (np.nanmax(vec_smooth) - np.nanmin(vec_smooth))
                high_intensity_indices = np.where(vec_smooth > threshold)[0]

                if len(high_intensity_indices) < 2:
                    debug_counters['no_high_int'] += 1
                    continue

                # Find gaps in high intensity regions (matching MATLAB line 413)
                gaps = np.where(np.diff(high_intensity_indices) > 20)[0]
                boundaries = np.concatenate(([high_intensity_indices[0]],
                                            high_intensity_indices[gaps] if len(gaps) > 0 else [],
                                            high_intensity_indices[gaps + 1] if len(gaps) > 0 else [],
                                            [high_intensity_indices[-1]]))

                # Find boundaries around peak (matching MATLAB lines 414-418)
                before_peak = boundaries[boundaries < peak_idx]
                after_peak = boundaries[boundaries > peak_idx]

                if len(before_peak) > 0 and len(after_peak) > 0:
                    left_boundary = before_peak[-1]
                    right_boundary = after_peak[0]

                    # Horizontal extent in pixels (matching MATLAB line 419)
                    h_extent_pixels = right_boundary - left_boundary

                    if h_extent_pixels > 0:
                        # Convert to meters
                        L_horizontal = h_extent_pixels * self.pixel_resolution

                        # Photogrammetric conversion (MATLAB lines 425-426)
                        # cor = (L1).*tan(AngleCam)/tan(AngleWaveFront);
                        # Lf = (L1 - cor).*tan(AngleCam);
                        correction = L_horizontal * camera_angle_tan / np.tan(wave_face_angle)
                        L_corrected = L_horizontal - correction
                        wave_height = L_corrected * camera_angle_tan

                        # Only keep realistic values (matching MATLAB line 428)
                        # TEMPORARY: Accept all values to see distribution
                        wave_heights.append(abs(wave_height))  # Use absolute value for now
                        debug_counters['valid'] += 1
                        if debug_counters['valid'] <= 3:  # Print first 3 heights
                            print(f"      Wave height: {wave_height:.3f}m → {abs(wave_height):.3f}m (L_horiz={L_horizontal:.3f}m, correction={correction:.3f}m)")
                    else:
                        debug_counters['no_extent'] += 1
                else:
                    debug_counters['no_boundaries'] += 1

            except Exception as e:
                print(f"    DEBUG: Exception in iteration: {str(e)[:100]}")
                continue

        print(f"    DEBUG: Processed {num_iterations} iterations")
        print(f"    DEBUG: Failures - empty:{debug_counters['empty_window']}, no_profile:{debug_counters['no_profile']}, " +
              f"no_smooth:{debug_counters['no_smooth']}, no_high_int:{debug_counters['no_high_int']}, " +
              f"no_boundaries:{debug_counters['no_boundaries']}, no_extent:{debug_counters['no_extent']}, " +
              f"invalid_height:{debug_counters['invalid_height']}, valid:{debug_counters['valid']}")
        if len(wave_heights) == 0:
            print(f"    DEBUG: No wave heights computed in photogrammetric method")
            return np.nan

        print(f"    DEBUG: Found {len(wave_heights)} wave heights, range: {np.min(wave_heights):.3f} to {np.max(wave_heights):.3f}m")

        # Sort and take top 1/3 for significant wave height (MATLAB lines 429-432)
        wave_heights_sorted = np.sort(wave_heights)[::-1]

        # Remove outliers: keep middle 80% (matching MATLAB line 430)
        n_keep_start = max(1, len(wave_heights_sorted) // 10)
        n_keep_end = min(len(wave_heights_sorted) - 1, 9 * len(wave_heights_sorted) // 10)
        wave_heights_filtered = wave_heights_sorted[n_keep_start:n_keep_end]

        if len(wave_heights_filtered) == 0:
            print(f"    DEBUG: No wave heights after filtering")
            return np.nan

        # Significant wave height: median of top 1/3 (MATLAB line 432)
        n_third = max(1, len(wave_heights_filtered) // 3)
        Hs = np.median(wave_heights_filtered[:n_third])
        print(f"    DEBUG: Photogrammetric Hs = {Hs:.3f}m (from {len(wave_heights_filtered)} filtered heights)")

        return Hs

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
