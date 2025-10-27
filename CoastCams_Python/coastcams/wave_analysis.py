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
from .roller_detection import roller_properties_taller


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

        # Step 4: Compute wave height using RollerPropertiesTaller + BreakerHeight
        # This now properly implements the MATLAB workflow
        try:
            # Detect wave breaking events
            print(f"  Detecting wave breaking events...")
            PosX, PosT, Lw, B_processed, Breakstd, Breakmean1, Breakmean2 = roller_properties_taller(
                timestack, self.dt
            )

            print(f"  Found {len(PosX)} breaking events")

            if len(PosX) > 0 and not np.any(np.isnan(PosX)):
                # Compute wave height using photogrammetric method
                print(f"  Computing photogrammetric wave height for {len(PosX)} events...")

                # Create dx array (spatial resolution at each position)
                # MATLAB: dx is a vector that can vary with position
                # For uniform spacing: dx = constant array
                if len(cross_shore_positions) > 1:
                    # Compute dx from differences (MATLAB: X1 = cumsum(dx))
                    dx_array = np.diff(cross_shore_positions)
                    # Extend to same length by appending last value
                    dx_array = np.append(dx_array, dx_array[-1])
                else:
                    dx_array = np.array([self.pixel_resolution])

                hs, hm = self._breaker_height(
                    B_processed, PosT, PosX, Lw,
                    cross_shore_positions, dx_array
                )

                print(f"  Photogrammetric result: hs={hs}, hm={hm}")

                if not np.isnan(hs) and hs > 0:
                    results['mean_Hs'] = hs
                    results['mean_Hm'] = hm
                    print(f"  Wave height (photogrammetric): Hs = {hs:.3f}m, Hm = {hm:.3f}m")
                else:
                    # Fallback to time-series method
                    results['mean_Hs'] = Hs_from_ts * 2.0
                    results['mean_Hm'] = results['mean_Hs'] * 0.7
                    print(f"  Wave height (time series fallback - photogrammetric returned NaN): Hs = {results['mean_Hs']:.3f}m")
            else:
                # No breaking events detected, use fallback
                results['mean_Hs'] = Hs_from_ts * 2.0
                results['mean_Hm'] = results['mean_Hs'] * 0.7
                print(f"  Wave height (time series fallback): Hs = {results['mean_Hs']:.3f}m")

        except Exception as e:
            print(f"  Warning: Roller detection failed: {str(e)}")
            import traceback
            traceback.print_exc()
            results['mean_Hs'] = Hs_from_ts * 2.0
            results['mean_Hm'] = results['mean_Hs'] * 0.7
            print(f"  Wave height (time series fallback): Hs = {results['mean_Hs']:.3f}m")

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

    def _breaker_height(self, B: np.ndarray, PosT: np.ndarray, PosX: np.ndarray,
                       Lw: np.ndarray, X1: np.ndarray, dx: np.ndarray) -> Tuple[float, float]:
        """
        Calculate wave height at breaking point using photogrammetric method.

        Matches MATLAB BreakerHeight function (lines 392-455).

        Parameters
        ----------
        B : np.ndarray
            Preprocessed timestack (time × space)
        PosT : np.ndarray
            Temporal positions of breaking waves
        PosX : np.ndarray
            Spatial positions of breaking waves
        Lw : np.ndarray
            Roller lengths array
        X1 : np.ndarray
            Cross-shore positions in meters
        dx : np.ndarray
            Spatial resolution array (m/pixel) - can vary with position

        Returns
        -------
        hs : float
            Significant wave height (m)
        hm : float
            Mean wave height (m)
        """
        try:
            # Wave face angle at breaking (MATLAB line 403)
            AngleWaveFront = 35 * np.pi / 180

            # Camera viewing angle tangent for each breaking position (MATLAB line 402)
            # AngleCam = abs(z0./X1(PosX))
            AngleCam = np.abs(self.camera_height / X1[PosX.astype(int)])

            print(f"    DEBUG _breaker_height: B.shape={B.shape}, len(PosT)={len(PosT)}, len(PosX)={len(PosX)}")
            print(f"    DEBUG: camera_height={self.camera_height:.2f}m")
            print(f"    DEBUG: X1 range: {np.min(X1):.2f} to {np.max(X1):.2f}m")
            print(f"    DEBUG: PosX range: {np.min(PosX):.0f} to {np.max(PosX):.0f}")

            # Initialize wave height array
            L1 = []

            # Process each breaking event (MATLAB lines 408-423)
            nan_count = 0
            valid_count = 0
            for i in range(len(PosT)):
                try:
                    # Extract window around this breaking event (±25 frames, offshore from break point)
                    t_start = max(0, int(PosT[i]) - 25)
                    t_end = min(B.shape[0], int(PosT[i]) + 25)
                    x_start = max(0, int(PosX[i]) - 50)
                    x_end = B.shape[1]

                    window = B[t_start:t_end, x_start:x_end]

                    if window.size == 0:
                        L1.append(np.nan)
                        nan_count += 1
                        continue

                    # Compute max-min range along time for each spatial position
                    # MATLAB: vec = FilterMean(FilterMean(nanmax(B(...)) - nanmin(B(...)), 20), 5)
                    spatial_ranges = []
                    for x_idx in range(window.shape[1]):
                        ts = window[:, x_idx]
                        if np.sum(~np.isnan(ts)) > len(ts) / 2:
                            range_val = np.nanmax(ts) - np.nanmin(ts)
                            spatial_ranges.append(range_val)
                        else:
                            spatial_ranges.append(np.nan)

                    if len(spatial_ranges) == 0:
                        L1.append(np.nan)
                        continue

                    vec = np.array(spatial_ranges)

                    # Filter twice (MATLAB lines use FilterMean twice with windows 20 and 5)
                    from .roller_detection import filter_mean
                    vec = filter_mean(filter_mean(vec, 20), 5)

                    # Find peak in first third (MATLAB line 411)
                    peak_idx = np.argmax(vec[:len(vec) // 3])

                    # Find high intensity regions (MATLAB line 412)
                    threshold = np.nanmin(vec) + 0.5 * (np.nanmax(vec) - np.nanmin(vec))
                    high_idx = np.where(vec > threshold)[0]

                    if len(high_idx) < 2:
                        L1.append(np.nan)
                        continue

                    # Find gaps (MATLAB line 413)
                    gaps = np.where(np.diff(high_idx) > 20)[0]
                    boundaries = np.sort(np.concatenate([
                        [high_idx[0]],
                        high_idx[gaps] if len(gaps) > 0 else [],
                        high_idx[gaps + 1] if len(gaps) > 0 else [],
                        [high_idx[-1]]
                    ]))

                    # Find boundaries around peak (MATLAB lines 414-418)
                    before_peak = boundaries[boundaries < peak_idx]
                    after_peak = boundaries[boundaries > peak_idx]

                    if len(before_peak) > 0 and len(after_peak) > 0:
                        idl = before_peak[-1]
                        idp = after_peak[0]
                        b = idp - idl

                        # Convert pixel distance to meters (MATLAB line 419)
                        # MATLAB: L1(i) = abs(b).*dx(PosX(i))' - access dx at breaking position
                        pos_idx = int(PosX[i])
                        dx_at_pos = dx[pos_idx] if pos_idx < len(dx) else dx[-1]
                        L1.append(np.abs(b) * dx_at_pos)
                        valid_count += 1
                    else:
                        L1.append(np.nan)
                        nan_count += 1

                except:
                    L1.append(np.nan)
                    nan_count += 1

            print(f"    DEBUG: Processed {len(PosT)} events: valid={valid_count}, nan={nan_count}")

            L1 = np.array(L1)

            print(f"    DEBUG: L1 computed, len={len(L1)}, valid={np.sum(~np.isnan(L1))}")
            if len(L1) > 0:
                print(f"    DEBUG: L1 range: {np.nanmin(L1):.3f} to {np.nanmax(L1):.3f}m")

            # Photogrammetric correction (MATLAB lines 425-426)
            cor = L1 * AngleCam / np.tan(AngleWaveFront)
            Lf = (L1 - cor) * AngleCam

            print(f"    DEBUG: Lf computed, valid={np.sum((Lf > 0) & (~np.isnan(Lf)))}")
            if np.sum((Lf > 0) & (~np.isnan(Lf))) > 0:
                print(f"    DEBUG: Lf range (valid): {np.nanmin(Lf[Lf > 0]):.3f} to {np.nanmax(Lf[Lf > 0]):.3f}m")

            # Filter valid heights (MATLAB lines 428-430)
            ind = np.where((Lf > 0) & (~np.isnan(Lf)))[0]

            if len(ind) == 0:
                print(f"    DEBUG: No valid Lf values! Returning NaN")
                return np.nan, np.nan

            Lord = np.sort(Lf[ind])[::-1]  # Sort descending

            # Remove outliers: keep middle 80% (MATLAB line 430)
            start_idx = max(1, round(len(Lord) / 10))
            end_idx = min(len(Lord) - 1, round(9 * len(Lord) / 10))
            Lord = Lord[start_idx:end_idx]

            if len(Lord) == 0:
                return np.nan, np.nan

            # Significant wave height: median of top third (MATLAB line 432)
            n_third = max(1, round(len(Lord) / 3))
            hs = np.nanmedian(Lord[:n_third])

            # Mean wave height: median of all (MATLAB line 433)
            hm = np.nanmedian(Lord)

            # Alternative hm calculation (MATLAB lines 435-448)
            # Accumulate spatial profiles across all breaking events
            try:
                vecm = np.zeros(350)
                for i in range(len(PosT)):
                    try:
                        # Extract window for this event
                        t_start = max(0, int(PosT[i]) - 25)
                        t_end = min(B.shape[0], int(PosT[i]) + 25)
                        x_start = max(0, int(PosX[i]) - 100)
                        x_end = B.shape[1]

                        window = B[t_start:t_end, x_start:x_end]

                        if window.size > 0:
                            # Compute max-min range along time for each spatial position
                            spatial_ranges = []
                            for x_idx in range(window.shape[1]):
                                ts = window[:, x_idx]
                                if np.sum(~np.isnan(ts)) > len(ts) / 2:
                                    range_val = np.nanmax(ts) - np.nanmin(ts)
                                    spatial_ranges.append(range_val)
                                else:
                                    spatial_ranges.append(np.nan)

                            if len(spatial_ranges) > 0:
                                vec = np.array(spatial_ranges)

                                # Filter twice (MATLAB FilterMean with windows 20 and 5)
                                from .roller_detection import filter_mean
                                vec = filter_mean(filter_mean(vec, 20), 5)

                                # Accumulate (limited to 350 points)
                                n_pts = min(len(vec), 350, len(vecm))
                                vecm[:n_pts] += vec[:n_pts]
                    except:
                        continue

                # Find max and min of derivative
                if np.sum(vecm != 0) > 10:
                    diff_vecm = np.diff(vecm)
                    gtf = np.argmax(diff_vecm)  # Position of max increase
                    gtd = np.argmin(diff_vecm)  # Position of max decrease

                    # Compute L1 from distance between these points
                    pos_idx = int(np.nanmedian(PosX))
                    dx_at_pos = dx[pos_idx] if pos_idx < len(dx) else dx[-1]
                    L1_alt = dx_at_pos * abs(gtd - gtf)

                    # Photogrammetric correction (same as main calculation)
                    AngleCam_mean = np.nanmean(AngleCam)
                    cor_alt = L1_alt * AngleCam_mean / np.tan(AngleWaveFront)
                    hm_alt = (L1_alt - cor_alt) * AngleCam_mean

                    # Use alternative if it's valid and different from primary
                    if not np.isnan(hm_alt) and hm_alt > 0:
                        # Average the two methods for robustness
                        hm = (hm + hm_alt) / 2.0
            except:
                pass  # Keep primary hm if alternative fails

            return hs, hm

        except Exception as e:
            print(f"    BreakerHeight error: {str(e)[:100]}")
            return np.nan, np.nan

    def _compute_wave_height_photogrammetric_old(self, timestack: np.ndarray,
                                             cross_shore_positions: np.ndarray,
                                             break_idx: int) -> float:
        """
        OLD VERSION - DEPRECATED
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

    def linear_c(self, T: float, c: np.ndarray, precision: float = 0.01) -> np.ndarray:
        """
        Calculate water depth for linear waves using dispersion relation.

        Matches MATLAB LinearC function (lines 135-156).
        Iteratively solves: ω² = g*k*tanh(k*d) for depth d.

        Parameters
        ----------
        T : float
            Wave period (seconds)
        c : np.ndarray
            Phase speed / wave celerity (m/s)
        precision : float, optional
            Convergence precision (default: 0.01)

        Returns
        -------
        np.ndarray
            Water depth (meters) for each celerity value
        """
        df = np.zeros_like(c)

        for i in range(len(c)):
            if np.isnan(c[i]) or c[i] <= 0:
                df[i] = np.nan
                continue

            # Initialize
            w = 2 * np.pi / T  # Angular frequency
            k = w / c[i]  # Wave number
            g = self.gravity  # Gravity
            do = 1000  # Arbitrary large value
            d = c[i]**2 / g  # Initial guess

            # Newton-Raphson iteration
            ct = 0
            max_iter = 100
            while abs(do - d) > precision and ct < max_iter:
                ct += 1
                do = d
                # Dispersion relation: ω² = g*k*tanh(k*d)
                dispe = w**2 - g * k * np.tanh(k * d)
                # Derivative for Newton-Raphson
                fdispe = -g * (k**2) / (np.cosh(k * d)**2)
                # Update depth
                d = d - dispe / fdispe

            df[i] = d

        return df

    def radon_c_indiv(self, dt: float, dx: np.ndarray, In: np.ndarray, Tm: float) -> np.ndarray:
        """
        Compute wave celerity using Radon transform on windowed segments.

        Matches MATLAB RadonCIndiv_20140903 function (lines 457-540).

        Parameters
        ----------
        dt : float
            Temporal resolution (seconds)
        dx : np.ndarray
            Spatial resolution array (m/pixel)
        In : np.ndarray
            Input timestack (time × space)
        Tm : float
            Wave period (seconds)

        Returns
        -------
        np.ndarray
            Wave celerity for each spatial position (m/s)
        """
        from skimage.transform import radon
        from scipy.ndimage import uniform_filter1d

        # Spatial and temporal resolution degradation
        pdx = 2  # Spatial decimation factor
        pdt = 1  # Temporal decimation factor

        # Temporal and spatial frequency
        freq_t = 1.0 / dt
        if isinstance(dx, np.ndarray):
            freq_x = 1.0 / (dx * pdx)
        else:
            freq_x = 1.0 / (dx * pdx)

        # Degrade resolution and compute difference
        # MATLAB: M=detrend(abs(diff(In(1:pdt:length(In(:,1)),1:pdx:length(In(1,:))),1))')';
        In_decimated = In[::pdt, ::pdx]
        M = np.abs(np.diff(In_decimated, axis=0))
        M = signal.detrend(M, axis=0)

        # Window size in space
        Wx = max(10, round(M.shape[1] / 10))

        # Initialize output
        CRadmoyt = np.full(M.shape[1], np.nan)

        # Loop over spatial positions
        for ix in range(Wx + 1, M.shape[1] - Wx, Wx - 1):
            try:
                # Extract window
                MR = signal.detrend(M[:, ix - Wx:ix + Wx + 1], axis=0)

                nt = MR.shape[0]
                if nt < 10:
                    continue

                # Radon transform
                iang = np.arange(1, 181)
                R = radon(MR, theta=iang, circle=False)

                nr = R.shape[0]
                amp = nr / nt

                # Compute resolution for each angle
                k = nt
                trk = np.floor(MR.shape[0] / 2) - np.floor(0 * np.cos(np.deg2rad(iang)) +
                                                            (MR.shape[0] / 2 - k * amp) * np.sin(np.deg2rad(iang)))
                trk = trk - np.min(trk)
                trk[trk == 0] = 1  # Avoid division by zero
                res = (nt * dt) / (trk * 2)

                # Filter R2 based on wave period
                R2 = np.zeros_like(R)
                for i, angle in enumerate(iang):
                    if res[i] > 0:
                        # Smooth with period-dependent windows
                        win1 = max(1, round((Tm + 2) / res[i]))
                        win2 = max(1, round((Tm - 2) / res[i]))

                        # Apply smoothing
                        R_smooth1 = uniform_filter1d(R[:, i], size=int(win1), mode='nearest')
                        R_detrend = R[:, i] - R_smooth1
                        R2[:, i] = uniform_filter1d(R_detrend, size=int(win2), mode='nearest')
                    else:
                        R2[:, i] = R[:, i]

                # Find peak angle in first 90 degrees (incident waves)
                vec = np.nanstd(R2[round(R2.shape[0] / 4):round(3 * R2.shape[0] / 4), :90], axis=0)
                vec = uniform_filter1d(vec, size=10, mode='nearest')
                a2 = np.argmax(vec)
                frd = vec[a2]

                # Check if peak is significant
                if vec[a2] > 1.15 * vec[89]:
                    # Compute celerity from angle
                    # MATLAB: C2=(1/mean(freq_x))/(tand(90-a2)*(1/mean(freq_t)))
                    angle_deg = iang[a2]
                    if isinstance(freq_x, np.ndarray):
                        freq_x_local = np.nanmedian(freq_x[pdx * (ix - Wx):pdx * (ix + Wx + 1)])
                    else:
                        freq_x_local = freq_x

                    tan_angle = np.tan(np.deg2rad(90 - angle_deg))
                    if tan_angle != 0:
                        C2 = (1.0 / freq_x_local) / (tan_angle * (1.0 / freq_t))
                    else:
                        C2 = np.nan
                else:
                    C2 = np.nan

                # Assign to output range
                CRadmoyt[ix - Wx:ix + Wx + 1] = C2

            except Exception as e:
                continue

        # Interpolate to fill gaps
        notnan = np.where((~np.isnan(CRadmoyt)) & (CRadmoyt > 0))[0]

        if isinstance(dx, np.ndarray):
            CRadmoy = np.full_like(dx, np.nan, dtype=float)
        else:
            CRadmoy = np.full(In.shape[1], np.nan, dtype=float)

        if len(notnan) > 1:
            try:
                min_idx = np.min(notnan * pdx)
                max_idx = np.max(notnan * pdx)
                if max_idx < len(CRadmoy):
                    x_interp = np.arange(min_idx, max_idx + 1)
                    CRadmoy[min_idx:max_idx + 1] = np.interp(
                        x_interp,
                        pdx * notnan,
                        CRadmoyt[notnan]
                    )
            except:
                pass

        return CRadmoy

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
