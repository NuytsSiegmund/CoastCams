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
from .matlab_preprocessing import (RollerDetector, PhotogrammetricHeightCalculator)


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

        # Initialize MATLAB-style processors
        self.roller_detector = RollerDetector(dt=self.dt)
        self.height_calculator = PhotogrammetricHeightCalculator(
            camera_height=self.camera_height,
            pixel_resolution=self.pixel_resolution
        )

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

        # Step 3: Apply MATLAB-style FIR bandpass filtering
        # Matches Get_Periode function (lines 88-114)
        filtered_ts = self._apply_matlab_bandpass_filter(intensity_valid)

        if filtered_ts is None or len(filtered_ts) < 10:
            print(f"  Filtering failed, using unfiltered signal")
            filtered_ts = intensity_valid
            # Simple detrending
            x = np.arange(len(filtered_ts))
            p = np.polyfit(x, filtered_ts, 1)
            filtered_ts = filtered_ts - (p[0] * x + p[1])
            filtered_ts = filtered_ts - np.mean(filtered_ts)

        # Use Wave_Char-style analysis: compute Hs = 4 × std
        # This is from the 1D time series, matching MATLAB line 159
        Hs_from_ts = 4.0 * np.std(filtered_ts)

        # Compute periods using zero-crossing method on filtered signal
        periods = self._compute_wave_periods_from_timeseries(filtered_ts)

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

        # Step 4: Use full MATLAB photogrammetric pipeline
        # Detect breaking waves and compute heights
        try:
            print(f"  Running full MATLAB photogrammetric pipeline...")

            # Detect rollers and breaking waves
            roller_results = self.roller_detector.detect_rollers(timestack)

            PosX = roller_results['PosX']
            PosT = roller_results['PosT']
            Lw = roller_results['Lw']
            B = roller_results['B']

            if len(PosX) > 0:
                # Calculate wave heights using photogrammetric method
                Hs, Hm = self.height_calculator.calculate_wave_heights(
                    B, PosT, PosX, Lw, cross_shore_positions
                )

                if not np.isnan(Hs) and Hs > 0:
                    results['mean_Hs'] = Hs
                    results['mean_Hm'] = Hm
                    print(f"  Wave height (photogrammetric): Hs = {Hs:.3f}m, Hm = {Hm:.3f}m")
                else:
                    # Fallback to time series method
                    calibration_factor = 2.0
                    results['mean_Hs'] = Hs_from_ts * calibration_factor
                    results['mean_Hm'] = results['mean_Hs'] * 0.7
                    print(f"  Wave height (fallback - time series): Hs = {results['mean_Hs']:.3f}m")
            else:
                # No breaking waves detected, use fallback
                print(f"  No breaking waves detected, using fallback method")
                calibration_factor = 2.0
                results['mean_Hs'] = Hs_from_ts * calibration_factor
                results['mean_Hm'] = results['mean_Hs'] * 0.7
                print(f"  Wave height (fallback): Hs = {results['mean_Hs']:.3f}m")

        except Exception as e:
            print(f"  Photogrammetric pipeline error: {e}")
            # Fallback to time series method
            calibration_factor = 2.0
            results['mean_Hs'] = Hs_from_ts * calibration_factor
            results['mean_Hm'] = results['mean_Hs'] * 0.7
            print(f"  Wave height (error fallback): Hs = {results['mean_Hs']:.3f}m")

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
        Compute wave height using simplified photogrammetric method.
        Matches MATLAB BreakerHeight approach but simplified without individual wave detection.

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
            return np.nan

        x_distance = cross_shore_positions[break_idx]

        if x_distance <= 0 or x_distance < 1.0:  # Too close to camera
            return np.nan

        # Camera viewing angle (MATLAB line 402)
        camera_angle = np.arctan(self.camera_height / x_distance)

        # Simplified approach: analyze intensity variations in breaking region
        # Extract region around breaking location (±50 spatial positions)
        x_start = max(0, break_idx - 50)
        x_end = min(timestack.shape[1], break_idx + 50)

        breaking_region = timestack[:, x_start:x_end]  # (time, space_window)

        if breaking_region.size == 0:
            return np.nan

        # Compute std along spatial dimension for each time point
        spatial_std = np.nanstd(breaking_region, axis=1)  # (time,)

        # Take 90th percentile of spatial std (representative of wave activity)
        wave_intensity_range = np.percentile(spatial_std[~np.isnan(spatial_std)], 90) if np.sum(~np.isnan(spatial_std)) > 0 else np.nan

        if np.isnan(wave_intensity_range) or wave_intensity_range < 0.01:
            return np.nan

        # Convert intensity variation to horizontal wave extent
        # Empirical: std of 0.1 in normalized intensity ≈ 10-20 pixels of wave face
        # This accounts for the spatial variation caused by wave breaking
        horizontal_extent_pixels = wave_intensity_range * 150.0  # Empirical scaling

        if horizontal_extent_pixels < 1.0:
            return np.nan

        # Convert to meters
        L_horizontal = horizontal_extent_pixels * self.pixel_resolution

        # Photogrammetric conversion (MATLAB lines 425-426)
        # correction = L × tan(camera_angle) / tan(wave_face_angle)
        correction = L_horizontal * np.tan(camera_angle) / np.tan(wave_face_angle)
        L_corrected = L_horizontal - correction

        # Vertical wave height = L_corrected × tan(camera_angle)
        wave_height = L_corrected * np.tan(camera_angle)

        # Only keep realistic values
        if 0.1 < wave_height < 5.0:
            return wave_height
        else:
            return np.nan

    def _apply_matlab_bandpass_filter(self, timeseries: np.ndarray) -> np.ndarray:
        """
        Apply MATLAB-style FIR bandpass filtering.
        Matches Get_Periode function from WaveParameters_CoastCams.m (lines 88-114).

        Creates a bandpass filter for wave periods between 1.5s and 20s by:
        1. Low-pass filter at 1.5s cutoff → removes high-frequency noise
        2. High-pass filter at 20s cutoff → removes long-period trends

        Parameters
        ----------
        timeseries : np.ndarray
            Input time series

        Returns
        -------
        np.ndarray
            Filtered and detrended time series, or None if filtering fails
        """
        try:
            from scipy.signal import firwin, filtfilt

            # MATLAB: Tcoupure = 1.5 (cutoff period in seconds)
            # MATLAB: fr = 1/dt (sampling frequency)
            # MATLAB: Val = (1/Tcoupure)*2*(1/fr) = (1/1.5) * 2 * dt
            # MATLAB: fil = fir1(ord, Val, 'low')

            sampling_freq = 1.0 / self.dt  # Hz
            S = timeseries.copy()

            # Step 1: Low-pass filter (removes periods < 1.5s)
            T_cutoff_low = 1.5  # seconds
            f_cutoff_low = 1.0 / T_cutoff_low  # Hz
            numtaps = min(501, len(S) // 4)  # Filter order (reduce from MATLAB's 1000 to avoid filtfilt issues)

            # Ensure odd number of taps for FIR filter
            if numtaps % 2 == 0:
                numtaps += 1

            if numtaps < 3 or len(S) < numtaps * 3:  # filtfilt needs 3x filter length
                return None

            # Normalize cutoff frequency
            nyquist = sampling_freq / 2.0
            normalized_cutoff = f_cutoff_low / nyquist

            if normalized_cutoff >= 1.0:
                normalized_cutoff = 0.99

            # Create low-pass FIR filter
            fir_coeff_low = firwin(numtaps, normalized_cutoff, window='hamming')

            # Apply filter (filtfilt for zero-phase filtering, like MATLAB's approach)
            S_filtered = filtfilt(fir_coeff_low, [1.0], S)

            # Detrend (MATLAB line 101)
            x = np.arange(len(S_filtered))
            p = np.polyfit(x, S_filtered, 1)
            S_detrended = S_filtered - (p[0] * x + p[1])

            # Step 2: High-pass filter (removes periods > 20s)
            T_cutoff_high = 20.0  # seconds
            f_cutoff_high = 1.0 / T_cutoff_high  # Hz
            normalized_cutoff_high = f_cutoff_high / nyquist

            if normalized_cutoff_high >= 1.0:
                normalized_cutoff_high = 0.99

            # Create high-pass FIR filter
            fir_coeff_high = firwin(numtaps, normalized_cutoff_high, window='hamming', pass_zero=False)

            # Apply filter
            S_filtered2 = filtfilt(fir_coeff_high, [1.0], S_detrended)

            # Final detrend (MATLAB line 114)
            x = np.arange(len(S_filtered2))
            p = np.polyfit(x, S_filtered2, 1)
            S_final = S_filtered2 - (p[0] * x + p[1])

            # Remove mean (matching MATLAB preprocessing)
            S_final = S_final - np.mean(S_final)

            return S_final

        except Exception as e:
            print(f"  FIR filtering error: {e}")
            return None

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
