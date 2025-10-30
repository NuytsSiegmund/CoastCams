"""
Cross-correlation analysis for CoastCams.

This module computes wave properties from spatiotemporal patterns
using cross-correlation techniques.
"""

import numpy as np
from scipy import signal
from typing import Dict, Tuple, Optional, List


class CrossCorrelationAnalyzer:
    """
    Analyze wave properties using cross-correlation.

    Computes wave celerity, wavelength, and period from
    spatiotemporal correlation patterns in timestack data.
    """

    def __init__(self, config=None):
        """
        Initialize cross-correlation analyzer.

        Parameters
        ----------
        config : CoastCamsConfig, optional
            Configuration object
        """
        if config is not None:
            self.dt = 1.0 / config.acquisition_frequency
            self.dx = config.pixel_resolution
            self.correlation_spacing = config.correlation_spacing
            self.time_lag = config.time_lag
        else:
            # Default parameters
            self.dt = 0.5  # Time step (seconds)
            self.dx = 0.1  # Spatial resolution (meters/pixel)
            self.correlation_spacing = 100  # Pixel spacing
            self.time_lag = 1  # Frame lag

    def analyze_timestack(self, timestack: np.ndarray) -> Dict:
        """
        Perform cross-correlation analysis on timestack.

        Parameters
        ----------
        timestack : np.ndarray
            2D timestack array (space x time)

        Returns
        -------
        Dict
            Cross-correlation results including celerity and wavelength
        """
        results = {}

        # Compute spatial cross-correlation
        spatial_cc = self._compute_spatial_correlation(timestack)
        results['spatial_correlation'] = spatial_cc

        # Compute temporal cross-correlation
        temporal_cc = self._compute_temporal_correlation(timestack)
        results['temporal_correlation'] = temporal_cc

        # Extract wave properties
        wave_props = self._extract_wave_properties(timestack)
        results.update(wave_props)

        return results

    def _compute_spatial_correlation(self, timestack: np.ndarray) -> np.ndarray:
        """
        Compute spatial correlation for each time step.

        Parameters
        ----------
        timestack : np.ndarray
            2D timestack array

        Returns
        -------
        np.ndarray
            Spatial correlation array
        """
        num_positions = timestack.shape[0]
        num_times = timestack.shape[1]

        # Select positions at regular spacing
        positions = np.arange(0, num_positions, self.correlation_spacing)

        correlations = []

        for t in range(num_times):
            profile = timestack[:, t]

            # Compute autocorrelation
            autocorr = signal.correlate(profile, profile, mode='full')
            autocorr = autocorr[len(autocorr)//2:]  # Keep positive lags only

            correlations.append(autocorr)

        return np.array(correlations)

    def _compute_temporal_correlation(self, timestack: np.ndarray) -> np.ndarray:
        """
        Compute temporal correlation for each spatial position.

        Parameters
        ----------
        timestack : np.ndarray
            2D timestack array

        Returns
        -------
        np.ndarray
            Temporal correlation array
        """
        num_positions = timestack.shape[0]

        correlations = []

        for i in range(num_positions):
            timeseries = timestack[i, :]

            # Compute autocorrelation
            autocorr = signal.correlate(timeseries, timeseries, mode='full')
            autocorr = autocorr[len(autocorr)//2:]  # Keep positive lags only

            correlations.append(autocorr)

        return np.array(correlations)

    def _extract_wave_properties(self, timestack: np.ndarray) -> Dict:
        """
        Extract wave celerity and wavelength from timestack.

        Uses MATLAB-style algorithm: fixed time lag, variable spatial lag.

        Parameters
        ----------
        timestack : np.ndarray
            2D timestack array (space x time)

        Returns
        -------
        Dict
            Wave properties (celerity, wavelength, period)
        """
        results = {}

        # Use MATLAB-style cross-correlation (fixed time lag, variable spatial lag)
        celerities = self._compute_celerity_matlab_style(timestack)

        # Store results
        results['celerities'] = celerities

        # Compute mean values
        valid_celerities = celerities[~np.isnan(celerities)]
        if len(valid_celerities) > 0:
            results['mean_celerity'] = np.mean(valid_celerities)
            results['std_celerity'] = np.std(valid_celerities)
        else:
            results['mean_celerity'] = np.nan
            results['std_celerity'] = np.nan

        # Note: wavelengths and periods computed elsewhere
        results['wavelengths'] = np.array([])
        results['periods'] = np.array([])
        results['mean_wavelength'] = np.nan
        results['std_wavelength'] = np.nan
        results['mean_period'] = np.nan
        results['std_period'] = np.nan

        return results

    def _compute_celerity_matlab_style(self, timestack: np.ndarray) -> np.ndarray:
        """
        Compute wave celerity using MATLAB's CrossCorrelation_CoastCams algorithm.

        Algorithm:
        1. Use FIXED time lag (dpha = 1 second, converted to frames)
        2. For each spatial position, compute correlation at different spatial lags
        3. Find spatial lag that maximizes correlation
        4. Convert to velocity: spatial_lag / time_lag

        Parameters
        ----------
        timestack : np.ndarray
            2D timestack array (space x time) - transposed from MATLAB!
            MATLAB: A2(time, space)
            Python: timestack(space, time)

        Returns
        -------
        np.ndarray
            Array of celerities at each spatial position (in m/s)
        """
        # Input timestack is already (time x space) format from main.py
        # Example: (1680, 689) = (time points, spatial positions)
        A2 = timestack
        print(f"  Timestack shape for correlation: {A2.shape} (time x space)")

        # Detrend (remove linear trend from each column)
        from scipy.signal import detrend
        A2 = detrend(A2, axis=0, type='linear')

        nt, nc = A2.shape  # nt = time points, nc = spatial points
        print(f"  After detrend: nt={nt} time points, nc={nc} spatial points")

        # Parameters (matching MATLAB)
        dpha = 1.0  # Time lag in seconds
        dc = self.correlation_spacing  # Spatial window (100 pixels)
        dc = int(np.round(dc / 2) * 2)  # Ensure even

        n = int(np.floor(dpha / self.dt))  # Time lag in frames

        # Validate
        if n <= 0 or n >= nt:
            print(f"  Warning: Invalid time lag n={n} (nt={nt}, dt={self.dt})")
            return np.full(nc, np.nan)

        if nc < dc:
            print(f"  Warning: Spatial dimension {nc} < dc={dc}")
            return np.full(nc, np.nan)

        # Preallocate correlation matrix
        # R2M will store correlation profiles for each spatial position
        num_positions = nc - dc + 1
        R2M = np.zeros((num_positions, dc - 1))

        # For each spatial position
        pos_idx = 0
        for ic in range(int(dc/2), int(nc - dc/2)):
            R2 = np.zeros(dc - 1)

            # For each spatial lag
            for lc in range(1, dc):
                pos1 = ic - int(dc/2)
                pos2 = ic - int(dc/2) + lc - 1

                # Validate indices
                if pos1 >= 0 and pos2 < nc and pos1 < nc and pos2 >= 0:
                    # Get time series at two positions with time lag
                    # Position 1: later times (n+1:nt)
                    # Position 2: earlier times (1:nt-n)
                    ts1 = A2[n:nt, pos1]  # Python: 0-indexed, so n:nt is equivalent to n+1:nt in MATLAB
                    ts2 = A2[0:nt-n, pos2]

                    # Compute correlation coefficient
                    if len(ts1) > 1 and len(ts2) > 1 and np.std(ts1) > 0 and np.std(ts2) > 0:
                        corr = np.corrcoef(ts1, ts2)[0, 1]
                        R2[lc - 1] = corr
                    else:
                        R2[lc - 1] = np.nan
                else:
                    R2[lc - 1] = np.nan

            if pos_idx < num_positions:
                R2M[pos_idx, :] = R2
                pos_idx += 1

        # Debug: show sample correlation profile
        if R2M.shape[0] > 0:
            sample_profile = R2M[R2M.shape[0]//2, :]  # Middle position
            print(f"  Sample correlation profile (middle position):")
            print(f"    First 10 lags: {sample_profile[:10]}")
            print(f"    Last 10 lags: {sample_profile[-10:]}")

        # Find spatial lag with maximum correlation for each position
        spatial_lags = np.argmax(R2M, axis=1) + 1  # +1 because lc starts at 1

        # Debug: check correlation values
        max_corr = np.max(R2M, axis=1)
        mean_max_corr = np.mean(max_corr)
        print(f"  Cross-correlation debug:")
        print(f"    R2M shape: {R2M.shape}")
        print(f"    Max correlation values: mean={mean_max_corr:.3f}, range=[{np.min(max_corr):.3f}, {np.max(max_corr):.3f}]")
        print(f"    Spatial lags: mean={np.mean(spatial_lags):.1f}, range=[{np.min(spatial_lags)}, {np.max(spatial_lags)}]")
        print(f"    Time lag n={n} frames ({n*self.dt:.1f}s), dpha={dpha}s")

        # Convert to velocity: spatial_lag / time_lag
        # spatial_lag is in pixels, dpha is in seconds
        # Result is in pixels/second
        velocities_pixels_per_sec = spatial_lags / dpha

        # Convert from pixels/second to m/s
        # Method 1: Multiply by pixel resolution
        # Method 2 (MATLAB): Divide by 10 (equivalent when res=0.1)
        velocities_ms = velocities_pixels_per_sec * self.dx

        # Apply moving average (MATLAB line 242: movmean(..., 10))
        from scipy.ndimage import uniform_filter1d
        window = 10
        velocities_smooth = uniform_filter1d(velocities_ms, size=window, mode='nearest')

        # Interpolate to full spatial grid (R2M covers only central positions)
        # Pad with NaN for positions outside the correlation window
        celerities_full = np.full(nc, np.nan)
        start_idx = int(dc/2)
        end_idx = start_idx + len(velocities_smooth)
        if end_idx > nc:
            end_idx = nc
            velocities_smooth = velocities_smooth[:nc - start_idx]
        celerities_full[start_idx:end_idx] = velocities_smooth

        return celerities_full

    def compute_phase_velocity(self, timestack: np.ndarray,
                               position1: int, position2: int) -> float:
        """
        Compute phase velocity between two spatial positions.

        Parameters
        ----------
        timestack : np.ndarray
            2D timestack array
        position1 : int
            First spatial position (pixels)
        position2 : int
            Second spatial position (pixels)

        Returns
        -------
        float
            Phase velocity (m/s)
        """
        if position1 >= timestack.shape[0] or position2 >= timestack.shape[0]:
            return np.nan

        # Get time series at both positions
        ts1 = timestack[position1, :]
        ts2 = timestack[position2, :]

        # Compute cross-correlation
        cc = signal.correlate(ts2, ts1, mode='full')
        lags = signal.correlation_lags(len(ts2), len(ts1), mode='full')

        # Find peak
        peak_idx = np.argmax(cc)
        time_lag = lags[peak_idx] * self.dt

        if time_lag <= 0:
            return np.nan

        # Compute velocity
        distance = abs(position2 - position1) * self.dx
        velocity = distance / time_lag

        return velocity

    def compute_group_velocity(self, timestack: np.ndarray) -> Dict:
        """
        Compute wave group velocity from timestack.

        Parameters
        ----------
        timestack : np.ndarray
            2D timestack array

        Returns
        -------
        Dict
            Group velocity parameters
        """
        # Compute envelope of signal
        from scipy.signal import hilbert

        envelopes = []

        for i in range(timestack.shape[0]):
            ts = timestack[i, :]
            analytic = hilbert(ts)
            envelope = np.abs(analytic)
            envelopes.append(envelope)

        envelope_stack = np.array(envelopes)

        # Compute phase velocity on envelope
        group_velocities = []

        num_positions = envelope_stack.shape[0]
        positions = np.arange(0, num_positions - self.correlation_spacing,
                            self.correlation_spacing)

        for pos in positions:
            env1 = envelope_stack[pos, :]
            env2 = envelope_stack[pos + self.correlation_spacing, :]

            # Cross-correlate envelopes
            cc = signal.correlate(env2, env1, mode='full')
            lags = signal.correlation_lags(len(env2), len(env1), mode='full')

            peak_idx = np.argmax(cc)
            time_lag = lags[peak_idx] * self.dt

            if time_lag > 0:
                distance = self.correlation_spacing * self.dx
                velocity = distance / time_lag
                group_velocities.append(velocity)

        return {
            'group_velocities': np.array(group_velocities),
            'mean_group_velocity': np.mean(group_velocities) if len(group_velocities) > 0 else np.nan
        }

    def compute_2d_correlation_matrix(self, timestack: np.ndarray,
                                     max_lag: Optional[int] = None) -> np.ndarray:
        """
        Compute 2D space-time correlation matrix.

        Parameters
        ----------
        timestack : np.ndarray
            2D timestack array
        max_lag : int, optional
            Maximum lag to compute (default: None, use all)

        Returns
        -------
        np.ndarray
            2D correlation matrix
        """
        if max_lag is None:
            max_lag = min(timestack.shape[0], timestack.shape[1]) // 4

        # Normalize timestack
        ts_norm = (timestack - np.mean(timestack)) / (np.std(timestack) + 1e-10)

        # Compute 2D correlation using FFT
        from scipy.fft import fft2, ifft2

        fft_ts = fft2(ts_norm)
        correlation = np.real(ifft2(fft_ts * np.conj(fft_ts)))

        # Shift to center
        correlation = np.fft.fftshift(correlation)

        # Extract region around zero lag
        center_y, center_x = np.array(correlation.shape) // 2
        correlation_roi = correlation[
            center_y - max_lag:center_y + max_lag,
            center_x - max_lag:center_x + max_lag
        ]

        return correlation_roi

    def estimate_dispersion_relation(self, timestack: np.ndarray) -> Dict:
        """
        Estimate dispersion relation from timestack using 2D FFT.

        Parameters
        ----------
        timestack : np.ndarray
            2D timestack array

        Returns
        -------
        Dict
            Frequency-wavenumber spectrum and dispersion relation
        """
        from scipy.fft import fft2, fftfreq

        # Compute 2D FFT
        fft_result = fft2(timestack)
        spectrum = np.abs(fft_result) ** 2

        # Get frequency and wavenumber axes
        freq = fftfreq(timestack.shape[1], d=self.dt)
        wavenumber = fftfreq(timestack.shape[0], d=self.dx)

        # Shift to center
        spectrum_shifted = np.fft.fftshift(spectrum)
        freq_shifted = np.fft.fftshift(freq)
        k_shifted = np.fft.fftshift(wavenumber)

        return {
            'spectrum': spectrum_shifted,
            'frequencies': freq_shifted,
            'wavenumbers': k_shifted
        }
