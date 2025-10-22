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

        Parameters
        ----------
        timestack : np.ndarray
            2D timestack array

        Returns
        -------
        Dict
            Wave properties (celerity, wavelength, period)
        """
        results = {}

        # Compute phase velocity using cross-correlation at different lags
        celerities = []
        wavelengths = []
        periods = []

        num_positions = timestack.shape[0]
        positions = np.arange(0, num_positions - self.correlation_spacing,
                            self.correlation_spacing)

        for pos in positions:
            # Get two spatial locations
            ts1 = timestack[pos, :]
            ts2 = timestack[pos + self.correlation_spacing, :]

            # Compute cross-correlation
            cc = signal.correlate(ts2, ts1, mode='full')
            lags = signal.correlation_lags(len(ts2), len(ts1), mode='full')

            # Find peak correlation
            peak_idx = np.argmax(cc)
            time_lag = lags[peak_idx] * self.dt

            if time_lag > 0:
                # Compute celerity
                distance = self.correlation_spacing * self.dx
                celerity = distance / time_lag
                celerities.append(celerity)

                # Estimate period from autocorrelation
                autocorr = signal.correlate(ts1, ts1, mode='full')
                autocorr = autocorr[len(autocorr)//2:]

                # Find first peak after zero lag
                if len(autocorr) > 10:
                    peaks = signal.find_peaks(autocorr)[0]
                    if len(peaks) > 0:
                        period = peaks[0] * self.dt
                        periods.append(period)

                        # Compute wavelength
                        wavelength = celerity * period
                        wavelengths.append(wavelength)

        # Store results
        results['celerities'] = np.array(celerities)
        results['wavelengths'] = np.array(wavelengths)
        results['periods'] = np.array(periods)

        # Compute mean values
        if len(celerities) > 0:
            results['mean_celerity'] = np.mean(celerities)
            results['std_celerity'] = np.std(celerities)
        else:
            results['mean_celerity'] = np.nan
            results['std_celerity'] = np.nan

        if len(wavelengths) > 0:
            results['mean_wavelength'] = np.mean(wavelengths)
            results['std_wavelength'] = np.std(wavelengths)
        else:
            results['mean_wavelength'] = np.nan
            results['std_wavelength'] = np.nan

        if len(periods) > 0:
            results['mean_period'] = np.mean(periods)
            results['std_period'] = np.std(periods)
        else:
            results['mean_period'] = np.nan
            results['std_period'] = np.nan

        return results

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
