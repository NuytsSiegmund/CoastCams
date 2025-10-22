"""
Sea level anomaly calculations for CoastCams.

This module computes sea level anomalies from shoreline positions
and bathymetry data.
"""

import numpy as np
from typing import Dict, Optional, List
from datetime import datetime


class SeaLevelAnalyzer:
    """
    Compute sea level anomalies from coastal observations.

    Estimates deviations from mean sea level using shoreline positions,
    wave runup, and bathymetry data.
    """

    def __init__(self, config=None):
        """
        Initialize sea level analyzer.

        Parameters
        ----------
        config : CoastCamsConfig, optional
            Configuration object
        """
        if config is not None:
            self.pixel_resolution = config.pixel_resolution
            self.camera_height = config.camera_height
        else:
            # Default parameters
            self.pixel_resolution = 0.1  # meters/pixel
            self.camera_height = 27.240  # meters above MSL

        # Mean sea level reference (can be calibrated)
        self.mean_sea_level = 0.0

    def compute_sla_from_shoreline(self, shoreline_positions: np.ndarray,
                                   depths: np.ndarray,
                                   method: str = 'linear') -> np.ndarray:
        """
        Compute sea level anomaly from shoreline positions and depth.

        Parameters
        ----------
        shoreline_positions : np.ndarray
            Shoreline positions in pixels or meters
        depths : np.ndarray
            Corresponding water depths
        method : str, optional
            Calculation method: 'linear' or 'shallow' (default: 'linear')

        Returns
        -------
        np.ndarray
            Sea level anomalies (meters)
        """
        # Convert positions to meters if needed
        if np.nanmean(shoreline_positions) > 100:
            # Assume pixels
            positions_m = shoreline_positions * self.pixel_resolution
        else:
            positions_m = shoreline_positions

        # Compute mean shoreline position
        mean_position = np.nanmean(positions_m)

        # Compute anomalies
        sla = np.zeros(len(shoreline_positions))

        for i in range(len(shoreline_positions)):
            if not np.isnan(positions_m[i]) and not np.isnan(depths[i]):
                # Deviation from mean position
                position_deviation = positions_m[i] - mean_position

                # Estimate SLA using depth and position change
                # Simple approach: SLA ~ depth_change
                sla[i] = position_deviation * np.tan(np.arctan(0.02))  # Assume 2% slope
            else:
                sla[i] = np.nan

        return sla

    def compute_sla_from_depth(self, observed_depths: np.ndarray,
                              reference_depths: np.ndarray) -> np.ndarray:
        """
        Compute sea level anomaly from depth measurements.

        Parameters
        ----------
        observed_depths : np.ndarray
            Observed water depths
        reference_depths : np.ndarray
            Reference (mean) depths

        Returns
        -------
        np.ndarray
            Sea level anomalies (meters)
        """
        sla = observed_depths - reference_depths

        return sla

    def compute_runup(self, wave_height: float, wave_period: float,
                     beach_slope: float) -> float:
        """
        Estimate wave runup using empirical formula.

        Parameters
        ----------
        wave_height : float
            Significant wave height (meters)
        wave_period : float
            Peak wave period (seconds)
        beach_slope : float
            Beach slope (dimensionless)

        Returns
        -------
        float
            Wave runup (meters)
        """
        if wave_height <= 0 or wave_period <= 0 or beach_slope <= 0:
            return np.nan

        # Iribarren number (surf similarity parameter)
        wavelength_deep = 1.56 * wave_period ** 2  # Deep water wavelength
        iribarren = beach_slope / np.sqrt(wave_height / wavelength_deep)

        # Stockdon et al. (2006) formula
        # R_2% = 1.1 * (setup + 0.5 * sqrt(S_inc^2 + S_ig^2))
        setup = 0.35 * beach_slope * np.sqrt(wave_height * wavelength_deep)
        S_inc = 0.75 * beach_slope * np.sqrt(wave_height * wavelength_deep)
        S_ig = 0.06 * np.sqrt(wave_height * wavelength_deep)

        runup = 1.1 * (setup + 0.5 * np.sqrt(S_inc**2 + S_ig**2))

        return runup

    def compute_tidal_range(self, shoreline_positions: List[float],
                           timestamps: List[datetime],
                           return_components: bool = False) -> Dict:
        """
        Estimate tidal range from shoreline variations.

        Parameters
        ----------
        shoreline_positions : List[float]
            Time series of shoreline positions
        timestamps : List[datetime]
            Corresponding timestamps
        return_components : bool, optional
            Return tidal components (default: False)

        Returns
        -------
        Dict
            Tidal statistics and optionally components
        """
        positions = np.array(shoreline_positions)

        # Remove NaN values
        valid_mask = ~np.isnan(positions)
        valid_positions = positions[valid_mask]

        if len(valid_positions) < 2:
            return {'range': np.nan, 'mean': np.nan}

        # Compute tidal range as difference between max and min
        tidal_range = np.max(valid_positions) - np.min(valid_positions)
        mean_position = np.mean(valid_positions)

        results = {
            'range': tidal_range,
            'mean': mean_position,
            'max': np.max(valid_positions),
            'min': np.min(valid_positions)
        }

        if return_components:
            # Simple harmonic analysis (would need more sophisticated method for real tides)
            # This is a placeholder for tidal component extraction
            results['components'] = self._extract_tidal_components(
                valid_positions, np.array(timestamps)[valid_mask]
            )

        return results

    def _extract_tidal_components(self, positions: np.ndarray,
                                  timestamps: np.ndarray) -> Dict:
        """
        Extract tidal harmonic components (simplified).

        Parameters
        ----------
        positions : np.ndarray
            Shoreline positions
        timestamps : np.ndarray
            Timestamps

        Returns
        -------
        Dict
            Tidal components (M2, S2, etc.)
        """
        # This is a simplified placeholder
        # Real tidal analysis would use harmonic analysis (e.g., UTide)

        from scipy.signal import periodogram

        # Convert timestamps to hours
        if len(timestamps) > 0:
            time_hours = np.array([(t - timestamps[0]).total_seconds() / 3600
                                  for t in timestamps])

            # Compute periodogram
            freqs, power = periodogram(positions, fs=1.0/np.mean(np.diff(time_hours)))

            # Find dominant periods
            dominant_idx = np.argsort(power)[-3:]  # Top 3 components
            dominant_periods = 1.0 / freqs[dominant_idx]

            return {
                'dominant_periods_hours': dominant_periods,
                'frequencies': freqs[dominant_idx],
                'power': power[dominant_idx]
            }

        return {}

    def compute_relative_tidal_range(self, shoreline_positions: np.ndarray,
                                    wave_heights: np.ndarray) -> np.ndarray:
        """
        Compute relative tidal range (RTR = tidal_range / wave_height).

        Parameters
        ----------
        shoreline_positions : np.ndarray
            Shoreline positions over time
        wave_heights : np.ndarray
            Significant wave heights

        Returns
        -------
        np.ndarray
            Relative tidal range values
        """
        # Estimate tidal component from low-pass filtered shoreline
        from scipy.signal import butter, filtfilt

        # Remove short-period variations (< 1 hour period)
        # Assuming positions are at regular intervals
        cutoff_freq = 1.0 / 3600  # 1 hour cutoff

        # This is simplified - actual implementation would need proper time series
        if len(shoreline_positions) > 10:
            b, a = butter(3, 0.1)  # Low-pass filter
            tidal_component = filtfilt(b, a, shoreline_positions)

            tidal_range = np.max(tidal_component) - np.min(tidal_component)
        else:
            tidal_range = np.max(shoreline_positions) - np.min(shoreline_positions)

        # Compute RTR
        mean_wave_height = np.nanmean(wave_heights)

        if mean_wave_height > 0:
            rtr = tidal_range / mean_wave_height
        else:
            rtr = np.nan

        return np.full(len(shoreline_positions), rtr)

    def calibrate_mean_sea_level(self, reference_positions: np.ndarray,
                                 reference_depths: np.ndarray):
        """
        Calibrate mean sea level from reference measurements.

        Parameters
        ----------
        reference_positions : np.ndarray
            Reference shoreline positions (meters)
        reference_depths : np.ndarray
            Reference water depths (meters)
        """
        # Compute mean values
        valid_mask = ~np.isnan(reference_positions) & ~np.isnan(reference_depths)

        if np.any(valid_mask):
            self.mean_sea_level = np.mean(reference_depths[valid_mask])

    def compute_storm_surge(self, observed_sea_level: np.ndarray,
                           predicted_tide: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Compute storm surge as residual from predicted tide.

        Parameters
        ----------
        observed_sea_level : np.ndarray
            Observed sea level values
        predicted_tide : np.ndarray, optional
            Predicted tidal values (if None, uses low-pass filtered signal)

        Returns
        -------
        np.ndarray
            Storm surge residuals
        """
        if predicted_tide is None:
            # Estimate tide as low-frequency component
            from scipy.signal import butter, filtfilt

            if len(observed_sea_level) > 10:
                b, a = butter(3, 0.1)
                predicted_tide = filtfilt(b, a, observed_sea_level)
            else:
                # Not enough data
                return np.zeros(len(observed_sea_level))

        # Storm surge is residual
        storm_surge = observed_sea_level - predicted_tide

        return storm_surge

    def compute_sla_statistics(self, sla_values: np.ndarray,
                              timestamps: Optional[List[datetime]] = None) -> Dict:
        """
        Compute statistical properties of sea level anomalies.

        Parameters
        ----------
        sla_values : np.ndarray
            Sea level anomaly values
        timestamps : List[datetime], optional
            Timestamps for time-based statistics

        Returns
        -------
        Dict
            SLA statistics
        """
        valid_sla = sla_values[~np.isnan(sla_values)]

        if len(valid_sla) == 0:
            return {
                'mean': np.nan,
                'std': np.nan,
                'min': np.nan,
                'max': np.nan,
                'range': np.nan
            }

        stats = {
            'mean': np.mean(valid_sla),
            'std': np.std(valid_sla),
            'min': np.min(valid_sla),
            'max': np.max(valid_sla),
            'range': np.max(valid_sla) - np.min(valid_sla),
            'median': np.median(valid_sla)
        }

        if timestamps is not None and len(timestamps) == len(sla_values):
            # Add time-based statistics
            stats['duration_hours'] = (timestamps[-1] - timestamps[0]).total_seconds() / 3600

        return stats
