"""
Bathymetry estimation for CoastCams.

This module estimates water depth profiles using linear wave theory
and wave dispersion relationships.
"""

import numpy as np
from scipy import optimize
from typing import Dict, Tuple, Optional, List
from .utils import linear_wave_celerity, calculate_depth_from_celerity


class BathymetryEstimator:
    """
    Estimate bathymetry (water depth) from wave parameters.

    Uses linear wave theory and dispersion relations to
    estimate depth profiles from observed wave characteristics.
    """

    def __init__(self, config=None):
        """
        Initialize bathymetry estimator.

        Parameters
        ----------
        config : CoastCamsConfig, optional
            Configuration object
        """
        if config is not None:
            self.gravity = config.gravity
            self.pixel_resolution = config.pixel_resolution
        else:
            # Default parameters
            self.gravity = 9.81  # m/s^2
            self.pixel_resolution = 0.1  # meters/pixel

    def estimate_depth_from_waves(self, wave_period: float,
                                  wave_celerity: float,
                                  method: str = 'linear') -> float:
        """
        Estimate water depth from wave period and celerity.

        Parameters
        ----------
        wave_period : float
            Wave period (seconds)
        wave_celerity : float
            Wave celerity (m/s)
        method : str, optional
            Estimation method: 'linear' or 'shallow' (default: 'linear')

        Returns
        -------
        float
            Estimated water depth (meters)
        """
        if wave_period <= 0 or wave_celerity <= 0:
            return np.nan

        if method == 'linear':
            # Use full linear dispersion relation
            depth = calculate_depth_from_celerity(wave_celerity, wave_period,
                                                 self.gravity)
        elif method == 'shallow':
            # Use shallow water approximation: C = sqrt(g*h)
            depth = (wave_celerity ** 2) / self.gravity
        else:
            raise ValueError(f"Unknown method: {method}")

        return depth

    def estimate_depth_profile(self, wave_periods: np.ndarray,
                              wave_celerities: np.ndarray,
                              cross_shore_positions: np.ndarray,
                              method: str = 'linear') -> Dict:
        """
        Estimate depth profile along cross-shore transect.

        Parameters
        ----------
        wave_periods : np.ndarray
            Wave periods at each position
        wave_celerities : np.ndarray
            Wave celerities at each position
        cross_shore_positions : np.ndarray
            Cross-shore positions (pixels or meters)
        method : str, optional
            Estimation method (default: 'linear')

        Returns
        -------
        Dict
            Depth profile and statistics
        """
        depths = np.zeros(len(wave_periods))

        for i in range(len(wave_periods)):
            if not np.isnan(wave_periods[i]) and not np.isnan(wave_celerities[i]):
                depths[i] = self.estimate_depth_from_waves(
                    wave_periods[i], wave_celerities[i], method
                )
            else:
                depths[i] = np.nan

        # Filter outliers and smooth
        depths_filtered = self._filter_depths(depths)
        depths_smoothed = self._smooth_depths(depths_filtered)

        results = {
            'depths': depths,
            'depths_filtered': depths_filtered,
            'depths_smoothed': depths_smoothed,
            'cross_shore_positions': cross_shore_positions,
            'mean_depth': np.nanmean(depths_filtered),
            'max_depth': np.nanmax(depths_filtered),
            'min_depth': np.nanmin(depths_filtered)
        }

        return results

    def _filter_depths(self, depths: np.ndarray,
                      max_depth: float = 50.0,
                      min_depth: float = 0.1) -> np.ndarray:
        """
        Filter unrealistic depth values.

        Parameters
        ----------
        depths : np.ndarray
            Raw depth estimates
        max_depth : float, optional
            Maximum realistic depth (default: 50 m)
        min_depth : float, optional
            Minimum realistic depth (default: 0.1 m)

        Returns
        -------
        np.ndarray
            Filtered depths
        """
        filtered = depths.copy()

        # Remove unrealistic values
        filtered[filtered < min_depth] = np.nan
        filtered[depths > max_depth] = np.nan

        # Remove statistical outliers
        valid_depths = filtered[~np.isnan(filtered)]

        if len(valid_depths) > 0:
            mean_depth = np.mean(valid_depths)
            std_depth = np.std(valid_depths)

            outlier_mask = np.abs(filtered - mean_depth) > (3 * std_depth)
            filtered[outlier_mask] = np.nan

        return filtered

    def _smooth_depths(self, depths: np.ndarray,
                      window_size: int = 5) -> np.ndarray:
        """
        Smooth depth profile using moving average.

        Parameters
        ----------
        depths : np.ndarray
            Depth values
        window_size : int, optional
            Smoothing window size (default: 5)

        Returns
        -------
        np.ndarray
            Smoothed depths
        """
        from scipy.ndimage import uniform_filter1d

        # Handle NaN values
        valid_mask = ~np.isnan(depths)

        if not np.any(valid_mask):
            return depths

        # Interpolate NaNs for smoothing
        x = np.arange(len(depths))
        valid_x = x[valid_mask]
        valid_depths = depths[valid_mask]

        if len(valid_x) < 2:
            return depths

        # Interpolate
        interpolated = np.interp(x, valid_x, valid_depths)

        # Smooth
        smoothed = uniform_filter1d(interpolated, size=window_size)

        # Restore NaN where original was NaN
        smoothed[~valid_mask] = np.nan

        return smoothed

    def compute_depth_gradient(self, depths: np.ndarray,
                              dx: float) -> np.ndarray:
        """
        Compute depth gradient (slope).

        Parameters
        ----------
        depths : np.ndarray
            Depth values
        dx : float
            Spatial spacing (meters)

        Returns
        -------
        np.ndarray
            Depth gradient (dimensionless)
        """
        gradient = np.gradient(depths, dx)

        return gradient

    def estimate_breaking_depth(self, wave_height: float,
                               gamma: float = 0.78) -> float:
        """
        Estimate depth at wave breaking location.

        Uses breaker index gamma = H/h at breaking.

        Parameters
        ----------
        wave_height : float
            Breaking wave height (meters)
        gamma : float, optional
            Breaker index (default: 0.78)

        Returns
        -------
        float
            Breaking depth (meters)
        """
        return wave_height / gamma

    def refine_depth_with_breaking(self, initial_depths: np.ndarray,
                                   wave_heights: np.ndarray,
                                   break_location: int,
                                   gamma: float = 0.78) -> np.ndarray:
        """
        Refine depth estimates using breaking wave constraint.

        Parameters
        ----------
        initial_depths : np.ndarray
            Initial depth estimates
        wave_heights : np.ndarray
            Wave heights
        break_location : int
            Index of breaking location
        gamma : float, optional
            Breaker index (default: 0.78)

        Returns
        -------
        np.ndarray
            Refined depth estimates
        """
        refined = initial_depths.copy()

        if break_location < 0 or break_location >= len(wave_heights):
            return refined

        # Estimate depth at breaking
        breaking_wave_height = wave_heights[break_location]
        breaking_depth = self.estimate_breaking_depth(breaking_wave_height, gamma)

        # Adjust depths to match breaking constraint
        if not np.isnan(initial_depths[break_location]):
            depth_offset = breaking_depth - initial_depths[break_location]
            refined = initial_depths + depth_offset

        return refined

    def compute_shore_slope(self, depths: np.ndarray,
                           cross_shore_positions: np.ndarray,
                           region: str = 'nearshore') -> float:
        """
        Compute shore slope in specified region.

        Parameters
        ----------
        depths : np.ndarray
            Depth profile
        cross_shore_positions : np.ndarray
            Cross-shore positions (meters)
        region : str, optional
            Region to compute slope: 'nearshore', 'offshore', or 'all'
            (default: 'nearshore')

        Returns
        -------
        float
            Shore slope (dimensionless)
        """
        # Filter valid data
        valid_mask = ~np.isnan(depths)
        valid_depths = depths[valid_mask]
        valid_positions = cross_shore_positions[valid_mask]

        if len(valid_depths) < 2:
            return np.nan

        # Select region
        if region == 'nearshore':
            # Use shallow region (depth < 5m)
            region_mask = valid_depths < 5.0
        elif region == 'offshore':
            # Use deeper region (depth > 5m)
            region_mask = valid_depths > 5.0
        else:
            # Use all data
            region_mask = np.ones(len(valid_depths), dtype=bool)

        if not np.any(region_mask):
            return np.nan

        # Fit linear slope
        from scipy.stats import linregress

        slope, intercept, r_value, p_value, std_err = linregress(
            valid_positions[region_mask],
            valid_depths[region_mask]
        )

        return slope

    def estimate_depth_from_shoaling(self, offshore_wave_height: float,
                                    local_wave_height: float,
                                    offshore_depth: float) -> float:
        """
        Estimate local depth from wave shoaling.

        Parameters
        ----------
        offshore_wave_height : float
            Wave height in deep water (meters)
        local_wave_height : float
            Local wave height (meters)
        offshore_depth : float
            Offshore water depth (meters)

        Returns
        -------
        float
            Estimated local depth (meters)
        """
        # Shoaling coefficient K_s = H_local / H_offshore = sqrt(C_offshore / C_local)
        # For linear waves: K_s = sqrt(tanh(kh_offshore) / tanh(kh_local))

        if offshore_wave_height <= 0 or local_wave_height <= 0:
            return np.nan

        K_s = local_wave_height / offshore_wave_height

        # This is a simplified approach - full solution requires iteration
        # Using shallow water approximation: K_s ~ (h_offshore / h_local)^0.25

        local_depth = offshore_depth / (K_s ** 4)

        return local_depth

    def compute_depth_statistics(self, depths: np.ndarray) -> Dict:
        """
        Compute statistical properties of depth profile.

        Parameters
        ----------
        depths : np.ndarray
            Depth values

        Returns
        -------
        Dict
            Depth statistics
        """
        valid_depths = depths[~np.isnan(depths)]

        if len(valid_depths) == 0:
            return {
                'mean': np.nan,
                'std': np.nan,
                'min': np.nan,
                'max': np.nan,
                'median': np.nan,
                'range': np.nan
            }

        return {
            'mean': np.mean(valid_depths),
            'std': np.std(valid_depths),
            'min': np.min(valid_depths),
            'max': np.max(valid_depths),
            'median': np.median(valid_depths),
            'range': np.max(valid_depths) - np.min(valid_depths)
        }
