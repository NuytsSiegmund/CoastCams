"""
Shoreline detection for CoastCams analysis.

This module provides multiple methods for detecting the shoreline
position in timestack images.
"""

import numpy as np
import cv2
from scipy import ndimage, signal
from typing import Tuple, Optional, List
from .utils import smooth2d


class ShorelineDetector:
    """
    Detect shoreline positions using multiple methods.

    Provides three detection methods:
    1. Grayscale intensity method (for rocky platforms)
    2. Red-minus-blue channel method
    3. Color convergence method
    """

    def __init__(self, method: int = 1, threshold: int = 30, config=None):
        """
        Initialize shoreline detector.

        Parameters
        ----------
        method : int, optional
            Detection method (1, 2, or 3) (default: 1)
        threshold : int, optional
            Intensity threshold for detection (default: 30)
        config : CoastCamsConfig, optional
            Configuration object
        """
        if config is not None:
            self.method = config.shoreline_method
            self.threshold = config.shoreline_threshold
        else:
            self.method = method
            self.threshold = threshold

    def detect(self, img: np.ndarray) -> Optional[np.ndarray]:
        """
        Detect shoreline position using configured method.

        Parameters
        ----------
        img : np.ndarray
            Input image (RGB or grayscale)

        Returns
        -------
        Optional[np.ndarray]
            Array of shoreline positions (one per column) or None if detection fails
        """
        if self.method == 1:
            return self.detect_grayscale(img)
        elif self.method == 2:
            return self.detect_red_minus_blue(img)
        elif self.method == 3:
            return self.detect_color_convergence(img)
        else:
            raise ValueError(f"Unknown detection method: {self.method}")

    def detect_grayscale(self, img: np.ndarray) -> np.ndarray:
        """
        Detect shoreline using grayscale intensity method.

        Best for rocky platforms with clear intensity contrast.

        Parameters
        ----------
        img : np.ndarray
            Input image (RGB or grayscale)

        Returns
        -------
        np.ndarray
            Shoreline position for each column (in pixels)
        """
        # Convert to grayscale if needed
        if len(img.shape) == 3:
            gray = cv2.cvtColor(img.astype(np.uint8), cv2.COLOR_RGB2GRAY)
        else:
            gray = img.astype(np.uint8)

        # Compute gradient magnitude
        gradient_y = np.abs(np.diff(gray.astype(np.float32), axis=0))

        # Find maximum gradient position in each column
        shoreline_positions = np.argmax(gradient_y, axis=0)

        # Apply smoothing to remove outliers
        shoreline_smooth = self._smooth_shoreline(shoreline_positions)

        # Convert to float to allow NaN values
        shoreline_smooth = shoreline_smooth.astype(np.float64)

        # Filter by threshold
        max_gradients = np.max(gradient_y, axis=0)
        shoreline_smooth[max_gradients < self.threshold] = np.nan

        return shoreline_smooth

    def detect_red_minus_blue(self, img: np.ndarray) -> np.ndarray:
        """
        Detect shoreline using red-minus-blue channel method.

        Uses color difference between water (blue) and land (red/brown).

        Parameters
        ----------
        img : np.ndarray
            Input RGB image

        Returns
        -------
        np.ndarray
            Shoreline position for each column (in pixels)
        """
        if len(img.shape) != 3 or img.shape[2] != 3:
            raise ValueError("Red-minus-blue method requires RGB image")

        # Extract red and blue channels
        red = img[:, :, 0].astype(np.float32)
        blue = img[:, :, 2].astype(np.float32)

        # Compute red-minus-blue difference
        rmb = red - blue

        # Compute gradient
        gradient_y = np.abs(np.diff(rmb, axis=0))

        # Find maximum gradient position in each column
        shoreline_positions = np.argmax(gradient_y, axis=0)

        # Apply smoothing
        shoreline_smooth = self._smooth_shoreline(shoreline_positions)

        # Convert to float to allow NaN values
        shoreline_smooth = shoreline_smooth.astype(np.float64)

        # Filter by threshold
        max_gradients = np.max(gradient_y, axis=0)
        shoreline_smooth[max_gradients < self.threshold] = np.nan

        return shoreline_smooth

    def detect_color_convergence(self, img: np.ndarray) -> np.ndarray:
        """
        Detect shoreline using color convergence method.

        Analyzes convergence of all color channels.

        Parameters
        ----------
        img : np.ndarray
            Input RGB image

        Returns
        -------
        np.ndarray
            Shoreline position for each column (in pixels)
        """
        if len(img.shape) != 3 or img.shape[2] != 3:
            raise ValueError("Color convergence method requires RGB image")

        # Compute gradients for each channel
        gradients = []
        for i in range(3):
            channel = img[:, :, i].astype(np.float32)
            gradient_y = np.abs(np.diff(channel, axis=0))
            gradients.append(gradient_y)

        # Combine gradients
        combined_gradient = np.mean(gradients, axis=0)

        # Find maximum gradient position in each column
        shoreline_positions = np.argmax(combined_gradient, axis=0)

        # Apply smoothing
        shoreline_smooth = self._smooth_shoreline(shoreline_positions)

        # Convert to float to allow NaN values
        shoreline_smooth = shoreline_smooth.astype(np.float64)

        # Filter by threshold
        max_gradients = np.max(combined_gradient, axis=0)
        shoreline_smooth[max_gradients < self.threshold] = np.nan

        return shoreline_smooth

    def _smooth_shoreline(self, positions: np.ndarray,
                         window_size: int = 15) -> np.ndarray:
        """
        Smooth shoreline positions using moving median filter.

        Parameters
        ----------
        positions : np.ndarray
            Raw shoreline positions
        window_size : int, optional
            Size of smoothing window (default: 15)

        Returns
        -------
        np.ndarray
            Smoothed shoreline positions
        """
        # Apply median filter
        smoothed = ndimage.median_filter(positions, size=window_size)

        return smoothed

    def detect_multiple_images(self, images: List[np.ndarray]) -> List[np.ndarray]:
        """
        Detect shoreline in multiple images.

        Parameters
        ----------
        images : List[np.ndarray]
            List of images

        Returns
        -------
        List[np.ndarray]
            List of shoreline position arrays
        """
        shorelines = []

        for img in images:
            shoreline = self.detect(img)
            shorelines.append(shoreline)

        return shorelines

    def get_average_position(self, shorelines: List[np.ndarray]) -> float:
        """
        Compute average shoreline position.

        Parameters
        ----------
        shorelines : List[np.ndarray]
            List of shoreline position arrays

        Returns
        -------
        float
            Average shoreline position (pixels)
        """
        all_positions = []

        for shoreline in shorelines:
            valid_positions = shoreline[~np.isnan(shoreline)]
            if len(valid_positions) > 0:
                all_positions.extend(valid_positions)

        if len(all_positions) == 0:
            return 0.0

        return np.mean(all_positions)

    def get_shoreline_variation(self, shorelines: List[np.ndarray]) -> Tuple[float, float]:
        """
        Compute shoreline variation statistics.

        Parameters
        ----------
        shorelines : List[np.ndarray]
            List of shoreline position arrays

        Returns
        -------
        Tuple[float, float]
            Mean and standard deviation of shoreline positions
        """
        all_positions = []

        for shoreline in shorelines:
            valid_positions = shoreline[~np.isnan(shoreline)]
            if len(valid_positions) > 0:
                all_positions.extend(valid_positions)

        if len(all_positions) == 0:
            return 0.0, 0.0

        return np.mean(all_positions), np.std(all_positions)

    def pixels_to_meters(self, pixel_position: float,
                        pixel_resolution: float) -> float:
        """
        Convert shoreline position from pixels to meters.

        Parameters
        ----------
        pixel_position : float
            Position in pixels
        pixel_resolution : float
            Resolution in meters/pixel

        Returns
        -------
        float
            Position in meters
        """
        return pixel_position * pixel_resolution

    def filter_outliers(self, shoreline: np.ndarray,
                       std_threshold: float = 3.0) -> np.ndarray:
        """
        Remove outliers from shoreline positions.

        Parameters
        ----------
        shoreline : np.ndarray
            Shoreline position array
        std_threshold : float, optional
            Number of standard deviations for outlier threshold (default: 3.0)

        Returns
        -------
        np.ndarray
            Filtered shoreline positions
        """
        filtered = shoreline.copy()

        # Get valid positions
        valid_mask = ~np.isnan(filtered)
        valid_positions = filtered[valid_mask]

        if len(valid_positions) == 0:
            return filtered

        # Compute statistics
        mean_pos = np.mean(valid_positions)
        std_pos = np.std(valid_positions)

        # Mark outliers as NaN
        outlier_mask = np.abs(filtered - mean_pos) > (std_threshold * std_pos)
        filtered[outlier_mask] = np.nan

        return filtered

    def interpolate_gaps(self, shoreline: np.ndarray) -> np.ndarray:
        """
        Interpolate missing shoreline positions.

        Parameters
        ----------
        shoreline : np.ndarray
            Shoreline position array with possible NaN values

        Returns
        -------
        np.ndarray
            Interpolated shoreline positions
        """
        interpolated = shoreline.copy()

        # Find NaN positions
        nan_mask = np.isnan(interpolated)

        if not np.any(nan_mask):
            return interpolated

        # Get indices
        x = np.arange(len(interpolated))
        valid_x = x[~nan_mask]
        valid_y = interpolated[~nan_mask]

        if len(valid_x) < 2:
            # Not enough points for interpolation
            return interpolated

        # Interpolate
        interpolated[nan_mask] = np.interp(x[nan_mask], valid_x, valid_y)

        return interpolated
