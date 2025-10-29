"""
Image preprocessing for CoastCams analysis.

This module provides filtering and denoising operations
for timestack images before wave analysis.
"""

import numpy as np
import cv2
from scipy import signal, ndimage
from typing import Tuple, Optional
from .utils import smooth2d, detrend_2d, butter_lowpass_filter, butter_highpass_filter


class ImagePreprocessor:
    """
    Preprocess timestack images for wave analysis.

    Applies filtering, smoothing, and detrending operations
    to reduce noise and enhance wave signals.
    """

    def __init__(self, config=None):
        """
        Initialize preprocessor.

        Parameters
        ----------
        config : CoastCamsConfig, optional
            Configuration object with processing parameters
        """
        if config is not None:
            self.lowpass_cutoff = config.lowpass_cutoff
            self.highpass_cutoff = config.highpass_cutoff
            self.filter_order = config.filter_order
            self.smoothing_window = config.smoothing_window
            self.acquisition_frequency = config.acquisition_frequency
        else:
            # Default parameters
            self.lowpass_cutoff = 0.05  # Hz
            self.highpass_cutoff = 0.005  # Hz
            self.filter_order = 5
            self.smoothing_window = 3
            self.acquisition_frequency = 2.0  # Hz

    def preprocess_image(self, img: np.ndarray,
                        apply_smoothing: bool = True,
                        apply_detrend: bool = True) -> np.ndarray:
        """
        Apply full preprocessing pipeline to image.

        Parameters
        ----------
        img : np.ndarray
            Input image (can be grayscale or color)
        apply_smoothing : bool, optional
            Apply spatial smoothing (default: True)
        apply_detrend : bool, optional
            Apply detrending (default: True)

        Returns
        -------
        np.ndarray
            Preprocessed image
        """
        processed = img.copy()

        # Convert to float for processing
        if processed.dtype == np.uint8:
            processed = processed.astype(np.float32)

        # Apply smoothing if requested
        if apply_smoothing:
            processed = self.smooth_image(processed)

        # Apply detrending if requested
        if apply_detrend:
            processed = self.detrend_image(processed)

        return processed

    def smooth_image(self, img: np.ndarray) -> np.ndarray:
        """
        Apply 2D spatial smoothing to image.

        Parameters
        ----------
        img : np.ndarray
            Input image

        Returns
        -------
        np.ndarray
            Smoothed image
        """
        if len(img.shape) == 3:
            # Process each channel separately
            smoothed = np.zeros_like(img)
            for i in range(img.shape[2]):
                smoothed[:, :, i] = smooth2d(img[:, :, i], self.smoothing_window)
        else:
            # Grayscale image
            smoothed = smooth2d(img, self.smoothing_window)

        return smoothed

    def detrend_image(self, img: np.ndarray, axis: int = 1) -> np.ndarray:
        """
        Remove linear trend from image along specified axis.

        Parameters
        ----------
        img : np.ndarray
            Input image
        axis : int, optional
            Axis along which to detrend (default: 1, horizontal)

        Returns
        -------
        np.ndarray
            Detrended image
        """
        if len(img.shape) == 3:
            # Process each channel separately
            detrended = np.zeros_like(img)
            for i in range(img.shape[2]):
                detrended[:, :, i] = detrend_2d(img[:, :, i], axis=axis)
        else:
            # Grayscale image
            detrended = detrend_2d(img, axis=axis)

        return detrended

    def filter_temporal(self, timeseries: np.ndarray,
                       filter_type: str = 'both') -> np.ndarray:
        """
        Apply temporal filtering to time series data.

        Parameters
        ----------
        timeseries : np.ndarray
            1D or 2D time series data (time axis should be last)
        filter_type : str, optional
            Type of filter: 'low', 'high', or 'both' (default: 'both')

        Returns
        -------
        np.ndarray
            Filtered time series
        """
        filtered = timeseries.copy()

        if filter_type in ['low', 'both']:
            # Apply low-pass filter
            if len(timeseries.shape) == 1:
                filtered = butter_lowpass_filter(
                    filtered, self.lowpass_cutoff,
                    self.acquisition_frequency, self.filter_order
                )
            else:
                # Apply along each row
                for i in range(filtered.shape[0]):
                    filtered[i, :] = butter_lowpass_filter(
                        filtered[i, :], self.lowpass_cutoff,
                        self.acquisition_frequency, self.filter_order
                    )

        if filter_type in ['high', 'both']:
            # Apply high-pass filter
            if len(timeseries.shape) == 1:
                filtered = butter_highpass_filter(
                    filtered, self.highpass_cutoff,
                    self.acquisition_frequency, self.filter_order
                )
            else:
                # Apply along each row
                for i in range(filtered.shape[0]):
                    filtered[i, :] = butter_highpass_filter(
                        filtered[i, :], self.highpass_cutoff,
                        self.acquisition_frequency, self.filter_order
                    )

        return filtered

    def extract_blue_channel(self, img: np.ndarray,
                            normalize: bool = True) -> np.ndarray:
        """
        Extract and optionally normalize blue channel.

        Blue channel often provides best contrast for water surface.

        Parameters
        ----------
        img : np.ndarray
            RGB image
        normalize : bool, optional
            Normalize to [0, 1] range (default: True)

        Returns
        -------
        np.ndarray
            Blue channel
        """
        if len(img.shape) != 3 or img.shape[2] != 3:
            raise ValueError("Image must be RGB (3 channels)")

        blue = img[:, :, 2].astype(np.float32)

        if normalize and blue.max() > 1.0:
            blue = blue / 255.0

        return blue

    def extract_red_minus_blue(self, img: np.ndarray,
                               normalize: bool = True) -> np.ndarray:
        """
        Extract red-minus-blue channel.

        Useful for shoreline detection and water/land contrast.

        Parameters
        ----------
        img : np.ndarray
            RGB image
        normalize : bool, optional
            Normalize to [0, 1] range (default: True)

        Returns
        -------
        np.ndarray
            Red-minus-blue difference
        """
        if len(img.shape) != 3 or img.shape[2] != 3:
            raise ValueError("Image must be RGB (3 channels)")

        red = img[:, :, 0].astype(np.float32)
        blue = img[:, :, 2].astype(np.float32)

        if normalize and red.max() > 1.0:
            red = red / 255.0
            blue = blue / 255.0

        red_minus_blue = red - blue

        return red_minus_blue

    def create_timestack_array(self, images: list,
                              cross_shore_range: Tuple[int, int]) -> np.ndarray:
        """
        Create 2D timestack array from list of images.

        Parameters
        ----------
        images : list
            List of image arrays
        cross_shore_range : Tuple[int, int]
            Range of cross-shore pixels to extract (min, max)

        Returns
        -------
        np.ndarray
            2D timestack array (space x time)
        """
        min_pixel, max_pixel = cross_shore_range

        timestacks = []

        for img in images:
            # Extract cross-shore slice
            if len(img.shape) == 3:
                # Use blue channel for color images
                blue = self.extract_blue_channel(img, normalize=True)
                slice_data = np.mean(blue[:, min_pixel:max_pixel], axis=0)
            else:
                # Grayscale image
                slice_data = np.mean(img[:, min_pixel:max_pixel], axis=0)

            timestacks.append(slice_data)

        # Stack into 2D array
        timestack_array = np.array(timestacks).T

        return timestack_array

    def enhance_waves(self, timestack: np.ndarray) -> np.ndarray:
        """
        Enhance wave signals in timestack data.

        Applies detrending and filtering to emphasize wave patterns.

        Parameters
        ----------
        timestack : np.ndarray
            2D timestack array (space x time)

        Returns
        -------
        np.ndarray
            Enhanced timestack
        """
        enhanced = timestack.copy()

        # Remove spatial mean
        enhanced = enhanced - np.mean(enhanced, axis=1, keepdims=True)

        # Remove temporal mean
        enhanced = enhanced - np.mean(enhanced, axis=0, keepdims=True)

        # Apply temporal filtering
        enhanced = self.filter_temporal(enhanced, filter_type='both')

        # Apply spatial smoothing
        enhanced = smooth2d(enhanced, window_size=self.smoothing_window)

        return enhanced

    def normalize_image(self, img: np.ndarray,
                       method: str = 'minmax') -> np.ndarray:
        """
        Normalize image intensity values.

        Parameters
        ----------
        img : np.ndarray
            Input image
        method : str, optional
            Normalization method: 'minmax' or 'zscore' (default: 'minmax')

        Returns
        -------
        np.ndarray
            Normalized image
        """
        if method == 'minmax':
            # Normalize to [0, 1]
            img_min = np.min(img)
            img_max = np.max(img)
            if img_max - img_min > 0:
                normalized = (img - img_min) / (img_max - img_min)
            else:
                normalized = img - img_min
        elif method == 'zscore':
            # Z-score normalization
            mean = np.mean(img)
            std = np.std(img)
            if std > 0:
                normalized = (img - mean) / std
            else:
                normalized = img - mean
        else:
            raise ValueError(f"Unknown normalization method: {method}")

        return normalized

    def remove_mean(self, data: np.ndarray, axis: Optional[int] = None) -> np.ndarray:
        """
        Remove mean from data along specified axis.

        Parameters
        ----------
        data : np.ndarray
            Input data
        axis : int, optional
            Axis along which to compute mean (default: None, remove global mean)

        Returns
        -------
        np.ndarray
            Data with mean removed
        """
        if axis is None:
            return data - np.mean(data)
        else:
            return data - np.mean(data, axis=axis, keepdims=True)
