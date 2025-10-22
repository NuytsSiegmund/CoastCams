"""
Image loading and management for CoastCams.

This module handles loading timestack images, extracting metadata,
and preparing images for analysis.
"""

import os
import re
import cv2
import numpy as np
from datetime import datetime
from typing import List, Tuple, Optional, Dict
from pathlib import Path
from tqdm import tqdm


class ImageLoader:
    """
    Load and manage timestack images for CoastCams analysis.

    This class handles automatic image discovery, metadata extraction,
    and image preprocessing (rotation, format conversion).
    """

    def __init__(self, input_dir: str, rotation_angle: int = 270, verbose: bool = True):
        """
        Initialize image loader.

        Parameters
        ----------
        input_dir : str
            Directory containing timestack images
        rotation_angle : int, optional
            Angle to rotate images (default: 270)
        verbose : bool, optional
            Print progress messages (default: True)
        """
        self.input_dir = input_dir
        self.rotation_angle = rotation_angle
        self.verbose = verbose

        self.image_files = []
        self.timestamps = []
        self.metadata = []

    def discover_images(self, pattern: str = "*.jpeg") -> int:
        """
        Automatically discover timestack images in input directory.

        Parameters
        ----------
        pattern : str, optional
            File pattern to match (default: '*.jpeg')

        Returns
        -------
        int
            Number of images found
        """
        # Find all matching image files
        image_paths = sorted(Path(self.input_dir).glob(pattern))

        self.image_files = []
        self.timestamps = []
        self.metadata = []

        if self.verbose:
            print(f"Discovering images in {self.input_dir}...")

        for img_path in image_paths:
            filename = img_path.name

            # Extract timestamp from filename
            timestamp = self._extract_timestamp(filename)

            if timestamp:
                self.image_files.append(str(img_path))
                self.timestamps.append(timestamp)

                # Store metadata
                self.metadata.append({
                    'filename': filename,
                    'path': str(img_path),
                    'timestamp': timestamp,
                    'datetime': timestamp
                })

        if self.verbose:
            print(f"Found {len(self.image_files)} images")
            if len(self.image_files) > 0:
                print(f"Time range: {self.timestamps[0]} to {self.timestamps[-1]}")

        return len(self.image_files)

    def _extract_timestamp(self, filename: str) -> Optional[datetime]:
        """
        Extract timestamp from filename.

        Expected format: S_X_YYYYMMDDHHSS.jpeg or similar

        Parameters
        ----------
        filename : str
            Image filename

        Returns
        -------
        Optional[datetime]
            Extracted datetime or None if parsing fails
        """
        # Try to extract 12-digit timestamp (YYYYMMDDHHSS)
        match = re.search(r'(\d{12})', filename)

        if match:
            timestamp_str = match.group(1)
            try:
                # Parse as YYYYMMDDHHSS
                dt = datetime.strptime(timestamp_str, '%Y%m%d%H%M')
                return dt
            except ValueError:
                pass

        # Try alternative format: YYYYMMDD_HHSS
        match = re.search(r'(\d{8})_(\d{4})', filename)
        if match:
            date_str = match.group(1)
            time_str = match.group(2)
            timestamp_str = date_str + time_str
            try:
                dt = datetime.strptime(timestamp_str, '%Y%m%d%H%M')
                return dt
            except ValueError:
                pass

        if self.verbose:
            print(f"Warning: Could not extract timestamp from {filename}")

        return None

    def load_image(self, index: int, as_grayscale: bool = False) -> Optional[np.ndarray]:
        """
        Load a single image by index.

        Parameters
        ----------
        index : int
            Index of image to load
        as_grayscale : bool, optional
            Load as grayscale image (default: False)

        Returns
        -------
        Optional[np.ndarray]
            Loaded image array or None if loading fails
        """
        if index < 0 or index >= len(self.image_files):
            print(f"Error: Index {index} out of range")
            return None

        img_path = self.image_files[index]

        # Load image
        if as_grayscale:
            img = cv2.imread(img_path, cv2.IMREAD_GRAYSCALE)
        else:
            img = cv2.imread(img_path, cv2.IMREAD_COLOR)
            # Convert BGR to RGB
            img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)

        if img is None:
            print(f"Error: Failed to load {img_path}")
            return None

        # Apply rotation if specified
        if self.rotation_angle != 0:
            img = self._rotate_image(img, self.rotation_angle)

        return img

    def load_all_images(self, as_grayscale: bool = False,
                       max_images: Optional[int] = None) -> List[np.ndarray]:
        """
        Load all discovered images.

        Parameters
        ----------
        as_grayscale : bool, optional
            Load as grayscale images (default: False)
        max_images : int, optional
            Maximum number of images to load (default: None, load all)

        Returns
        -------
        List[np.ndarray]
            List of loaded image arrays
        """
        images = []

        num_to_load = len(self.image_files)
        if max_images is not None:
            num_to_load = min(num_to_load, max_images)

        if self.verbose:
            print(f"Loading {num_to_load} images...")

        iterator = range(num_to_load)
        if self.verbose:
            iterator = tqdm(iterator, desc="Loading images")

        for i in iterator:
            img = self.load_image(i, as_grayscale=as_grayscale)
            if img is not None:
                images.append(img)

        return images

    def _rotate_image(self, img: np.ndarray, angle: int) -> np.ndarray:
        """
        Rotate image by specified angle.

        Parameters
        ----------
        img : np.ndarray
            Input image
        angle : int
            Rotation angle in degrees (90, 180, 270)

        Returns
        -------
        np.ndarray
            Rotated image
        """
        if angle == 90:
            return cv2.rotate(img, cv2.ROTATE_90_CLOCKWISE)
        elif angle == 180:
            return cv2.rotate(img, cv2.ROTATE_180)
        elif angle == 270:
            return cv2.rotate(img, cv2.ROTATE_90_COUNTERCLOCKWISE)
        else:
            # Use general rotation for arbitrary angles
            height, width = img.shape[:2]
            center = (width // 2, height // 2)
            rotation_matrix = cv2.getRotationMatrix2D(center, angle, 1.0)
            return cv2.warpAffine(img, rotation_matrix, (width, height))

    def get_image_info(self, index: int) -> Dict:
        """
        Get metadata for a specific image.

        Parameters
        ----------
        index : int
            Image index

        Returns
        -------
        Dict
            Image metadata
        """
        if index < 0 or index >= len(self.metadata):
            return {}

        return self.metadata[index]

    def get_time_range(self) -> Tuple[Optional[datetime], Optional[datetime]]:
        """
        Get time range of loaded images.

        Returns
        -------
        Tuple[Optional[datetime], Optional[datetime]]
            Start and end datetime
        """
        if len(self.timestamps) == 0:
            return None, None

        return self.timestamps[0], self.timestamps[-1]

    def filter_by_time_range(self, start_time: datetime,
                            end_time: datetime) -> 'ImageLoader':
        """
        Filter images by time range.

        Parameters
        ----------
        start_time : datetime
            Start of time range
        end_time : datetime
            End of time range

        Returns
        -------
        ImageLoader
            New ImageLoader with filtered images
        """
        filtered_loader = ImageLoader(self.input_dir, self.rotation_angle, self.verbose)

        for i, timestamp in enumerate(self.timestamps):
            if start_time <= timestamp <= end_time:
                filtered_loader.image_files.append(self.image_files[i])
                filtered_loader.timestamps.append(self.timestamps[i])
                filtered_loader.metadata.append(self.metadata[i])

        if self.verbose:
            print(f"Filtered to {len(filtered_loader.image_files)} images")

        return filtered_loader

    def get_summary(self) -> str:
        """
        Get summary of loaded images.

        Returns
        -------
        str
            Summary string
        """
        if len(self.image_files) == 0:
            return "No images loaded"

        summary = f"ImageLoader Summary:\n"
        summary += f"  Total images: {len(self.image_files)}\n"
        summary += f"  Directory: {self.input_dir}\n"

        if len(self.timestamps) > 0:
            start_time, end_time = self.get_time_range()
            duration = end_time - start_time
            summary += f"  Time range: {start_time} to {end_time}\n"
            summary += f"  Duration: {duration}\n"

            # Calculate average time interval
            if len(self.timestamps) > 1:
                intervals = [(self.timestamps[i+1] - self.timestamps[i]).total_seconds()
                           for i in range(len(self.timestamps)-1)]
                avg_interval = np.mean(intervals)
                summary += f"  Average interval: {avg_interval:.1f} seconds\n"

        return summary
