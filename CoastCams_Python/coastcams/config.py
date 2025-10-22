"""
Configuration management for CoastCams.

This module handles loading and managing configuration parameters
for CoastCams analysis, with smart defaults for non-expert users.
"""

import os
import yaml
from typing import Dict, Any, Optional
from pathlib import Path


class CoastCamsConfig:
    """
    Configuration manager for CoastCams analysis.

    This class provides a user-friendly interface for configuring
    CoastCams analysis with sensible defaults for non-experts.
    """

    def __init__(self, config_file: Optional[str] = None):
        """
        Initialize configuration with defaults or from file.

        Parameters
        ----------
        config_file : str, optional
            Path to YAML configuration file. If None, uses defaults.
        """
        # Set default configuration
        self._set_defaults()

        # Load from file if provided
        if config_file and os.path.exists(config_file):
            self.load_from_file(config_file)

    def _set_defaults(self):
        """Set default configuration parameters for typical coastal monitoring."""

        # === Camera Parameters ===
        self.camera_height = 27.240  # Camera height above mean sea level (meters)
        self.pixel_resolution = 0.1  # Pixel resolution (meters/pixel)
        self.rotation_angle = 270  # Image rotation angle (degrees)
        self.acquisition_frequency = 2.0  # Images per second (Hz)

        # === Image Processing Parameters ===
        self.cross_shore_width = 1600  # Number of pixels for cross-shore analysis
        self.min_cross_shore = 1  # Minimum pixel for analysis
        self.max_cross_shore = None  # Maximum pixel (None = auto)

        # === Shoreline Detection ===
        self.shoreline_method = 1  # 1=Grayscale, 2=Red-Blue, 3=Color convergence
        self.shoreline_threshold = 30  # Intensity threshold for detection

        # === Wave Analysis Parameters ===
        self.time_lag = 1  # Time lag for phase correlation (frames)
        self.correlation_spacing = 100  # Spacing between correlation windows (pixels)
        self.wave_period_min = 4.0  # Minimum expected wave period (seconds)
        self.wave_period_max = 25.0  # Maximum expected wave period (seconds)

        # === Physical Constants ===
        self.gravity = 9.81  # Gravitational acceleration (m/s^2)
        self.water_density = 1025  # Seawater density (kg/m^3)

        # === Filtering Parameters ===
        self.lowpass_cutoff = 0.05  # Low-pass filter cutoff frequency (Hz)
        self.highpass_cutoff = 0.005  # High-pass filter cutoff frequency (Hz)
        self.filter_order = 5  # Butterworth filter order
        self.smoothing_window = 3  # 2D smoothing window size

        # === Output Parameters ===
        self.output_interval = 15  # Resampling interval for output (minutes)
        self.save_plots = True  # Save visualization plots
        self.save_intermediate = False  # Save intermediate results
        self.output_format = 'csv'  # Output format: 'csv', 'excel', or 'both'

        # === File Paths (auto-detected if not specified) ===
        self.input_dir = 'Timestacks'  # Input directory for timestack images
        self.output_dir = 'Output'  # Output directory for results
        self.coordinate_file = None  # Optional coordinate transformation file

        # === Visualization Options ===
        self.plot_shoreline = True
        self.plot_waves = True
        self.plot_bathymetry = True
        self.plot_timeseries = True

        # === Advanced Options ===
        self.use_coordinate_transform = False  # Apply coordinate transformation
        self.parallel_processing = False  # Use multiprocessing (future feature)
        self.verbose = True  # Print progress messages

    def load_from_file(self, config_file: str):
        """
        Load configuration from YAML file.

        Parameters
        ----------
        config_file : str
            Path to YAML configuration file
        """
        with open(config_file, 'r') as f:
            config_data = yaml.safe_load(f)

        # Update configuration with file values
        for key, value in config_data.items():
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                print(f"Warning: Unknown configuration parameter '{key}' ignored")

    def save_to_file(self, config_file: str):
        """
        Save current configuration to YAML file.

        Parameters
        ----------
        config_file : str
            Path to save configuration file
        """
        config_dict = self.to_dict()

        with open(config_file, 'w') as f:
            yaml.dump(config_dict, f, default_flow_style=False, sort_keys=False)

        print(f"Configuration saved to {config_file}")

    def to_dict(self) -> Dict[str, Any]:
        """
        Convert configuration to dictionary.

        Returns
        -------
        Dict[str, Any]
            Configuration as dictionary
        """
        # Get all non-private attributes
        config_dict = {
            key: value for key, value in self.__dict__.items()
            if not key.startswith('_')
        }

        return config_dict

    def validate(self) -> bool:
        """
        Validate configuration parameters.

        Returns
        -------
        bool
            True if configuration is valid, False otherwise
        """
        errors = []

        # Validate camera parameters
        if self.camera_height <= 0:
            errors.append("Camera height must be positive")

        if self.pixel_resolution <= 0:
            errors.append("Pixel resolution must be positive")

        # Validate wave parameters
        if self.wave_period_min >= self.wave_period_max:
            errors.append("Minimum wave period must be less than maximum")

        # Validate shoreline method
        if self.shoreline_method not in [1, 2, 3]:
            errors.append("Shoreline method must be 1, 2, or 3")

        # Validate file paths
        if not os.path.exists(self.input_dir):
            errors.append(f"Input directory '{self.input_dir}' does not exist")

        # Print errors if any
        if errors:
            print("Configuration validation errors:")
            for error in errors:
                print(f"  - {error}")
            return False

        return True

    def auto_detect_paths(self, base_dir: Optional[str] = None):
        """
        Auto-detect input and output directories.

        Parameters
        ----------
        base_dir : str, optional
            Base directory to search in. If None, uses current directory.
        """
        if base_dir is None:
            base_dir = os.getcwd()

        # Look for Timestacks directory
        timestacks_dir = os.path.join(base_dir, 'Timestacks')
        if os.path.exists(timestacks_dir):
            self.input_dir = timestacks_dir

        # Create Output directory if it doesn't exist
        output_dir = os.path.join(base_dir, 'Output')
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        self.output_dir = output_dir

        # Look for coordinate file
        coord_file = os.path.join(base_dir, 'UserScripts', 'Coordinates.txt')
        if os.path.exists(coord_file):
            self.coordinate_file = coord_file
            self.use_coordinate_transform = True

    def print_summary(self):
        """Print a summary of current configuration."""
        print("\n" + "="*60)
        print("CoastCams Configuration Summary")
        print("="*60)

        print("\nCamera Parameters:")
        print(f"  Height above MSL: {self.camera_height:.3f} m")
        print(f"  Pixel resolution: {self.pixel_resolution:.3f} m/pixel")
        print(f"  Rotation angle: {self.rotation_angle}°")
        print(f"  Acquisition frequency: {self.acquisition_frequency} Hz")

        print("\nWave Analysis:")
        print(f"  Period range: {self.wave_period_min}-{self.wave_period_max} s")
        print(f"  Correlation spacing: {self.correlation_spacing} pixels")
        print(f"  Cross-shore width: {self.cross_shore_width} pixels")

        print("\nShoreline Detection:")
        methods = {1: "Grayscale", 2: "Red-Blue", 3: "Color Convergence"}
        print(f"  Method: {methods[self.shoreline_method]}")
        print(f"  Threshold: {self.shoreline_threshold}")

        print("\nOutput Settings:")
        print(f"  Input directory: {self.input_dir}")
        print(f"  Output directory: {self.output_dir}")
        print(f"  Output interval: {self.output_interval} minutes")
        print(f"  Output format: {self.output_format}")
        print(f"  Save plots: {self.save_plots}")

        print("="*60 + "\n")

    @classmethod
    def create_template(cls, output_file: str = "config_template.yaml"):
        """
        Create a template configuration file.

        Parameters
        ----------
        output_file : str, optional
            Path to save template file (default: 'config_template.yaml')
        """
        config = cls()
        config.save_to_file(output_file)
        print(f"Template configuration created at {output_file}")
        print("Edit this file to customize your analysis settings.")
