"""
CoastCams - Python Coastal Wave Analysis Toolbox

A comprehensive Python port of the MATLAB CoastCams toolbox for analyzing
coastal wave dynamics and nearshore morphology from land-based video monitoring.

Author: Ported to Python
License: GNU General Public License v3
Original: Nuyts et al. (2023), Environmental Modelling & Software, 168, 105800
"""

__version__ = "1.0.0"
__author__ = "CoastCams Development Team"

from .config import CoastCamsConfig
from .image_loader import ImageLoader
from .preprocessing import ImagePreprocessor
from .shoreline import ShorelineDetector
from .wave_analysis import WaveAnalyzer
from .bathymetry import BathymetryEstimator
from .visualize import CoastCamsVisualizer

__all__ = [
    'CoastCamsConfig',
    'ImageLoader',
    'ImagePreprocessor',
    'ShorelineDetector',
    'WaveAnalyzer',
    'BathymetryEstimator',
    'CoastCamsVisualizer',
]
