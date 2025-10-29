# CoastCams - Python Edition

**Automated Coastal Wave Analysis Toolbox**

A complete Python port of the MATLAB CoastCams toolbox for analyzing coastal wave dynamics, bathymetry, and nearshore morphology from land-based video monitoring systems.

---

## Overview

CoastCams-Python is a fully automated, modular, and user-friendly tool designed for non-experts to analyze coastal processes. It processes oblique orthorectified timestack images to extract:

- **Wave parameters**: height, period, celerity, energy
- **Bathymetry**: underwater depth profiles
- **Shoreline positions**: water-land boundary detection
- **Sea level anomalies**: deviations from mean sea level
- **Morphological changes**: temporal evolution of coastal features

---

## Key Features

✅ **Fully Automated** - Minimal user input required
✅ **Modular Design** - Run complete workflow or specific tasks
✅ **Non-Expert Friendly** - Smart defaults and clear documentation
✅ **Multiple Analysis Methods** - 3 shoreline detection methods, 2 bathymetry approaches
✅ **Comprehensive Visualization** - Automatic generation of publication-quality plots
✅ **Flexible Configuration** - YAML-based or programmatic configuration

---

## Installation

### Prerequisites

- Python 3.7 or higher
- pip package manager

### Quick Install

```bash
# Clone the repository
cd CoastCams

# Install dependencies
pip install -r requirements.txt
```

### Manual Dependency Installation

```bash
pip install numpy opencv-python scipy scikit-image pandas matplotlib seaborn PyYAML tqdm pywavelets
```

---

## Quick Start

### 1. Run Complete Analysis (Easiest)

```bash
# Run with default settings
python main.py

# Run with custom configuration
python main.py --config my_config.yaml

# Create a template configuration file
python main.py --create-config
```

### 2. Run Specific Tasks

CoastCams provides example scripts for focused analyses:

```bash
# Shoreline detection only
python examples/run_shoreline_only.py

# Wave analysis only
python examples/run_wave_analysis_only.py

# Bathymetry estimation only
python examples/run_bathymetry_only.py

# Full workflow
python examples/run_full_workflow.py
```

---

## Directory Structure

```
CoastCams/
├── coastcams/                  # Main Python package
│   ├── __init__.py            # Package initialization
│   ├── config.py              # Configuration management
│   ├── image_loader.py        # Image I/O and metadata
│   ├── preprocessing.py       # Image filtering and denoising
│   ├── shoreline.py           # Shoreline detection (3 methods)
│   ├── wave_analysis.py       # Wave parameter extraction
│   ├── cross_correlation.py   # Cross-correlation analysis
│   ├── bathymetry.py          # Depth estimation
│   ├── sea_level.py           # Sea level anomaly calculation
│   ├── visualize.py           # Plotting and visualization
│   └── utils.py               # Helper functions
├── examples/                   # Example scripts
│   ├── run_full_workflow.py
│   ├── run_shoreline_only.py
│   ├── run_wave_analysis_only.py
│   └── run_bathymetry_only.py
├── Timestacks/                 # Input directory (timestack images)
├── Output/                     # Output directory (results and plots)
├── main.py                     # Main workflow script
├── requirements.txt            # Python dependencies
├── config.yaml                 # Default configuration file
└── README_PYTHON.md           # This file
```

---

## Configuration

### Creating a Configuration File

```bash
python main.py --create-config
```

This creates `config_template.yaml` with all available options.

### Configuration Options

```yaml
# Camera Parameters
camera_height: 27.240          # Height above mean sea level (m)
pixel_resolution: 0.1          # Meters per pixel
rotation_angle: 270            # Image rotation (degrees)
acquisition_frequency: 2.0     # Images per second (Hz)

# Shoreline Detection
shoreline_method: 1            # 1=Grayscale, 2=Red-Blue, 3=Color convergence
shoreline_threshold: 30        # Intensity threshold

# Wave Analysis
wave_period_min: 4.0           # Minimum wave period (s)
wave_period_max: 25.0          # Maximum wave period (s)
correlation_spacing: 100       # Correlation window spacing (pixels)

# Output Settings
output_interval: 15            # Resampling interval (minutes)
save_plots: true               # Save visualization plots
output_format: 'csv'           # 'csv', 'excel', or 'both'

# File Paths
input_dir: 'Timestacks'        # Input image directory
output_dir: 'Output'           # Output directory
```

### Programmatic Configuration

```python
from coastcams.config import CoastCamsConfig

# Create config
config = CoastCamsConfig()

# Customize settings
config.shoreline_method = 2  # Use red-minus-blue method
config.wave_period_min = 5.0
config.wave_period_max = 20.0
config.save_plots = True

# Save to file
config.save_to_file('my_config.yaml')
```

---

## Usage Examples

### Example 1: Full Automated Analysis

```python
from main import CoastCamsWorkflow

# Create workflow with default settings
workflow = CoastCamsWorkflow()

# Run complete analysis
workflow.run_full_analysis()
```

### Example 2: Shoreline Detection Only

```python
from coastcams import ImageLoader, ShorelineDetector, CoastCamsConfig

# Initialize
config = CoastCamsConfig()
loader = ImageLoader('Timestacks', rotation_angle=270)
detector = ShorelineDetector(method=1, threshold=30)

# Load image
loader.discover_images()
img = loader.load_image(0)

# Detect shoreline
shoreline = detector.detect(img)

# Get statistics
mean_pos, std_pos = detector.get_shoreline_variation([shoreline])
print(f"Mean shoreline position: {mean_pos:.1f} pixels")
```

### Example 3: Wave Parameter Analysis

```python
from coastcams import (ImageLoader, ImagePreprocessor,
                       WaveAnalyzer, CoastCamsConfig)
import numpy as np

# Initialize
config = CoastCamsConfig()
loader = ImageLoader('Timestacks')
preprocessor = ImagePreprocessor(config)
analyzer = WaveAnalyzer(config)

# Load and process
loader.discover_images()
images = loader.load_all_images()

# Create timestack
timestack = preprocessor.create_timestack_array(
    images, cross_shore_range=(0, 1600)
)

# Analyze waves
positions = np.arange(timestack.shape[0]) * config.pixel_resolution
results = analyzer.analyze_timestack(timestack, positions)

print(f"Mean wave height: {results['mean_Hs']:.2f} m")
print(f"Mean wave period: {results['mean_Tm']:.2f} s")
```

### Example 4: Bathymetry Estimation

```python
from coastcams import BathymetryEstimator, CoastCamsConfig
import numpy as np

# Initialize
config = CoastCamsConfig()
estimator = BathymetryEstimator(config)

# Example wave parameters
wave_periods = np.array([8.0, 9.0, 10.0, 11.0])  # seconds
wave_celerities = np.array([10.0, 9.5, 9.0, 8.5])  # m/s
positions = np.array([0, 50, 100, 150])  # meters

# Estimate depth profile
bathymetry = estimator.estimate_depth_profile(
    wave_periods, wave_celerities, positions, method='linear'
)

print(f"Mean depth: {bathymetry['mean_depth']:.2f} m")
print(f"Depth range: {bathymetry['min_depth']:.2f} - {bathymetry['max_depth']:.2f} m")
```

---

## Analysis Workflow

The complete CoastCams workflow consists of 7 steps:

1. **Image Loading** - Automatically discover and load timestack images
2. **Shoreline Detection** - Detect water-land boundary using selected method
3. **Image Preprocessing** - Filter, denoise, and create timestack arrays
4. **Wave Analysis** - Extract wave height, period, and energy
5. **Cross-Correlation** - Compute wave celerity and wavelength
6. **Bathymetry Estimation** - Estimate water depth profiles
7. **Sea Level Analysis** - Compute sea level anomalies

---

## Shoreline Detection Methods

### Method 1: Grayscale Intensity
- Best for: Rocky platforms with clear intensity contrast
- Uses: Gradient of grayscale intensity

### Method 2: Red-Minus-Blue
- Best for: Sandy beaches with water-land color difference
- Uses: Difference between red and blue color channels

### Method 3: Color Convergence
- Best for: Mixed environments
- Uses: Combined gradient from all color channels

---

## Output Files

### CSV/Excel Data Files
- `coastcams_results_YYYYMMDD_HHMMSS.csv` - Complete analysis results
- `analysis_summary.txt` - Text summary report

### Visualization Plots
- `timestack.png` - Enhanced timestack image
- `shoreline_method_X.png` - Shoreline detection results
- `wave_parameters.png` - Time series of wave parameters
- `wave_profiles.png` - Cross-shore wave profiles
- `bathymetry.png` - Depth profile
- `bathymetry_comparison.png` - Method comparison
- `comprehensive_analysis.png` - Multi-panel overview

---

## Advanced Features

### Custom Preprocessing

```python
from coastcams import ImagePreprocessor

preprocessor = ImagePreprocessor()
preprocessor.lowpass_cutoff = 0.1  # Hz
preprocessor.highpass_cutoff = 0.01  # Hz
preprocessor.smoothing_window = 5

# Apply custom preprocessing
processed = preprocessor.preprocess_image(img)
```

### Wave Group Analysis

```python
from coastcams import WaveAnalyzer

analyzer = WaveAnalyzer()

# Analyze wave grouping
timeseries = timestack[100, :]  # At specific location
group_analysis = analyzer.analyze_wave_groups(timeseries)

print(f"Number of wave groups: {group_analysis['num_groups']}")
print(f"Mean group period: {group_analysis['mean_group_period']:.1f} s")
```

### Spectral Analysis

```python
from coastcams import WaveAnalyzer

analyzer = WaveAnalyzer()

# Compute spectral moments
timeseries = timestack[100, :]
moments = analyzer.compute_spectral_moments(timeseries)

print(f"Spectral Hs: {moments['significant_height']:.2f} m")
```

---

## Troubleshooting

### Problem: No images found
**Solution**: Check that images are in the `Timestacks/` directory and match the expected naming pattern (`S_X_YYYYMMDDHHSS.jpeg`)

### Problem: Shoreline detection fails
**Solution**: Try different detection methods (1, 2, or 3) or adjust the threshold parameter

### Problem: Unrealistic depth values
**Solution**: Check wave parameter estimates, adjust pixel resolution, or use different bathymetry method

### Problem: Import errors
**Solution**: Ensure all dependencies are installed: `pip install -r requirements.txt`

---

## Performance Tips

- **Large datasets**: Process images in batches or use `max_images` parameter
- **Memory issues**: Reduce `cross_shore_width` or process fewer images at once
- **Speed**: Disable plotting (`save_plots: false`) for faster processing

---

## Comparison with MATLAB Version

| Feature | MATLAB | Python |
|---------|--------|--------|
| Core functionality | ✅ | ✅ |
| Shoreline detection | ✅ | ✅ (3 methods) |
| Wave analysis | ✅ | ✅ |
| Bathymetry estimation | ✅ | ✅ |
| Visualization | ✅ | ✅ (Enhanced) |
| Configuration | Script-based | YAML + Programmatic |
| Modularity | Limited | High |
| User-friendliness | Expert | Non-expert |
| Automation | Manual | Fully automated |

---

## Contributing

Contributions are welcome! Areas for improvement:
- Additional shoreline detection methods
- Machine learning integration
- Real-time processing capabilities
- GPU acceleration
- Web interface

---

## Citation

If you use CoastCams in your research, please cite:

```
Nuyts, S., et al. (2023). CoastCams: A MATLAB toolbox for coastal monitoring.
Environmental Modelling & Software, 168, 105800.
```

---

## License

GNU General Public License v3.0

---

## Support

For questions, issues, or feature requests:
- Check the examples in `examples/` directory
- Review configuration options in `config.yaml`
- Examine module documentation in source files

---

## Acknowledgments

This Python port maintains compatibility with the original MATLAB CoastCams toolbox while adding:
- Modular architecture for flexible usage
- Automated workflows for non-experts
- Enhanced visualization capabilities
- Modern Python scientific computing stack

---

**Version**: 1.0.0
**Last Updated**: 2024

**Happy Coastal Monitoring!** 🌊
