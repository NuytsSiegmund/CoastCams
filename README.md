# CoastCams

CoastCams is an open-source MATLAB toolbox designed to quantify key wave parameters, mean water levels, bathymetry, and morphology in the nearshore environment using oblique orthorectified timestack images from land-based coastal monitoring systems.

## Key Features

- Unified and simplified method for coastal monitoring
- Accessible to coastal managers, engineers, and scientists
- User-friendly approach to identify key drivers in the coastal zone
- Enhanced bathymetry measurements for improved depth estimation
- Comprehensive set of functions for diverse coastal analyses
- Automated error handling and logging
- Improved image processing and analysis workflows

## Python Version Available

**NEW!** A complete Python port of CoastCams is now available in the [`CoastCams_Python/`](CoastCams_Python/) directory.

The Python version offers:
- ✅ **Fully automated workflow** - minimal user input required
- ✅ **Modular design** - run complete analysis or specific tasks
- ✅ **Non-expert friendly** - smart defaults and comprehensive documentation
- ✅ **No MATLAB required** - uses open-source Python libraries
- ✅ **Enhanced features** - additional analysis methods and visualizations

Quick start:
```bash
cd CoastCams_Python
pip install -r requirements.txt
python main.py
```

See [`CoastCams_Python/README.md`](CoastCams_Python/README.md) for full documentation.

## Citation

When using CoastCams, please cite:

Nuyts, S., Almar, R., Morichon, D., Dealbera, S., Abalia, A., Muñoz, J. M., Abessolo, G. O., & Regard, V. (2023). CoastCams: A MATLAB toolbox making accessible estimations of nearshore processes, mean water levels, and morphology from timestack images. Environmental Modelling & Software, 168, 105800. https://doi.org/10.1016/j.envsoft.2023.105800

## Installation

1. Clone the repository:
   ```
   git clone https://github.com/NuytsSiegmund/CoastCams.git
   ```
2. Add the CoastCams directory and its subfolders to your MATLAB path.

## Usage

### Overview

CoastCams analyzes georectified timestack images from coastal video cameras. The creation of timestack images and georectification can be achieved with the [Quantitative Coastal Imaging Toolbox](https://github.com/Coastal-Imaging-Research-Network/CIRN-Quantitative-Coastal-Imaging-Toolbox).

### Configuration

Open `UserScripts/S01_AnalysisTimestackImages.m` and configure the following sections:

1. **Paths**: Set paths for input images and output results.
2. **Image Selection**: Specify the naming convention for your timestack images.
3. **Processing Parameters**: Set camera and image processing parameters.
4. **Analysis Options**: Choose shoreline detection methods and plotting options.

### Key Parameters

- `dt`: Camera acquisition frequency (e.g., 2 images per second)
- `H_camera`: Camera height above MSL in meters
- `res`: Pixel size on timestack image in meters
- `rotation`: Image rotation to ensure waves come from the top-left corner
- `dur`: Duration of the timestack images in minutes
- `ShoreMethod`: Shoreline detection method (1, 2, or 3)
- `compare_bathymetry_option`: Enable comparison of bathymetry methods

### Running the Analysis

Execute `S01_AnalysisTimestackImages.m` in MATLAB. The script will:

1. Process all timestack images in the specified directory
2. Calculate wave parameters, water levels, and bathymetry
3. Generate plots and save results

## New Features

### Enhanced Bathymetry Measurements

The updated version includes improved algorithms for bathymetry estimation:

- Calculates bathymetry for each timestep
- Computes mean bathymetry across all timesteps
- Generates a bathymetry profile plot
- Saves bathymetry data for further analysis

### Additional Functions

- `WaveParameters_CoastCams`: Computes various wave parameters
- `CrossCorrelation_CoastCams`: Performs cross-correlation analysis
- `LinearC`: Calculates depth using linear wave theory
- `plot_coastcams_main`: Generates main plots for analysis results

### Error Handling and Logging

Improved error handling with a dedicated log file for tracking issues during processing.

## Outputs

- Timetable of processed data (`WP`)
- Daily averaged data (`WP_daily`)
- Saved workspace with all variables
- Text file with processed data
- Various plots saved as PNG files
- Error log
- Analysis summary report

## Contributing and Issues
Having a problem? Post an issue in the [Issues Page](https://github.com/NuytsSiegmund/CoastCams/issues)

If you're willing to contribute: 

1. Fork the repository. A fork is a copy on which you can make your changes.
2. Create a new branch on your fork
3. Commit your changes and push them to your branch
4. When the branch is ready to be merged, create a Pull Request (how to make a clean pull request explained [here](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request))



