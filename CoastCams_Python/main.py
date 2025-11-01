"""
CoastCams - Python Implementation
Coastal Wave Analysis from Timestack Images

This script computes nearshore processes from timestack images:
- Wave height, wave period, wave celerity
- Mean water levels and sea level anomaly
- Morphology (shoreline positions)

Based on MATLAB script S01_AnalysisTimestackImages.m
"""

import os
import sys
import numpy as np
from datetime import datetime
from pathlib import Path
from typing import List, Dict, Optional
import warnings

# Import CoastCams modules
from coastcams.config import CoastCamsConfig
from coastcams.image_loader import ImageLoader
from coastcams.shoreline import ShorelineDetector
from coastcams.preprocessing import ImagePreprocessor
from coastcams.wave_analysis import WaveAnalyzer
from coastcams.cross_correlation import CrossCorrelationAnalyzer
from coastcams.bathymetry import BathymetryEstimator
from coastcams.visualize import CoastCamsVisualizer
from coastcams.matlab_preprocessing import MATLABPreprocessor

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore', category=RuntimeWarning)


def smooth2(data: np.ndarray, Nr: int, Nc: int) -> np.ndarray:
    """
    2D smoothing function matching MATLAB's smooth2.

    Parameters
    ----------
    data : np.ndarray
        2D array to smooth
    Nr : int
        Number of rows for smoothing window
    Nc : int
        Number of columns for smoothing window

    Returns
    -------
    np.ndarray
        Smoothed 2D array
    """
    from scipy.ndimage import uniform_filter1d

    if data.ndim != 2:
        raise ValueError("Input must be 2D array")

    result = np.copy(data)

    # Column-wise smoothing (spatial dimension)
    if Nc > 0:
        for i in range(result.shape[0]):
            result[i, :] = uniform_filter1d(result[i, :], size=Nc, mode='nearest')

    # Row-wise smoothing (time dimension) - minimal with Nr=1
    if Nr > 1:
        for j in range(result.shape[1]):
            result[:, j] = uniform_filter1d(result[:, j], size=Nr, mode='nearest')

    return result


def movmean(data: np.ndarray, window: int) -> np.ndarray:
    """
    Moving average matching MATLAB's movmean.

    Parameters
    ----------
    data : np.ndarray
        Input array
    window : int
        Window size for moving average

    Returns
    -------
    np.ndarray
        Smoothed array
    """
    from scipy.ndimage import uniform_filter1d

    if data.ndim == 1:
        return uniform_filter1d(data, size=window, mode='nearest')
    else:
        # Apply along last dimension
        result = np.zeros_like(data)
        if data.ndim == 2:
            for i in range(data.shape[0]):
                result[i, :] = uniform_filter1d(data[i, :], size=window, mode='nearest')
        return result


def main():
    """
    Main analysis function following MATLAB S01_AnalysisTimestackImages.m workflow.
    """

    print("\n" + "="*70)
    print("CoastCams - Coastal Wave Analysis")
    print("="*70 + "\n")

    # =====================================================================
    # A: Configuration Setup
    # =====================================================================

    # Load configuration
    config_path = Path(__file__).parent / 'config.yaml'
    if config_path.exists():
        config = CoastCamsConfig(str(config_path))
    else:
        print(f"Configuration file not found at {config_path}, using defaults")
        config = CoastCamsConfig()

    # Set up paths
    repo_path = Path(__file__).parent.parent
    img_path = repo_path / 'Timestacks'

    # Resolve output directory path
    out_path = Path(config.output_dir)
    if not out_path.is_absolute():
        # If relative, resolve relative to repo path
        out_path = repo_path / config.output_dir

    # Create output directory
    try:
        out_path.mkdir(parents=True, exist_ok=True)
    except PermissionError as e:
        # Fall back to a safe location if we don't have permissions
        import tempfile
        out_path = Path(tempfile.gettempdir()) / 'CoastCams_Output'
        out_path.mkdir(parents=True, exist_ok=True)
        print(f"Warning: Could not create output directory at {config.output_dir}")
        print(f"Using temporary directory: {out_path}")

    # Calculate dt from acquisition frequency
    dt = 1.0 / config.acquisition_frequency

    # Set icmax if not specified
    if config.max_cross_shore is None:
        config.max_cross_shore = 680  # Default from MATLAB

    print("\n" + "="*60)
    print("CoastCams Configuration Summary")
    print("="*60 + "\n")
    print(f"Camera Parameters:")
    print(f"  Height above MSL: {config.camera_height:.3f} m")
    print(f"  Pixel resolution: {config.pixel_resolution:.3f} m/pixel")
    print(f"  Rotation angle: {config.rotation_angle}°")
    print(f"  Acquisition frequency: {config.acquisition_frequency} Hz")
    print(f"\nWave Analysis:")
    print(f"  Period range: {config.wave_period_min}-{config.wave_period_max} s")
    print(f"  Correlation spacing: {config.correlation_spacing} pixels")
    print(f"  Cross-shore width: {config.cross_shore_width} pixels")
    print(f"\nShoreline Detection:")
    print(f"  Method: {config.shoreline_method}")
    print(f"  Threshold: {config.shoreline_threshold}")
    print(f"\nOutput Settings:")
    print(f"  Input directory: {img_path}")
    print(f"  Output directory: {out_path}")
    print(f"  Output interval: {config.output_interval} minutes")
    print(f"  Output format: {config.output_format}")
    print(f"  Save plots: {config.save_plots}")
    print("="*60 + "\n")

    # =====================================================================
    # B: Initialize Result Storage (MATLAB line 86)
    # =====================================================================

    # Initialize lists/arrays to store results from each timestack
    Time_TS = []              # Timestamps
    WaveCelerity = []         # Wave celerity matrix (timestacks × spatial points)
    WaveLength = []           # Wavelength matrix
    Hs_TS = []                # Significant wave height per timestack
    Tp_TS = []                # Peak period per timestack
    Tm_TS = []                # Mean period per timestack
    WaveEnergy = []           # Wave energy
    RollerLength = []         # Roller length
    BreakpointDepth = []      # Breakpoint depth
    BreakpointLocation = []   # Breakpoint location
    WaterDepth = []           # Water depth from photogrammetry (depth_s)
    WaterDepth_L = []         # Water depth from linear wave theory
    ShorePosition = []        # Shoreline positions
    Stack_av = []             # Average timestack intensities

    # =====================================================================
    # C: Load Images
    # =====================================================================

    print("[1/7] Loading timestack images...")
    loader = ImageLoader(str(img_path), rotation_angle=config.rotation_angle)

    # Discover images in directory
    num_images = loader.discover_images(pattern="S_*.jpeg")

    if num_images == 0:
        print(f"Error: No images found in {img_path}")
        return

    # Print summary
    print(loader.get_summary())

    # Initialize analyzers
    shoreline_detector = ShorelineDetector(
        method=config.shoreline_method,
        threshold=config.shoreline_threshold
    )

    preprocessor = ImagePreprocessor(config=config)

    wave_analyzer = WaveAnalyzer(config=config)

    correlation_analyzer = CrossCorrelationAnalyzer(config=config)

    matlab_preprocessor = MATLABPreprocessor(dt=dt, use_radon=True)

    bathymetry_estimator = BathymetryEstimator(config=config)

    visualizer = CoastCamsVisualizer(
        output_dir=str(out_path),
        save_plots=True
    )

    print(f"\nNOTE: Processing {len(loader.image_files)} timestack images individually")
    print("(Each image is a space x time array analyzed independently)\n")

    # =====================================================================
    # D: Process Each Timestack (MATLAB line 108-300)
    # =====================================================================

    dc = config.correlation_spacing

    for i, (timestamp, image_path) in enumerate(zip(loader.timestamps, loader.image_files)):
        print("="*70)
        print(f"Image {i+1}/{len(loader.image_files)}: {timestamp}")
        print("="*70)

        # G1: Extract datetime (MATLAB line 110)
        Time_TS.append(timestamp)

        # G2: Load timestack (MATLAB line 116)
        try:
            timestack = loader.load_image(i)
            print(f"  Loaded timestack shape: {timestack.shape}")
        except Exception as e:
            print(f"  Error loading image: {e}")
            # Append NaN values for this timestack
            _append_nan_results(i, WaveCelerity, WaveLength, Hs_TS, Tp_TS, Tm_TS,
                               WaveEnergy, RollerLength, BreakpointDepth,
                               BreakpointLocation, WaterDepth, WaterDepth_L,
                               ShorePosition, Stack_av)
            continue

        # G3: Extract shoreline position (MATLAB line 126)
        print("[2/7] Detecting shoreline...")
        try:
            shoreline = shoreline_detector.detect(timestack)
            shoreline_pos = np.nanmean(shoreline) * config.pixel_resolution
            ShorePosition.append(shoreline_pos)
        except Exception as e:
            print(f"  Warning: Shoreline detection failed: {e}")
            ShorePosition.append(np.nan)

        # G4-G5: Preprocessing (MATLAB line 169-184)
        print("[3/7] Preprocessing timestack...")
        try:
            # Convert to grayscale if needed and preprocess
            if len(timestack.shape) == 3:
                # Use blue channel (MATLAB line 175: S0(:,:,3))
                timestack_gray = timestack[:, :, 2].astype(np.float64)
            else:
                timestack_gray = timestack.astype(np.float64)

            preprocessed = preprocessor.preprocess_image(timestack_gray)
            print(f"  Preprocessed timestack shape: {preprocessed.shape}")
        except Exception as e:
            print(f"  Error in preprocessing: {e}")
            import traceback
            traceback.print_exc()
            _append_nan_results(i, WaveCelerity, WaveLength, Hs_TS, Tp_TS, Tm_TS,
                               WaveEnergy, RollerLength, BreakpointDepth,
                               BreakpointLocation, WaterDepth, WaterDepth_L,
                               ShorePosition, Stack_av)
            continue

        # G6: Compute Wave Parameters (MATLAB line 210)
        # This is equivalent to calling WaveParameters_CoastCams in MATLAB
        print("[4/7] Analyzing wave parameters...")
        try:
            # Create cross-shore positions array
            cross_shore_positions = np.arange(preprocessed.shape[1]) * config.pixel_resolution

            wave_results = wave_analyzer.analyze_timestack(
                preprocessed,
                cross_shore_positions=cross_shore_positions
            )

            # Extract results
            hs = wave_results.get('mean_Hs', np.nan)
            hm = wave_results.get('mean_Hm', np.nan)
            Tm = wave_results.get('mean_Tm', np.nan)
            Tp = wave_results.get('mean_Tp', np.nan)
            roller_length = wave_results.get('roller_length', np.nan)
            break_location = wave_results.get('break_location', np.nan)
            break_depth = wave_results.get('break_depth', np.nan)

            # depth_s from photogrammetry (MATLAB line 210, 218)
            depth_s = wave_results.get('depth_profile', None)

            Hs_TS.append(hs)
            Tm_TS.append(Tm)
            Tp_TS.append(Tp)
            RollerLength.append(roller_length)
            BreakpointLocation.append(break_location + (dc / 2))  # MATLAB line 229
            BreakpointDepth.append(break_depth)

        except Exception as e:
            print(f"  Error in wave analysis: {e}")
            import traceback
            traceback.print_exc()
            Hs_TS.append(np.nan)
            Tm_TS.append(np.nan)
            Tp_TS.append(np.nan)
            RollerLength.append(np.nan)
            BreakpointLocation.append(np.nan)
            BreakpointDepth.append(np.nan)
            depth_s = None

        # Store depth_s with movmean smoothing (MATLAB line 218)
        if depth_s is not None and not np.all(np.isnan(depth_s)):
            # MATLAB divides by 10 and applies movmean(10)
            # This suggests depth_s might be in different units (cm?)
            depth_s_smoothed = movmean(depth_s, 10)
            WaterDepth.append(depth_s_smoothed)
        else:
            # Store NaN array of appropriate size
            WaterDepth.append(np.full(config.max_cross_shore - dc, np.nan))

        # G7: Cross-correlation calculations (MATLAB line 234)
        print("[5/7] Performing cross-correlation...")
        try:
            corr_results = correlation_analyzer.analyze_timestack(preprocessed)

            # Extract celerity and wavelength (MATLAB line 234)
            # Keys are 'Cf1' and 'WLe1' (raw MATLAB outputs)
            Cf1 = corr_results.get('Cf1', np.array([]))
            WLe1 = corr_results.get('WLe1', np.array([]))

            if len(Cf1) > 0:
                # Apply movmean smoothing (MATLAB line 242-243)
                # MATLAB divides by 10 - likely converting units
                Cf1_smoothed = movmean(Cf1, 10)
                WLe1_smoothed = movmean(WLe1, 10)

                WaveCelerity.append(Cf1_smoothed)
                WaveLength.append(WLe1_smoothed)
            else:
                WaveCelerity.append(np.full(config.max_cross_shore - dc, np.nan))
                WaveLength.append(np.full(config.max_cross_shore - dc, np.nan))

        except Exception as e:
            print(f"  Error in cross-correlation: {e}")
            import traceback
            traceback.print_exc()
            WaveCelerity.append(np.full(config.max_cross_shore - dc, np.nan))
            WaveLength.append(np.full(config.max_cross_shore - dc, np.nan))

        # G8: Calculate depth using linear wave theory (MATLAB line 247)
        print("[6/7] Calculating water depth from linear wave theory...")
        try:
            # Get the celerity array we just calculated
            if len(WaveCelerity) > 0 and not np.all(np.isnan(WaveCelerity[-1])):
                celerity_array = WaveCelerity[-1]

                # Call LinearC (bathymetry estimator) - MATLAB line 247
                depths_linear = bathymetry_estimator.estimate_depth_linear_wave_theory(
                    peak_period=Tp if not np.isnan(Tp) else 10.0,
                    celerities=celerity_array
                )

                # Apply movmean smoothing (MATLAB line 252)
                depths_linear_smoothed = movmean(depths_linear, 10)
                WaterDepth_L.append(depths_linear_smoothed)
            else:
                WaterDepth_L.append(np.full(config.max_cross_shore - dc, np.nan))
        except Exception as e:
            print(f"  Error in depth calculation: {e}")
            WaterDepth_L.append(np.full(config.max_cross_shore - dc, np.nan))

        # G10: Average timestack for visualization (MATLAB line 281)
        try:
            avg_profile = np.mean(preprocessed, axis=0)
            Stack_av.append(avg_profile)
        except Exception as e:
            print(f"  Error calculating average timestack: {e}")
            Stack_av.append(np.full(preprocessed.shape[1], np.nan))

        print("[7/7] Timestack analysis complete")
        print(f"  → Wave Height: {hs:.3f} m")
        print(f"  → Wave Period (Tm): {Tm:.2f} s")
        print(f"  → Wave Period (Tp): {Tp:.2f} s")
        if not np.isnan(break_location):
            print(f"  → Break Location: {break_location:.2f} m")
        print()

    # =====================================================================
    # E: Convert Lists to Matrices
    # =====================================================================

    print("\n" + "="*70)
    print("Post-Processing: Building Matrices and Calculating SLA")
    print("="*70 + "\n")

    # Convert lists to numpy arrays
    # Ensure all rows have same length by padding with NaN
    def pad_to_matrix(list_of_arrays):
        """Convert list of arrays to 2D matrix, padding with NaN."""
        if len(list_of_arrays) == 0:
            return np.array([[]])

        # Find maximum length
        max_len = max(len(arr) if hasattr(arr, '__len__') and not isinstance(arr, (int, float))
                      else 1 for arr in list_of_arrays)

        # Pad each array
        matrix = []
        for arr in list_of_arrays:
            if not hasattr(arr, '__len__') or isinstance(arr, (int, float)):
                # Scalar value
                padded = np.full(max_len, arr)
            else:
                arr = np.array(arr)
                if len(arr) < max_len:
                    padded = np.pad(arr, (0, max_len - len(arr)),
                                   constant_values=np.nan)
                else:
                    padded = arr[:max_len]
            matrix.append(padded)

        return np.array(matrix)

    WaveCelerity = pad_to_matrix(WaveCelerity)
    WaveLength = pad_to_matrix(WaveLength)
    WaterDepth = pad_to_matrix(WaterDepth)
    WaterDepth_L = pad_to_matrix(WaterDepth_L)
    Stack_av = pad_to_matrix(Stack_av)

    print(f"WaveCelerity matrix: {WaveCelerity.shape} (timestacks × spatial points)")
    print(f"WaterDepth_L matrix: {WaterDepth_L.shape} (timestacks × spatial points)")
    print(f"WaterDepth (photogrammetry) matrix: {WaterDepth.shape}")

    # =====================================================================
    # F: Calculate Sea Level Anomaly (MATLAB line 256-272)
    # =====================================================================

    print("\nCalculating Sea Level Anomaly (SLA)...")

    # Determine smoothing parameters (MATLAB line 257-259)
    rows, cols = WaveCelerity.shape
    Nr = min(rows, 1)   # Minimal row smoothing
    Nc = min(cols, 30)  # Spatial smoothing

    try:
        # Smooth WaveCelerity and WaterDepth_L (MATLAB line 260, 264)
        print(f"  Applying smooth2 with Nr={Nr}, Nc={Nc}...")
        Csmooth_S = smooth2(WaveCelerity, Nr, Nc)
        Csmooth_L = smooth2(WaterDepth_L, Nr, Nc)

        # Calculate SLA (MATLAB line 266-267)
        # SLA = smoothed - mean(smoothed)
        sLimit_size = min(Csmooth_S.shape[1], Csmooth_L.shape[1])

        SLA_S = Csmooth_S[:, :sLimit_size] - np.nanmean(Csmooth_S[:, :sLimit_size])
        SLA_L = Csmooth_L[:, :sLimit_size] - np.nanmean(Csmooth_L[:, :sLimit_size])

        print(f"  SLA_S matrix: {SLA_S.shape}")
        print(f"  SLA_L matrix: {SLA_L.shape}")

        # Calculate per-timestack averages (for CSV export)
        sla_shallow_values = np.nanmean(SLA_S, axis=1)  # Average across space
        sla_linear_values = np.nanmean(SLA_L, axis=1)   # Average across space

        print(f"  SLA_S range: {np.nanmin(SLA_S):.3f} to {np.nanmax(SLA_S):.3f} m/s")
        print(f"  SLA_L range: {np.nanmin(SLA_L):.3f} to {np.nanmax(SLA_L):.3f} m")

    except Exception as e:
        print(f"  Error calculating SLA: {e}")
        import traceback
        traceback.print_exc()
        SLA_S = np.full_like(WaveCelerity, np.nan)
        SLA_L = np.full_like(WaterDepth_L, np.nan)
        sla_shallow_values = np.full(len(Time_TS), np.nan)
        sla_linear_values = np.full(len(Time_TS), np.nan)

    # =====================================================================
    # G: Calculate Relative Tidal Range (RTR) (MATLAB line 275-278)
    # =====================================================================

    print("\nCalculating Relative Tidal Range (RTR)...")

    try:
        Level_TS = np.nanmean(WaveCelerity, axis=1)  # Average water level per timestack
        nRTR_thresh = Level_TS - np.nanmin(Level_TS)
        nRTR_thresh[nRTR_thresh < 0.2] = np.nan  # Threshold to avoid divide by zero
        RTR = np.array(Hs_TS) / nRTR_thresh

        print(f"  RTR: mean={np.nanmean(RTR):.3f}, range=[{np.nanmin(RTR):.3f}, {np.nanmax(RTR):.3f}]")
    except Exception as e:
        print(f"  Error calculating RTR: {e}")
        RTR = np.full(len(Time_TS), np.nan)

    # =====================================================================
    # H: Calculate Bathymetry (MATLAB line 310-366)
    # =====================================================================

    print("\nCalculating Bathymetry...")

    try:
        # Bathymetry = WaterDepth - SLA_S (MATLAB line 338)
        # Need to match dimensions
        min_rows = min(WaterDepth.shape[0], SLA_S.shape[0])
        min_cols = min(WaterDepth.shape[1], SLA_S.shape[1])

        WaterDepth_adj = WaterDepth[:min_rows, :min_cols]
        SLA_S_adj = SLA_S[:min_rows, :min_cols]

        Bathymetry_full = WaterDepth_adj - SLA_S_adj

        # Mean bathymetry across all timesteps (MATLAB line 341)
        Bathymetry = np.nanmean(Bathymetry_full, axis=0)

        print(f"  Bathymetry profile: {len(Bathymetry)} points")
        print(f"  Depth range: {np.nanmin(Bathymetry):.2f} to {np.nanmax(Bathymetry):.2f} m")

        # Create cross-shore distance array
        x_bathymetry = np.arange(len(Bathymetry)) * config.pixel_resolution

    except Exception as e:
        print(f"  Error calculating bathymetry: {e}")
        import traceback
        traceback.print_exc()
        Bathymetry = np.full(WaterDepth.shape[1], np.nan)
        x_bathymetry = np.arange(len(Bathymetry)) * config.pixel_resolution

    # =====================================================================
    # I: Calculate per-timestack averaged values for CSV export
    # =====================================================================

    print("\nCalculating averaged values for CSV export...")

    # Water depth per timestack (spatially averaged)
    water_depth_values = np.nanmean(WaterDepth_L, axis=1)

    print(f"  Water depths (per timestack): {len(water_depth_values)} values")
    print(f"  Range: {np.nanmin(water_depth_values):.2f} to {np.nanmax(water_depth_values):.2f} m")

    # =====================================================================
    # J: Store Results in Dictionary
    # =====================================================================

    results = {
        'timestamps': Time_TS,
        'wave_heights_timeseries': np.array(Hs_TS),
        'wave_periods_timeseries': np.array(Tm_TS),
        'wave_periods_peak_timeseries': np.array(Tp_TS),
        'wavelengths': WaveLength,
        'celerities': WaveCelerity,
        'roller_lengths': np.array(RollerLength),
        'breakpoint_locations': np.array(BreakpointLocation),
        'breakpoint_depths': np.array(BreakpointDepth),
        'shoreline_positions': np.array(ShorePosition),
        'depths': water_depth_values,  # Per-timestack averaged depth
        'sla_values': sla_linear_values,  # SLA_L (linear wave theory)
        'sla_shallow_values': sla_shallow_values,  # SLA_S (shallow water)
        'rtr_values': RTR,
        # Matrices for visualization
        'sla_matrix_linear': SLA_L,
        'sla_matrix_shallow': SLA_S,
        'average_timestack': Stack_av,
        # Bathymetry
        'bathymetry': Bathymetry,
        'bathymetry_x': x_bathymetry,
        # Aggregate statistics
        'mean_Hs': np.nanmean(Hs_TS),
        'mean_Tm': np.nanmean(Tm_TS),
        'mean_Tp': np.nanmean(Tp_TS),
        'mean_depth': np.nanmean(water_depth_values),
    }

    # =====================================================================
    # K: Export Results
    # =====================================================================

    print("\n" + "="*70)
    print("Exporting Results")
    print("="*70 + "\n")

    # Export to CSV
    timestamp_str = datetime.now().strftime('%Y%m%d_%H%M%S')
    csv_filename = f'coastcams_results_{timestamp_str}.csv'
    csv_path = out_path / csv_filename

    try:
        import pandas as pd

        # Create DataFrame with per-timestack values
        df = pd.DataFrame({
            'Timestamp': Time_TS,
            'Hs': Hs_TS,
            'Tm': Tm_TS,
            'Tp': Tp_TS,
            'WaterDepth': water_depth_values,
            'SLA_L': sla_linear_values,
            'SLA_S': sla_shallow_values,
            'RTR': RTR,
            'BreakpointLocation': BreakpointLocation,
            'BreakpointDepth': BreakpointDepth,
            'ShorelinePosition': ShorePosition,
            'RollerLength': RollerLength,
        })

        df.to_csv(csv_path, index=False)
        print(f"Results exported to: {csv_path}")

    except Exception as e:
        print(f"Error exporting to CSV: {e}")
        import traceback
        traceback.print_exc()

    # Create summary report
    summary_path = out_path / 'analysis_summary.txt'
    try:
        with open(summary_path, 'w') as f:
            f.write("="*70 + "\n")
            f.write("CoastCams Analysis Summary Report\n")
            f.write("="*70 + "\n\n")
            f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            f.write(f"Number of Images Processed: {len(Time_TS)}\n")
            f.write(f"Time Range: {Time_TS[0]} to {Time_TS[-1]}\n\n")
            f.write("Analysis Results:\n")
            f.write("-"*70 + "\n")
            f.write(f"  Mean Significant Wave Height: {results['mean_Hs']:.3f} m\n")
            f.write(f"  Mean Wave Period (Tm): {results['mean_Tm']:.3f} s\n")
            f.write(f"  Mean Peak Period (Tp): {results['mean_Tp']:.3f} s\n")
            f.write(f"  Mean Water Depth: {results['mean_depth']:.3f} m\n")
            f.write("\n" + "="*70 + "\n")

        print(f"Summary report saved to {summary_path}")

    except Exception as e:
        print(f"Error creating summary report: {e}")

    # =====================================================================
    # L: Create Visualizations (MATLAB line 437)
    # =====================================================================

    print("\nCreating visualizations...")

    try:
        # MATLAB-style summary plot
        visualizer.plot_matlab_style_summary(
            timestamps=Time_TS,
            average_timestack=Stack_av,
            sla_matrix=SLA_L,
            wave_heights=np.array(Hs_TS),
            wave_periods=np.array(Tp_TS),
            water_levels=water_depth_values,
            rotation=config.rotation_angle,
            filename='coastcams_matlab_summary.png'
        )

    except Exception as e:
        print(f"Error creating visualizations: {e}")
        import traceback
        traceback.print_exc()

    # =====================================================================
    # M: Completion
    # =====================================================================

    print("\n" + "="*70)
    print("CoastCams Analysis Completed Successfully")
    print("="*70)
    print(f"\nResults saved in: {out_path}")
    print(f"CSV file: {csv_filename}")
    print(f"Summary report: analysis_summary.txt\n")

    return results


def _append_nan_results(i, WaveCelerity, WaveLength, Hs_TS, Tp_TS, Tm_TS,
                       WaveEnergy, RollerLength, BreakpointDepth,
                       BreakpointLocation, WaterDepth, WaterDepth_L,
                       ShorePosition, Stack_av):
    """Helper function to append NaN values when image processing fails."""
    Hs_TS.append(np.nan)
    Tp_TS.append(np.nan)
    Tm_TS.append(np.nan)
    RollerLength.append(np.nan)
    BreakpointLocation.append(np.nan)
    BreakpointDepth.append(np.nan)
    ShorePosition.append(np.nan)

    # Append NaN arrays for matrix results
    nan_array = np.full(600, np.nan)  # Default size
    WaveCelerity.append(nan_array)
    WaveLength.append(nan_array)
    WaterDepth.append(nan_array)
    WaterDepth_L.append(nan_array)
    Stack_av.append(nan_array)


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print("\n\nAnalysis interrupted by user.")
        sys.exit(1)
    except Exception as e:
        print(f"\n\nFatal error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
