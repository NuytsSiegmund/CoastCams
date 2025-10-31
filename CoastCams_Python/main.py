#!/usr/bin/env python3
"""
CoastCams - Main Workflow Script

Fully automated coastal wave analysis from timestack images.

This script provides a complete, user-friendly workflow for analyzing
coastal wave dynamics, bathymetry, and shoreline positions.

Usage:
    python main.py                    # Run with default settings
    python main.py --config my_config.yaml    # Run with custom config
    python main.py --help             # Show help
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import cv2
from datetime import datetime
from pathlib import Path
from tqdm import tqdm

# Import CoastCams modules
from coastcams.config import CoastCamsConfig
from coastcams.image_loader import ImageLoader
from coastcams.preprocessing import ImagePreprocessor
from coastcams.shoreline import ShorelineDetector
from coastcams.wave_analysis import WaveAnalyzer
from coastcams.cross_correlation import CrossCorrelationAnalyzer
from coastcams.bathymetry import BathymetryEstimator
from coastcams.sea_level import SeaLevelAnalyzer
from coastcams.visualize import CoastCamsVisualizer


class CoastCamsWorkflow:
    """
    Complete CoastCams analysis workflow.

    This class orchestrates the entire analysis pipeline from
    image loading to result export and visualization.
    """

    def __init__(self, config_file: str = None):
        """
        Initialize workflow.

        Parameters
        ----------
        config_file : str, optional
            Path to configuration file
        """
        # Load configuration
        self.config = CoastCamsConfig(config_file)

        # Auto-detect paths if not specified
        self.config.auto_detect_paths()

        # Validate configuration
        if not self.config.validate():
            print("Configuration validation failed. Please check your settings.")
            sys.exit(1)

        # Initialize modules (resolve paths relative to config file)
        self.image_loader = ImageLoader(
            self.config.resolve_path(self.config.input_dir),
            self.config.rotation_angle,
            self.config.verbose
        )

        self.preprocessor = ImagePreprocessor(self.config)
        self.shoreline_detector = ShorelineDetector(config=self.config)
        self.wave_analyzer = WaveAnalyzer(self.config)
        self.correlation_analyzer = CrossCorrelationAnalyzer(self.config)
        self.bathymetry_estimator = BathymetryEstimator(self.config)
        self.sea_level_analyzer = SeaLevelAnalyzer(self.config)
        self.visualizer = CoastCamsVisualizer(
            self.config.resolve_path(self.config.output_dir),
            self.config.save_plots
        )

        # Storage for results
        self.results = {}

    def run_full_analysis(self):
        """
        Run complete CoastCams analysis workflow.

        This is the main method that executes all analysis steps.
        """
        print("\n" + "="*70)
        print("CoastCams - Coastal Wave Analysis")
        print("="*70 + "\n")

        # Print configuration summary
        self.config.print_summary()

        # Step 1: Load images
        print("\n[1/7] Loading timestack images...")
        num_images = self.image_loader.discover_images()

        if num_images == 0:
            print("Error: No images found. Please check input directory.")
            return

        print(self.image_loader.get_summary())
        print(f"\nNOTE: Processing {num_images} timestack images individually")
        print("(Each image is a space x time array analyzed independently)\n")

        # Process each timestack image individually (matching MATLAB workflow)
        all_results = []

        for img_idx in range(num_images):
            print(f"\n{'='*70}")
            print(f"Image {img_idx+1}/{num_images}: {self.image_loader.timestamps[img_idx]}")
            print(f"{'='*70}")

            # Load single image (this image IS a timestack: space x time)
            img = self.image_loader.load_image(img_idx)
            print(f"  Loaded timestack shape: {img.shape}")

            # Step 2: Detect shoreline for this timestack
            print("[2/7] Detecting shoreline...")
            shoreline = self.shoreline_detector.detect(img)

            # Step 3: Preprocess timestack
            print("[3/7] Preprocessing timestack...")
            timestack = self._preprocess_single_timestack(img)
            print(f"  Preprocessed timestack shape: {timestack.shape}")

            # Compute average timestack (MATLAB line 281: mean(Timestack_rot, 1, 'omitnan'))
            # timestack shape is (time, space), so mean along axis=0 (time) gives spatial profile
            avg_timestack_profile = np.nanmean(timestack, axis=0)  # Shape: (space,)

            # Step 4: Analyze wave parameters
            print("[4/7] Analyzing wave parameters...")
            wave_results = self._analyze_single_timestack(timestack, shoreline)

            # Step 5: Perform cross-correlation analysis (G7 in MATLAB)
            print("[5/7] Performing cross-correlation...")
            correlation_results = self.correlation_analyzer.analyze_timestack(timestack)

            # Extract raw Cf1 and WLe1 from cross-correlation (MATLAB outputs)
            # NOTE: In Python's cross_correlation.py, Cf1 is already:
            # 1. Converted to m/s (multiplied by pixel resolution)
            # 2. Smoothed with movmean(10)
            # So Cf1 is NOT divided by 10, and movmean is already applied
            Cf1 = correlation_results.get('Cf1', None)  # Already in m/s, already smoothed
            WLe1 = correlation_results.get('WLe1', None)  # Raw wavelengths

            # For MATLAB compatibility:
            # MATLAB divides Cf1 by 10 because their units are different
            # Python already has correct units (m/s), so we use Cf1 directly as WaveCelerity
            if Cf1 is not None and len(Cf1) > 0:
                WaveCelerity = Cf1  # Already in m/s and smoothed
                num_valid = np.sum(~np.isnan(WaveCelerity))
                print(f"  WaveCelerity array: {len(WaveCelerity)} spatial points ({num_valid} valid)")
                if num_valid > 0:
                    print(f"    Range: [{np.nanmin(WaveCelerity):.2f}, {np.nanmax(WaveCelerity):.2f}] m/s")
                else:
                    print(f"    Warning: All WaveCelerity values are NaN")
            else:
                WaveCelerity = None
                print(f"  Warning: No Cf1 available")

            if WLe1 is not None and len(WLe1) > 0:
                WaveLength = WLe1  # Already computed
                num_valid = np.sum(~np.isnan(WaveLength))
                print(f"  WaveLength array: {len(WaveLength)} spatial points ({num_valid} valid)")
                if num_valid > 0:
                    print(f"    Range: [{np.nanmin(WaveLength):.2f}, {np.nanmax(WaveLength):.2f}] m")
                else:
                    print(f"    Warning: All WaveLength values are NaN")
            else:
                WaveLength = None
                print(f"  Warning: No WLe1 available")

            # Step 6: Calculate WaterDepth_L using LinearC (G8 in MATLAB)
            # MATLAB line 247: [df] = LinearC(Tp, WaveCelerity(i,:), 0.01)
            # Tp is SCALAR, WaveCelerity is ARRAY
            print("[6/7] Calculating water depth from linear wave theory...")

            # Get Tp (peak period) - scalar
            Tp = wave_results.get('mean_Tp', wave_results.get('mean_Tm', 8.0))
            print(f"  Using peak period Tp = {Tp:.2f}s (scalar)")

            if WaveCelerity is not None:
                # Calculate depth for each spatial position using LinearC
                from coastcams.utils import calculate_depth_from_celerity

                WaterDepth_L = []
                num_valid = 0
                for c in WaveCelerity:
                    if not np.isnan(c) and c > 0:
                        depth = calculate_depth_from_celerity(c, Tp, precision=0.01)
                        WaterDepth_L.append(depth)
                        if not np.isnan(depth):
                            num_valid += 1
                    else:
                        WaterDepth_L.append(np.nan)

                WaterDepth_L = np.array(WaterDepth_L)

                # Don't apply movmean here - it will be applied in _aggregate_all_results
                # MATLAB: WaterDepth_L = [WaterDepth_L; movmean(df, 10)];
                # The movmean is applied to df before stacking, but we'll handle it in aggregation

                print(f"  WaterDepth_L array: {len(WaterDepth_L)} spatial points")
                print(f"    Valid depths: {num_valid}/{len(WaterDepth_L)}")
                if num_valid > 0:
                    print(f"    Depth range: [{np.nanmin(WaterDepth_L):.2f}, {np.nanmax(WaterDepth_L):.2f}] m")
                else:
                    print(f"    Warning: All depth values are NaN")

                # Mean depth for this timestack
                mean_depth = np.nanmean(WaterDepth_L)
                if not np.isnan(mean_depth):
                    print(f"  Mean water depth: {mean_depth:.2f} m")
                else:
                    print(f"  Warning: Mean water depth is NaN")
            else:
                WaterDepth_L = None
                mean_depth = np.nan
                print(f"  Warning: Cannot calculate WaterDepth_L (no WaveCelerity)")

            # For backward compatibility, create bathymetry_results structure
            bathymetry_results = {
                'mean_depth': mean_depth,
                'depths_filtered': WaterDepth_L,
                'cross_shore_positions': wave_results.get('cross_shore_positions', np.array([]))
            }

            # Debug NaN values
            mean_depth = bathymetry_results.get('mean_depth', np.nan)
            if np.isnan(mean_depth):
                print(f"  ⚠️  WARNING: mean_depth is NaN!")
                print(f"      - Celerities available: {len(correlation_results.get('celerities', []))} values")
                print(f"      - Mean celerity: {correlation_results.get('mean_celerity', np.nan):.3f}")
                print(f"      - Mean wave period: {wave_results.get('mean_Tm', np.nan):.2f} s")

            # Step 7: Compute sea level anomaly
            print("[7/7] Computing sea level anomaly...")
            sla = self._compute_sla_single(shoreline, bathymetry_results)

            # Calculate additional parameters to match MATLAB
            Hs = wave_results.get('mean_Hs', np.nan)
            Tm = wave_results.get('mean_Tm', np.nan)

            # Calculate scalars for CSV from arrays
            # MATLAB aggregates at the end: nanmean(WaveCelerity(i,:))
            celerity_mean = np.nanmean(WaveCelerity) if WaveCelerity is not None else np.nan
            wavelength_mean = np.nanmean(WaveLength) if WaveLength is not None else np.nan
            water_depth_mean = np.nanmean(WaterDepth_L) if WaterDepth_L is not None else np.nan

            # Calculate wave energy: E = Hs^2 / 2
            wave_energy = (Hs ** 2) / 2.0 if not np.isnan(Hs) else np.nan

            # Calculate shallow water depth: depth_s = C^2 / g (using mean celerity)
            depth_shallow = (celerity_mean ** 2) / 9.81 if not np.isnan(celerity_mean) else np.nan

            # Extract roller/breaking parameters from wave analysis
            roller_length = wave_results.get('roller_length', np.nan)
            breakpoint_location = wave_results.get('break_location', np.nan)

            # Calculate breaking depth using shallow water approximation at breaking location
            # MATLAB line 82-83: BD = depth_s(round(nanmean(PosX)))
            # where depth_s = C^2 / 9.81
            breakpoint_depth = np.nan
            if not np.isnan(breakpoint_location) and not np.isnan(depth_shallow):
                # Use the shallow water depth at this location as the breaking depth
                breakpoint_depth = depth_shallow
            elif not np.isnan(breakpoint_location) and WaterDepth_L is not None:
                # Try to get depth from WaterDepth_L profile at breaking location
                positions = bathymetry_results.get('cross_shore_positions', [])
                if len(positions) > 0 and len(WaterDepth_L) > 0:
                    # Find closest position to breakpoint
                    idx = np.argmin(np.abs(np.array(positions) - breakpoint_location))
                    if idx < len(WaterDepth_L):
                        breakpoint_depth = WaterDepth_L[idx]

            # Store results for this image
            # IMPORTANT: Store ARRAYS (not scalars) for matrix building
            result = {
                'timestamp': self.image_loader.timestamps[img_idx],
                'shoreline_mean': np.nanmean(shoreline) if shoreline is not None else np.nan,
                'wave_height': Hs,
                'wave_period_mean': Tm,
                'wave_period_peak': Tp,  # Already defined above
                # Store both arrays and scalars
                'wave_celerity_array': WaveCelerity,  # ARRAY for matrix building
                'wave_celerity': celerity_mean,  # Scalar for CSV
                'wavelength_array': WaveLength,  # ARRAY for matrix building
                'wavelength': wavelength_mean,  # Scalar for CSV
                'water_depth_array': WaterDepth_L,  # ARRAY for matrix building (WaterDepth_L)
                'water_depth': water_depth_mean,  # Scalar for CSV
                'wave_energy': wave_energy,
                'depth_shallow': depth_shallow,
                'roller_length': roller_length,
                'breakpoint_location': breakpoint_location,
                'breakpoint_depth': breakpoint_depth,
                'sla': sla,
                'rtr': np.nan,  # Will calculate after all timestacks processed
                # Store depth profile for backward compatibility
                'depth_profile': WaterDepth_L,
                'cross_shore_positions': bathymetry_results.get('cross_shore_positions', None),
                # Store average timestack profile for visualization
                'avg_timestack_profile': avg_timestack_profile,
                # Store raw Cf1 for SLA_S calculation
                'celerity_array_raw': Cf1,  # Raw Cf1 (before movmean transformation)
            }
            all_results.append(result)

            print(f"  → Wave Height: {result['wave_height']:.3f} m" if not np.isnan(result['wave_height']) else "  → Wave Height: N/A")
            print(f"  → Wave Period (Tm): {result['wave_period_mean']:.2f} s" if not np.isnan(result['wave_period_mean']) else "  → Wave Period: N/A")
            print(f"  → Wavelength: {result['wavelength']:.2f} m" if not np.isnan(result['wavelength']) else "  → Wavelength: N/A")
            print(f"  → Celerity: {result['wave_celerity']:.3f} m/s" if not np.isnan(result['wave_celerity']) else "  → Celerity: N/A")
            print(f"  → Water Depth: {result['water_depth']:.2f} m" if not np.isnan(result['water_depth']) else "  → Water Depth: N/A")

        # Aggregate all results
        self._aggregate_all_results(all_results)

        # Export results
        print("\nExporting results...")
        self._export_results()

        # Create visualizations
        print("\nCreating visualizations...")
        self._create_visualizations()

        print("\n" + "="*70)
        print("Analysis complete!")
        print(f"Results saved to: {self.config.output_dir}")
        print("="*70 + "\n")

    def _preprocess_single_timestack(self, img):
        """
        Preprocess a single timestack image.

        Parameters
        ----------
        img : np.ndarray
            Input timestack image (time x space x channels), e.g., (1680, 689, 3)

        Returns
        -------
        np.ndarray
            Preprocessed timestack (time x space), e.g., (1680, 689)
        """
        # Convert to grayscale if needed
        if len(img.shape) == 3:
            timestack = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY).astype(np.float32) / 255.0
        else:
            timestack = img.astype(np.float32) / 255.0

        # Enhance wave signals
        timestack_enhanced = self.preprocessor.enhance_waves(timestack)

        return timestack_enhanced

    def _analyze_single_timestack(self, timestack, shoreline):
        """
        Analyze wave parameters from a single timestack.

        Parameters
        ----------
        timestack : np.ndarray
            Preprocessed timestack (time x space), e.g., (1680, 689)
        shoreline : np.ndarray
            Detected shoreline positions

        Returns
        -------
        dict
            Wave analysis results
        """
        # Create cross-shore position array (spatial dimension)
        # Timestack shape: (time, space) = (1680, 689)
        num_positions = timestack.shape[1]  # Use SPACE dimension (689), not time (1680)
        cross_shore_positions = np.arange(num_positions) * self.config.pixel_resolution

        # Analyze this timestack
        wave_results = self.wave_analyzer.analyze_timestack(
            timestack, cross_shore_positions, shorelines=[shoreline] if shoreline is not None else None
        )

        return wave_results

    def _compute_sla_single(self, shoreline, bathymetry_results):
        """
        Compute sea level anomaly for a single timestack.

        Parameters
        ----------
        shoreline : np.ndarray
            Detected shoreline positions
        bathymetry_results : dict
            Bathymetry estimation results

        Returns
        -------
        float
            Sea level anomaly (meters)
        """
        if shoreline is not None:
            valid_pos = shoreline[~np.isnan(shoreline)]
            if len(valid_pos) > 0:
                avg_shoreline = np.mean(valid_pos)
            else:
                return np.nan
        else:
            return np.nan

        # Get depth
        mean_depth = bathymetry_results.get('mean_depth', 2.0)
        if np.isnan(mean_depth):
            mean_depth = 2.0  # Default

        # Compute SLA
        sla = self.sea_level_analyzer.compute_sla_from_shoreline(
            np.array([avg_shoreline]), np.array([mean_depth])
        )

        return sla[0] if len(sla) > 0 else np.nan

    def _aggregate_all_results(self, all_results):
        """
        Aggregate results from all individual timstacks.

        Parameters
        ----------
        all_results : list
            List of result dictionaries, one per timestack image
        """
        from scipy.ndimage import uniform_filter1d

        # Extract arrays for each parameter
        timestamps = [r['timestamp'] for r in all_results]
        shorelines = [r['shoreline_mean'] for r in all_results]
        wave_heights = [r['wave_height'] for r in all_results]
        wave_periods_mean = [r['wave_period_mean'] for r in all_results]
        wave_periods_peak = [r['wave_period_peak'] for r in all_results]
        wave_celerities = [r['wave_celerity'] for r in all_results]
        wavelengths = [r['wavelength'] for r in all_results]
        wave_energies = [r['wave_energy'] for r in all_results]
        depths_shallow = [r['depth_shallow'] for r in all_results]
        roller_lengths = [r['roller_length'] for r in all_results]
        breakpoint_locations = [r['breakpoint_location'] for r in all_results]
        breakpoint_depths = [r['breakpoint_depth'] for r in all_results]

        # ===== MATLAB Depth Smoothing Approach =====
        # In MATLAB (S01_AnalysisTimestackImages.m, lines 252-264):
        # 1. WaterDepth_L = [WaterDepth_L; movmean(df(:,1:numel(sLimit)), 10)];
        # 2. Csmooth_L = smooth2(WaterDepth_L, Nr=1, Nc=30);
        # 3. SLA_S = Csmooth_L - nanmean(Csmooth_L);

        # Extract raw depth profiles (not pre-smoothed)
        depth_profiles_list = []
        cross_shore_pos_list = []

        for r in all_results:
            if r.get('depth_profile') is not None:
                # Get the filtered depths (not yet smoothed across timestacks)
                depth_profiles_list.append(r['depth_profile'])
                cross_shore_pos_list.append(r.get('cross_shore_positions'))

        if len(depth_profiles_list) > 0:
            # Stack into 2D matrix: (num_timestacks × num_spatial_positions)
            min_length = min(len(profile) for profile in depth_profiles_list)
            depth_profiles_trimmed = [profile[:min_length] for profile in depth_profiles_list]
            depth_matrix = np.array(depth_profiles_trimmed)  # Shape: (timestacks, space)

            reference_positions = cross_shore_pos_list[0][:min_length]

            print(f"Depth matrix shape: {depth_matrix.shape} (timestacks × spatial points)")

            # Step 1: Apply 10-point moving average to each row (each timestack)
            # This matches MATLAB's movmean(df, 10)
            depth_matrix_ma = np.zeros_like(depth_matrix)
            for i in range(depth_matrix.shape[0]):
                # Use uniform_filter1d for moving average (size=10)
                depth_matrix_ma[i, :] = uniform_filter1d(
                    depth_matrix[i, :], size=10, mode='nearest'
                )

            # Step 2: Apply 2D smoothing (smooth2 with Nr=1, Nc=30)
            # Nr=1 means minimal smoothing across timestacks
            # Nc=30 means 30-point smoothing spatially
            from scipy.ndimage import uniform_filter

            # smooth2 with (Nr=1, Nc=30) means window (2*1+1, 2*30+1) = (3, 61)
            # But MATLAB's smooth2 is more like (1, 30) based on the parameters
            # Let's apply spatial smoothing with window size 30
            depth_matrix_smooth = np.zeros_like(depth_matrix_ma)
            for i in range(depth_matrix_ma.shape[0]):
                depth_matrix_smooth[i, :] = uniform_filter1d(
                    depth_matrix_ma[i, :], size=30, mode='nearest'
                )

            print(f"Applied MATLAB smoothing: movmean(10) + smooth2(Nr=1, Nc=30)")

            # CRITICAL FIX: Compute SLA from the 2D matrix FIRST (matches MATLAB line 266-267)
            # MATLAB: SLA_L = Csmooth_L - nanmean(Csmooth_L)
            # This preserves the spatial structure during SLA calculation
            sla_matrix = depth_matrix_smooth - np.nanmean(depth_matrix_smooth)
            print(f"SLA matrix shape: {sla_matrix.shape}")
            print(f"  SLA 2D range: {np.nanmin(sla_matrix):.3f} to {np.nanmax(sla_matrix):.3f} m")

            # Extract per-timestack SLA values by averaging across space
            # This gives one SLA value per timestack for the CSV
            slas = np.nanmean(sla_matrix, axis=1)  # Shape: (num_timestacks,)

            # Also compute mean depth per timestack (for CSV output)
            water_depths_array = np.nanmean(depth_matrix_smooth, axis=1)

            # Compute averaged bathymetry profile across all timestacks
            averaged_depth_profile = np.nanmean(depth_matrix_smooth, axis=0)
            bathymetry_positions = reference_positions

            # Calculate beach slope from bathymetry gradient
            nearshore_length = max(10, min_length // 5)
            x_nearshore = bathymetry_positions[:nearshore_length]
            y_nearshore = averaged_depth_profile[:nearshore_length]

            valid_mask = ~np.isnan(y_nearshore)
            if np.sum(valid_mask) >= 2:
                slope_fit = np.polyfit(x_nearshore[valid_mask], y_nearshore[valid_mask], 1)
                beach_slope = abs(slope_fit[0])
            else:
                beach_slope = 0.02

            print(f"Averaged bathymetry profile: {len(averaged_depth_profile)} points")
            print(f"  Depth range: {np.nanmin(averaged_depth_profile):.2f} to {np.nanmax(averaged_depth_profile):.2f} m")
            print(f"  Beach slope from bathymetry: {beach_slope:.4f} ({beach_slope*100:.2f}%)")
            print(f"Per-timestack SLA values (spatially averaged): {slas}")
            print(f"  SLA 1D range: {np.nanmin(slas):.3f} to {np.nanmax(slas):.3f} m")
        else:
            # Fallback if no depth profiles available
            print("⚠️  WARNING: No depth profiles available, using fallback SLA calculation")
            print("   This may produce less accurate SLA values!")
            water_depths_array = np.array([r['water_depth'] for r in all_results])
            averaged_depth_profile = np.array([])
            bathymetry_positions = np.array([])
            beach_slope = 0.02
            sla_matrix = None

            # Compute SLA from mean depths (less accurate fallback)
            mean_water_depth = np.nanmean(water_depths_array)
            slas = water_depths_array - mean_water_depth

        print(f"Water depths (per timestack, spatially averaged): {water_depths_array}")
        print(f"Mean water depth across all: {np.nanmean(water_depths_array):.3f} m")

        # Calculate RTR (Relative Tidal Range) - MATLAB lines 274-278
        # Level_TS = nanmean(WaveCelerity, 2);
        # nRTR_thresh = Level_TS - min(Level_TS);
        # nRTR_thresh(nRTR_thresh < 0.2) = NaN; % threshold to avoid dividing by 0
        # RTR = hs./nRTR_thresh;
        water_levels = np.array(wave_celerities)  # Use celerity as water level proxy
        min_water_level = np.nanmin(water_levels)
        nRTR_thresh = water_levels - min_water_level
        nRTR_thresh[nRTR_thresh < 0.2] = np.nan  # Avoid division by near-zero
        rtrs = np.array(wave_heights) / nRTR_thresh
        print(f"RTR (Relative Tidal Range): mean={np.nanmean(rtrs):.3f}, range=[{np.nanmin(rtrs):.3f}, {np.nanmax(rtrs):.3f}]")

        # ===== Calculate SLA_S from WaveCelerity Matrix (MATLAB line 256-267) =====
        # MATLAB:
        # Csmooth_S = smooth2(WaveCelerity, Nr=1, Nc=30);
        # SLA_S = Csmooth_S - nanmean(Csmooth_S);
        wave_celerity_arrays = []
        for r in all_results:
            if r.get('wave_celerity_array') is not None:
                wave_celerity_arrays.append(r['wave_celerity_array'])

        if len(wave_celerity_arrays) > 0:
            # Build WaveCelerity matrix
            min_len = min(len(arr) for arr in wave_celerity_arrays)
            WaveCelerity_matrix = np.array([arr[:min_len] for arr in wave_celerity_arrays])
            print(f"\nBuilding WaveCelerity matrix for SLA_S:")
            print(f"  Matrix shape: {WaveCelerity_matrix.shape} (timestacks × spatial points)")
            print(f"  Celerity range: [{np.nanmin(WaveCelerity_matrix):.2f}, {np.nanmax(WaveCelerity_matrix):.2f}] m/s")

            # Apply smooth2 with Nc=30 (spatial smoothing)
            # MATLAB's smooth2(WaveCelerity, Nr=1, Nc=30) means:
            # - Nr=1: minimal smoothing across timestacks (3 rows)
            # - Nc=30: 30-point smoothing spatially (61 columns)
            # For simplicity, we'll apply spatial smoothing only (size=30)
            Csmooth_S = np.zeros_like(WaveCelerity_matrix)
            for i in range(WaveCelerity_matrix.shape[0]):
                Csmooth_S[i, :] = uniform_filter1d(
                    WaveCelerity_matrix[i, :], size=30, mode='nearest'
                )

            # Calculate SLA_S from smoothed matrix
            SLA_S_matrix = Csmooth_S - np.nanmean(Csmooth_S)
            print(f"  SLA_S 2D range: {np.nanmin(SLA_S_matrix):.3f} to {np.nanmax(SLA_S_matrix):.3f} m/s")

            # Extract per-timestack SLA_S values by averaging across space
            sla_shallow = np.nanmean(SLA_S_matrix, axis=1)  # Shape: (timestacks,)
            print(f"  SLA_S (spatially averaged): {sla_shallow}")
            print(f"  SLA_S 1D range: {np.nanmin(sla_shallow):.3f} to {np.nanmax(sla_shallow):.3f} m/s")
            print(f"SLA_S (shallow water): mean={np.nanmean(sla_shallow):.3f} m")
        else:
            sla_shallow = np.full(len(all_results), np.nan)

        # Build average timestack matrix (MATLAB line 281: Stack_av)
        # Extract average timestack profiles from all results
        avg_timestack_profiles = []
        for r in all_results:
            if 'avg_timestack_profile' in r:
                avg_timestack_profiles.append(r['avg_timestack_profile'])

        if len(avg_timestack_profiles) > 0:
            # Stack into 2D matrix: (num_timestacks × spatial_pixels)
            min_length = min(len(profile) for profile in avg_timestack_profiles)
            avg_timestack_trimmed = [profile[:min_length] for profile in avg_timestack_profiles]
            average_timestack = np.array(avg_timestack_trimmed)  # Shape: (timestacks, space)
            print(f"Average timestack matrix shape: {average_timestack.shape}")
        else:
            average_timestack = None
            print("Warning: No average timestack profiles available")

        # Store in results dictionary
        self.results = {
            'timestamps': timestamps,
            'shoreline_positions': np.array(shorelines),
            'wave_heights_timeseries': np.array(wave_heights),
            'wave_periods_timeseries': np.array(wave_periods_mean),
            'wave_periods_peak_timeseries': np.array(wave_periods_peak),
            'celerities': np.array(wave_celerities),
            'wavelengths': np.array(wavelengths),
            'wave_energies': np.array(wave_energies),
            'depths': water_depths_array,
            'depths_shallow': np.array(depths_shallow),
            'roller_lengths': np.array(roller_lengths),
            'breakpoint_locations': np.array(breakpoint_locations),
            'breakpoint_depths': np.array(breakpoint_depths),
            'sla_values': slas,  # SLA_L (linear wave theory)
            'sla_shallow_values': sla_shallow,  # SLA_S (shallow water)
            'rtr_values': rtrs,  # Relative Tidal Range
            'beach_slope': beach_slope,
            # Averaged bathymetry profile
            'depths_smoothed': averaged_depth_profile,
            'cross_shore_positions': bathymetry_positions,
            # Aggregate statistics
            'mean_Hs': np.nanmean(wave_heights),
            'mean_Tm': np.nanmean(wave_periods_mean),
            'mean_Tp': np.nanmean(wave_periods_peak),
            'mean_celerity': np.nanmean(wave_celerities),
            'mean_depth': np.nanmean(water_depths_array),
            # For MATLAB-style visualization
            'average_timestack': average_timestack,
            'sla_matrix': sla_matrix,
            'sla_shallow_matrix': sla_shallow_matrix if 'sla_shallow_matrix' in locals() else None,
        }

    def _detect_shorelines(self):
        """Detect shorelines in all images."""
        images = self.image_loader.load_all_images()
        shorelines = []

        for i, img in enumerate(tqdm(images, desc="Detecting shorelines")):
            shoreline = self.shoreline_detector.detect(img)
            shorelines.append(shoreline)

        return shorelines

    def _create_timestack(self):
        """Create timestack array from images."""
        images = self.image_loader.load_all_images()

        # Determine cross-shore range
        if self.config.max_cross_shore is None:
            max_pixel = min(img.shape[1] for img in images)
        else:
            max_pixel = self.config.max_cross_shore

        cross_shore_range = (self.config.min_cross_shore, max_pixel)

        # Create timestack
        timestack = self.preprocessor.create_timestack_array(images, cross_shore_range)

        # Enhance wave signals
        timestack_enhanced = self.preprocessor.enhance_waves(timestack)

        return timestack_enhanced

    def _analyze_waves(self, timestack, shorelines):
        """Analyze wave parameters from timestack."""
        # Create cross-shore position array
        num_positions = timestack.shape[0]
        cross_shore_positions = np.arange(num_positions) * self.config.pixel_resolution

        # Analyze timestack with shoreline data for improved wave height estimation
        wave_results = self.wave_analyzer.analyze_timestack(
            timestack, cross_shore_positions, shorelines=shorelines
        )

        return wave_results

    def _analyze_correlations(self, timestack):
        """Perform cross-correlation analysis."""
        correlation_results = self.correlation_analyzer.analyze_timestack(timestack)

        return correlation_results

    def _estimate_bathymetry(self, wave_results, correlation_results):
        """Estimate bathymetry from wave data across full spatial domain."""
        # Extract full spatial grid
        cross_shore_positions_full = wave_results.get('cross_shore_positions', np.array([]))
        num_positions = len(cross_shore_positions_full)

        # Get celerities from MATLAB-style correlation analysis
        if 'celerities' in correlation_results and len(correlation_results['celerities']) > 0:
            celerities_full = correlation_results['celerities']

            # MATLAB-style correlation returns full array (with NaN at edges)
            # If length doesn't match, pad or trim
            if len(celerities_full) < num_positions:
                # Pad with NaN
                padded = np.full(num_positions, np.nan)
                padded[:len(celerities_full)] = celerities_full
                celerities_full = padded
            elif len(celerities_full) > num_positions:
                # Trim
                celerities_full = celerities_full[:num_positions]

            # Count valid celerities
            num_valid = np.sum(~np.isnan(celerities_full))
            print(f"  Celerities: {num_valid}/{num_positions} valid values")

            # Interpolate NaN values if we have enough valid points
            if num_valid >= 2:
                from scipy.interpolate import interp1d
                valid_mask = ~np.isnan(celerities_full)
                valid_indices = np.where(valid_mask)[0]
                valid_celerities = celerities_full[valid_mask]

                # Interpolate (but don't extrapolate too far)
                interp_func = interp1d(valid_indices, valid_celerities,
                                      kind='linear', bounds_error=False,
                                      fill_value=(valid_celerities[0], valid_celerities[-1]))
                all_indices = np.arange(num_positions)
                celerities_full = interp_func(all_indices)
                celerities_full = np.abs(celerities_full)  # Ensure positive
            elif num_valid > 0:
                # Use mean of valid values
                mean_celerity = np.nanmean(celerities_full)
                print(f"  Warning: Only {num_valid} valid celerities, using mean={mean_celerity:.2f} m/s")
                celerities_full = np.full(num_positions, mean_celerity)
            else:
                # No valid celerities - use default
                print(f"  Warning: No valid celerities found, using default=5.0 m/s")
                celerities_full = np.full(num_positions, 5.0)
        else:
            # Estimate celerities from wave periods if correlation failed
            print(f"  Warning: No celerities from correlation, using default=5.0 m/s")
            celerities_full = np.full(num_positions, 5.0)

        # Use PEAK wave period for all positions (matching MATLAB line 247)
        # MATLAB: [df] = LinearC(Tp, WaveCelerity(i,:), 0.01);
        # Tp is the PEAK period (scalar), not mean period
        peak_period = wave_results.get('mean_Tp', wave_results.get('mean_Tm', 8.0))
        print(f"  Using peak period Tp = {peak_period:.2f}s for depth calculation")
        wave_periods_full = np.full(num_positions, peak_period)

        # Estimate depth profile across full spatial domain
        bathymetry_results = self.bathymetry_estimator.estimate_depth_profile(
            wave_periods_full, celerities_full, cross_shore_positions_full
        )

        return bathymetry_results

    def _compute_sea_level(self, shorelines, bathymetry_results):
        """Compute sea level anomalies."""
        # Get average shoreline positions
        avg_shorelines = []

        for shoreline in shorelines:
            if shoreline is not None:
                valid_pos = shoreline[~np.isnan(shoreline)]
                if len(valid_pos) > 0:
                    avg_shorelines.append(np.mean(valid_pos))
                else:
                    avg_shorelines.append(np.nan)
            else:
                avg_shorelines.append(np.nan)

        avg_shorelines = np.array(avg_shorelines)

        # Get depths (use mean depth for all positions)
        if 'depths_smoothed' in bathymetry_results:
            depths = bathymetry_results['depths_smoothed']
            mean_depth = np.nanmean(depths)
        else:
            mean_depth = 2.0  # Default

        depths_array = np.full(len(avg_shorelines), mean_depth)

        # Compute SLA
        sla = self.sea_level_analyzer.compute_sla_from_shoreline(
            avg_shorelines, depths_array
        )

        sla_results = {
            'sla_values': sla,
            'shoreline_positions': avg_shorelines,
            'statistics': self.sea_level_analyzer.compute_sla_statistics(sla)
        }

        return sla_results

    def _aggregate_results(self, shorelines, wave_results, correlation_results,
                          bathymetry_results, sla_results, timestack):
        """Aggregate all results into main results dictionary."""
        self.results = {
            'shorelines': shorelines,
            'timestack': timestack,
            **wave_results,
            **correlation_results,
            **bathymetry_results,
            **sla_results
        }

        # Add timestamps
        self.results['timestamps'] = self.image_loader.timestamps

    def _export_results(self):
        """Export results to files matching MATLAB timetable format."""
        timestamps = self.results.get('timestamps', [])

        if len(timestamps) == 0:
            print("Warning: No timestamps available for export")
            return

        # Create DataFrame matching MATLAB timetable columns:
        # {'BreakPointLocation', 'BreakpointDepth', 'Hs', 'WaveEnergy', 'RollerLength',
        #  'WaveCelerity', 'Tp', 'Tm', 'WaveLength', 'WaterDepth', 'ShorelinePosition',
        #  'SLA_S', 'SLA_L', 'RTR', 'Bathymetry'}
        data = {
            'Timestamp': timestamps,
            'BreakPointLocation': self.results.get('breakpoint_locations', []),
            'BreakpointDepth': self.results.get('breakpoint_depths', []),
            'Hs': self.results.get('wave_heights_timeseries', []),
            'WaveEnergy': self.results.get('wave_energies', []),
            'RollerLength': self.results.get('roller_lengths', []),
            'WaveCelerity': self.results.get('celerities', []),
            'Tp': self.results.get('wave_periods_peak_timeseries', []),
            'Tm': self.results.get('wave_periods_timeseries', []),
            'WaveLength': self.results.get('wavelengths', []),
            'WaterDepth': self.results.get('depths', []),
            'ShorelinePosition': self.results.get('shoreline_positions', []),
            'SLA_S': self.results.get('sla_shallow_values', []),
            'SLA_L': self.results.get('sla_values', []),
            'RTR': self.results.get('rtr_values', []),
            'Bathymetry': self.results.get('depths', []),  # Same as WaterDepth
        }

        df = pd.DataFrame(data)

        # Set Timestamp as index to create timetable (like MATLAB)
        df.set_index('Timestamp', inplace=True)

        # Export to CSV
        output_file = os.path.join(
            self.config.output_dir,
            f"coastcams_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
        )

        df.to_csv(output_file, index=False)
        print(f"Results exported to: {output_file}")

        # Create summary report
        self.visualizer.create_summary_report(self.results, self.config)

    def _create_visualizations(self):
        """Create all visualizations."""
        print("\n=== Creating Visualizations ===")
        timestamps = self.results.get('timestamps', [])
        print(f"Number of timestamps: {len(timestamps)}")

        # Plot MATLAB-style summary matching plot_coastcams_main(Time_TS, Stack_av, SLA_S, Hs_TS, Tp_TS, rotation)
        if ('average_timestack' in self.results and
            'wave_heights_timeseries' in self.results and
            len(timestamps) > 0):
            print("Plotting MATLAB-style summary...")
            # Use SLA_S (shallow water SLA) to match MATLAB
            sla_for_plot = self.results.get('sla_shallow_matrix', self.results.get('sla_matrix', None))
            # Get water levels (mean sea level / water depth) for plotting
            water_levels = self.results.get('depths', None)
            self.visualizer.plot_matlab_style_summary(
                timestamps=timestamps,
                average_timestack=self.results['average_timestack'],
                sla_matrix=sla_for_plot,
                wave_heights=self.results['wave_heights_timeseries'],
                wave_periods=self.results['wave_periods_timeseries'],
                water_levels=water_levels,
                rotation=self.config.rotation_angle
            )

        # Plot timestack
        if 'timestack' in self.results:
            print("Plotting timestack...")
            self.visualizer.plot_timestack(self.results['timestack'])

        # Plot wave parameters
        if 'wave_heights_timeseries' in self.results and len(timestamps) > 0:
            print("Plotting wave parameters...")
            wave_heights = self.results['wave_heights_timeseries']
            wave_periods = self.results['wave_periods_timeseries']
            wave_celerities = self.results.get('celerities', None)

            self.visualizer.plot_wave_parameters(timestamps, wave_heights, wave_periods, wave_celerities)

        # Plot bathymetry
        if 'depths_smoothed' in self.results and 'cross_shore_positions' in self.results:
            print("Plotting bathymetry...")
            self.visualizer.plot_bathymetry(
                self.results['cross_shore_positions'],
                self.results['depths_smoothed']
            )

        # Plot comprehensive analysis
        print("Plotting comprehensive analysis...")
        self.visualizer.plot_comprehensive_analysis(self.results, timestamps)

        print("=== Visualization Complete ===\n")


def main():
    """Main entry point for CoastCams analysis."""
    parser = argparse.ArgumentParser(
        description='CoastCams - Automated Coastal Wave Analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        '--config', '-c',
        type=str,
        default=None,
        help='Path to configuration file (YAML)'
    )

    parser.add_argument(
        '--create-config',
        action='store_true',
        help='Create a template configuration file'
    )

    args = parser.parse_args()

    # Create template config if requested
    if args.create_config:
        CoastCamsConfig.create_template()
        return

    # Run analysis
    try:
        # Use config.yaml by default if no config specified
        config_file = args.config
        if config_file is None:
            default_config = os.path.join(os.path.dirname(__file__), 'config.yaml')
            if os.path.exists(default_config):
                config_file = default_config

        workflow = CoastCamsWorkflow(config_file)
        workflow.run_full_analysis()
    except KeyboardInterrupt:
        print("\n\nAnalysis interrupted by user.")
        sys.exit(0)
    except Exception as e:
        print(f"\nError during analysis: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
