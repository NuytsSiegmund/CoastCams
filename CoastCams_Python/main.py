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

            # Step 4: Analyze wave parameters
            print("[4/7] Analyzing wave parameters...")
            wave_results = self._analyze_single_timestack(timestack, shoreline)

            # Step 5: Perform cross-correlation analysis
            print("[5/7] Performing cross-correlation...")
            correlation_results = self.correlation_analyzer.analyze_timestack(timestack)

            # Step 6: Estimate bathymetry
            print("[6/7] Estimating bathymetry...")
            bathymetry_results = self._estimate_bathymetry(wave_results, correlation_results)

            # Step 7: Compute sea level anomaly
            print("[7/7] Computing sea level anomaly...")
            sla = self._compute_sla_single(shoreline, bathymetry_results)

            # Store results for this image
            result = {
                'timestamp': self.image_loader.timestamps[img_idx],
                'shoreline_mean': np.nanmean(shoreline) if shoreline is not None else np.nan,
                'wave_height': wave_results.get('mean_Hs', np.nan),
                'wave_period': wave_results.get('mean_Tm', np.nan),
                'wave_celerity': correlation_results.get('mean_celerity', np.nan),
                'water_depth': bathymetry_results.get('mean_depth', np.nan),
                'sla': sla,
            }
            all_results.append(result)

            print(f"  → Wave Height: {result['wave_height']:.3f} m" if not np.isnan(result['wave_height']) else "  → Wave Height: N/A")
            print(f"  → Wave Period: {result['wave_period']:.2f} s" if not np.isnan(result['wave_period']) else "  → Wave Period: N/A")
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
        # Create cross-shore position array
        num_positions = timestack.shape[0]
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
        # Extract arrays for each parameter
        timestamps = [r['timestamp'] for r in all_results]
        shorelines = [r['shoreline_mean'] for r in all_results]
        wave_heights = [r['wave_height'] for r in all_results]
        wave_periods = [r['wave_period'] for r in all_results]
        wave_celerities = [r['wave_celerity'] for r in all_results]
        water_depths = [r['water_depth'] for r in all_results]
        slas = [r['sla'] for r in all_results]

        # Store in results dictionary
        self.results = {
            'timestamps': timestamps,
            'shoreline_positions': np.array(shorelines),
            'wave_heights_timeseries': np.array(wave_heights),
            'wave_periods_timeseries': np.array(wave_periods),
            'celerities': np.array(wave_celerities),
            'depths': np.array(water_depths),
            'sla_values': np.array(slas),
            # Aggregate statistics
            'mean_Hs': np.nanmean(wave_heights),
            'mean_Tm': np.nanmean(wave_periods),
            'mean_celerity': np.nanmean(wave_celerities),
            'mean_depth': np.nanmean(water_depths),
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
        """Estimate bathymetry from wave data."""
        # Extract wave parameters
        wave_periods = wave_results.get('wave_periods', np.array([]))
        cross_shore_positions = wave_results.get('cross_shore_positions', np.array([]))

        # Get celerities from correlation analysis
        if 'celerities' in correlation_results:
            celerities = correlation_results['celerities']

            # Match lengths
            min_len = min(len(wave_periods), len(celerities))
            wave_periods = wave_periods[:min_len]
            celerities = celerities[:min_len]
            positions = cross_shore_positions[:min_len]
        else:
            # Estimate celerities from wave periods
            positions = cross_shore_positions
            celerities = np.array([5.0] * len(wave_periods))  # Default estimate

        # Estimate depth profile
        bathymetry_results = self.bathymetry_estimator.estimate_depth_profile(
            wave_periods, celerities, positions
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
        """Export results to files."""
        timestamps = self.results.get('timestamps', [])

        if len(timestamps) == 0:
            print("Warning: No timestamps available for export")
            return

        # Create DataFrame
        data = {
            'Timestamp': timestamps,
        }

        # Add shoreline positions
        if 'shoreline_positions' in self.results:
            data['ShorelinePosition_pixels'] = self.results['shoreline_positions']
            data['ShorelinePosition_m'] = (self.results['shoreline_positions'] *
                                          self.config.pixel_resolution)

        # Add wave parameters (per-image values)
        if 'wave_heights_timeseries' in self.results:
            data['WaveHeight_m'] = self.results['wave_heights_timeseries']

        if 'wave_periods_timeseries' in self.results:
            data['WavePeriod_s'] = self.results['wave_periods_timeseries']

        if 'celerities' in self.results:
            data['WaveCelerity_m_s'] = self.results['celerities']

        if 'depths' in self.results:
            data['WaterDepth_m'] = self.results['depths']

        # Add SLA
        if 'sla_values' in self.results:
            data['SeaLevelAnomaly_m'] = self.results['sla_values']

        df = pd.DataFrame(data)

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

        # Plot timestack
        if 'timestack' in self.results:
            print("Plotting timestack...")
            self.visualizer.plot_timestack(self.results['timestack'])

        # Plot wave parameters
        if 'mean_Hs' in self.results and len(timestamps) > 0:
            print("Plotting wave parameters...")
            wave_heights = np.full(len(timestamps), self.results['mean_Hs'])
            wave_periods = np.full(len(timestamps), self.results.get('mean_Tm', np.nan))

            self.visualizer.plot_wave_parameters(timestamps, wave_heights, wave_periods)

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
