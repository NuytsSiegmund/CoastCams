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

        # Step 2: Detect shorelines
        print("\n[2/7] Detecting shorelines...")
        shorelines = self._detect_shorelines()

        # Step 3: Preprocess images and create timestack
        print("\n[3/7] Preprocessing images and creating timestack...")
        timestack = self._create_timestack()

        # Step 4: Analyze wave parameters
        print("\n[4/7] Analyzing wave parameters...")
        wave_results = self._analyze_waves(timestack, shorelines)

        # Step 5: Perform cross-correlation analysis
        print("\n[5/7] Performing cross-correlation analysis...")
        correlation_results = self._analyze_correlations(timestack)

        # Step 6: Estimate bathymetry
        print("\n[6/7] Estimating bathymetry...")
        bathymetry_results = self._estimate_bathymetry(wave_results, correlation_results)

        # Step 7: Compute sea level anomalies
        print("\n[7/7] Computing sea level anomalies...")
        sla_results = self._compute_sea_level(shorelines, bathymetry_results)

        # Aggregate all results
        self._aggregate_results(shorelines, wave_results, correlation_results,
                               bathymetry_results, sla_results, timestack)

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

        # Add wave parameters
        # Use time-varying wave heights if available
        if 'wave_heights_timeseries' in self.results:
            wave_heights_ts = self.results['wave_heights_timeseries']
            # Ensure lengths match
            if len(wave_heights_ts) == len(timestamps):
                data['WaveHeight_m'] = wave_heights_ts
            else:
                print(f"Warning: Wave heights length ({len(wave_heights_ts)}) != timestamps ({len(timestamps)})")
                data['WaveHeight_m'] = [self.results.get('mean_Hs', np.nan)] * len(timestamps)
        elif 'mean_Hs' in self.results:
            data['WaveHeight_m'] = [self.results['mean_Hs']] * len(timestamps)

        # Wave period - currently computed as aggregate, could be per-timestep in future
        if 'mean_Tm' in self.results and not np.isnan(self.results['mean_Tm']):
            data['WavePeriod_s'] = [self.results['mean_Tm']] * len(timestamps)
        else:
            data['WavePeriod_s'] = [np.nan] * len(timestamps)

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
