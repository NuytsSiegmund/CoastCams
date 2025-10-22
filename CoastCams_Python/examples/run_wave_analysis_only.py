#!/usr/bin/env python3
"""
Example: Wave Analysis Only

This script demonstrates how to run only wave parameter analysis
without full workflow or bathymetry estimation.
"""

import sys
import os
import numpy as np

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from coastcams.config import CoastCamsConfig
from coastcams.image_loader import ImageLoader
from coastcams.preprocessing import ImagePreprocessor
from coastcams.wave_analysis import WaveAnalyzer
from coastcams.visualize import CoastCamsVisualizer


def main():
    """Run wave analysis example."""

    print("="*70)
    print("CoastCams - Wave Analysis Example")
    print("="*70)

    # Initialize configuration from config.yaml
    config_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'config.yaml')
    config = CoastCamsConfig(config_path)

    # Initialize modules
    loader = ImageLoader(config.input_dir, config.rotation_angle)
    preprocessor = ImagePreprocessor(config)
    wave_analyzer = WaveAnalyzer(config)
    visualizer = CoastCamsVisualizer(config.output_dir)

    # Load images
    print("\nLoading images...")
    num_images = loader.discover_images()

    if num_images == 0:
        print("Error: No images found!")
        return

    print(f"Found {num_images} images")
    images = loader.load_all_images()

    # Create timestack
    print("\nCreating timestack...")
    cross_shore_range = (config.min_cross_shore, config.cross_shore_width)
    timestack = preprocessor.create_timestack_array(images, cross_shore_range)

    print(f"Timestack shape: {timestack.shape}")
    print(f"  Spatial points: {timestack.shape[0]}")
    print(f"  Time points: {timestack.shape[1]}")

    # Enhance wave signals
    print("\nEnhancing wave signals...")
    timestack_enhanced = preprocessor.enhance_waves(timestack)

    # Visualize timestack
    visualizer.plot_timestack(
        timestack_enhanced,
        timestamps=loader.timestamps,
        title='Enhanced Timestack'
    )

    # Analyze wave parameters
    print("\nAnalyzing wave parameters...")
    cross_shore_positions = np.arange(timestack.shape[0]) * config.pixel_resolution

    wave_results = wave_analyzer.analyze_timestack(
        timestack_enhanced,
        cross_shore_positions
    )

    # Display results
    print("\n" + "="*70)
    print("Wave Analysis Results")
    print("="*70)

    if 'mean_Hs' in wave_results:
        print(f"\nSignificant Wave Height:")
        print(f"  Mean Hs: {wave_results['mean_Hs']:.3f} m")
        print(f"  Max Hs: {wave_results['max_Hs']:.3f} m")

    if 'mean_Tm' in wave_results:
        print(f"\nWave Period:")
        print(f"  Mean Tm: {wave_results['mean_Tm']:.3f} s")

    # Analyze individual time series at specific location
    print("\n\nDetailed Analysis at Breaking Point:")

    # Find breaking point (max wave height location)
    if 'wave_heights' in wave_results:
        wave_heights = wave_results['wave_heights']
        break_idx = np.nanargmax(wave_heights)
        break_position = cross_shore_positions[break_idx]

        print(f"  Breaking location: {break_position:.2f} m")
        print(f"  Breaking wave height: {wave_heights[break_idx]:.3f} m")

        # Analyze timeseries at this location
        timeseries = timestack_enhanced[break_idx, :]

        # Compute spectral properties
        spectral_moments = wave_analyzer.compute_spectral_moments(timeseries)

        print(f"\nSpectral Analysis:")
        print(f"  Spectral Hs: {spectral_moments['significant_height']:.3f} m")
        print(f"  Zeroth moment (m0): {spectral_moments['m0']:.6f}")

        # Analyze wave groups
        group_analysis = wave_analyzer.analyze_wave_groups(timeseries)

        print(f"\nWave Grouping:")
        print(f"  Number of groups: {group_analysis['num_groups']}")
        if not np.isnan(group_analysis['mean_group_period']):
            print(f"  Mean group period: {group_analysis['mean_group_period']:.2f} s")

    # Create visualization
    if 'wave_heights' in wave_results:
        print("\nCreating wave parameter plots...")

        # Create plots along cross-shore
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 1, figsize=(12, 8))

        # Wave height profile
        axes[0].plot(cross_shore_positions, wave_results['wave_heights'],
                    'b-', linewidth=2)
        axes[0].set_xlabel('Cross-shore Distance (m)')
        axes[0].set_ylabel('Wave Height (m)')
        axes[0].set_title('Wave Height Profile')
        axes[0].grid(True, alpha=0.3)

        # Wave period profile
        axes[1].plot(cross_shore_positions, wave_results['wave_periods'],
                    'g-', linewidth=2)
        axes[1].set_xlabel('Cross-shore Distance (m)')
        axes[1].set_ylabel('Wave Period (s)')
        axes[1].set_title('Wave Period Profile')
        axes[1].grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(os.path.join(config.output_dir, 'wave_profiles.png'),
                   dpi=300, bbox_inches='tight')
        plt.close()

    print("\n" + "="*70)
    print("Wave analysis complete!")
    print(f"Results saved to: {config.output_dir}")
    print("="*70)


if __name__ == '__main__':
    main()
