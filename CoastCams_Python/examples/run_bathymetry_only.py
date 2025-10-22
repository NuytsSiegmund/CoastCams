#!/usr/bin/env python3
"""
Example: Bathymetry Estimation Only

This script demonstrates how to estimate bathymetry (water depth)
from wave parameters.
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
from coastcams.cross_correlation import CrossCorrelationAnalyzer
from coastcams.bathymetry import BathymetryEstimator
from coastcams.visualize import CoastCamsVisualizer


def main():
    """Run bathymetry estimation example."""

    print("="*70)
    print("CoastCams - Bathymetry Estimation Example")
    print("="*70)

    # Initialize configuration from config.yaml
    config_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'config.yaml')
    config = CoastCamsConfig(config_path)

    # Initialize modules
    loader = ImageLoader(config.input_dir, config.rotation_angle, verbose=False)
    preprocessor = ImagePreprocessor(config)
    wave_analyzer = WaveAnalyzer(config)
    correlation_analyzer = CrossCorrelationAnalyzer(config)
    bathymetry_estimator = BathymetryEstimator(config)
    visualizer = CoastCamsVisualizer(config.output_dir)

    # Load and process images
    print("\nLoading and processing images...")
    num_images = loader.discover_images()

    if num_images == 0:
        print("Error: No images found!")
        return

    images = loader.load_all_images()

    # Create timestack
    cross_shore_range = (config.min_cross_shore, config.cross_shore_width)
    timestack = preprocessor.create_timestack_array(images, cross_shore_range)
    timestack_enhanced = preprocessor.enhance_waves(timestack)

    # Analyze waves and correlations
    print("Analyzing wave parameters...")
    cross_shore_positions = np.arange(timestack.shape[0]) * config.pixel_resolution

    wave_results = wave_analyzer.analyze_timestack(timestack_enhanced,
                                                   cross_shore_positions)

    print("Performing cross-correlation analysis...")
    correlation_results = correlation_analyzer.analyze_timestack(timestack_enhanced)

    # Estimate bathymetry
    print("\nEstimating bathymetry...")

    wave_periods = wave_results.get('wave_periods', np.array([]))

    if 'celerities' in correlation_results:
        celerities = correlation_results['celerities']

        # Match lengths
        min_len = min(len(wave_periods), len(celerities))
        wave_periods = wave_periods[:min_len]
        celerities = celerities[:min_len]
        positions = cross_shore_positions[:min_len]
    else:
        # Use default celerities if correlation failed
        print("  Warning: Using estimated celerities")
        positions = cross_shore_positions
        celerities = np.array([5.0] * len(wave_periods))

    # Estimate depth profile using linear wave theory
    print("  Using linear wave theory...")
    bathymetry_linear = bathymetry_estimator.estimate_depth_profile(
        wave_periods, celerities, positions, method='linear'
    )

    # Estimate using shallow water approximation
    print("  Using shallow water approximation...")
    bathymetry_shallow = bathymetry_estimator.estimate_depth_profile(
        wave_periods, celerities, positions, method='shallow'
    )

    # Display results
    print("\n" + "="*70)
    print("Bathymetry Results")
    print("="*70)

    print("\nLinear Wave Theory:")
    stats_linear = bathymetry_estimator.compute_depth_statistics(
        bathymetry_linear['depths_smoothed']
    )
    print(f"  Mean depth: {stats_linear['mean']:.2f} m")
    print(f"  Max depth: {stats_linear['max']:.2f} m")
    print(f"  Min depth: {stats_linear['min']:.2f} m")
    print(f"  Depth range: {stats_linear['range']:.2f} m")

    print("\nShallow Water Approximation:")
    stats_shallow = bathymetry_estimator.compute_depth_statistics(
        bathymetry_shallow['depths_smoothed']
    )
    print(f"  Mean depth: {stats_shallow['mean']:.2f} m")
    print(f"  Max depth: {stats_shallow['max']:.2f} m")
    print(f"  Min depth: {stats_shallow['min']:.2f} m")

    # Compute shore slope
    print("\nShore Slope Analysis:")
    slope_nearshore = bathymetry_estimator.compute_shore_slope(
        bathymetry_linear['depths_smoothed'],
        positions,
        region='nearshore'
    )

    if not np.isnan(slope_nearshore):
        print(f"  Nearshore slope: {slope_nearshore:.4f} (1:{abs(1/slope_nearshore):.1f})")

    # Visualize bathymetry
    print("\nCreating bathymetry visualizations...")

    # Plot both methods for comparison
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 1, figsize=(12, 10))

    # Linear wave theory
    axes[0].plot(positions, -bathymetry_linear['depths_smoothed'],
                'b-', linewidth=2, label='Estimated')
    axes[0].axhline(y=0, color='cyan', linestyle='--', linewidth=1.5,
                   alpha=0.7, label='Water Surface')
    axes[0].set_xlabel('Cross-shore Distance (m)')
    axes[0].set_ylabel('Elevation (m)')
    axes[0].set_title('Bathymetry - Linear Wave Theory')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    axes[0].invert_yaxis()

    # Shallow water approximation
    axes[1].plot(positions, -bathymetry_shallow['depths_smoothed'],
                'g-', linewidth=2, label='Estimated')
    axes[1].axhline(y=0, color='cyan', linestyle='--', linewidth=1.5,
                   alpha=0.7, label='Water Surface')
    axes[1].set_xlabel('Cross-shore Distance (m)')
    axes[1].set_ylabel('Elevation (m)')
    axes[1].set_title('Bathymetry - Shallow Water Approximation')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    axes[1].invert_yaxis()

    plt.tight_layout()
    plt.savefig(os.path.join(config.output_dir, 'bathymetry_comparison.png'),
               dpi=300, bbox_inches='tight')
    plt.close()

    # Create detailed bathymetry plot
    visualizer.plot_bathymetry(
        positions,
        bathymetry_linear['depths_smoothed'],
        title='Estimated Bathymetry Profile',
        filename='bathymetry_detailed.png'
    )

    print("\n" + "="*70)
    print("Bathymetry estimation complete!")
    print(f"Results saved to: {config.output_dir}")
    print("="*70)


if __name__ == '__main__':
    main()
