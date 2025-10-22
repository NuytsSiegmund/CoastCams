#!/usr/bin/env python3
"""
Example: Shoreline Detection Only

This script demonstrates how to run only shoreline detection
without performing full wave analysis.
"""

import sys
import os
import numpy as np

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from coastcams.config import CoastCamsConfig
from coastcams.image_loader import ImageLoader
from coastcams.shoreline import ShorelineDetector
from coastcams.visualize import CoastCamsVisualizer


def main():
    """Run shoreline detection example."""

    print("="*70)
    print("CoastCams - Shoreline Detection Example")
    print("="*70)

    # Initialize configuration from config.yaml
    config_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'config.yaml')
    config = CoastCamsConfig(config_path)

    print("\nConfiguration:")
    print(f"  Input directory: {config.input_dir}")
    print(f"  Shoreline method: {config.shoreline_method}")
    print(f"  Threshold: {config.shoreline_threshold}")

    # Initialize image loader
    print("\nLoading images...")
    loader = ImageLoader(config.input_dir, config.rotation_angle)
    num_images = loader.discover_images()

    if num_images == 0:
        print("Error: No images found!")
        return

    print(f"Found {num_images} images")

    # Initialize shoreline detector
    detector = ShorelineDetector(config=config)

    # Test different methods
    methods = {
        1: "Grayscale",
        2: "Red-minus-Blue",
        3: "Color Convergence"
    }

    print("\nTesting different shoreline detection methods...")

    for method_id, method_name in methods.items():
        print(f"\n  Method {method_id}: {method_name}")

        # Update detector method
        detector.method = method_id

        # Load first image
        img = loader.load_image(0)

        if img is None:
            print("    Error loading image")
            continue

        # Detect shoreline
        shoreline = detector.detect(img)

        if shoreline is None:
            print("    Detection failed")
            continue

        # Compute statistics
        valid_positions = shoreline[~np.isnan(shoreline)]

        if len(valid_positions) > 0:
            print(f"    Mean position: {np.mean(valid_positions):.1f} pixels")
            print(f"    Std deviation: {np.std(valid_positions):.1f} pixels")
            print(f"    Valid points: {len(valid_positions)}/{len(shoreline)}")

            # Visualize result
            visualizer = CoastCamsVisualizer(config.output_dir)
            visualizer.plot_shoreline(
                img, shoreline,
                title=f'Shoreline Detection - Method {method_id} ({method_name})',
                filename=f'shoreline_method_{method_id}.png'
            )

    # Detect shorelines in all images
    print(f"\n\nDetecting shorelines in all {num_images} images...")

    detector.method = config.shoreline_method  # Reset to config method

    all_shorelines = []
    for i in range(num_images):
        img = loader.load_image(i)
        if img is not None:
            shoreline = detector.detect(img)
            all_shorelines.append(shoreline)

    # Compute statistics across all images
    if len(all_shorelines) > 0:
        mean_pos, std_pos = detector.get_shoreline_variation(all_shorelines)

        print(f"\nShoreline Statistics (all images):")
        print(f"  Mean position: {mean_pos:.1f} pixels")
        print(f"  Std deviation: {std_pos:.1f} pixels")
        print(f"  Position in meters: {mean_pos * config.pixel_resolution:.2f} m")

    print("\n" + "="*70)
    print("Shoreline detection complete!")
    print(f"Results saved to: {config.output_dir}")
    print("="*70)


if __name__ == '__main__':
    main()
