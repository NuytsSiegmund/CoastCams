"""
Quick test of Radon transform on single image.
"""
import numpy as np
import cv2
from coastcams.matlab_preprocessing import RollerDetector, PhotogrammetricHeightCalculator
from coastcams.wave_analysis import WaveAnalyzer
from coastcams.config import CoastCamsConfig

# Create default config
config = CoastCamsConfig()
config.pixel_resolution = 0.1
config.camera_height = 27.24
config.acquisition_frequency = 2.0

# Load a single image
img_path = '../Timestacks/S_1_202110130745.jpeg'
img = cv2.imread(img_path)

print(f"Loaded image (raw): {img.shape}")

# Rotate 270 degrees (matching config.yaml rotation_angle)
img = cv2.rotate(img, cv2.ROTATE_90_COUNTERCLOCKWISE)

print(f"Loaded image (rotated): {img.shape}")

# Convert to grayscale and normalize
if len(img.shape) == 3:
    timestack = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY).astype(np.float32) / 255.0
else:
    timestack = img.astype(np.float32) / 255.0

print(f"Timestack shape: {timestack.shape}")

# Create analyzer
analyzer = WaveAnalyzer(config)

# Create cross-shore positions
cross_shore_positions = np.arange(timestack.shape[1]) * config.pixel_resolution

print(f"\nCross-shore positions:")
print(f"  Range: {cross_shore_positions.min():.1f} - {cross_shore_positions.max():.1f} m")
print(f"  Shape: {cross_shore_positions.shape}")

print("\n" + "="*70)
print("Testing FULL PHOTOGRAMMETRIC PIPELINE with Radon transform")
print("="*70)

# Detect rollers
print("\n1. Detecting breaking waves...")
roller_results = analyzer.roller_detector.detect_rollers(timestack)

PosX = roller_results['PosX']
PosT = roller_results['PosT']
Lw = roller_results['Lw']
B = roller_results['B']

print(f"   Detected {len(PosX)} breaking waves")
print(f"   PosX range: {PosX.min()}-{PosX.max() if len(PosX) > 0 else 'N/A'}")
print(f"   PosT range: {PosT.min()}-{PosT.max() if len(PosT) > 0 else 'N/A'}")

if len(PosX) > 0:
    # Calculate wave heights
    print("\n2. Calculating photogrammetric wave heights...")
    Hs, Hm = analyzer.height_calculator.calculate_wave_heights(
        B, PosT, PosX, Lw, cross_shore_positions
    )

    print(f"\nRESULTS:")
    print(f"  Significant wave height (Hs): {Hs:.3f} m")
    print(f"  Mean wave height (Hm): {Hm:.3f} m")

    if not np.isnan(Hs) and Hs > 0:
        print("\n✓ SUCCESS: Photogrammetric method worked!")
    else:
        print("\n✗ FAILED: Photogrammetric method returned NaN")
else:
    print("\n✗ No breaking waves detected")
