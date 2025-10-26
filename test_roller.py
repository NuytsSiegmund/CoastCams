#!/usr/bin/env python
"""Quick test of roller detection and wave height calculation."""

import sys
sys.path.insert(0, '/home/user/CoastCams/CoastCams_Python')

import numpy as np
from PIL import Image
from coastcams.roller_detection import roller_properties_taller
from coastcams.wave_analysis import WaveAnalyzer
from coastcams.config import CoastCamsConfig

# Load test image
img = Image.open('/home/user/CoastCams/Timestacks/S_1_202110130745.jpeg')
img = np.rot90(np.array(img), k=3)
timestack = img.astype(np.float32) / 255.0

print(f'Timestack shape: {timestack.shape}')

# Remove mean from each column (simple preprocessing)
timestack_preprocessed = timestack[:, :, 0].copy()
for i in range(timestack_preprocessed.shape[1]):
    timestack_preprocessed[:, i] -= np.mean(timestack_preprocessed[:, i])

print(f'Preprocessed shape: {timestack_preprocessed.shape}')

# Run roller detection
print('\nRunning roller detection...')
PosX, PosT, Lw, B, Breakstd, Breakmean1, Breakmean2 = roller_properties_taller(timestack_preprocessed, 0.5)
print(f'  Found {len(PosX)} breaking events')
print(f'  PosX range: {np.min(PosX)}-{np.max(PosX)}')
print(f'  PosT range: {np.min(PosT)}-{np.max(PosT)}')

# Create cross-shore positions
cross_shore_positions = np.arange(timestack_preprocessed.shape[1]) * 0.1

# Test BreakerHeight
print('\nTesting BreakerHeight...')
config = CoastCamsConfig()
analyzer = WaveAnalyzer(config)

try:
    hs, hm = analyzer._breaker_height(
        B, PosT, PosX, Lw,
        cross_shore_positions, 0.1
    )
    print(f'  Result: hs={hs:.3f}m, hm={hm:.3f}m')
except Exception as e:
    print(f'  Error: {e}')
    import traceback
    traceback.print_exc()
