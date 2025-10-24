# CoastCams Python - Issues and Fixes

## Date: 2025-10-24

## Critical Issues Identified

### 1. Cross-Shore Position Dimension Bug (FIXED)
**Location**: `main.py` lines 224-225 and 340-341

**Problem**: Using `timestack.shape[0]` (time dimension) instead of `timestack.shape[1]` (space dimension) to create cross-shore positions.

```python
# WRONG:
num_positions = timestack.shape[0]  # This is TIME (1680 frames)

# CORRECT:
num_positions = timestack.shape[1]  # This is SPACE (689 pixels)
```

**Impact**: Caused all downstream calculations to use wrong spatial coordinates.

**Status**: ✅ FIXED in `main.py`

---

### 2. Photogrammetric Wave Height Method is Fundamentally Broken
**Location**: `coastcams/wave_analysis.py` lines 334-481

**Problem**: The MATLAB `BreakerHeight` function requires outputs from `RollerPropertiesTaller`:
- `PosT`: Temporal positions of individual breaking wave events
- `PosX`: Spatial positions of breaking waves
- `Lw`: Roller lengths

The Python version attempts to reconstruct this information from the timestack alone, but:
1. The calculation produces negative wave heights (correction > L_horizontal)
2. Taking absolute value gives heights 10x too large (50m vs MATLAB's 5.5m)
3. The method requires individual wave event detection that isn't implemented

**Root Cause**: Missing prerequisite - need to implement `RollerPropertiesTaller` equivalent first.

**Current Workaround**: Method returns NaN and falls back to time-series method.

**Status**: ⚠️ REQUIRES MAJOR REFACTORING - Beyond scope of simple fix

---

### 3. Wave Analysis Time vs Space Confusion
**Location**: `coastcams/wave_analysis.py` lines 394-410

**Problem**: Original code computed max-min along SPATIAL dimension for each TIME step, but should compute along TIME dimension for each SPATIAL position.

**Status**: ✅ FIXED - Now correctly processes spatial profiles

---

### 4. Angle vs Tangent Confusion
**Location**: `coastcams/wave_analysis.py` line 371

**Problem**: MATLAB's `AngleCam` variable is actually `tan(angle)`, not the angle itself. Python code was using `arctan()` then `tan()` again, applying the transformation twice.

**Status**: ✅ FIXED - Now uses `camera_angle_tan` directly

---

## Comparison: MATLAB vs Python Results

### MATLAB Output:
- Mean Significant Wave Height: **5.48 m**
- Mean Peak Wave Period: **12.39 s**

### Python Output (Current):
- Mean Significant Wave Height: **50.74 m** (with broken photogrammetric)
- Mean Wave Period: **10.69 s**

### Python Output (Fallback method only):
- Mean Significant Wave Height: **~1.08 m** (5x too small)
- Mean Wave Period: **10.69 s** (reasonably close)

---

## Required Next Steps

### Immediate (Simple Fixes):
1. ✅ Fix cross-shore dimension bug
2. ✅ Fix time/space confusion in photogrammetric method
3. ✅ Fix angle/tangent confusion
4. ⚠️ Disable broken photogrammetric method (return NaN)
5. ⚠️ Remove debug output

### Medium Term (Requires Implementation):
1. Implement `RollerPropertiesTaller` equivalent:
   - Detect individual wave breaking events
   - Track roller properties
   - Extract temporal and spatial positions
2. Fix photogrammetric method with proper inputs
3. Validate against MATLAB results

### Alternative Approach:
Instead of implementing photogrammetric method, calibrate the time-series method to match MATLAB using physical principles (not arbitrary scaling factor).

---

## Files Modified:
- `CoastCams_Python/main.py` - Fixed cross-shore dimension bug
- `CoastCams_Python/coastcams/wave_analysis.py` - Multiple fixes to photogrammetric method

---

## Testing:
- Processed 16 timestack images
- All images load correctly
- Cross-shore positions now correct (689 spatial points)
- Photogrammetric method processes 163 iterations per image
- Wave period calculation appears correct (~10-12s vs MATLAB ~12s)
- Wave height still requires proper implementation

---

## Recommendation:
The photogrammetric wave height method requires significant additional work (RollerPropertiesTaller implementation). For a working solution, consider:

**Option A**: Implement roller detection as prerequisite
**Option B**: Use alternative wave height estimation method based on validated physical principles
**Option C**: Focus on other analysis outputs that are working correctly (periods, celer, bathymetry)

