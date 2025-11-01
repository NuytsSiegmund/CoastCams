# Critical Bug Fix: movmean NaN Handling

## Problem Identified

**Root Cause:** The `movmean()` function was using scipy's `uniform_filter1d`, which propagates NaN values throughout the smoothing window. This caused ALL smoothed values to become NaN when the input contained any NaN values.

## Evidence from Debug Output

```
DEBUG: Cf1 shape: (689,)
DEBUG: Cf1 non-NaN values: 590/689          ✅ Correlation produces valid data
DEBUG: Cf1 range: [0.10, 5.12]              ✅ Reasonable values
DEBUG: Stored celerity array with 0 non-NaN values  ❌ ALL DATA LOST after smoothing!
```

**Analysis:**
- Cross-correlation produces 590 valid values out of 689 (14% NaN)
- After applying movmean with window=10, result has 0 valid values (100% NaN)
- With 99 NaN values distributed across 689 positions, a window size of 10 means nearly every window contains at least one NaN
- scipy's `uniform_filter1d` propagates NaN: if ANY value in the window is NaN, the output is NaN

## Impact

This bug caused a cascade of empty values:

```
Cf1 (correlation)
    ↓ movmean → ALL NaN
WaveCelerity matrix = ALL NaN
    ↓ LinearC(Tp, WaveCelerity)
WaterDepth_L = ALL NaN (can't calculate from NaN celerity)
    ↓ smooth2, then -mean
SLA_L = ALL NaN (mean of NaN = NaN)
    ↓ WaterDepth - SLA_S
Bathymetry = ALL NaN
```

**Result:** Empty CSV columns for:
- WaterDepth
- SLA_L
- SLA_S
- Bathymetry
- BreakpointDepth

## Solution

Replaced scipy's `uniform_filter1d` with pandas' `rolling().mean()`:

```python
def movmean(data: np.ndarray, window: int) -> np.ndarray:
    """Moving average matching MATLAB's movmean with NaN handling."""
    import pandas as pd

    series = pd.Series(data)
    result = series.rolling(window=window, center=True, min_periods=1).mean()
    return result.values
```

**Key improvements:**
1. **NaN-aware:** Computes mean of non-NaN values in each window
2. **center=True:** Matches MATLAB's centered window behavior
3. **min_periods=1:** Produces output even if only 1 valid value in window

## Expected Results

With the fix:
- Input: 590/689 non-NaN values (14% NaN)
- Output: ~590/689 non-NaN values (preserves valid data)
- WaveCelerity matrix: Populated with valid celerity values
- WaterDepth_L: Calculated from valid celerities
- SLA_L, SLA_S: Calculated from valid depth matrices
- Bathymetry: Calculated from valid depth and SLA values

## Additional Fix: Wave Height Sanity Check

Added bounds checking for absurd wave heights (e.g., 191m, 199m):

```python
if not np.isnan(hs) and hs > 20.0:
    print(f"WARNING: Absurd wave height {hs:.1f}m detected, setting to NaN")
    hs = np.nan
```

This prevents WaveAnalyzer errors from corrupting the CSV output.

## Files Modified

- `CoastCams_Python/main.py`:
  - Commit eceafdb: Fixed movmean NaN handling
  - Commit 45b325b: Added wave height sanity check

## Testing

Run the analysis and verify:

1. **Debug output shows preservation of valid values:**
   ```
   DEBUG: Cf1 non-NaN values: 590/689
   DEBUG: After /10 division: 590/689 non-NaN values
   DEBUG: After smoothing: ~590/689 non-NaN values  ← Should preserve data!
   ```

2. **CSV contains valid values in previously empty columns:**
   - WaterDepth: Non-empty depth values
   - SLA_L: Sea level anomaly values
   - SLA_S: Shallow water SLA values
   - Bathymetry: Calculated bathymetry

3. **No absurd wave heights:**
   - All Hs values < 20m (typical: 0.2-2.0m)
   - Warnings logged for any filtered values

## MATLAB Comparison

MATLAB's `movmean()` has implicit NaN-aware behavior in newer versions, or the MATLAB implementation may produce fewer NaN values from cross-correlation. The key is that MATLAB's workflow successfully produces valid smoothed celerity values, which our scipy implementation was destroying.
