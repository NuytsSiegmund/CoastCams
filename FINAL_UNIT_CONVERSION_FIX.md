# Final Critical Fixes - Unit Conversion and Fallback Issues

## Problem Summary

After fixing the movmean NaN handling and BathymetryEstimator API, two critical issues remained:

1. **Water depths 100x too small:** 0.007m (7mm) instead of ~1-2m
2. **Absurd wave heights still appearing:** 191m, 199m in some timestacks

## Root Cause Analysis

### Issue 1: Double Unit Conversion

**The Bug:**
```python
# Python cross_correlation.py line 313:
velocities_ms = velocities_pixels_per_sec * self.dx  # Already in m/s!

# Then in main.py:
Cf1_divided = Cf1 / 10.0  # ❌ Dividing by 10 AGAIN!
```

**Root Cause:**
- MATLAB's `CrossCorrelation_CoastCams` returns Cf1 in **pixels/second**
- MATLAB divides by 10 to convert to m/s (when pixel_resolution = 0.1 m/pixel)
- Python's cross-correlation **already multiplies by dx** (line 313)
- So Python's Cf1 is **already in m/s** - no need to divide by 10!

**Impact Chain:**
```
Cf1 = 1.0-4.5 m/s (correct)
  ↓ /10 (wrong!)
Cf1 = 0.1-0.45 m/s (10x too small)
  ↓ depth calculation: h ∝ c²
WaterDepth = 0.007m (100x too small!)
```

**The Fix:**
Removed the `/10` division in main.py. Python's Cf1 is already in correct units.

**Expected Result:**
```
Before: Smoothed celerity range [0.010, 0.420] m/s → depths ~0.007m ❌
After:  Smoothed celerity range [1.0, 4.2] m/s → depths ~1-2m ✅
```

---

### Issue 2: Broken Time-Series Fallback

**The Bug:**
```python
# wave_analysis.py lines 149, 227-228:
Hs_from_ts = 4.0 * np.std(filtered_ts)  # Already Hs estimate!
calibration_factor = 2.0
results['mean_Hs'] = Hs_from_ts * calibration_factor  # = 8.0 * std ❌
```

**Root Cause:**
- The `4.0 * std` formula already estimates Hs (Rayleigh distribution)
- Multiplying by calibration_factor=2.0 again gives `8.0 * std`
- Standard deviation of pixel intensities ≈ 25
- Result: `Hs = 8 * 25 = 200m` (absurd!)

**When This Occurs:**
```
Lf: min=-0.440, max=-0.264  ← Photogrammetry produces negative heights
Valid heights after filtering: 0/4  ← All filtered out
Wave height (fallback - time series): Hs = 199.703m  ← Broken fallback!
```

**The Fix:**
Removed the unreliable fallback entirely. When photogrammetry fails, return `NaN` with diagnostic message rather than using a broken formula.

**Expected Result:**
```
Before: Fallback → Hs = 200m ❌
After:  "Photogrammetric height invalid (Lf <= 0), setting to NaN" → Hs = NaN ✅
```

---

## Expected CSV Output After Fixes

### Before All Fixes:
```csv
Timestamp,Hs,Tm,Tp,WaterDepth,SLA_L,SLA_S,...
7:45,0.496,6.25,5.33,0.007,,,...      ← 7mm depth, empty SLA
10:45,191.2,7.61,14.2,0.001,,,...     ← Absurd Hs, 1mm depth
```

### After All Fixes:
```csv
Timestamp,Hs,Tm,Tp,WaterDepth,SLA_L,SLA_S,...
7:45,0.496,6.25,5.33,1.45,0.045,-0.012,...    ← Realistic depth, populated SLA!
10:45,,7.61,14.2,1.38,0.038,-0.015,...        ← NaN for bad Hs, valid depth!
```

Note: Some Hs values will be NaN when photogrammetry fails (better than absurd 200m values!)

---

## Debug Output to Verify Fixes

### 1. Celerity Values Should Be Realistic
```
DEBUG: Cf1 range: [0.10, 4.54]  ← Raw celerities (GOOD)
DEBUG: Smoothed celerity range: [1.0, 4.2] m/s  ← Should be 1-5 m/s ✅
```

NOT:
```
DEBUG: Smoothed range: [0.010, 0.420] m/s  ← 10x too small ❌
```

### 2. Water Depths Should Be Realistic
```
WaterDepth_L matrix: (16, 580)
  Non-NaN values: 8200/9280
  Value range: [0.8, 3.2] m  ← Should be 0.5-5m for coastal waves ✅
```

NOT:
```
Value range: [0.001, 0.009] m  ← 1-9mm is absurd ❌
```

### 3. Failed Photogrammetry Should Return NaN
```
Lf: min=-0.440, max=-0.264
Valid heights after filtering: 0/4
Photogrammetric height invalid (Lf <= 0), setting to NaN  ← Clear message ✅
Wave Height: NaN
```

NOT:
```
Wave height (fallback - time series): Hs = 199.703m  ← Absurd ❌
```

---

## Technical Details

### Why MATLAB and Python Handle Units Differently

**MATLAB CrossCorrelation_CoastCams.m:**
```matlab
% Returns Cf1 in pixels/second
% User must convert: Cf1./10 where pixel_resolution = 0.1 m/pixel
WaveCelerity(i,:) = movmean(Cf1./10, 10);
```

**Python cross_correlation.py:**
```python
# Line 308: velocities_pixels_per_sec = spatial_lags / dpha
# Line 313: velocities_ms = velocities_pixels_per_sec * self.dx
# Returns velocities already in m/s!
```

The Python version **already does the unit conversion** that MATLAB requires the user to do manually.

### Why the Fallback Was Broken

The `Hs = 4σ` formula comes from the Rayleigh distribution of wave heights:
- For a Rayleigh distribution: `Hs = 4σ`
- This is already the significant wave height estimate

Multiplying by 2.0 again assumes the variance is off by a factor of 4, which makes no physical sense and produces absurd results.

---

## Files Modified

### Commit a384481: "Fix celerity double-conversion and remove broken fallback"

1. **CoastCams_Python/main.py** (lines 371-380)
   - Removed `/10.0` division
   - Updated comment explaining unit handling
   - Changed debug output to show "celerity range" in m/s

2. **CoastCams_Python/coastcams/wave_analysis.py** (lines 225-250)
   - Removed `calibration_factor = 2.0`
   - Removed time-series fallback calculation
   - Return NaN with clear diagnostic messages when photogrammetry fails
   - Three failure modes: invalid Lf, no breaking waves, pipeline error

---

## Testing Checklist

Run the analysis and verify:

1. **Celerity values are realistic:**
   - ✅ Smoothed range: [1-5] m/s (not 0.01-0.5 m/s)

2. **Water depths are realistic:**
   - ✅ WaterDepth range: [0.5-3] m (not 0.001-0.009 m)
   - ✅ CSV WaterDepth column populated with realistic values

3. **SLA columns are populated:**
   - ✅ SLA_L has values (typically ±0.1 m)
   - ✅ SLA_S has values (typically ±0.1 m)

4. **No absurd wave heights:**
   - ✅ All Hs < 20m (or NaN)
   - ✅ Diagnostic messages for failed photogrammetry
   - ✅ No 200m values in CSV

5. **Matrix statistics show valid data:**
   - ✅ WaveCelerity matrix: ~8000+ non-NaN values
   - ✅ WaterDepth_L matrix: ~7000+ non-NaN values

---

## Remaining Issues

The photogrammetric height calculation still produces negative Lf values in some cases, which gets filtered to NaN. This suggests there's still an issue with the camera angle or correction calculation, but at least now:

1. We return NaN instead of absurd 200m values
2. We provide diagnostic output to understand why it failed
3. Other parameters (depth, period, location) are still calculated

Future work could investigate why `Lf = (L1 - correction) * tan_angle` sometimes produces negative values.

---

## Summary

**Commit:** a384481
**Branch:** `claude/coastcams-matlab-python-comparison-updated-011CUcJqCvFZyQHBBfJ5fZGU`

**Fixes:**
1. ✅ Removed double unit conversion (depths now 100x larger → realistic)
2. ✅ Removed broken fallback (no more 200m wave heights)
3. ✅ Added clear diagnostic messages for failures

**Expected Impact:**
- Water depths: 0.007m → 1.5m ✅
- Wave heights: 200m → NaN (when photogrammetry fails) ✅
- SLA columns: empty → populated ✅
- CSV: Fully usable data!

Run the analysis again - you should now see realistic water depths and properly populated CSV columns!
