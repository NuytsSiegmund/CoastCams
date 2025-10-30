# MATLAB vs Python CoastCams - Critical Bug Fixes

## Summary

Fixed critical bugs in Python's cross-correlation implementation that were causing:
- Wave celerities of 0.1 m/s instead of realistic 4-8 m/s
- NaN depth values for 13 out of 16 timestacks
- NaN SLA values throughout the analysis

## Root Causes Identified

### 1. Wrong Cross-Correlation Algorithm
**MATLAB Approach** (`CrossCorrelation_CoastCams.m` lines 28-81):
- Uses **FIXED time lag**: `n = floor(dpha/dt) = floor(1/0.5) = 2 frames` (1 second)
- **Searches for spatial lag**: Computes correlation at different spatial separations (1 to dc-1 pixels)
- Finds the spatial lag that maximizes correlation
- Converts to velocity: `spatial_lag / time_lag`

**Python Original Approach** (WRONG):
- Used **FIXED spatial separation** (100 pixels)
- **Searched for time lag** that maximizes correlation
- Could find spurious correlations at incorrect time lags
- Failed for weak or noisy signals

**Why This Matters**:
- MATLAB's approach constrains the search to a physically meaningful time scale (wave period)
- Python's approach could find correlations at any time lag, including noise correlations
- For the test data, Python was consistently finding maximum correlation at spatial lag=1 pixel, suggesting the transpose bug was making all positions look identical

### 2. Incorrect Timestack Transpose
**Issue**: Python code was transposing the timestack from (time, space) to (space, time), effectively swapping dimensions

**Evidence from main.py**:
```python
# Line 206: "Input timestack image (time x space x channels), e.g., (1680, 689, 3)"
# Line 211: "Preprocessed timestack (time x space), e.g., (1680, 689)"
# Line 231: "timestack (time x space), e.g., (1680, 689)"
```

**cross_correlation.py Original Code** (WRONG):
```python
A2 = timestack.T  # Transposed (1680, 689) to (689, 1680) = WRONG!
```

**Result**:
- With swapped dimensions, the correlation was comparing time series that were actually spatial profiles
- This caused all spatial lags to have maximum correlation at lag=1
- Celerity = 1 pixel / 1 second = 1 pixel/s × 0.1 m/pixel = 0.1 m/s

### 3. Unit Conversion
**MATLAB** (`S01_AnalysisTimestackImages.m` line 242):
```matlab
WaveCelerity(i,:) = movmean(Cf1./10, 10)
```
- `Cf1` is in pixels/second
- Division by 10 converts to m/s (since pixel resolution = 0.1 m/pixel)
- Applies moving average with window=10

**Python Fixed**:
```python
velocities_pixels_per_sec = spatial_lags / dpha  # pixels/second
velocities_ms = velocities_pixels_per_sec * self.dx  # m/s (dx = 0.1)
velocities_smooth = uniform_filter1d(velocities_ms, size=10)  # moving average
```

## Files Modified

### `CoastCams_Python/coastcams/cross_correlation.py`

1. **Rewrote `_extract_wave_properties()` method**:
   - Calls new `_compute_celerity_matlab_style()` method
   - Returns full array of celerities (with NaN at edges where correlation window doesn't fit)

2. **Added `_compute_celerity_matlab_style()` method** (lines 177-288):
   - Implements MATLAB's exact algorithm
   - Fixed time lag: `n = floor(dpha/dt)` frames
   - Variable spatial lag: loops from 1 to dc-1 pixels
   - Computes correlation coefficient at each spatial lag
   - Finds maximum correlation -> spatial lag
   - Converts to velocity with proper units
   - Applies moving average (window=10)

3. **Removed incorrect transpose**:
   - Changed from `A2 = timestack.T` to `A2 = timestack`
   - Added documentation clarifying shape convention

### `CoastCams_Python/main.py`

1. **Updated `_estimate_bathymetry()` method** (lines 500-562):
   - Now handles full celerity array from MATLAB-style correlation
   - Interpolates NaN values at edges if enough valid points exist
   - Falls back to mean value if only few valid celerities
   - Added detailed logging of celerity validity

## MATLAB Algorithm Details

### CrossCorrelation_CoastCams.m Key Lines:

```matlab
Line 28:  n = floor(dpha/dt);  % Time lag in frames (e.g., 1s / 0.5s = 2 frames)
Line 42:  if ic - dc/2 > 0 && ic - dc/2 + lc - 1 <= nc
Line 43:    R7 = corrcoef(A2(n+1:nt, ic-dc/2), A2(1:(nt-n), ic-dc/2+lc-1));
Line 44:    R2(lc) = R7(2,1);
```

This correlates:
- Signal at position `ic-dc/2` at later times (`n+1:nt`)
- Signal at position `ic-dc/2+lc-1` at earlier times (`1:nt-n`)

Looking for onshore wave propagation: wave appears first at offshore position, then closer position.

```matlab
Line 80:  [~, R2M] = max(R2M');  % Find spatial lag with max correlation
Line 81:  R2M = R2M ./ dpha;     % Convert to pixels/second
```

### S01_AnalysisTimestackImages.m Key Lines:

```matlab
Line 234:  [Cf1, WLe1, Tp1, Hs1, RM] = CrossCorrelation_CoastCams(S2, dpha, dt, dc);
Line 242:  WaveCelerity(i,1:size(Cf1,2)) = movmean(Cf1./10, 10);
Line 247:  [df] = LinearC(Tp, WaveCelerity(i,:), 0.01);
Line 252:  WaterDepth_L = [WaterDepth_L; movmean(df(:,1:numel(sLimit)), 10)];
```

## Testing Results (Before Fix)

From debug output with transpose bug:
```
Sample correlation profile (middle position):
  First 10 lags: [0.99279904 0.9736769  0.79937607 0.53350991 ...]
  Spatial lags: mean=1.0, range=[1, 1]  <-- ALL lags are 1!
  Mean celerity: 0.100 m/s  <-- Way too low!
  Time lag n=2 frames (1.0s), dpha=1.0s
```

Maximum correlation at lag=1 indicates the transpose bug made all spatial positions look identical.

## Expected Results (After Fix)

With correct dimensions and algorithm:
- Spatial lags should vary from ~20-80 pixels (for waves moving 2-8 m/s)
- Celerities should be 2-8 m/s (typical coastal waves)
- Depth calculations should succeed for most/all timestacks
- SLA should have no NaN values (matching MATLAB's R²=0.98)

## Next Steps

1. ✅ Test corrected cross-correlation implementation
2. Verify realistic celerity values (4-8 m/s range)
3. Confirm depth calculation produces valid values for most timestacks
4. Verify SLA calculation has no NaN values
5. Compare Python SLA output to MATLAB for same dataset

## References

- MATLAB code: `/home/user/CoastCams/UserScripts/S01_AnalysisTimestackImages.m`
- MATLAB correlation: `/home/user/CoastCams/UserScripts/Subfunctions/CrossCorrelation_CoastCams.m`
- MATLAB depth solver: `/home/user/CoastCams/UserScripts/Subfunctions/LinearC.m`
