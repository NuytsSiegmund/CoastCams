# Diagnosis: Why Columns Are Still Empty

## The Problem

After the complete rewrite, the CSV still has empty columns:
- WaterDepth: **EMPTY**
- SLA_L: **EMPTY**
- SLA_S: **EMPTY**
- BreakpointDepth: **EMPTY**
- Rows 13-14 have absurd wave heights (191m, 199m)

## Root Cause Analysis

### Issue 1: Cross-Correlation Not Producing Valid Arrays

The workflow in MATLAB is:
```matlab
% Line 234: Get Cf1 and WLe1 from cross-correlation
[Cf1, WLe1, Tp1, Hs1, RM] = CrossCorrelation_CoastCams(S2, dpha, dt, dc);

% Line 242-243: Store in matrices with smoothing
WaveCelerity(i,1:size(Cf1,2)) = movmean(Cf1./10, 10);  % NOTE: /10 division!
WaveLength(i,1:size(WLe1,2)) = movmean(WLe1,10);
```

In our Python code, we're calling:
```python
corr_results = correlation_analyzer.analyze_timestack(preprocessed)
Cf1 = corr_results.get('Cf1', np.array([]))
WLe1 = corr_results.get('WLe1', np.array([]))
```

**Problem**: The `analyze_timestack()` method expects certain input format and may not be producing the Cf1/WLe1 arrays that match MATLAB's output.

### Issue 2: Timestack Shape Mismatch

MATLAB CrossCorrelation_CoastCams expects:
- Input: `S2` which is preprocessed timestack
- Format: The code has `A2(time, space)` according to comments

Our Python code passes:
- `preprocessed` which is `(time, space)` from `ImagePreprocessor`

But looking at the cross_correlation.py code:
```python
# Line 222-224 in cross_correlation.py:
# Input timestack is already (time x space) format from main.py
# Example: (1680, 689) = (time points, spatial positions)
A2 = timestack
```

This **should** be correct, but we need to verify the correlation is actually computing.

### Issue 3: Missing Division by 10

MATLAB does: `Cf1./10` before storing in WaveCelerity.

Our code doesn't divide by 10:
```python
Cf1_smoothed = movmean(Cf1, 10)  # Missing /10!
```

This could cause unit issues downstream.

### Issue 4: Linear Depth Calculation

MATLAB (line 247):
```matlab
[df] = LinearC(Tp, WaveCelerity(i,:), 0.01);
```

Our Python code:
```python
depths_linear = bathymetry_estimator.estimate_depth_linear_wave_theory(
    peak_period=Tp if not np.isnan(Tp) else 10.0,
    celerities=celerity_array
)
```

If `celerity_array` is all NaN from the correlation step, then `depths_linear` will be all NaN.

## Immediate Fixes Needed

### Fix 1: Add /10 Division to WaveCelerity

In main.py around line 361:
```python
# OLD:
Cf1_smoothed = movmean(Cf1, 10)

# NEW:
Cf1_smoothed = movmean(Cf1 / 10.0, 10)  # Match MATLAB's /10
```

### Fix 2: Add Debug Output to See What Correlation Returns

Add after line 351:
```python
corr_results = correlation_analyzer.analyze_timestack(preprocessed)
print(f"  DEBUG: Cf1 shape: {corr_results.get('Cf1', np.array([])).shape}")
print(f"  DEBUG: Cf1 has {np.sum(~np.isnan(corr_results.get('Cf1', np.array([]))))} non-NaN values")
```

### Fix 3: Check If Correlation Analyzer Needs Different Input

The CrossCorrelationAnalyzer might not be implementing the MATLAB algorithm correctly. We may need to:

1. Call the MATLAB-style method directly
2. Or re-implement CrossCorrelation_CoastCams from scratch
3. Or use a simpler correlation approach that we know works

### Fix 4: Verify Matrix Building Logic

After converting lists to matrices (line 427-431), add debug:
```python
WaveCelerity = pad_to_matrix(WaveCelerity)
print(f"DEBUG: WaveCelerity matrix built: {WaveCelerity.shape}")
print(f"DEBUG: Non-NaN values in WaveCelerity: {np.sum(~np.isnan(WaveCelerity))}/{WaveCelerity.size}")
print(f"DEBUG: Sample values: {WaveCelerity[0, :10]}")  # First row, first 10 values
```

## Next Steps

1. **Run with debug output** to see what correlation actually returns
2. **Add the /10 division** to match MATLAB
3. **Check if WaveCelerity matrix has valid values** before smooth2
4. **If still NaN**, investigate the CrossCorrelationAnalyzer implementation

## Alternative Approach

If the CrossCorrelationAnalyzer is too complex and not working, we could:

1. Implement a simpler MATLAB-matching correlation function directly in main.py
2. Or import and call the exact MATLAB algorithm parameters
3. Or use a known-working correlation method from the existing codebase

The key is to get valid Cf1 arrays so that:
- WaveCelerity matrix is populated
- WaterDepth_L can be calculated
- SLA can be calculated
- All columns in CSV are filled

---

**Created**: 2025-11-01
**Status**: Active Investigation
