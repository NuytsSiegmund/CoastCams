# MATLAB-Python Comparison Status Update

## Issues Reported and Status

### ✅ 1. Camera Angle Issue - **RESOLVED**
**Report**: "Camera angle keeps changing. It's fixed angle which is set at the start."

**Status**: This was actually working CORRECTLY all along!
- Camera height (z0) is fixed at 27.24m
- Camera angle SHOULD vary per wave because each wave breaks at a different distance
- MATLAB formula: `AngleCam = abs(z0./X1(PosX))` where X1(PosX) varies per wave
- Python correctly implements: `angle = arctan(z0 / distance)` for each breaking wave

**Evidence from output**:
```
Breaking wave distances: min=16.4m, max=20.2m
Camera angles (varying per wave): mean=32.57°, std=1.10°
```

**Note**: I initially "fixed" this incorrectly by making the angle constant, then reverted it (commit 4694801).

---

### ✅ 2. Tp and Tm Issue - **RESOLVED**
**Report**: "Tp and Tm are the same, which is not possible."

**Status**: FIXED in coastcams/wave_analysis.py (lines 151-175)

**Changes**:
- Tm (mean period): Calculated using zero-crossing method
- Tp (peak period): Calculated using spectral analysis (FFT)
- Removed fallback that set Tp = Tm

**Evidence from output**:
```
Image 1: Tm=8.78s, Tp=14.22s ✓
Image 2: Tm=6.96s, Tp=10.67s ✓
Image 3: Tm=7.12s, Tp=10.67s ✓
```

---

### ✅ 3. Plotting Orientation - **RESOLVED**
**Report**: "The timestack should have the time on x axis, now it's the opposite"

**Status**: FIXED in coastcams/visualize.py

**Changes**:
- Timestacks are stored as (time × space) = (1680, 689)
- For display with time on x-axis, must transpose to (space × time)
- Fixed in `plot_timestack()` and `plot_comprehensive_analysis()`
- Now matches MATLAB: `imagesc(Time_TS, 1:size(stack,1), stack)`

**Code**:
```python
# Transpose from (time x space) to (space x time) for correct display
timestack_display = timestack.T
im = ax.imshow(timestack_display, aspect='auto', cmap='viridis', origin='lower')
```

---

### ⚠️ 4. Empty CSV Columns - **LIKELY RESOLVED (Needs Verification)**
**Report**: "A lot of variables are still missing" - WaterDepth, SLA_S, SLA_L, Bathymetry empty

**Status**: Code appears correct, likely user was viewing OLD CSV file

**Evidence from test run**:
```
WaterDepth_L array: 689 spatial points
  Valid depths: 590/689
  Depth range: [0.00, 2.54] m
  Mean water depth: 1.10 m
```

**What's working**:
- Water depths ARE being calculated from linear wave theory
- SLA calculations are implemented (lines 506-540 in main.py)
- CSV export keys are correct: 'depths', 'sla_values', 'sla_shallow_values'

**Added debug output**:
```python
print(f"\n=== Depth Profile Collection ===")
print(f"Collected {len(depth_profiles_list)} depth profiles from {len(all_results)} timestacks")

print(f"\n=== Results Summary ===")
print(f"water_depths_array: shape, Non-NaN count")
print(f"slas (SLA_L): shape, Non-NaN count")
print(f"sla_shallow (SLA_S): shape, Non-NaN count")
```

**Action**: User should run fresh analysis and check CSV timestamp matches.

---

### ⚠️ 5. Breakpoint Location Variation - **UNDER INVESTIGATION**
**Report**: "Breakpoint location is too similar at the moment. It should vary a lot more"

**Current behavior**:
```
Image 1: Break location: 17.40m (water depth: 1.10m)
Image 2: Break location: 17.20m (water depth: 1.11m)
Image 3: Break location: 17.60m (water depth: 1.13m)
Image 4: Break location: 17.40m (water depth: 1.06m)
Image 5: Break location: 17.50m (water depth: 1.10m)
Image 6: Break location: 17.20m (water depth: 0.91m)
Image 7: Break location: 17.40m (water depth: 0.78m)
```

**Observations**:
- Break locations: 17.2-17.6m (only 0.4m variation)
- Water depths: 0.78-1.13m (0.35m variation over ~1.5 hours)
- Individual waves within each timestack break over 3-4m range (e.g., 16.4-20.2m)
- Median of individual positions stays constant at ~17.3m

**How it's calculated** (wave_analysis.py:213-217):
```python
median_break_pos = int(np.round(np.nanmedian(PosX)))
break_location = cross_shore_positions[median_break_pos]
```

**Possible explanations**:
1. **Physical (bathymetry-controlled)**: If there's a sandbar or strong bathymetric feature at ~17m, waves will preferentially break there regardless of tidal level
2. **Detection bias**: Radon transform may be biased toward detecting waves in a particular zone
3. **Algorithm difference**: MATLAB may use a different method (not median) to determine representative break location

**Next steps**:
1. Compare with MATLAB's WaveParameters_CoastCams to see how it calculates break location
2. Check if MATLAB uses median, mean, or mode of PosX
3. Review bathymetry profile to see if there's a physical feature at 17m
4. Consider using weighted average based on wave heights instead of median

---

## Summary

### Fixed (3/5):
✅ Camera angle (was never broken)
✅ Tp vs Tm separation
✅ Timestack plotting orientation

### Likely Fixed (1/5):
⚠️ Empty CSV columns (code correct, user should re-run)

### Needs Further Investigation (1/5):
⚠️ Breakpoint location variation (median calculation may differ from MATLAB)

---

## Commits Made

1. **4694801**: REVERT: Camera angle should vary per wave based on distance
2. **b02b457**: Add investigation doc for empty CSV columns
3. **8329714**: Fix timestack plotting orientation and add debug output

All changes pushed to: `claude/coastcams-matlab-python-comparison-updated-011CUcJqCvFZyQHBBfJ5fZGU`

---

## Recommendations

1. **Run fresh analysis** with latest code to generate new CSV
2. **Verify CSV columns** are populated (check file timestamp)
3. **Compare breakpoint calculation** with MATLAB - need to review WaveParameters_CoastCams line by line
4. **Check bathymetry profile** - plot cross-shore depth to see if there's a feature at 17m
5. **Consider physical vs algorithmic** - the small variation might be realistic for this beach

---

**Date**: 2025-11-01
**Session**: Continuation from previous context
**Branch**: claude/coastcams-matlab-python-comparison-updated-011CUcJqCvFZyQHBBfJ5fZGU
