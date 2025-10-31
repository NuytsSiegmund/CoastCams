# Python Restructure Plan to Match MATLAB

## Current Problem

Python doesn't follow MATLAB's exact workflow sequence and misses critical outputs.

## MATLAB Workflow (Lines 210-280)

```matlab
% For each timestack image:

1. WaveParameters_CoastCams() returns:
   - C (celerity from Radon)
   - depth_s (shallow water depth = C²/g)
   - depth_l (linear wave depth from Radon celerity)
   - hs, hm, Tm, Tp (wave parameters)
   - Er, rolL (energy, roller length)
   - BreakLocs, BreakDepth (breaking parameters)

2. Store: WaterDepth(i,:) = movmean(depth_s./10, 10)

3. CrossCorrelation_CoastCams(S2, ...) returns:
   - Cf1 (celerity from cross-correlation)
   - WLe1 (wavelength)
   - Tp1, Hs1 (wave parameters)

4. Store: WaveCelerity(i,:) = movmean(Cf1./10, 10)

5. Store: WaveLength(i,:) = movmean(WLe1, 10)

6. LinearC(Tp, WaveCelerity(i,:), 0.01) returns:
   - df (depth from CrossCorr celerity)

7. Store: WaterDepth_L = [WaterDepth_L; movmean(df, 10)]

8. Calculate SLA_S from WaveCelerity matrix
9. Calculate SLA_L from WaterDepth_L matrix
10. Calculate RTR
```

## Key Insight: TWO DIFFERENT CELERITIES

**MATLAB has TWO separate celerity calculations:**

1. **From WaveParameters_CoastCams** (Radon-based):
   - `C` = celerity from RadonCIndiv
   - Used to calculate `depth_s` and `depth_l`
   - Stored in `WaterDepth` variable

2. **From CrossCorrelation_CoastCams**:
   - `Cf1` = celerity from cross-correlation
   - Used to calculate `WaveCelerity`
   - Used in LinearC to calculate `WaterDepth_L`
   - Used to calculate `SLA_S`

## Python Current Issues

1. ❌ Only uses CrossCorrelation celerity
2. ❌ Doesn't capture WaveParameters_CoastCams celerity (C)
3. ❌ Doesn't calculate depth_s from WaveParameters
4. ❌ Missing proper wavelength (WLe1) from CrossCorrelation
5. ❌ SLA_S calculation is empty because celerity matrix has NaN

## Required Changes

### 1. Modify wave_analysis.py

Add method to return ALL WaveParameters_CoastCams outputs:
```python
def analyze_timestack(self, timestack, cross_shore_positions):
    """Match MATLAB WaveParameters_CoastCams exactly."""
    return {
        'C': radon_celerities,           # NEW
        'depth_s': shallow_depths,        # NEW
        'depth_l': linear_depths,         # NEW
        'hs': Hs,
        'hm': Hm,
        'Tm': Tm,
        'Tp': Tp,
        'Er': energy_array,               # NEW (full array)
        'rolL': roller_length_array,      # NEW (full array)
        'nbwave': n_waves,                # NEW
        'Breakstd': break_std,            # NEW
        'BreakLocs': break_locations,
        'BreakDepth': break_depth
    }
```

### 2. Modify cross_correlation.py

Return ALL CrossCorrelation outputs:
```python
def analyze_timestack(self, timestack):
    """Match MATLAB CrossCorrelation_CoastCams exactly."""
    return {
        'Cf1': celerities_array,          # Raw celerities
        'WLe1': wavelengths_array,        # Raw wavelengths
        'Tp1': periods_array,             # Periods
        'Hs1': heights_array,             # Heights
        'RM': correlation_matrix          # Full correlation matrix
    }
```

### 3. Restructure main.py

Follow MATLAB sequence exactly:
```python
for img_idx in range(num_images):
    # G2: Load image
    img = load_image(img_idx)

    # G3: Shoreline
    shoreline = detect_shoreline(img)

    # G4: Image formatting
    dx, dxXX, sLimit = setup_coordinates(nc)

    # G5: Preprocessing
    S2 = preprocess_image(S1)

    # G6: WaveParameters_CoastCams
    wave_params = analyze_timestack(timestack)
    C = wave_params['C']
    depth_s = wave_params['depth_s']
    depth_l = wave_params['depth_l']
    hs = wave_params['hs']
    # ... store all outputs

    # Store depth from Radon method
    WaterDepth[i,:] = movmean(depth_s/10, 10)

    # G7: CrossCorrelation
    corr_results = cross_correlation(S2)
    Cf1 = corr_results['Cf1']
    WLe1 = corr_results['WLe1']

    # Store celerity from cross-correlation
    WaveCelerity[i,:] = movmean(Cf1/10, 10)
    WaveLength[i,:] = movmean(WLe1, 10)

    # G8: LinearC using CrossCorr celerity
    df = LinearC(Tp, WaveCelerity[i,:], 0.01)
    WaterDepth_L[i,:] = movmean(df, 10)

    # G9: SLA calculations (after loop)
    # ... calculate SLA_S and SLA_L from matrices
```

## Output Variables

After loop, aggregate into timetable with ONE ROW per timestack:

```python
{
    'BreakPointLocation': scalar per timestack,
    'BreakpointDepth': scalar per timestack,
    'Hs': scalar per timestack,
    'WaveEnergy': nanmean of Er array,
    'RollerLength': nanmean of rolL array,
    'WaveCelerity': nanmean of WaveCelerity array,
    'Tp': scalar per timestack,
    'Tm': scalar per timestack,
    'WaveLength': nanmean of WaveLength array,
    'WaterDepth': nanmean of WaterDepth_L array,
    'ShorelinePosition': scalar per timestack,
    'SLA_S': nanmean of SLA_S matrix row,
    'SLA_L': nanmean of SLA_L matrix row,
    'RTR': scalar per timestack,
    'Bathymetry': nanmean(WaterDepth - SLA_S, axis=0)
}
```

## Implementation Priority

1. ✅ First: Get CrossCorrelation returning Cf1, WLe1 properly
2. ✅ Second: Store WaveCelerity and WaveLength matrices
3. ✅ Third: Calculate WaterDepth_L using LinearC(Tp, WaveCelerity)
4. ✅ Fourth: Calculate SLA_S from WaveCelerity matrix
5. ✅ Fifth: Calculate SLA_L from WaterDepth_L matrix
6. ⚠️ Sixth: Extract depth_s from WaveParameters (if needed)

This will ensure all MATLAB outputs are captured correctly.
