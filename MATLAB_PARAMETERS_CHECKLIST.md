# CoastCams: MATLAB vs Python Parameters Checklist

## All MATLAB Variables Calculated ✅

| MATLAB Variable | Python Equivalent | Status | Notes |
|----------------|------------------|--------|-------|
| `Img_date` | `timestamps` | ✅ Calculated | From image filenames |
| `WaveCelerity` | `celerities` | ✅ Calculated | From MATLAB-style cross-correlation |
| `df` / `WaterDepth_L` | `depths` | ✅ Calculated | From LinearC (linear wave theory) |
| `WaveLength` | `wavelengths` | ✅ Calculated | L = C × T |
| `Hs_TS` | `wave_heights_timeseries` | ✅ Calculated | Significant wave height |
| `Tp_TS` | `wave_periods_peak_timeseries` | ✅ Calculated | Peak period |
| `Tm_TS` | `wave_periods_timeseries` | ✅ Calculated | Mean period |
| `WaveEnergy` | `wave_energies` | ✅ Calculated | E = Hs²/2 |
| `RollerLength` | `roller_lengths` | ✅ Calculated | From roller detection |
| `BreakpointDepth` | `breakpoint_depths` | ✅ Calculated | From breaking wave analysis |
| `BreakpointLocation` | `breakpoint_locations` | ✅ Calculated | Breaking wave position |
| `Cf1` | `celerities` | ✅ Calculated | Same as WaveCelerity |
| `WLe1` | `wavelengths` | ✅ Calculated | Same as WaveLength |
| `Depth_S` | `depths_shallow` | ✅ Calculated | C²/g (shallow water approx) |
| `ShorePosition` | `shoreline_positions` | ✅ Calculated | From shoreline detection |
| `Stack_av` | `average_timestack` | ✅ Calculated | Average timestack matrix |
| `BL_Coordinates` | Part of roller detection | ✅ Calculated | Breaking line coordinates |
| `Time_TS` | `timestamps` | ✅ Calculated | Time array |
| `SLA_S` | `sla_shallow_values` | ✅ Calculated | Shallow water SLA |
| `SLA_L` | `sla_values` | ✅ Calculated | Linear wave theory SLA |
| `RTR` | `rtr_values` | ✅ Calculated | Relative Tidal Range |

## Timetable Columns ✅

Python DataFrame with datetime index matching MATLAB timetable format:

```python
Columns = {
    'BreakPointLocation': breakpoint_locations,
    'BreakpointDepth': breakpoint_depths,
    'Hs': wave_heights,
    'WaveEnergy': wave_energies,
    'RollerLength': roller_lengths,
    'WaveCelerity': celerities,
    'Tp': peak_periods,
    'Tm': mean_periods,
    'WaveLength': wavelengths,
    'WaterDepth': depths,
    'ShorelinePosition': shoreline_positions,
    'SLA_S': sla_shallow_values,
    'SLA_L': sla_values,
    'RTR': rtr_values,
    'Bathymetry': depths  # Same as WaterDepth
}
```

**Output Format**: CSV file with Timestamp as index (pandas DataFrame with `set_index('Timestamp')`)

## Plots Generated ✅

### MATLAB-style Summary Plot
Function signature matching MATLAB:
```matlab
% MATLAB
plot_coastcams_main(Time_TS, Stack_av, SLA_S, Hs_TS, Tp_TS, rotation)
```

```python
# Python equivalent
visualizer.plot_matlab_style_summary(
    timestamps=timestamps,           # Time_TS
    average_timestack=Stack_av,      # Stack_av
    sla_matrix=SLA_S,               # SLA_S (shallow water)
    wave_heights=Hs_TS,             # Hs_TS
    wave_periods=Tp_TS,             # Tp_TS (using Tm if Tp unavailable)
    rotation=rotation               # rotation angle
)
```

**4-Panel Layout**:
1. Average timestack (grayscale, rotated)
2. SLA heatmap (jet colormap)
3. Wave height time series (red line)
4. Wave period time series (blue line)

## Key Calculations

### 1. Wavelength
```python
wavelength = celerity * mean_period
```

### 2. Wave Energy
```python
wave_energy = Hs² / 2
```
(Simplified form of E = (1/16) × ρ × g × Hs²)

### 3. Shallow Water Depth
```python
depth_shallow = celerity² / 9.81
```

### 4. RTR (Relative Tidal Range)
```python
water_levels = celerities
min_water_level = nanmin(water_levels)
nRTR_thresh = water_levels - min_water_level
nRTR_thresh[nRTR_thresh < 0.2] = NaN  # Avoid division by near-zero
RTR = Hs / nRTR_thresh
```

### 5. SLA_S (Shallow Water)
```python
# Smooth celerity matrix spatially (Nc=30)
celerity_matrix_smooth = uniform_filter1d(celerity_matrix, size=30, axis=1)
# Calculate SLA from smoothed matrix
SLA_S_matrix = celerity_matrix_smooth - nanmean(celerity_matrix_smooth)
# Average across space for time series
SLA_S = nanmean(SLA_S_matrix, axis=1)
```

### 6. SLA_L (Linear Wave Theory)
```python
# Apply movmean(10) then smooth2(Nr=1, Nc=30) to depth matrix
depth_matrix_ma = uniform_filter1d(depth_matrix, size=10, axis=1)
depth_matrix_smooth = uniform_filter1d(depth_matrix_ma, size=30, axis=1)
# Calculate SLA from smoothed matrix
SLA_L_matrix = depth_matrix_smooth - nanmean(depth_matrix_smooth)
# Average across space for time series
SLA_L = nanmean(SLA_L_matrix, axis=1)
```

## Implementation Files

- **`main.py`**: Main workflow with all parameter calculations
- **`cross_correlation.py`**: MATLAB-style cross-correlation for celerities
- **`bathymetry.py`**: Depth estimation using linear wave theory
- **`wave_analysis.py`**: Wave height, period, and breaking parameters
- **`matlab_preprocessing.py`**: Roller detection and photogrammetric heights
- **`visualize.py`**: MATLAB-style summary plots

## Output Files

1. **CSV File**: `coastcams_results_YYYYMMDD_HHMMSS.csv`
   - Contains all timetable columns
   - Timestamp index
   - All parameters per timestack

2. **Plots**:
   - `coastcams_matlab_summary.png`: 4-panel MATLAB-style plot
   - Additional bathymetry, wave parameter plots

## Validation

To validate against MATLAB:
1. Run both versions on same dataset
2. Compare CSV outputs column by column
3. Check SLA_L achieves R² = 0.98 vs pressure sensor data
4. Verify no NaN values in SLA (critical requirement)
