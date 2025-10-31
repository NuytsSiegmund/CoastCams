# Critical Fixes to Match MATLAB Exactly

## Summary of Remaining Issues

From your data showing empty fields, the problems are:

1. **SLA_S is empty** - WaveCelerity matrix not being built correctly
2. **WaterDepth has only 4/16 values** - LinearC needs the WaveCelerity ARRAY, not scalar
3. **Some WaveEnergy/RollerLength empty** - Need full arrays from WaveParameters

## Root Cause: Python Aggregates Too Early

**MATLAB stores MATRICES during the loop, aggregates AFTER:**

```matlab
for i = 1:length(Img)
    % Store FULL ARRAYS (not scalars)
    WaveCelerity(i,1:size(Cf1,2)) = movmean(Cf1./10, 10);
    WaveLength(i,1:size(WLe1,2)) = movmean(WLe1,10);
    WaterDepth_L(i,:) = movmean(df, 10);
end

% AFTER loop: aggregate to scalars for timetable
for i = 1:num_rows
    BreakpointDepth_scalar(i) = BreakpointDepth(i);  % Already scalar
    WaveEnergy_scalar(i) = nanmean(WaveEnergy(i,:)); % Mean of array
    WaveCelerity_scalar(i) = nanmean(WaveCelerity(i,:)); % Mean of array
end
```

**Python currently:**
```python
# Wrong: stores scalars immediately
result = {
    'wave_celerity': correlation_results.get('mean_celerity', np.nan),  # ❌ SCALAR
    ...
}
```

## Fix 1: Store Full Arrays in Results

In `main.py`, change storage to match MATLAB:

```python
# Current (WRONG):
result = {
    'wave_celerity': correlation_results.get('mean_celerity', np.nan),  # Scalar
}

# Fixed (CORRECT):
result = {
    'wave_celerity_array': correlation_results.get('celerities', None),  # Full array
    'wave_celerity_mean': np.nanmean(celerities) if celerities is not None else np.nan,
    'wavelength_array': correlation_results.get('wavelengths', None),  # NEW
    'depth_array': bathymetry_results.get('depths_filtered', None),  # Already stored
}
```

## Fix 2: Calculate WaterDepth_L from WaveCelerity Array

**Current Python bug:**
```python
# Using scalar mean - WRONG
peak_period = wave_results.get('mean_Tp', 8.0)  # Scalar
wave_periods_full = np.full(num_positions, peak_period)
bathymetry_results = estimate_depth_profile(wave_periods_full, celerities_full)
```

**MATLAB correct approach:**
```matlab
% Line 247: Use Tp (scalar) and WaveCelerity array
[df] = LinearC(Tp, WaveCelerity(i,:), 0.01);
% Tp is SCALAR, WaveCelerity(i,:) is ARRAY
```

**Fixed Python:**
```python
# Get the full celerity array from cross-correlation
Cf1 = correlation_results.get('celerities', None)  # Full array

# Apply MATLAB transformation: movmean(Cf1/10, 10)
from scipy.ndimage import uniform_filter1d
WaveCelerity = uniform_filter1d(Cf1 / 10.0, size=10, mode='nearest')

# Use LinearC with SCALAR Tp and ARRAY WaveCelerity
Tp = wave_results.get('mean_Tp', wave_results.get('mean_Tm', 8.0))
depths = []
for c in WaveCelerity:
    if not np.isnan(c):
        d = calculate_depth_from_celerity(c, Tp, precision=0.01)
        depths.append(d)
    else:
        depths.append(np.nan)

# Apply movmean again
WaterDepth_L = uniform_filter1d(np.array(depths), size=10, mode='nearest')
```

## Fix 3: Build Matrices for SLA Calculation

**Current Python issue:** SLA_S is calculated but the celerity matrix has too many NaN values.

**MATLAB approach (lines 256-267):**
```matlab
% Build WaveCelerity matrix: rows=timestacks, cols=spatial positions
% WaveCelerity(i, 1:size(Cf1,2)) stores each timestack

% After all timestacks processed:
Csmooth_S = smooth2(WaveCelerity, Nr=1, Nc=30);
SLA_S = Csmooth_S - nanmean(Csmooth_S);  % 2D matrix
% Then for timetable: take nanmean(SLA_S(i,:)) for each row
```

**Fixed Python in _aggregate_all_results():**
```python
# Build WaveCelerity matrix
wave_celerity_matrix = []
for r in all_results:
    if r.get('wave_celerity_array') is not None:
        wave_celerity_matrix.append(r['wave_celerity_array'])

if len(wave_celerity_matrix) > 0:
    # Ensure same length
    min_len = min(len(arr) for arr in wave_celerity_matrix)
    wave_celerity_trimmed = [arr[:min_len] for arr in wave_celerity_matrix]
    WaveCelerity = np.array(wave_celerity_trimmed)  # Shape: (timestacks, space)

    # Apply smooth2 with Nr=1, Nc=30
    Csmooth_S = smooth2(WaveCelerity, Nr=1, Nc=30)

    # Calculate SLA_S
    SLA_S = Csmooth_S - np.nanmean(Csmooth_S)

    # For CSV: take spatial mean for each timestack
    sla_s_values = np.nanmean(SLA_S, axis=1)  # Shape: (timestacks,)
```

## Fix 4: WaveLength Calculation

**MATLAB (line 243):**
```matlab
WaveLength(i,1:size(WLe1,2)) = movmean(WLe1,10);
```

**Python needs to:**
```python
# Get WLe1 from cross-correlation
WLe1 = correlation_results.get('wavelengths_raw', None)  # NEW output needed

# Apply movmean
WaveLength = uniform_filter1d(WLe1, size=10, mode='nearest')

# Store in result
result['wavelength_array'] = WaveLength
result['wavelength_mean'] = np.nanmean(WaveLength)
```

## Implementation Steps

### Step 1: Modify cross_correlation.py

Add raw outputs BEFORE any processing:

```python
def _compute_celerity_matlab_style(self, timestack):
    # ... existing code ...

    # Return RAW celerities (Cf1 equivalent)
    celerities_full[start_idx:end_idx] = velocities_ms  # BEFORE smoothing

    return celerities_full  # Raw Cf1
```

Then in `_extract_wave_properties()`:

```python
def _extract_wave_properties(self, timestack):
    # Get raw Cf1
    Cf1 = self._compute_celerity_matlab_style(timestack)

    # Also compute wavelengths (WLe1)
    WLe1 = self._compute_wavelength(timestack)  # NEW method needed

    return {
        'Cf1': Cf1,  # Raw celerities
        'WLe1': WLe1,  # Raw wavelengths
        'celerities': Cf1,  # For backward compatibility
    }
```

### Step 2: Modify main.py loop

```python
for img_idx in range(num_images):
    # ... existing code ...

    # Get cross-correlation results
    correlation_results = self.correlation_analyzer.analyze_timestack(timestack)
    Cf1 = correlation_results.get('Cf1', None)
    WLe1 = correlation_results.get('WLe1', None)

    # Apply MATLAB transformations
    if Cf1 is not None:
        WaveCelerity = uniform_filter1d(Cf1 / 10.0, size=10, mode='nearest')
    else:
        WaveCelerity = None

    if WLe1 is not None:
        WaveLength = uniform_filter1d(WLe1, size=10, mode='nearest')
    else:
        WaveLength = None

    # Calculate WaterDepth_L using Tp and WaveCelerity array
    Tp = wave_results.get('mean_Tp', wave_results.get('mean_Tm', 8.0))
    if WaveCelerity is not None:
        depths = np.array([
            calculate_depth_from_celerity(c, Tp, precision=0.01)
            for c in WaveCelerity
        ])
        WaterDepth_L = uniform_filter1d(depths, size=10, mode='nearest')
    else:
        WaterDepth_L = None

    # Store ARRAYS (not scalars)
    result = {
        'timestamp': self.image_loader.timestamps[img_idx],
        'wave_celerity_array': WaveCelerity,  # Full array
        'wavelength_array': WaveLength,  # Full array
        'water_depth_array': WaterDepth_L,  # Full array
        # ... other fields
    }
```

### Step 3: Aggregate correctly in _aggregate_all_results()

```python
# Build matrices
WaveCelerity_matrix = []
WaveLength_matrix = []
WaterDepth_L_matrix = []

for r in all_results:
    if r.get('wave_celerity_array') is not None:
        WaveCelerity_matrix.append(r['wave_celerity_array'])
    if r.get('wavelength_array') is not None:
        WaveLength_matrix.append(r['wavelength_array'])
    if r.get('water_depth_array') is not None:
        WaterDepth_L_matrix.append(r['water_depth_array'])

# Convert to numpy arrays (pad/trim to same length)
# ... processing ...

# Calculate SLA_S from WaveCelerity_matrix
Csmooth_S = smooth2(WaveCelerity_matrix, Nr=1, Nc=30)
SLA_S = Csmooth_S - np.nanmean(Csmooth_S)

# Calculate SLA_L from WaterDepth_L_matrix
Csmooth_L = smooth2(WaterDepth_L_matrix, Nr=1, Nc=30)
SLA_L = Csmooth_L - np.nanmean(Csmooth_L)

# For CSV: take spatial means
wave_celerity_scalars = np.nanmean(WaveCelerity_matrix, axis=1)
wavelength_scalars = np.nanmean(WaveLength_matrix, axis=1)
water_depth_scalars = np.nanmean(WaterDepth_L_matrix, axis=1)
sla_s_scalars = np.nanmean(SLA_S, axis=1)
sla_l_scalars = np.nanmean(SLA_L, axis=1)
```

This matches MATLAB's exact workflow!
