# Investigation: Empty CSV Columns (WaterDepth, SLA_S, SLA_L, Bathymetry)

## Problem Statement

User reports that several CSV columns are still empty:
- `WaterDepth`: EMPTY
- `SLA_S`: EMPTY
- `SLA_L`: EMPTY
- `Bathymetry`: EMPTY

Yet the data shows that many other columns ARE populated correctly (Hs, Tp, Tm, WaveCelerity, etc.)

## Root Causes Identified

### 1. Camera Angle Issue (NOW FIXED)

**Problem**: I incorrectly "fixed" the camera angle to be constant, but MATLAB correctly calculates a different angle for each wave:

```matlab
AngleCam = abs(z0./X1(PosX));
```

- `z0` (camera height) is FIXED
- `X1(PosX)` (distance of each breaking wave) VARIES
- Therefore `AngleCam` varies for each wave based on where it breaks

**Fix**: Reverted to calculating angle for each wave individually based on its distance from camera.

### 2. Tp Calculation (NEEDS INVESTIGATION)

**MATLAB Code**:
```matlab
Pv = pwelch(detrend(S-mean(S)));
[valmx, indmx] = max(Pv);
Tp = 1./(indmx./(dt*length(Pv)));
```

**Current Python**: Uses scipy's welch with different parameters. Need to verify this matches MATLAB's pwelch exactly.

### 3. CSV Export Keys (NEEDS VERIFICATION)

Looking at `main.py` line 852:
```python
'WaterDepth': self.results.get('depths', []),
```

And lines 854-855:
```python
'SLA_S': self.results.get('sla_shallow_values', []),
'SLA_L': self.results.get('sla_values', []),
```

These keys SHOULD be populated in `_aggregate_all_results` (lines 634, 639-640):
```python
'depths': water_depths_array,
'sla_values': slas,  # SLA_L
'sla_shallow_values': sla_shallow,  # SLA_S
```

**Hypothesis**: Perhaps `water_depths_array`, `slas`, or `sla_shallow` are all NaN arrays?

## What Data IS Being Calculated

From the CSV provided by user:
- ✅ BreakPointLocation: 17.4, 17.2, 17.6, etc. (POPULATED)
- ✅ BreakpointDepth: 0.668, 0.656, 0.651, etc. (POPULATED)
- ✅ Hs: 0.479, 0.377, 0.364, etc. (POPULATED)
- ✅ WaveEnergy: 0.115, 0.071, 0.066, etc. (POPULATED)
- ✅ RollerLength: 1.764, 1.822, 1.735, etc. (POPULATED)
- ✅ WaveCelerity: 2.561, 2.537, 2.527, etc. (POPULATED)
- ✅ Tp: 14.22, 10.67, 10.67, etc. (POPULATED and DIFFERENT from Tm!)
- ✅ Tm: 8.78, 6.96, 7.12, etc. (POPULATED)
- ❌ WaterDepth: EMPTY
- ❌ SLA_S: EMPTY
- ❌ SLA_L: EMPTY
- ❌ Bathymetry: EMPTY

## Investigation Steps

### Step 1: Check if `water_depths_array` is calculated

In `_aggregate_all_results`, line 430:
```python
water_depths_array = np.nanmean(depth_matrix_smooth, axis=1)
```

This depends on `depth_profiles_list` being populated. Each result in `all_results` needs:
```python
r.get('depth_profile') is not None
```

Which is stored in line 291 of main loop:
```python
'depth_profile': WaterDepth_L,
```

Where `WaterDepth_L` is calculated in lines 189-210.

**Verification needed**: Are the `WaterDepth_L` arrays being stored correctly in `all_results`?

### Step 2: Check if `slas` is calculated

Lines 419-427 calculate SLA:
```python
sla_matrix = depth_matrix_smooth - np.nanmean(depth_matrix_smooth)
slas = np.nanmean(sla_matrix, axis=1)
```

This depends on the same `depth_matrix` as water_depths_array.

### Step 3: Check if `sla_shallow` is calculated

Lines 561-587 calculate SLA_S from WaveCelerity matrix:
```python
if len(wave_celerity_arrays) > 0:
    WaveCelerity_matrix = np.array([arr[:min_len] for arr in wave_celerity_arrays])
    # ... smoothing ...
    sla_shallow = np.nanmean(SLA_S_matrix, axis=1)
else:
    sla_shallow = np.full(len(all_results), np.nan)
```

This depends on `wave_celerity_arrays` being populated from:
```python
r.get('wave_celerity_array') is not None
```

## Likely Issue

**Theory**: The issue is probably NOT in the calculation or export code, but in the CSV file being viewed.

The user might be looking at an OLD CSV file that was generated BEFORE the recent fixes. The newest fixes include:
1. Matrix building workflow (commit bf5395f)
2. Tp/Tm separation (commit 123677e)
3. Camera angle revert (commit 4694801)

## Next Steps

1. **Run a fresh analysis** with the latest code
2. **Check the OUTPUT timestamp** to ensure it's from the latest run
3. **Add debug output** to verify arrays are being stored:
   - Print `len(depth_profiles_list)` in _aggregate_all_results
   - Print `len(wave_celerity_arrays)`
   - Print shapes of `water_depths_array`, `slas`, `sla_shallow`
4. **Verify CSV export** is reading from correct keys

## MATLAB Reference

From the MATLAB code provided, the key outputs are:
```matlab
function [C,depth_s,depth_l,hs,hm,Tm,Tp,Er,rolL,nbwave,...
          Breakstd,Breakmean1,Breakmean2,BreakLocs, BreakDepth]= ...
          WaveParameters_CoastCams(Img,dt,CoordCam,dx,dur)
```

Where:
- `depth_l` = from LinearC(Tp, C, 0.01) - linear wave theory depth
- `depth_s` = C^2/9.81 - shallow water depth
- Both are ARRAYS along the spatial dimension
- Final CSV should have scalar values (mean/nanmean of arrays)

## Breakpoint Location Variation

User mentioned breakpoint locations are "too similar" (all around 17-18m). Looking at MATLAB:
```matlab
BreakLocs = str2num(B(B~=' '));
```

Where:
```matlab
B = num2str(round(X1((round(nanmean(PosX))))));
```

So `BreakLocs` is based on `X1` (cross-shore coordinate) at position `round(nanmean(PosX))`.

The variation comes from:
1. Different `PosX` (breaking positions in pixels) for each timestack
2. Different `X1` (cross-shore coordinates) mapping

Python needs to ensure `PosX` is varying appropriately based on wave conditions.
