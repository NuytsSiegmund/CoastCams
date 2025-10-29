# CoastCams MATLAB vs Python Comparison Report

**Date:** October 29, 2025
**Purpose:** Identify discrepancies and issues between MATLAB and Python implementations

---

## Executive Summary

The Python port of CoastCams provides a **simplified, modular architecture** but is **missing several critical analysis features** present in the MATLAB version. The Python version implements basic wave analysis but lacks advanced functionalities including:

- Radon transform-based wave separation
- Roller property detection and breaking wave analysis
- Complex wave height calculations using camera geometry
- Several specialized filtering and preprocessing methods
- Detailed breakpoint location detection

**Recommendation:** The Python version requires significant enhancements to achieve feature parity with MATLAB.

---

## 1. Architecture Comparison

### MATLAB Implementation
- **Type:** Monolithic script with function library
- **Main file:** `S01_AnalysisTimestackImages.m` (463 lines)
- **Core module:** `WaveParameters_CoastCams.m` (2,053 lines - contains most analysis logic)
- **Approach:** Sequential processing with extensive try-catch error handling
- **Code organization:** Single main script + 10 supporting functions in Subfunctions/

### Python Implementation
- **Type:** Object-oriented, modular package
- **Main file:** `main.py` (390 lines)
- **Architecture:** Separate classes for each component (11 modules)
- **Approach:** Class-based with workflow orchestration
- **Code organization:** Package structure with separate modules for each feature

**Assessment:** ✅ Python has better code organization but ❌ at the cost of missing functionality

---

## 2. Major Missing Features in Python

### 2.1 Wave Parameter Analysis

#### MATLAB Has (Python Missing):
1. **Radon Transform Processing**
   - `RadonCIndiv_20140903()` - Advanced celerity computation using Radon transform
   - `RadonSeparationmodif()` - Separates incident and reflected wave components
   - `FiltreRadon()` - Radon-based filtering between specific angles
   - **Location:** UserScripts/Subfunctions/WaveParameters_CoastCams.m:457-540, 1982-2046

2. **Roller Properties Detection**
   - `RollerPropertiesTaller()` - Detects breaking wave rollers from timestack
   - Computes roller position, length, and temporal extent
   - Uses image preprocessing with Radon separation
   - **Location:** UserScripts/Subfunctions/WaveParameters_CoastCams.m:283-390

3. **Advanced Wave Height Calculation**
   - `BreakerHeight()` - Calculates wave height using camera geometry
   - Accounts for camera angle and wave face angle at breaking
   - Applies geometric corrections based on camera height and position
   - **Location:** UserScripts/Subfunctions/WaveParameters_CoastCams.m:392-455

4. **Breaking Point Analysis**
   - Detects breakpoint location using roller standard deviation
   - Computes breaking depth at the identified location
   - Provides multiple breakpoint metrics (mean, std)
   - **Location:** UserScripts/Subfunctions/WaveParameters_CoastCams.m:38-84

#### Python Has (Simplified):
- Basic wave height from zero-crossing method (wave_analysis.py:129-159)
- Simple wave period computation (wave_analysis.py:161-189)
- Basic spectral analysis using Welch method (wave_analysis.py:191-227)

**Impact:** ⚠️ **CRITICAL** - Python cannot perform breaking wave analysis or accurate wave height estimation

---

### 2.2 Image Preprocessing

#### MATLAB Implementation:
```matlab
function [B]=ImagePreProcessing_20090121Taller(A,icmin,icmax,dt,resc,methodPreT)
```
- **6 different preprocessing methods** (methodPreT = 0-5+)
- Method 1: FIR filters (order=1000) with specific cutoffs (1.5s low, 20s high)
- Method 2: Adaptive filtering based on wave characteristics
- Method 3: Multi-stage smoothing with detrending
- Method 4: Zero-averaging with ideal filter
- Method 5: Hilbert-based normalization
- **Location:** UserScripts/Subfunctions/ImagePreProcessing_CoastCams.m:1-241

#### Python Implementation:
```python
class ImagePreprocessor:
```
- **Single preprocessing approach** using Butterworth filters
- Low-pass: 0.05 Hz (config)
- High-pass: 0.005 Hz (config)
- Filter order: 5 (much lower than MATLAB's 1000)
- **Location:** CoastCams_Python/coastcams/preprocessing.py

**Issues:**
1. ❌ Python lacks multiple preprocessing methods
2. ❌ Different filter types (FIR vs Butterworth) may produce different results
3. ❌ Much lower filter order may affect signal quality
4. ❌ Missing Radon-based preprocessing

---

### 2.3 Cross-Correlation Analysis

#### MATLAB Implementation:
```matlab
function [R2M, L2M, T2M, Hs, RM] = CrossCorrelation_Coastcams(A2, dpha, dt, dc)
```
- Computes correlation matrix R2M(nc-dc+1, dc-1) for all spatial positions
- Extracts wavelength (L2M) for each position
- Computes wave period (T2M) and significant wave height (Hs) per position
- Uses `Wave_Char()` function to extract detailed wave characteristics
- **Location:** UserScripts/Subfunctions/CrossCorrelation_CoastCams.m

#### Python Implementation:
```python
class CrossCorrelationAnalyzer:
```
- Simplified correlation using scipy.signal.correlate
- Computes mean celerity and wavelength
- Missing detailed spatial correlation matrix
- Does not integrate with Wave_Char equivalent
- **Location:** CoastCams_Python/coastcams/cross_correlation.py

**Issues:**
1. ❌ Python doesn't compute full correlation matrix (R2M)
2. ❌ Missing per-position wavelength extraction
3. ❌ No integration with detailed wave characteristic analysis
4. ⚠️ Different correlation approach may yield different results

---

### 2.4 Utility Functions

#### MATLAB Has (Python Missing):
1. **Wave_Char()** - Comprehensive wave characteristic extraction
   - Multiple methods (meth = 1 or 2)
   - Computes hs, htiers, hrms, trms, Hmax, h, Tp, t
   - **Referenced in:** WaveParameters_CoastCams.m:124, CrossCorrelation_CoastCams.m:56,69

2. **smooth2()** - 2D smoothing with edge handling
   - **Location:** UserScripts/Subfunctions/smooth2.m:1-3

3. **lmax()** - Local maxima detection with filtering
   - **Location:** UserScripts/Subfunctions/lmax.m:1-153
   - Also in ImagePreProcessing_CoastCams.m:90-153

4. **FilterMean()** - Mean filtering utility
   - **Location:** UserScripts/Subfunctions/FilterMean.m

5. **LinearC()** - Depth calculation using linear wave theory
   - Iterative solver for dispersion relation
   - **Location:** UserScripts/Subfunctions/LinearC.m:1-19
   - Also in WaveParameters_CoastCams.m:135-156

#### Python Has (Partial):
- `utils.py` contains some utilities but missing many MATLAB functions
- Has `local_maxima()` but may not match lmax() exactly
- Has `smooth2d()` but implementation may differ
- Missing `Wave_Char()` entirely

**Impact:** ⚠️ Missing utility functions means many analysis steps cannot be performed

---

## 3. Algorithm Differences

### 3.1 Wave Period Calculation

#### MATLAB (WaveParameters_CoastCams.m:87-133):
```matlab
function [trms,Tp]=Get_Periode(S,dt)
```
- Uses FIR filters (order 1000) for band-pass filtering
- Applies Welch PSD estimation
- Uses Wave_Char() for detailed period extraction
- Returns both mean period (trms) and peak period (Tp)

#### Python (wave_analysis.py:161-227):
```python
def _compute_wave_periods(self, timeseries: np.ndarray) -> List[float]:
```
- Uses zero-crossing method
- Welch method for peak period only
- No FIR filtering
- Simpler implementation

**Consequence:** May produce different period values

---

### 3.2 Wave Celerity Calculation

#### MATLAB:
```matlab
C = smooth(abs(RadonCIndiv_20140903(dt, dx, B, Tm)),3)';
```
- Uses Radon transform to compute celerity
- Accounts for wave period in the Radon processing
- Applies smoothing with window size 3

#### Python:
```python
def compute_wave_celerity(self, timestack: np.ndarray, dx: float) -> np.ndarray:
```
- Uses simple cross-correlation between successive time steps
- No Radon transform
- Different algorithm entirely

**Consequence:** ⚠️ **MAJOR DIFFERENCE** - Celerity values will likely be significantly different

---

### 3.3 Bathymetry Estimation

#### MATLAB (S01_AnalysisTimestackImages.m:310-380):
```matlab
Bathymetry_full = WaterDepth_adjusted - SLA_S_adjusted;
Bathymetry = nanmean(Bathymetry_full, 1);
```
- Uses two depth calculation methods:
  1. Shallow water: `depth_s = C.^2/9.81`
  2. Linear wave theory: `LinearC(Tp,C,0.01)`
- Subtracts sea level anomaly (SLA)
- Averages across all timesteps

#### Python (bathymetry.py):
```python
class BathymetryEstimator:
```
- Uses dispersion relation solver
- Iterative Newton-Raphson method
- May not match MATLAB exactly

**Assessment:** ⚠️ Both use linear wave theory but implementation details differ

---

## 4. Workflow Differences

### MATLAB Workflow:
1. Load image → Rotate
2. Extract shoreline (shoreline_position.m)
3. Image formatting and preprocessing
4. **Call WaveParameters_CoastCams** (does most analysis)
   - Detects rollers and breaking waves
   - Computes wave heights using camera geometry
   - Estimates celerity with Radon transform
   - Calculates depths
5. Cross-correlation (CrossCorrelation_CoastCams.m)
6. Bathymetry calculation
7. Sea level anomaly
8. Output generation

### Python Workflow:
1. Load images → Rotate
2. Detect shorelines (shoreline.py)
3. Create timestack
4. Analyze waves (wave_analysis.py) - **simplified**
5. Cross-correlation (cross_correlation.py) - **simplified**
6. Estimate bathymetry
7. Compute sea level anomaly
8. Export results

**Key Difference:** Python skips the complex wave parameter analysis done by `WaveParameters_CoastCams`

---

## 5. Configuration Differences

### MATLAB Parameters:
```matlab
dt = 1/2;           % 2 Hz acquisition
H_camera = 27.240;  % Camera height
res = 0.1;          % 0.1 m/pixel
rotation = 270;
dur = 14;           % 14 minutes
Nlim = 1600;        % Timestack width
dpha = 1;           % Time lag for correlation
dc = 100;           % Correlation spacing
```

### Python Parameters (config.yaml):
```yaml
acquisition_frequency: 2.0
camera_height: 27.24
pixel_resolution: 0.1
rotation_angle: 270
cross_shore_width: 1600
time_lag: 1
correlation_spacing: 100
```

**Assessment:** ✅ Configuration values match between versions

---

## 6. Input/Output Differences

### MATLAB Outputs:
- Timetable with 15 variables:
  - BreakPointLocation, BreakpointDepth
  - Hs, WaveEnergy, RollerLength
  - WaveCelerity, Tp, Tm, WaveLength
  - WaterDepth, ShorelinePosition
  - SLA_S, SLA_L, RTR, Bathymetry
- Saves .mat workspace file
- Exports to .txt file
- Generates multiple plots

### Python Outputs:
- CSV file with subset of parameters:
  - Timestamp
  - ShorelinePosition
  - WaveHeight, WavePeriod
  - SeaLevelAnomaly
- Missing: BreakpointLocation, BreakpointDepth, RollerLength, WaveEnergy
- Generates plots

**Issues:** ❌ Python outputs are incomplete compared to MATLAB

---

## 7. Specific Code Issues Found

### Issue 1: Missing Radon Transform Implementation
- **File:** wave_analysis.py
- **Problem:** No Radon transform functions
- **Impact:** Cannot separate incident/reflected waves or compute accurate celerity
- **MATLAB equivalent:** WaveParameters_CoastCams.m:457-540, 1982-2046

### Issue 2: Incomplete Wave Height Calculation
- **File:** wave_analysis.py:129-159
- **Problem:** Uses simple zero-crossing, doesn't account for camera geometry
- **Impact:** Wave heights will be less accurate, especially near breaking
- **MATLAB equivalent:** WaveParameters_CoastCams.m:392-455

### Issue 3: Missing Breaking Wave Detection
- **File:** wave_analysis.py (entire module)
- **Problem:** No roller detection or breakpoint analysis
- **Impact:** Cannot identify where waves break or compute breaking parameters
- **MATLAB equivalent:** WaveParameters_CoastCams.m:283-390

### Issue 4: Simplified Cross-Correlation
- **File:** cross_correlation.py:42-70
- **Problem:** Doesn't compute full spatial correlation matrix
- **Impact:** Less detailed spatial wave analysis
- **MATLAB equivalent:** CrossCorrelation_CoastCams.m:36-50

### Issue 5: Missing Wave_Char Function
- **File:** utils.py
- **Problem:** Critical wave characteristic extraction function not implemented
- **Impact:** Cannot extract detailed wave statistics
- **MATLAB equivalent:** Called in WaveParameters_CoastCams.m:124, Wave_Char.m

### Issue 6: Different Filter Types
- **File:** preprocessing.py:134-183
- **Problem:** Uses Butterworth instead of FIR filters
- **Impact:** Frequency response differs, may affect wave detection
- **MATLAB equivalent:** ImagePreProcessing_CoastCams.m:16-32

### Issue 7: Missing LinearC Iterative Solver
- **File:** utils.py
- **Problem:** May have simplified linear wave theory implementation
- **Impact:** Depth calculations may be less accurate
- **MATLAB equivalent:** WaveParameters_CoastCams.m:135-156

### Issue 8: Incomplete Output Variables
- **File:** main.py:279-321
- **Problem:** Only exports subset of MATLAB variables
- **Impact:** Users lose access to important parameters
- **MATLAB equivalent:** S01_AnalysisTimestackImages.m:381-418

---

## 8. Performance Considerations

### MATLAB:
- Runs in single-threaded mode
- Extensive use of for loops
- Heavy matrix operations (Radon transform)
- **Estimated processing time:** 5-10 minutes for 16 images (based on code complexity)

### Python:
- Simpler algorithms → faster execution
- Uses NumPy vectorization
- Parallel processing option (config: parallel_processing: false)
- **Estimated processing time:** 1-3 minutes for 16 images

**Trade-off:** Python is faster but less comprehensive

---

## 9. Testing and Validation

### MATLAB:
- Extensive error handling with try-catch blocks
- Input validation in shoreline_position.m
- Warning messages for data issues
- Error logging to file (S01_AnalysisTimestackImages.m:420-434)

### Python:
- Some error handling
- Configuration validation
- Try-catch in main workflow
- No error logging file

**Assessment:** ⚠️ MATLAB has more robust error handling

---

## 10. Recommendations

### Critical (Must Fix):
1. ❌ **Implement Radon Transform Processing**
   - Port `RadonCIndiv_20140903`, `RadonSeparationmodif`, `FiltreRadon`
   - Required for accurate wave celerity calculation

2. ❌ **Implement Roller Detection**
   - Port `RollerPropertiesTaller` function
   - Critical for breaking wave analysis

3. ❌ **Implement Wave_Char Function**
   - Port complete wave characteristic extraction
   - Used by multiple modules

4. ❌ **Add Breaking Wave Height Calculation**
   - Port `BreakerHeight` with camera geometry
   - Necessary for accurate wave height estimates

5. ❌ **Complete Cross-Correlation Implementation**
   - Add full spatial correlation matrix computation
   - Match MATLAB's detailed approach

### Important (Should Fix):
6. ⚠️ **Add Multiple Preprocessing Methods**
   - Port all 6 preprocessing methods from MATLAB
   - Allow user to select method

7. ⚠️ **Implement FIR Filters**
   - Add FIR filter option to match MATLAB exactly
   - Important for signal processing consistency

8. ⚠️ **Add Complete Output Variables**
   - Export all parameters that MATLAB exports
   - Include RollerLength, WaveEnergy, BreakpointLocation, etc.

### Nice to Have:
9. ✓ Add error logging file
10. ✓ Add more comprehensive input validation
11. ✓ Improve progress reporting for long-running analysis

---

## 11. Summary Table

| Feature | MATLAB | Python | Status |
|---------|---------|---------|---------|
| Radon Transform | ✅ Full implementation | ❌ Missing | **Critical** |
| Roller Detection | ✅ Implemented | ❌ Missing | **Critical** |
| Wave_Char Function | ✅ Comprehensive | ❌ Missing | **Critical** |
| Breaking Wave Analysis | ✅ With geometry | ❌ Simple only | **Critical** |
| Cross-Correlation Matrix | ✅ Full spatial | ⚠️ Simplified | **Major** |
| Preprocessing Methods | ✅ 6 methods | ⚠️ 1 method | **Major** |
| Filter Type | ✅ FIR (order 1000) | ⚠️ Butterworth (order 5) | **Major** |
| Wave Height Calculation | ✅ Geometry-based | ⚠️ Zero-crossing only | **Major** |
| Wave Period Calculation | ✅ Multiple methods | ⚠️ Basic methods | **Minor** |
| Bathymetry Estimation | ✅ Implemented | ✅ Implemented | **OK** |
| Shoreline Detection | ✅ 3 methods | ✅ 3 methods | **OK** |
| Sea Level Anomaly | ✅ Implemented | ✅ Implemented | **OK** |
| Configuration Management | ✅ In-script | ✅ YAML file | **OK (Better)** |
| Code Organization | ⚠️ Monolithic | ✅ Modular | **OK (Better)** |
| Output Variables | ✅ 15 variables | ⚠️ 5 variables | **Major** |
| Error Handling | ✅ Extensive | ⚠️ Basic | **Minor** |
| Documentation | ⚠️ Comments | ✅ Docstrings | **OK (Better)** |

---

## 12. Conclusion

The Python version of CoastCams provides a **clean, modular architecture** with better code organization than MATLAB. However, it is currently a **simplified implementation** that lacks many critical analysis features:

**Missing Critical Functionality (30-40% feature parity):**
- No Radon transform processing
- No roller/breaking wave detection
- Simplified wave height calculation
- Incomplete cross-correlation analysis
- Missing many utility functions

**Implementation Differences:**
- Different filtering approaches (Butterworth vs FIR)
- Different wave analysis algorithms
- Reduced output variables

**Python Advantages:**
- Better code organization and modularity
- Configuration file support
- Better documentation (docstrings)
- Potentially faster (but less comprehensive)

**Recommendation:**
The Python version requires **significant development work** to achieve functional parity with the MATLAB version. The missing features are not minor - they represent core analysis capabilities of CoastCams. A staged development approach is recommended:

1. **Phase 1 (Critical):** Implement Radon transform, roller detection, Wave_Char
2. **Phase 2 (Major):** Complete cross-correlation, add preprocessing methods
3. **Phase 3 (Polish):** Match all output variables, improve error handling

**Estimated effort:** 2-4 weeks of development for Phase 1, 1-2 weeks each for Phases 2-3.

---

## 13. Detailed Function Mapping

### Functions to Port from MATLAB to Python:

| MATLAB Function | Location | Lines | Python Module | Status |
|----------------|----------|-------|---------------|--------|
| WaveParameters_CoastCams | WaveParameters_CoastCams.m | 1-2054 | wave_analysis.py | ⚠️ Partially implemented |
| RollerPropertiesTaller | WaveParameters_CoastCams.m:283 | 108 | wave_analysis.py | ❌ Missing |
| BreakerHeight | WaveParameters_CoastCams.m:392 | 64 | wave_analysis.py | ❌ Missing |
| RadonCIndiv_20140903 | WaveParameters_CoastCams.m:457 | 84 | cross_correlation.py | ❌ Missing |
| Get_Periode | WaveParameters_CoastCams.m:87 | 47 | wave_analysis.py | ⚠️ Different implementation |
| LinearC | WaveParameters_CoastCams.m:135 | 22 | utils.py | ⚠️ Check implementation |
| ImagePreProcessing_20090121Taller | WaveParameters_CoastCams.m:158 | 47 | preprocessing.py | ❌ Missing methods |
| RadonSeparationmodif | WaveParameters_CoastCams.m:1982 | 25 | preprocessing.py | ❌ Missing |
| FiltreRadon | WaveParameters_CoastCams.m:2009 | 37 | preprocessing.py | ❌ Missing |
| smooth2 | smooth2.m | 4 | utils.py | ✅ Implemented |
| lmax | lmax.m | 153 | utils.py | ⚠️ Check implementation |
| Wave_Char | Wave_Char.m | 201 | wave_analysis.py | ❌ Missing |
| FilterMean | FilterMean.m | 57 | utils.py | ❌ Missing |
| CrossCorrelation_CoastCams | CrossCorrelation_CoastCams.m | 94 | cross_correlation.py | ⚠️ Simplified |
| shoreline_position | shoreline_position.m | 165 | shoreline.py | ✅ Implemented |

---

**End of Report**
