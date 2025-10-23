"""
MATLAB-style preprocessing for CoastCams timestacks.

This module implements the full MATLAB preprocessing pipeline from
WaveParameters_CoastCams.m, including:
- ImagePreProcessing_20090121Taller: FIR filtering + Radon transform
- RollerPropertiesTaller: Wave and roller detection
- BreakerHeight: Photogrammetric wave height calculation
"""

import numpy as np
from scipy import signal
from scipy.ndimage import uniform_filter1d
from typing import Tuple, Optional, Dict, List
from skimage.transform import radon, iradon


class MATLABPreprocessor:
    """
    MATLAB-style preprocessing for timestack images.

    Implements the preprocessing pipeline from WaveParameters_CoastCams.m
    """

    def __init__(self, dt: float = 0.5, use_radon: bool = False):
        """
        Initialize preprocessor.

        Parameters
        ----------
        dt : float
            Temporal resolution (sampling period) in seconds
        use_radon : bool
            Whether to apply Radon transform (slow but more accurate)
        """
        self.dt = dt
        self.use_radon = use_radon

    def preprocess_timestack(self, timestack: np.ndarray) -> np.ndarray:
        """
        Preprocess timestack with FIR filtering and Radon transform.

        Matches MATLAB ImagePreProcessing_20090121Taller function (lines 158-204).

        Parameters
        ----------
        timestack : np.ndarray
            Input timestack (time × space), e.g., (1680, 689)

        Returns
        -------
        np.ndarray
            Preprocessed timestack with incident waves separated
        """
        nt, nc = timestack.shape

        # Spatial resolution for sampling (MATLAB line 163)
        # Process every ~1% of spatial points
        resc = max(1, round(nc / 100))

        print(f"  Preprocessing timestack: {nt}×{nc}, sampling every {resc} spatial points")

        # Initialize output
        B = np.zeros_like(timestack)

        # Loop over spatial grid (MATLAB lines 169-189)
        for ic in range(0, nc, resc):
            # Extract time series at this spatial location
            time_series = timestack[:, ic].copy()

            # Apply bandpass filtering (1.5s - 20s)
            filtered_series = self._apply_fir_bandpass(time_series)

            if filtered_series is not None:
                B[:, ic] = filtered_series

        # Interpolate missing values (MATLAB lines 191-197)
        for it in range(nt):
            row = B[it, :]
            valid_indices = np.where(np.abs(row) > 1e-10)[0]

            if len(valid_indices) > 1:
                # Interpolate to fill gaps
                B[it, :] = np.interp(np.arange(nc), valid_indices, row[valid_indices])

        # Replace NaNs with zeros (MATLAB line 199)
        B[np.isnan(B)] = 0.0

        # Apply Radon transform to separate incident/reflected waves (MATLAB line 202)
        if self.use_radon:
            print("    Applying Radon transform (this may take a while)...")
            B_separated = self._radon_separation(B)
            return B_separated
        else:
            print("    Skipping Radon transform (use_radon=False)")
            return B

    def _apply_fir_bandpass(self, time_series: np.ndarray) -> Optional[np.ndarray]:
        """
        Apply FIR bandpass filter (1.5s low-pass, 20s high-pass).

        Matches MATLAB ImagePreProcessing lines 170-186.

        Parameters
        ----------
        time_series : np.ndarray
            Input time series

        Returns
        -------
        np.ndarray or None
            Filtered time series
        """
        try:
            from scipy.signal import firwin, lfilter

            S = time_series.copy()

            # Low-pass filter (Tcoupure = 1.5s, MATLAB line 171)
            Tcoupure = 1.5
            fr = 1.0 / self.dt
            Val = (1.0 / Tcoupure) * 2.0 * (1.0 / fr)
            ord = 1000

            # Create FIR filter
            fil = firwin(ord + 1, Val, window='hamming')

            # Apply filter via convolution (MATLAB line 176)
            kk1 = np.convolve(fil, S, mode='full')

            # Remove edge effects (MATLAB line 177)
            edge = ord // 2
            if len(kk1) > 2 * edge:
                S_filtered = kk1[edge:len(kk1) - edge]

                # Detrend (MATLAB line 177)
                S_filtered = signal.detrend(S_filtered)

                # Ensure same length as input
                if len(S_filtered) > len(S):
                    S_filtered = S_filtered[:len(S)]
                elif len(S_filtered) < len(S):
                    # Pad if needed
                    pad_len = len(S) - len(S_filtered)
                    S_filtered = np.pad(S_filtered, (0, pad_len), mode='edge')
            else:
                return None

            # High-pass filter (Tcoupure = 20s, MATLAB line 179)
            Tcoupure = 20.0
            Val = (1.0 / Tcoupure) * 2.0 * (1.0 / fr)

            # Create high-pass FIR filter
            fil_high = firwin(ord + 1, Val, window='hamming', pass_zero=False)

            # Apply filter
            kk1 = np.convolve(fil_high, S_filtered, mode='full')

            # Remove edge effects
            if len(kk1) > 2 * edge:
                y = kk1[edge:len(kk1) - edge]

                # Detrend
                y = signal.detrend(y)

                # Ensure same length
                if len(y) > len(S):
                    y = y[:len(S)]
                elif len(y) < len(S):
                    pad_len = len(S) - len(y)
                    y = np.pad(y, (0, pad_len), mode='edge')
            else:
                return S_filtered

            return y

        except Exception as e:
            print(f"    FIR filtering failed: {e}")
            return None

    def _radon_separation(self, timestack: np.ndarray) -> np.ndarray:
        """
        Apply Radon transform to separate incident/reflected waves.

        Matches MATLAB RadonSeparationmodif function (lines 1982-2007).

        Parameters
        ----------
        timestack : np.ndarray
            Preprocessed timestack (time × space)

        Returns
        -------
        np.ndarray
            Timestack with incident waves separated
        """
        nt, nx = timestack.shape

        # Check dimensions (MATLAB line 1994)
        if nx > nt:
            print(f"    Warning: nx > nt ({nx} > {nt}), skipping Radon separation")
            return timestack

        # Pre-treatment: remove mean from each spatial column (MATLAB lines 1999-2001)
        M = timestack.copy()
        for i in range(nx):
            M[:, i] = M[:, i] - np.nanmean(M[:, i])

        # Apply Radon filter (MATLAB line 2003)
        # Filter angles 1-89 degrees (incident waves from offshore)
        Sin = self._filtre_radon(M, 1, 89)

        return Sin

    def _filtre_radon(self, image: np.ndarray, angle_min: int, angle_max: int) -> np.ndarray:
        """
        Filter image using Radon transform.

        Matches MATLAB FiltreRadon function (lines 2009-2045).

        Parameters
        ----------
        image : np.ndarray
            Input image (time × space)
        angle_min : int
            Minimum angle in degrees
        angle_max : int
            Maximum angle in degrees

        Returns
        -------
        np.ndarray
            Filtered image
        """
        c1, c2 = image.shape

        # Perform Radon transform (MATLAB line 2023)
        theta = np.arange(0, 180)
        R = radon(image, theta=theta, circle=False)

        # Reconstruct with specified angle range (MATLAB line 2026)
        theta_filtered = np.arange(angle_min, angle_max + 1)

        # Extract the filtered sinogram
        R_filtered = R[:, angle_min:angle_max + 1]

        # Inverse Radon transform (MATLAB line 2026)
        try:
            I = iradon(R_filtered, theta=theta_filtered, filter_name='hann', circle=False, output_size=max(c1, c2))

            # The reconstructed image is square with size = max(c1, c2)
            # We need to extract the original shape (c1, c2)

            # Get the center of the reconstructed image
            I_size = I.shape[0]
            center = I_size // 2

            # Create output image
            S2 = np.full_like(image, np.nan, dtype=np.float64)

            # Extract the region corresponding to original shape
            if c1 <= c2:
                # More columns: extract centered rows
                row_start = center - c1 // 2
                row_end = row_start + c1
                col_start = center - c2 // 2
                col_end = col_start + c2

                if row_end <= I_size and col_end <= I_size:
                    S2[:, :] = I[row_start:row_end, col_start:col_end]
                else:
                    # Fallback: return original
                    return image
            else:
                # More rows: similar extraction
                row_start = center - c1 // 2
                row_end = row_start + c1
                col_start = center - c2 // 2
                col_end = col_start + c2

                if row_end <= I_size and col_end <= I_size:
                    S2[:, :] = I[row_start:row_end, col_start:col_end]
                else:
                    return image

            # Multiply by 0.5 (MATLAB line 2044)
            S2 = S2 * 0.5

            return S2

        except Exception as e:
            print(f"    Radon reconstruction failed: {e}, returning original")
            return image


class RollerDetector:
    """
    Detect breaking waves and roller properties.

    Implements MATLAB RollerPropertiesTaller function (lines 283-390).
    """

    def __init__(self, dt: float = 0.5):
        """
        Initialize roller detector.

        Parameters
        ----------
        dt : float
            Temporal resolution in seconds
        """
        self.dt = dt
        self.preprocessor = MATLABPreprocessor(dt)

    def detect_rollers(self, timestack: np.ndarray) -> Dict:
        """
        Detect breaking waves and compute roller properties.

        Matches MATLAB RollerPropertiesTaller function.

        Parameters
        ----------
        timestack : np.ndarray
            Raw timestack (time × space)

        Returns
        -------
        Dict
            Contains:
            - PosX: Breaking positions (spatial indices)
            - PosT: Breaking times (temporal indices)
            - Lw: Roller lengths
            - B: Preprocessed timestack
            - Breakstd: Standard deviation profile
        """
        I = timestack.astype(np.float64)
        nt, nc = I.shape

        print(f"  Detecting rollers in {nt}×{nc} timestack...")

        # Preprocess the image (MATLAB line 292)
        B = self.preprocessor.preprocess_timestack(I)

        # Compute std of preprocessed timestack along time (MATLAB line 295)
        # Use middle 90% of time to avoid edge effects
        t_start = round(nt / 20)
        t_end = round(19 * nt / 20)

        Breakstd_raw = np.nanstd(B[t_start:t_end, :], axis=0)  # std along time for each spatial position

        # Apply spatial smoothing (MATLAB line 295)
        Breakstd = filter_mean(filter_mean(Breakstd_raw, round(nc / 50)), round(nc / 30))

        # Compute threshold (MATLAB line 298)
        thresh = 2.0 * np.nanmax(Breakstd) / 3.0

        # Normalize Breakstd (MATLAB lines 301-305)
        Breakstd0 = (Breakstd - np.nanmin(Breakstd)) / (np.nanmax(Breakstd) - np.nanmin(Breakstd))

        # Additional breaking indicators from raw intensities
        Breakmean1 = filter_mean(np.nanmax(I, axis=0), round(nc / 50)) - filter_mean(np.nanmean(I, axis=0), round(nc / 50))
        Breakmean1 = (Breakmean1 - np.nanmin(Breakmean1)) / np.nanmax(Breakmean1 - np.nanmin(Breakmean1) + 1e-10)

        # Weight std profile (MATLAB lines 308-311)
        Breakstd_diff = -np.diff(np.concatenate(([Breakstd[0]], Breakstd)))
        Breakstd_norm = (Breakstd_diff - np.nanmin(Breakstd_diff)) / (np.nanmax(Breakstd_diff) - np.nanmin(Breakstd_diff) + 1e-10)

        # Apply weighting to preprocessed timestack
        B_weighted = B * Breakstd_norm[np.newaxis, :]

        # Set edges to zero (MATLAB line 311)
        B_weighted[0:5, :] = 0
        B_weighted[-5:, :] = 0

        # Create binary mask for breaking regions (MATLAB lines 313-316)
        B2 = np.zeros_like(B_weighted)
        B2[B_weighted > thresh] = 1
        B2[0:5, :] = 0
        B2[-5:, :] = 0

        # Detect rollers in each time frame (MATLAB lines 319-337)
        Lw = np.full_like(B2, np.nan)

        for i in range(nt - 1, -1, -1):  # Loop backwards like MATLAB
            ind = np.where(B2[i, :] == 1)[0]

            if len(ind) > 2:
                # Find gaps (multiple rollers)
                fw = np.where(np.diff(ind) > 10)[0]
                nbw = len(fw) + 1

                # Get boundaries
                if len(fw) > 0:
                    indw = sorted(np.concatenate(([ind[0]], ind[fw], ind[fw + 1], [ind[-1]])))
                else:
                    indw = [ind[0], ind[-1]]

                # Compute roller lengths
                for v in range(nbw):
                    start_idx = indw[2 * v]
                    end_idx = indw[2 * v + 1] if 2 * v + 1 < len(indw) else indw[-1]
                    mid_idx = round((start_idx + end_idx) / 2)
                    roller_length = end_idx - start_idx + 1

                    if 0 <= mid_idx < nc:
                        Lw[i, mid_idx] = roller_length

        Lw[Lw == 0] = np.nan

        # Find deflection points (breaking wave fronts) (MATLAB lines 341-348)
        Def = np.full(nt, np.nan)
        for i in range(nt):
            try:
                first_roller = np.nanargmax(Lw[i, :] > 0)
                if Lw[i, first_roller] > 0:
                    Def[i] = first_roller
            except:
                pass

        # Refine positions (MATLAB lines 352-374)
        oi = np.where(~np.isnan(Def))[0]

        if len(oi) > 1:
            # Find jumps in breaking position
            oo = np.where(np.abs(np.diff(Def[oi])) > 10)[0]
            oo = np.concatenate(([0], oo + 1))

            PosX = Def[oi[oo]]
            PosT = oi[oo]

            try:
                PosX = np.round(PosX).astype(int)

                # Threshold based on combined breaking indicators (MATLAB lines 363-369)
                vec = 0.5 * (Breakstd0 + Breakstd_norm)
                seuil = np.nanmin(vec) + 2 * (np.nanmax(vec) - np.nanmin(vec)) / 3

                # Filter by threshold
                valid_indices = PosX < len(vec)
                PosX_filtered = PosX[valid_indices]
                PosT_filtered = PosT[valid_indices]

                if len(PosX_filtered) > 0:
                    io = vec[PosX_filtered] > seuil
                    PosX = PosX_filtered[io]
                    PosT = PosT_filtered[io]

            except Exception as e:
                print(f"    Warning in position refinement: {e}")
                PosX = np.array([])
                PosT = np.array([])
        else:
            PosX = np.array([])
            PosT = np.array([])

        print(f"    Detected {len(PosX)} breaking waves")

        return {
            'PosX': PosX,
            'PosT': PosT,
            'Lw': Lw,
            'B': B,
            'Breakstd': Breakstd,
            'Breakstd0': Breakstd0,
        }


class PhotogrammetricHeightCalculator:
    """
    Calculate wave heights using photogrammetric method.

    Implements MATLAB BreakerHeight function (lines 392-455).
    """

    def __init__(self, camera_height: float = 27.24, pixel_resolution: float = 0.1):
        """
        Initialize calculator.

        Parameters
        ----------
        camera_height : float
            Camera height above mean sea level in meters
        pixel_resolution : float
            Spatial resolution in meters per pixel
        """
        self.camera_height = camera_height
        self.pixel_resolution = pixel_resolution

    def calculate_wave_heights(self, B: np.ndarray, PosT: np.ndarray, PosX: np.ndarray,
                               Lw: np.ndarray, cross_shore_positions: np.ndarray) -> Tuple[float, float]:
        """
        Calculate significant and mean wave heights using photogrammetric method.

        Matches MATLAB BreakerHeight function (lines 392-455).

        Parameters
        ----------
        B : np.ndarray
            Preprocessed timestack (time × space)
        PosT : np.ndarray
            Breaking times (temporal indices)
        PosX : np.ndarray
            Breaking positions (spatial indices)
        Lw : np.ndarray
            Roller lengths array
        cross_shore_positions : np.ndarray
            Cross-shore positions in meters

        Returns
        -------
        Tuple[float, float]
            (Hs, Hm) - Significant and mean wave heights in meters
        """
        try:
            # Wave face angle at breaking (MATLAB line 403)
            wave_face_angle = np.deg2rad(35.0)

            # Camera angles for each breaking position (MATLAB line 402)
            AngleCam = np.abs(self.camera_height / cross_shore_positions[PosX])

            # Smooth roller lengths (MATLAB line 406)
            Lwi = np.nanmax(Lw, axis=0)
            valid_Lw = ~np.isnan(Lwi) & (Lwi > 0)

            if np.sum(valid_Lw) < 2:
                return np.nan, np.nan

            # Smooth with moving average
            Lwi_smooth = uniform_filter1d(Lwi, size=min(80, len(Lwi) // 10), mode='nearest')

            L1 = []

            # Analyze each breaking wave (MATLAB lines 408-423)
            for i in range(len(PosT)):
                try:
                    t = PosT[i]
                    x = PosX[i]

                    # Extract window around breaking event (MATLAB line 410)
                    t_start = max(0, t - 25)
                    t_end = min(B.shape[0], t + 25)
                    x_start = max(0, x - 50)
                    x_end = B.shape[1]

                    # Compute spatial range for each time in window
                    window = B[t_start:t_end, x_start:x_end]

                    if window.size == 0:
                        continue

                    # Max - min along space for each time
                    spatial_range = np.nanmax(window, axis=1) - np.nanmin(window, axis=1)

                    # Smooth (MATLAB line 410)
                    vec = filter_mean(filter_mean(spatial_range, 20), 5)

                    if len(vec) < 3:
                        continue

                    # Find peak in first third (MATLAB line 411)
                    first_third_end = max(1, len(vec) // 3)
                    ngt = np.nanmax(vec[:first_third_end])
                    gft = np.nanargmax(vec[:first_third_end])

                    # Find extent of high intensity region (MATLAB lines 412-418)
                    threshold = np.nanmin(vec) + 0.5 * (np.nanmax(vec) - np.nanmin(vec))
                    ii = np.where(vec > threshold)[0]

                    if len(ii) < 2:
                        continue

                    # Find boundaries with gaps
                    gaps = np.where(np.diff(ii) > 20)[0]

                    if len(gaps) > 0:
                        id_arr = np.concatenate(([ii[0]], ii[gaps], ii[gaps + 1], [ii[-1]]))
                    else:
                        id_arr = np.array([ii[0], ii[-1]])

                    # Find boundaries around peak
                    idl = id_arr[id_arr < gft]
                    idp = id_arr[id_arr > gft]

                    if len(idl) > 0 and len(idp) > 0:
                        idl_end = idl[-1]
                        idp_start = idp[0]

                        # Horizontal extent in pixels (MATLAB line 419)
                        b = idp_start - idl_end
                        L1_val = abs(b) * self.pixel_resolution
                        L1.append(L1_val)

                except Exception as e:
                    continue

            if len(L1) == 0:
                return np.nan, np.nan

            L1 = np.array(L1)

            # Photogrammetric conversion (MATLAB lines 425-426)
            median_angle = np.nanmean(AngleCam)
            correction = L1 * np.tan(median_angle) / np.tan(wave_face_angle)
            Lf = (L1 - correction) * np.tan(median_angle)

            # Filter valid heights (MATLAB line 428)
            ind = np.where((Lf > 0) & ~np.isnan(Lf))[0]

            if len(ind) == 0:
                return np.nan, np.nan

            # Sort and remove outliers (MATLAB line 429-430)
            Lord = np.sort(Lf[ind])[::-1]  # Descending

            # Keep middle 80%
            start_idx = max(0, round(len(Lord) / 10))
            end_idx = min(len(Lord), round(9 * len(Lord) / 10))
            Lord_filtered = Lord[start_idx:end_idx]

            if len(Lord_filtered) == 0:
                return np.nan, np.nan

            # Significant wave height: median of top 1/3 (MATLAB line 432)
            top_third_count = max(1, round(len(Lord_filtered) / 3))
            hs = np.nanmedian(Lord_filtered[:top_third_count])

            # Mean wave height: median of all (MATLAB line 433)
            hm = np.nanmedian(Lord_filtered)

            return hs, hm

        except Exception as e:
            print(f"    BreakerHeight calculation failed: {e}")
            return np.nan, np.nan


def filter_mean(signal: np.ndarray, width: int) -> np.ndarray:
    """
    Filter data using running mean filter.

    Matches MATLAB FilterMean function (lines 1-58).

    Parameters
    ----------
    signal : np.ndarray
        Input signal (1D array)
    width : int
        Filter width

    Returns
    -------
    np.ndarray
        Filtered signal
    """
    if signal.ndim > 1:
        # Flatten if needed
        signal = signal.flatten()

    la = len(signal)

    if la < (2 * width):
        print(f"Warning: signal too small for filtering (length={la}, width={width})")
        return signal

    # Pad signal (MATLAB lines 24-26)
    fsig2 = np.zeros(la + 2 * width)
    fsig2[0:width] = np.mean(signal[0:width])
    fsig2[width:width + la] = signal
    fsig2[width + la:] = np.mean(signal[la - width:la])

    # Apply running mean (MATLAB lines 30-33)
    fsig = np.zeros(la)
    for kk in range(width, width + la):
        k = kk - width
        fsig[k] = np.mean(fsig2[kk - width:kk + width + 1])

    # Replace NaNs with original values (MATLAB line 57)
    nan_mask = np.isnan(fsig)
    fsig[nan_mask] = signal[nan_mask]

    return fsig
