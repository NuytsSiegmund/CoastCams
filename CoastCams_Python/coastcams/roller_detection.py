"""
Roller detection and wave breaking event identification.

This module implements the RollerPropertiesTaller function from MATLAB,
which detects individual wave breaking events in timestack images.
"""

import numpy as np
from scipy import signal
from scipy.ndimage import uniform_filter1d
from typing import Tuple, Optional


def filter_mean(a: np.ndarray, b: int) -> np.ndarray:
    """
    Filter data using running mean filter.

    Matches MATLAB FilterMean.m function.

    Parameters
    ----------
    a : np.ndarray
        Input vector or 1D array
    b : int
        Filter width (half-window size)

    Returns
    -------
    np.ndarray
        Filtered signal
    """
    # Handle 1D input
    a = np.atleast_1d(a).flatten()
    la = len(a)

    if la < (2 * b):
        print('Warning: vector too small for filtering')
        return a

    # Pad signal with edge values
    fsig2 = np.zeros(la + 2 * b)
    fsig2[0:b] = np.mean(a[0:b])
    fsig2[b:b + la] = a
    fsig2[b + la:] = np.mean(a[la - b:la])

    # Apply running mean filter
    fsig = np.zeros(la)
    for k in range(la):
        kk = k + b
        fsig[k] = np.mean(fsig2[kk - b:kk + b + 1])

    # Replace NaN values with original
    nan_mask = np.isnan(fsig)
    fsig[nan_mask] = a[nan_mask]

    return fsig


def radon_separation(M: np.ndarray) -> np.ndarray:
    """
    Separate incident wave component using Radon transform.

    Matches MATLAB RadonSeparationmodif and FiltreRadon functions.

    Parameters
    ----------
    M : np.ndarray
        Timestack (time × space)

    Returns
    -------
    np.ndarray
        Incident component (filtered between angles 1-89 degrees)
    """
    from skimage.transform import radon, iradon

    nt, nx = M.shape

    if nx > nt:
        # If dimensions are swapped, return as-is
        return M

    # Pre-treatment: demean each spatial column (MATLAB line 2000)
    Sin = M.copy()
    for i in range(nx):
        Sin[:, i] = Sin[:, i] - np.nanmean(Sin[:, i])

    try:
        # Apply Radon filtering to separate incident component
        # MATLAB: [Sin]=FiltreRadon(M(1:nt,1:nx),1,89)

        # Radon transform with angles 0-179 (MATLAB line 2023)
        theta = np.arange(180)
        R = radon(Sin, theta=theta, circle=False)

        # Filter to keep only angles 1-89 (incident waves, MATLAB line 2026)
        Lm = 1
        Lx = 89
        R_filtered = R[:, Lm:Lx+1]
        theta_filtered = theta[Lm:Lx+1]

        # Inverse Radon transform (MATLAB line 2026)
        I = iradon(R_filtered, theta=theta_filtered, filter_name='hann', output_size=nt)

        # Extract center region to match original size (MATLAB lines 2029-2040)
        mx = I.shape[0]
        mn = min(nt, nx)

        if mx > mn:
            ind_start = round(mx/2 - mn/2)
            ind_end = ind_start + mn
            I_cropped = I[:, ind_start:ind_end]
        else:
            I_cropped = I

        # Handle size mismatch
        S2 = np.full_like(Sin, np.nan)

        if I_cropped.shape[1] <= nx:
            # Center the result
            c1, c2 = nt, nx
            if c2 > c1:
                start_col = (c2 - I_cropped.shape[1]) // 2
                S2[:min(nt, I_cropped.shape[0]), start_col:start_col + I_cropped.shape[1]] = I_cropped[:min(nt, I_cropped.shape[0]), :]
            else:
                S2[:min(nt, I_cropped.shape[0]), :min(nx, I_cropped.shape[1])] = I_cropped[:min(nt, I_cropped.shape[0]), :min(nx, I_cropped.shape[1])]

        # Multiply by 0.5 (MATLAB line 2044)
        S2 = S2 * 0.5

        # Replace NaN with original demeaned values
        S2[np.isnan(S2)] = Sin[np.isnan(S2)]

        return S2

    except Exception as e:
        # If Radon transform fails, return demeaned version
        print(f"Warning: Radon transform failed ({str(e)}), using demeaned data")
        return Sin


def image_preprocessing(A: np.ndarray, dt: float) -> np.ndarray:
    """
    Preprocess timestack image with bandpass filtering.

    Matches MATLAB ImagePreProcessing_20090121Taller function.

    Parameters
    ----------
    A : np.ndarray
        Input timestack (time × space)
    dt : float
        Temporal resolution (seconds)

    Returns
    -------
    np.ndarray
        Preprocessed timestack
    """
    nt, nc = A.shape[:2]

    # Spatial resolution for processing
    resc = max(1, round(nc / 100))

    B = np.zeros_like(A[:, :, 0] if A.ndim == 3 else A)

    # Design FIR filters to match MATLAB exactly
    # MATLAB: ord = 1000, fir1(ord, Val, 'low'/'high')
    fr = 1.0 / dt

    # Low-pass filter (cutoff = 1.5s)
    Tcoupure_low = 1.5
    Val_low = (1.0 / Tcoupure_low) * 2 * (1.0 / fr)
    ord = 1000

    # Ensure normalized frequency is valid
    if Val_low >= 1.0:
        Val_low = 0.99

    # Design FIR filters matching MATLAB fir1()
    try:
        fil_low = signal.firwin(ord + 1, Val_low, window='hamming')
    except:
        # Fallback if filter design fails
        fil_low = signal.firwin(min(ord, nt-1), Val_low, window='hamming')

    # High-pass filter (cutoff = 20s)
    Tcoupure_high = 20
    Val_high = (1.0 / Tcoupure_high) * 2 * (1.0 / fr)

    if Val_high >= 1.0:
        Val_high = 0.99

    try:
        fil_high = signal.firwin(ord + 1, Val_high, pass_zero=False, window='hamming')
    except:
        fil_high = signal.firwin(min(ord, nt-1), Val_high, pass_zero=False, window='hamming')

    # Process every resc-th spatial column
    for ic in range(0, nc, resc):
        # Extract time series at this spatial position
        if A.ndim == 3:
            ts = A[:, ic, 0]
        else:
            ts = A[:, ic]

        # Low-pass filter using convolution (matching MATLAB)
        kk1 = np.convolve(fil_low, ts, mode='full')
        # Trim edges like MATLAB: kk1((max(ord)/2)+1 : length(kk1)-(max(ord)/2))
        edge_trim = ord // 2
        S = signal.detrend(kk1[edge_trim:len(kk1)-edge_trim])

        # High-pass filter using convolution
        kk2 = np.convolve(fil_high, S, mode='full')
        # Trim edges
        y = signal.detrend(kk2[edge_trim:len(kk2)-edge_trim])

        B[:, ic] = y

    # Interpolate across spatial dimension
    ii = np.where((np.abs(np.mean(B, axis=0)) > 0) & (~np.isnan(np.mean(B, axis=0))))[0]

    if len(ii) > 1:
        for irt in range(nt):
            try:
                B[irt, :] = np.interp(np.arange(nc), ii, B[irt, ii])
            except:
                pass

    # Replace NaN with 0
    B[np.isnan(B)] = 0

    # Apply Radon separation to extract incident wave component (MATLAB line 202)
    B = radon_separation(B)

    return B


def roller_properties_taller(S: np.ndarray, dt: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray,
                                                                  np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Detect wave breaking events and roller properties.

    Matches MATLAB RollerPropertiesTaller function.

    Parameters
    ----------
    S : np.ndarray
        Input timestack (time × space), single channel
    dt : float
        Temporal resolution (seconds)

    Returns
    -------
    PosX : np.ndarray
        Spatial positions of breaking waves (pixel indices)
    PosT : np.ndarray
        Temporal positions of breaking waves (frame indices)
    Lw : np.ndarray
        Roller lengths array (time × space)
    B : np.ndarray
        Preprocessed timestack
    Breakstd : np.ndarray
        Standard deviation profile
    Breakmean1 : np.ndarray
        Mean profile 1
    Breakmean2 : np.ndarray
        Mean profile 2
    """
    # Convert to double
    I = S.astype(np.float64)

    # Get dimensions (note: MATLAB transposes, so we adjust)
    if I.ndim == 3:
        nt, nc, nv = I.shape
        I = I[:, :, 0]  # Use first channel
    else:
        nt, nc = I.shape

    # Perform preprocessing
    B = image_preprocessing(I, dt)

    # Compute standard deviation profile
    # Use middle 90% of time to avoid edge effects
    t_start = max(0, round(nt / 20))
    t_end = min(nt, round(19 * nt / 20))
    Breakstd_raw = np.nanstd(B[t_start:t_end, :], axis=0)

    # Smooth the standard deviation
    Breakstd = filter_mean(filter_mean(Breakstd_raw, round(nc / 50)), round(nc / 30))

    # Compute threshold
    thresh = 2.0 * np.nanmax(Breakstd) / 3.0

    # Normalize standard deviation
    Breakstd0 = (Breakstd - np.nanmin(Breakstd)) / np.nanmax(Breakstd - np.nanmin(Breakstd))

    # Compute additional breaking metrics
    Breakmean1 = filter_mean(np.nanmax(I, axis=0), round(nc / 50)) - filter_mean(np.nanmean(I, axis=0), round(nc / 50))
    Breakmean1 = (Breakmean1 - np.nanmin(Breakmean1)) / np.nanmax(Breakmean1 - np.nanmin(Breakmean1))

    Breakmean2 = filter_mean(np.nanmean(I, axis=0), round(nc / 50))
    Breakmean2 = (Breakmean2 - np.nanmin(Breakmean2)) / np.nanmax(Breakmean2 - np.nanmin(Breakmean2))

    # Compute spatial derivative and normalize
    Breakstd_diff = -np.diff(np.concatenate([[Breakstd[0]], Breakstd]))
    Breakstd_norm = (Breakstd_diff - np.nanmin(Breakstd_diff)) / np.nanmax(Breakstd_diff - np.nanmin(Breakstd_diff))

    # Weight the preprocessed image by standard deviation
    B = B * Breakstd_norm[np.newaxis, :]

    # Zero out edges
    B[0:5, :] = 0
    B[nt - 5:nt, :] = 0

    # Create binary mask for rollers
    B2 = np.zeros_like(B)
    B2[B > thresh] = 1
    B2[0:5, :] = 0
    B2[nt - 5:nt, :] = 0

    # Detect rollers in each frame
    Lw = np.full((nt, nc), np.nan)

    for i in range(nt - 1, -1, -1):  # Iterate backwards like MATLAB
        # Find roller pixels in this frame
        ind = np.where(B2[i, :] == 1)[0]

        if len(ind) > 2:
            # Detect multiple rollers (gaps > 10 pixels)
            gaps = np.where(np.diff(ind) > 10)[0]
            nbw = len(gaps) + 1

            # Compute boundaries of each roller
            indw = np.sort(np.concatenate([[ind[0]],
                                          ind[gaps] if len(gaps) > 0 else [],
                                          ind[gaps + 1] if len(gaps) > 0 else [],
                                          [ind[-1]]]))

            # Calculate length of each roller
            for v in range(nbw):
                start_idx = indw[2 * v]
                end_idx = indw[2 * v + 1]
                mid_pos = round(np.mean([start_idx, end_idx]))
                roller_length = end_idx - start_idx + 1

                if 0 <= mid_pos < nc:
                    Lw[i, mid_pos] = roller_length

    # Replace 0 with NaN
    Lw[Lw == 0] = np.nan

    # Find deflection points (first roller in each frame)
    Def = np.full(nt, np.nan)
    for i in range(nt):
        valid = np.where(Lw[i, :] > 0)[0]
        if len(valid) > 0:
            Def[i] = valid[0]

    # Refine positions
    oi = np.where(~np.isnan(Def))[0]

    if len(oi) < 2:
        # No valid breaking events
        return np.array([np.nan]), np.array([np.nan]), Lw, B, Breakstd, Breakmean1, Breakmean2

    # Find jumps in deflection position
    oo_idx = np.where(np.abs(np.diff(Def[oi])) > 10)[0]
    oo = oi[np.concatenate([[0], oo_idx + 1])]

    PosX = Def[oo]
    PosT = oo

    try:
        PosX = np.round(PosX).astype(int)

        # Compute threshold for filtering
        vec = 0.5 * (Breakstd0 + Breakstd_norm)
        seuil = np.nanmin(vec) + 2 * (np.nanmax(vec) - np.nanmin(vec)) / 3.0

        # Filter by threshold
        valid_mask = np.zeros(len(PosX), dtype=bool)
        for idx in range(len(PosX)):
            if 0 <= PosX[idx] < len(vec):
                if vec[PosX[idx]] > seuil:
                    valid_mask[idx] = True

        PosX = PosX[valid_mask]
        PosT = PosT[valid_mask]

    except:
        PosX = np.array([np.nan])
        PosT = np.array([np.nan])

    # Remove overlapping events
    to_remove = np.zeros(len(PosT), dtype=bool)

    for i in range(len(PosX)):
        for k in range(len(PosX)):
            if i != k and not to_remove[k]:
                try:
                    if 0 <= PosX[i] < nc and 0 <= PosX[k] < nc:
                        # Check trajectory between events
                        ix = np.round(np.linspace(PosX[i], PosX[k], 20)).astype(int)
                        it = np.round(np.linspace(PosT[i], PosT[k], 20)).astype(int)

                        # Clip to valid ranges
                        ix = np.clip(ix, 0, nc - 1)
                        it = np.clip(it, 0, nt - 1)

                        # Get diagonal values
                        diag_vals = B[it, ix]

                        # If trajectory is mostly negative and later in time, remove it
                        if np.sum(diag_vals < 0) < len(ix) / 5 and PosT[i] < PosT[k]:
                            to_remove[k] = True
                except:
                    pass

    # Remove marked events
    PosX = PosX[~to_remove]
    PosT = PosT[~to_remove]

    return PosX, PosT, Lw, B, Breakstd, Breakmean1, Breakmean2
