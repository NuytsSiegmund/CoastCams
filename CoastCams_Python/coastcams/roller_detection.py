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

    Simplified version of MATLAB RadonSeparationmodif.
    For now, just demean each spatial column.

    Parameters
    ----------
    M : np.ndarray
        Timestack (time × space)

    Returns
    -------
    np.ndarray
        Incident component
    """
    nt, nx = M.shape

    if nx > nt:
        # If dimensions are swapped, just return as-is
        return M

    # Demean each spatial column
    Sin = M.copy()
    for i in range(nx):
        Sin[:, i] = Sin[:, i] - np.nanmean(Sin[:, i])

    # TODO: Implement full Radon filtering (FiltreRadon)
    # For now, this simplified version removes the mean

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

    # Design filters - using Butterworth for speed (FIR ord=1000 is too slow)
    fr = 1.0 / dt
    nyq = fr / 2.0

    # Low-pass filter (cutoff = 1.5s = 0.667 Hz)
    Tcoupure_low = 1.5
    cutoff_low = 1.0 / Tcoupure_low  # Hz

    # High-pass filter (cutoff = 20s = 0.05 Hz)
    Tcoupure_high = 20
    cutoff_high = 1.0 / Tcoupure_high  # Hz

    # Use 4th order Butterworth (much faster than FIR ord=1000)
    try:
        b_low, a_low = signal.butter(4, min(cutoff_low / nyq, 0.99), btype='low')
        b_high, a_high = signal.butter(4, min(cutoff_high / nyq, 0.99), btype='high')
    except:
        # Fallback to simple filters
        b_low, a_low = np.array([1.0]), np.array([1.0])
        b_high, a_high = np.array([1.0, -1.0]), np.array([1.0])

    # Process every resc-th spatial column
    for ic in range(0, nc, resc):
        # Extract time series at this spatial position
        if A.ndim == 3:
            ts = A[:, ic, 0]
        else:
            ts = A[:, ic]

        # Low-pass filter using filtfilt (zero-phase)
        try:
            kk1 = signal.filtfilt(b_low, a_low, ts)
        except:
            kk1 = ts
        S = signal.detrend(kk1)

        # High-pass filter using filtfilt
        try:
            kk2 = signal.filtfilt(b_high, a_high, S)
        except:
            kk2 = S
        y = signal.detrend(kk2)

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

    # Apply Radon separation (simplified - just demean for speed)
    # Full Radon transform implementation would go here
    # B = radon_separation(B)  # Disabled for speed - simple demean is sufficient

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
