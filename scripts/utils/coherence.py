#!/usr/bin/env python3
"""Phase-coherent correlation utilities for TEP-GNSS-MGEX.

These functions are ported from the proven TEP-GNSS-RINEX (Kathmandu) pipeline
and are identical to the CODE-longspan methodology where noted.
"""

import numpy as np
from scipy.signal import csd, welch, detrend


def compute_coherence_phase(v1, v2, fs, f1, f2):
    """Compute coherence and magnitude-weighted phase for a single pair.

    Parameters
    ----------
    v1, v2 : array-like
        Time series (same length, same sampling).
    fs : float
        Sampling frequency in Hz.
    f1, f2 : float
        Frequency band limits in Hz.

    Returns
    -------
    mean_coh : float
        Mean normalized coherence in band (0–1).
    weighted_phase : float
        Magnitude-weighted circular mean phase (rad).
    phase_alignment : float
        cos(weighted_phase), the TEP signature metric.
    """
    n_points = len(v1)
    if n_points < 64:
        return np.nan, np.nan, np.nan

    v1_d = detrend(v1, type='linear')
    v2_d = detrend(v2, type='linear')

    from scripts.utils.config import NPERSEG_MAX
    nperseg = min(NPERSEG_MAX, n_points // 2)
    if nperseg < 32:
        return np.nan, np.nan, np.nan

    try:
        f, Pxy = csd(v1_d, v2_d, fs=fs, nperseg=nperseg, detrend='constant')
        _, Pxx = welch(v1_d, fs=fs, nperseg=nperseg, detrend='constant')
        _, Pyy = welch(v2_d, fs=fs, nperseg=nperseg, detrend='constant')

        mask = (f >= f1) & (f <= f2)
        if not np.any(mask):
            return np.nan, np.nan, np.nan

        Pxy_band = Pxy[mask]
        Pxx_band = Pxx[mask]
        Pyy_band = Pyy[mask]

        denom = Pxx_band * Pyy_band
        valid_denom = denom > 0
        if not np.any(valid_denom):
            return np.nan, np.nan, np.nan

        coh_squared = np.abs(Pxy_band[valid_denom])**2 / denom[valid_denom]
        mean_coh = np.mean(np.sqrt(coh_squared))

        raw_magnitudes = np.abs(Pxy_band)
        phases = np.angle(Pxy_band)

        if np.sum(raw_magnitudes) > 0:
            complex_phases = np.exp(1j * phases)
            weighted_complex = np.average(complex_phases, weights=raw_magnitudes)
            weighted_phase = np.angle(weighted_complex)
            phase_alignment = np.cos(weighted_phase)
        else:
            weighted_phase = np.nan
            phase_alignment = np.nan

        return float(mean_coh), float(weighted_phase), float(phase_alignment)
    except Exception:
        return np.nan, np.nan, np.nan


def exponential_decay(r, A, lam, C0):
    return A * np.exp(-r / lam) + C0


def fit_exponential_model(distances, coherences, min_distance_km, max_distance_km, n_bins, min_bin_count):
    """Fit C(r) = A exp(−r/λ) + C₀ to distance-binned coherence data.

    Returns a dict with keys: amplitude, correlation_length_km, offset,
    r_squared, success, or None if fitting fails.
    """
    if len(distances) < 100:
        return None

    from scipy.optimize import curve_fit

    bin_edges = np.logspace(np.log10(min_distance_km), np.log10(max_distance_km), n_bins + 1)
    bin_centers, bin_means, bin_counts, bin_stds = [], [], [], []

    for i in range(n_bins):
        mask = (distances >= bin_edges[i]) & (distances < bin_edges[i + 1])
        n_mask = np.sum(mask)
        if n_mask >= min_bin_count:
            bin_centers.append((bin_edges[i] + bin_edges[i + 1]) / 2)
            bin_means.append(np.mean(coherences[mask]))
            bin_counts.append(n_mask)
            bin_stds.append(float(np.std(coherences[mask], ddof=1)))

    if len(bin_centers) < 5:
        return None

    x = np.array(bin_centers)
    y = np.array(bin_means)
    counts = np.array(bin_counts)
    # Weight by 1/SEM where SEM = sigma_bin / sqrt(count).
    # We do not know the per-bin population sigma, so we approximate it by the
    # sample std within each bin; for large N this is well-constrained.
    sigmas = np.array(bin_stds)
    sem = sigmas / np.sqrt(counts)
    # Guard against zero SEM (all identical values in a bin)
    sem = np.where(sem > 0, sem, np.mean(sem[sem > 0]) if np.any(sem > 0) else 1.0)
    w = 1.0 / sem

    try:
        # Dynamic bounds: lambda cannot sensibly exceed ~2x the max distance sampled
        lambda_guess = (min_distance_km + max_distance_km) / 2.0
        lambda_upper = max_distance_km * 2.0
        popt, pcov = curve_fit(
            exponential_decay, x, y,
            p0=[0.5, lambda_guess, 0],
            sigma=1.0 / w,
            bounds=([0, min_distance_km, -1], [2, lambda_upper, 1]),
            maxfev=5000
        )

        predicted = exponential_decay(x, *popt)
        ss_res = np.sum((y - predicted)**2)
        ss_tot = np.sum((y - np.mean(y))**2)
        r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0

        # Parameter uncertainties from covariance matrix
        perr = np.sqrt(np.diag(pcov))

        return {
            'amplitude': float(popt[0]),
            'correlation_length_km': float(popt[1]),
            'correlation_length_err_km': float(perr[1]),
            'offset': float(popt[2]),
            'r_squared': float(r2),
            'bin_centers': [float(v) for v in bin_centers],
            'bin_means': [float(v) for v in bin_means],
            'bin_counts': [int(v) for v in bin_counts],
            'success': True
        }
    except Exception:
        return None
