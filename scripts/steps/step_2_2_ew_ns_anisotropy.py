#!/usr/bin/env python3
"""Step 2.2: EW > NS anisotropy test.

Frozen prediction: east–west phase_alignment exceeds north–south.
Test: azimuth-sector matched station pairs, EW vs NS correlation-length ratio.

STATISTICAL RIGOR AND MULTIPLE-TESTING CORRECTION
-------------------------------------------------
This step performs ONE primary directional test per constellation:
  - H0: lambda_EW / lambda_NS <= 1 (null, no anisotropy)
  - H1: lambda_EW / lambda_NS > 1   (TEP prediction)

The p_value reported for the full-range ratio is the primary test.
The longitude-difference-filtered subset (lon_diff < 30deg) is an
EXPLORATORY secondary test motivated by the RINEX-era finding that
ionospheric decorrelation inverts the ratio at long baselines.
Its p-value is reported but NOT treated as an independent replication;
instead it provides a physical-consistency check.

No Bonferroni/FDR correction is applied to the primary full-range test
because it is a single directional hypothesis. The four constellations
share identical pair tables (combined-clock architecture), so they are
NOT independent tests; replication is assessed via bootstrap CI, not via
a family-wise correction across constellations.

NOTE: This test was NOT formally pre-registered (no pre-registration
document exists for TEP-GNSS-MGEX). The directional hypothesis was
motivated by the TEP framework, but the specific analysis choices
(phase_alignment metric, Welch CSD, 10-degree CMB grid) were refined
through exploratory development.
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

import json
import numpy as np
from scipy import stats
from scripts.utils.logger import TEPLogger, set_step_logger, print_status
from scripts.utils.config import (
    OUTPUTS_DIR, LOGS_DIR, DATA_DIR, CONSTELLATIONS, MIN_ANISOTROPY_RATIO,
    MIN_DISTANCE_KM, MAX_DISTANCE_KM, N_BINS, MIN_BIN_COUNT, N_BOOTSTRAP,
    MIN_PAIRS_FIT, MIN_PAIRS_ANALYSIS, MIN_PAIRS_BOOTSTRAP, SIGNIFICANCE_ALPHA,
    SECTOR_WIDTH_DEG, LON_DIFF_THRESHOLD_DEG, BOOTSTRAP_SEED
)
from scripts.utils.coherence import fit_exponential_model
from scripts.utils.geospatial import _sector_mask

PROCESSED_DIR = DATA_DIR / "processed"

LOGS_DIR.mkdir(parents=True, exist_ok=True)
logger = TEPLogger("step_2_2", log_file_path=LOGS_DIR / "step_2_2_anisotropy.log")
set_step_logger(logger)


def _load_station_longitudes():
    """Load station coordinates (ECEF XYZ) and return {sta: lon_deg}."""
    coord_file = PROCESSED_DIR / "station_coordinates.json"
    if not coord_file.exists():
        return {}
    with open(coord_file) as f:
        coords = json.load(f)
    lon_map = {}
    for sta, (x, y, z) in coords.items():
        lon = float(np.degrees(np.arctan2(float(y), float(x))))
        lon_map[sta] = lon
    return lon_map


def _compute_lon_diff(sta1_list, sta2_list, lon_map):
    """Compute absolute longitude difference for each pair, handling wrap-around."""
    lons1 = np.array([lon_map.get(s, np.nan) for s in sta1_list])
    lons2 = np.array([lon_map.get(s, np.nan) for s in sta2_list])
    diff = np.abs((lons1 - lons2 + 180.0) % 360.0 - 180.0)
    return diff


def _bootstrap_ratio(dist, coh, az, sector_centers, sector_names,
                     ew_sectors, ns_sectors, n_boot, seed=BOOTSTRAP_SEED, max_sample=None):
    """Bootstrap the λ_EW/λ_NS ratio by resampling pairs with replacement.

    For large datasets each resample draws at most ``max_sample`` pairs (with
    replacement). With ~40 log-spaced distance bins this still leaves thousands
    of pairs per bin, so the exponential fit is stable while runtime stays
    bounded.

    Returns dict with ratio_mean, ratio_std, ci_low, ci_high (95%), p_value
    (one-sided P(ratio <= 1)), and n_boot_valid. Returns None if too few
    valid bootstrap samples were obtained.
    """
    rng = np.random.default_rng(seed)
    n = len(dist)
    draw = n if max_sample is None else min(n, max_sample)
    if max_sample is not None and n > max_sample:
        logger.warning(f"Bootstrap subsampling {max_sample}/{n} pairs for stability")
    ratios = []
    for _ in range(n_boot):
        idx = rng.choice(n, size=draw, replace=True)
        b_dist, b_coh, b_az = dist[idx], coh[idx], az[idx]
        b_lambdas = {}
        for sector, center in zip(sector_names, sector_centers):
            mask = _sector_mask(b_az, center, width=SECTOR_WIDTH_DEG)
            if np.sum(mask) < MIN_PAIRS_FIT:
                continue
            fit = fit_exponential_model(
                b_dist[mask], b_coh[mask],
                MIN_DISTANCE_KM, MAX_DISTANCE_KM, N_BINS, MIN_BIN_COUNT
            )
            if fit and fit.get('success'):
                b_lambdas[sector] = fit['correlation_length_km']
        b_ew = [b_lambdas[s] for s in ew_sectors if s in b_lambdas]
        b_ns = [b_lambdas[s] for s in ns_sectors if s in b_lambdas]
        if len(b_ew) >= 1 and len(b_ns) >= 1:
            ns_mean = np.mean(b_ns)
            if ns_mean > 0:
                ratios.append(float(np.mean(b_ew) / ns_mean))
    if len(ratios) < max(10, n_boot // 10):
        return None
    arr = np.array(ratios)
    return {
        "ratio_mean": float(np.mean(arr)),
        "ratio_std": float(np.std(arr)),
        "ci_low": float(np.percentile(arr, 2.5)),
        "ci_high": float(np.percentile(arr, 97.5)),
        "p_value": float(np.mean(arr <= 1.0)),
        "n_boot_valid": int(len(arr)),
    }


def _fit_lambda_weighted(dist, coh, weights, min_d, max_d, n_bins, min_count):
    """Weighted exponential fit returning λ only. Bin means are weighted by the
    per-pair weights; the effective per-bin count is the summed weight."""
    from scipy.optimize import curve_fit
    bin_edges = np.logspace(np.log10(min_d), np.log10(max_d), n_bins + 1)
    x, y, w = [], [], []
    for i in range(n_bins):
        m = (dist >= bin_edges[i]) & (dist < bin_edges[i + 1])
        if not np.any(m):
            continue
        wm = weights[m]
        sw = wm.sum()
        if sw < min_count:
            continue
        x.append((bin_edges[i] + bin_edges[i + 1]) / 2.0)
        y.append(float(np.sum(wm * coh[m]) / sw))
        w.append(np.sqrt(sw))
    if len(x) < 5:
        return None
    x, y, w = np.array(x), np.array(y), np.array(w)
    try:
        lambda_guess = (min_d + max_d) / 2.0
        lambda_upper = max_d * 2.0
        popt, _ = curve_fit(
            lambda r, A, lam, C0: A * np.exp(-r / lam) + C0, x, y,
            p0=[0.5, lambda_guess, 0], sigma=1.0 / w,
            bounds=([0, min_d, -1], [2, lambda_upper, 1]), maxfev=5000
        )
        return float(popt[1])
    except Exception:
        return None


def _bootstrap_ratio_clustered(dist, coh, az, s1_idx, s2_idx, n_stations,
                               sector_centers, sector_names, ew_sectors,
                               ns_sectors, n_boot, seed=BOOTSTRAP_SEED, max_sample=None):
    """Station-clustered (dyadic vertex) bootstrap of the λ_EW/λ_NS ratio.

    Each iteration resamples the ``n_stations`` stations with replacement and
    weights every pair by the product of its two endpoints' resample
    multiplicities. This respects the fact that pairs sharing a station are not
    independent, so the resulting CI is wider and more honest than a naive
    pair-resampling bootstrap. Returns the same dict shape as ``_bootstrap_ratio``.
    """
    rng = np.random.default_rng(seed)
    n = len(dist)
    if max_sample is not None and n > max_sample:
        logger.warning(f"Station-clustered bootstrap subsampling {max_sample}/{n} pairs for stability")
        sel = rng.choice(n, size=max_sample, replace=False)
        dist, coh, az = dist[sel], coh[sel], az[sel]
        s1_idx, s2_idx = s1_idx[sel], s2_idx[sel]
    ratios = []
    for _ in range(n_boot):
        draw = rng.integers(0, n_stations, size=n_stations)
        counts = np.bincount(draw, minlength=n_stations).astype(np.float64)
        pw = counts[s1_idx] * counts[s2_idx]
        nz = pw > 0
        if np.sum(nz) < MIN_PAIRS_BOOTSTRAP:
            continue
        d, c, a, w = dist[nz], coh[nz], az[nz], pw[nz]
        lambdas = {}
        for sector, center in zip(sector_names, sector_centers):
            mask = _sector_mask(a, center, width=SECTOR_WIDTH_DEG)
            if np.sum(mask) < MIN_PAIRS_FIT:
                continue
            lam = _fit_lambda_weighted(
                d[mask], c[mask], w[mask],
                MIN_DISTANCE_KM, MAX_DISTANCE_KM, N_BINS, MIN_BIN_COUNT
            )
            if lam is not None:
                lambdas[sector] = lam
        b_ew = [lambdas[s] for s in ew_sectors if s in lambdas]
        b_ns = [lambdas[s] for s in ns_sectors if s in lambdas]
        if len(b_ew) >= 1 and len(b_ns) >= 1:
            ns_mean = np.mean(b_ns)
            if ns_mean > 0:
                ratios.append(float(np.mean(b_ew) / ns_mean))
    if len(ratios) < max(10, n_boot // 10):
        return None
    arr = np.array(ratios)
    return {
        "method": "station_clustered",
        "ratio_mean": float(np.mean(arr)),
        "ratio_std": float(np.std(arr)),
        "ci_low": float(np.percentile(arr, 2.5)),
        "ci_high": float(np.percentile(arr, 97.5)),
        "p_value": float(np.mean(arr <= 1.0)),
        "n_boot_valid": int(len(arr)),
    }


def _azimuthal_profile(dist, coh, az, window=45, step=15):
    """Continuous λ(θ) profile: fit λ in overlapping azimuth windows."""
    profile = []
    for center in range(0, 360, step):
        mask = _sector_mask(az, center, width=window)
        if np.sum(mask) < MIN_PAIRS_FIT:
            continue
        fit = fit_exponential_model(
            dist[mask], coh[mask],
            MIN_DISTANCE_KM, MAX_DISTANCE_KM, N_BINS, MIN_BIN_COUNT
        )
        if fit and fit.get('success'):
            profile.append({
                "azimuth_deg": int(center),
                "lambda_km": float(fit['correlation_length_km']),
                "lambda_err_km": float(fit.get('correlation_length_err_km', 0.0)),
                "r_squared": float(fit['r_squared']),
                "n_pairs": int(np.sum(mask)),
            })
    return profile


def _longitude_stratified(dist, coh, az, lon_diff, sector_names, sector_centers,
                          ew_sectors, ns_sectors):
    """EW/NS ratio as a function of station-pair longitude difference.

    Demonstrates whether the anisotropy degrades as pairs span larger
    longitude separations (i.e. larger local-solar-time differences), the
    signature of ionospheric rather than geometric/TEP origin.
    """
    ld_bins = [(0, 15), (15, 30), (30, 60), (60, 120), (120, 180)]
    out = {}
    for low, high in ld_bins:
        key = f"{low}-{high}deg"
        m = (lon_diff >= low) & (lon_diff < high) & np.isfinite(lon_diff)
        if np.sum(m) < MIN_PAIRS_ANALYSIS // 2:
            out[key] = {"status": "insufficient_data", "n_pairs": int(np.sum(m))}
            continue
        d, c, a = dist[m], coh[m], az[m]
        ew_mask = _sector_mask(a, 90) | _sector_mask(a, 270)
        ns_mask = _sector_mask(a, 0) | _sector_mask(a, 180)
        ew_mean = float(np.mean(c[ew_mask])) if np.any(ew_mask) else np.nan
        ns_mean = float(np.mean(c[ns_mask])) if np.any(ns_mask) else np.nan
        ratio = float(ew_mean / ns_mean) if (np.isfinite(ew_mean) and np.isfinite(ns_mean) and ns_mean != 0) else None
        out[key] = {
            "method": "mean_phase_alignment",
            "ew_mean": ew_mean,
            "ns_mean": ns_mean,
            "ratio": ratio,
            "n_pairs": int(np.sum(m)),
        }
    return out


class Step22EWNSAnisotropy:
    def run(self):
        print_status("Step 2.2: EW/NS Anisotropy", "INFO")
        import copy
        results = {}
        _cached = None  # MGEX is a single solution; compute once, replicate across constellations
        for const_name in CONSTELLATIONS:
            if _cached is not None:
                results[const_name] = copy.deepcopy(_cached)
                continue
            pair_file = OUTPUTS_DIR / "step_2_0_mgex_pairs.json"
            if not pair_file.exists():
                results[const_name] = {"status": "no_data", "note": f"Pair file not found: {pair_file}"}
                continue

            with open(pair_file) as f:
                records = json.load(f)
            if not records:
                results[const_name] = {"status": "no_data", "note": "Empty pair records"}
                continue

            # Use phase_alignment (cos of magnitude-weighted CSD phase) — the TEP signature metric.
            # The old code incorrectly used "coherence" (mean spectral magnitude).
            dist = np.array([r["distance_km"] for r in records])
            coh = np.array([r["phase_alignment"] for r in records])
            az = np.array([r["azimuth_deg"] for r in records])

            sta1_all = np.array([r["sta1"] for r in records])
            sta2_all = np.array([r["sta2"] for r in records])

            # Filter invalid values (failed CSD computations)
            valid = np.isfinite(dist) & np.isfinite(coh) & np.isfinite(az)
            dist, coh, az = dist[valid], coh[valid], az[valid]
            sta1_arr = sta1_all[valid]
            sta2_arr = sta2_all[valid]

            if len(coh) < MIN_PAIRS_ANALYSIS:
                results[const_name] = {"status": "insufficient_data", "n_valid": int(len(coh))}
                continue

            # Integer station codes for the station-clustered (dyadic) bootstrap
            uniq_sta = np.unique(np.concatenate([sta1_arr, sta2_arr]))
            sta_code = {s: i for i, s in enumerate(uniq_sta)}
            n_stations = len(uniq_sta)
            s1_idx = np.array([sta_code[s] for s in sta1_arr])
            s2_idx = np.array([sta_code[s] for s in sta2_arr])

            # 8-sector classification (matches original TEP-GNSS methodology)
            sector_names = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW']
            sector_centers = [0, 45, 90, 135, 180, 225, 270, 315]

            sector_lambdas = {}
            sector_fits = {}
            for sector, center in zip(sector_names, sector_centers):
                mask = _sector_mask(az, center, width=SECTOR_WIDTH_DEG)
                n_sector = int(np.sum(mask))
                if n_sector < MIN_PAIRS_FIT:
                    continue
                fit = fit_exponential_model(
                    dist[mask], coh[mask],
                    MIN_DISTANCE_KM, MAX_DISTANCE_KM, N_BINS, MIN_BIN_COUNT
                )
                if fit and fit.get('success'):
                    sector_lambdas[sector] = float(fit['correlation_length_km'])
                    # Strip heavy arrays from JSON output
                    sector_fits[sector] = {
                        k: v for k, v in fit.items()
                        if k not in ('bin_centers', 'bin_means', 'bin_counts')
                    }

            if len(sector_lambdas) < 4:
                results[const_name] = {
                    "status": "insufficient_sectors",
                    "n_sectors": len(sector_lambdas),
                    "sector_lambdas": sector_lambdas
                }
                continue

            # λ_EW from E and W sectors; λ_NS from N and S sectors
            ew_sectors = ['E', 'W']
            ns_sectors = ['N', 'S']
            ew_lambdas = [sector_lambdas[s] for s in ew_sectors if s in sector_lambdas]
            ns_lambdas = [sector_lambdas[s] for s in ns_sectors if s in sector_lambdas]

            if len(ew_lambdas) < 1 or len(ns_lambdas) < 1:
                results[const_name] = {
                    "status": "insufficient_ew_ns",
                    "ew_n": len(ew_lambdas),
                    "ns_n": len(ns_lambdas),
                    "sector_lambdas": sector_lambdas
                }
                continue

            ew_lambda = float(np.mean(ew_lambdas))
            ns_lambda = float(np.mean(ns_lambdas))
            ratio = float(ew_lambda / ns_lambda) if ns_lambda > 0 else None

            # Distance-stratified anisotropy analysis
            # RINEX finding: <500 km → TEP dominates (ratio > 1)
            #               >1000 km → ionospheric decorrelation dominates (ratio < 1)
            # For short distances λ fitting is impossible (no long-baseline pairs to
            # constrain the exponential).  Use mean phase_alignment as proxy instead.
            dist_bins = {
                "short": (0, 500),
                "medium": (500, 1000),
                "long": (1000, float('inf'))
            }
            dist_stratified = {}
            for bin_name, (low, high) in dist_bins.items():
                if high == float('inf'):
                    dmask = (dist >= low)
                else:
                    dmask = (dist >= low) & (dist < high)
                if np.sum(dmask) < MIN_PAIRS_ANALYSIS // 2:
                    dist_stratified[bin_name] = {"status": "insufficient_data"}
                    continue

                # For short/medium distances: use mean phase_alignment as coherence proxy
                # For long distances: fit λ per sector when enough data spans full range
                if high != float('inf') and high <= 1000:
                    ew_mask = dmask & (_sector_mask(az, 90) | _sector_mask(az, 270))
                    ns_mask = dmask & (_sector_mask(az, 0) | _sector_mask(az, 180))
                    ew_mean = float(np.mean(coh[ew_mask])) if np.sum(ew_mask) > 0 else np.nan
                    ns_mean = float(np.mean(coh[ns_mask])) if np.sum(ns_mask) > 0 else np.nan
                    if np.isfinite(ew_mean) and np.isfinite(ns_mean) and ns_mean != 0:
                        bin_ratio = float(ew_mean / ns_mean)
                        dist_stratified[bin_name] = {
                            "method": "mean_phase_alignment",
                            "ew_mean": ew_mean,
                            "ns_mean": ns_mean,
                            "ratio": bin_ratio,
                            "n_ew": int(np.sum(ew_mask)),
                            "n_ns": int(np.sum(ns_mask))
                        }
                    else:
                        dist_stratified[bin_name] = {"status": "insufficient_data"}
                else:
                    # Long distances: fit λ per sector
                    bin_sector_lambdas = {}
                    for sector, center in zip(sector_names, sector_centers):
                        mask = dmask & _sector_mask(az, center, width=SECTOR_WIDTH_DEG)
                        if np.sum(mask) < MIN_PAIRS_FIT // 2:
                            continue
                        fit = fit_exponential_model(
                            dist[mask], coh[mask],
                            MIN_DISTANCE_KM, MAX_DISTANCE_KM, N_BINS, MIN_BIN_COUNT
                        )
                        if fit and fit.get('success'):
                            bin_sector_lambdas[sector] = float(fit['correlation_length_km'])
                    b_ew = [bin_sector_lambdas[s] for s in ew_sectors if s in bin_sector_lambdas]
                    b_ns = [bin_sector_lambdas[s] for s in ns_sectors if s in bin_sector_lambdas]
                    if len(b_ew) >= 1 and len(b_ns) >= 1:
                        bin_ratio = float(np.mean(b_ew) / np.mean(b_ns))
                        dist_stratified[bin_name] = {
                            "method": "lambda_fit",
                            "ew_lambda_km": float(np.mean(b_ew)),
                            "ns_lambda_km": float(np.mean(b_ns)),
                            "ratio": bin_ratio,
                            "n_sectors": len(bin_sector_lambdas),
                            "sector_lambdas": bin_sector_lambdas
                        }
                    else:
                        dist_stratified[bin_name] = {"status": "insufficient_sectors"}

            # --- Longitude-difference filtered analysis ---
            # RINEX finding: lon_diff < 30° isolates TEP signal by keeping
            # stations at similar local solar times, reducing ionospheric
            # decorrelation that inverts the anisotropy at long baselines.
            lon_map = _load_station_longitudes()
            lon_diff = _compute_lon_diff(sta1_arr, sta2_arr, lon_map)
            ld_mask = lon_diff < LON_DIFF_THRESHOLD_DEG
            n_ld = int(np.sum(ld_mask))
            ld_result = {"n_pairs": n_ld}
            if n_ld >= MIN_PAIRS_ANALYSIS:
                ld_dist = dist[ld_mask]
                ld_coh = coh[ld_mask]
                ld_az = az[ld_mask]
                ld_sector_lambdas = {}
                for sector, center in zip(sector_names, sector_centers):
                    mask = _sector_mask(ld_az, center, width=SECTOR_WIDTH_DEG)
                    if np.sum(mask) < MIN_PAIRS_FIT:
                        continue
                    fit = fit_exponential_model(
                        ld_dist[mask], ld_coh[mask],
                        MIN_DISTANCE_KM, MAX_DISTANCE_KM, N_BINS, MIN_BIN_COUNT
                    )
                    if fit and fit.get('success'):
                        ld_sector_lambdas[sector] = float(fit['correlation_length_km'])
                ld_ew = [ld_sector_lambdas[s] for s in ew_sectors if s in ld_sector_lambdas]
                ld_ns = [ld_sector_lambdas[s] for s in ns_sectors if s in ld_sector_lambdas]
                if len(ld_ew) >= 1 and len(ld_ns) >= 1:
                    ld_ew_lambda = float(np.mean(ld_ew))
                    ld_ns_lambda = float(np.mean(ld_ns))
                    ld_ratio = float(ld_ew_lambda / ld_ns_lambda) if ld_ns_lambda > 0 else None
                    ld_result.update({
                        "ew_lambda_km": ld_ew_lambda,
                        "ns_lambda_km": ld_ns_lambda,
                        "ratio": ld_ratio,
                        "n_sectors": len(ld_sector_lambdas),
                        "sector_lambdas": ld_sector_lambdas
                    })
                    # Bootstrap CI + significance for the headline filtered ratio.
                    # Pair-resampling bootstrap (legacy, CI likely too tight):
                    ld_boot = _bootstrap_ratio(
                        ld_dist, ld_coh, ld_az, sector_centers, sector_names,
                        ew_sectors, ns_sectors, N_BOOTSTRAP
                    )
                    if ld_boot is not None:
                        ld_result["bootstrap"] = ld_boot
                    # Station-clustered bootstrap (honest CI accounting for shared stations):
                    ld_boot_c = _bootstrap_ratio_clustered(
                        ld_dist, ld_coh, ld_az, s1_idx[ld_mask], s2_idx[ld_mask],
                        n_stations, sector_centers, sector_names,
                        ew_sectors, ns_sectors, N_BOOTSTRAP
                    )
                    if ld_boot_c is not None:
                        ld_result["bootstrap_station_clustered"] = ld_boot_c
                    boot_str = (
                        f", 95%CI=[{ld_boot_c['ci_low']:.2f},{ld_boot_c['ci_high']:.2f}], "
                        f"p={ld_boot_c['p_value']:.4f} (clustered)" if ld_boot_c else ""
                    )
                    print_status(
                        f"    {const_name} (lon_diff<30°): λ_EW={ld_ew_lambda:.0f}km, "
                        f"λ_NS={ld_ns_lambda:.0f}km, ratio={ld_ratio:.3f}{boot_str}",
                        "INFO"
                    )
                else:
                    ld_result["status"] = "insufficient_ew_ns"
            else:
                ld_result["status"] = "insufficient_data"

            # Bootstrap significance for the full-range ratio: resample pairs and refit
            main_boot = _bootstrap_ratio(
                dist, coh, az, sector_centers, sector_names,
                ew_sectors, ns_sectors, N_BOOTSTRAP
            )
            main_boot_c = _bootstrap_ratio_clustered(
                dist, coh, az, s1_idx, s2_idx, n_stations,
                sector_centers, sector_names, ew_sectors, ns_sectors, N_BOOTSTRAP
            )
            if main_boot_c is not None:
                p_value = main_boot_c["p_value"]
                ratio_err = main_boot_c["ratio_std"]
            elif main_boot is not None:
                p_value = main_boot["p_value"]
                ratio_err = main_boot["ratio_std"]
            else:
                p_value = 1.0
                ratio_err = None

            # Continuous azimuthal λ(θ) profile and longitude-stratified ratio
            az_profile = _azimuthal_profile(dist, coh, az)
            lon_strat = _longitude_stratified(
                dist, coh, az, lon_diff, sector_names, sector_centers,
                ew_sectors, ns_sectors
            )

            effect_size_met = ratio is not None and ratio >= MIN_ANISOTROPY_RATIO
            significance_met = p_value < SIGNIFICANCE_ALPHA

            results[const_name] = {
                "ew_lambda_km": ew_lambda,
                "ns_lambda_km": ns_lambda,
                "ew_ns_ratio": ratio,
                "ew_ns_ratio_err": ratio_err,
                "p_value": p_value,
                "sector_lambdas": sector_lambdas,
                "sector_fits": sector_fits,
                "distance_stratified": dist_stratified,
                "lon_diff_filtered": ld_result,
                "longitude_stratified": lon_strat,
                "azimuthal_profile": az_profile,
                "bootstrap_full_range": main_boot,
                "bootstrap_full_range_clustered": main_boot_c,
                "n_pairs": int(len(coh)),
                "n_stations": int(n_stations),
                "min_ratio_threshold": MIN_ANISOTROPY_RATIO,
                "effect_size_met": bool(effect_size_met),
                "significance_met": bool(significance_met),
                "status": "success",
                "replicated": bool(effect_size_met and significance_met)
            }
            _cached = results[const_name]
            rep_str = "REPLICATED" if (effect_size_met and significance_met) else "NOT_REPLICATED"
            print_status(
                f"    {const_name}: λ_EW={ew_lambda:.0f}km, λ_NS={ns_lambda:.0f}km, "
                f"ratio={ratio:.3f}, p={p_value:.4f} [{rep_str}]",
                "INFO"
            )
            for bin_name, vals in dist_stratified.items():
                if "ratio" in vals:
                    if vals.get("method") == "mean_phase_alignment":
                        print_status(
                            f"      {bin_name}: ew_mean={vals['ew_mean']:.4f}, "
                            f"ns_mean={vals['ns_mean']:.4f}, ratio={vals['ratio']:.3f}",
                            "INFO"
                        )
                    else:
                        print_status(
                            f"      {bin_name}: λ_EW={vals['ew_lambda_km']:.0f}km, "
                            f"λ_NS={vals['ns_lambda_km']:.0f}km, ratio={vals['ratio']:.3f}",
                            "INFO"
                        )

        out_file = OUTPUTS_DIR / "step_2_2_ew_ns_anisotropy.json"
        with open(out_file, 'w') as f:
            json.dump(results, f, indent=2)
        print_status(f"Anisotropy results written to {out_file}", "INFO")
        return {"status": "success", "results": str(out_file)}
