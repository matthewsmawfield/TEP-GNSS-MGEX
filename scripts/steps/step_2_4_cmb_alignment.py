#!/usr/bin/env python3
"""Step 2.4: CMB-frame alignment via orbital-velocity geometry.

Frozen prediction: the EW/NS phase-alignment anisotropy correlates with
Earth's orbital velocity projected onto the CMB dipole direction.

Methodology (from TEP-GNSS-II Cairo, revised):
    1. Aggregate pair phase_alignment values by day and by month into
       EW−NS mean phase-alignment differences (the difference metric
       used in step_2_3, not the unstable ew/ns ratio).
    2. For each day/month, compute Earth orbital velocity vector in
       equatorial coordinates.
    3. For a proposed background direction (RA, Dec), compute predicted
       anisotropy = dot(v_orb, n_background).
    4. Correlate predicted vs observed EW−NS differences (Pearson r).
    5. Grid-search background directions; evaluate alignment using the
       direction (max +r or max −r) closest to the CMB dipole.
    6. Permutation null test with full grid search per permutation.

The difference metric pa_diff = ew_mean − ns_mean is used instead of
the ratio ew/ns because:
  (a) phase_alignment = cos(weighted_phase) ranges from −1 to +1, so
      the ratio diverges when the NS mean is near zero;
  (b) the ns > 0 filter required by the ratio metric silently drops
      ~36% of days, introducing a non-linear selection bias;
  (c) step_2_3 (orbital coupling) uses pa_diff as its supplementary,
      model-free metric, so this scan is a directional characterization
      of the same observable rather than an independent confirmation.

Monthly aggregation uses complete chronological calendar months, preserving
the 2025--2026 temporal sequence rather than pooling months from different
years.  It reduces daily noise but cannot establish a preferred sky direction
from this short (16-month) interval. Daily results are retained as a
supplementary diagnostic.
"""

import sys
from pathlib import Path
from datetime import datetime, timedelta

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

import json
import numpy as np
from scipy import stats
from scripts.utils.logger import TEPLogger, set_step_logger, print_status
from scripts.utils.config import (
    OUTPUTS_DIR, LOGS_DIR, CONSTELLATIONS,
    CMB_DIPOLE_RA_DEG, CMB_DIPOLE_DEC_DEG,
    LOOK_ELSEWHERE_FACTOR, CMB_ALIGNMENT_THRESHOLD_DEG, SIGNIFICANCE_ALPHA,
    CMB_GRID_RA_N, CMB_GRID_DEC_N, CMB_PERMUTATION_N, CMB_SECTOR_WIDTHS,
    MIN_DAYS_CMB, MIN_DAYS_PER_MONTH, BOOTSTRAP_SEED, SECTOR_WIDTH_DEG,
    SENSITIVITY_MIN_DISTANCES
)
from scripts.utils.geospatial import _sector_mask, _earth_orbital_velocity_equatorial

LOGS_DIR.mkdir(parents=True, exist_ok=True)
logger = TEPLogger("step_2_4", log_file_path=LOGS_DIR / "step_2_4_cmb_alignment.log")
set_step_logger(logger)


def angular_separation(ra1, dec1, ra2, dec2):
    """Great-circle angular separation in degrees (spherical law of cosines)."""
    ra1, dec1, ra2, dec2 = map(np.radians, [ra1, dec1, ra2, dec2])
    cos_sep = np.sin(dec1) * np.sin(dec2) + np.cos(dec1) * np.cos(dec2) * np.cos(ra1 - ra2)
    return float(np.degrees(np.arccos(np.clip(cos_sep, -1, 1))))


def _direction_unit_vector(ra_deg, dec_deg):
    """Convert RA/Dec to unit vector in equatorial coordinates."""
    ra, dec = np.radians(ra_deg), np.radians(dec_deg)
    return np.array([
        np.cos(dec) * np.cos(ra),
        np.cos(dec) * np.sin(ra),
        np.sin(dec)
    ])


def _month_from_day_key(day_key):
    """Convert YYYYDOY to a chronological YYYY-MM calendar-month key."""
    try:
        year = int(day_key[:4])
        doy = int(day_key[4:])
        dt = datetime(year, 1, 1) + timedelta(days=doy - 1)
        return dt.strftime("%Y-%m")
    except Exception:
        return None


def _ew_ns_difference(records, sector_width=SECTOR_WIDTH_DEG):
    """Aggregate pair records by day into EW−NS phase_alignment differences.

    Uses the same 8-sector definitions as step_2_2 and step_2_3:
    EW = E (center 90°) + W (center 270°)
    NS = N (center 0°) + S (center 180°)

    The difference metric pa_diff = ew_mean − ns_mean requires no sign
    filter on the NS mean, unlike the ratio ew/ns which diverges when
    ns is near zero.

    Returns: dict {day_key: {"ew": float, "ns": float, "diff": float}}
    """
    day_data = {}
    for r in records:
        day = r.get("day", "")
        if not day:
            continue
        if day not in day_data:
            day_data[day] = {"ew": [], "ns": []}
        az = r.get("azimuth_deg", 0)
        # Use phase_alignment (TEP signature metric)
        coh = r.get("phase_alignment", np.nan)
        if not np.isfinite(coh):
            continue
        # Match step_2_2 sector definitions: EW = E + W, NS = N + S
        if _sector_mask(az, 90, width=sector_width) or _sector_mask(az, 270, width=sector_width):
            day_data[day]["ew"].append(coh)
        elif _sector_mask(az, 0, width=sector_width) or _sector_mask(az, 180, width=sector_width):
            day_data[day]["ns"].append(coh)

    diffs = {}
    for day, vals in day_data.items():
        ew = np.mean(vals["ew"]) if vals["ew"] else np.nan
        ns = np.mean(vals["ns"]) if vals["ns"] else np.nan
        if np.isfinite(ew) and np.isfinite(ns):
            diffs[day] = {"ew": float(ew), "ns": float(ns), "diff": float(ew - ns)}
    return diffs


def _monthly_from_daily(day_diffs):
    """Aggregate daily EW−NS differences into monthly means.

    Returns: dict {YYYY-MM: {"diff": float, "n_days": int}}.  Partial
    months are excluded so all observations have comparable support.
    """
    month_vals = {}
    for day, vals in day_diffs.items():
        m = _month_from_day_key(day)
        if m is None:
            continue
        if m not in month_vals:
            month_vals[m] = []
        month_vals[m].append(vals["diff"])

    monthly = {}
    for m, vals in month_vals.items():
        if len(vals) >= MIN_DAYS_PER_MONTH:
            monthly[m] = {"diff": float(np.mean(vals)), "n_days": len(vals)}
    return monthly


def _predicted_anisotropy(obs_values, day_keys, ra_bg, dec_bg):
    """Compute Pearson r between predicted (velocity·background) and observed.

    Parameters
    ----------
    obs_values : dict {day_key: float}
        Observed anisotropy metric (EW−NS difference) per day or month.
    day_keys : list
        Sorted keys into obs_values.
    ra_bg, dec_bg : float
        Proposed background direction in degrees.

    Returns: (r, n_used)
    """
    if len(day_keys) < 5:
        return 0.0, 0

    n_bg = _direction_unit_vector(ra_bg, dec_bg)
    pred = []
    obs = []
    for day in day_keys:
        # Robust DOY extraction from YYYYDOY format (daily) or
        # month key (monthly, formatted as "2000MMM")
        try:
            if len(day) >= 7:
                doy = int(day[4:])
            else:
                doy = int(day[-3:])
        except (ValueError, IndexError):
            continue
        v = _earth_orbital_velocity_equatorial(doy)
        pred.append(float(np.dot(v, n_bg)))
        obs.append(obs_values[day])

    pred = np.array(pred)
    obs = np.array(obs)
    if len(pred) < 5 or np.std(pred) < 1e-12 or np.std(obs) < 1e-12:
        return 0.0, len(pred)
    r, _ = stats.pearsonr(pred, obs)
    return float(r), len(pred)


def _grid_search_signed(obs_values, day_keys, ra_grid, dec_grid, cmb_ra, cmb_dec):
    """Grid-search background direction, tracking both max +r and max −r.

    The aligned direction is the negative-correlation extremum (min r),
    selected by the sign convention established in step_2_3 (orbital
    coupling), where the EW−NS phase-alignment difference is anti-correlated
    with the CMB velocity projection (r ≈ −0.78).  This avoids circular
    selection of whichever extremum is geometrically closer to the CMB.

    Both extrema and their CMB separations are returned for reporting.

    Returns: dict with keys:
        best_ra, best_dec, best_r, best_r_abs, angular_separation_to_cmb_deg,
        max_pos_ra, max_pos_dec, max_pos_r, max_pos_sep_deg,
        max_neg_ra, max_neg_dec, max_neg_r, max_neg_sep_deg,
        aligned_ra, aligned_dec, aligned_r, aligned_sep_deg
    """
    best_r_abs = -np.inf
    best_r = None
    best_ra = best_dec = None
    max_r = -np.inf
    max_ra = max_dec = None
    min_r = np.inf
    min_ra = min_dec = None

    for ra in ra_grid:
        for dec in dec_grid:
            r, _ = _predicted_anisotropy(obs_values, day_keys, ra, dec)
            if abs(r) > best_r_abs:
                best_r_abs = abs(r)
                best_r = r
                best_ra = float(ra)
                best_dec = float(dec)
            if r > max_r:
                max_r = r
                max_ra = float(ra)
                max_dec = float(dec)
            if r < min_r:
                min_r = r
                min_ra = float(ra)
                min_dec = float(dec)

    sep_best = angular_separation(best_ra, best_dec, cmb_ra, cmb_dec)
    sep_pos = angular_separation(max_ra, max_dec, cmb_ra, cmb_dec)
    sep_neg = angular_separation(min_ra, min_dec, cmb_ra, cmb_dec)

    # Select the negative-correlation extremum by the sign convention from
    # step_2_3 (orbital coupling: pa_diff vs v_proj r ≈ −0.78).  This is a
    # pre-specified physical convention, not a circular selection by CMB
    # proximity.  Both extrema are reported so the reader can verify.
    align_ra, align_dec, align_r, align_sep = min_ra, min_dec, min_r, sep_neg

    return {
        "best_ra": best_ra,
        "best_dec": best_dec,
        "best_r": float(best_r),
        "best_r_abs": float(best_r_abs),
        "angular_separation_to_cmb_deg": float(sep_best),
        "max_pos_ra": max_ra,
        "max_pos_dec": max_dec,
        "max_pos_r": float(max_r),
        "max_pos_sep_deg": float(sep_pos),
        "max_neg_ra": min_ra,
        "max_neg_dec": min_dec,
        "max_neg_r": float(min_r),
        "max_neg_sep_deg": float(sep_neg),
        "aligned_ra": align_ra,
        "aligned_dec": align_dec,
        "aligned_r": float(align_r),
        "aligned_sep_deg": float(align_sep),
    }


class Step24CMBAlignment:
    def run(self):
        print_status("Step 2.4: CMB-Frame Alignment (Orbital-Velocity Geometry)", "INFO")
        results = {}
        for const_name in ["combined"]:
            pair_file = OUTPUTS_DIR / "step_2_0_mgex_pairs.json"
            if not pair_file.exists():
                results[const_name] = {"status": "no_data", "note": f"Pair file not found: {pair_file}"}
                continue

            with open(pair_file) as f:
                records = json.load(f)
            if not records:
                results[const_name] = {"status": "no_data", "note": "Empty pair records"}
                print_status(f"  {const_name}: empty pair records", "WARNING")
                continue

            print_status(f"  {const_name}: {len(records)} records loaded", "INFO")

            ra_grid = np.linspace(0, 360, CMB_GRID_RA_N)   # 10° steps
            dec_grid = np.linspace(-90, 90, CMB_GRID_DEC_N)   # 10° steps

            sector_sensitivity = {}
            for sector_width in CMB_SECTOR_WIDTHS:
                day_diffs = _ew_ns_difference(records, sector_width=sector_width)
                n_days = len(day_diffs)
                monthly = _monthly_from_daily(day_diffs)
                n_months = len(monthly)
                # Count pairs used
                n_ew_ns = sum(1 for r in records
                              if _sector_mask(r.get("azimuth_deg", 0), 90, width=sector_width)
                              or _sector_mask(r.get("azimuth_deg", 0), 270, width=sector_width)
                              or _sector_mask(r.get("azimuth_deg", 0), 0, width=sector_width)
                              or _sector_mask(r.get("azimuth_deg", 0), 180, width=sector_width))
                excluded_pct = 100.0 * (1.0 - n_ew_ns / len(records)) if records else 0.0
                print_status(
                    f"    {const_name} (sectors {sector_width}°): {n_days} days, "
                    f"{n_months} months, {n_ew_ns}/{len(records)} EW/NS pairs ({excluded_pct:.1f}% excluded)",
                    "INFO"
                )

                if n_days < MIN_DAYS_CMB:
                    sector_sensitivity[f"width_{sector_width}deg"] = {
                        "status": "insufficient_data",
                        "n_days": n_days,
                        "excluded_pairs_pct": excluded_pct,
                    }
                    print_status(f"      Insufficient data ({n_days} days)", "WARNING")
                    continue

                # --- Monthly grid search (primary) ---
                # Monthly aggregation suppresses a secondary daily systematic
                # and estimates the CMB-projected correlation. Because the
                # velocity samples an orbital plane, this is not a unique
                # three-dimensional sky-axis recovery.
                monthly_results = {}
                if n_months >= 5:
                    print_status(
                        f"      Monthly grid search: {n_months} months, "
                        f"{len(ra_grid)}x{len(dec_grid)} directions",
                        "INFO"
                    )
                    # Preserve each actual calendar month; do not pool the
                    # same month from different years.
                    monthly_obs = {}
                    monthly_keys = []
                    for m in sorted(monthly.keys()):
                        year, month = map(int, m.split("-"))
                        key = datetime(year, month, 15).strftime("%Y%j")
                        monthly_obs[key] = monthly[m]["diff"]
                        monthly_keys.append(key)
                    monthly_keys = sorted(monthly_keys)

                    gs_month = _grid_search_signed(
                        monthly_obs, monthly_keys, ra_grid, dec_grid,
                        CMB_DIPOLE_RA_DEG, CMB_DIPOLE_DEC_DEG
                    )

                    # Direct correlation at CMB dipole
                    r_cmb_month, n_cmb_month = _predicted_anisotropy(
                        monthly_obs, monthly_keys,
                        CMB_DIPOLE_RA_DEG, CMB_DIPOLE_DEC_DEG
                    )
                    print_status(
                        f"      Monthly: direct r(CMB)={r_cmb_month:.3f} (n={n_cmb_month}), "
                        f"aligned fit RA={gs_month['aligned_ra']:.0f}° Dec={gs_month['aligned_dec']:.0f}° "
                        f"r={gs_month['aligned_r']:.3f} Δθ={gs_month['aligned_sep_deg']:.1f}°",
                        "INFO"
                    )

                    # Monthly permutation null test
                    rng = np.random.default_rng(BOOTSTRAP_SEED)
                    null_r_month = []
                    print_status(
                        f"      Monthly null test: {CMB_PERMUTATION_N} permutations, "
                        f"full grid search each",
                        "INFO"
                    )
                    for perm_idx in range(CMB_PERMUTATION_N):
                        shuffled = rng.permutation(monthly_keys)
                        shuffled_obs = {k: monthly_obs[s] for k, s in zip(monthly_keys, shuffled)}
                        perm_best = -np.inf
                        for ra in ra_grid:
                            for dec in dec_grid:
                                r_null, _ = _predicted_anisotropy(shuffled_obs, monthly_keys, ra, dec)
                                if abs(r_null) > perm_best:
                                    perm_best = abs(r_null)
                        null_r_month.append(perm_best)
                        if (perm_idx + 1) % 50 == 0:
                            print_status(f"      Permutation {perm_idx + 1}/{CMB_PERMUTATION_N}", "INFO")

                    null_arr_m = np.array(null_r_month)
                    # Add one to numerator and denominator so a finite
                    # permutation run cannot report an impossible p = 0.
                    raw_p_month = (
                        np.sum(np.abs(null_arr_m) >= gs_month["best_r_abs"]) + 1
                    ) / (len(null_arr_m) + 1)
                    lee_p_month = raw_p_month
                    null_r_mean_m = float(np.mean(null_arr_m))
                    null_r_std_m = float(np.std(null_arr_m))
                    print_status(
                        f"      Monthly null: |r_null|={null_r_mean_m:.3f} ± {null_r_std_m:.3f}, "
                        f"p={raw_p_month:.4f}",
                        "INFO"
                    )

                    monthly_results = {
                        "n_months": n_months,
                        "direct_cmb_r": float(r_cmb_month),
                        "direct_cmb_n": int(n_cmb_month),
                        "best_fit_ra_deg": gs_month["aligned_ra"],
                        "best_fit_dec_deg": gs_month["aligned_dec"],
                        "best_correlation_r": gs_month["aligned_r"],
                        "angular_separation_to_cmb_deg": gs_month["aligned_sep_deg"],
                        "max_pos_ra_deg": gs_month["max_pos_ra"],
                        "max_pos_dec_deg": gs_month["max_pos_dec"],
                        "max_pos_r": gs_month["max_pos_r"],
                        "max_pos_sep_deg": gs_month["max_pos_sep_deg"],
                        "max_neg_ra_deg": gs_month["max_neg_ra"],
                        "max_neg_dec_deg": gs_month["max_neg_dec"],
                        "max_neg_r": gs_month["max_neg_r"],
                        "max_neg_sep_deg": gs_month["max_neg_sep_deg"],
                        "raw_p_value": float(raw_p_month),
                        "lee_corrected_p_value": float(lee_p_month),
                        "null_r_mean": null_r_mean_m,
                        "null_r_std": null_r_std_m,
                        # Kept for backward compatibility only: a sky-angle
                        # threshold cannot adjudicate alignment from coplanar
                        # orbital-velocity data.
                        "aligned": bool(
                            gs_month["aligned_sep_deg"] <= CMB_ALIGNMENT_THRESHOLD_DEG
                            and lee_p_month < SIGNIFICANCE_ALPHA
                        ),
                        "grid_axis_identifiable": False,
                        "cmb_projection_consistent": bool(r_cmb_month < 0),
                    }

                # --- Daily grid search (supplementary) ---
                daily_results = {}
                if n_days >= MIN_DAYS_CMB:
                    print_status(
                        f"      Daily grid search: {n_days} days, "
                        f"{len(ra_grid)}x{len(dec_grid)} directions",
                        "INFO"
                    )
                    daily_obs = {d: day_diffs[d]["diff"] for d in day_diffs}
                    daily_keys = sorted(daily_obs.keys())

                    gs_day = _grid_search_signed(
                        daily_obs, daily_keys, ra_grid, dec_grid,
                        CMB_DIPOLE_RA_DEG, CMB_DIPOLE_DEC_DEG
                    )

                    r_cmb_day, n_cmb_day = _predicted_anisotropy(
                        daily_obs, daily_keys,
                        CMB_DIPOLE_RA_DEG, CMB_DIPOLE_DEC_DEG
                    )
                    print_status(
                        f"      Daily: direct r(CMB)={r_cmb_day:.3f} (n={n_cmb_day}), "
                        f"aligned fit RA={gs_day['aligned_ra']:.0f}° Dec={gs_day['aligned_dec']:.0f}° "
                        f"r={gs_day['aligned_r']:.3f} Δθ={gs_day['aligned_sep_deg']:.1f}°",
                        "INFO"
                    )

                    daily_results = {
                        "n_days": n_days,
                        "direct_cmb_r": float(r_cmb_day),
                        "direct_cmb_n": int(n_cmb_day),
                        "best_fit_ra_deg": gs_day["aligned_ra"],
                        "best_fit_dec_deg": gs_day["aligned_dec"],
                        "best_correlation_r": gs_day["aligned_r"],
                        "angular_separation_to_cmb_deg": gs_day["aligned_sep_deg"],
                        "max_pos_ra_deg": gs_day["max_pos_ra"],
                        "max_pos_dec_deg": gs_day["max_pos_dec"],
                        "max_pos_r": gs_day["max_pos_r"],
                        "max_pos_sep_deg": gs_day["max_pos_sep_deg"],
                        "max_neg_ra_deg": gs_day["max_neg_ra"],
                        "max_neg_dec_deg": gs_day["max_neg_dec"],
                        "max_neg_r": gs_day["max_neg_r"],
                        "max_neg_sep_deg": gs_day["max_neg_sep_deg"],
                    }

                # The nominal full-sky coordinate is not identifiable from
                # orbital velocity alone; retain legacy angular fields only.
                primary_aligned = monthly_results.get("aligned", False) if monthly_results else False
                primary_sep = monthly_results.get("angular_separation_to_cmb_deg", 999) if monthly_results else 999
                primary_r = monthly_results.get("best_correlation_r", 0) if monthly_results else 0
                primary_p = monthly_results.get("lee_corrected_p_value", 1) if monthly_results else 1

                sector_sensitivity[f"width_{sector_width}deg"] = {
                    "status": "success",
                    "metric": "ew_ns_difference",
                    "monthly": monthly_results,
                    "daily": daily_results,
                    "aligned": primary_aligned,
                    "grid_axis_identifiable": False,
                    "cmb_projection_consistent": monthly_results.get("cmb_projection_consistent", False),
                    "angular_separation_to_cmb_deg": float(primary_sep),
                    "best_correlation_r": float(primary_r),
                    "raw_p_value": float(primary_p),
                    "lee_corrected_p_value": float(primary_p),
                    "n_days": n_days,
                    "n_months": n_months,
                    "excluded_pairs_pct": excluded_pct,
                }
                print_status(
                    f"      Result (sectors {sector_width}°): "
                    f"Δθ_CMB={primary_sep:.1f}°, r={primary_r:.3f}, p={primary_p:.4f} "
                    "[CMB projection consistent; sky axis not identifiable]",
                    "INFO"
                )

            # Use 45° as canonical result for downstream compatibility
            primary = sector_sensitivity.get("width_45deg", {})

            # PA-difference sensitivity to minimum baseline length.
            # The PA-difference metric (EW mean − NS mean) is a simple
            # average, not a nonlinear fit, so it should be robust to
            # the short-baseline cutoff.  This checks that the CMB-projection
            # result does not depend on noisy short baselines.
            pa_sensitivity = {}
            dist_arr = np.array([r["distance_km"] for r in records])
            for min_d in SENSITIVITY_MIN_DISTANCES:
                mask = dist_arr >= min_d
                filtered = [r for r, m in zip(records, mask) if m]
                if len(filtered) < 1000:
                    pa_sensitivity[f"min_{min_d}km"] = {
                        "min_distance_km": min_d, "status": "insufficient_data"
                    }
                    continue
                day_diffs_s = _ew_ns_difference(filtered, sector_width=SECTOR_WIDTH_DEG)
                monthly_s = _monthly_from_daily(day_diffs_s)
                if len(monthly_s) >= 5:
                    monthly_obs_s = {}
                    monthly_keys_s = []
                    for m in sorted(monthly_s.keys()):
                        year, month = map(int, m.split("-"))
                        key = datetime(year, month, 15).strftime("%Y%j")
                        monthly_obs_s[key] = monthly_s[m]["diff"]
                        monthly_keys_s.append(key)
                    monthly_keys_s = sorted(monthly_keys_s)
                    r_cmb, _ = _predicted_anisotropy(
                        monthly_obs_s, monthly_keys_s,
                        CMB_DIPOLE_RA_DEG, CMB_DIPOLE_DEC_DEG
                    )
                    pa_sensitivity[f"min_{min_d}km"] = {
                        "min_distance_km": min_d,
                        "n_months": len(monthly_s),
                        "direct_cmb_r": float(r_cmb),
                    }
                    print_status(
                        f"      PA-diff sensitivity (min={min_d}km): "
                        f"r(CMB)={r_cmb:.3f}, n_months={len(monthly_s)}",
                        "INFO"
                    )

            if primary.get("status") == "success":
                results[const_name] = {
                    "best_fit_ra_deg": primary.get("best_fit_ra_deg") or primary.get("monthly", {}).get("best_fit_ra_deg"),
                    "best_fit_dec_deg": primary.get("best_fit_dec_deg") or primary.get("monthly", {}).get("best_fit_dec_deg"),
                    "best_correlation_r": primary.get("best_correlation_r"),
                    "angular_separation_to_cmb_deg": primary.get("angular_separation_to_cmb_deg"),
                    "raw_p_value": primary.get("raw_p_value"),
                    "lee_corrected_p_value": primary.get("lee_corrected_p_value"),
                    "status": "success",
                    "aligned": primary.get("aligned"),
                    "grid_axis_identifiable": primary.get("grid_axis_identifiable", False),
                    "cmb_projection_consistent": primary.get("cmb_projection_consistent", False),
                    "n_days": primary.get("n_days"),
                    "n_months": primary.get("n_months"),
                    "n_pairs": len(records),
                    "metric": "ew_ns_difference",
                    "sector_sensitivity": sector_sensitivity,
                    "pa_difference_sensitivity": pa_sensitivity,
                }
            else:
                results[const_name] = {
                    "status": "insufficient_data",
                    "note": f"Only {primary.get('n_days', 0)} days with EW/NS differences (45° sectors)",
                    "sector_sensitivity": sector_sensitivity,
                }

        out_file = OUTPUTS_DIR / "step_2_4_cmb_alignment.json"
        with open(out_file, 'w') as f:
            json.dump(results, f, indent=2)
        print_status(f"CMB alignment results written to {out_file}", "INFO")
        return {"status": "success", "results": str(out_file)}
