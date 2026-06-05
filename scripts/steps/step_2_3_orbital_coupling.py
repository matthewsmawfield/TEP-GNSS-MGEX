#!/usr/bin/env python3
"""Step 2.3: Orbital-velocity coupling.

Frozen prediction: monthly λ / EW–NS ratio correlates with Earth orbital speed.
Test: stratify pair records by month, fit λ and EW/NS per month, regress
against Earth orbital speed.
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
    OUTPUTS_DIR, LOGS_DIR, CONSTELLATIONS, N_ORBITAL_TESTS,
    CMB_DIPOLE_RA_DEG, CMB_DIPOLE_DEC_DEG, SIGNIFICANCE_ALPHA,
    MIN_DISTANCE_KM, MAX_DISTANCE_KM, N_BINS, MIN_BIN_COUNT,
    MIN_PAIRS_FIT, MIN_MONTHLY_RECORDS
)
from scripts.utils.coherence import fit_exponential_model
from scripts.utils.geospatial import _sector_mask, _earth_orbital_velocity_equatorial

LOGS_DIR.mkdir(parents=True, exist_ok=True)
logger = TEPLogger("step_2_3", log_file_path=LOGS_DIR / "step_2_3_orbital_coupling.log")
set_step_logger(logger)


def _month_from_day_key(day_key):
    """Convert YYYYDOY string to month integer (1-12)."""
    try:
        year = int(day_key[:4])
        doy = int(day_key[4:])
        dt = datetime(year, 1, 1) + timedelta(days=doy - 1)
        return dt.month
    except Exception:
        return None


def _direction_unit_vector(ra_deg, dec_deg):
    """Convert RA/Dec to unit vector in equatorial coordinates."""
    ra, dec = np.radians(ra_deg), np.radians(dec_deg)
    return np.array([
        np.cos(dec) * np.cos(ra),
        np.cos(dec) * np.sin(ra),
        np.sin(dec)
    ])


def earth_orbital_speed_kms(month):
    """Approximate Earth orbital speed variation by month (sinusoidal model).
    Kept for backward compatibility; the pipeline now uses velocity projection.
    Peak ~30.3 km/s in January, minimum ~29.3 km/s in July.
    """
    return 29.78 + 0.5 * np.cos(2 * np.pi * (month - 1) / 12)


def _monthly_velocity_projection(month, cmb_ra=167.0, cmb_dec=-7.0):
    """Earth orbital velocity projected onto the CMB dipole direction (km/s).

    Positive when Earth moves toward the CMB dipole; negative when moving away.
    The TEP prediction is that the anisotropy correlates with this projection,
    not with scalar speed.
    """
    # Representative DOY for each month (15th of month)
    doy = int((datetime(2000, month, 15) - datetime(2000, 1, 1)).days + 1)
    v_orb = _earth_orbital_velocity_equatorial(doy)
    n_cmb = _direction_unit_vector(cmb_ra, cmb_dec)
    return float(np.dot(v_orb, n_cmb))


class Step23OrbitalCoupling:
    def run(self):
        print_status("Step 2.3: Orbital-Velocity Coupling", "INFO")
        results = {}
        for const_name in CONSTELLATIONS:
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

            # Group records by month (from day_key YYYYDOY)
            month_records = {}
            for r in records:
                month = _month_from_day_key(r.get("day", ""))
                if month is None:
                    continue
                if month not in month_records:
                    month_records[month] = []
                month_records[month].append(r)

            print_status(f"  {const_name}: {len(month_records)} months with data", "INFO")

            months = []
            speeds_scalar = []
            speeds_proj = []
            lambdas = []
            ew_ns_ratios = []
            mean_pa_diffs = []  # Supplementary metric: direct EW − NS mean phase_alignment

            for m in sorted(month_records.keys()):
                recs = month_records[m]
                if len(recs) < MIN_MONTHLY_RECORDS:
                    continue
                dist = np.array([r["distance_km"] for r in recs])
                coh = np.array([r["phase_alignment"] for r in recs])
                fit = fit_exponential_model(dist, coh, MIN_DISTANCE_KM, MAX_DISTANCE_KM, N_BINS, MIN_BIN_COUNT)
                if fit and fit["success"]:
                    # EW/NS ratio for this month — fit λ separately for EW and NS sectors
                    # This matches the step_2_2 methodology: ratio = λ_EW / λ_NS
                    az_month = np.array([r.get("azimuth_deg", 0) for r in recs])
                    ew_mask = _sector_mask(az_month, 90) | _sector_mask(az_month, 270)
                    ns_mask = _sector_mask(az_month, 0) | _sector_mask(az_month, 180)

                    fit_ew = fit_exponential_model(
                        dist[ew_mask], coh[ew_mask],
                        MIN_DISTANCE_KM, MAX_DISTANCE_KM, N_BINS, MIN_BIN_COUNT
                    ) if np.sum(ew_mask) >= MIN_PAIRS_FIT else None
                    fit_ns = fit_exponential_model(
                        dist[ns_mask], coh[ns_mask],
                        MIN_DISTANCE_KM, MAX_DISTANCE_KM, N_BINS, MIN_BIN_COUNT
                    ) if np.sum(ns_mask) >= MIN_PAIRS_FIT else None

                    ew_lambda = fit_ew["correlation_length_km"] if fit_ew and fit_ew["success"] else None
                    ns_lambda = fit_ns["correlation_length_km"] if fit_ns and fit_ns["success"] else None
                    ratio = float(ew_lambda / ns_lambda) if (ew_lambda is not None and ns_lambda is not None and ns_lambda > 0) else np.nan

                    # Supplementary metric: direct mean PA difference (no exponential fit)
                    ew_mean_pa = float(np.mean(coh[ew_mask])) if np.sum(ew_mask) > 0 else np.nan
                    ns_mean_pa = float(np.mean(coh[ns_mask])) if np.sum(ns_mask) > 0 else np.nan
                    pa_diff = float(ew_mean_pa - ns_mean_pa) if np.isfinite(ew_mean_pa) and np.isfinite(ns_mean_pa) else np.nan

                    months.append(m)
                    # Original test: scalar orbital speed (Paper 1-3 methodology)
                    speed_scalar = earth_orbital_speed_kms(m)
                    # Audit addition: velocity-vector projection onto CMB dipole
                    v_proj = _monthly_velocity_projection(m, CMB_DIPOLE_RA_DEG, CMB_DIPOLE_DEC_DEG)
                    speeds_scalar.append(speed_scalar)
                    speeds_proj.append(v_proj)
                    lambdas.append(fit["correlation_length_km"])
                    ew_ns_ratios.append(ratio)
                    mean_pa_diffs.append(pa_diff)
                    ratio_str = f"{ratio:.2f}" if np.isfinite(ratio) else "N/A"
                    pa_diff_str = f"{pa_diff:.4f}" if np.isfinite(pa_diff) else "N/A"
                    print_status(
                        f"    Month {m:02d}: λ={fit['correlation_length_km']:.0f}km, "
                        f"EW/NS ratio={ratio_str}, PA_diff={pa_diff_str}, "
                        f"speed={speed_scalar:.2f}km/s, v_proj={v_proj:+.2f}km/s",
                        "INFO"
                    )
                else:
                    print_status(f"    Month {m:02d}: fit failed (n={len(recs)})", "WARNING")

            lambdas_arr = np.array(lambdas)
            speeds_scalar_arr = np.array(speeds_scalar)
            speeds_proj_arr = np.array(speeds_proj)
            ratios_arr = np.array(ew_ns_ratios)
            pa_diffs_arr = np.array(mean_pa_diffs)
            valid = ~np.isnan(lambdas_arr) & ~np.isnan(ratios_arr)
            valid_pa = ~np.isnan(pa_diffs_arr)
            n_valid = int(np.sum(valid))
            n_valid_pa = int(np.sum(valid_pa))
            print_status(f"  {const_name}: {n_valid} valid monthly fits, {n_valid_pa} valid PA_diff values", "INFO")

            if n_valid >= 3 or n_valid_pa >= 3:
                # Scalar speed correlations (original Paper 1-3 test)
                r_lambda_s, p_lambda_s = stats.pearsonr(speeds_scalar_arr[valid], lambdas_arr[valid]) if n_valid >= 3 else (np.nan, np.nan)
                r_ratio_s, p_ratio_s = stats.pearsonr(speeds_scalar_arr[valid], ratios_arr[valid]) if n_valid >= 3 else (np.nan, np.nan)
                # Velocity projection correlations (audit addition)
                r_lambda_p, p_lambda_p = stats.pearsonr(speeds_proj_arr[valid], lambdas_arr[valid]) if n_valid >= 3 else (np.nan, np.nan)
                r_ratio_p, p_ratio_p = stats.pearsonr(speeds_proj_arr[valid], ratios_arr[valid]) if n_valid >= 3 else (np.nan, np.nan)
                # Supplementary: PA-diff correlations (short-baseline robust metric)
                r_pa_s, p_pa_s = stats.pearsonr(speeds_scalar_arr[valid_pa], pa_diffs_arr[valid_pa]) if n_valid_pa >= 3 else (np.nan, np.nan)
                r_pa_p, p_pa_p = stats.pearsonr(speeds_proj_arr[valid_pa], pa_diffs_arr[valid_pa]) if n_valid_pa >= 3 else (np.nan, np.nan)
                # Bonferroni correction: two correlations tested per predictor for λ/ratio
                p_lambda_s_bonf = min(float(p_lambda_s) * N_ORBITAL_TESTS, 1.0) if n_valid >= 3 else np.nan
                p_ratio_s_bonf = min(float(p_ratio_s) * N_ORBITAL_TESTS, 1.0) if n_valid >= 3 else np.nan
                p_lambda_p_bonf = min(float(p_lambda_p) * N_ORBITAL_TESTS, 1.0) if n_valid >= 3 else np.nan
                p_ratio_p_bonf = min(float(p_ratio_p) * N_ORBITAL_TESTS, 1.0) if n_valid >= 3 else np.nan
                # PA-diff is a single supplementary metric per predictor; no Bonferroni needed
                # (it is not one of multiple tests for the same hypothesis)
                results[const_name] = {
                    "orbital_speed_kms": float(np.mean(speeds_scalar_arr[valid_pa])),
                    "velocity_projection_kms": float(np.mean(speeds_proj_arr[valid_pa])),
                    # Scalar speed (original methodology)
                    "lambda_vs_speed_r": float(r_lambda_s),
                    "lambda_vs_speed_p": float(p_lambda_s),
                    "lambda_vs_speed_p_bonferroni": p_lambda_s_bonf,
                    "ew_ns_ratio_vs_speed_r": float(r_ratio_s),
                    "ew_ns_ratio_vs_speed_p": float(p_ratio_s),
                    "ew_ns_ratio_vs_speed_p_bonferroni": p_ratio_s_bonf,
                    # Velocity projection (audit addition)
                    "lambda_vs_vproj_r": float(r_lambda_p),
                    "lambda_vs_vproj_p": float(p_lambda_p),
                    "lambda_vs_vproj_p_bonferroni": p_lambda_p_bonf,
                    "ew_ns_ratio_vs_vproj_r": float(r_ratio_p),
                    "ew_ns_ratio_vs_vproj_p": float(p_ratio_p),
                    "ew_ns_ratio_vs_vproj_p_bonferroni": p_ratio_p_bonf,
                    # Supplementary PA-diff (short-baseline robust metric)
                    "pa_diff_vs_speed_r": float(r_pa_s),
                    "pa_diff_vs_speed_p": float(p_pa_s),
                    "pa_diff_vs_vproj_r": float(r_pa_p),
                    "pa_diff_vs_vproj_p": float(p_pa_p),
                    "n_tests_corrected": N_ORBITAL_TESTS,
                    "significant_after_correction": bool(
                        (n_valid >= 3 and (p_lambda_s_bonf < SIGNIFICANCE_ALPHA or p_ratio_s_bonf < SIGNIFICANCE_ALPHA or
                                           p_lambda_p_bonf < SIGNIFICANCE_ALPHA or p_ratio_p_bonf < SIGNIFICANCE_ALPHA)) or
                        (n_valid_pa >= 3 and (p_pa_s < SIGNIFICANCE_ALPHA or p_pa_p < SIGNIFICANCE_ALPHA))
                    ),
                    "status": "success",
                    "n_months": int(np.sum(valid)),
                    "n_months_pa_diff": int(np.sum(valid_pa)),
                    "months": [int(m) for m in months]
                }
                print_status(
                    f"    {const_name}: λ–speed r={r_lambda_s:.3f} (p={p_lambda_s:.3f}, bonf={p_lambda_s_bonf:.3f}), "
                    f"ratio–speed r={r_ratio_s:.3f} (p={p_ratio_s:.3f}, bonf={p_ratio_s_bonf:.3f}) | "
                    f"λ–v_proj r={r_lambda_p:.3f} (p={p_lambda_p:.3f}, bonf={p_lambda_p_bonf:.3f}), "
                    f"ratio–v_proj r={r_ratio_p:.3f} (p={p_ratio_p:.3f}, bonf={p_ratio_p_bonf:.3f}) | "
                    f"PA_diff–speed r={r_pa_s:.3f} (p={p_pa_s:.4f}), "
                    f"PA_diff–v_proj r={r_pa_p:.3f} (p={p_pa_p:.4f})",
                    "INFO"
                )
            else:
                results[const_name] = {"status": "insufficient_data", "note": "Too few valid monthly fits"}

        out_file = OUTPUTS_DIR / "step_2_3_orbital_coupling.json"
        with open(out_file, 'w') as f:
            json.dump(results, f, indent=2)
        print_status(f"Orbital coupling results written to {out_file}", "INFO")
        return {"status": "success", "results": str(out_file)}
