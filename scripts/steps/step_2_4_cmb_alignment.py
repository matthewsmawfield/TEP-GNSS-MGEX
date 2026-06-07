#!/usr/bin/env python3
"""Step 2.4: CMB-frame alignment via orbital-velocity geometry.

Frozen prediction: daily EW/NS anisotropy ratio correlates with Earth's
orbital velocity projected onto the CMB dipole direction.

Methodology (from TEP-GNSS-II Cairo):
    1. Aggregate pair phase_alignment values by day → daily EW/NS ratio
    2. For each day, compute Earth orbital velocity vector in equatorial coords
    3. For a proposed background direction (RA, Dec), compute predicted
       anisotropy = dot(v_orb, n_background)
    4. Correlate predicted vs observed daily ratios
    5. Grid-search background directions; compare best-fit to CMB dipole
    6. Permutation null test + look-elsewhere correction
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
    OUTPUTS_DIR, LOGS_DIR, CONSTELLATIONS,
    CMB_DIPOLE_RA_DEG, CMB_DIPOLE_DEC_DEG,
    LOOK_ELSEWHERE_FACTOR, CMB_ALIGNMENT_THRESHOLD_DEG, SIGNIFICANCE_ALPHA,
    CMB_GRID_RA_N, CMB_GRID_DEC_N, CMB_PERMUTATION_N, CMB_SECTOR_WIDTHS,
    MIN_DAYS_CMB, BOOTSTRAP_SEED, SECTOR_WIDTH_DEG
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


def _daily_ew_ns_ratios(records, sector_width=SECTOR_WIDTH_DEG):
    """Aggregate pair records by day into EW/NS phase_alignment ratios.

    Uses the same 8-sector definitions as step_2_2 and step_2_3:
    EW = E (center 90°) + W (center 270°)
    NS = N (center 0°) + S (center 180°)

    Parameters
    ----------
    records : list of dict
        Pair records with 'day', 'azimuth_deg', 'phase_alignment' keys.
    sector_width : int
        Full width of azimuth sector in degrees (default 45°).

    Returns: dict {day_key: {"ew": float, "ns": float, "ratio": float}}
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

    ratios = {}
    for day, vals in day_data.items():
        ew = np.mean(vals["ew"]) if vals["ew"] else np.nan
        ns = np.mean(vals["ns"]) if vals["ns"] else np.nan
        if np.isfinite(ew) and np.isfinite(ns) and ns > 0:
            ratios[day] = {"ew": float(ew), "ns": float(ns), "ratio": float(ew / ns)}
    return ratios


def _predicted_anisotropy(day_ratios, ra_bg, dec_bg):
    """Compute Pearson r between predicted (velocity·background) and observed EW/NS ratio."""
    days = sorted(day_ratios.keys())
    if len(days) < 5:
        return 0.0, days

    n_bg = _direction_unit_vector(ra_bg, dec_bg)
    pred = []
    obs = []
    for day in days:
        # Robust DOY extraction from YYYYDOY format
        try:
            doy = int(day[4:]) if len(day) >= 7 else int(day[-3:])
        except (ValueError, IndexError):
            continue
        v = _earth_orbital_velocity_equatorial(doy)
        pred.append(float(np.dot(v, n_bg)))
        obs.append(day_ratios[day]["ratio"])

    pred = np.array(pred)
    obs = np.array(obs)
    if np.std(pred) < 1e-12 or np.std(obs) < 1e-12:
        return 0.0, days
    r, _ = stats.pearsonr(pred, obs)
    return float(r), days


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

            # Compute with default 45° sectors and also 60° to expose sensitivity
            sector_sensitivity = {}
            for sector_width in CMB_SECTOR_WIDTHS:
                day_ratios = _daily_ew_ns_ratios(records, sector_width=sector_width)
                n_days = len(day_ratios)
                # Count pairs used
                n_ew_ns = sum(1 for r in records
                              if _sector_mask(r.get("azimuth_deg", 0), 90, width=sector_width)
                              or _sector_mask(r.get("azimuth_deg", 0), 270, width=sector_width)
                              or _sector_mask(r.get("azimuth_deg", 0), 0, width=sector_width)
                              or _sector_mask(r.get("azimuth_deg", 0), 180, width=sector_width))
                excluded_pct = 100.0 * (1.0 - n_ew_ns / len(records)) if records else 0.0
                print_status(
                    f"    {const_name} (sectors {sector_width}°): {n_days} days, "
                    f"{n_ew_ns}/{len(records)} EW/NS pairs ({excluded_pct:.1f}% excluded)",
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

                # Grid search over background directions
                ra_grid = np.linspace(0, 360, CMB_GRID_RA_N)   # 10° steps
                dec_grid = np.linspace(-90, 90, CMB_GRID_DEC_N)   # 10° steps
                print_status(
                    f"      Grid search: {len(ra_grid)}x{len(dec_grid)} directions",
                    "INFO"
                )
                best_r_abs = -np.inf
                best_r = None
                best_ra = best_dec = None

                for ra in ra_grid:
                    for dec in dec_grid:
                        r, _ = _predicted_anisotropy(day_ratios, ra, dec)
                        if abs(r) > best_r_abs:
                            best_r_abs = abs(r)
                            best_r = r
                            best_ra = float(ra)
                            best_dec = float(dec)

                sep = angular_separation(best_ra, best_dec, CMB_DIPOLE_RA_DEG, CMB_DIPOLE_DEC_DEG)

                # Null: shuffle day labels
                rng = np.random.default_rng(BOOTSTRAP_SEED)
                null_r = []
                days_list = list(day_ratios.keys())
                for _ in range(CMB_PERMUTATION_N):
                    shuffled_days = rng.permutation(days_list)
                    shuffled_ratios = {d: day_ratios[sd] for d, sd in zip(days_list, shuffled_days)}
                    r_null, _ = _predicted_anisotropy(shuffled_ratios, best_ra, best_dec)
                    null_r.append(r_null)

                # Two-tailed: test whether |null| >= |best|
                null_arr = np.array(null_r)
                raw_p = np.mean(np.abs(null_arr) >= np.abs(best_r))
                lee_p = min(raw_p * LOOK_ELSEWHERE_FACTOR, 1.0)
                null_r_mean = float(np.mean(null_arr))
                null_r_std = float(np.std(null_arr))
                print_status(
                    f"      Null test: |r_null| = {null_r_mean:.3f} ± {null_r_std:.3f} "
                    f"(n=200), raw p={raw_p:.4f}, LEE p={lee_p:.4f}",
                    "INFO"
                )

                sector_sensitivity[f"width_{sector_width}deg"] = {
                    "status": "success",
                    "best_fit_ra_deg": best_ra,
                    "best_fit_dec_deg": best_dec,
                    "best_correlation_r": float(best_r),
                    "angular_separation_to_cmb_deg": float(sep),
                    "raw_p_value": float(raw_p),
                    "lee_corrected_p_value": float(lee_p),
                    "aligned": bool(sep <= CMB_ALIGNMENT_THRESHOLD_DEG and lee_p < SIGNIFICANCE_ALPHA),
                    "n_days": n_days,
                    "excluded_pairs_pct": excluded_pct,
                }
                print_status(
                    f"      Best fit: RA={best_ra:.1f}°, Dec={best_dec:.1f}°, r={best_r:.3f}, "
                    f"Δθ_CMB={sep:.1f}° [{'ALIGNED' if sep <= CMB_ALIGNMENT_THRESHOLD_DEG else 'NOT ALIGNED'}]",
                    "SUCCESS" if sep <= CMB_ALIGNMENT_THRESHOLD_DEG else "INFO"
                )

            # Use 45° as canonical result for downstream compatibility
            primary = sector_sensitivity.get("width_45deg", {})
            if primary.get("status") == "success":
                results[const_name] = {
                    "best_fit_ra_deg": primary.get("best_fit_ra_deg"),
                    "best_fit_dec_deg": primary.get("best_fit_dec_deg"),
                    "best_correlation_r": primary.get("best_correlation_r"),
                    "angular_separation_to_cmb_deg": primary.get("angular_separation_to_cmb_deg"),
                    "raw_p_value": primary.get("raw_p_value"),
                    "lee_corrected_p_value": primary.get("lee_corrected_p_value"),
                    "status": "success",
                    "aligned": primary.get("aligned"),
                    "n_days": primary.get("n_days"),
                    "n_pairs": len(records),
                    "sector_sensitivity": sector_sensitivity,
                }
            else:
                results[const_name] = {
                    "status": "insufficient_data",
                    "note": f"Only {primary.get('n_days', 0)} days with EW/NS ratios (45° sectors)",
                    "sector_sensitivity": sector_sensitivity,
                }

        out_file = OUTPUTS_DIR / "step_2_4_cmb_alignment.json"
        with open(out_file, 'w') as f:
            json.dump(results, f, indent=2)
        print_status(f"CMB alignment results written to {out_file}", "INFO")
        return {"status": "success", "results": str(out_file)}
