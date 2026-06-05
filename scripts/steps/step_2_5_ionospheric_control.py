#!/usr/bin/env python3
"""Step 2.5: Ionospheric control.

Frozen prediction: TEP signal is not ionospheric.
Test: Kp stratification and storm-day exclusion using authentic space weather data.

Aligned with TEP-GNSS-II methodology:
- Fetches Kp from GFZ Potsdam (historical) or NOAA SWPC (recent)
- Stratifies pair records by daily Kp into quiet / active / storm-excluded subsets
- Fits lambda for each subset; TEP prediction: signal persists in quiet conditions.
"""

import sys
from pathlib import Path
from datetime import date

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

import json
import numpy as np
from scripts.utils.logger import TEPLogger, set_step_logger, print_status
from scripts.utils.config import OUTPUTS_DIR, LOGS_DIR, CONSTELLATIONS, KP_QUIET, KP_ACTIVE, IONO_PERSISTENCE_THRESHOLD, MIN_DISTANCE_KM, MAX_DISTANCE_KM, N_BINS, MIN_BIN_COUNT
from scripts.utils.coherence import fit_exponential_model
from scripts.utils.space_weather_data import get_kp_data_for_range, _day_key_to_date

LOGS_DIR.mkdir(parents=True, exist_ok=True)
logger = TEPLogger("step_2_5", log_file_path=LOGS_DIR / "step_2_5_ionospheric_control.log")
set_step_logger(logger)


def _fit_lambda(records):
    """Fit correlation length from a list of pair records."""
    if not records:
        return None
    dist = np.array([r["distance_km"] for r in records])
    # Use phase_alignment (TEP signature metric) instead of coherence (isotropic amplitude)
    coh = np.array([r["phase_alignment"] for r in records])
    fit = fit_exponential_model(dist, coh, MIN_DISTANCE_KM, MAX_DISTANCE_KM, N_BINS, MIN_BIN_COUNT)
    return fit["correlation_length_km"] if fit and fit["success"] else None


def _fit_lambda_by_kp(records, kp_map, kp_threshold):
    """Fit lambda for quiet records with Kp <= threshold."""
    filtered = [r for r in records if kp_map.get(r["day"], 99) <= kp_threshold]
    return _fit_lambda(filtered), len(filtered)


def _fit_lambda_by_kp_min(records, kp_map, kp_threshold):
    """Fit lambda for active records with Kp >= threshold.

    Note: default of -1 ensures days with no Kp coverage are excluded
    from the active subset (they cannot be confirmed as active)."""
    filtered = [r for r in records if kp_map.get(r["day"], -1) >= kp_threshold]
    return _fit_lambda(filtered), len(filtered)


def _fit_lambda_excluding_storms(records, kp_map, storm_threshold):
    """Fit lambda excluding records from storm days (Kp >= threshold)."""
    filtered = [r for r in records if kp_map.get(r["day"], 99) < storm_threshold]
    return _fit_lambda(filtered), len(filtered)


class Step25IonosphericControl:
    def run(self):
        print_status("Step 2.5: Ionospheric Control", "INFO")
        results = {}

        pair_file = OUTPUTS_DIR / "step_2_0_mgex_pairs.json"
        if not pair_file.exists():
            print_status(f"Pair file not found: {pair_file}", "ERROR")
            for c in CONSTELLATIONS:
                results[c] = {"status": "no_data", "note": f"Pair file not found: {pair_file}"}
            out_file = OUTPUTS_DIR / "step_2_5_ionospheric_control.json"
            with open(out_file, 'w') as f:
                json.dump(results, f, indent=2)
            return {"status": "no_data", "results": str(out_file)}

        with open(pair_file) as f:
            records = json.load(f)

        if not records:
            for c in CONSTELLATIONS:
                results[c] = {"status": "no_data", "note": "Empty pair records"}
            out_file = OUTPUTS_DIR / "step_2_5_ionospheric_control.json"
            with open(out_file, 'w') as f:
                json.dump(results, f, indent=2)
            return {"status": "no_data", "results": str(out_file)}

        # Determine date range from pair records
        days = sorted(set(r["day"] for r in records))
        start_date = _day_key_to_date(days[0])
        end_date = _day_key_to_date(days[-1])
        print_status(f"  Loaded {len(records)} pair records over {len(days)} days ({start_date} to {end_date})", "INFO")

        # Fetch authentic Kp data
        kp_map = get_kp_data_for_range(start_date, end_date)
        has_kp = bool(kp_map)
        n_kp_matched = sum(1 for d in days if d in kp_map)
        print_status(f"  Kp coverage: {n_kp_matched}/{len(days)} days matched", "INFO" if has_kp else "WARNING")

        # Base lambda (all records)
        base_lambda = _fit_lambda(records)
        n_total = len(records)
        base_str = f"{base_lambda:.0f}" if base_lambda is not None else "N/A"
        print_status(f"  Base lambda (all records): {base_str} km (n={n_total})", "INFO" if base_lambda is not None else "WARNING")

        for const_name in CONSTELLATIONS:
            entry = {
                "base_lambda_km": base_lambda,
                "base_n_pairs": n_total,
                "ionofree_lambda_km": base_lambda,  # MGEX CLK is iono-free by construction
                "status": "success" if base_lambda is not None else "no_data",
            }

            if has_kp and base_lambda is not None:
                # Quiet: Kp <= KP_QUIET (typically 2)
                quiet_lambda, quiet_n = _fit_lambda_by_kp(records, kp_map, KP_QUIET)
                # Active: Kp >= KP_ACTIVE (typically 5)
                active_lambda, active_n = _fit_lambda_by_kp_min(records, kp_map, KP_ACTIVE)
                # Storm-excluded: Kp < KP_ACTIVE
                storm_excl_lambda, storm_n = _fit_lambda_excluding_storms(records, kp_map, KP_ACTIVE)

                entry.update({
                    "quiet_kp_lambda_km": quiet_lambda,
                    "quiet_n_pairs": quiet_n,
                    "active_kp_lambda_km": active_lambda,
                    "active_n_pairs": active_n,
                    "storm_excluded_lambda_km": storm_excl_lambda,
                    "storm_excluded_n_pairs": storm_n,
                    "kp_quiet_threshold": KP_QUIET,
                    "kp_active_threshold": KP_ACTIVE,
                    "kp_coverage_days": n_kp_matched,
                    "total_days": len(days),
                })

                # TEP prediction: signal persists in quiet conditions
                quiet_persist = (
                    quiet_lambda is not None
                    and base_lambda is not None
                    and base_lambda != 0
                    and abs(quiet_lambda - base_lambda) / base_lambda < IONO_PERSISTENCE_THRESHOLD
                )
                entry["quiet_persistence"] = bool(quiet_persist)

                quiet_str = f"{quiet_lambda:.0f}" if quiet_lambda is not None else "N/A"
                active_str = f"{active_lambda:.0f}" if active_lambda is not None else "N/A"
                storm_str = f"{storm_excl_lambda:.0f}" if storm_excl_lambda is not None else "N/A"
                print_status(
                    f"    {const_name}: base={base_lambda:.0f}km, "
                    f"quiet(Kp<={KP_QUIET})={quiet_str}km (n={quiet_n}), "
                    f"active(Kp>={KP_ACTIVE})={active_str}km (n={active_n}), "
                    f"storm-excl={storm_str}km (n={storm_n})",
                    "INFO"
                )
                persist_str = "PASS" if quiet_persist else "FAIL"
                if quiet_lambda is not None and base_lambda is not None and base_lambda != 0:
                    persist_msg = (
                        f"      Persistence check: |quiet-base|/base = "
                        f"{abs(quiet_lambda - base_lambda) / base_lambda:.2%} [{persist_str}]"
                    )
                else:
                    persist_msg = f"      Persistence check: N/A [{persist_str}]"
                print_status(persist_msg, "INFO")
            else:
                entry["note"] = "Kp data unavailable; stratification skipped"
                base_str_else = f"{base_lambda:.0f}" if base_lambda is not None else "N/A"
                print_status(f"    {const_name}: base={base_str_else}km (Kp data unavailable)", "WARNING")

            results[const_name] = entry

        out_file = OUTPUTS_DIR / "step_2_5_ionospheric_control.json"
        with open(out_file, 'w') as f:
            json.dump(results, f, indent=2)
        print_status(f"Ionospheric control results written to {out_file}", "INFO")
        return {"status": "success", "results": str(out_file)}
