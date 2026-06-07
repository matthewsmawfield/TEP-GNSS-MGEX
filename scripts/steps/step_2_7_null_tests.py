#!/usr/bin/env python3
"""Step 2.7: Null control tests.

Tests: temporal shuffle, spatial shuffle, phase randomization, solar rotation.
Expected null for all control tests.
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

import json
import numpy as np
from scripts.utils.logger import TEPLogger, set_step_logger, print_status
from scripts.utils.config import (
    OUTPUTS_DIR, LOGS_DIR, CONSTELLATIONS,
    MIN_DISTANCE_KM, MAX_DISTANCE_KM, N_BINS, MIN_BIN_COUNT,
    NULL_TEST_R2_THRESHOLD, SIGNIFICANCE_ALPHA, NULL_TEST_SEED
)
from scripts.utils.coherence import fit_exponential_model
from scripts.utils.geospatial import _sector_mask

LOGS_DIR.mkdir(parents=True, exist_ok=True)
logger = TEPLogger("step_2_7", log_file_path=LOGS_DIR / "step_2_7_null_tests.log")
set_step_logger(logger)


class Step27NullTests:
    def run(self):
        print_status("Step 2.7: Null Control Tests", "INFO")
        results = {}
        rng = np.random.default_rng(NULL_TEST_SEED)

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

            dist = np.array([r["distance_km"] for r in records])
            # Use phase_alignment (TEP signature metric) instead of coherence (isotropic amplitude)
            coh = np.array([r["phase_alignment"] for r in records])
            az = np.array([r["azimuth_deg"] for r in records])

            def _run_null(label, dist_arr, coh_arr):
                fit = fit_exponential_model(
                    dist_arr, coh_arr, MIN_DISTANCE_KM, MAX_DISTANCE_KM, N_BINS, MIN_BIN_COUNT
                )
                if fit and fit.get("success"):
                    null_pass = fit["r_squared"] < NULL_TEST_R2_THRESHOLD
                    status = "PASS" if null_pass else "FAIL"
                    print_status(
                        f"    {label}: λ={fit['correlation_length_km']:.0f}km, R²={fit['r_squared']:.3f} [{status}]",
                        "INFO" if null_pass else "WARNING"
                    )
                else:
                    print_status(f"    {label}: no fit (null_pass=True)", "INFO")
                return fit

            # Temporal shuffle: randomize phase_alignment values independently of distance
            coh_temporal = rng.permutation(coh)
            fit_temporal = _run_null("Temporal shuffle", dist, coh_temporal)

            # Spatial shuffle: randomize distances independently of phase_alignment
            dist_spatial = rng.permutation(dist)
            fit_spatial = _run_null("Spatial shuffle", dist_spatial, coh)

            # Phase randomization: randomize the sign of each phase_alignment value.
            # This is equivalent to adding either 0 or π to each pair's phase angle,
            # which properly destroys coherent phase structure while preserving the
            # marginal amplitude distribution.
            coh_phase = coh * np.where(rng.random(len(coh)) < 0.5, -1, 1)
            fit_phase = _run_null("Phase randomization", dist, coh_phase)

            # Solar rotation test: rotate azimuth by 90° then shuffle values
            # between EW sectors (E+W) and all remaining pairs.  If the
            # anisotropy is tied to a fixed spatial frame, rotating the azimuth
            # labels and shuffling destroys the coherent structure.
            az_solar = (az + 90) % 360
            ew_mask_rot = _sector_mask(az_solar, 90) | _sector_mask(az_solar, 270)
            # All non-EW pairs (including diagonal sectors) receive the rest of
            # the shuffled values; this preserves total pair count.
            non_ew_mask_rot = ~ew_mask_rot
            n_ew_rot = np.sum(ew_mask_rot)
            n_non_ew_rot = np.sum(non_ew_mask_rot)
            coh_shuffled = coh.copy()
            rng.shuffle(coh_shuffled)
            coh_solar = np.empty_like(coh)
            coh_solar[ew_mask_rot] = coh_shuffled[:n_ew_rot]
            coh_solar[non_ew_mask_rot] = coh_shuffled[n_ew_rot:]
            fit_solar = _run_null("Solar rotation null", dist, coh_solar)

            def _extract(fit):
                return {
                    "lambda_km": fit["correlation_length_km"] if fit and fit.get("success") else None,
                    "R2": fit["r_squared"] if fit and fit.get("success") else None,
                    "null_pass": (fit is None or not fit.get("success") or fit["r_squared"] < NULL_TEST_R2_THRESHOLD)
                }

            n_pass = sum([
                _extract(fit_temporal)["null_pass"],
                _extract(fit_spatial)["null_pass"],
                _extract(fit_phase)["null_pass"],
                _extract(fit_solar)["null_pass"],
            ])
            results[const_name] = {
                "temporal_shuffle": _extract(fit_temporal),
                "spatial_shuffle": _extract(fit_spatial),
                "phase_random": _extract(fit_phase),
                "solar_rotation": _extract(fit_solar),
                "null_tests_passed": n_pass,
                "total_null_tests": 4,
                "status": "success",
                "n_pairs": len(records)
            }
            print_status(
                f"    {const_name}: {n_pass}/4 null tests passed",
                "SUCCESS" if n_pass == 4 else "WARNING"
            )

        out_file = OUTPUTS_DIR / "step_2_7_null_tests.json"
        with open(out_file, 'w') as f:
            json.dump(results, f, indent=2)
        print_status(f"Null test results written to {out_file}", "INFO")
        return {"status": "success", "results": str(out_file)}
