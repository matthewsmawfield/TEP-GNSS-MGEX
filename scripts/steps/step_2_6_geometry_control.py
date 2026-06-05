#!/usr/bin/env python3
"""Step 2.6: Geometry control.

Frozen prediction: TEP signal is not network geometry.
Test: hemisphere-balanced and distance-matched station subsets.
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
    GEOMETRY_SAMPLE_SIZE, MIN_PAIRS_FIT, HEMISPHERE_STRATA_N, BOOTSTRAP_SEED
)
from scripts.utils.coherence import fit_exponential_model

LOGS_DIR.mkdir(parents=True, exist_ok=True)
logger = TEPLogger("step_2_6", log_file_path=LOGS_DIR / "step_2_6_geometry_control.log")
set_step_logger(logger)


class Step26GeometryControl:
    def run(self):
        print_status("Step 2.6: Geometry Control", "INFO")
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
                continue

            dist = np.array([r["distance_km"] for r in records])
            # Use phase_alignment (TEP signature metric) instead of coherence (isotropic amplitude)
            coh = np.array([r["phase_alignment"] for r in records])
            mid_lat = np.array([r["mid_lat"] for r in records])

            # Hemisphere split
            north_mask = mid_lat >= 0
            south_mask = mid_lat < 0
            print_status(
                f"  {const_name}: {np.sum(north_mask)} N pairs, {np.sum(south_mask)} S pairs",
                "INFO"
            )

            def _fit(mask):
                if np.sum(mask) < MIN_PAIRS_FIT:
                    return None
                return fit_exponential_model(
                    dist[mask], coh[mask],
                    MIN_DISTANCE_KM, MAX_DISTANCE_KM, N_BINS, MIN_BIN_COUNT
                )

            fit_north = _fit(north_mask)
            fit_south = _fit(south_mask)
            if fit_north and fit_north.get("success"):
                print_status(
                    f"    N(raw): λ={fit_north['correlation_length_km']:.0f}km, R²={fit_north['r_squared']:.3f}",
                    "INFO"
                )
            if fit_south and fit_south.get("success"):
                print_status(
                    f"    S(raw): λ={fit_south['correlation_length_km']:.0f}km, R²={fit_south['r_squared']:.3f}",
                    "INFO"
                )

            # Distance-stratified balanced hemisphere analysis
            # The IGS/MGEX network is N-biased (214 vs 119 stations).  To avoid
            # geometry artifacts we downsample the dominant hemisphere within each
            # distance decile so both hemispheres have identical distance distributions.
            north_idx = np.where(north_mask)[0]
            south_idx = np.where(south_mask)[0]

            n_strata = HEMISPHERE_STRATA_N
            dist_range = dist[north_mask | south_mask]
            if len(dist_range) > 0:
                stratum_edges = np.percentile(
                    dist_range, np.linspace(0, 100, n_strata + 1)
                )
            else:
                stratum_edges = np.linspace(MIN_DISTANCE_KM, MAX_DISTANCE_KM, n_strata + 1)

            rng = np.random.default_rng(BOOTSTRAP_SEED)
            balanced_north_idx = []
            balanced_south_idx = []
            for i in range(n_strata):
                low = stratum_edges[i]
                high = stratum_edges[i + 1]
                if i == n_strata - 1:
                    north_in = north_idx[(dist[north_idx] >= low) & (dist[north_idx] <= high)]
                    south_in = south_idx[(dist[south_idx] >= low) & (dist[south_idx] <= high)]
                else:
                    north_in = north_idx[(dist[north_idx] >= low) & (dist[north_idx] < high)]
                    south_in = south_idx[(dist[south_idx] >= low) & (dist[south_idx] < high)]
                n_keep = min(len(north_in), len(south_in))
                if n_keep > 0:
                    balanced_north_idx.extend(
                        rng.choice(north_in, size=n_keep, replace=False)
                    )
                    balanced_south_idx.extend(
                        rng.choice(south_in, size=n_keep, replace=False)
                    )

            print_status(
                f"  {const_name}: {len(balanced_north_idx)} balanced N, {len(balanced_south_idx)} balanced S pairs",
                "INFO"
            )

            balanced_north_mask = np.zeros(len(dist), dtype=bool)
            balanced_south_mask = np.zeros(len(dist), dtype=bool)
            balanced_north_mask[balanced_north_idx] = True
            balanced_south_mask[balanced_south_idx] = True

            fit_north_balanced = _fit(balanced_north_mask)
            fit_south_balanced = _fit(balanced_south_mask)
            if fit_north_balanced and fit_north_balanced.get("success"):
                print_status(
                    f"    N(bal): λ={fit_north_balanced['correlation_length_km']:.0f}km, R²={fit_north_balanced['r_squared']:.3f}",
                    "INFO"
                )
            if fit_south_balanced and fit_south_balanced.get("success"):
                print_status(
                    f"    S(bal): λ={fit_south_balanced['correlation_length_km']:.0f}km, R²={fit_south_balanced['r_squared']:.3f}",
                    "INFO"
                )

            # Distance-matched subset: sample a random subset preserving distance distribution
            # 50k points needed for stable fitting over 40 log-spaced bins; 2k was far too few.
            rng = np.random.default_rng(BOOTSTRAP_SEED)
            n_sample = min(GEOMETRY_SAMPLE_SIZE, len(dist))
            sample_idx = rng.choice(len(dist), size=n_sample, replace=False)
            fit_matched = fit_exponential_model(
                dist[sample_idx], coh[sample_idx],
                MIN_DISTANCE_KM, MAX_DISTANCE_KM, N_BINS, MIN_BIN_COUNT
            )
            if fit_matched and fit_matched.get("success"):
                print_status(
                    f"    Distance-matched (n={n_sample}): λ={fit_matched['correlation_length_km']:.0f}km, R²={fit_matched['r_squared']:.3f}",
                    "INFO"
                )

            results[const_name] = {
                "northern_hemisphere_lambda_km": fit_north["correlation_length_km"] if fit_north and fit_north["success"] else None,
                "northern_hemisphere_R2": fit_north["r_squared"] if fit_north and fit_north["success"] else None,
                "southern_hemisphere_lambda_km": fit_south["correlation_length_km"] if fit_south and fit_south["success"] else None,
                "southern_hemisphere_R2": fit_south["r_squared"] if fit_south and fit_south["success"] else None,
                "northern_hemisphere_balanced_lambda_km": fit_north_balanced["correlation_length_km"] if fit_north_balanced and fit_north_balanced["success"] else None,
                "northern_hemisphere_balanced_R2": fit_north_balanced["r_squared"] if fit_north_balanced and fit_north_balanced["success"] else None,
                "southern_hemisphere_balanced_lambda_km": fit_south_balanced["correlation_length_km"] if fit_south_balanced and fit_south_balanced["success"] else None,
                "southern_hemisphere_balanced_R2": fit_south_balanced["r_squared"] if fit_south_balanced and fit_south_balanced["success"] else None,
                "distance_matched_lambda_km": fit_matched["correlation_length_km"] if fit_matched and fit_matched["success"] else None,
                "distance_matched_R2": fit_matched["r_squared"] if fit_matched and fit_matched["success"] else None,
                "status": "success",
                "n_pairs": len(records)
            }

        out_file = OUTPUTS_DIR / "step_2_6_geometry_control.json"
        with open(out_file, 'w') as f:
            json.dump(results, f, indent=2)
        print_status(f"Geometry control results written to {out_file}", "INFO")
        return {"status": "success", "results": str(out_file)}
