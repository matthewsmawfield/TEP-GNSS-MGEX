#!/usr/bin/env python3
"""Step 2.1: Correlation length λ_T fitting.

Frozen prediction: λ_T within 1,000–4,000 km (MGEX-appropriate range).
Test: Fit C(r) = A exp(−r/λ) + C₀ on distance bins for each constellation.
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

import json
import numpy as np
from scripts.utils.logger import TEPLogger, set_step_logger, print_status
from scripts.utils.config import (
    OUTPUTS_DIR, FIGURES_DIR, LOGS_DIR, CONSTELLATIONS,
    MIN_DISTANCE_KM, MAX_DISTANCE_KM, N_BINS, MIN_BIN_COUNT,
    LAMBDA_MIN_KM, LAMBDA_MAX_KM
)
from scripts.utils.coherence import fit_exponential_model

LOGS_DIR.mkdir(parents=True, exist_ok=True)
logger = TEPLogger("step_2_1", log_file_path=LOGS_DIR / "step_2_1_correlation.log")
set_step_logger(logger)


class Step21CorrelationLength:
    def run(self):
        print_status("Step 2.1: Correlation Length Fitting", "INFO")
        results = {}

        pair_file = OUTPUTS_DIR / "step_2_0_mgex_pairs.json"
        if not pair_file.exists():
            print_status(f"Pair file not found: {pair_file}", "ERROR")
            results["combined"] = {"status": "no_data", "note": f"Pair file not found: {pair_file}"}
            out_file = OUTPUTS_DIR / "step_2_1_correlation_length.json"
            with open(out_file, 'w') as f:
                json.dump(results, f, indent=2)
            return {"status": "no_data", "results": str(out_file)}

        with open(pair_file) as f:
            records = json.load(f)
        if not records:
            print_status("Empty pair records", "WARNING")
            results["combined"] = {"status": "no_data", "note": "Empty pair records"}
            out_file = OUTPUTS_DIR / "step_2_1_correlation_length.json"
            with open(out_file, 'w') as f:
                json.dump(results, f, indent=2)
            return {"status": "no_data", "results": str(out_file)}

        distances = np.array([r["distance_km"] for r in records])
        # Use phase_alignment (cos of magnitude-weighted CSD phase) — the TEP signature metric.
        # coherence (mean spectral magnitude) is isotropic and lacks directional phase structure.
        coherences = np.array([r["phase_alignment"] for r in records])

        const_name = "combined"
        print_status(f"  Fitting {const_name}...", "INFO")

        fit = fit_exponential_model(
            distances, coherences,
            MIN_DISTANCE_KM, MAX_DISTANCE_KM, N_BINS, MIN_BIN_COUNT
        )

        if fit and fit["success"]:
            results[const_name] = {
                "lambda_km": fit["correlation_length_km"],
                "lambda_err_km": fit["correlation_length_err_km"],
                "A": fit["amplitude"],
                "C0": fit["offset"],
                "R2": fit["r_squared"],
                "status": "success",
                "in_range": bool(LAMBDA_MIN_KM <= fit["correlation_length_km"] <= LAMBDA_MAX_KM),
                "n_pairs": len(records)
            }
            print_status(f"    {const_name}: λ = {fit['correlation_length_km']:.0f} ± {fit['correlation_length_err_km']:.0f} km, R² = {fit['r_squared']:.3f}", "INFO")
        else:
            results[const_name] = {
                "status": "fit_failed",
                "note": "Exponential fit did not converge",
                "n_pairs": len(records)
            }

        out_file = OUTPUTS_DIR / "step_2_1_correlation_length.json"
        with open(out_file, 'w') as f:
            json.dump(results, f, indent=2)
        print_status(f"Correlation length results written to {out_file}", "INFO")
        return {"status": "success", "results": str(out_file)}
