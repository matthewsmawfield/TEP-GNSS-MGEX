#!/usr/bin/env python3
"""TEP-GNSS-MGEX Analysis Pipeline Master Script
===============================================
Orchestrates the full held-out multi-constellation replication pipeline.

Usage:
    python scripts/run_all.py
    python scripts/run_all.py --start-step 1
    python scripts/run_all.py --skip-steps 3,4

Author: Matthew Lukin Smawfield
Date: May 2026
Version: v0.1-Suva
License: CC-BY-4.0
"""

import sys
import time
import json
import argparse
import os
from pathlib import Path
from datetime import datetime, timezone

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.logger import TEPLogger, set_step_logger, print_status, check_memory_usage
from scripts.utils.config import LOGS_DIR, OUTPUTS_DIR, FIGURES_DIR

# Import steps
from scripts.steps.step_1_0_data_acquisition import Step10DataAcquisition
from scripts.steps.step_1_1_station_metadata import Step11StationMetadata
from scripts.steps.step_2_0_mgex_clock_correlation import Step20MGEXClockCorrelation
from scripts.steps.step_2_1_correlation_length import Step21CorrelationLength
from scripts.steps.step_2_2_ew_ns_anisotropy import Step22EWNSAnisotropy
from scripts.steps.step_2_3_orbital_coupling import Step23OrbitalCoupling
from scripts.steps.step_2_4_cmb_alignment import Step24CMBAlignment
from scripts.steps.step_2_5_ionospheric_control import Step25IonosphericControl
from scripts.steps.step_2_6_geometry_control import Step26GeometryControl
from scripts.steps.step_2_8_satellite_clock_analysis import Step28SatelliteClockAnalysis
from scripts.steps.step_2_7_null_tests import Step27NullTests
from scripts.steps.step_3_0_cross_constellation import Step30CrossConstellation

STEP_REGISTRY = {
    0: ("1_0_data_acquisition", Step10DataAcquisition, "Data Acquisition (NASA CDDIS / IGS)"),
    1: ("1_1_station_metadata", Step11StationMetadata, "Station Metadata & Coordinates"),
    2: ("2_0_mgex", Step20MGEXClockCorrelation, "MGEX Clock Correlation Analysis (Multi-GNSS)"),
    3: ("2_1_correlation", Step21CorrelationLength, "Correlation Length λ_T Fitting"),
    4: ("2_2_anisotropy", Step22EWNSAnisotropy, "EW > NS Anisotropy Test"),
    5: ("2_3_orbital", Step23OrbitalCoupling, "Orbital-Velocity Coupling"),
    6: ("2_4_cmb", Step24CMBAlignment, "CMB-Frame Alignment Scan"),
    7: ("2_5_iono", Step25IonosphericControl, "Ionospheric Control (Kp / storm exclusion)"),
    8: ("2_6_geometry", Step26GeometryControl, "Geometry Control (hemisphere / distance matched)"),
    9: ("2_8_satellite", Step28SatelliteClockAnalysis, "Satellite-Clock Per-Constellation Analysis (SP3)"),
    10: ("2_7_null", Step27NullTests, "Null Control Tests"),
    11: ("3_0_cross", Step30CrossConstellation, "MGEX Validation Synthesis"),
}

# Map step numbers to their module names for logger lookup
STEP_MODULE_MAP = {
    0: "scripts.steps.step_1_0_data_acquisition",
    1: "scripts.steps.step_1_1_station_metadata",
    2: "scripts.steps.step_2_0_mgex_clock_correlation",
    3: "scripts.steps.step_2_1_correlation_length",
    4: "scripts.steps.step_2_2_ew_ns_anisotropy",
    5: "scripts.steps.step_2_3_orbital_coupling",
    6: "scripts.steps.step_2_4_cmb_alignment",
    7: "scripts.steps.step_2_5_ionospheric_control",
    8: "scripts.steps.step_2_6_geometry_control",
    9: "scripts.steps.step_2_8_satellite_clock_analysis",
    10: "scripts.steps.step_2_7_null_tests",
    11: "scripts.steps.step_3_0_cross_constellation",
}


def clear_pipeline_logs():
    """Clear all step log files and the pipeline master log.

    Must be called AFTER all step modules are imported so their loggers exist.
    Uses TEPLogger.clear_log() to properly close and reopen file handlers.
    """
    # Clear pipeline master log
    master_path = LOGS_DIR / "pipeline_master.log"
    try:
        with open(master_path, 'w') as f:
            f.write("")
    except Exception:
        pass

    # Clear each step's logger via its module-level TEPLogger instance
    for module_name in STEP_MODULE_MAP.values():
        if module_name in sys.modules:
            mod = sys.modules[module_name]
            logger = getattr(mod, 'logger', None)
            if logger and hasattr(logger, 'clear_log'):
                logger.clear_log()


def run_pipeline(args: argparse.Namespace) -> dict:
    LOGS_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    # Clear all log files (AFTER imports so step loggers exist)
    clear_pipeline_logs()

    level = "DEBUG" if args.verbose else "INFO"
    pipeline_logger = TEPLogger(
        "pipeline_master",
        log_file_path=LOGS_DIR / "pipeline_master.log",
        level=level
    )
    set_step_logger(pipeline_logger)

    print_status("=" * 70, "TITLE")
    print_status("TEP-GNSS-MGEX ANALYSIS PIPELINE", "TITLE")
    print_status("Paper 4: Held-Out Multi-Constellation Replication (Suva)", "TITLE")
    print_status("=" * 70, "TITLE")
    print_status(f"Project Root: {PROJECT_ROOT}", "INFO")
    print_status(f"Started: {datetime.now(timezone.utc).isoformat()}", "INFO")
    print_status(f"Python: {sys.version.split()[0]}", "INFO")
    print_status(f"Verbose: {args.verbose}", "INFO")
    print_status("")

    start_step = args.start_step if args.start_step is not None else 0
    stop_step = args.stop_step if args.stop_step is not None else 11
    skip_steps = set(args.skip_steps) if args.skip_steps else set()

    results = {
        "pipeline_start": datetime.now(timezone.utc).isoformat(),
        "steps_completed": [],
        "steps_failed": [],
        "steps_skipped": [],
        "execution_times": {},
        "status": "RUNNING"
    }

    total_start_time = time.time()

    for step_num in range(start_step, stop_step + 1):
        if step_num in skip_steps:
            print_status(f">>> STEP {step_num:02d}: SKIPPED", "WARNING")
            results["steps_skipped"].append(step_num)
            continue
        if step_num not in STEP_REGISTRY:
            print_status(f">>> STEP {step_num:02d}: NOT FOUND", "ERROR")
            results["steps_failed"].append(step_num)
            continue

        step_id, StepClass, step_desc = STEP_REGISTRY[step_num]
        print_status(f">>> STEP {step_num:02d}: {step_desc.upper()}", "TITLE")
        print_status("")
        check_memory_usage(f"before step {step_num}")

        # Redirect logging to this step's individual log file
        step_module_name = STEP_MODULE_MAP.get(step_num)
        step_module_logger = None
        if step_module_name and step_module_name in sys.modules:
            step_module = sys.modules[step_module_name]
            step_module_logger = getattr(step_module, "logger", None)
        if step_module_logger:
            set_step_logger(step_module_logger)

        step_start_time = time.time()
        try:
            step_instance = StepClass()
            step_result = step_instance.run()
            results["steps_completed"].append(step_num)
            step_time = time.time() - step_start_time
            results["execution_times"][step_id] = round(step_time, 2)
            print_status(f"\nStep {step_num:02d} completed in {step_time:.1f}s", "SUCCESS")
        except Exception as e:
            results["steps_failed"].append(step_num)
            step_time = time.time() - step_start_time
            results["execution_times"][step_id] = round(step_time, 2)
            print_status(f"\nStep {step_num:02d} failed after {step_time:.1f}s", "ERROR")
            print_status(f"Error: {e}", "ERROR")
            if not args.continue_on_error:
                print_status("\nPipeline halted. Use --continue-on-error to proceed.", "WARNING")
                results["status"] = "FAILED"
                break
        finally:
            # Restore pipeline master logger for orchestration messages
            set_step_logger(pipeline_logger)

        check_memory_usage(f"after step {step_num}")
        print_status("")
        print_status("-" * 70, "TITLE")
        print_status("")

    total_time = time.time() - total_start_time
    results["pipeline_end"] = datetime.now(timezone.utc).isoformat()
    results["total_time_seconds"] = round(total_time, 2)

    print_status("=" * 70, "TITLE")
    print_status("PIPELINE EXECUTION SUMMARY", "TITLE")
    print_status("=" * 70, "TITLE")
    print_status("")
    print_status(f"Completed: {len(results['steps_completed'])}", "INFO")
    print_status(f"Skipped:   {len(results['steps_skipped'])}", "INFO")
    print_status(f"Failed:    {len(results['steps_failed'])}", "INFO")
    print_status(f"Total time: {total_time:.1f}s", "INFO")
    print_status("")

    if results["steps_failed"]:
        results["status"] = "COMPLETED_WITH_ERRORS" if args.continue_on_error else "FAILED"
        print_status(f"Status: {results['status']}", "WARNING")
    else:
        results["status"] = "SUCCESS"
        print_status("Status: ALL STEPS COMPLETED SUCCESSFULLY", "SUCCESS")

    results_file = OUTPUTS_DIR / "pipeline_results.json"
    with open(results_file, 'w') as f:
        json.dump(results, f, indent=2)
    print_status(f"Pipeline results saved to: {results_file}", "INFO")
    return results


def main():
    parser = argparse.ArgumentParser(
        description="TEP-GNSS-MGEX Pipeline - Held-Out Multi-Constellation Replication",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Pipeline Steps:
  0  - Data Acquisition
  1  - Station Metadata
  2  - MGEX Clock Correlation
  3  - Correlation Length
  4  - EW/NS Anisotropy
  5  - Orbital-Velocity Coupling
  6  - CMB-Frame Alignment
  7  - Ionospheric Control
  8  - Geometry Control
  9  - Satellite-Clock Per-Constellation Analysis (SP3)
  10 - Null Tests
  11 - Cross-Constellation Synthesis
        """
    )
    parser.add_argument("--start-step", type=int, choices=range(0, 12), help="First step (0-11)")
    parser.add_argument("--stop-step", type=int, choices=range(0, 12), help="Last step (0-11)")
    parser.add_argument("--skip-steps", type=lambda s: [int(x) for x in s.split(",")], help="Comma-separated skip list")
    parser.add_argument("--continue-on-error", action="store_true", help="Continue on step failure")
    parser.add_argument("--list-steps", action="store_true", help="List steps and exit")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose (DEBUG) logging and per-step memory reporting")
    args = parser.parse_args()

    if args.list_steps:
        for num, (step_id, _, desc) in STEP_REGISTRY.items():
            print(f"  {num:02d}. {step_id:20s} - {desc}")
        return

    try:
        results = run_pipeline(args)
        if results["status"] == "SUCCESS":
            sys.exit(0)
        elif results["status"] == "COMPLETED_WITH_ERRORS":
            sys.exit(1)
        else:
            sys.exit(2)
    except KeyboardInterrupt:
        print("\n\nPipeline interrupted by user.")
        sys.exit(130)
    except Exception as e:
        print(f"\n\nPipeline crashed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(3)


if __name__ == "__main__":
    main()
