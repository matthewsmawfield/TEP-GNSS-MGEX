#!/usr/bin/env python3
"""Step 2.8: Per-constellation satellite-clock correlation analysis.

Uses MGEX SP3 orbit files (satellite positions) and CLK files (satellite clocks)
to compute within-constellation phase coherence as a function of satellite
separation.  Per-epoch common-mode removal is applied per constellation to
suppress shared system-time references.

Frozen prediction: each constellation independently shows a correlation length
in the thousands-of-km range when satellite clocks are common-mode removed.
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

import json
import gzip
import time
import numpy as np
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp

from scripts.utils.logger import TEPLogger, set_step_logger, print_status
from scripts.utils.config import (
    OUTPUTS_DIR, LOGS_DIR, RAW_DIR, DATA_START, DATA_END, CONSTELLATIONS,
    MIN_DISTANCE_KM, MAX_DISTANCE_KM, N_BINS, MIN_BIN_COUNT, FS_HZ, F1_HZ, F2_HZ,
    MAX_SATELLITE_DISTANCE_KM, MIN_EPOCHS, MIN_SATELLITE_PAIRS,
    MIN_CONSTELLATION_PAIRS
)
from scripts.utils.coherence import compute_coherence_phase, fit_exponential_model
from scripts.utils.geospatial import _safe_datetime_with_leap_second

LOGS_DIR.mkdir(parents=True, exist_ok=True)
logger = TEPLogger("step_2_8", log_file_path=LOGS_DIR / "step_2_8_satellite_clock.log")
set_step_logger(logger)

# ---------------------------------------------------------------------------
# SP3-c parser
# ---------------------------------------------------------------------------

def parse_sp3_file(sp3_path):
    """Parse an SP3-c file, returning {epoch_str: {sv_id: (x,y,z)}}.

    Epochs are formatted as 'YYYY-MM-DD HH:MM:SS' to match CLK epochs.
    Coordinates are in km.
    """
    positions = defaultdict(dict)
    if str(sp3_path).endswith('.gz'):
        opener = lambda: gzip.open(sp3_path, 'rt')
    else:
        opener = lambda: open(sp3_path, 'r')

    try:
        with opener() as f:
            lines = f.readlines()
    except Exception as e:
        print_status(f"  SP3 parse error {sp3_path.name}: {e}", "WARNING")
        return {}

    epoch_dt = None
    for line in lines:
        if not line.strip():
            continue
        if line[0] == '*':
            # Epoch line: *  YYYY  MM  DD  HH  MM  SS.SSSSSSSSS
            parts = line[1:].split()
            if len(parts) >= 6:
                try:
                    y, m, d = int(parts[0]), int(parts[1]), int(parts[2])
                    hr, mn = int(parts[3]), int(parts[4])
                    sec = float(parts[5])
                    epoch_dt = _safe_datetime_with_leap_second(y, m, d, hr, mn, sec)
                except (ValueError, IndexError):
                    epoch_dt = None
        elif line[0] == 'P' and epoch_dt is not None:
            # Position record: P{sv_id} {X} {Y} {Z} ...
            sv_id = line[1:4].strip()
            try:
                x = float(line[4:18].strip())
                y = float(line[18:32].strip())
                z = float(line[32:46].strip())
                # SP3 position units vary: standard is km, some files use 0.1 mm.
                # Use physical plausibility (satellite orbits ~7,000-42,000 km from
                # Earth's centre) to auto-detect and avoid mis-scaling.
                mag = np.sqrt(x*x + y*y + z*z)
                if mag > 5e7:  # 0.1 mm units (e.g. 2e8)
                    x, y, z = x / 1e7, y / 1e7, z / 1e7
                elif mag > 5e3:  # km units (standard, ~2e4)
                    pass
                else:  # Possibly metres or 1.0 mm
                    x, y, z = x / 1e3, y / 1e3, z / 1e3
                epoch_str = epoch_dt.strftime('%Y-%m-%d %H:%M:%S')
                positions[epoch_str][sv_id] = (x, y, z)
            except (ValueError, IndexError):
                continue
    return dict(positions)


# ---------------------------------------------------------------------------
# CLK AS (satellite-clock) parser
# ---------------------------------------------------------------------------

def parse_clk_as_records(clk_path):
    """Parse AS (satellite-clock) records from a CLK file.

    Returns: {sv_id: [(epoch_str, clock_offset_seconds)]}
    """
    records = defaultdict(list)
    opener = gzip.open if str(clk_path).endswith('.gz') else open
    mode = 'rt' if str(clk_path).endswith('.gz') else 'r'

    try:
        with opener(clk_path, mode) as f:
            for line in f:
                if not line.startswith('AS '):
                    continue
                parts = line.split()
                if len(parts) < 10:
                    continue
                sv_id = parts[1]
                try:
                    y, m, d = int(parts[2]), int(parts[3]), int(parts[4])
                    hr, mn = int(parts[5]), int(parts[6])
                    sec = float(parts[7])
                    # Index 8 is number-of-values (usually 2: bias + drift)
                    # Index 9 is clock bias in seconds
                    clock_s = float(parts[9])
                    epoch_dt = _safe_datetime_with_leap_second(y, m, d, hr, mn, sec)
                    epoch_str = epoch_dt.strftime('%Y-%m-%d %H:%M:%S')
                    records[sv_id].append((epoch_str, clock_s))
                except (ValueError, IndexError):
                    continue
    except Exception as e:
        print_status(f"  CLK AS parse error {clk_path.name}: {e}", "WARNING")

    return {sv: np.array(recs, dtype=object) for sv, recs in records.items()}


# ---------------------------------------------------------------------------
# Constellation mapping
# ---------------------------------------------------------------------------

CONSTELLATION_MAP = {
    'G': 'GPS', 'R': 'GLONASS', 'E': 'Galileo', 'C': 'BeiDou', 'J': 'QZSS'
}


def sv_to_constellation(sv_id):
    """Map satellite ID (e.g. 'G01', 'C06') to constellation name."""
    prefix = sv_id[0].upper()
    return CONSTELLATION_MAP.get(prefix)


# ---------------------------------------------------------------------------
# Per-day processing
# ---------------------------------------------------------------------------

def process_day(day_key, sp3_positions, clk_as_records):
    """Process one day: match SP3 positions with CLK clocks, remove common mode,
    compute satellite-pair phase coherence.

    Returns list of pair dicts with keys:
        sv1, sv2, constellation, epoch, separation_km, phase_alignment
    """
    pairs = []
    fs = FS_HZ
    f1, f2 = F1_HZ, F2_HZ

    # Group CLK records by constellation
    const_clocks = defaultdict(lambda: defaultdict(list))  # const -> epoch -> [(sv, clock_ns)]
    for sv_id, recs in clk_as_records.items():
        const = sv_to_constellation(sv_id)
        if not const:
            continue
        for epoch_str, clock_s in recs:
            const_clocks[const][epoch_str].append((sv_id, clock_s * 1e9))  # convert to ns

    # Remove per-epoch common mode and compute coherence for each constellation
    for const, epoch_data in const_clocks.items():
        # Need at least 2 satellites per epoch to form pairs
        valid_epochs = []
        for epoch_str, sv_clocks in epoch_data.items():
            if len(sv_clocks) < 2:
                continue
            # Check if we have SP3 positions for all these sats at this epoch
            sp3_epoch = sp3_positions.get(epoch_str, {})
            matched = []
            for sv_id, clock_ns in sv_clocks:
                if sv_id in sp3_epoch:
                    matched.append((sv_id, clock_ns, sp3_epoch[sv_id]))
            if len(matched) < 2:
                continue

            # Robust common-mode removal: subtract per-epoch median across
            # visible satellites to suppress shared system-time reference.
            clocks = np.array([c for _, c, _ in matched])
            median_clock = float(np.median(clocks))
            cleaned = [(sv, c - median_clock, pos) for sv, c, pos in matched]
            valid_epochs.append((epoch_str, cleaned))

        if len(valid_epochs) < 2:
            continue

        # Build time series per satellite (aligned epochs)
        all_epochs = sorted({ep for ep, _ in valid_epochs})
        epoch_to_idx = {ep: i for i, ep in enumerate(all_epochs)}
        sv_series = defaultdict(lambda: ([], []))  # sv -> (epoch_indices, clock_values)
        for epoch_str, cleaned in valid_epochs:
            ep_idx = epoch_to_idx[epoch_str]
            for sv, clock_ns, pos in cleaned:
                sv_series[sv][0].append(ep_idx)
                sv_series[sv][1].append(clock_ns)

        # Form pairs and compute coherence
        sv_list = list(sv_series.keys())
        for i in range(len(sv_list)):
            for j in range(i + 1, len(sv_list)):
                sv1, sv2 = sv_list[i], sv_list[j]
                idx1, vals1 = sv_series[sv1]
                idx2, vals2 = sv_series[sv2]

                # Find common epochs
                common_idx = np.intersect1d(idx1, idx2)
                if len(common_idx) < MIN_EPOCHS:  # Need enough points for CSD
                    continue

                idx_to_val1 = {k: v for k, v in zip(idx1, vals1)}
                idx_to_val2 = {k: v for k, v in zip(idx2, vals2)}
                v1 = np.array([idx_to_val1[k] for k in common_idx])
                v2 = np.array([idx_to_val2[k] for k in common_idx])

                # Compute phase coherence
                mean_coh, weighted_phase, phase_alignment = compute_coherence_phase(v1, v2, fs, f1, f2)
                if not np.isfinite(phase_alignment):
                    continue

                # Compute mean separation from SP3 positions across all common epochs
                separations = []
                for ep_idx in common_idx:
                    ep_str = all_epochs[ep_idx]
                    pos1 = sp3_positions[ep_str].get(sv1)
                    pos2 = sp3_positions[ep_str].get(sv2)
                    if pos1 is None or pos2 is None:
                        continue
                    sep = np.sqrt(
                        (pos1[0] - pos2[0])**2 +
                        (pos1[1] - pos2[1])**2 +
                        (pos1[2] - pos2[2])**2
                    )
                    separations.append(sep)
                if len(separations) < 1:
                    continue
                separation_km = float(np.median(separations))
                if separation_km < MIN_DISTANCE_KM or separation_km > MAX_SATELLITE_DISTANCE_KM:
                    continue

                pairs.append({
                    "sv1": sv1,
                    "sv2": sv2,
                    "constellation": const,
                    "day": day_key,
                    "separation_km": float(separation_km),
                    "phase_alignment": float(phase_alignment),
                    "n_epochs": int(len(common_idx)),
                })

    return pairs


# ---------------------------------------------------------------------------
# Main step class
# ---------------------------------------------------------------------------

class Step28SatelliteClockAnalysis:
    def run(self, max_workers=None):
        print_status("Step 2.8: Satellite-Clock Per-Constellation Analysis", "INFO")
        print_status("=" * 60, "INFO")

        if max_workers is None:
            max_workers = max(1, mp.cpu_count() - 1)

        # -------------------------------------------------------------------
        # Find available days (need both CLK and SP3)
        # -------------------------------------------------------------------
        clk_dir = RAW_DIR / "clk"
        sp3_dir = RAW_DIR / "sp3"

        if not clk_dir.exists() or not sp3_dir.exists():
            print_status("CLK or SP3 directories not found. Run data acquisition first.", "ERROR")
            return {"status": "failed", "reason": "missing_data"}

        # Map day_key -> (clk_path, sp3_path)
        day_files = {}
        for clk_path in sorted(clk_dir.glob("*.clk")):
            day_key = clk_path.stem.split("_")[0]
            # Find matching SP3 (same day, any AC)
            sp3_matches = list(sp3_dir.glob(f"{day_key}_*.sp3"))
            if sp3_matches:
                day_files[day_key] = (clk_path, sp3_matches[0])

        print_status(f"Days with both CLK and SP3: {len(day_files)}", "INFO")
        if len(day_files) < 10:
            print_status("Too few matched days for satellite analysis.", "ERROR")
            return {"status": "failed", "reason": "insufficient_matched_days"}

        # -------------------------------------------------------------------
        # Process days in parallel
        # -------------------------------------------------------------------
        print_status(f"Processing {len(day_files)} days with {max_workers} workers...", "INFO")

        all_pairs = []
        completed = 0
        t0 = time.time()

        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {}
            for day_key, (clk_path, sp3_path) in day_files.items():
                future = executor.submit(_worker_process_day, day_key, str(sp3_path), str(clk_path))
                futures[future] = day_key

            for future in as_completed(futures):
                day_key = futures[future]
                completed += 1
                try:
                    day_pairs = future.result()
                    all_pairs.extend(day_pairs)
                    if completed % 10 == 0 or completed == len(day_files):
                        elapsed = time.time() - t0
                        rate = completed / elapsed if elapsed > 0 else 0
                        eta = (len(day_files) - completed) / rate if rate > 0 else 0
                        print_status(
                            f"  [{completed}/{len(day_files)}] {day_key} | "
                            f"Pairs: {len(day_pairs)} | Rate: {rate:.1f} days/s | "
                            f"ETA: {eta/60:.1f} min",
                            "PROCESS"
                        )
                except Exception as e:
                    print_status(f"  Day {day_key} failed: {e}", "WARNING")

        print_status(f"Total satellite-pair records: {len(all_pairs)}", "INFO")

        if len(all_pairs) < MIN_SATELLITE_PAIRS:
            print_status("Insufficient satellite pairs for analysis.", "ERROR")
            return {"status": "failed", "reason": "insufficient_pairs"}

        # -------------------------------------------------------------------
        # Save raw pair records
        # -------------------------------------------------------------------
        pair_file = OUTPUTS_DIR / "step_2_8_satellite_pairs.json"
        with open(pair_file, 'w') as f:
            json.dump(all_pairs, f, indent=2)
        print_status(f"Satellite pair records written to {pair_file}", "INFO")

        # -------------------------------------------------------------------
        # Per-constellation exponential fit
        # -------------------------------------------------------------------
        const_results = {}
        for const in CONSTELLATIONS.keys():
            const_pairs = [p for p in all_pairs if p['constellation'] == const]
            if len(const_pairs) < MIN_CONSTELLATION_PAIRS:
                const_results[const] = {"status": "insufficient_pairs", "n_pairs": len(const_pairs)}
                continue

            sep = np.array([p['separation_km'] for p in const_pairs])
            pa = np.array([p['phase_alignment'] for p in const_pairs])
            valid = np.isfinite(sep) & np.isfinite(pa)
            sep, pa = sep[valid], pa[valid]

            if len(sep) < MIN_CONSTELLATION_PAIRS:
                const_results[const] = {"status": "insufficient_valid", "n_pairs": len(sep)}
                continue

            # Satellite separations can reach ~40 000 km (antipodal GEO/MEO),
            # so we use a larger upper bound than the ground-station default.
            fit = fit_exponential_model(sep, pa, MIN_DISTANCE_KM, MAX_SATELLITE_DISTANCE_KM, N_BINS, MIN_BIN_COUNT)
            if fit and fit.get('success'):
                # Strip heavy arrays from JSON
                fit_out = {k: v for k, v in fit.items() if k not in ('bin_centers', 'bin_means', 'bin_counts')}
                boundary_conv = bool(fit['correlation_length_km'] >= MAX_SATELLITE_DISTANCE_KM * 1.9)
                const_results[const] = {
                    "status": "success",
                    "lambda_km": fit['correlation_length_km'],
                    "lambda_err_km": fit.get('correlation_length_err_km', 0.0),
                    "amplitude": fit['amplitude'],
                    "offset": fit['offset'],
                    "r_squared": fit['r_squared'],
                    "n_pairs": int(len(sep)),
                    "boundary_convergence": boundary_conv,
                    "fit": fit_out,
                }
                print_status(
                    f"  {const}: λ = {fit['correlation_length_km']:.0f} ± "
                    f"{fit.get('correlation_length_err_km', 0):.0f} km, "
                    f"R² = {fit['r_squared']:.3f}, n = {len(sep)}",
                    "INFO"
                )
            else:
                const_results[const] = {"status": "fit_failed", "n_pairs": int(len(sep))}
                print_status(f"  {const}: Fit failed (n={len(sep)})", "WARNING")

        # -------------------------------------------------------------------
        # Save results
        # -------------------------------------------------------------------
        out_file = OUTPUTS_DIR / "step_2_8_satellite_clock_analysis.json"
        with open(out_file, 'w') as f:
            json.dump(const_results, f, indent=2)
        print_status(f"Satellite-clock analysis written to {out_file}", "INFO")

        return {"status": "success", "results": str(out_file), "n_pairs": len(all_pairs)}


# ---------------------------------------------------------------------------
# Worker function (must be top-level for pickling)
# ---------------------------------------------------------------------------

def _worker_process_day(day_key, sp3_path, clk_path):
    """Worker function for parallel day processing."""
    sp3_positions = parse_sp3_file(sp3_path)
    clk_as_records = parse_clk_as_records(clk_path)
    return process_day(day_key, sp3_positions, clk_as_records)


if __name__ == "__main__":
    step = Step28SatelliteClockAnalysis()
    step.run()
