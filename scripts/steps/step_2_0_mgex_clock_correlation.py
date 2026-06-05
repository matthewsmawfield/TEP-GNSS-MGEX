#!/usr/bin/env python3
"""Step 2.0: MGEX multi-GNSS receiver clock correlation analysis.

Parses MGEX CLK files (AR records) to extract per-station receiver clock
offsets, builds daily time series, and computes phase-coherent pair
coherences for all valid station pairs.

MGEX CLK files contain a single multi-GNSS solution (GPS+GLO+GAL+BDS).
The same pair table is replicated for each constellation name in the
pipeline to maintain downstream compatibility.
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

import json
import os
import time
import numpy as np
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from scipy.signal import detrend
from scripts.utils.logger import TEPLogger, set_step_logger, print_status
from scripts.utils.config import (
    RAW_DIR, PROCESSED_DIR, OUTPUTS_DIR, FIGURES_DIR, LOGS_DIR,
    CONSTELLATIONS, F1_HZ, F2_HZ, FS_HZ,
    MIN_DISTANCE_KM, MAX_DISTANCE_KM, MIN_EPOCHS,
    MAX_STATION_STD_NS, MIN_STATION_STD_NS
)
from scripts.utils.geospatial import haversine_km, calculate_azimuth, ecef_to_lla, _safe_datetime_with_leap_second
from scripts.utils.coherence import compute_coherence_phase

LOGS_DIR.mkdir(parents=True, exist_ok=True)
logger = TEPLogger("step_2_0", log_file_path=LOGS_DIR / "step_2_0_mgex_clock.log")
set_step_logger(logger)


def _parse_clk_ar_records(clk_path):
    """Parse AR (receiver clock) records from a CLK file.

    Format: AR STATION YYYY MM DD HH MM SS.SSSSSS N CLOCK_OFFSET SIGMA
    Returns DataFrame-like dict: {station: [(timestamp, offset_sec), ...]}
    """
    records = defaultdict(list)
    try:
        with open(clk_path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if not line.startswith("AR "):
                    continue
                parts = line.split()
                if len(parts) < 10:
                    continue
                try:
                    station = parts[1].upper()
                    year, month, day = int(parts[2]), int(parts[3]), int(parts[4])
                    hour, minute = int(parts[5]), int(parts[6])
                    second = float(parts[7])
                    clock_offset = float(parts[9])  # seconds
                    ts = _safe_datetime_with_leap_second(year, month, day, hour, minute, second)
                    records[station].append((ts, clock_offset))
                except (ValueError, IndexError):
                    continue
    except Exception as e:
        logger.warning(f"Error reading {clk_path.name}: {e}")
    return records


class Step20MGEXClockCorrelation:
    def run(self):
        print_status("Step 2.0: MGEX Clock Correlation Analysis [OPTIMIZED v2]", "INFO")
        OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)
        FIGURES_DIR.mkdir(parents=True, exist_ok=True)
        PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

        # -------------------------------------------------------------------
        # Load station coordinates
        # -------------------------------------------------------------------
        coords_file = PROCESSED_DIR / "station_coordinates.json"
        station_coords = {}
        if coords_file.exists():
            with open(coords_file) as f:
                ecef_coords = json.load(f)
            for sta, xyz in ecef_coords.items():
                if len(xyz) == 3:
                    lat, lon = ecef_to_lla(xyz[0], xyz[1], xyz[2])
                    station_coords[sta] = {"lat": lat, "lon": lon}
            print_status(f"  Loaded {len(station_coords)} station coordinates", "INFO")
        else:
            print_status("  Warning: station_coordinates.json not found", "WARNING")
            return {"status": "failed", "reason": "missing_coordinates"}

        # -------------------------------------------------------------------
        # Parse CLK files
        # -------------------------------------------------------------------
        clk_dir = RAW_DIR / "clk"
        if not clk_dir.exists():
            print_status("  CLK directory not found", "ERROR")
            return {"status": "failed", "reason": "no_clk_files"}

        clk_files = sorted(clk_dir.glob("*.clk"))
        if not clk_files:
            print_status("  No CLK files found", "ERROR")
            return {"status": "failed", "reason": "no_clk_files"}

        print_status(f"  Parsing {len(clk_files)} CLK files...", "INFO")

        # Aggregate records by day -> station
        day_station_data = defaultdict(lambda: defaultdict(list))
        for i, clk_path in enumerate(clk_files, 1):
            records = _parse_clk_ar_records(clk_path)
            for station, entries in records.items():
                if not entries:
                    continue
                # Group by actual day (YYYYDOY) in case a file spans midnight
                for e in entries:
                    day_key = e[0].strftime("%Y%j")
                    day_station_data[day_key][station].append(e)
            if i % 50 == 0 or i == len(clk_files):
                print_status(f"    Parsed {i}/{len(clk_files)} CLK files | {len(day_station_data)} days so far", "PROCESS")

        print_status(f"  Found {len(day_station_data)} days with data", "INFO")

        # -------------------------------------------------------------------
        # Save per-station daily NPZ (downstream compatibility)
        # -------------------------------------------------------------------
        npz_dir = PROCESSED_DIR / "mgex"
        npz_dir.mkdir(parents=True, exist_ok=True)

        total_npz = 0
        for day_key, station_entries in day_station_data.items():
            for station, entries in station_entries.items():
                if len(entries) < MIN_EPOCHS:
                    continue
                # Sort by timestamp
                entries_sorted = sorted(entries, key=lambda x: x[0])
                timestamps = np.array([e[0].timestamp() for e in entries_sorted])
                # Clock offset in seconds -> nanoseconds
                clock_bias_ns = np.array([e[1] * 1e9 for e in entries_sorted])
                npz_path = npz_dir / f"{station}_{day_key}.npz"
                np.savez(npz_path,
                         timestamps=timestamps,
                         clock_bias_ns=clock_bias_ns,
                         pos_jitter=np.zeros_like(timestamps))
                total_npz += 1

        print_status(f"  Saved {total_npz} station-day NPZ files", "INFO")

        # -------------------------------------------------------------------
        # Compute pair coherences
        # -------------------------------------------------------------------
        pair_records = self._compute_pairs(day_station_data, station_coords)

        # Save the same pair table for each constellation entry
        # MGEX is a single multi-GNSS solution; write one pair file for all constellations
        table_file = OUTPUTS_DIR / "step_2_0_mgex_pairs.json"
        with open(table_file, "w") as f:
            json.dump(pair_records, f, indent=2)

        out_file = OUTPUTS_DIR / "step_2_0_pair_coherences.json"
        n_stations = len(station_coords)
        with open(out_file, "w") as f:
            json.dump({c: {"status": "success", "n_pairs": len(pair_records), "n_stations": n_stations}
                       for c in CONSTELLATIONS}, f, indent=2)

        print_status(f"Pair coherence: {len(pair_records)} records written ({n_stations} stations)", "INFO")
        return {"status": "success", "results": str(out_file), "n_pairs": len(pair_records), "n_stations": n_stations}

    def _compute_pairs(self, day_station_data, station_coords, max_workers=None):
        """Compute phase-coherent pair coherences across all days (parallelized).

        Uses a single ThreadPoolExecutor with day-level parallelism for efficiency,
        avoiding the overhead of creating 474 separate executors.
        """
        if max_workers is None:
            max_workers = min(os.cpu_count() or 4, 16)

        # Pre-build all day tasks
        day_tasks = []
        for day_key, station_entries in day_station_data.items():
            valid_stations = [
                sta for sta in station_entries
                if sta in station_coords and len(station_entries[sta]) >= MIN_EPOCHS
            ]
            if len(valid_stations) < 2:
                continue

            # Fast path: shared timestamps within a day
            first_ts_raw = [e[0].timestamp() for e in station_entries[valid_stations[0]]]
            first_ts = set(first_ts_raw)
            if len(first_ts) < len(first_ts_raw):
                logger.warning(
                    f"Duplicate timestamps in {day_key} station {valid_stations[0]}: "
                    f"{len(first_ts_raw) - len(first_ts)} duplicates dropped"
                )
            all_match = all(
                {e[0].timestamp() for e in station_entries[sta]} == first_ts
                for sta in valid_stations[1:]
            )

            if all_match:
                common_ts = sorted(first_ts)
                sta_series = {}
                for sta in valid_stations:
                    entry_dict = {e[0].timestamp(): e[1] * 1e9 for e in station_entries[sta]}
                    sta_series[sta] = np.array([entry_dict[ts] for ts in common_ts])
            else:
                timestamps_by_sta = {
                    sta: np.array(sorted({e[0].timestamp() for e in station_entries[sta]}))
                    for sta in valid_stations
                }
                common_ts = timestamps_by_sta[valid_stations[0]]
                for sta in valid_stations[1:]:
                    common_ts = np.intersect1d(common_ts, timestamps_by_sta[sta])
                    if len(common_ts) < MIN_EPOCHS:
                        break
                if len(common_ts) < MIN_EPOCHS:
                    continue
                sta_series = {}
                for sta in valid_stations:
                    entry_dict = {e[0].timestamp(): e[1] * 1e9 for e in station_entries[sta]}
                    sta_series[sta] = np.array([entry_dict[ts] for ts in common_ts])

            if len(sta_series) < 2:
                continue

            # Robust common-mode removal:
            # 1. Detrend each station to remove receiver-specific drift.
            # 2. Filter out stations with excessive residual variation
            #    (malformed clocks, missing data blocks, etc.).
            # 3. Filter out flat-line stations with negligible variation
            #    (stable clocks held constant by AC; CSD phases are pure noise).
            # 4. Subtract per-epoch median across good stations to remove
            #    shared AC time-scale reference without outlier contamination.
            # Use centralized thresholds from config.py

            sta_series_detrended = {}
            outlier_stations = []
            flatline_stations = []
            for sta, series in sta_series.items():
                detrended = detrend(series, type='linear')
                std = float(np.std(detrended))
                if std > MAX_STATION_STD_NS:
                    outlier_stations.append((sta, std))
                elif std < MIN_STATION_STD_NS:
                    flatline_stations.append((sta, std))
                else:
                    sta_series_detrended[sta] = detrended

            if outlier_stations:
                logger.warning(
                    f"{day_key}: filtered {len(outlier_stations)} outlier stations "
                    f"(>{MAX_STATION_STD_NS:.0f} ns std): "
                    + ", ".join(f"{s}({std:.0f}ns)" for s, std in outlier_stations[:5])
                    + ("..." if len(outlier_stations) > 5 else "")
                )
            if flatline_stations:
                logger.warning(
                    f"{day_key}: filtered {len(flatline_stations)} flat-line stations "
                    f"(<{MIN_STATION_STD_NS:.0f} ns std): "
                    + ", ".join(f"{s}({std:.3f}ns)" for s, std in flatline_stations[:5])
                    + ("..." if len(flatline_stations) > 5 else "")
                )

            sta_series = sta_series_detrended
            if len(sta_series) < 2:
                continue

            stacked = np.vstack([sta_series[sta] for sta in sta_series])
            epoch_median = np.median(stacked, axis=0)
            for sta in sta_series:
                sta_series[sta] = sta_series[sta] - epoch_median

            pair_args = []
            stations_list = list(sta_series.keys())
            for i in range(len(stations_list)):
                for j in range(i + 1, len(stations_list)):
                    s1, s2 = stations_list[i], stations_list[j]
                    lat1, lon1 = station_coords[s1]["lat"], station_coords[s1]["lon"]
                    lat2, lon2 = station_coords[s2]["lat"], station_coords[s2]["lon"]
                    dist = float(haversine_km(lat1, lon1, lat2, lon2))
                    if dist < MIN_DISTANCE_KM or dist > MAX_DISTANCE_KM:
                        continue
                    az = float(calculate_azimuth(lat1, lon1, lat2, lon2))
                    mid_lat = float((lat1 + lat2) / 2)
                    pair_args.append((s1, s2, dist, az, mid_lat, sta_series[s1], sta_series[s2], day_key))

            if pair_args:
                day_tasks.append((day_key, pair_args))

        # Single executor: process entire days in parallel across workers
        pair_records = []
        t0 = time.time()
        completed_days = 0
        print_status(f"  Computing pairs for {len(day_tasks)} days using {max_workers} workers...", "INFO")

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_day = {
                executor.submit(self._compute_day_pairs, pair_args, day_key): (day_key, len(pair_args))
                for day_key, pair_args in day_tasks
            }

            for future in as_completed(future_to_day):
                day_key, n_expected = future_to_day[future]
                day_results = future.result()
                pair_records.extend(day_results)
                completed_days += 1
                if completed_days % 10 == 0 or completed_days == len(day_tasks):
                    elapsed = time.time() - t0
                    rate = completed_days / elapsed if elapsed > 0 else 0
                    eta = (len(day_tasks) - completed_days) / rate if rate > 0 else 0
                    print_status(
                        f"  [{completed_days}/{len(day_tasks)}] Day {day_key}: "
                        f"{len(day_results)}/{n_expected} pairs | "
                        f"Rate: {rate:.1f} days/s | ETA: {eta/60:.1f} min",
                        "PROCESS"
                    )

        print_status(f"  Total pairs computed: {len(pair_records)} in {time.time()-t0:.1f}s", "INFO")
        return pair_records

    def _compute_day_pairs(self, pair_args, day_key):
        """Compute all pairs for a single day (runs inside executor worker)."""
        day_results = []
        for args in pair_args:
            result = self._compute_one_pair(args)
            if result:
                day_results.append(result)
        return day_results

    def _compute_one_pair(self, args):
        """Compute coherence for a single station pair."""
        s1, s2, dist, az, mid_lat, v1, v2, day_key = args
        mean_coh, weighted_phase, phase_alignment = compute_coherence_phase(
            v1, v2, FS_HZ, F1_HZ, F2_HZ
        )
        if np.isfinite(mean_coh) and np.isfinite(phase_alignment):
            return {
                "sta1": s1, "sta2": s2,
                "distance_km": dist,
                "azimuth_deg": az,
                "mid_lat": mid_lat,
                "metric": "clock_bias_ns",
                "coherence": float(mean_coh),
                "phase_alignment": float(phase_alignment),
                "day": day_key,
                "n_points": len(v1)
            }
        return None
