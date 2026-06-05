#!/usr/bin/env python3
"""Step 1.1: Extract station coordinates from MGEX CLK file headers.

CLK files embed precise ECEF (XYZ) coordinates in the RINEX header section
under the "SOLN STA NAME / NUM" record. These are extracted once and cached
as station_coordinates.json for downstream pair-distance computation.
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

import json
from scripts.utils.logger import TEPLogger, set_step_logger, print_status
from scripts.utils.config import RAW_DIR, PROCESSED_DIR, LOGS_DIR, COORD_MM_THRESHOLD_METRES

LOGS_DIR.mkdir(parents=True, exist_ok=True)
logger = TEPLogger("step_1_1", log_file_path=LOGS_DIR / "step_1_1_station_metadata.log")
set_step_logger(logger)


def _parse_clk_header_coords(clk_path):
    """Parse station ECEF coordinates from a single CLK file header.

    Returns dict: {station_name: [x_m, y_m, z_m], ...}
    """
    coords = {}
    try:
        with open(clk_path, "r", encoding="utf-8", errors="ignore") as f:
            in_header = True
            for line in f:
                if not in_header:
                    break
                if "END OF HEADER" in line:
                    in_header = False
                    break
                if "SOLN STA NAME / NUM" in line:
                    content = line.split("SOLN STA NAME / NUM")[0]
                    parts = content.split()
                    if len(parts) >= 5:
                        name = parts[0].strip().upper()
                        try:
                            x = float(parts[2])
                            y = float(parts[3])
                            z = float(parts[4])
                            # CLK headers store coords in millimetres; convert to metres
                            mag = (x**2 + y**2 + z**2) ** 0.5
                            if mag > COORD_MM_THRESHOLD_METRES:  # > 100 Mm → units are mm
                                x /= 1000.0
                                y /= 1000.0
                                z /= 1000.0
                            coords[name] = [round(x, 4), round(y, 4), round(z, 4)]
                        except ValueError:
                            continue
    except Exception as e:
        logger.warning(f"Error reading {clk_path.name}: {e}")
    return coords


class Step11StationMetadata:
    def run(self):
        print_status("Step 1.1: Station Metadata from CLK Headers", "INFO")
        PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

        clk_dir = RAW_DIR / "clk"
        if not clk_dir.exists():
            print_status(f"CLK directory not found: {clk_dir}", "ERROR")
            return {"status": "failed", "reason": "no_clk_files"}

        clk_files = sorted(clk_dir.glob("*.clk"))
        if not clk_files:
            print_status("No CLK files found", "ERROR")
            return {"status": "failed", "reason": "no_clk_files"}

        print_status(f"Scanning {len(clk_files)} CLK files for coordinates...", "INFO")

        # Scan in reverse chronological order so first-seen (latest) coords are kept
        stations = {}
        for i, clk_path in enumerate(reversed(clk_files)):
            file_coords = _parse_clk_header_coords(clk_path)
            for name, xyz in file_coords.items():
                if name not in stations:
                    stations[name] = xyz
            if (i + 1) % 50 == 0 or i == 0:
                print_status(f"  Scanned {i+1}/{len(clk_files)} files, {len(stations)} stations so far", "INFO")

        print_status(f"Extracted coordinates for {len(stations)} unique stations", "INFO")

        # Save JSON for downstream
        coords_file = PROCESSED_DIR / "station_coordinates.json"
        with open(coords_file, "w") as f:
            json.dump(stations, f, indent=2)
        print_status(f"Coordinates written to {coords_file}", "INFO")

        return {
            "status": "success",
            "n_stations": len(stations),
            "coordinates_file": str(coords_file),
        }
