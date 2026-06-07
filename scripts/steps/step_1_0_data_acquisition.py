#!/usr/bin/env python3
"""Step 1.0: Download MGEX multi-GNSS receiver clock (CLK) and orbit (SP3) files.

Authentication:
    ~/.netrc entry:
        machine urs.earthdata.nasa.gov
            login <your-username>
            password <your-password>
    Or environment variables:
        CDDIS_USER=<your-username>
        CDDIS_PASS=<your-password>

CLK Sources (CDDIS, priority order):
    1. COD0MGXFIN — CODE MGEX final (GPS+GLO+GAL+BDS+QZSS)
    2. WUM0MGXFIN — Wuhan MGEX final
    3. GRG0MGXFIN — CNES MGEX final
    4. COD0OPSFIN — CODE OPS final (fallback)

SP3 Sources (CODE MGEX orbits, same AC as CLK):
    1. COD0MGXFIN — CODE MGEX final orbits
    2. WUM0MGXFIN — Wuhan MGEX final orbits (fallback)

File formats:
    CLK: {AC}0MGXFIN_YYYYDDD0000_01D_30S_CLK.CLK.gz  (~3–13 MB/day)
    SP3: {AC}0MGXFIN_YYYYDDD0000_01D_05M_ORB.SP3.gz   (~1–3 MB/day)
    Location: https://cddis.nasa.gov/archive/gnss/products/{gps_week}/

Output structure:
    data/raw/clk/{year}{doy:03d}_{ac}.clk   -> Decompressed CLK files
    data/raw/sp3/{year}{doy:03d}_{ac}.sp3  -> Decompressed SP3 files
    data/external/data_provenance.json       -> Download metadata
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

import os
import json
import gzip
import hashlib
import threading
import time as _time
from datetime import datetime, timedelta, timezone
from collections import defaultdict

import requests
from concurrent.futures import ThreadPoolExecutor, as_completed
from scripts.utils.logger import TEPLogger, set_step_logger, print_status, check_memory_usage
from scripts.utils.config import RAW_DIR, LOGS_DIR, EXTERNAL_DIR, DATA_START, DATA_END

LOGS_DIR.mkdir(parents=True, exist_ok=True)
logger = TEPLogger("step_1_0", log_file_path=LOGS_DIR / "step_1_0_data_acquisition.log")
set_step_logger(logger)

# Thread-safe lock for concurrent hash-registry writes
_hash_lock = threading.Lock()

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
CDDIS_BASE = "https://cddis.nasa.gov/archive/gnss/products"

# Analysis window, configured in config.py (DATA_START/DATA_END).
# Deliberately disjoint from the original TEP-GNSS-RINEX training period.
START_DATE = datetime(*DATA_START)
END_DATE = datetime(*DATA_END)

# AC priority: best multi-GNSS coverage first
AC_PRIORITY = [
    {"ac": "COD", "name": "CODE MGEX", "pattern": "{ac}0MGXFIN_{year}{doy:03d}0000_01D_30S_CLK.CLK.gz"},
    {"ac": "WUM", "name": "Wuhan MGEX", "pattern": "{ac}0MGXFIN_{year}{doy:03d}0000_01D_30S_CLK.CLK.gz"},
    {"ac": "GRG", "name": "CNES MGEX", "pattern": "{ac}0MGXFIN_{year}{doy:03d}0000_01D_30S_CLK.CLK.gz"},
    {"ac": "COD", "name": "CODE OPS", "pattern": "{ac}0OPSFIN_{year}{doy:03d}0000_01D_30S_CLK.CLK.gz"},
]

MAX_DOWNLOAD_WORKERS = 4
MAX_RETRIES = 3
RETRY_DELAY_SEC = 5


# ---------------------------------------------------------------------------
# Authentication
# ---------------------------------------------------------------------------
def get_auth():
    """Get NASA Earthdata credentials from ~/.netrc or environment."""
    try:
        import netrc
        auth = netrc.netrc().authenticators("urs.earthdata.nasa.gov")
        if auth:
            return (auth[0], auth[2])
    except Exception:
        pass
    user = os.getenv("CDDIS_USER")
    passwd = os.getenv("CDDIS_PASS")
    if user and passwd:
        return (user, passwd)
    return None


# ---------------------------------------------------------------------------
# GPS week / day-of-week calculation
# ---------------------------------------------------------------------------
def date_to_gps_week_dow(date):
    """Convert datetime to GPS week and day-of-week (0=Sunday)."""
    gps_epoch = datetime(1980, 1, 6)
    days = (date - gps_epoch).days
    week = days // 7
    dow = days % 7
    return week, dow


# ---------------------------------------------------------------------------
# MGEX CLK file download
# ---------------------------------------------------------------------------
def _compute_file_hash(path):
    """Compute SHA-256 hex digest of a file."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def _verify_hash(path, expected_hash):
    """Verify file SHA-256 matches expected_hash. Returns bool."""
    if not expected_hash:
        return False
    return _compute_file_hash(path) == expected_hash


def download_clk_file(year, doy, auth):
    """Download one day's MGEX CLK file, trying ACs in priority order.

    Returns: Path to decompressed .clk file, or None if all ACs failed.
    """
    clk_dir = RAW_DIR / "clk"
    clk_dir.mkdir(parents=True, exist_ok=True)

    date = datetime(year, 1, 1) + timedelta(days=doy - 1)
    gps_week, dow = date_to_gps_week_dow(date)

    # Load existing hash registry for integrity verification on re-runs
    hash_registry_path = EXTERNAL_DIR / "file_hashes.json"
    hash_registry = {}
    if hash_registry_path.exists():
        try:
            with open(hash_registry_path) as f:
                hash_registry = json.load(f)
        except Exception:
            pass

    # Check if any AC file already exists for this day and verify hash
    for ac_info in AC_PRIORITY:
        ac = ac_info["ac"]
        candidate = clk_dir / f"{year}{doy:03d}_{ac}.clk"
        if candidate.exists() and candidate.stat().st_size > 1000:
            key = str(candidate.relative_to(PROJECT_ROOT))
            expected = hash_registry.get(key, {}).get("sha256")
            if expected and _verify_hash(candidate, expected):
                return candidate
            elif not expected:
                # First run: register hash silently
                return candidate
            else:
                # Hash mismatch: corrupted or modified; re-download
                candidate.unlink(missing_ok=True)
                logger.warning(f"Hash mismatch for {candidate.name}; re-downloading")

    session = requests.Session()
    session.auth = auth

    # Try each AC in priority order
    for ac_info in AC_PRIORITY:
        ac = ac_info["ac"]
        filename = ac_info["pattern"].format(ac=ac, year=year, doy=doy)
        output_path = clk_dir / f"{year}{doy:03d}_{ac}.clk"
        compressed_path = clk_dir / f"{year}{doy:03d}_{ac}.clk.gz"

        # URL candidates: standard week dir, then _IGS20 suffix dir
        week_dirs = [str(gps_week), f"{gps_week}_IGS20"]
        urls = [f"{CDDIS_BASE}/{w}/{filename}" for w in week_dirs]

        for url in urls:
            for attempt in range(MAX_RETRIES):
                try:
                    resp = session.get(url, timeout=60, stream=True)
                    if resp.status_code != 200:
                        continue

                    # Download
                    with open(compressed_path, "wb") as f:
                        for chunk in resp.iter_content(chunk_size=8192):
                            if chunk:
                                f.write(chunk)

                    # Decompress
                    with open(compressed_path, "rb") as gz_f:
                        content = gzip.decompress(gz_f.read())
                    with open(output_path, "wb") as f:
                        f.write(content)

                    compressed_path.unlink(missing_ok=True)

                    # Validate decompressed header to catch truncated/corrupt files
                    valid_header = False
                    if output_path.exists() and output_path.stat().st_size > 1000:
                        try:
                            with open(output_path, "rb") as f:
                                header = f.read(300).decode("ascii", errors="ignore")
                            # CLK files are RINEX CLK; expect version/type header
                            if "RINEX VERSION" in header or "PGM / RUN BY / DATE" in header:
                                valid_header = True
                        except Exception:
                            pass
                    if valid_header:
                        # Register SHA-256 for provenance and re-run verification
                        file_hash = _compute_file_hash(output_path)
                        key = str(output_path.relative_to(PROJECT_ROOT))
                        with _hash_lock:
                            # Re-read fresh copy to avoid clobbering other threads
                            fresh_registry = {}
                            if hash_registry_path.exists():
                                try:
                                    with open(hash_registry_path) as f:
                                        fresh_registry = json.load(f)
                                except Exception:
                                    pass
                            fresh_registry[key] = {
                                "sha256": file_hash,
                                "downloaded_at": datetime.now(timezone.utc).isoformat(),
                                "size_bytes": output_path.stat().st_size,
                            }
                            with open(hash_registry_path, "w") as f:
                                json.dump(fresh_registry, f, indent=2)
                        return output_path
                    else:
                        output_path.unlink(missing_ok=True)
                        raise Exception("decompressed file empty or invalid header")

                except requests.RequestException:
                    if attempt < MAX_RETRIES - 1:
                        _time.sleep(RETRY_DELAY_SEC)
                    continue
                except Exception:
                    if attempt < MAX_RETRIES - 1:
                        _time.sleep(RETRY_DELAY_SEC)
                    continue

    return None


# ---------------------------------------------------------------------------
# SP3 orbit file download
# ---------------------------------------------------------------------------
SP3_AC_PRIORITY = [
    {"ac": "COD", "name": "CODE MGEX", "pattern": "{ac}0MGXFIN_{year}{doy:03d}0000_01D_05M_ORB.SP3.gz"},
    {"ac": "WUM", "name": "Wuhan MGEX", "pattern": "{ac}0MGXFIN_{year}{doy:03d}0000_01D_05M_ORB.SP3.gz"},
    {"ac": "COD", "name": "CODE MGEX (15m)", "pattern": "{ac}0MGXFIN_{year}{doy:03d}0000_01D_15M_ORB.SP3.gz"},
]


def download_sp3_file(year, doy, auth):
    """Download one day's MGEX SP3 orbit file, trying ACs in priority order.

    Returns: Path to decompressed .sp3 file, or None if all ACs failed.
    """
    sp3_dir = RAW_DIR / "sp3"
    sp3_dir.mkdir(parents=True, exist_ok=True)

    date = datetime(year, 1, 1) + timedelta(days=doy - 1)
    gps_week, dow = date_to_gps_week_dow(date)

    # Load existing hash registry for integrity verification on re-runs
    hash_registry_path = EXTERNAL_DIR / "file_hashes.json"
    hash_registry = {}
    if hash_registry_path.exists():
        try:
            with open(hash_registry_path) as f:
                hash_registry = json.load(f)
        except Exception:
            pass

    # Check if any AC file already exists for this day and verify hash
    for ac_info in SP3_AC_PRIORITY:
        ac = ac_info["ac"]
        candidate = sp3_dir / f"{year}{doy:03d}_{ac}.sp3"
        if candidate.exists() and candidate.stat().st_size > 1000:
            key = str(candidate.relative_to(PROJECT_ROOT))
            expected = hash_registry.get(key, {}).get("sha256")
            if expected and _verify_hash(candidate, expected):
                return candidate
            elif not expected:
                return candidate
            else:
                candidate.unlink(missing_ok=True)
                logger.warning(f"Hash mismatch for {candidate.name}; re-downloading")

    session = requests.Session()
    session.auth = auth

    for ac_info in SP3_AC_PRIORITY:
        ac = ac_info["ac"]
        filename = ac_info["pattern"].format(ac=ac, year=year, doy=doy)
        output_path = sp3_dir / f"{year}{doy:03d}_{ac}.sp3"
        compressed_path = sp3_dir / f"{year}{doy:03d}_{ac}.sp3.gz"

        week_dirs = [str(gps_week), f"{gps_week}_IGS20"]
        urls = [f"{CDDIS_BASE}/{w}/{filename}" for w in week_dirs]

        for url in urls:
            for attempt in range(MAX_RETRIES):
                try:
                    resp = session.get(url, timeout=60, stream=True)
                    if resp.status_code != 200:
                        continue
                    with open(compressed_path, "wb") as f:
                        for chunk in resp.iter_content(chunk_size=8192):
                            if chunk:
                                f.write(chunk)
                    with open(compressed_path, "rb") as gz_f:
                        content = gzip.decompress(gz_f.read())
                    with open(output_path, "wb") as f:
                        f.write(content)
                    compressed_path.unlink(missing_ok=True)
                    # Validate decompressed header to catch truncated/corrupt files
                    valid_header = False
                    if output_path.exists() and output_path.stat().st_size > 1000:
                        try:
                            with open(output_path, "rb") as f:
                                header = f.read(100).decode("ascii", errors="ignore")
                            # SP3 files start with #c or #a (format identifier)
                            if header.startswith("#") or "SP3" in header.upper():
                                valid_header = True
                        except Exception:
                            pass
                    if valid_header:
                        # Register SHA-256 for provenance and re-run verification
                        file_hash = _compute_file_hash(output_path)
                        key = str(output_path.relative_to(PROJECT_ROOT))
                        with _hash_lock:
                            # Re-read fresh copy to avoid clobbering other threads
                            fresh_registry = {}
                            if hash_registry_path.exists():
                                try:
                                    with open(hash_registry_path) as f:
                                        fresh_registry = json.load(f)
                                except Exception:
                                    pass
                            fresh_registry[key] = {
                                "sha256": file_hash,
                                "downloaded_at": datetime.now(timezone.utc).isoformat(),
                                "size_bytes": output_path.stat().st_size,
                            }
                            with open(hash_registry_path, "w") as f:
                                json.dump(fresh_registry, f, indent=2)
                        return output_path
                    else:
                        output_path.unlink(missing_ok=True)
                        raise Exception("decompressed file empty or invalid header")
                except requests.RequestException:
                    if attempt < MAX_RETRIES - 1:
                        _time.sleep(RETRY_DELAY_SEC)
                    continue
                except Exception:
                    if attempt < MAX_RETRIES - 1:
                        _time.sleep(RETRY_DELAY_SEC)
                    continue
    return None


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------
class Step10DataAcquisition:
    def run(self, max_workers=MAX_DOWNLOAD_WORKERS):
        print_status("Step 1.0: MGEX CLK Data Acquisition (NASA CDDIS)", "INFO")
        print_status("=" * 60, "INFO")

        auth = get_auth()
        if not auth:
            msg = (
                "No CDDIS credentials found.\n"
                "  Option 1: Add to ~/.netrc:\n"
                "    machine urs.earthdata.nasa.gov\n"
                "      login <username>\n"
                "      password <password>\n"
                "  Option 2: Set env vars CDDIS_USER and CDDIS_PASS\n"
                "Register at: https://urs.earthdata.nasa.gov/"
            )
            print_status(msg, "ERROR")
            return {"status": "failed", "reason": "missing_credentials"}

        dates = []
        current = START_DATE
        while current <= END_DATE:
            dates.append(current)
            current += timedelta(days=1)

        print_status(f"Date range: {START_DATE.date()} to {END_DATE.date()} ({len(dates)} days)", "INFO")
        print_status(f"CDDIS user: {auth[0]}", "INFO")
        print_status(f"Priority ACs: COD (MGEX) -> WUM -> GRG -> COD (OPS)", "INFO")

        RAW_DIR.mkdir(parents=True, exist_ok=True)
        (RAW_DIR / "clk").mkdir(parents=True, exist_ok=True)
        (RAW_DIR / "sp3").mkdir(parents=True, exist_ok=True)
        EXTERNAL_DIR.mkdir(parents=True, exist_ok=True)

        # -------------------------------------------------------------------
        # Download CLK files
        # -------------------------------------------------------------------
        clk_results = defaultdict(lambda: {"downloaded": 0, "failed": 0, "by_ac": defaultdict(int)})
        print_status("Downloading MGEX CLK files...", "INFO")
        print_status(f"  Workers: {max_workers} | Retries: {MAX_RETRIES} | Total days: {len(dates)}", "INFO")

        completed = 0
        clk_success = 0
        clk_fail = 0
        t0 = _time.time()

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {}
            for dt in dates:
                year, doy = dt.year, dt.timetuple().tm_yday
                future = executor.submit(download_clk_file, year, doy, auth)
                futures[future] = (year, doy, "clk")

            for future in as_completed(futures):
                year, doy, ftype = futures[future]
                completed += 1
                try:
                    clk_path = future.result()
                    if clk_path:
                        ac_code = clk_path.stem.split("_")[-1]
                        clk_results[year]["downloaded"] += 1
                        clk_results[year]["by_ac"][ac_code] += 1
                        clk_success += 1
                        if completed % 10 == 0 or completed == len(dates):
                            elapsed = _time.time() - t0
                            rate = completed / elapsed if elapsed > 0 else 0
                            eta = (len(dates) - completed) / rate if rate > 0 else 0
                            print_status(
                                f"  [{completed}/{len(dates)}] {year}/{doy:03d} CLK OK ({ac_code}) | "
                                f"Rate: {rate:.1f} files/s | ETA: {eta/60:.1f} min",
                                "PROCESS"
                            )
                    else:
                        clk_results[year]["failed"] += 1
                        clk_fail += 1
                        if completed % 10 == 0 or completed == len(dates):
                            print_status(
                                f"  [{completed}/{len(dates)}] {year}/{doy:03d} CLK FAILED | "
                                f"Success: {clk_success} | Failed: {clk_fail}",
                                "WARNING"
                            )
                except Exception as e:
                    clk_results[year]["failed"] += 1
                    clk_fail += 1
                    logger.warning(f"CLK download error {year}/{doy:03d}: {e}")

        total_clk_ok = sum(v["downloaded"] for v in clk_results.values())
        total_clk_fail = sum(v["failed"] for v in clk_results.values())
        clk_ac_breakdown = defaultdict(int)
        for yr_data in clk_results.values():
            for ac, cnt in yr_data["by_ac"].items():
                clk_ac_breakdown[ac] += cnt

        elapsed_clk = _time.time() - t0
        print_status(f"  CLK files: {total_clk_ok} downloaded, {total_clk_fail} failed", "INFO")
        print_status(f"  Elapsed: {elapsed_clk/60:.1f} min", "INFO")
        for ac, cnt in sorted(clk_ac_breakdown.items()):
            print_status(f"    {ac}: {cnt} days", "INFO")

        # -------------------------------------------------------------------
        # Download SP3 files
        # -------------------------------------------------------------------
        sp3_results = defaultdict(lambda: {"downloaded": 0, "failed": 0, "by_ac": defaultdict(int)})
        print_status("Downloading MGEX SP3 orbit files...", "INFO")
        print_status(f"  Workers: {max_workers} | Retries: {MAX_RETRIES} | Total days: {len(dates)}", "INFO")

        completed_sp3 = 0
        sp3_success = 0
        sp3_fail = 0
        t0_sp3 = _time.time()

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {}
            for dt in dates:
                year, doy = dt.year, dt.timetuple().tm_yday
                future = executor.submit(download_sp3_file, year, doy, auth)
                futures[future] = (year, doy, "sp3")

            for future in as_completed(futures):
                year, doy, ftype = futures[future]
                completed_sp3 += 1
                try:
                    sp3_path = future.result()
                    if sp3_path:
                        ac_code = sp3_path.stem.split("_")[-1]
                        sp3_results[year]["downloaded"] += 1
                        sp3_results[year]["by_ac"][ac_code] += 1
                        sp3_success += 1
                        if completed_sp3 % 10 == 0 or completed_sp3 == len(dates):
                            elapsed = _time.time() - t0_sp3
                            rate = completed_sp3 / elapsed if elapsed > 0 else 0
                            eta = (len(dates) - completed_sp3) / rate if rate > 0 else 0
                            print_status(
                                f"  [{completed_sp3}/{len(dates)}] {year}/{doy:03d} SP3 OK ({ac_code}) | "
                                f"Rate: {rate:.1f} files/s | ETA: {eta/60:.1f} min",
                                "PROCESS"
                            )
                    else:
                        sp3_results[year]["failed"] += 1
                        sp3_fail += 1
                        if completed_sp3 % 10 == 0 or completed_sp3 == len(dates):
                            print_status(
                                f"  [{completed_sp3}/{len(dates)}] {year}/{doy:03d} SP3 FAILED | "
                                f"Success: {sp3_success} | Failed: {sp3_fail}",
                                "WARNING"
                            )
                except Exception as e:
                    sp3_results[year]["failed"] += 1
                    sp3_fail += 1
                    logger.warning(f"SP3 download error {year}/{doy:03d}: {e}")

        total_sp3_ok = sum(v["downloaded"] for v in sp3_results.values())
        total_sp3_fail = sum(v["failed"] for v in sp3_results.values())
        sp3_ac_breakdown = defaultdict(int)
        for yr_data in sp3_results.values():
            for ac, cnt in yr_data["by_ac"].items():
                sp3_ac_breakdown[ac] += cnt

        elapsed_sp3 = _time.time() - t0_sp3
        print_status(f"  SP3 files: {total_sp3_ok} downloaded, {total_sp3_fail} failed", "INFO")
        print_status(f"  Elapsed: {elapsed_sp3/60:.1f} min", "INFO")
        for ac, cnt in sorted(sp3_ac_breakdown.items()):
            print_status(f"    {ac}: {cnt} days", "INFO")

        check_memory_usage("after data acquisition")

        # -------------------------------------------------------------------
        # Provenance
        # -------------------------------------------------------------------
        provenance = {
            "source": "NASA CDDIS / IGS MGEX",
            "base_url": CDDIS_BASE,
            "start_date": START_DATE.isoformat(),
            "end_date": END_DATE.isoformat(),
            "n_days_requested": len(dates),
            "clk_ac_priority": [a["name"] for a in AC_PRIORITY],
            "sp3_ac_priority": [a["name"] for a in SP3_AC_PRIORITY],
            "auth_method": "netrc" if (Path.home() / ".netrc").exists() else "env",
            "auth_user": auth[0],
            "clk_results": {
                str(y): {
                    "downloaded": v["downloaded"],
                    "failed": v["failed"],
                    "by_ac": dict(v["by_ac"]),
                }
                for y, v in clk_results.items()
            },
            "sp3_results": {
                str(y): {
                    "downloaded": v["downloaded"],
                    "failed": v["failed"],
                    "by_ac": dict(v["by_ac"]),
                }
                for y, v in sp3_results.items()
            },
            "clk_ac_breakdown": dict(clk_ac_breakdown),
            "sp3_ac_breakdown": dict(sp3_ac_breakdown),
            "hash_registry": str(EXTERNAL_DIR / "file_hashes.json"),
            "downloaded_at": datetime.now(timezone.utc).isoformat(),
        }
        prov_file = EXTERNAL_DIR / "data_provenance.json"
        with open(prov_file, "w") as f:
            json.dump(provenance, f, indent=2)
        print_status(f"Provenance written to {prov_file}", "INFO")
        print_status(f"SHA-256 hashes written to {EXTERNAL_DIR / 'file_hashes.json'}", "INFO")

        print_status("=" * 60, "INFO")
        print_status("Data acquisition complete.", "INFO")
        return {
            "status": "success",
            "days_requested": len(dates),
            "clk_files": total_clk_ok,
            "clk_failed": total_clk_fail,
            "sp3_files": total_sp3_ok,
            "sp3_failed": total_sp3_fail,
            "clk_ac_breakdown": dict(clk_ac_breakdown),
            "sp3_ac_breakdown": dict(sp3_ac_breakdown),
            "provenance": str(prov_file),
        }


if __name__ == "__main__":
    step = Step10DataAcquisition()
    step.run()
