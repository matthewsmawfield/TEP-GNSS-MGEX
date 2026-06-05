#!/usr/bin/env python3
"""TEP-GNSS-MGEX Space Weather Data Utilities

Lightweight Kp/Ap index fetcher from official sources:
- GFZ Potsdam (historical Kp, 1932-present)
- NOAA SWPC (recent Kp, last 30 days rolling)

Aligned with TEP-GNSS-II methodology but without pandas dependency.

Author: Matthew Lukin Smawfield
License: CC-BY-4.0
"""

import urllib.request
import urllib.error
import ssl
import json
from datetime import datetime, date
from pathlib import Path

# Anchor to project root for logger import
import sys
_PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_PROJECT_ROOT))

from scripts.utils.logger import print_status


def _day_key_to_date(day_key: str) -> date:
    """Convert YYYYDOY string to Python date."""
    year = int(day_key[:4])
    doy = int(day_key[4:])
    return date.fromordinal(date(year, 1, 1).toordinal() + doy - 1)


def _date_to_day_key(d: date) -> str:
    """Convert Python date to YYYYDOY string."""
    doy = (d - date(d.year, 1, 1)).days + 1
    return f"{d.year}{doy:03d}"


def fetch_kp_gfz(start_date: date, end_date: date) -> dict:
    """Fetch historical Kp/Ap from GFZ Potsdam. Returns {YYYYDOY: kp_float}."""
    url = "https://www-app3.gfz-potsdam.de/kp_index/Kp_ap_since_1932.txt"
    print_status(f"Fetching Kp from GFZ Potsdam ({start_date} to {end_date})...", "INFO")

    kp_map = {}
    try:
        ctx = ssl.create_default_context()
        with urllib.request.urlopen(url, context=ctx, timeout=60) as resp:
            content = resp.read().decode("utf-8").splitlines()

        # GFZ format: YYYY MM DD hh.h hh._m days days_m Kp ap D
        # Kp is at index 7 (0-based), ap at index 8. 8 intervals/day.
        # Missing data indicated by -1.000 for Kp.
        daily_kp = {}
        for line in content:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 10:
                continue
            try:
                year = int(parts[0])
                month = int(parts[1])
                day = int(parts[2])
                kp_val = float(parts[7])
                if kp_val < 0:  # missing data flag
                    continue
                dt = date(year, month, day)
                if start_date <= dt <= end_date:
                    dk = _date_to_day_key(dt)
                    daily_kp.setdefault(dk, []).append(kp_val)
            except (ValueError, IndexError):
                continue

        for dk, vals in daily_kp.items():
            kp_map[dk] = round(sum(vals) / len(vals), 1)

        print_status(f"GFZ Kp: {len(kp_map)} days fetched", "SUCCESS")
        return kp_map

    except Exception as e:
        print_status(f"GFZ Kp fetch failed: {e}", "WARNING")
        return {}


def fetch_kp_noaa(start_date: date, end_date: date) -> dict:
    """Fetch recent Kp from NOAA SWPC API. Returns {YYYYDOY: kp_float}."""
    api_url = "https://services.swpc.noaa.gov/json/planetary_k_index_1m.json"
    print_status(f"Fetching Kp from NOAA SWPC...", "INFO")

    kp_map = {}
    try:
        ctx = ssl.create_default_context()
        with urllib.request.urlopen(api_url, context=ctx, timeout=30) as resp:
            data = json.loads(resp.read().decode("utf-8"))

        # NOAA returns 3-hourly Kp values; aggregate to daily mean
        daily_kp = {}
        for record in data:
            try:
                date_str = record["time_tag"][:10]  # YYYY-MM-DD
                dt = datetime.strptime(date_str, "%Y-%m-%d").date()
                if start_date <= dt <= end_date:
                    kp_val = float(record.get("kp_index", 2.0))
                    day_key = _date_to_day_key(dt)
                    daily_kp.setdefault(day_key, []).append(kp_val)
            except (KeyError, ValueError, TypeError):
                continue

        for day_key, vals in daily_kp.items():
            kp_map[day_key] = round(sum(vals) / len(vals), 1)

        print_status(f"NOAA Kp: {len(kp_map)} days fetched", "SUCCESS")
        return kp_map

    except Exception as e:
        print_status(f"NOAA Kp fetch failed: {e}", "WARNING")
        return {}


def get_kp_data_for_range(start_date: date, end_date: date) -> dict:
    """Fetch Kp data for a date range, trying GFZ first then NOAA.

    Returns dict: {YYYYDOY string: Kp float value}
    """
    kp_map = fetch_kp_gfz(start_date, end_date)
    if not kp_map:
        kp_map = fetch_kp_noaa(start_date, end_date)
    return kp_map
