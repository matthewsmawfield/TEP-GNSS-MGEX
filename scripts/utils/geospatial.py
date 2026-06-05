#!/usr/bin/env python3
"""Geospatial utilities for TEP-GNSS-MGEX."""

import numpy as np
from datetime import datetime, timedelta


def _safe_datetime_with_leap_second(y, m, d, hr, mn, sec):
    """Build a datetime robustly, handling leap seconds (sec >= 60).

    Leap-second epochs (sec == 60) are accepted by rolling forward into the
    next minute so that valid data are not silently lost.
    """
    sec_whole = int(sec)
    usec = int((sec % 1) * 1e6)
    dt = datetime(y, m, d, hr, mn, 0)
    dt += timedelta(seconds=sec_whole, microseconds=usec)
    return dt


def _sector_mask(azimuths, center, width=45):
    """Return boolean mask for azimuths within ±width/2 of center.

    Wraps correctly around 0°/360°.  Used consistently across anisotropy,
    orbital-coupling, CMB-alignment and null-test steps.
    """
    half = width / 2
    diff = np.abs((azimuths - center + 180) % 360 - 180)
    return diff <= half


def _earth_orbital_velocity_equatorial(doy):
    """Return Earth orbital velocity vector (km/s) in equatorial J2000 coordinates.

    Model: circular orbit, speed varies sinusoidally (perihelion ~day 3).
    """
    speed = 29.78 + 0.5 * np.cos(2 * np.pi * (doy - 3) / 365.25)
    # Ecliptic coords: x toward vernal equinox, y toward summer solstice
    prog = 2 * np.pi * (doy - 80) / 365.25
    # Earth is at ecliptic longitude (prog + π), not prog (Sun's longitude).
    # Velocity for CCW orbit at angle θ: v_x = -speed*sin(θ), v_y = speed*cos(θ)
    # With θ = prog + π: v_x = speed*sin(prog), v_y = -speed*cos(prog)
    v_ecl_x = speed * np.sin(prog)
    v_ecl_y = -speed * np.cos(prog)
    v_ecl_z = 0.0
    # Ecliptic -> Equatorial (J2000 obliquity 23.4367 deg)
    eps = np.radians(23.4367)
    v_eq_x = v_ecl_x
    v_eq_y = v_ecl_y * np.cos(eps) - v_ecl_z * np.sin(eps)
    v_eq_z = v_ecl_y * np.sin(eps) + v_ecl_z * np.cos(eps)
    return np.array([v_eq_x, v_eq_y, v_eq_z])


def ecef_to_lla(x, y, z):
    """Convert ECEF coordinates to latitude/longitude (WGS84)."""
    a = 6378137.0
    f = 1 / 298.257223563
    e2 = f * (2 - f)
    lon = np.arctan2(y, x)
    p = np.sqrt(x**2 + y**2)
    lat = np.arctan2(z, p * (1 - e2))
    for _ in range(5):
        N = a / np.sqrt(1 - e2 * np.sin(lat)**2)
        lat = np.arctan2(z + e2 * N * np.sin(lat), p)
    return np.degrees(lat), np.degrees(lon)


def haversine_km(lat1, lon1, lat2, lon2):
    """Great-circle distance between two points in kilometres."""
    R = 6371.0
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
    return R * 2 * np.arcsin(np.sqrt(np.clip(a, 0.0, 1.0)))


def calculate_azimuth(lat1, lon1, lat2, lon2):
    """Forward azimuth from point 1 to point 2 (degrees, 0° = north)."""
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlon = lon2 - lon1
    y = np.sin(dlon) * np.cos(lat2)
    x = np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(dlon)
    azimuth = np.degrees(np.arctan2(y, x))
    return (azimuth + 360) % 360
