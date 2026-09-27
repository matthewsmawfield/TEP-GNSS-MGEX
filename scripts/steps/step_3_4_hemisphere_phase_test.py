#!/usr/bin/env python3
"""
TEP-GNSS-MGEX Analysis - STEP 3.4: Hemisphere-Split Phase Test
==============================================================

Closes issue 14-8 (corpus issues register): the distance-decay structure
is episodic (two of four quarterly windows return R^2 <= 0.02), and the
question is whether the modulation of the coherent channel is driven by a
global heliocentric mechanism (Earth orbital velocity, TEP) or by
hemisphere-local seasonal drivers (insolation, ionospheric and
tropospheric cycles).

The two classes of driver predict opposite hemisphere behaviour:

  * Global kinematic driver: the same annual modulation phase appears in
    the Northern- and Southern-hemisphere sub-networks (in-phase), and
    each hemisphere's monthly series correlates with Earth's orbital
    velocity with the same sign.
  * Local-seasonal driver: insolation reverses sign between hemispheres,
    so a meteorological or ionospheric seasonal driver produces
    ANTI-phased annual modulation (phase difference ~180 deg), or at
    minimum no common phase.

Design (per hemisphere class, using true station latitudes):

  1. Every pair is classified intra-NH (both endpoints lat >= 0) or
     intra-SH (both endpoints lat < 0).  Cross-hemisphere pairs are
     reported as a reference class.
  2. Per chronological calendar month, two observables are computed per
     class:
       a. the EW - NS phase-alignment contrast (pa_diff) -- the channel
          step 2.3 found modulated by the orbital-velocity projection
          (combined r = -0.77);
       b. the mean phase alignment (coherence amplitude).
  3. An annual sinusoid is fitted to each hemisphere's monthly series;
     the phase difference between hemispheres is the discriminator.
  4. Each hemisphere's monthly series is regressed against Earth's
     orbital speed and against the CMB-frame velocity projection, with
     sign consistency required for a single global driver.

Uncertainties: sinusoid coefficient errors are propagated to a phase
uncertainty by the delta method; the NH-SH phase difference carries the
two phase errors in quadrature.  The month-level NH vs SH correlation of
the pa_diff series is reported with a circular-shift permutation p-value.

Inputs:
    results/outputs/step_2_0_mgex_pairs.json    pair-day records
    data/processed/station_coordinates.json     station ECEF (m)

Output:
    results/outputs/step_3_4_hemisphere_phase_test.json
"""

import gc
import json
import math
import sys
from datetime import datetime
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np

from scripts.utils.logger import TEPLogger, set_step_logger, print_status
from scripts.utils.config import (
    OUTPUTS_DIR, LOGS_DIR,
    MIN_MONTHLY_RECORDS, MIN_DAYS_PER_MONTH,
)
from scripts.utils.geospatial import _sector_mask, _earth_orbital_velocity_equatorial

LOGS_DIR.mkdir(parents=True, exist_ok=True)
logger = TEPLogger("step_3_4", log_file_path=LOGS_DIR / "step_3_4_hemisphere_phase_test.log")
set_step_logger(logger)

SEED = 20260927
N_PERM = 20000


def ecef_to_lat_deg(xyz):
    x, y, z = xyz
    r = (x * y + y * y + z * z) ** 0.5 if False else (x * x + y * y + z * z) ** 0.5
    return float(np.degrees(np.arcsin(z / r))) if r > 0 else 0.0


def month_from_day_key(day_key):
    """YYYYDOY -> chronological 'YYYY-MM' calendar month key."""
    try:
        year = int(day_key[:4])
        doy = int(day_key[4:])
        if doy < 1 or doy > 366:
            return None
        d = datetime(year, 1, 1).toordinal() + doy - 1
        dd = datetime.fromordinal(d)
        return f"{dd.year:04d}-{dd.month:02d}"
    except Exception:
        return None


def month_mid_ordinal(month_key):
    y, m = map(int, month_key.split("-"))
    return datetime(y, m, 15).toordinal()


def month_mid_doy(month_key):
    y, m = map(int, month_key.split("-"))
    return datetime(y, m, 15).timetuple().tm_yday


def fit_annual_sinusoid(t_ord, y):
    """OLS fit y ~ 1 + sin(2 pi t / 365.25) + cos(2 pi t / 365.25).

    Returns amplitude, phase (deg, convention: phase of maximum
    elongation day-of-year = phase_deg/360*365.25), coefficient SEs and an
    F-test p-value against the constant model.
    """
    t = np.asarray(t_ord, float)
    y = np.asarray(y, float)
    w = 2 * np.pi * t / 365.25
    X = np.column_stack([np.ones_like(t), np.sin(w), np.cos(w)])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = max(len(y) - 3, 1)
    sigma2 = float(resid @ resid / dof)
    cov = sigma2 * np.linalg.pinv(X.T @ X)
    se = np.sqrt(np.maximum(np.diag(cov), 0.0))
    amp = float(math.hypot(beta[1], beta[2]))
    phase = float(math.atan2(beta[1], beta[2]))
    phase_deg = (math.degrees(phase) + 360.0) % 360.0
    # delta-method phase SE: phase = atan2(b1, b2); dphase ~ |db| / amp
    phase_se = float(math.sqrt(max(se[1] ** 2 + se[2] ** 2, 0.0)) / amp) if amp > 0 else float("nan")
    amp_se = float(se[1]) if len(se) > 1 else float("nan")
    rss1 = float(resid @ resid)
    rss0 = float(((y - y.mean()) ** 2).sum())
    F = ((rss0 - rss1) / 2.0) / (rss1 / dof) if rss1 > 0 else float("inf")
    try:
        from scipy.stats import f as _f
        p = float(_f.sf(F, 2, dof))
    except Exception:
        p = float("nan")
    return {
        "amplitude": amp,
        "amplitude_se": amp_se,
        "phase_deg": phase_deg,
        "phase_se_deg": math.degrees(phase_se) if math.isfinite(phase_se) else float("nan"),
        "peak_doy": phase_deg / 360.0 * 365.25,
        "F": float(F),
        "p_value": p,
        "n_months": int(len(y)),
    }


def regress_on_driver(month_keys, y, cmb_ra=167.0, cmb_dec=-7.0):
    """Correlate the monthly series with orbital speed and CMB-projected
    orbital velocity, with a circular-shift permutation p-value."""
    y = np.asarray(y, float)
    n = len(y)
    doys = np.array([month_mid_doy(m) for m in month_keys])
    speeds = np.array([29.78 + 0.5 * np.cos(2 * np.pi * (d - 3) / 365.25) for d in doys])
    ra, dec = np.radians(cmb_ra), np.radians(cmb_dec)
    cmb = np.array([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])
    vproj = np.array([float(np.dot(_earth_orbital_velocity_equatorial(d), cmb)) for d in doys])

    def _corr(x):
        sx = (x - x.mean()) / (x.std() or 1.0)
        sy = (y - y.mean()) / (y.std() or 1.0)
        return float(np.mean(sx * sy))

    def _perm_p(x):
        """Circular-shift null over ALL n-1 distinct rotations (identical
        convention to step 2.3).  With n monthly points the null is
        discrete, so the attainable floor is 1/n."""
        r_obs = abs(_corr(x))
        cnt = 0
        for shift in range(1, n):
            sx = (np.roll(x, shift) - x.mean()) / (x.std() or 1.0)
            sy = (y - y.mean()) / (y.std() or 1.0)
            if abs(float(np.mean(sx * sy))) >= r_obs:
                cnt += 1
        return (cnt + 1) / n

    out = {
        "r_vs_orbital_speed": _corr(speeds),
        "perm_p_vs_orbital_speed": _perm_p(speeds),
        "r_vs_cmb_vproj": _corr(vproj),
        "perm_p_vs_cmb_vproj": _perm_p(vproj),
        "n_months": n,
        "perm_method": "all n-1 circular shifts; discrete floor 1/n",
    }
    # parametric Pearson p as a complementary estimate
    try:
        from scipy.stats import t as _t
        import math as _m
        for key, x in (("speed", speeds), ("cmb_vproj", vproj)):
            r = abs(_corr(x))
            if n > 2 and r < 1:
                tstat = r * _m.sqrt((n - 2) / max(1 - r * r, 1e-12))
                out[f"pearson_p_vs_{key}"] = float(2 * _t.sf(tstat, n - 2))
    except Exception:
        pass
    return out


def monthly_series(month_key_arr, day_ord_arr, value, mask):
    """Aggregate per-chronological-month mean of value[mask].

    month_key_arr holds the 'YYYY-MM' key for each record; day_ord_arr the
    day ordinal (used to enforce a minimum spread of observing days per
    month).  Returns (month_keys, means, ns, day_counts)."""
    agg = {}
    mk = month_key_arr[mask]
    do = day_ord_arr[mask]
    va = value[mask]
    for m, o, v in zip(mk, do, va):
        if m is None:
            continue
        cell = agg.get(m)
        if cell is None:
            cell = agg[m] = [0.0, 0, set()]
        cell[0] += v; cell[1] += 1; cell[2].add(o)
    month_keys, means, ns, dc = [], [], [], []
    for m in sorted(agg):
        s, c, days = agg[m]
        if c < MIN_MONTHLY_RECORDS or len(days) < MIN_DAYS_PER_MONTH:
            continue
        month_keys.append(m)
        means.append(s / c)
        ns.append(c)
        dc.append(len(days))
    return month_keys, np.array(means), np.array(ns), np.array(dc)


class Step34HemispherePhaseTest:
    def run(self):
        return _run_impl()


def _run_impl():
    print_status("Step 3.4: Hemisphere-Split Phase Test", "INFO")
    results = {"step": "3.4", "name": "hemisphere_split_phase_test"}

    pair_file = OUTPUTS_DIR / "step_2_0_mgex_pairs.json"
    coord_file = PROJECT_ROOT / "data" / "processed" / "station_coordinates.json"
    if not pair_file.exists() or not coord_file.exists():
        print_status("Missing pair records or station coordinates", "ERROR")
        results["status"] = "no_data"
        out = OUTPUTS_DIR / "step_3_4_hemisphere_phase_test.json"
        json.dump(results, open(out, "w"), indent=2)
        return {"status": "no_data"}

    coords_raw = json.load(open(coord_file))
    lat = {k: ecef_to_lat_deg(v) for k, v in coords_raw.items()}
    print_status(f"  Station latitudes: {len(lat)} "
                 f"(NH {sum(1 for v in lat.values() if v >= 0)}, "
                 f"SH {sum(1 for v in lat.values() if v < 0)})", "INFO")

    records = json.load(open(pair_file))
    n = len(records)
    print_status(f"  Loaded {n} pair-day records", "INFO")

    pa = np.empty(n)
    az = np.empty(n)
    lat1 = np.empty(n)
    lat2 = np.empty(n)
    day_ord = np.empty(n, dtype=np.int64)
    month_key = np.empty(n, dtype=object)
    missing = 0
    for i, r in enumerate(records):
        l1 = lat.get(r["sta1"]); l2 = lat.get(r["sta2"])
        if l1 is None or l2 is None:
            missing += 1
            l1 = l2 = np.nan
        pa[i] = r["phase_alignment"]
        az[i] = r.get("azimuth_deg", 0.0)
        lat1[i] = l1; lat2[i] = l2
        dk = r.get("day", "")
        try:
            year = int(dk[:4]); doy = int(dk[4:])
            o = datetime(year, 1, 1).toordinal() + doy - 1 if 1 <= doy <= 366 else -1
        except Exception:
            o = -1
        day_ord[i] = o
        month_key[i] = month_from_day_key(dk)
    del records
    gc.collect()

    ok = ~np.isnan(lat1) & (day_ord > 0)
    pa, az, lat1, lat2, day_ord, month_key = (
        pa[ok], az[ok], lat1[ok], lat2[ok], day_ord[ok], month_key[ok])

    cls_nh = (lat1 >= 0) & (lat2 >= 0)
    cls_sh = (lat1 < 0) & (lat2 < 0)
    cls_x = ~(cls_nh | cls_sh)
    print_status(f"  Classes: intra-NH {int(cls_nh.sum())}, intra-SH {int(cls_sh.sum())}, "
                 f"cross-hemisphere {int(cls_x.sum())}", "INFO")

    # EW / NS azimuth sectors (paired sectors at 0/180 and 90/270 deg)
    ew = _sector_mask(az, 90) | _sector_mask(az, 270)
    ns = _sector_mask(az, 0) | _sector_mask(az, 180)
    # signed sector: +1 EW, -1 NS, 0 outside
    pa_signed = np.where(ew, pa, np.where(ns, -pa, np.nan))

    out_classes = {}
    for name, cmask in [("intra_NH", cls_nh), ("intra_SH", cls_sh), ("cross_hemi", cls_x)]:
        # pa_diff channel: monthly mean of signed PA over EW/NS sectors
        sector_ok = cmask & ~np.isnan(pa_signed)
        mk_d, pa_diff, n_d, dc_d = monthly_series(month_key, day_ord, np.nan_to_num(pa_signed), sector_ok)

        # amplitude channel: monthly mean PA over all sectors
        mk_a, mean_pa, n_a, dc_a = monthly_series(month_key, day_ord, pa, cmask)

        entry = {
            "n_pair_days": int(cmask.sum()),
            "pa_diff_monthly": {"months": mk_d, "values": [float(v) for v in pa_diff],
                                "n_pairs": [int(x) for x in n_d], "n_days": [int(x) for x in dc_d]},
            "mean_pa_monthly": {"months": mk_a, "values": [float(v) for v in mean_pa],
                                "n_pairs": [int(x) for x in n_a], "n_days": [int(x) for x in dc_a]},
        }
        if len(pa_diff) >= 8:
            entry["pa_diff_sinusoid"] = fit_annual_sinusoid(
                [month_mid_ordinal(m) for m in mk_d], pa_diff)
            entry["pa_diff_driver"] = regress_on_driver(mk_d, pa_diff)
        if len(mean_pa) >= 8:
            entry["mean_pa_sinusoid"] = fit_annual_sinusoid(
                [month_mid_ordinal(m) for m in mk_a], mean_pa)
            entry["mean_pa_driver"] = regress_on_driver(mk_a, mean_pa)
        out_classes[name] = entry
        print_status(f"  {name}: {len(mk_d)} months usable for pa_diff series", "INFO")

    # Hemisphere phase comparison on the pa_diff channel
    comp = {}
    nh = out_classes.get("intra_NH", {}).get("pa_diff_sinusoid")
    sh = out_classes.get("intra_SH", {}).get("pa_diff_sinusoid")
    if nh and sh:
        dphi = (nh["phase_deg"] - sh["phase_deg"] + 540.0) % 360.0 - 180.0
        se = math.hypot(nh["phase_se_deg"], sh["phase_se_deg"])
        comp["pa_diff_phase_NH_deg"] = nh["phase_deg"]
        comp["pa_diff_phase_SH_deg"] = sh["phase_deg"]
        comp["pa_diff_phase_difference_deg"] = float(dphi)
        comp["pa_diff_phase_difference_se_deg"] = float(se)
        comp["in_phase"] = bool(abs(dphi) < 90.0)
        comp["anti_phase"] = bool(abs(dphi) > 150.0)

    # Direct NH-vs-SH month-level correlation on shared months
    nhm = out_classes.get("intra_NH", {}).get("pa_diff_monthly", {})
    shm = out_classes.get("intra_SH", {}).get("pa_diff_monthly", {})
    shared = sorted(set(nhm.get("months", [])) & set(shm.get("months", [])))
    if len(shared) >= 8:
        nhi = [nhm["months"].index(m) for m in shared]
        shi = [shm["months"].index(m) for m in shared]
        a = np.array([nhm["values"][i] for i in nhi])
        b = np.array([shm["values"][i] for i in shi])
        r = float(np.corrcoef(a, b)[0, 1])
        n = len(a)
        cnt = 0
        for shift in range(1, n):
            if abs(float(np.corrcoef(np.roll(a, shift), b)[0, 1])) >= abs(r):
                cnt += 1
        comp["nh_sh_monthly_corr"] = {
            "r": r,
            "perm_p": (cnt + 1) / n,
            "n_shared_months": n,
            "perm_method": "all n-1 circular shifts; discrete floor 1/n",
        }

    # Driver-sign consistency: a global driver predicts the same sign of
    # the orbital-velocity correlation in both hemispheres.
    nh_d = out_classes.get("intra_NH", {}).get("pa_diff_driver", {})
    sh_d = out_classes.get("intra_SH", {}).get("pa_diff_driver", {})
    if nh_d and sh_d:
        comp["speed_corr_sign_consistent"] = bool(
            np.sign(nh_d["r_vs_orbital_speed"]) == np.sign(sh_d["r_vs_orbital_speed"]))
        comp["vproj_corr_sign_consistent"] = bool(
            np.sign(nh_d["r_vs_cmb_vproj"]) == np.sign(sh_d["r_vs_cmb_vproj"]))

    results["classes"] = out_classes
    results["comparison"] = comp

    if comp.get("anti_phase"):
        results["verdict"] = ("NH and SH pa_diff modulations are anti-phased -- consistent with "
                              "a hemisphere-local seasonal driver, not a global kinematic coupling.")
    elif comp.get("in_phase") and comp.get("vproj_corr_sign_consistent"):
        results["verdict"] = ("NH and SH pa_diff modulations share a common annual phase and a "
                              "consistent orbital-velocity sign -- consistent with a global "
                              "heliocentric driver, not hemisphere-local seasons.")
    else:
        results["verdict"] = ("Phase structure does not cleanly discriminate: inspect per-hemisphere "
                              "sinusoid significance and driver correlations.")

    out = OUTPUTS_DIR / "step_3_4_hemisphere_phase_test.json"
    with open(out, "w") as f:
        json.dump(results, f, indent=2)
    print_status(f"  Wrote {out}", "INFO")
    print_status(f"  Verdict: {results['verdict']}", "INFO")
    return {"status": "success", "results": str(out)}


if __name__ == "__main__":
    Step34HemispherePhaseTest().run()
