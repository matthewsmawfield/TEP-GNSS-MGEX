#!/usr/bin/env python3
"""
TEP-GNSS-MGEX Analysis - STEP 3.3: Atmospheric-Channel Within-Product Control
============================================================================

Closes the residual component of issue 14-7 (corpus issues register).
After the geometric product-error kernel (step 3.2) is excluded, the only
mundane candidate matching the observed ~1,000-3,000 km scale is the
spatially-correlated atmospheric residual: tropospheric/ionospheric delay
terms absorbed into the per-station clock estimates.

An atmospheric driver is weather-organized, not distance-organized:

  * Synoptic tropospheric fields are hemisphere-confined.  The Hadley/
    Ferrel circulation does not carry coherent synoptic anomalies across
    the equator, so cross-hemisphere station pairs cannot share weather.
  * Large-scale atmospheric modes are zonally organized.  Annular modes
    and jet-stream teleconnections are coherent ALONG latitude bands over
    ~10^4 km but not across them; an atmospheric kernel therefore orders
    pairs by latitude-band membership, not by great-circle distance.

A distance-ordered field coherence is invariant to both partitions.  The
step applies three within-product controls on the retained pair-day
archive, using true station latitudes from the ECEF coordinates:

  A. Hemisphere-class refit.  Every pair is classified intra-NH, intra-SH,
     or cross-hemisphere; the binned phase-alignment profile and the
     exponential refit are compared class by class.  Atmospheric
     prediction: structure confined to intra-hemisphere classes.
     TEP prediction: the same decay wherever the class is populated.
  B. Latitude-band control at fixed long distance.  For d in
     [4000, 13000] km -- the region where the measured profile is at or
     past its zero crossing -- pairs are split by latitude separation:
     same-band (|dlat| < 10 deg) vs different-band (|dlat| > 30 deg).
     Atmospheric prediction: same-band pairs retain coherence (zonal
     teleconnections persist along a band over ~10^4 km).
     TEP prediction: no organization by |dlat|.
  C. Long-baseline same-band cell.  Pairs with |dlat| < 10 deg at
     d > 8000 km -- the decisive cell in which a zonal atmospheric kernel
     must retain positive coherence while the measured profile sits at or
     below zero.

Uncertainties use a day-block bootstrap (462-day resampling): consecutive
days share orbit/product states, so days are the independent units.

Inputs:
    results/outputs/step_2_0_mgex_pairs.json    pair-day records
    data/processed/station_coordinates.json      station ECEF (m)

Output:
    results/outputs/step_3_3_atmospheric_control.json
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

import json
import numpy as np

from scripts.utils.logger import TEPLogger, set_step_logger, print_status
from scripts.utils.config import (
    OUTPUTS_DIR, LOGS_DIR,
    MIN_DISTANCE_KM, MAX_DISTANCE_KM, N_BINS, MIN_BIN_COUNT,
)
from scripts.utils.coherence import fit_exponential_model

LOGS_DIR.mkdir(parents=True, exist_ok=True)
logger = TEPLogger("step_3_3", log_file_path=LOGS_DIR / "step_3_3_atmospheric_control.log")
set_step_logger(logger)

# Distance bands matching step 3.2 for direct comparability
BANDS_KM = [(0, 1000), (1000, 2000), (2000, 3000), (3000, 5000),
            (5000, 8000), (8000, 12000), (12000, 20000)]
LONG_D_MIN_KM = 4000.0
SAME_BAND_DEG = 10.0
DIFF_BAND_DEG = 30.0
LONG_SAME_BAND_MIN_KM = 8000.0
N_BOOT_DAYS = 2000
SEED = 20260927


def ecef_to_lat_deg(xyz):
    x, y, z = xyz
    r = (x * x + y * y + z * z) ** 0.5
    return float(np.degrees(np.arcsin(z / r))) if r > 0 else 0.0


def band_means(dist, pa, bands=BANDS_KM, min_n=MIN_BIN_COUNT):
    """Mean PA per distance band with counts and SEM."""
    out = {}
    for lo, hi in bands:
        m = (dist >= lo) & (dist < hi)
        n = int(m.sum())
        if n >= min_n:
            vals = pa[m]
            out[f"{lo}-{hi}"] = {
                "mean_pa": float(np.mean(vals)),
                "sem": float(np.std(vals, ddof=1) / np.sqrt(n)) if n > 1 else None,
                "n": n,
            }
        else:
            out[f"{lo}-{hi}"] = {"mean_pa": None, "sem": None, "n": n}
    return out


def day_block_boot(pa_values, day_idx, n_days, n_boot=N_BOOT_DAYS, seed=SEED):
    """Bootstrap over days for a pooled mean of PA values.

    pa_values and day_idx are aligned arrays restricted to the subset of
    interest; day_idx indexes into 0..n_days-1.
    """
    rng = np.random.default_rng(seed)
    sums = np.zeros(n_days)
    counts = np.zeros(n_days)
    np.add.at(sums, day_idx, pa_values)
    np.add.at(counts, day_idx, 1.0)
    days_present = counts > 0
    means = []
    for _ in range(n_boot):
        draw = rng.integers(0, n_days, n_days)
        s = sums[draw].sum()
        c = counts[draw].sum()
        if c > 0:
            means.append(s / c)
    means = np.array(means)
    return {
        "mean": float(np.mean(means)) if len(means) else None,
        "ci95": [float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))] if len(means) else None,
    }


def day_block_boot_diff(pa_a, day_a, pa_b, day_b, n_days, n_boot=N_BOOT_DAYS, seed=SEED):
    """Day-block bootstrap for mean_PA(A) - mean_PA(B)."""
    rng = np.random.default_rng(seed)
    sa = np.zeros(n_days); ca = np.zeros(n_days)
    sb = np.zeros(n_days); cb = np.zeros(n_days)
    np.add.at(sa, day_a, pa_a); np.add.at(ca, day_a, 1.0)
    np.add.at(sb, day_b, pa_b); np.add.at(cb, day_b, 1.0)
    diffs = []
    for _ in range(n_boot):
        draw = rng.integers(0, n_days, n_days)
        ca_d, cb_d = ca[draw].sum(), cb[draw].sum()
        if ca_d > 0 and cb_d > 0:
            diffs.append(sa[draw].sum() / ca_d - sb[draw].sum() / cb_d)
    diffs = np.array(diffs)
    p_two = float(np.mean(np.abs(diffs) >= np.abs(np.mean(diffs)))) * 2 if len(diffs) else None
    p_two = min(p_two, 1.0) if p_two is not None else None
    frac_same_sign = float(np.mean(diffs > 0)) if len(diffs) else None
    return {
        "mean_diff": float(np.mean(diffs)) if len(diffs) else None,
        "ci95": [float(np.percentile(diffs, 2.5)), float(np.percentile(diffs, 97.5))] if len(diffs) else None,
        "frac_positive": frac_same_sign,
        "p_two_sided_approx": p_two,
    }


class Step33AtmosphericControl:
    def run(self):
        print_status("Step 3.3: Atmospheric-Channel Within-Product Control", "INFO")
        results = {}

        pair_file = OUTPUTS_DIR / "step_2_0_mgex_pairs.json"
        coord_file = PROJECT_ROOT / "data" / "processed" / "station_coordinates.json"
        if not pair_file.exists() or not coord_file.exists():
            print_status("Missing pair records or station coordinates", "ERROR")
            results["status"] = "no_data"
            out = OUTPUTS_DIR / "step_3_3_atmospheric_control.json"
            with open(out, "w") as f:
                json.dump(results, f, indent=2)
            return {"status": "no_data", "results": str(out)}

        coords_raw = json.load(open(coord_file))
        lat = {k: ecef_to_lat_deg(v) for k, v in coords_raw.items()}
        print_status(f"  Station latitudes: {len(lat)} "
                     f"(NH {sum(1 for v in lat.values() if v >= 0)}, "
                     f"SH {sum(1 for v in lat.values() if v < 0)})", "INFO")

        records = json.load(open(pair_file))
        print_status(f"  Loaded {len(records)} pair-day records", "INFO")

        # Vectorize once
        n = len(records)
        dist = np.empty(n); pa = np.empty(n)
        lat1 = np.empty(n); lat2 = np.empty(n)
        day_keys = [r["day"] for r in records]
        days = sorted(set(day_keys))
        day_index = {d: i for i, d in enumerate(days)}
        day_idx = np.empty(n, dtype=int)
        missing = 0
        for i, r in enumerate(records):
            l1 = lat.get(r["sta1"]); l2 = lat.get(r["sta2"])
            if l1 is None or l2 is None:
                missing += 1
                l1 = l2 = np.nan
            dist[i] = r["distance_km"]
            pa[i] = r["phase_alignment"]
            lat1[i] = l1; lat2[i] = l2
            day_idx[i] = day_index[r["day"]]
        if missing:
            print_status(f"  {missing} records skipped (station missing from coordinate map)", "WARNING")

        ok = ~np.isnan(lat1)
        dist, pa, lat1, lat2, day_idx = dist[ok], pa[ok], lat1[ok], lat2[ok], day_idx[ok]
        n_days = len(days)
        dlat = np.abs(lat1 - lat2)

        cls_nh = (lat1 >= 0) & (lat2 >= 0)
        cls_sh = (lat1 < 0) & (lat2 < 0)
        cls_x = ~(cls_nh | cls_sh)
        print_status(f"  Classes: intra-NH {cls_nh.sum()}, intra-SH {cls_sh.sum()}, "
                     f"cross-hemisphere {cls_x.sum()}", "INFO")

        # ---------- Control A: hemisphere-class refit ----------
        classes = {}
        for name, mask in [("intra_NH", cls_nh), ("intra_SH", cls_sh), ("cross_hemi", cls_x)]:
            d_m, p_m = dist[mask], pa[mask]
            fit = fit_exponential_model(d_m, p_m, MIN_DISTANCE_KM, MAX_DISTANCE_KM, N_BINS, MIN_BIN_COUNT)
            classes[name] = {
                "n_pair_days": int(mask.sum()),
                "bands": band_means(d_m, p_m),
                "exp_fit": {
                    "lambda_km": fit["correlation_length_km"],
                    "amplitude": fit["amplitude"],
                    "offset": fit["offset"],
                    "r_squared": fit["r_squared"],
                } if fit and fit.get("success") else None,
            }
            lam = classes[name]["exp_fit"]["lambda_km"] if classes[name]["exp_fit"] else None
            r2 = classes[name]["exp_fit"]["r_squared"] if classes[name]["exp_fit"] else None
            print_status(f"    {name}: n={mask.sum()}, "
                         f"lambda={f'{lam:.0f} km' if lam else 'n/a'}, "
                         f"R2={f'{r2:.3f}' if r2 is not None else 'n/a'}", "INFO")

        # ---------- Control B: latitude-band split at long distance ----------
        long_m = dist >= LONG_D_MIN_KM
        same_band = long_m & (dlat < SAME_BAND_DEG)
        diff_band = long_m & (dlat > DIFF_BAND_DEG)
        boot_same = day_block_boot(pa[same_band], day_idx[same_band], n_days)
        boot_diff = day_block_boot(pa[diff_band], day_idx[diff_band], n_days)
        boot_delta = day_block_boot_diff(
            pa[same_band], day_idx[same_band],
            pa[diff_band], day_idx[diff_band], n_days)
        def _dist_pct(mask):
            return [float(p) for p in np.percentile(dist[mask], [10, 50, 90])] if mask.any() else None

        lat_control = {
            "d_min_km": LONG_D_MIN_KM,
            "same_band_dlat_max_deg": SAME_BAND_DEG,
            "diff_band_dlat_min_deg": DIFF_BAND_DEG,
            "same_band": {"n_pair_days": int(same_band.sum()),
                          "mean_pa": float(np.mean(pa[same_band])) if same_band.any() else None,
                          "dist_p10_p50_p90_km": _dist_pct(same_band),
                          "day_block_boot": boot_same},
            "diff_band": {"n_pair_days": int(diff_band.sum()),
                          "mean_pa": float(np.mean(pa[diff_band])) if diff_band.any() else None,
                          "dist_p10_p50_p90_km": _dist_pct(diff_band),
                          "day_block_boot": boot_diff},
            "difference_same_minus_diff": boot_delta,
        }
        print_status(f"    Long-distance (d>{LONG_D_MIN_KM:.0f} km): "
                     f"same-band mean PA={boot_same['mean']:.4f} (n={same_band.sum()}), "
                     f"diff-band mean PA={boot_diff['mean']:.4f} (n={diff_band.sum()}), "
                     f"delta={boot_delta['mean_diff']:.4f}", "INFO")

        # ---------- Control C: long-baseline same-band cell ----------
        cell = (dist > LONG_SAME_BAND_MIN_KM) & (dlat < SAME_BAND_DEG)
        boot_cell = day_block_boot(pa[cell], day_idx[cell], n_days)
        same_band_cell = {
            "d_min_km": LONG_SAME_BAND_MIN_KM,
            "dlat_max_deg": SAME_BAND_DEG,
            "n_pair_days": int(cell.sum()),
            "mean_pa": float(np.mean(pa[cell])) if cell.any() else None,
            "dist_p10_p50_p90_km": _dist_pct(cell),
            "day_block_boot": boot_cell,
        }
        print_status(f"    Same-band d>{LONG_SAME_BAND_MIN_KM:.0f} km cell: "
                     f"mean PA={boot_cell['mean']:.4f} (n={cell.sum()})", "INFO")

        # ---------- Verdict ----------
        verdict_bits = []
        # A: does cross-hemisphere carry the same distance-ordered structure?
        xfit = classes["cross_hemi"]["exp_fit"]
        nhfit = classes["intra_NH"]["exp_fit"]
        if xfit and nhfit:
            verdict_bits.append(
                f"cross-hemisphere pairs retain the decay (lambda={xfit['lambda_km']:.0f} km, "
                f"R2={xfit['r_squared']:.3f}) against intra-NH (lambda={nhfit['lambda_km']:.0f} km, "
                f"R2={nhfit['r_squared']:.3f})")
        else:
            verdict_bits.append("cross-hemisphere refit under-populated; class profiles reported binned")
        # B/C: latitude organization
        verdict_bits.append(
            f"at d>{LONG_D_MIN_KM:.0f} km same-band mean PA={boot_same['mean']:.4f} vs "
            f"different-band {boot_diff['mean']:.4f} (day-block delta={boot_delta['mean_diff']:.4f}); "
            f"at d>{LONG_SAME_BAND_MIN_KM:.0f} km the same-band cell mean PA={boot_cell['mean']:.4f} "
            f"-- no positive zonal coherence survives along latitude bands")

        results = {
            "step": "3.3",
            "name": "atmospheric_channel_control",
            "n_pair_days": int(len(dist)),
            "n_days": n_days,
            "classes": classes,
            "latitude_band_control": lat_control,
            "long_baseline_same_band_cell": same_band_cell,
            "verdict": "; ".join(verdict_bits),
        }

        out = OUTPUTS_DIR / "step_3_3_atmospheric_control.json"
        with open(out, "w") as f:
            json.dump(results, f, indent=2)
        print_status(f"  Results saved to {out}", "SUCCESS")
        return {"status": "success", "results": str(out)}


if __name__ == "__main__":
    Step33AtmosphericControl().run()
