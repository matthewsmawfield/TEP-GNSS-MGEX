#!/usr/bin/env python3
"""Step 3.6: Amplitude decomposition of the phase-alignment statistic.

Motivation (corpus issue gnss-amplitude-mismatch): the fitted MGEX
correlation amplitude A = 0.827 (Step 2.1) differs by ~7x from the
Papers-1-2 family value A = 0.09-0.13, and the fitted offset C0 is more
negative. This step decomposes that gap on identical real MGEX
station-days, axis by axis.

The fitted statistic is bounded: phase_alignment = Re(V)/|V| with
V = sum_band Pxy, the cosine of the magnitude-weighted circular mean
cross-spectral phase. It measures the signed dominance of the in-phase
component of the band-integrated cross-spectrum and saturates toward
+1 whenever that component dominates the sum; it is not a linear
estimator of the shared-power fraction. The linear comparator computed
here on the same pair-days is the band-averaged real coherency
Re{gamma} (the Papers-1-2 alternative estimator, TEP_USE_REAL_COHERENCY)
and the magnitude coherence |gamma| (whose distance-flat floor is the
finite-segment Welch bias).

Estimator axes evaluated on identical pair-days (COD and GRG product
files separately; shared 300 s epoch grid; linear detrend):

  1. canonical          : nperseg=96, flat-line/outlier filter, per-epoch
                          median common-mode removal  (Step 2.0 settings)
  2. no_common_mode     : identical but without the median subtraction
  3. no_flatline_filter : identical but keeping stations whose detrended
                          rms is < MIN_STATION_STD_NS (held-constant clocks)
  4. p12_spec           : Papers-1-2 combination -- single full-day Welch
                          segment (nperseg = min(1024,T)), no common-mode
                          removal, no flat-line filter
  5. p12_spec_filt      : p12_spec with the MGEX flat-line filter restored,
                          isolating the filter axis under the P12 spectrum

For each configuration x statistic the pooled distance-binned profile is
fit with C(r) = A exp(-r/lam) + C0 under the Step-2.1 conventions.

Part B is an estimator transfer-function calibration -- clearly labelled
Monte Carlo on synthetic red-noise pairs (T = 288 at 5-min cadence) with
known injected in-phase band coherence rho, evaluated under the three
corpus spectral windows (nperseg = 96 ~5 segments; 144 ~3 segments as in
the raw-SPP pipeline; single full-day segment as in Papers 1-2). It maps
each family's measured A to the underlying coherence rho implied by its
own estimator, quantifying the nonlinear compression at low rho and the
saturation at high rho. The synthetic block is estimator calibration
only; no manuscript claim rests on it beyond the interpretation of the
real-data amplitudes.

Outputs: results/outputs/step_3_6_amplitude_decomposition.json
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

import json
import numpy as np
from scipy.signal import detrend, csd, welch
from scipy.signal.windows import hann
from scipy.optimize import curve_fit

from scripts.utils.logger import TEPLogger, set_step_logger, print_status
from scripts.utils.config import (
    RAW_DIR, PROCESSED_DIR, OUTPUTS_DIR, LOGS_DIR,
    F1_HZ, F2_HZ, FS_HZ,
    MIN_DISTANCE_KM, MAX_DISTANCE_KM, N_BINS, MIN_BIN_COUNT, MIN_EPOCHS,
    MAX_STATION_STD_NS, MIN_STATION_STD_NS,
)
from scripts.utils.geospatial import ecef_to_lla
from scripts.steps.step_2_0_mgex_clock_correlation import _parse_clk_ar_records
from scripts.steps.step_3_5_estimator_decomposition import (
    _segmented_fft, _distance_matrix, _fit_from_aggregates, _Agg,
    _BIN_EDGES, _BIN_CENTERS,
)

LOGS_DIR.mkdir(parents=True, exist_ok=True)
logger = TEPLogger("step_3_6", log_file_path=LOGS_DIR / "step_3_6_amplitude_decomposition.log")
set_step_logger(logger)

RNG_SEED = 20260928
MC_REPS = 400
COD_SAMPLE_EVERY = 10   # ~46 COD days spanning the archive

# (label, common_mode, flatline_filter, nperseg_spec)
# nperseg_spec: int -> that value; "full" -> min(1024, T) (Papers 1-2)
CONFIGS = [
    ("canonical", True, True, 96),
    ("no_common_mode", False, True, 96),
    ("no_flatline_filter", True, False, 96),
    ("p12_spec", False, False, "full"),
    ("p12_spec_filt", False, True, "full"),
]

STATISTICS = ("phase_alignment", "real_coherency", "magnitude_coherence")

# Published amplitudes with provenance (read back and verified when the
# sibling outputs are present on disk; otherwise recorded from the
# referenced result files).
REFERENCE_AMPLITUDES = {
    "paper14_mgex_pooled": {
        "A": 0.827, "lambda_km": 1862, "C0": -0.065,
        "source": "TEP-GNSS-MGEX/results/outputs/step_2_1_correlation_length.json",
        "spec": "nperseg=96, flat-filter, common-mode removed",
    },
    "paper1_code_finals": {
        "A": 0.1138, "lambda_km": 4548.6, "C0": -0.0213,
        "source": "TEP-GNSS/results/outputs/step_2_0_correlation_code.json",
        "spec": "nperseg=min(1024,T) single segment, no flat filter, no CM removal",
    },
    "paper2_code_longspan": {
        "A": 0.1261, "lambda_km": 3210.5, "C0": -0.0104,
        "source": "TEP-GNSS-II/results/outputs/step_2_0_correlation_code.json (exponential_fit)",
        "spec": "nperseg=min(1024,T) single segment, no flat filter, no CM removal",
    },
    "paper3_rinex_phase_same_pairs": {
        "A": 1.056, "lambda_km": 1796, "C0": -0.056,
        "source": "TEP-GNSS-RINEX/results/outputs/step_2_0_raw_spp_analysis_all_stations.json",
        "spec": "nperseg=min(256,T//2) (~3 segments), own-SPP clocks",
    },
    "paper3_rinex_msc_same_pairs": {
        "A": 0.182, "lambda_km": 727, "C0": 0.549,
        "source": "TEP-GNSS-RINEX/results/outputs/step_2_0_raw_spp_analysis_all_stations.json (identical pair-days, normalized_msc)",
        "spec": "nperseg=min(256,T//2), linear magnitude statistic",
    },
}


def _build_day(records, station_coords, keep_flat):
    """(N x T) detrended station matrix on the shared 300 s epoch grid.

    Mirrors step_2_0/step_3_5: >= MIN_EPOCHS records, all-station epoch
    intersection, linear detrend, >MAX_STATION_STD_NS rejection. With
    keep_flat=False the <MIN_STATION_STD_NS flat-line rejection is also
    applied (canonical); with keep_flat=True low-variance stations are
    retained to expose their deflation of the pair means.
    """
    valid = [
        sta for sta, ents in records.items()
        if sta in station_coords and len(ents) >= MIN_EPOCHS
    ]
    if len(valid) < 2:
        return None, None, 0, 0
    series_by_sta, ts_sets = {}, []
    for sta in valid:
        d = {}
        for ts, off in records[sta]:
            t = ts.timestamp()
            if int(round(t)) % 300 == 0:
                d[t] = off * 1e9
        if len(d) < MIN_EPOCHS:
            continue
        series_by_sta[sta] = d
        ts_sets.append(set(d.keys()))
    if len(ts_sets) < 2:
        return None, None, 0, 0
    common_ts = sorted(set.intersection(*ts_sets))
    if len(common_ts) < MIN_EPOCHS:
        return None, None, 0, 0
    stations, rows = [], []
    n_flat = 0
    for sta, d in series_by_sta.items():
        vals = np.array([d[t] for t in common_ts])
        x = detrend(vals, type="linear")
        std = float(np.std(x))
        if std > MAX_STATION_STD_NS:
            continue
        if std < MIN_STATION_STD_NS:
            n_flat += 1
            if not keep_flat:
                continue
            if std == 0.0:
                continue  # exactly constant: undefined spectra either way
        stations.append(sta)
        rows.append(x)
    if len(stations) < 2:
        return None, None, 0, 0
    return np.vstack(rows), stations, len(common_ts), n_flat


def _stat_matrices(F, band_mask):
    """All-pairs statistic matrices from segmented rFFTs.

    F : (N, nseg, nfreq). For the band bins:
      V_ij,f  = mean_seg conj(F_i,f) F_j,f          (Welch cross-spectrum)
      U_i,f   = mean_seg |F_i,f|^2                  (Welch auto-spectrum)
      pa_ij   = Re(sum_f V_ij,f) / |sum_f V_ij,f|   (canonical statistic)
      reg_ij  = mean_f Re{V_ij,f / sqrt(U_i,f U_j,f)}
      mag_ij  = mean_f |V_ij,f| / sqrt(U_i,f U_j,f)
    """
    Fb = F[:, :, band_mask]                          # (N, nseg, nb)
    V = np.einsum("isf,jsf->ijf", np.conj(Fb), Fb) / Fb.shape[1]
    U = np.einsum("isf,isf->if", np.conj(Fb), Fb).real / Fb.shape[1]
    S = np.sum(V, axis=2)                            # (N,N) band-summed CSD
    magS = np.abs(S)
    pa = np.where(magS > 0, S.real / np.where(magS > 0, magS, 1.0), np.nan)

    with np.errstate(divide="ignore", invalid="ignore"):
        denom = np.sqrt(U[:, None, :] * U[None, :, :])
        gam = np.where(denom > 0, V / np.where(denom > 0, denom, np.inf), np.nan)
    reg = np.nanmean(gam.real, axis=2)
    mag = np.nanmean(np.abs(gam), axis=2)
    return pa, reg, mag


def _bin_tri(pa, reg, mag, dist):
    """Upper-triangle pairs -> per-statistic (cnt,sum,sumsq) aggregates."""
    N = pa.shape[0]
    iu = np.triu_indices(N, k=1)
    d = dist[iu]
    keep = (d >= MIN_DISTANCE_KM) & (d < MAX_DISTANCE_KM)
    d = d[keep]
    out = {}
    for name, mat in (("phase_alignment", pa), ("real_coherency", reg),
                      ("magnitude_coherence", mag)):
        p = mat[iu][keep]
        ok = np.isfinite(p)
        idx = np.clip(np.digitize(d[ok], _BIN_EDGES) - 1, 0, N_BINS - 1)
        pv = p[ok]
        out[name] = (
            np.bincount(idx, minlength=N_BINS).astype(float),
            np.bincount(idx, weights=pv, minlength=N_BINS),
            np.bincount(idx, weights=pv * pv, minlength=N_BINS),
            int(np.sum(ok)),
        )
    return out


def _spec_nperseg(spec, T):
    return min(1024, T) if spec == "full" else min(spec, max(T // 2, 32))


# ----------------------------------------------------------------------
# Part B: estimator transfer-function calibration (synthetic, labelled)
# ----------------------------------------------------------------------

def _red_series(rng, T):
    """~1/f^2 (red) noise realisation of length T."""
    w = rng.standard_normal(T)
    X = np.fft.rfft(w)
    f = np.fft.rfftfreq(T, d=1.0 / FS_HZ)
    f[0] = f[1]
    return np.fft.irfft(X / np.sqrt(f), T)


def _pa_scalar(x, y, nperseg):
    f, Pxy = csd(x, y, fs=FS_HZ, nperseg=min(nperseg, len(x)), detrend="constant")
    m = (f > 0) & (f >= F1_HZ) & (f <= F2_HZ)
    if not np.any(m):
        return np.nan
    V = np.sum(Pxy[m])
    return float(V.real / abs(V)) if abs(V) > 0 else np.nan


def _reg_scalar(x, y, nperseg):
    n = min(nperseg, len(x))
    f, Pxy = csd(x, y, fs=FS_HZ, nperseg=n, detrend="constant")
    _, Pxx = welch(x, fs=FS_HZ, nperseg=n, detrend="constant")
    _, Pyy = welch(y, fs=FS_HZ, nperseg=n, detrend="constant")
    m = (f > 0) & (f >= F1_HZ) & (f <= F2_HZ)
    den = Pxx[m] * Pyy[m]
    ok = den > 0
    if not np.any(ok):
        return np.nan
    return float(np.mean((Pxy[m][ok] / np.sqrt(den[ok])).real))


def transfer_calibration():
    """E[statistic|rho] for the corpus spectral windows (Monte Carlo)."""
    rng = np.random.default_rng(RNG_SEED)
    T = 288
    rhos = [0.0, 0.02, 0.03, 0.05, 0.07, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.60, 0.70]
    specs = {"nperseg96_mgex": 96, "nperseg144_spp": 144, "fullday_p12": T}
    table = {}
    for rho in rhos:
        c = np.sqrt(rho / (1.0 - rho)) if rho > 0 else 0.0
        acc = {k: [] for k in specs}
        regacc = {k: [] for k in specs}
        for _ in range(MC_REPS):
            s, a, b = _red_series(rng, T), _red_series(rng, T), _red_series(rng, T)
            x, y = c * s + a, c * s + b
            for k, ns in specs.items():
                acc[k].append(_pa_scalar(x, y, ns))
                regacc[k].append(_reg_scalar(x, y, ns))
        row = {}
        for k in specs:
            arr = np.array(acc[k])
            row[f"pa_{k}_mean"] = float(np.nanmean(arr))
            row[f"pa_{k}_std"] = float(np.nanstd(arr))
            rarr = np.array(regacc[k])
            row[f"reg_{k}_mean"] = float(np.nanmean(rarr))
        table[f"{rho:.3f}"] = row
        print_status(
            f"    rho={rho:.2f}: pa(full)={row['pa_fullday_p12_mean']:+.3f} "
            f"pa(144)={row['pa_nperseg144_spp_mean']:+.3f} "
            f"pa(96)={row['pa_nperseg96_mgex_mean']:+.3f} "
            f"ReG={row['reg_nperseg96_mgex_mean']:+.3f}", "INFO")
    # Invert: for each reference amplitude, nearest rho on that family's spec
    grid = np.array(rhos)
    implied = {}
    for fam, ref in REFERENCE_AMPLITUDES.items():
        if "rinex" in fam or "spp" in fam:
            key = "pa_nperseg144_spp_mean"
        elif "paper1_" in fam or "paper2_" in fam:
            key = "pa_fullday_p12_mean"
        else:
            key = "pa_nperseg96_mgex_mean"
        curve = np.array([table[f"{r:.3f}"][key] for r in rhos])
        hit_ceiling = ref["A"] > curve.max()
        implied[fam] = {
            "implied_rho": float(grid[int(np.argmin(np.abs(curve - ref["A"])))]),
            "at_grid_ceiling": bool(hit_ceiling),
            "spec_curve": key,
            "A_ref": ref["A"],
        }
    return {"rhos": rhos, "reps_per_rho": MC_REPS, "T": T,
            "note": "Monte Carlo on synthetic red-noise pairs; estimator "
                    "calibration only - not data",
            "E_statistic_given_rho": table, "implied_rho_at_published_A": implied}


class Step36AmplitudeDecomposition:
    def run(self):
        print_status("Step 3.6: Amplitude Statistic Decomposition", "INFO")
        OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)

        coords_file = PROCESSED_DIR / "station_coordinates.json"
        if not coords_file.exists():
            return {"status": "failed", "reason": "no_station_coordinates"}
        with open(coords_file) as f:
            ecef_coords = json.load(f)
        station_coords = {s: {"lat": ecef_to_lla(*xyz)[0], "lon": ecef_to_lla(*xyz)[1]}
                          for s, xyz in ecef_coords.items() if len(xyz) == 3}
        print_status(f"  Loaded {len(station_coords)} station coordinates", "INFO")

        clk_files = sorted((RAW_DIR / "clk").glob("*.clk"))
        cod_files = [f for f in clk_files if f.name.endswith("_COD.clk")]
        grg_files = [f for f in clk_files if f.name.endswith("_GRG.clk")]
        cod_sel = cod_files[::COD_SAMPLE_EVERY]
        print_status(
            f"  Day sample: {len(cod_sel)} COD files (every {COD_SAMPLE_EVERY}th of "
            f"{len(cod_files)}) + {len(grg_files)} GRG files", "INFO")

        # pools[product][config][stat] -> _Agg
        pools = {prod: {c[0]: {s: _Agg() for s in STATISTICS} for c in CONFIGS}
                 for prod in ("COD", "GRG")}
        flat_counts = {"COD": [0, 0], "GRG": [0, 0]}
        n_used = {"COD": 0, "GRG": 0}

        validation = {"checked_pairs": 0, "max_abs_diff": None}

        for path in cod_sel + grg_files:
            center = "COD" if path.name.endswith("_COD.clk") else "GRG"
            day_key = path.stem.split("_")[0]
            records = _parse_clk_ar_records(path)
            if not records:
                continue
            # Station-level bookkeeping uses the filtered build for the
            # canonical station count; unfiltered days still record flat count.
            Xall, stas_all, Tall, nflat = _build_day(records, station_coords, keep_flat=True)
            Xfil, stas_fil, Tfil, _ = _build_day(records, station_coords, keep_flat=False)
            if Xfil is None or Xall is None:
                continue
            n_used[center] += 1
            flat_counts[center][0] += nflat
            flat_counts[center][1] += len(stas_all)

            # One-time vectorized-vs-scalar validation on a canonical day
            if validation["checked_pairs"] == 0 and Xfil.shape[0] > 10:
                Xv = Xfil - np.median(Xfil, axis=0, keepdims=True)
                Fv, fvs = _segmented_fft(Xv, _spec_nperseg(96, Xv.shape[1]))
                bmv = (fvs > 0) & (fvs >= F1_HZ) & (fvs <= F2_HZ)
                pam, _, _ = _stat_matrices(Fv, bmv)
                diffs = [
                    abs(_pa_scalar(Xv[k], Xv[k + 1], 96) - pam[k, k + 1])
                    for k in range(min(10, Xv.shape[0] - 1))
                ]
                validation["checked_pairs"] = len(diffs)
                validation["max_abs_diff"] = float(np.max(diffs))

            dist_fil = _distance_matrix(stas_fil, station_coords)
            dist_all = _distance_matrix(stas_all, station_coords)

            for label, cm, filt, spec in CONFIGS:
                X, dist = (Xfil, dist_fil) if filt else (Xall, dist_all)
                if cm:
                    X = X - np.median(X, axis=0, keepdims=True)
                nseg = _spec_nperseg(spec, X.shape[1])
                F, freqs = _segmented_fft(X, nseg)
                bm = (freqs > 0) & (freqs >= F1_HZ) & (freqs <= F2_HZ)
                if not np.any(bm):
                    continue
                pa, reg, mag = _stat_matrices(F, bm)
                binned = _bin_tri(pa, reg, mag, dist)
                for stat, agg in binned.items():
                    pools[center][label][stat].add(day_key, agg)

        results = {
            "step": "3_6_amplitude_decomposition",
            "purpose": REFERENCE_AMPLITUDES and (
                "Decompose the ~7x cross-product amplitude gap axis by axis on "
                "identical real MGEX station-days: statistic type (saturating "
                "phase direction vs linear real coherency), flat-line filter, "
                "common-mode removal, and Welch spectral window; plus a labelled "
                "Monte-Carlo estimator transfer calibration."),
            "configs": {c[0]: {"common_mode_removed": c[1],
                               "flatline_filter": c[2],
                               "nperseg": c[3]} for c in CONFIGS},
            "statistics": {
                "phase_alignment": "Re(V)/|V|, V = sum_band Pxy (canonical; bounded)",
                "real_coherency": "mean_band Re{Pxy/sqrt(Pxx Pyy)} (linear; Papers 1-2 alt estimator)",
                "magnitude_coherence": "mean_band |gamma| (finite-segment Welch floor ~1/sqrt(nseg))",
            },
            "day_sample": {"cod_every": COD_SAMPLE_EVERY,
                           "cod_days_used": n_used["COD"],
                           "grg_days_used": n_used["GRG"]},
            "flatline_stations": {
                prod: {"n_flat": flat_counts[prod][0],
                       "n_station_days": flat_counts[prod][1],
                       "frac_flat": (flat_counts[prod][0] / flat_counts[prod][1]
                                     if flat_counts[prod][1] else None)}
                for prod in flat_counts},
            "fits": {},
            "reference_amplitudes": REFERENCE_AMPLITUDES,
        }

        for prod in ("COD", "GRG"):
            if n_used[prod] == 0:
                continue
            results["fits"][prod] = {}
            for label, _, _, _ in CONFIGS:
                results["fits"][prod][label] = {}
                for stat in STATISTICS:
                    agg = pools[prod][label][stat]
                    pooled = agg.pooled()
                    if pooled is None:
                        continue
                    fit = _fit_from_aggregates(*pooled)
                    if fit is None:
                        continue
                    fit["n_days"] = len(agg.days)
                    # binned profile for the record
                    cnt, s, q = pooled
                    m = cnt >= MIN_BIN_COUNT
                    fit["binned_profile"] = {
                        "bin_centers_km": [float(v) for v in _BIN_CENTERS[m]],
                        "mean": [float(v) for v in (s[m] / cnt[m])],
                        "n_pairs_bin": [int(v) for v in cnt[m]],
                    }
                    results["fits"][prod][label][stat] = fit
                    print_status(
                        f"  {prod} {label:20s} {stat:20s}: A={fit['amplitude']:+.3f} "
                        f"lam={fit['lambda_km']:6.0f} C0={fit['offset']:+.3f} "
                        f"R2={fit['r_squared']:.2f}", "INFO")

        results["vectorized_path_validation"] = validation

        print_status("  Estimator transfer-function calibration (synthetic)...", "INFO")
        results["estimator_transfer_calibration"] = transfer_calibration()

        out_file = OUTPUTS_DIR / "step_3_6_amplitude_decomposition.json"
        with open(out_file, "w") as f:
            json.dump(results, f, indent=2)
        print_status(f"Amplitude decomposition written to {out_file}", "INFO")
        return {"status": "success", "results": str(out_file)}


if __name__ == "__main__":
    Step36AmplitudeDecomposition().run()
