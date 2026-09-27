#!/usr/bin/env python3
"""Step 3.5: Estimator decomposition of the MGEX correlation scale.

Motivation (issue gnss-lambda-instability): the pooled MGEX correlation
length (lambda = 1,862 +/- 155 km, Paper 14) lies below the IGS-final
family (3,330-4,549 km, Papers 1-2), which those papers identify with
the geometric saturation radius R_T(M_earth) ~ 4,150 km. Both pipelines
measure the same statistic - the cosine of the magnitude-weighted
circular mean cross-spectral phase - but differ in documented estimator
settings:

  1. Spectral resolution: Papers 1-2 use nperseg = min(1024, n_points),
     i.e. a single full-day Welch segment at the 5-min common cadence
     (df ~ 11.6 uHz, first resolved bin 11.6 uHz). MGEX uses
     nperseg = 96 (df ~ 34.7 uHz, first resolved bin 34.7 uHz), so the
     two lowest-frequency (longest-range) bins are absent.
  2. Common mode: MGEX subtracts the per-epoch network median before
     pairing; Papers 1-2 do not.
  3. Product family: IGS final GPS-only receiver clocks vs CODE/GRG
     MGEX combined multi-GNSS receiver clocks.

On identical real MGEX station-days this step evaluates a factorial
decomposition of (1) and (2): all-pairs phase-alignment matrices are
computed under {nperseg 96, full-day} x {common-mode on, off}, with an
additional band-floor variant (>=34.7 uHz bins only on the fine grid)
separating spectral resolution from band content. A per-sub-band lambda
profile on the full-day matrix reproduces the Paper-1 multiband
ordering test on MGEX data. Because the station-days are fixed, any
shift in the fitted lambda is attributable to the estimator axis alone.

A cadence diagnostic on GRG days checks a mixed-cadence subtlety: GRG
products are natively 30-s sampled, so the configured FS_HZ = 1/300
mislabels the frequency grid on GRG station-days (the claimed 10-500 uHz
band is truly ~0.1-5 mHz there). The diagnostic evaluates GRG days
three ways - native cadence under the configured fs, native cadence
under the true fs, and 5-min decimation under the configured fs -
bounding the contribution of those days to the pooled scale.

Note on data structure: unlike step_2_0 (which pools all products into
calendar-day station lists), this step builds one station matrix per
product-day file, keeping the estimator comparison free of
cross-product timestamp mixing.

Outputs: results/outputs/step_3_5_estimator_decomposition.json
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

import json
import numpy as np
from scipy.signal import detrend, csd
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
from scripts.utils.coherence import exponential_decay, compute_coherence_phase
from scripts.steps.step_2_0_mgex_clock_correlation import _parse_clk_ar_records

LOGS_DIR.mkdir(parents=True, exist_ok=True)
logger = TEPLogger("step_3_5", log_file_path=LOGS_DIR / "step_3_5_estimator_decomposition.log")
set_step_logger(logger)

# MGEX effective band floor: first resolved Welch bin at nperseg=96,
# df = FS_HZ/96 ~ 34.7 uHz.
F_MGEX_FLOOR_HZ = FS_HZ / 96.0

# Paper-1 named sub-bands (microHz) for the per-band lambda profile.
PAPER1_BANDS_UHZ = {
    "10-20": (10.0, 20.0),
    "20-30": (20.0, 30.0),
    "30-40": (30.0, 40.0),
    "40-50": (40.0, 50.0),
    "50-75": (50.0, 75.0),
    "75-100": (75.0, 100.0),
    "100-200": (100.0, 200.0),
    "200-350": (200.0, 350.0),
    "350-500": (350.0, 500.0),
}

N_DAY_BOOTSTRAP = 200
EARTH_R_KM = 6371.0
RNG_SEED = 20260927

_BIN_EDGES = np.logspace(np.log10(MIN_DISTANCE_KM), np.log10(MAX_DISTANCE_KM), N_BINS + 1)
_BIN_CENTERS = (_BIN_EDGES[:-1] + _BIN_EDGES[1:]) / 2.0


def _segmented_fft(X, nperseg):
    """Hann-windowed, mean-detrended segmented rFFT (scipy csd conventions)."""
    T = X.shape[1]
    nperseg = min(nperseg, T)
    step = max(1, nperseg - nperseg // 2)
    win = hann(nperseg, sym=False)
    segs = []
    for s in range(0, T - nperseg + 1, step):
        seg = X[:, s:s + nperseg].astype(float)
        seg = seg - seg.mean(axis=1, keepdims=True)
        segs.append(seg * win)
    F = np.stack([np.fft.rfft(sg, axis=1) for sg in segs], axis=1)  # (N, nseg, nfreq)
    freqs = np.fft.rfftfreq(nperseg, d=1.0 / FS_HZ)
    return F, freqs


def _pa_matrix_from_F(F, band_mask):
    """All-pairs phase-alignment matrix from segmented FFTs.

    phase_alignment = cos(angle(sum_f Pxy_f)) = Re(V)/|V| with
    V_ij = sum_{f in band} mean_seg S_i(seg,f) * conj(S_j(seg,f)).
    This is identical to the magnitude-weighted circular mean phase used
    by compute_coherence_phase and by Papers 1-2.
    """
    Fb = F[:, :, band_mask]
    V = np.einsum("isf,jsf->ij", Fb, np.conj(Fb)) / Fb.shape[1]
    mag = np.abs(V)
    pa = np.where(mag > 0, V.real / np.where(mag > 0, mag, 1.0), np.nan)
    return pa


def _pair_dist(s1, s2, station_coords):
    """Haversine distance (km) between two stations; None if unmapped."""
    if s1 not in station_coords or s2 not in station_coords:
        return None
    lat1, lon1 = np.deg2rad(station_coords[s1]["lat"]), np.deg2rad(station_coords[s1]["lon"])
    lat2, lon2 = np.deg2rad(station_coords[s2]["lat"]), np.deg2rad(station_coords[s2]["lon"])
    a = np.sin((lat2 - lat1) / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2
    return float(EARTH_R_KM * 2 * np.arcsin(np.sqrt(np.clip(a, 0.0, 1.0))))


def _bin_pair_records(recs):
    """List of (dist_km, pa) -> (count, sum, sumsq) per log-distance bin."""
    cnt = np.zeros(N_BINS)
    s = np.zeros(N_BINS)
    q = np.zeros(N_BINS)
    for d, p in recs:
        if not np.isfinite(p) or d < MIN_DISTANCE_KM or d >= MAX_DISTANCE_KM:
            continue
        idx = int(np.clip(np.digitize(d, _BIN_EDGES) - 1, 0, N_BINS - 1))
        cnt[idx] += 1
        s[idx] += p
        q[idx] += p * p
    return cnt, s, q


def _distance_matrix(stations, station_coords):
    lats = np.deg2rad([station_coords[s]["lat"] for s in stations])
    lons = np.deg2rad([station_coords[s]["lon"] for s in stations])
    dlat = lats[:, None] - lats[None, :]
    dlon = lons[:, None] - lons[None, :]
    a = np.sin(dlat / 2) ** 2 + np.cos(lats[:, None]) * np.cos(lats[None, :]) * np.sin(dlon / 2) ** 2
    return EARTH_R_KM * 2 * np.arcsin(np.sqrt(np.clip(a, 0.0, 1.0)))


def _bin_pairs(pa, dist):
    """Upper-triangle pairs -> (count, sum, sumsq) per log-distance bin."""
    N = pa.shape[0]
    iu = np.triu_indices(N, k=1)
    d = dist[iu]
    p = pa[iu]
    keep = np.isfinite(p) & (d >= MIN_DISTANCE_KM) & (d < MAX_DISTANCE_KM)
    d = d[keep]
    p = p[keep]
    if len(d) == 0:
        return None
    idx = np.clip(np.digitize(d, _BIN_EDGES) - 1, 0, N_BINS - 1)
    cnt = np.bincount(idx, minlength=N_BINS).astype(float)
    s = np.bincount(idx, weights=p, minlength=N_BINS)
    q = np.bincount(idx, weights=p * p, minlength=N_BINS)
    return cnt, s, q


def _fit_from_aggregates(counts, sums, sumsq):
    """Fit C(r) = A exp(-r/lam) + C0 on bin aggregates (step_2_1 conventions)."""
    mask = counts >= MIN_BIN_COUNT
    if np.sum(mask) < 5:
        return None
    n = counts[mask]
    means = sums[mask] / n
    var = (sumsq[mask] - n * means**2) / np.maximum(n - 1, 1)
    stds = np.sqrt(np.maximum(var, 0.0))
    sem = stds / np.sqrt(n)
    good = sem > 0
    fallback = np.mean(sem[good]) if np.any(good) else 1.0
    sem = np.where(sem > 0, sem, fallback)
    x = _BIN_CENTERS[mask]
    try:
        lam0 = (MIN_DISTANCE_KM + MAX_DISTANCE_KM) / 2.0
        lam_hi = MAX_DISTANCE_KM * 2.0
        popt, pcov = curve_fit(
            exponential_decay, x, means,
            p0=[0.5, lam0, 0.0],
            sigma=sem,
            bounds=([0, MIN_DISTANCE_KM, -1], [2, lam_hi, 1]),
            maxfev=5000,
        )
        pred = exponential_decay(x, *popt)
        ss_res = float(np.sum((means - pred) ** 2))
        ss_tot = float(np.sum((means - np.mean(means)) ** 2))
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
        perr = np.sqrt(np.diag(pcov))
        return {
            "lambda_km": float(popt[1]),
            "lambda_err_km": float(perr[1]),
            "amplitude": float(popt[0]),
            "offset": float(popt[2]),
            "r_squared": float(r2),
            "n_bins_used": int(np.sum(mask)),
            "n_pairs": int(np.sum(counts)),
        }
    except Exception:
        return None


class _Agg:
    """Per-day binned aggregate store: day -> (cnt, sum, sumsq)."""

    def __init__(self):
        self.days = {}

    def add(self, day_key, agg):
        if agg is not None:
            self.days[day_key] = agg

    def pooled(self):
        if not self.days:
            return None
        cnt = np.sum([v[0] for v in self.days.values()], axis=0)
        s = np.sum([v[1] for v in self.days.values()], axis=0)
        q = np.sum([v[2] for v in self.days.values()], axis=0)
        return cnt, s, q

    def bootstrap(self, rng, n_boot=N_DAY_BOOTSTRAP):
        days = list(self.days.keys())
        if len(days) < 3:
            return {}
        lam = []
        for _ in range(n_boot):
            pick = rng.choice(days, size=len(days), replace=True)
            cnt = np.sum([self.days[k][0] for k in pick], axis=0)
            s = np.sum([self.days[k][1] for k in pick], axis=0)
            q = np.sum([self.days[k][2] for k in pick], axis=0)
            fit = _fit_from_aggregates(cnt, s, q)
            if fit:
                lam.append(fit["lambda_km"])
        if len(lam) < 20:
            return {}
        lam = np.array(lam)
        return {
            "sigma_lambda_km": float(np.std(lam, ddof=1)),
            "ci95_km": [float(np.percentile(lam, 2.5)), float(np.percentile(lam, 97.5))],
            "median_km": float(np.median(lam)),
            "n_boot_used": int(len(lam)),
        }


def _build_day_matrix(records, station_coords, decimate_5min):
    """Assemble aligned (N x T) station matrix for one product-day file.

    Mirrors step_2_0 filtering: stations with >= MIN_EPOCHS records and
    known coordinates; common-epoch intersection; linear detrend;
    outlier (>MAX_STATION_STD_NS) and flat-line (<MIN_STATION_STD_NS)
    rejection. With decimate_5min the per-station timestamp sets are
    restricted to the 5-min UTC boundary before intersection (30-s
    products such as GRG are thinned to the canonical cadence).
    """
    valid = [
        sta for sta, ents in records.items()
        if sta in station_coords and len(ents) >= MIN_EPOCHS
    ]
    if len(valid) < 2:
        return None, None, 0

    series_by_sta = {}
    ts_sets = []
    for sta in valid:
        d = {}
        for ts, off in records[sta]:
            d[ts.timestamp()] = off * 1e9  # seconds -> ns
        ts_arr = np.array(sorted(d.keys()))
        if decimate_5min:
            keep = (np.round(ts_arr).astype(np.int64) % 300) == 0
            ts_arr = ts_arr[keep]
        series_by_sta[sta] = d
        ts_sets.append(set(ts_arr.tolist()))

    common_ts = sorted(set.intersection(*ts_sets))
    if len(common_ts) < MIN_EPOCHS:
        return None, None, 0

    stations, rows = [], []
    for sta in valid:
        d = series_by_sta[sta]
        vals = np.array([d[t] for t in common_ts])
        x = detrend(vals, type="linear")
        std = float(np.std(x))
        if std > MAX_STATION_STD_NS or std < MIN_STATION_STD_NS:
            continue
        stations.append(sta)
        rows.append(x)

    if len(stations) < 2:
        return None, None, 0
    return np.vstack(rows), stations, len(common_ts)


def _pair_pa_scipy(x, y, fs, nperseg, f_lo, f_hi):
    """Per-pair phase alignment via scipy.signal.csd (reference path)."""
    f, Pxy = csd(x, y, fs=fs, nperseg=min(nperseg, len(x)), detrend="constant")
    mask = (f > 0) & (f >= f_lo) & (f <= f_hi)
    if not np.any(mask):
        return np.nan
    V = np.sum(Pxy[mask])
    return float(V.real / abs(V)) if abs(V) > 0 else np.nan


class Step35EstimatorDecomposition:
    def run(self):
        coords_file = PROCESSED_DIR / "station_coordinates.json"
        if not coords_file.exists():
            return {"status": "failed", "reason": "no_station_coordinates"}
        with open(coords_file) as f:
            ecef_coords = json.load(f)
        station_coords = {}
        for sta, xyz in ecef_coords.items():
            lat, lon = ecef_to_lla(*xyz)
            station_coords[sta] = {"lat": lat, "lon": lon}
        print_status(f"  Loaded {len(station_coords)} station coordinates", "INFO")

        clk_files = sorted((RAW_DIR / "clk").glob("*.clk"))
        if not clk_files:
            return {"status": "failed", "reason": "no_clk_files"}

        cod_files = [f for f in clk_files if f.name.endswith("_COD.clk")]
        grg_files = [f for f in clk_files if f.name.endswith("_GRG.clk")]
        # COD every ~8th file spans the full interval; all GRG files are
        # kept for the cadence diagnostic.
        cod_sel = cod_files[::8]
        print_status(
            f"  Day sample: {len(cod_sel)} COD files (every 8th of {len(cod_files)}) "
            f"+ {len(grg_files)} GRG files", "INFO")

        rng = np.random.default_rng(RNG_SEED)

        # (label, common_mode, nperseg_mode, f_lo)
        #   nperseg_mode: 96 -> MGEX canonical; 'full' -> min(1024, T) (Paper 1)
        configs = [
            ("mgex_canonical", True, 96, F1_HZ),
            ("mgex_no_cm", False, 96, F1_HZ),
            ("fullday_cm_fullband", True, "full", F1_HZ),
            ("fullday_nocm_fullband", False, "full", F1_HZ),
            ("fullday_cm_35uhz", True, "full", F_MGEX_FLOOR_HZ),
            ("fullday_nocm_35uhz", False, "full", F_MGEX_FLOOR_HZ),
        ]
        pools = {c[0]: _Agg() for c in configs}
        band_pools = {"fullday_cm": defaultdict_agg(), "fullday_nocm": defaultdict_agg()}
        grg_diag = {k: _Agg() for k in
                    ("native_cfgfs_cm", "native_truefs_cm", "decim_cfgfs_cm",
                     "pair_aligned_5min")}

        validation = {"checked_pairs": 0, "max_abs_diff": None,
                      "pearson_r": None}
        val_done = False

        n_used = {"COD": 0, "GRG": 0}
        for path in cod_sel + grg_files:
            center = "COD" if path.name.endswith("_COD.clk") else "GRG"
            day_key = path.stem.split("_")[0]
            records = _parse_clk_ar_records(path)
            if not records:
                continue

            # Canonical-cadence path: 5-min-decimated common grid.
            X, stations, T = _build_day_matrix(records, station_coords, decimate_5min=True)
            if X is not None:
                n_used[center] += 1
                X_cm = X - np.median(X, axis=0, keepdims=True)
                dist = _distance_matrix(stations, station_coords)

                F96, f96 = _segmented_fft(X, min(96, T // 2))
                F96_cm, _ = _segmented_fft(X_cm, min(96, T // 2))
                Ffull, ffull = _segmented_fft(X, min(1024, T))
                Ffull_cm, _ = _segmented_fft(X_cm, min(1024, T))

                masks96 = {
                    F1_HZ: (f96 > 0) & (f96 >= F1_HZ) & (f96 <= F2_HZ),
                }
                masksfull = {
                    F1_HZ: (ffull > 0) & (ffull >= F1_HZ) & (ffull <= F2_HZ),
                    F_MGEX_FLOOR_HZ: (ffull > 0) & (ffull >= F_MGEX_FLOOR_HZ) & (ffull <= F2_HZ),
                }
                F_by_cfg = {
                    ("mgex_canonical"): (F96_cm, masks96[F1_HZ]),
                    ("mgex_no_cm"): (F96, masks96[F1_HZ]),
                    ("fullday_cm_fullband"): (Ffull_cm, masksfull[F1_HZ]),
                    ("fullday_nocm_fullband"): (Ffull, masksfull[F1_HZ]),
                    ("fullday_cm_35uhz"): (Ffull_cm, masksfull[F_MGEX_FLOOR_HZ]),
                    ("fullday_nocm_35uhz"): (Ffull, masksfull[F_MGEX_FLOOR_HZ]),
                }
                pa_canonical = None
                for label, _, _, _ in configs:
                    Fmat, bm = F_by_cfg[label]
                    if not np.any(bm):
                        continue
                    pa = _pa_matrix_from_F(Fmat, bm)
                    if label == "mgex_canonical":
                        pa_canonical = pa
                    pools[label].add(day_key, _bin_pairs(pa, dist))

                # Per-sub-band lambda profile on the full-day matrix.
                for cm_label, Fmat in (("fullday_cm", Ffull_cm), ("fullday_nocm", Ffull)):
                    for bname, (lo_uhz, hi_uhz) in PAPER1_BANDS_UHZ.items():
                        bm = (ffull > 0) & (ffull >= lo_uhz * 1e-6) & (ffull <= hi_uhz * 1e-6)
                        if not np.any(bm):
                            continue
                        pa_b = _pa_matrix_from_F(Fmat, bm)
                        band_pools[cm_label][bname].add(day_key, _bin_pairs(pa_b, dist))

                # One-time numerical validation of the vectorized path.
                if not val_done and center == "COD" and pa_canonical is not None:
                    N = X_cm.shape[0]
                    refv, newv = [], []
                    nperseg_ref = min(96, X_cm.shape[1] // 2)
                    for _ in range(600):
                        a, b = rng.integers(0, N, 2)
                        if a == b:
                            continue
                        v_ref = _pair_pa_scipy(X_cm[a], X_cm[b], FS_HZ, nperseg_ref, F1_HZ, F2_HZ)
                        v_new = pa_canonical[a, b]
                        if np.isfinite(v_ref) and np.isfinite(v_new):
                            refv.append(v_ref)
                            newv.append(v_new)
                    if refv:
                        diffs = np.abs(np.array(refv) - np.array(newv))
                        validation["checked_pairs"] = int(len(diffs))
                        validation["max_abs_diff"] = float(np.max(diffs))
                        validation["pearson_r"] = float(np.corrcoef(refv, newv)[0, 1])
                        val_done = True

            # GRG cadence diagnostic on the native (30-s) grid.
            if center == "GRG":
                Xn, stations_n, Tn = _build_day_matrix(records, station_coords, decimate_5min=False)
                if Xn is not None and len(stations_n) >= 2:
                    Xn_cm = Xn - np.median(Xn, axis=0, keepdims=True)
                    dist_n = _distance_matrix(stations_n, station_coords)
                    nperseg96 = min(96, Tn // 2)
                    if nperseg96 >= 32:
                        Fn, fn = _segmented_fft(Xn_cm, nperseg96)
                        # Buggy path: configured fs=1/300 on 30-s data.
                        bm_bug = (fn > 0) & (fn >= F1_HZ) & (fn <= F2_HZ)
                        # fn built with FS_HZ already; for true fs the grid is 10x.
                        pa_bug = _pa_matrix_from_F(Fn, bm_bug)
                        grg_diag["native_cfgfs_cm"].add(day_key, _bin_pairs(pa_bug, dist_n))
                        # True fs: bins sit 10x higher; only ~347 uHz bin in band.
                        bm_true = (fn * 10.0 > 0) & (fn * 10.0 >= F1_HZ) & (fn * 10.0 <= F2_HZ)
                        if np.any(bm_true):
                            pa_true = _pa_matrix_from_F(Fn, bm_true)
                            grg_diag["native_truefs_cm"].add(day_key, _bin_pairs(pa_true, dist_n))
                # Corrected path: the decimated canonical config already
                # accumulated this day in pools["mgex_canonical"].
                if day_key in pools["mgex_canonical"].days:
                    grg_diag["decim_cfgfs_cm"].days[day_key] = pools["mgex_canonical"].days[day_key]

                # Corrected per-center estimate on pairwise 5-min-intersected
                # series: the step_3_1 GRG_indep estimator repaired — per-pair
                # timestamp intersection on the 300 s UTC grid (≥ MIN_EPOCHS
                # shared epochs) with the true cadence labels. Sparse
                # asynchronous coverage keeps many days below the all-station
                # common-grid minimum, but pairs intersected individually
                # retain far more usable content.
                sta_series = {}
                for sta, ents in records.items():
                    if sta not in station_coords or len(ents) < MIN_EPOCHS:
                        continue
                    d = {}
                    for ts, off in ents:
                        t = ts.timestamp()
                        if round(t) % 300 == 0:
                            d[t] = off * 1e9
                    if len(d) >= MIN_EPOCHS:
                        sta_series[sta] = d
                if len(sta_series) >= 2:
                    keys = sorted(sta_series)
                    day_recs = []
                    for i in range(len(keys)):
                        a = sta_series[keys[i]]
                        for j in range(i + 1, len(keys)):
                            b = sta_series[keys[j]]
                            common = sorted(set(a) & set(b))
                            if len(common) < MIN_EPOCHS:
                                continue
                            x = np.array([a[t] for t in common])
                            y = np.array([b[t] for t in common])
                            _, _, pa_ij = compute_coherence_phase(
                                x, y, FS_HZ, F1_HZ, F2_HZ)
                            if not np.isfinite(pa_ij):
                                continue
                            d_km = _pair_dist(keys[i], keys[j], station_coords)
                            if d_km is None:
                                continue
                            day_recs.append((d_km, pa_ij))
                    if day_recs:
                        grg_diag["pair_aligned_5min"].add(
                            day_key, _bin_pair_records(day_recs))

            total = n_used["COD"] + n_used["GRG"]
            if total % 10 == 0:
                print_status(f"  Processed {n_used['COD']} COD + {n_used['GRG']} GRG days", "PROCESS")

        print_status(f"  Usable days: {n_used['COD']} COD, {n_used['GRG']} GRG", "INFO")

        results = {}
        for label, cm, nps, flo in configs:
            agg = pools[label].pooled()
            if agg is None:
                continue
            fit = _fit_from_aggregates(*agg)
            if not fit:
                continue
            boot = pools[label].bootstrap(rng)
            fit["bootstrap"] = boot
            fit["n_days"] = len(pools[label].days)
            results[label] = fit
            print_status(
                f"  {label:24s} lambda = {fit['lambda_km']:8.1f} km  "
                f"R2 = {fit['r_squared']:.3f}  ({fit['n_days']} days, {fit['n_pairs']:,} pairs)",
                "INFO")

        band_results = {}
        for cm_label, bands in band_pools.items():
            band_results[cm_label] = {}
            for bname, aggobj in bands.items():
                agg = aggobj.pooled()
                if agg is None:
                    continue
                fit = _fit_from_aggregates(*agg)
                if fit:
                    fit["n_days"] = len(aggobj.days)
                    band_results[cm_label][bname] = fit

        grg_results = {}
        for label, aggobj in grg_diag.items():
            agg = aggobj.pooled()
            if agg is None:
                continue
            fit = _fit_from_aggregates(*agg)
            grg_results[label] = {
                "lambda_km": fit["lambda_km"] if fit else None,
                "r_squared": fit["r_squared"] if fit else None,
                "n_pairs": int(np.sum(agg[0])),
                "n_days": len(aggobj.days),
            }
            if fit:
                print_status(
                    f"  GRG {label:20s} lambda = {fit['lambda_km']:8.1f} km  "
                    f"R2 = {fit['r_squared']:.3f}", "INFO")

        # Pooled-lambda impact of the cadence handling on this day sample:
        # as-published = COD days canonical + GRG-only days on the buggy
        # native-cadence grid; fixed = same days with GRG decimated to 5 min.
        cadence_impact = {}
        cod_days = {k: v for k, v in pools["mgex_canonical"].days.items()
                    if k in {p.stem.split("_")[0] for p in cod_sel}}
        for tag, grg_source in (("as_published", grg_diag["native_cfgfs_cm"]),
                                ("cadence_fixed", grg_diag["decim_cfgfs_cm"])):
            merged = dict(cod_days)
            merged.update(grg_source.days)
            if not merged:
                continue
            cnt = np.sum([v[0] for v in merged.values()], axis=0)
            s = np.sum([v[1] for v in merged.values()], axis=0)
            q = np.sum([v[2] for v in merged.values()], axis=0)
            fit = _fit_from_aggregates(cnt, s, q)
            if fit:
                fit["n_days"] = len(merged)
                cadence_impact[tag] = fit
                print_status(
                    f"  pooled {tag:16s} lambda = {fit['lambda_km']:8.1f} km "
                    f"({len(merged)} days)", "INFO")

        out = {
            "step": "3_5_estimator_decomposition",
            "purpose": "Isolate estimator-conditioning contributions (Welch "
                       "segment length / band floor, per-epoch common-mode "
                       "removal) to the fitted correlation length on identical "
                       "real MGEX station-days; includes a GRG native-cadence "
                       "frequency-axis diagnostic.",
            "configs": {
                label: {"common_mode_removed": cm,
                        "nperseg": ("min(1024,T)" if nps == "full" else nps),
                        "f_lo_hz": flo, "f_hi_hz": F2_HZ}
                for label, cm, nps, flo in configs
            },
            "day_sample": {
                "cod_files_every": 8,
                "cod_days_used": n_used["COD"],
                "grg_days_used": n_used["GRG"],
            },
            "pooled_fits": results,
            "per_band_lambda": band_results,
            "grg_cadence_diagnostic": grg_results,
            "grg_cadence_pooled_impact": cadence_impact,
            "vectorized_path_validation": validation,
            "reference_values": {
                "mgex_canonical_pooled_lambda_km": 1861.6,
                "igs_final_family_lambda_km": [3329.8, 3394.6, 3763.6, 4548.6],
                "r_t_earth_paper6_km": 4150.0,
            },
        }

        out_file = OUTPUTS_DIR / "step_3_5_estimator_decomposition.json"
        with open(out_file, "w") as f:
            json.dump(out, f, indent=2)
        print_status(f"  Results written to {out_file}", "SUCCESS")

        return {"status": "success", "results": str(out_file),
                "n_cod_days": n_used["COD"], "n_grg_days": n_used["GRG"]}


def defaultdict_agg():
    from collections import defaultdict
    return defaultdict(_Agg)


def main():
    step = Step35EstimatorDecomposition()
    result = step.run()
    print(json.dumps({"status": result["status"]}, indent=2))


if __name__ == "__main__":
    main()
