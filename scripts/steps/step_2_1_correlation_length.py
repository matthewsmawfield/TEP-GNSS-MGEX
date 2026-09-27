#!/usr/bin/env python3
"""Step 2.1: Correlation length λ_T fitting.

Frozen prediction: λ_T within 1,000–4,000 km (MGEX-appropriate range).
Test: Fit C(r) = A exp(−r/λ) + C₀ on distance bins for each constellation.
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

import json
import numpy as np
from datetime import date
from scripts.utils.logger import TEPLogger, set_step_logger, print_status
from scripts.utils.config import (
    OUTPUTS_DIR, FIGURES_DIR, LOGS_DIR, CONSTELLATIONS,
    MIN_DISTANCE_KM, MAX_DISTANCE_KM, N_BINS, MIN_BIN_COUNT,
    LAMBDA_MIN_KM, LAMBDA_MAX_KM, SENSITIVITY_MIN_DISTANCES
)
from scripts.utils.coherence import fit_exponential_model, exponential_decay

LOGS_DIR.mkdir(parents=True, exist_ok=True)
logger = TEPLogger("step_2_1", log_file_path=LOGS_DIR / "step_2_1_correlation.log")
set_step_logger(logger)

N_BOOTSTRAP = 300
BOOTSTRAP_SEED = 20260927


def _binned_fit_from_aggregates(counts, sums, sumsq, bin_centers,
                                min_distance_km, max_distance_km, min_bin_count):
    """Refit C(r) = A exp(-r/lam) + C0 from per-bin aggregates.

    Reproduces fit_exponential_model exactly: same log-spaced bin centres,
    same min-count masking, same SEM weighting, same curve_fit call.
    """
    from scipy.optimize import curve_fit

    mask = counts >= min_bin_count
    if np.sum(mask) < 5:
        return None
    n = counts[mask]
    means = sums[mask] / n
    # ddof=1 sample variance from aggregates
    var = (sumsq[mask] - n * means**2) / np.maximum(n - 1, 1)
    stds = np.sqrt(np.maximum(var, 0.0))
    sem = stds / np.sqrt(n)
    good = sem > 0
    fallback = np.mean(sem[good]) if np.any(good) else 1.0
    sem = np.where(sem > 0, sem, fallback)

    x = bin_centers[mask]
    try:
        lambda_guess = (min_distance_km + max_distance_km) / 2.0
        lambda_upper = max_distance_km * 2.0
        popt, _ = curve_fit(
            exponential_decay, x, means,
            p0=[0.5, lambda_guess, 0],
            sigma=sem,
            bounds=([0, min_distance_km, -1], [2, lambda_upper, 1]),
            maxfev=5000
        )
        return float(popt[1])
    except Exception:
        return None


def _binned_profile(records, min_distance_km, max_distance_km, n_bins, min_bin_count):
    """Return (bin_centers, means, sem, counts) for the qualifying log bins."""
    distances = np.array([r["distance_km"] for r in records])
    y = np.array([r["phase_alignment"] for r in records])
    in_range = (distances >= min_distance_km) & (distances < max_distance_km)
    distances, y = distances[in_range], y[in_range]
    bin_edges = np.logspace(np.log10(min_distance_km), np.log10(max_distance_km), n_bins + 1)
    bin_idx = np.clip(np.digitize(distances, bin_edges) - 1, 0, n_bins - 1)
    cnt = np.bincount(bin_idx, minlength=n_bins)
    s = np.bincount(bin_idx, weights=y, minlength=n_bins)
    q = np.bincount(bin_idx, weights=y * y, minlength=n_bins)
    mask = cnt >= min_bin_count
    n = cnt[mask].astype(float)
    centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0
    centers = centers[mask]
    means = s[mask] / n
    var = np.maximum((q[mask] - n * means**2) / np.maximum(n - 1, 1), 0.0)
    sem = np.sqrt(var) / np.sqrt(n)
    sem = np.where(sem > 0, sem, np.mean(sem[sem > 0]) if np.any(sem > 0) else 1.0)
    return centers, means, sem, n


def fit_morphology_comparison(records, min_distance_km, max_distance_km,
                              n_bins, min_bin_count):
    """Compare decay morphologies on identical binned data and SEM weights.

    Candidates: the exponential (primary estimator), a stretched exponential
    with free exponent (nests the exponential and tests whether the data
    depart from it), the algebraic suppressor family
    A/[1+(r/r_t)^k]+C0 with k free and at the derived Paper-0 values k=4
    (pairwise exponent) and k=4/3 (comparable-mass transition slope), and a
    Gaussian. AIC selects among them.
    """
    from scipy.optimize import curve_fit

    x, mean, sem, _ = _binned_profile(records, min_distance_km, max_distance_km,
                                      n_bins, min_bin_count)
    lam_hi = max_distance_km * 2.0
    lam0 = (min_distance_km + max_distance_km) / 2.0

    models = {
        "exponential": (lambda r, A, l, C0: A * np.exp(-r / l) + C0,
                        [0.5, lam0, 0.0], ([0, min_distance_km, -1], [2, lam_hi, 1]),
                        ["A", "lambda_km", "C0"]),
        "stretched_exponential": (lambda r, A, l, b, C0: A * np.exp(-(r / l) ** b) + C0,
                                  [0.8, 3000.0, 0.8, -0.05], ([0, min_distance_km, 0.1, -1], [2, 50000.0, 4.0, 1]),
                                  ["A", "lambda_km", "stretch_b", "C0"]),
        "algebraic_free_k": (lambda r, A, rt, k, C0: A / (1.0 + (r / rt) ** k) + C0,
                             [0.8, 2000.0, 2.0, -0.05], ([0, min_distance_km, 0.1, -1], [2, 50000.0, 8.0, 1]),
                             ["A", "r_t_km", "k", "C0"]),
        "algebraic_k4": (lambda r, A, rt, C0: A / (1.0 + (r / rt) ** 4) + C0,
                         [0.8, 2000.0, -0.05], ([0, min_distance_km, -1], [2, 50000.0, 1]),
                         ["A", "r_t_km", "C0"]),
        "algebraic_k43": (lambda r, A, rt, C0: A / (1.0 + (r / rt) ** (4.0 / 3.0)) + C0,
                          [0.8, 2000.0, -0.05], ([0, min_distance_km, -1], [2, 50000.0, 1]),
                          ["A", "r_t_km", "C0"]),
        "gaussian": (lambda r, A, l, C0: A * np.exp(-(r / l) ** 2) + C0,
                     [0.8, 4000.0, 0.0], ([0, min_distance_km, -1], [2, lam_hi, 1]),
                     ["A", "lambda_km", "C0"]),
    }

    out = {}
    for name, (fn, p0, bounds, pnames) in models.items():
        try:
            popt, _ = curve_fit(fn, x, mean, p0=p0, sigma=sem, bounds=bounds, maxfev=20000)
            pred = fn(x, *popt)
            chi2 = float(np.sum(((mean - pred) / sem) ** 2))
            dof = max(len(x) - len(popt), 1)
            ss_res = float(np.sum((mean - pred) ** 2))
            ss_tot = float(np.sum((mean - np.mean(mean)) ** 2))
            out[name] = {
                "params": {pn: float(pv) for pn, pv in zip(pnames, popt)},
                "chi2": chi2,
                "dof": int(dof),
                "chi2_red": chi2 / dof,
                "aic": chi2 + 2 * len(popt),
                "r_squared": float(1 - ss_res / ss_tot) if ss_tot > 0 else 0.0,
                "success": True,
            }
        except Exception:
            out[name] = {"success": False}
    ok = {k: v for k, v in out.items() if v.get("success")}
    if ok:
        out["aic_winner"] = min(ok, key=lambda k: ok[k]["aic"])
    out["n_bins"] = int(len(x))
    return out


def cluster_robust_lambda_uncertainty(records, min_distance_km, max_distance_km,
                                      n_bins, min_bin_count,
                                      n_bootstrap=N_BOOTSTRAP, seed=BOOTSTRAP_SEED):
    """Cluster-robust sigma_lambda via day/month block bootstrap + month jackknife.

    Pair-day records are not independent: they share stations within each day
    and share orbit/product states across consecutive days. The formal
    covariance error of the binned fit therefore understates the true scale
    uncertainty. Resampling whole epochs (days or calendar months) preserves
    all within-epoch clustering and yields the honest epoch-level spread.
    """
    distances = np.array([r["distance_km"] for r in records])
    y = np.array([r["phase_alignment"] for r in records])
    days = np.array([r["day"] for r in records])

    in_range = (distances >= min_distance_km) & (distances < max_distance_km)
    distances, y, days = distances[in_range], y[in_range], days[in_range]

    bin_edges = np.logspace(np.log10(min_distance_km), np.log10(max_distance_km), n_bins + 1)
    bin_idx = np.clip(np.digitize(distances, bin_edges) - 1, 0, n_bins - 1)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0

    unique_days, day_of = np.unique(days, return_inverse=True)
    n_days = len(unique_days)

    # day x bin aggregates
    flat = day_of * n_bins + bin_idx
    C_db = np.bincount(flat, minlength=n_days * n_bins).reshape(n_days, n_bins)
    S_db = np.bincount(flat, weights=y, minlength=n_days * n_bins).reshape(n_days, n_bins)
    Q_db = np.bincount(flat, weights=y * y, minlength=n_days * n_bins).reshape(n_days, n_bins)

    def fit_with_multiplicity(m):
        with np.errstate(all='ignore'):
            return _binned_fit_from_aggregates(
                m @ C_db, m @ S_db, m @ Q_db, bin_centers,
                min_distance_km, max_distance_km, min_bin_count)

    # Calendar-month label per unique day (day strings are YYYYDDD)
    def _month(dstr):
        d = date(int(dstr[:4]), 1, 1).toordinal() + int(dstr[4:]) - 1
        dd = date.fromordinal(d)
        return dd.year * 12 + dd.month

    months_of_day = np.array([_month(d) for d in unique_days])
    unique_months, month_of_day = np.unique(months_of_day, return_inverse=True)
    n_months = len(unique_months)

    rng = np.random.default_rng(seed)
    one = np.ones(n_days)

    day_boot = []
    for _ in range(n_bootstrap):
        m = np.bincount(rng.integers(0, n_days, n_days), minlength=n_days).astype(float)
        lam = fit_with_multiplicity(m)
        if lam is not None:
            day_boot.append(lam)

    month_boot = []
    for _ in range(n_bootstrap):
        mm = np.bincount(rng.integers(0, n_months, n_months), minlength=n_months).astype(float)
        lam = fit_with_multiplicity(mm[month_of_day])
        if lam is not None:
            month_boot.append(lam)

    jack = []
    for j in range(n_months):
        m = one.copy()
        m[month_of_day == j] = 0.0
        lam = fit_with_multiplicity(m)
        if lam is not None:
            jack.append(lam)

    lam_full = fit_with_multiplicity(one)

    out = {
        "n_days": n_days,
        "n_months": n_months,
        "n_bootstrap": n_bootstrap,
        "seed": seed,
        "lambda_reaggregated_km": lam_full,
    }
    if day_boot:
        arr = np.array(day_boot)
        out["day_block_bootstrap"] = {
            "n_fits": len(arr),
            "sigma_lambda_km": float(np.std(arr, ddof=1)),
            "ci95_km": [float(np.percentile(arr, 2.5)), float(np.percentile(arr, 97.5))],
            "median_km": float(np.median(arr)),
        }
    if month_boot:
        arr = np.array(month_boot)
        out["month_block_bootstrap"] = {
            "n_fits": len(arr),
            "sigma_lambda_km": float(np.std(arr, ddof=1)),
            "ci95_km": [float(np.percentile(arr, 2.5)), float(np.percentile(arr, 97.5))],
            "median_km": float(np.median(arr)),
        }
    if len(jack) >= 3:
        arr = np.array(jack)
        out["month_jackknife"] = {
            "n_fits": len(arr),
            "lambda_range_km": [float(arr.min()), float(arr.max())],
            "sigma_lambda_km": float(np.sqrt((len(arr) - 1) / len(arr) * np.sum((arr - arr.mean())**2))),
        }
    return out


class Step21CorrelationLength:
    def run(self):
        print_status("Step 2.1: Correlation Length Fitting", "INFO")
        results = {}

        pair_file = OUTPUTS_DIR / "step_2_0_mgex_pairs.json"
        if not pair_file.exists():
            print_status(f"Pair file not found: {pair_file}", "ERROR")
            results["combined"] = {"status": "no_data", "note": f"Pair file not found: {pair_file}"}
            out_file = OUTPUTS_DIR / "step_2_1_correlation_length.json"
            with open(out_file, 'w') as f:
                json.dump(results, f, indent=2)
            return {"status": "no_data", "results": str(out_file)}

        with open(pair_file) as f:
            records = json.load(f)
        if not records:
            print_status("Empty pair records", "WARNING")
            results["combined"] = {"status": "no_data", "note": "Empty pair records"}
            out_file = OUTPUTS_DIR / "step_2_1_correlation_length.json"
            with open(out_file, 'w') as f:
                json.dump(results, f, indent=2)
            return {"status": "no_data", "results": str(out_file)}

        distances = np.array([r["distance_km"] for r in records])
        # Use phase_alignment (cos of magnitude-weighted CSD phase) — the TEP signature metric.
        # coherence (mean spectral magnitude) is isotropic and lacks directional phase structure.
        coherences = np.array([r["phase_alignment"] for r in records])

        const_name = "combined"
        print_status(f"  Fitting {const_name}...", "INFO")

        fit = fit_exponential_model(
            distances, coherences,
            MIN_DISTANCE_KM, MAX_DISTANCE_KM, N_BINS, MIN_BIN_COUNT
        )

        if fit and fit["success"]:
            # Sensitivity analysis: λ at multiple minimum distances
            sensitivity = {}
            for min_d in SENSITIVITY_MIN_DISTANCES:
                s_fit = fit_exponential_model(
                    distances, coherences,
                    min_d, MAX_DISTANCE_KM, N_BINS, MIN_BIN_COUNT
                )
                if s_fit and s_fit.get("success"):
                    sensitivity[f"min_{min_d}km"] = {
                        "lambda_km": s_fit["correlation_length_km"],
                        "lambda_err_km": s_fit["correlation_length_err_km"],
                        "R2": s_fit["r_squared"],
                        "min_distance_km": min_d,
                        # The lambda lower bound equals min_d: flag boundary
                        # contact so a pinned value is not read as a
                        # measurement interior to the allowed range.
                        "at_lambda_lower_bound": bool(
                            s_fit["correlation_length_km"] <= min_d * 1.001
                        ),
                    }
            # Cluster-robust uncertainty: pair-days share stations within days
            # and orbit/product states across days, so the formal covariance
            # error understates the real scale scatter. Resample whole epochs.
            print_status("    Computing cluster-robust λ uncertainty (day/month block bootstrap, month jackknife)...", "INFO")
            cluster = cluster_robust_lambda_uncertainty(
                records, MIN_DISTANCE_KM, MAX_DISTANCE_KM, N_BINS, MIN_BIN_COUNT
            )
            lam_reagg = cluster.get("lambda_reaggregated_km")
            if lam_reagg is not None and abs(lam_reagg - fit["correlation_length_km"]) > 1.0:
                print_status(f"    WARNING: re-aggregated fit λ = {lam_reagg:.1f} km differs from primary fit ({fit['correlation_length_km']:.1f} km)", "WARNING")
            for key in ("day_block_bootstrap", "month_block_bootstrap", "month_jackknife"):
                if key in cluster:
                    print_status(f"      {key}: σ_λ = {cluster[key]['sigma_lambda_km']:.0f} km (n_fits={cluster[key]['n_fits']})", "INFO")

            # Morphology comparison on identical bins/weights: the chosen
            # exponential against the derived algebraic suppressor family
            # (k free, k=4, k=4/3), a stretched exponential (nests the
            # exponential via free b), and a Gaussian.
            morph = fit_morphology_comparison(
                records, MIN_DISTANCE_KM, MAX_DISTANCE_KM, N_BINS, MIN_BIN_COUNT
            )
            if morph.get("aic_winner"):
                for mname, mres in morph.items():
                    if isinstance(mres, dict) and mres.get("success"):
                        print_status(
                            f"      morphology[{mname}]: chi2={mres['chi2']:.1f} dof={mres['dof']} "
                            f"AIC={mres['aic']:.1f} R2={mres['r_squared']:.3f} params={mres['params']}",
                            "INFO"
                        )
                print_status(f"      morphology AIC winner: {morph['aic_winner']}", "INFO")

            results[const_name] = {
                "lambda_km": fit["correlation_length_km"],
                "lambda_err_km": fit["correlation_length_err_km"],
                "A": fit["amplitude"],
                "C0": fit["offset"],
                "R2": fit["r_squared"],
                "sensitivity_analysis": sensitivity,
                "cluster_uncertainty": cluster,
                "morphology_comparison": morph,
                "status": "success",
                "in_range": bool(LAMBDA_MIN_KM <= fit["correlation_length_km"] <= LAMBDA_MAX_KM),
                "n_pairs": len(records)
            }
            print_status(f"    {const_name}: λ = {fit['correlation_length_km']:.0f} ± {fit['correlation_length_err_km']:.0f} km, R² = {fit['r_squared']:.3f}", "INFO")
            print_status(f"    Sensitivity (λ vs min distance):", "INFO")
            for key, vals in sensitivity.items():
                print_status(
                    f"      {key}: λ={vals['lambda_km']:.0f} km, R²={vals['R2']:.3f}",
                    "INFO"
                )
        else:
            results[const_name] = {
                "status": "fit_failed",
                "note": "Exponential fit did not converge",
                "n_pairs": len(records)
            }

        out_file = OUTPUTS_DIR / "step_2_1_correlation_length.json"
        with open(out_file, 'w') as f:
            json.dump(results, f, indent=2)
        print_status(f"Correlation length results written to {out_file}", "INFO")
        return {"status": "success", "results": str(out_file)}
