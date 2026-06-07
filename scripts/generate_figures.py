#!/usr/bin/env python3
"""TEP-GNSS-MGEX Figure Generation
=================================
Generates publication-quality figures from pipeline JSON outputs.

Usage:
    python scripts/generate_figures.py
    python scripts/generate_figures.py --constellation GPS

Author: Matthew Lukin Smawfield
Date: May 2026
License: CC-BY-4.0
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

import json
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scripts.utils.logger import TEPLogger, print_status
from scripts.utils.config import OUTPUTS_DIR, FIGURES_DIR, CONSTELLATIONS, LAMBDA_MIN_KM, LAMBDA_MAX_KM

FIGURES_DIR.mkdir(parents=True, exist_ok=True)


def load_json(filename):
    path = OUTPUTS_DIR / filename
    if not path.exists():
        return {}
    with open(path) as f:
        return json.load(f)


def plot_correlation_length(constellations=None):
    """Figure 1: Correlation length."""
    data = load_json("step_2_1_correlation_length.json")
    consts = constellations or (list(CONSTELLATIONS.keys()) if any(c in data for c in CONSTELLATIONS) else ["combined"])

    fig, ax = plt.subplots(figsize=(8, 5))
    lambdas = []
    labels = []
    for c in consts:
        if c in data and data[c].get("lambda_km") is not None:
            lambdas.append(data[c]["lambda_km"])
            labels.append(c)

    if lambdas:
        x = np.arange(len(labels))
        ax.bar(x, lambdas, color="#2c3e50", alpha=0.8)
        ax.axhline(LAMBDA_MIN_KM, color="green", linestyle="--", alpha=0.5, label="Predicted range")
        ax.axhline(LAMBDA_MAX_KM, color="green", linestyle="--", alpha=0.5)
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.set_ylabel("Correlation length λ (km)")
        ax.set_title("Figure 1: Correlation length per constellation")
        ax.legend()
    else:
        ax.text(0.5, 0.5, "No data available\nRun step 2.1 first", ha="center", va="center", transform=ax.transAxes)
        ax.set_title("Figure 1: Correlation length per constellation (pending)")

    fig.tight_layout()
    out_path = FIGURES_DIR / "fig_1_correlation_length.png"
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print_status(f"Saved {out_path}")


def plot_ew_ns_anisotropy(constellations=None):
    """Figure 2: EW vs NS correlation length."""
    data = load_json("step_2_2_ew_ns_anisotropy.json")
    consts = constellations or (list(CONSTELLATIONS.keys()) if any(c in data for c in CONSTELLATIONS) else ["combined"])

    fig, ax = plt.subplots(figsize=(8, 5))
    ew_vals = []
    ns_vals = []
    labels = []
    for c in consts:
        if c in data and data[c].get("ew_lambda_km") is not None:
            ew_vals.append(data[c]["ew_lambda_km"])
            ns_vals.append(data[c]["ns_lambda_km"])
            labels.append(c)

    if ew_vals:
        x = np.arange(len(labels))
        width = 0.35
        ax.bar(x - width/2, ew_vals, width, label="EW", color="#2980b9")
        ax.bar(x + width/2, ns_vals, width, label="NS", color="#e74c3c")
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.set_ylabel("Correlation length λ (km)")
        ax.set_title("Figure 2: EW vs NS correlation length")
        ax.legend()
    else:
        ax.text(0.5, 0.5, "No data available\nRun step 2.2 first", ha="center", va="center", transform=ax.transAxes)
        ax.set_title("Figure 2: EW vs NS correlation length (pending)")

    fig.tight_layout()
    out_path = FIGURES_DIR / "fig_2_ew_ns_anisotropy.png"
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print_status(f"Saved {out_path}")


def plot_orbital_coupling(constellations=None):
    """Figure 3: Orbital-velocity coupling.

    Displays all three primary orbital metrics:
      - λ vs orbital speed (original Paper 1-3 metric)
      - EW/NS ratio vs orbital speed
      - PA_diff vs orbital speed (short-baseline robust supplementary metric)
    """
    data = load_json("step_2_3_orbital_coupling.json")
    consts = constellations or (list(CONSTELLATIONS.keys()) if any(c in data for c in CONSTELLATIONS) else ["combined"])

    fig, ax = plt.subplots(figsize=(9, 5))

    metrics = [
        ("λ vs v_orb", "lambda_vs_speed_r"),
        ("EW/NS ratio vs v_orb", "ew_ns_ratio_vs_speed_r"),
        ("PA_diff vs v_orb", "pa_diff_vs_speed_r"),
    ]
    n_metrics = len(metrics)
    width = 0.25
    has_data = False

    for i, c in enumerate(consts):
        if c not in data:
            continue
        for j, (label, key) in enumerate(metrics):
            r = data[c].get(key)
            if r is not None and np.isfinite(r):
                has_data = True
                x = i + (j - (n_metrics - 1) / 2) * width
                color = "#27ae60" if r > 0 else "#c0392b"
                ax.bar(x, r, width, color=color, alpha=0.8, label=label if i == 0 else "")

    if has_data:
        ax.set_xticks(np.arange(len(consts)))
        ax.set_xticklabels(consts)
        ax.set_ylabel("Pearson r")
        ax.set_title("Figure 3: Orbital-velocity coupling")
        ax.axhline(0, color="black", linewidth=0.5)
        ax.legend(loc="best")
    else:
        ax.text(0.5, 0.5, "No data available\nRun step 2.3 first", ha="center", va="center", transform=ax.transAxes)
        ax.set_title("Figure 3: Orbital-velocity coupling (pending)")

    fig.tight_layout()
    out_path = FIGURES_DIR / "fig_3_orbital_coupling.png"
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print_status(f"Saved {out_path}")


def plot_cmb_alignment(constellations=None):
    """Figure 4: CMB alignment angular separation."""
    data = load_json("step_2_4_cmb_alignment.json")
    consts = constellations or (list(CONSTELLATIONS.keys()) if any(c in data for c in CONSTELLATIONS) else ["combined"])

    fig, ax = plt.subplots(figsize=(8, 5))
    labels = []
    sep_vals = []
    for c in consts:
        if c in data and data[c].get("angular_separation_to_cmb_deg") is not None:
            sep_vals.append(data[c]["angular_separation_to_cmb_deg"])
            labels.append(c)

    if sep_vals:
        x = np.arange(len(labels))
        ax.bar(x, sep_vals, color="#8e44ad", alpha=0.8)
        ax.axhline(30, color="green", linestyle="--", alpha=0.5, label="30° threshold")
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.set_ylabel("Angular separation from CMB dipole (°)")
        ax.set_title("Figure 4: CMB-frame alignment")
        ax.legend()
    else:
        ax.text(0.5, 0.5, "No data available\nRun step 2.4 first", ha="center", va="center", transform=ax.transAxes)
        ax.set_title("Figure 4: CMB-frame alignment (pending)")

    fig.tight_layout()
    out_path = FIGURES_DIR / "fig_4_cmb_alignment.png"
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print_status(f"Saved {out_path}")


def plot_null_tests(constellations=None):
    """Figure 5: Null control test R² values."""
    data = load_json("step_2_7_null_tests.json")
    consts = constellations or (list(CONSTELLATIONS.keys()) if any(c in data for c in CONSTELLATIONS) else ["combined"])

    fig, ax = plt.subplots(figsize=(10, 5))
    tests = ["temporal_shuffle", "spatial_shuffle", "phase_random", "solar_rotation"]
    test_labels = ["Temporal\nshuffle", "Spatial\nshuffle", "Phase\nrandom", "Solar\nrotation"]
    x = np.arange(len(tests))
    width = 0.2

    for i, c in enumerate(consts):
        if c in data and data[c].get("status") == "success":
            r2_vals = []
            for t in tests:
                r2 = data[c].get(t, {}).get("R2")
                r2_vals.append(r2 if r2 is not None else 0)
            ax.bar(x + (i - len(consts)/2 + 0.5) * width, r2_vals, width, label=c, alpha=0.8)

    ax.axhline(0.3, color="red", linestyle="--", alpha=0.5, label="Null threshold")
    ax.set_xticks(x)
    ax.set_xticklabels(test_labels)
    ax.set_ylabel("R²")
    ax.set_title("Figure 5: Null control tests")
    ax.legend()
    fig.tight_layout()
    out_path = FIGURES_DIR / "fig_5_null_tests.png"
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print_status(f"Saved {out_path}")


def plot_cross_constellation():
    """Figure 6: MGEX validation synthesis."""
    data = load_json("step_3_0_cross_constellation.json")
    if not data:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.text(0.5, 0.5, "No data available\nRun step 3.0 first", ha="center", va="center", transform=ax.transAxes)
        ax.set_title("Figure 6: MGEX validation synthesis (pending)")
        fig.tight_layout()
        out_path = FIGURES_DIR / "fig_6_cross_constellation.png"
        fig.savefig(out_path, dpi=300)
        plt.close(fig)
        print_status(f"Saved {out_path}")
        return

    fig, ax = plt.subplots(figsize=(8, 5))
    checks = data.get("single_solution_checks", {})
    categories = ["λ in range", "Iono persist", "Anisotropy", "Orbital", "CMB aligned"]
    values = [
        1 if checks.get("lambda_in_range") else 0,
        1 if checks.get("ionospheric_persistence") else 0,
        1 if checks.get("anisotropy_replicated") else 0,
        1 if checks.get("orbital_coupling_detected") else 0,
        1 if checks.get("cmb_aligned") else 0,
    ]
    colors = ["#27ae60" if v else "#c0392b" for v in values]
    ax.bar(range(len(categories)), values, color=colors, alpha=0.8)
    ax.set_xticks(range(len(categories)))
    ax.set_xticklabels(categories)
    ax.set_ylim(0, 1.2)
    ax.set_ylabel("Pass / Fail")
    ax.set_title("Figure 6: MGEX validation synthesis")
    verdict = data.get("overall_replication_verdict", "unknown")
    n_pass = data.get("checks_passed", 0)
    n_total = data.get("checks_total", len(categories))
    ax.text(0.5, 0.95, f"Verdict: {verdict} ({n_pass}/{n_total})", ha="center", transform=ax.transAxes, fontsize=12, fontweight="bold")
    fig.tight_layout()
    out_path = FIGURES_DIR / "fig_6_cross_constellation.png"
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print_status(f"Saved {out_path}")


def main():
    parser = argparse.ArgumentParser(description="Generate TEP-GNSS-MGEX figures")
    parser.add_argument("--constellation", choices=list(CONSTELLATIONS.keys()), help="Filter to one constellation")
    args = parser.parse_args()

    consts = [args.constellation] if args.constellation else None

    print_status("=" * 50)
    print_status("TEP-GNSS-MGEX Figure Generation")
    print_status("=" * 50)

    plot_correlation_length(consts)
    plot_ew_ns_anisotropy(consts)
    plot_orbital_coupling(consts)
    plot_cmb_alignment(consts)
    plot_null_tests(consts)
    plot_cross_constellation()

    print_status("\nAll figures generated.")


if __name__ == "__main__":
    main()
