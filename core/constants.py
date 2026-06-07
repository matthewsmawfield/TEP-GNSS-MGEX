#!/usr/bin/env python3
"""
TEP Core Constants
================

Canonical physical and phenomenological parameters for the Temporal Equivalence
Principle (TEP) framework.  All TEP papers should import from this module to ensure
consistency across the corpus.  Human-readable registry: parameter_registry.yaml
in this directory.  Do not duplicate these values in project scripts.

Version: TEP v0.9 (Jakarta)
"""

import numpy as np

VERSION = "0.9"
VERSION_CODENAME = "Jakarta"
VERSION_STRING = f"TEP v{VERSION} ({VERSION_CODENAME})"

# =============================================================================
# PHYSICAL CONSTANTS (CODATA 2018)
# =============================================================================
G_NEWTON = 6.67430e-11          # m^3 kg^-1 s^-2
C_LIGHT = 299792458.0            # m s^-1
M_PLANCK = 2.176434e-8           # kg (Planck mass m_P = sqrt(hbar*c/G))
M_SUN = 1.98847e30               # kg
M_EARTH = 5.972e24               # kg
R_EARTH = 6.371e6                # m
R_SUN = 6.96e8                   # m
MPC_TO_M = 3.08567758e22         # m

# =============================================================================
# TEP UNIVERSAL PARAMETERS
# =============================================================================

# Conformal coupling strength.
# phi is dimensionless (measured in reduced-Planck-mass units:
# phi = phi_tilde / M_pl), so the conformal factor is:
#     A(phi) = exp(beta_A * phi)
# with no further M_pl normalization in the code.
BETA_A = -1.0                 # Dimensionless conformal coupling (locked lab-scale convention)

# Phenomenological screening coefficient in the TEP-SPIN tanh ansatz.
# This is NOT the fundamental conformal coupling (BETA_A). It is a
# calibrated parameter of the density-dependent screening model.
BETA_SPIN = 0.01                 # Dimensionless; Paper 24

# Solar-system PPN bound on conformal coupling from Cassini time-delay test.
BETA_CASSINI_MAX = 0.0034        # Bertotti et al. 2003

# Phenomenological saturation scale for Temporal Topology screening.
# When local proximity approaches rho_c (observationally proxied by density),
# the scalar field saturates and A(phi) -> 1, suppressing TEP effects.
RHO_C = 20.0                     # g cm^-3

# Coherence length for lab-scale scalar field
LAB_COHERENCE_LENGTH_M = 50000.0  # 50 km crustal column

# Reference mass scale for geometric coupling beta_geom
M_REF = 1.0e18                   # kg (threshold mass where phi_mass ~ beta_geom)

# Temporal Topology coherence length — canonical value for all forward analysis
# (25-year multi-center GNSS baseline, Papers 1–2/6). Do not substitute short-run
# verification estimates (e.g. Paper 14 MGEX ~1396 km on ~1 yr span).
SCREENING_LENGTH_KM = 4200.0

# MGEX held-out verification only (Paper 14; TEP-GNSS-MGEX step_2_1). Not used in
# NIST/UCD forward models or FEM dimensional normalization.
LAMBDA_T_MGEX_KM = 1396.19
LAMBDA_T_MGEX_ERR_KM = 90.19
LAMBDA_T_MGEX_R2 = 0.486

# Multi-center GNSS exponential fits (Paper 1; TEP-GNSS step_2_0_correlation_analysis_summary.json)
GNSS_LAMBDA_T_LONGSPAN_CODE_KM = 4201
GNSS_LAMBDA_T_LONGSPAN_CODE_ERR_KM = 1967
GNSS_LAMBDA_T_EXPONENTIAL_BY_CENTER = {
    "CODE": {"lambda_km": 4549, "ci_low_km": 1198, "ci_high_km": 5918},
    "IGS": {"lambda_km": 3764, "ci_low_km": 3197, "ci_high_km": 4871},
    "ESA": {"lambda_km": 3330, "ci_low_km": 2532, "ci_high_km": 3984},
}

# Lab-scale coupling constants (TEP-NIST Paper 21)
# alpha_log sign is fixed by the TEP field equation in the (+,-,-,-) metric
# signature: nabla_mu[K(phi) nabla^mu phi] = -alpha(phi) T with alpha = beta_A/M_Pl < 0.
# For non-relativistic dust T = +rho, the static limit gives nabla^2 phi ~ +|alpha| rho,
# so phi decreases with increasing density: dphi/drho < 0.  Since the
# phenomenological ansatz is phi_rho = alpha_log * ln(rho/rho_c), this requires
# alpha_log < 0.  The magnitude |7.66e-3| was determined from the requirement
# that the TEP model reproduce the correct order of magnitude for laboratory
# metrology shifts.
ALPHA_LOG = -7.66e-3             # Density-sector coupling (negative by field-equation sign)
BETA_GEOM = 1.50e-4              # Mass-sector geometric coupling
