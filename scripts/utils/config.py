#!/usr/bin/env python3
"""TEP-GNSS-MGEX configuration and paths."""

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]

DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"
EXTERNAL_DIR = DATA_DIR / "external"

RESULTS_DIR = PROJECT_ROOT / "results"
OUTPUTS_DIR = RESULTS_DIR / "outputs"
FIGURES_DIR = RESULTS_DIR / "figures"
LOGS_DIR = PROJECT_ROOT / "logs"

# Frozen analysis parameters (must not change after data inspection)
F1_HZ = 1e-5       # 10 µHz (28 h period)
F2_HZ = 5e-4       # 500 µHz (33 min period)
SAMPLING_PERIOD_SEC = 300.0  # MGEX CLK: 5-minute sampling
FS_HZ = 1.0 / SAMPLING_PERIOD_SEC

# Distance binning
MIN_DISTANCE_KM = 50
MAX_DISTANCE_KM = 13000
N_BINS = 40
MIN_BIN_COUNT = 10

# Station-quality filtering in step_2_0 (outlier and flat-line detection)
MAX_STATION_STD_NS = 1000.0   # ~1 µs residual threshold for malformed clocks
MIN_STATION_STD_NS = 1.0      # flat-line threshold for clocks held constant by AC

# Satellite separation upper bound (antipodal MEO/GEO satellites can reach ~40 000 km)
MAX_SATELLITE_DISTANCE_KM = 50000

# Constellation definitions for MGEX pipelines
CONSTELLATIONS = {
    "GPS": {"sys": "G", "freqs": ["L1", "L2"], "orbital_radius_km": 20200},
    "GLONASS": {"sys": "R", "freqs": ["L1", "L2"], "orbital_radius_km": 19100},
    "Galileo": {"sys": "E", "freqs": ["E1", "E5a"], "orbital_radius_km": 23222},
    "BeiDou": {"sys": "C", "freqs": ["B1", "B2"], "orbital_radius_km": 21500},
}

# CMB dipole direction (Planck 2018, J2000)
CMB_DIPOLE_RA_DEG = 167.0
CMB_DIPOLE_DEC_DEG = -7.0
CMB_DIPOLE_SPEED_KM_S = 369.0

# Kp thresholds for ionospheric stratification
KP_QUIET = 2
KP_ACTIVE = 5

# Maximum relative difference between quiet-Kp and base λ for "signal persists"
IONO_PERSISTENCE_THRESHOLD = 0.3

# Look-elsewhere correction factor for full-sky scan
LOOK_ELSEWHERE_FACTOR = 50

# CMB alignment angular threshold (degrees)
CMB_ALIGNMENT_THRESHOLD_DEG = 30

# Significance threshold for all hypothesis tests
SIGNIFICANCE_ALPHA = 0.05

# Minimum physically meaningful EW/NS anisotropy ratio (prevents
# statistically significant but negligible effects at large N)
MIN_ANISOTROPY_RATIO = 1.05

# Expected correlation length range for MGEX phase_alignment data.
# The original TEP-GNSS (Cairo) reported 3000–5000 km using the
# isotropic coherence metric.  MGEX multi-GNSS combined clocks
# and the phase_alignment metric produce shorter λ (~1000–3500 km).
LAMBDA_MIN_KM = 1000
LAMBDA_MAX_KM = 4000

# Minimum epochs per station-day for reliable coherence computation.
# Must be >= 64 because compute_coherence_phase requires n_points >= 64 for
# Welch cross-spectral density with nperseg >= 32.
MIN_EPOCHS = 64

# Bootstrap resampling iterations for anisotropy significance testing.
# 30 was too few for stable p-values/CIs on the headline ratios.
N_BOOTSTRAP = 500

# Minimum record / pair counts for reliable statistical estimates
MIN_PAIRS_FIT = 100           # minimum pairs for a single exponential fit
MIN_PAIRS_ANALYSIS = 1000     # minimum pairs for a full analysis step
MIN_PAIRS_BOOTSTRAP = 1000    # minimum weighted pairs for bootstrap iterations
MIN_MONTHLY_RECORDS = 100     # minimum records for a monthly orbital-coupling slice
MIN_SATELLITE_PAIRS = 1000    # minimum total satellite pairs for analysis
MIN_CONSTELLATION_PAIRS = 100 # minimum per-constellation satellite pairs

# Welch CSD maximum segment length (samples); actual nperseg is min(NPERSEG_MAX, n_points//2).
# Reduced from 256 to 96 for 5-minute MGEX sampling: a 288-point day yields 5 Welch
# segments instead of 3, cutting CSD variance by ~40 % while keeping ~35 µHz bins.
NPERSEG_MAX = 96

# Geometry-control distance-matched random sample size.
# 50 k points needed for stable fitting over 40 log-spaced bins.
GEOMETRY_SAMPLE_SIZE = 50000

# Number of independent orbital-coupling correlations tested (λ-vs-speed and
# EW/NS-ratio-vs-speed). Used for Bonferroni multiple-comparisons correction.
N_ORBITAL_TESTS = 2

# Azimuth sector width for EW/NS anisotropy tests (degrees)
SECTOR_WIDTH_DEG = 45

# Longitude-difference threshold for ionospheric-decorrelation filtering (degrees)
LON_DIFF_THRESHOLD_DEG = 30.0

# CMB full-sky grid dimensions for axis scan
CMB_GRID_RA_N = 37    # 10° steps from 0–360
CMB_GRID_DEC_N = 19   # 10° steps from −90–90

# CMB permutation null-test iterations
CMB_PERMUTATION_N = 200

# CMB analysis: sector widths to test for sensitivity (degrees)
CMB_SECTOR_WIDTHS = [45, 60]

# Minimum days with valid ratios for CMB analysis
MIN_DAYS_CMB = 5

# Hemisphere-balancing strata count for geometry control
HEMISPHERE_STRATA_N = 10

# Fixed random seeds for reproducibility
BOOTSTRAP_SEED = 42
NULL_TEST_SEED = 99

# Coordinate mm→m conversion threshold in CLK headers (metres)
# CLK headers store coords in millimetres; Earth radius ~6,371 km,
# so 100 Mm is a safe upper bound for detecting mm-scale units.
COORD_MM_THRESHOLD_METRES = 100_000_000

# Null-test R² threshold (structure is "negligible" below this)
NULL_TEST_R2_THRESHOLD = 0.3

# MGEX analysis window (YYYY, M, D). Deliberately disjoint from the
# original TEP-GNSS-RINEX training period. Extend END forward only as newly
# finalized CLK products become available; do NOT extend backward.
DATA_START = (2025, 1, 1)
DATA_END = (2026, 5, 1)
