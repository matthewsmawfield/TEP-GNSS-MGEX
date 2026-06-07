# TEP-GNSS-MGEX: MGEX Multi-GNSS Clock Replication

## Paper 14: Global Time Echoes

**Building on Papers 1–3:**
- Paper 1: GPS-only analysis (3 centers: CODE, IGS, ESA)
- Paper 2: 25-year long-span CODE GPS analysis
- Paper 3: Raw RINEX / SPP consistency test
- **Paper 14 (this project): Analysis of the public MGEX combined multi-GNSS receiver-clock product**

## Scientific Objective

This paper analyses the MGEX combined multi-GNSS receiver-clock product (CODE `COD0MGXFIN`) for the window 2025-01-01 to 2026-05-01 — an epoch disjoint from Papers 1–3 and a data product independent of the GPS PPP family. The independence is of **epoch and data product**, not of constellation: MGEX clock files contain a single combined solution per station, so the four systems cannot be separated. Per-constellation and raw-RINEX independence are already established in Papers 1–3.

The analysis evaluates the same signatures examined in earlier papers:

1. **Correlation length**: λ in the thousands-of-km range
2. **East–west anisotropy**: EW > NS for azimuth-sector matched station pairs
3. **Orbital-velocity coupling**: monthly λ / EW–NS ratio vs Earth orbital speed
4. **CMB-frame alignment**: full-sky anisotropy scan with look-elsewhere correction
5. **Ionospheric independence**: Kp stratification, storm-day exclusion
6. **Network-geometry robustness**: hemisphere-balanced and distance-matched station subsets

## Headline Result

Results: λ = 1396 ± 90 km (R² = 0.486, 1.75M pairs, 256 stations); signal persists on quiet days (Kp ≤ 2) and strengthens during active conditions (Kp ≥ 5, smaller sample); passes all four null tests. The full-range anisotropy is modest (ratio 1.23, p = 0.48) but the longitude-matched subset shows east–west excess (ratio 2.28, pair-bootstrap p = 0.002, station-clustered p = 0.244). Traditional monthly λ and EW/NS ratio metrics do not correlate with orbital velocity (Bonferroni p > 0.5), but a supplementary PA-difference metric shows coupling (r = −0.670, p = 0.017). The CMB-frame test detects a significant anisotropy axis (LEE p < 0.0001) at RA = 60°, Dec = −60° that lies 92° from the CMB dipole, consistent with an ionospheric origin rather than a TEP/CMB-frame effect. A supplementary satellite-clock analysis finds the exponential model actively rejected (R² < 0 for all constellations), consistent with the single-reference-time nature of the MGEX combined solution.

## Data Source

**Public MGEX combined multi-GNSS receiver-clock product from NASA CDDIS / IGS:**
- Product: CODE `COD0MGXFIN` daily CLK files (WUM, GRG as day-level fallbacks)
- Observable: combined multi-GNSS receiver-clock offset per station (AR records)
- Sampling: 5-minute (300 s)
- Stations: 256 globally distributed
- Time span: 2025-01-01 to 2026-05-01 (473 usable days)
- No raw RINEX observations and no positioning are used

## Methodology

Receiver-clock offsets are read from the AR records of each CLK file, converted to nanoseconds, linearly detrended, and analysed via the magnitude-weighted cross-spectral phase (phase-alignment metric) in the band [10 µHz, 500 µHz]. The product is a single combined solution, so the analysis is computed once.

| Signature | Test |
|---|---|
| λ in thousands of km | Fit C(r) = A exp(−r/λ) + C₀ on 40 log distance bins (50–13,000 km) |
| EW > NS | Azimuth-sector matched pairs; longitude-matched subset; 500× bootstrap |
| Orbital-velocity coupling | Monthly λ / EW–NS ratio vs Earth orbital velocity projected onto CMB dipole (Bonferroni-corrected); supplementary PA-difference metric (monthly mean EW − NS phase alignment) |
| CMB-frame alignment | Full-sky 37×19 axis scan with look-elsewhere correction |
| Not ionospheric | Kp stratification, storm-day exclusion |
| Not network geometry | Hemisphere-balanced and distance-matched station subsets |

## Repository Structure

```
TEP-GNSS-MGEX/
├── data/
│   ├── raw/clk/          # MGEX combined CLK files (gitignored)
│   ├── processed/        # NPZ clock time series, station metadata
│   └── external/         # Kp indices, data provenance
├── scripts/
│   ├── steps/
│   │   ├── step_1_0_data_acquisition.py
│   │   ├── step_1_1_station_metadata.py
│   │   ├── step_2_0_mgex_clock_correlation.py
│   │   ├── step_2_1_correlation_length.py
│   │   ├── step_2_2_ew_ns_anisotropy.py
│   │   ├── step_2_3_orbital_coupling.py
│   │   ├── step_2_4_cmb_alignment.py
│   │   ├── step_2_5_ionospheric_control.py
│   │   ├── step_2_6_geometry_control.py
│   │   ├── step_2_7_null_tests.py
│   │   ├── step_2_8_satellite_clock_analysis.py
│   │   └── step_3_0_cross_constellation.py
│   ├── utils/
│   └── run_all.py
├── results/
│   ├── outputs/          # JSON analysis results
│   └── figures/          # PNG/PDF plots
├── site/                 # Manuscript generation
│   ├── components/
│   ├── public/
│   └── build.js
├── logs/
├── requirements.txt
├── CITATION.cff
├── VERSION.json
└── README.md
```

## Quick Start

```bash
pip install -r requirements.txt
python scripts/run_all.py
```

Build manuscript:
```bash
cd site && npm install && npm run build
```

## License

CC-BY-4.0

## Citation

DOI: [10.5281/zenodo.20572727](https://doi.org/10.5281/zenodo.20572727)

## Abstract

The Temporal Equivalence Principle (TEP) predicts correlated phase-coherent disturbances in GNSS timekeeping with a spatial correlation length of order thousands of kilometres, an east–west anisotropy exceeding the north–south counterpart, coupling to Earth orbital velocity, and a preferred axis near the CMB rest frame. Papers 1 and 2 established these signatures in GPS-only precise point positioning (PPP) products from three analysis centres spanning 2000–2025, and Paper 3 reproduced them in raw RINEX single-point positioning (SPP), demonstrating independence from analysis-centre orbit and clock models. This paper presents an analysis of an independent data product: the public MGEX combined multi-GNSS receiver-clock solution (CODE COD0MGXFIN) distributed by NASA CDDIS, for the window 2025-01-01 to 2026-05-01. Receiver-clock offsets for 256 globally distributed stations are read directly from the daily 5-minute CLK files; no positioning is performed. Because MGEX clock files contain a single combined multi-GNSS solution per station, this test provides product-type independence and uses a largely held-out 2025–2026 epoch relative to Papers 1–3, rather than providing per-constellation independence, which requires raw per-system processing and is not available from this combined-clock product. The analysis evaluates the same signatures examined in earlier papers: correlation length, azimuthal anisotropy, orbital-velocity coupling, CMB-frame alignment, ionospheric independence, and geometric robustness. The isotropic correlation length is λ = 1396 ± 90 km (R² = 0.486, 1.75 million pairs), shorter than the 3,000–5,000 km reported for GPS PPP and consistent with the different metric and combined-clock product. The signal persists on geomagnetically quiet days (Kp ≤ 2) and appears to strengthen during active conditions (Kp ≥ 5, smaller sample); all four null controls collapse to negligible structure (R² ≈ 0). The full-range anisotropy is modest (ratio 1.23, p = 0.48) but a suggestive east–west excess emerges in the longitude-matched subset (ratio 2.28, pair-bootstrap p = 0.002, 95% CI [1.20, 2.69]; spatially-clustered resampling p = 0.244, 95% CI [0.16, 14.12]). The primary monthly λ and EW/NS orbital-coupling tests are not significant. A supplementary short-baseline phase-alignment metric suggestively recovers orbital-velocity modulation (r = −0.670, p = 0.017). The CMB-frame test detects a significant anisotropy axis (LEE p < 0.0001) at RA = 60°, Dec = −60° that lies 92° from the CMB dipole, more consistent with ionospheric or product-geometry contamination than with stable CMB-frame alignment. A product-limited satellite-clock analysis using SP3 orbit geometry finds the exponential model actively rejected (R² < 0), consistent with the single-reference-time nature of the MGEX combined solution. Frozen predictions and outcomes (MGEX combined-clock product, exploratory): Correlation length: λ in the thousands-of-km range — recovered, λ = 1396 ± 323 km (R² = 0.501). Anisotropy: EW > NS — full-range ratio 1.23 (p = 0.48, not significant under station-clustered bootstrap); longitude-matched subset ratio 2.28 (pair-bootstrap p &lt; 0.001, station-clustered p = 0.292). Orbital coupling: monthly λ and EW/NS ratio vs Earth orbital velocity — not recovered via traditional metrics (Bonferroni p &gt; 0.5); supplementary PA-difference metric recovers coupling (r = −0.670, p = 0.017). CMB alignment: full-sky axis scan with look-elsewhere correction — anisotropy axis detected (LEE p &lt; 0.0001) at RA = 60°, Dec = −60° that lies 92° from the CMB dipole, consistent with an ionospheric rather than TEP origin. Ionospheric independence: Kp stratification and storm-day exclusion — signal persists on quiet days. Geometry robustness: hemisphere-balanced and distance-matched subsets — λ stable under distance matching.


See `CITATION.cff` and `CITATION.bib`.
