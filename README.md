# TEP-GNSS-MGEX: Held-Out Replication in the Public MGEX Combined-Clock Product

## Paper 14: Global Time Echoes IV

**Building on Papers 1–3:**
- Paper 1: GPS-only analysis (3 centers: CODE, IGS, ESA)
- Paper 2: 25-year long-span CODE GPS analysis
- Paper 3: Raw RINEX / SPP consistency test
- **Paper 14 (this project): Held-out replication in the public MGEX combined multi-GNSS receiver-clock product**

## Scientific Objective

This is a **strict held-out replication paper**, not a discovery paper. It analyses the MGEX combined multi-GNSS receiver-clock product (CODE `COD0MGXFIN`) for the held-out window 2025-01-01 to 2026-05-01 — an epoch disjoint from Papers 1–3 and a data product independent of the GPS PPP family. The independence is of **epoch and data product**, not of constellation: MGEX clock files contain a single combined solution per station, so the four systems cannot be separated. Per-constellation and raw-RINEX independence are already established in Papers 1–3.

All predictions are frozen before any data are inspected:

1. **Correlation length**: λ in the thousands-of-km range (pass band 1,500–4,000 km for the phase-alignment metric)
2. **East–west anisotropy**: EW > NS for azimuth-sector matched station pairs
3. **Orbital-velocity coupling**: monthly λ / EW–NS ratio vs Earth orbital speed
4. **CMB-frame alignment**: full-sky anisotropy scan with look-elsewhere correction
5. **Ionospheric independence**: Kp stratification, storm-day exclusion
6. **Network-geometry robustness**: hemisphere-balanced and distance-matched station subsets

## Headline Result

Four of five evaluated signatures recovered (verdict 4/5): λ = 1396 ± 90 km (R² = 0.486, 1.75M pairs, 333 stations); signal persists on quiet days and passes all four null tests. The full-range anisotropy is modest (ratio 1.23, p = 0.48) but the longitude-matched subset shows east–west excess (ratio 2.28, pair-bootstrap p = 0.004, station-clustered p = 0.244). Traditional monthly λ and EW/NS ratio metrics do not correlate with orbital velocity (Bonferroni p > 0.5), but a supplementary PA-difference metric recovers coupling (r = −0.670, p = 0.017). The CMB-frame test detects a significant anisotropy axis (LEE p < 0.0001) at RA = 60°, Dec = −60° that lies 92° from the CMB dipole, consistent with an ionospheric origin rather than a TEP/CMB-frame effect. A supplementary satellite-clock analysis finds no detectable spatial correlation, consistent with the single-reference-time nature of the MGEX combined solution.

## Data Source

**Public MGEX combined multi-GNSS receiver-clock product from NASA CDDIS / IGS:**
- Product: CODE `COD0MGXFIN` daily CLK files (WUM, GRG as day-level fallbacks)
- Observable: combined multi-GNSS receiver-clock offset per station (AR records)
- Sampling: 5-minute (300 s)
- Stations: 333 globally distributed
- Time span: 2025-01-01 to 2026-05-01 (473 usable days)
- No raw RINEX observations and no positioning are used

## Methodology

Receiver-clock offsets are read from the AR records of each CLK file, converted to nanoseconds, linearly detrended, and analysed via the magnitude-weighted cross-spectral phase (phase-alignment metric) in the band [10 µHz, 500 µHz]. Because the product is a single combined solution, the analysis is computed once; per-constellation entries are identical by construction and retained only for reporting compatibility.

| Frozen prediction | Test |
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

## Abstract

The Temporal Equivalence Principle (TEP) predicts correlated phase-coherent disturbances in GNSS timekeeping with a spatial correlation length of order thousands of kilometres, an east–west anisotropy exceeding the north–south counterpart, coupling to Earth orbital velocity, and a preferred axis near the CMB rest frame. Papers 1 and 2 established these signatures in GPS-only precise point positioning (PPP) products from three analysis centres spanning 2000–2025, and Paper 3 reproduced them in raw RINEX single-point positioning (SPP), demonstrating independence from analysis-centre orbit and clock models. This paper presents a held-out replication using a deliberately independent data product: the public MGEX combined multi-GNSS receiver-clock solution (CODE COD0MGXFIN) distributed by NASA CDDIS, for the held-out window 2025-01-01 to 2026-05-01. Receiver-clock offsets for 333 globally distributed stations are read directly from the daily 5-minute CLK files; no positioning is performed. Because MGEX clock files contain a single combined multi-GNSS solution per station, this test provides temporal and product-type independence from Papers 1–3 rather than per-constellation independence, which those papers already address. All six predictions—correlation length, azimuthal anisotropy, orbital-velocity coupling, CMB-frame alignment, ionospheric independence, and geometric robustness—were frozen before the held-out data were inspected. Four of the five evaluated signatures are recovered. The isotropic correlation length is λ = 1396 ± 90 km (R² = 0.486, 1.75 million pairs), shorter than the 3,000–5,000 km reported for GPS PPP and consistent with the different metric and combined-clock product. The signal persists on geomagnetically quiet days and all four null controls return λ consistent with no structure. The full-range anisotropy is modest (ratio 1.23, p = 0.52) but the predicted east–west excess emerges in the longitude-matched subset (ratio 2.28, p &lt; 0.001). The traditional monthly λ and EW/NS ratio do not correlate with orbital velocity (Bonferroni p &gt; 0.5), but a supplementary PA-difference metric that avoids exponential-fit noise recovers coupling (r = −0.670, p = 0.017). The CMB-frame test detects a significant anisotropy axis (LEE p &lt; 0.0001) at RA = 60°, Dec = −60° that lies 92° from the CMB dipole, consistent with an ionospheric origin. A supplementary satellite-clock analysis using SP3 orbit geometry finds no detectable spatial correlation, consistent with the single-reference-time nature of the MGEX combined solution. Frozen predictions and outcomes (MGEX combined-clock product, held out): Correlation length: λ in the thousands-of-km range — recovered, λ = 1396 ± 323 km (R² = 0.501). Anisotropy: EW > NS — full-range ratio 1.23 (p = 0.52, not significant under station-clustered bootstrap); longitude-matched subset ratio 2.28 (pair-bootstrap p &lt; 0.001, station-clustered p = 0.292). Orbital coupling: monthly λ and EW/NS ratio vs Earth orbital velocity — not recovered via traditional metrics (Bonferroni p &gt; 0.5); supplementary PA-difference metric recovers coupling (r = −0.670, p = 0.017). CMB alignment: full-sky axis scan with look-elsewhere correction — anisotropy axis detected (LEE p &lt; 0.0001) at RA = 60°, Dec = −60° that lies 92° from the CMB dipole, consistent with an ionospheric rather than TEP origin. Ionospheric independence: Kp stratification and storm-day exclusion — signal persists on quiet days. Geometry robustness: hemisphere-balanced and distance-matched subsets — λ stable under distance matching.


See `CITATION.cff` and `CITATION.bib`.
