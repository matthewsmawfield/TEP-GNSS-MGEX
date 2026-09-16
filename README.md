# Global Time Echoes: MGEX Multi-GNSS Clock Replication, 2025–2026

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20572726.svg)](https://doi.org/10.5281/zenodo.20572726)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

![Global Time Echoes: MGEX Multi-GNSS Clock Replication, 2025–2026](site/public/image.webp)

**Author:** Matthew Lukin Smawfield  
**Version:** v0.2 (Suva)  
**First published:** 7 June 2026 · **Last updated:** 16 September 2026  
**Status:** Preprint  
**DOI:** [10.5281/zenodo.20572726](https://doi.org/10.5281/zenodo.20572726)  
**Website:** [https://mlsmawfield.com/tep/gnss-mgex](https://mlsmawfield.com/tep/gnss-mgex)  
**Paper Series:** TEP Series: Paper 14 (MGEX Multi-GNSS)

## Abstract

The Temporal Equivalence Principle (TEP) predicts correlated phase-coherent disturbances in GNSS timekeeping with a spatial correlation length of order thousands of kilometres, an east–west anisotropy, coupling to Earth orbital velocity, and a preferred axis near the CMB rest frame. Papers 1 and 2 examined these signatures in GPS precise point-positioning products from three analysis centres spanning 2000–2025; Paper 3 used raw RINEX single-point positioning.

This paper is a product-type independence stress test using the public CODE COD0MGXFIN combined multi-GNSS receiver-clock solution distributed by NASA CDDIS from 2025-01-01 to 2026-05-01. Receiver-clock offsets for 418 globally distributed stations are read directly from daily 5-minute CLK files. The combined solution provides no per-constellation independence, but supplies a distinct product family and a largely held-out 2025–2026 epoch.

The isotropic correlation length is λ = 1862 ± 112 km (R² = 0.699; 1.83 million pair-day observations). It is stable at 1,819–2,000 km across minimum-baseline cutoffs of 50–2,000 km and persists on geomagnetically quiet days (λ = 1708 km). The full-range anisotropy is not robust (ratio 1.29, station-clustered p = 0.49); the longitude-matched subset is suggestive by pair bootstrap but not spatially clustered resampling. The primary monthly λ and EW/NS orbital metrics are null (Bonferroni p > 0.5). A model-free PA-difference statistic gives associations with orbital speed (r = −0.716, p = 0.0018) and CMB-dipole velocity projection (r = −0.767, p = 0.0005) across 16 chronological full months. The physical hypothesis and sign were motivated by earlier TEP-GNSS papers, making this cross-product corroboration; the operationalisation remains exploratory within MGEX because it was selected after the primary estimators were unstable. A 37 × 19 directional grid finds strong within-grid modulation (200-permutation plus-one p = 0.005). Its nominal maximum is 82.7° from the CMB dipole, but this does not contradict the CMB result: orbital velocity constrains only an axis’s projection onto the orbital plane, and that projection is highly degenerate. The direct CMB-projection association is therefore consistent but not a unique three-dimensional axis localization. All four null controls pass under a robust analysis. A supplementary satellite-clock analysis does not recover the exponential structure, consistent with the shared reference-time architecture of the combined MGEX product.

## Overview

This project tests the public CODE `COD0MGXFIN` combined multi-GNSS receiver-clock product over 2025-01-01 to 2026-05-01. It extends Papers 1–3 to a different clock-product family and a largely held-out epoch. The MGEX clock files are a single combined solution per station, so this is not an independent per-constellation analysis.

The isotropic correlation length is λ = 1862 ± 112 km (R² = 0.699; 1.83M pair-day observations; 418 stations) and remains stable over the tested baseline cutoffs. It persists on quiet geomagnetic days; the active-day estimate is too imprecise to support a claimed change. The full-range anisotropy is not robust (EW/NS ratio 1.29, station-clustered p = 0.49).

Primary monthly λ and EW/NS orbital metrics are null after correction. A model-free PA-difference statistic is strongly associated with orbital speed (r = −0.716, p = 0.0018) and CMB-projected orbital velocity (r = −0.767, p = 0.0005) across 16 chronological full months. The expected sign was motivated by previous TEP-GNSS papers, making the result corroborative across products; the statistic remains exploratory within MGEX because it was chosen after the primary monthly fits were unstable.

A 37 × 19 directional grid detects strong within-grid modulation (plus-one permutation p = 0.005 from 200 shuffles). Its nominal maximum, RA = 230°, Dec = +60°, is 82.7° from the CMB dipole, but this does not contradict CMB consistency: Earth-orbital velocity identifies only an axis’s projection onto the orbital plane. The direct CMB-projection association is strong, while a unique three-dimensional axis is not identifiable from this design. The supplementary satellite-clock analysis does not recover the exponential structure, consistent with the combined product’s shared reference-time architecture.

## Data and methods

- **Product:** CODE `COD0MGXFIN` daily CLK files, with WUM and GRG fallbacks
- **Observable:** combined multi-GNSS receiver-clock offsets (AR records)
- **Sampling:** 5 minutes
- **Analysis:** detrended offsets; magnitude-weighted cross-spectral phase
  alignment over 10–500 µHz; station-pair distance fits and directional tests

| Signature | Test |
| --- | --- |
| Correlation length | Fit $$C(r)=A\exp(-r/\lambda)+C_0$$ across 40 log-distance bins |
| Directional anisotropy | EW/NS sector comparisons and bootstrap controls |
| Orbital coupling | Chronological monthly primary metrics plus PA-difference sensitivity analysis |
| CMB alignment | 37 × 19 full-sky scan with complete-grid permutations |
| Environmental and geometry controls | Kp stratification, storm exclusion, and distance matching |

## The TEP Research Program

| Paper | Repository | Title | DOI |
|-------|-----------|-------|-----|
| **Paper 0** | [TEP](https://github.com/matthewsmawfield/TEP) | Temporal Equivalence Principle: Dynamic Time & Emergent Light Speed | [10.5281/zenodo.16921911](https://doi.org/10.5281/zenodo.16921911) |
| **Paper 1** | [TEP-GNSS](https://github.com/matthewsmawfield/TEP-GNSS) | Global Time Echoes: Distance-Structured Correlations in GNSS Clocks | [10.5281/zenodo.17127229](https://doi.org/10.5281/zenodo.17127229) |
| **Paper 2** | [TEP-GNSS-II](https://github.com/matthewsmawfield/TEP-GNSS-II) | Global Time Echoes: 25-Year Temporal Evolution | [10.5281/zenodo.17517141](https://doi.org/10.5281/zenodo.17517141) |
| **Paper 3** | [TEP-GNSS-RINEX](https://github.com/matthewsmawfield/TEP-GNSS-RINEX) | Global Time Echoes: Raw RINEX Validation of Distance-Structured Correlations in GNSS Clocks | [10.5281/zenodo.17860166](https://doi.org/10.5281/zenodo.17860166) |
| **Paper 4** | [TEP-GL](https://github.com/matthewsmawfield/TEP-GL) | Temporal-Spatial Coupling in Gravitational Lensing: A Reinterpretation of Dark Matter Observations | [10.5281/zenodo.17982540](https://doi.org/10.5281/zenodo.17982540) |
| **Paper 5** | [TEP-GTE](https://github.com/matthewsmawfield/TEP-GTE) | Global Time Echoes: Empirical Validation of the Temporal Equivalence Principle | [10.5281/zenodo.18004832](https://doi.org/10.5281/zenodo.18004832) |
| **Paper 6** | [TEP-UCD](https://github.com/matthewsmawfield/TEP-UCD) | Universal Critical Density: Unifying Atomic, Galactic, and Compact Object Scales | [10.5281/zenodo.18064365](https://doi.org/10.5281/zenodo.18064365) |
| **Paper 7** | [TEP-RBH](https://github.com/matthewsmawfield/TEP-RBH) | The Soliton Wake: A Runaway Black Hole as a Gravitational Soliton | [10.5281/zenodo.18059250](https://doi.org/10.5281/zenodo.18059250) |
| **Paper 8** | [TEP-SLR](https://github.com/matthewsmawfield/TEP-SLR) | Global Time Echoes: Optical-Domain Consistency Test via Satellite Laser Ranging | [10.5281/zenodo.18064581](https://doi.org/10.5281/zenodo.18064581) |
| **Paper 9** | [TEP-EXP](https://github.com/matthewsmawfield/TEP-EXP) | What Do Precision Tests of General Relativity Actually Measure? | [10.5281/zenodo.18109760](https://doi.org/10.5281/zenodo.18109760) |
| **Paper 10** | [TEP-COS](https://github.com/matthewsmawfield/TEP-COS) | Suppressed Density Scaling in Globular Cluster Pulsars | [10.5281/zenodo.18165798](https://doi.org/10.5281/zenodo.18165798) |
| **Paper 11** | [TEP-H0](https://github.com/matthewsmawfield/TEP-H0) | The Cepheid Bias: Resolving the Hubble Tension | [10.5281/zenodo.18209702](https://doi.org/10.5281/zenodo.18209702) |
| **Paper 12** | [TEP-JWST](https://github.com/matthewsmawfield/TEP-JWST) | A Unified Resolution to the JWST High-Redshift Anomalies | [10.5281/zenodo.19000827](https://doi.org/10.5281/zenodo.19000827) |
| **Paper 13** | [TEP-WB](https://github.com/matthewsmawfield/TEP-WB) | Temporal Shear Recovery in Gaia DR3 Wide Binaries | [10.5281/zenodo.19102061](https://doi.org/10.5281/zenodo.19102061) |
| **Paper 14** | **TEP-GNSS-MGEX** (This repo) | MGEX Multi-GNSS Clock Replication | [10.5281/zenodo.20572726](https://doi.org/10.5281/zenodo.20572726) |
| **Paper 15** | [TEP-EFA](https://github.com/matthewsmawfield/TEP-EFA) | Temporal Shear in the Earth Flyby Anomaly | [10.5281/zenodo.19454862](https://doi.org/10.5281/zenodo.19454862) |
| **Paper 16** | [TEP-J0437](https://github.com/matthewsmawfield/TEP-J0437) | Synchronization Holonomy in Pulsar Scintillation | [10.5281/zenodo.19454620](https://doi.org/10.5281/zenodo.19454620) |
| **Paper 17** | [TEP-LLR](https://github.com/matthewsmawfield/TEP-LLR) | Lunar Laser Ranging and the Nordtvedt Effect | [10.5281/zenodo.19446029](https://doi.org/10.5281/zenodo.19446029) |
| **Paper 18** | [TEP-HC](https://github.com/matthewsmawfield/TEP-HC) | EFT Mapping and Acoustic Peak Constraints via hi_class | [10.5281/zenodo.20572722](https://doi.org/10.5281/zenodo.20572722) |
| **Paper 19** | [TEP-LENS](https://github.com/matthewsmawfield/TEP-LENS) | Blind-Prediction Residual Test in Multiply-Imaged Supernovae | — |
| **Paper 26** | [TEP-C0](https://github.com/matthewsmawfield/TEP-C0) | A Covariant Alternative to Cosmic Expansion | [10.5281/zenodo.20370143](https://doi.org/10.5281/zenodo.20370143) |
| **Paper 28** | [TEP-BH](https://github.com/matthewsmawfield/TEP-BH) | Black Holes and the Temporal Horizon | [10.5281/zenodo.21677826](https://doi.org/10.5281/zenodo.21677826) |

## Repository Structure

```text
TEP-GNSS-MGEX/
├── data/                    # Raw and processed MGEX clock data (not versioned)
├── logs/                    # Execution logs
├── results/
│   ├── outputs/             # Derived JSON results
│   └── figures/             # Generated figures
├── scripts/
│   ├── steps/               # Sequential analysis pipeline
│   └── utils/               # Shared configuration and utilities
├── site/
│   └── components/          # HTML source of truth for manuscript
└── README.md
```

## Reproduction Pipeline

To reproduce the analysis and generate the manuscript:

1. **Install Dependencies:**
   ```bash
   pip install -r requirements.txt
   npm install --prefix site
   ```

2. **Run Analysis Pipeline:**
   ```bash
   python scripts/run_all.py
   python scripts/generate_figures.py
   ```

3. **Build Manuscript and Site:**
   ```bash
   npm run build --prefix site
   ```

Downloading from CDDIS requires authenticated access. Raw data, processed data, and generated analysis outputs are excluded from version control.

## Citation

```bibtex
@article{smawfield2026tepgnssmgex,
  title={Global Time Echoes: MGEX Multi-GNSS Clock Replication, 2025--2026},
  author={Smawfield, Matthew Lukin},
  journal={Zenodo},
  year={2026},
  doi={10.5281/zenodo.20572726},
  note={Preprint v0.2 (Suva)}
}
```

See `CITATION.cff` and `CITATION.bib`.

## License

Creative Commons Attribution 4.0 International (CC-BY-4.0).

## Open Science Statement

This is an open research preprint repository. Manuscript sources, pipeline code, and derived outputs are provided to support transparent inspection and independent reproduction.

## Contact

Email: matthewsmawfield@gmail.com  
ORCID: [0009-0003-8219-3159](https://orcid.org/0009-0003-8219-3159)
