# Global Time Echoes: MGEX Multi-GNSS Clock Replication, 2025–2026

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20572726.svg)](https://doi.org/10.5281/zenodo.20572726)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

![Global Time Echoes: MGEX Multi-GNSS Clock Replication, 2025–2026](site/public/image.webp)

**Author:** Matthew Lukin Smawfield  
**Version:** v0.3 (Suva)  
**First published:** 7 June 2026 · **Last updated:** 30 September 2026  
**Status:** Preprint  
**DOI:** [10.5281/zenodo.20572726](https://doi.org/10.5281/zenodo.20572726)  
**Website:** [https://mlsmawfield.com/tep/gnss-mgex](https://mlsmawfield.com/tep/gnss-mgex)  
**Paper Series:** TEP Series: Paper 14 (MGEX Multi-GNSS)

## Abstract





The Temporal Equivalence Principle (TEP) predicts correlated phase-coherent disturbances in GNSS timekeeping with a spatial correlation length of order thousands of kilometres, an east–west anisotropy, coupling to Earth orbital velocity, and an annual-phase direction near Earth's aphelion velocity tangent. The directed CMB-apex template is anti-phase rather than detected as a rest frame. Papers 1 and 2 examined these signatures in GPS precise point-positioning products from three analysis centres spanning 2000–2025; Paper 3 used raw RINEX single-point positioning. This paper is a product-type independence stress test using the public CODE COD0MGXFIN combined multi-GNSS receiver-clock solution distributed by NASA CDDIS from 2025-01-01 to 2026-05-01. Receiver-clock offsets for 418 globally distributed stations are read directly from daily 5-minute CLK files. The combined solution provides no per-constellation independence, but supplies a distinct product family and a largely held-out 2025–2026 epoch. The isotropic correlation length is λ = 1862 km (R² = 0.699; 1.83 million pair-day observations) with a formal within-fit error of ± 112 km and a cluster-robust statistical uncertainty of ± 155 km obtained by resampling whole epochs (month-block bootstrap over 17 calendar months; day-block ± 86 km; leave-one-month-out jackknife ± 94 km). The scale is epoch- and product-dependent at the ~30–40% level: single-window realisations span 1,050–1,935 km, and matched-day estimates are 2,744 km for CODE against 1,933 km for the independently processed GRG product chain (same tracking network) — the latter landing within 4% of the pooled scale — so the defensible content is the thousands-of-kilometres coherence rather than the precise value. It is stable at 1,819–2,000 km across minimum-baseline cutoffs of 50–2,000 km and persists on geomagnetically quiet days ( λ = 1708 km). The full-range anisotropy is not robust (ratio 1.29, station-clustered p = 0.49); the longitude-matched subset is suggestive by pair bootstrap but not spatially clustered resampling. The primary monthly λ and EW/NS orbital metrics are null (Bonferroni p > 0.5). A model-free PA-difference statistic gives associations with orbital speed (r = −0.716, p = 0.0018) and CMB-dipole velocity projection (r = −0.767, p = 0.0005) across 16 chronological full months. Autocorrelation-aware corrections temper but do not erase the association: the effective-degrees-of-freedom values are p = 0.072 (speed) and p = 0.018 (v_proj), matched AR(1) surrogates give p = 0.019 and p = 0.003, the phase-preserving circular-shift null reaches its 1/16 resolution floor for the projection channel, and leave-one-month-out correlations stay within r = −0.66 to −0.86, so no single phase alignment dominates. The association also survives operational-asymmetry controls: a composition-matched EW−NS contrast over distance × latitude cells gives r = −0.65 (p = 0.006) for the projection channel, and a uniformly-tracked-station subset gives r = −0.74 (p = 0.001). Hemisphere-split checks return the same signed anticorrelation in northern, southern, and cross-hemisphere pair classes (r = −0.45, −0.88, −0.49) at a common annual phase (Δφ = 22° ± 44°), favouring a global heliocentric driver over hemisphere-local seasonal systematics. The physical hypothesis and sign were motivated by earlier TEP-GNSS papers, making this cross-product corroboration; the operationalisation remains exploratory within MGEX because it was selected after the primary estimators were unstable. A 37 × 19 directional grid finds strong within-grid modulation (200-permutation plus-one p = 0.005). Its nominal maximum is 82.7° from the CMB dipole. Orbital velocity constrains only an axis’s projection onto the orbital plane, and that projection is highly degenerate (the CMB and nominal-maximum projected predictors correlate at r = 0.984, and only 61 of 703 grid directions meet or exceed the direct-CMB |r|). The CMB-labelled projection is therefore one member of a broad family of annual templates, not evidence for a CMB rest frame or a unique three-dimensional axis. All four null controls pass under the robust analysis (minimum baseline 500 km); under the primary short-baseline parameters the solar-rotation reassignment null retains R² = 0.50, traced to a small-sample bin artifact at 100–152 km, so the robust analysis is reported as the clean result (the fitted scale is stable, 1,862 vs 1,819 km). The two leading mundane attribution channels are explicitly modelled and strongly disfavoured by the tested controls rather than assumed absent: a forward-modelled shared satellite-error kernel is flat through the signal band and decays only on the constellation-footprint scale, and the weather-organized partitions contradict the atmospheric signature — cross-hemisphere pairs retain the cleanest decay (λ = 2101 km) while same-latitude-band pairs beyond 8,000 km return −0.118 rather than the positive zonal coherence an atmospheric driver requires. A supplementary satellite-clock analysis does not recover the exponential structure, consistent with the shared reference-time architecture of the combined MGEX product: because satellite clocks are not independently estimated in that product, it tests only the receiver-clock channel, and the satellite-clock null is a product limitation rather than a prediction test.


## Overview

This project tests the public CODE `COD0MGXFIN` combined multi-GNSS receiver-clock product over 2025-01-01 to 2026-05-01. It extends Papers 1–3 to a different clock-product family and a largely held-out epoch. The MGEX clock files are a single combined solution per station, so this is not an independent per-constellation analysis.

The isotropic correlation length is λ = 1862 km (R² = 0.699; 1.83M pair-day observations; 418 stations; formal ± 112 km, cluster-robust ± 155 km) and remains stable over the tested baseline cutoffs. It persists on quiet geomagnetic days; the active-day estimate is too imprecise to support a claimed change. The full-range anisotropy is not robust (EW/NS ratio 1.29, station-clustered p = 0.49).

Primary monthly λ and EW/NS orbital metrics are null after correction. A model-free PA-difference statistic is strongly associated with orbital speed (r = −0.716, p = 0.0018) and CMB-projected orbital velocity (r = −0.767, p = 0.0005) across 16 chronological full months. The expected sign was motivated by previous TEP-GNSS papers, making the result corroborative across products; the statistic remains exploratory within MGEX because it was chosen after the primary monthly fits were unstable.

A 37 × 19 directional grid detects within-grid annual modulation (plus-one permutation p = 0.005 from 200 shuffles). Its nominal maximum, RA = 230°, Dec = +60°, is 82.7° from the CMB dipole. Earth-orbital velocity identifies only an axis’s projection onto the orbital plane, and 61 of 703 directions perform at least as well as the CMB-labelled template. This is annual directional structure, not a CMB-frame detection or a unique three-dimensional localization.

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
| Annual-phase direction | 37 × 19 full-sky template scan with complete-grid permutations and projection-degeneracy audit |
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
| **Paper 17** | [TEP-LLR](https://github.com/matthewsmawfield/TEP-LLR) | Lunar Laser Ranging and the Nordtvedt Effect | [10.5281/zenodo.19446028](https://doi.org/10.5281/zenodo.19446028) |
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
  title={Global Time Echoes: MGEX Multi-GNSS Clock Replication, 2025–2026},
  author={Smawfield, Matthew Lukin},
  journal={Zenodo},
  year={2026},
  doi={10.5281/zenodo.20572726},
  note={Preprint v0.3 (Suva)}
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
