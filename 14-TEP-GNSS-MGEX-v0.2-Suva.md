# Global Time Echoes: MGEX Multi-GNSS Clock Replication, 2025–2026
**Matthew Lukin Smawfield**
Version: v0.2 (Suva)
First published: 7 June 2026 · Last updated: 16 September 2026
DOI: 10.5281/zenodo.20572726

---

## Abstract

The Temporal Equivalence Principle (TEP) predicts correlated phase-coherent
disturbances in GNSS timekeeping with a spatial correlation length of order
thousands of kilometres, an east–west anisotropy, coupling to Earth orbital
velocity, and a preferred axis near the CMB rest frame. Papers 1 and 2 examined
these signatures in GPS precise point-positioning products from three analysis
centres spanning 2000–2025; Paper 3 used raw RINEX single-point positioning.

This paper is a product-type independence stress test using the public CODE
COD0MGXFIN combined multi-GNSS receiver-clock solution distributed by NASA
CDDIS from 2025-01-01 to 2026-05-01. Receiver-clock offsets for 418 globally
distributed stations are read directly from daily 5-minute CLK files. The
combined solution provides no per-constellation independence, but supplies a
distinct product family and a largely held-out 2025–2026 epoch.

The isotropic correlation length is $\lambda = 1862$ km (R² = 0.699;
1.83 million pair-day observations) with a formal within-fit error of ±112 km
and a cluster-robust statistical uncertainty of ±155 km obtained by resampling
whole epochs (month-block bootstrap over 17 calendar months; day-block ±86 km;
leave-one-month-out jackknife ±94 km). The scale is epoch- and
product-dependent at the ~30–40% level: single-window realisations span
1,050–1,935 km, and matched-day estimates are 2,744 km for CODE against
1,933 km for the software-independent GRG chain (same tracking network) — the latter landing within 4% of
the pooled scale — so the defensible content is the
thousands-of-kilometres coherence rather than the precise value. It is stable at 1,819–2,000 km across minimum-baseline
cutoffs of 50–2,000 km and persists on geomagnetically quiet days
($\lambda = 1708$ km). The full-range anisotropy is
not robust (ratio 1.29, station-clustered p = 0.49); the longitude-matched
subset is suggestive by pair bootstrap but not spatially clustered resampling.
The primary monthly $\lambda$ and EW/NS orbital metrics are null (Bonferroni
p > 0.5). A model-free PA-difference statistic gives associations with orbital
speed (r = −0.716, p = 0.0018) and CMB-dipole velocity projection
(r = −0.767, p = 0.0005) across 16 chronological full months.
Autocorrelation-aware corrections temper but do not erase the association:
the effective-degrees-of-freedom values are p = 0.072 (speed) and p = 0.018
(v_proj), matched AR(1) surrogates give p = 0.019 and p = 0.003, the
phase-preserving circular-shift null reaches its 1/16 resolution floor for
the projection channel, and leave-one-month-out correlations stay within
r = −0.66 to −0.86, so no single phase alignment dominates. The association
also survives operational-asymmetry controls: a composition-matched EW−NS
contrast over distance × latitude cells gives r = −0.65 (p = 0.006) for the
projection channel, and a uniformly-tracked-station subset gives r = −0.74
(p = 0.001). Hemisphere-split checks return the same signed
anticorrelation in northern, southern, and cross-hemisphere pair classes
(r = −0.45, −0.88, −0.49) at a common annual phase (Δφ = 22° ± 44°),
favouring a global heliocentric driver over hemisphere-local seasonal
systematics. The physical
hypothesis and sign were motivated by earlier TEP-GNSS papers, making this
cross-product corroboration; the operationalisation remains exploratory within
MGEX because it was selected after the primary estimators were unstable.
A 37 × 19 directional grid finds strong within-grid modulation (200-permutation
plus-one p = 0.005). Its nominal maximum is 82.7° from the CMB dipole, but this
does not contradict the CMB result: orbital velocity constrains only an axis’s
projection onto the orbital plane, and that projection is highly degenerate
(the CMB and nominal-maximum projected predictors correlate at r = 0.984, and
only 61 of 703 grid directions meet or exceed the direct-CMB |r|).
The direct CMB-projection association is therefore a directional consistency
check rather than an independent confirmation, and is not a unique
three-dimensional axis localization. All four null controls pass under the
robust analysis (minimum baseline 500 km); under the primary short-baseline
parameters the solar-rotation reassignment null retains R² = 0.50, traced to a
small-sample bin artifact at 100–152 km, so the robust analysis is reported as
the clean result (the fitted scale is stable, 1,862 vs 1,819 km). The two
leading mundane attribution channels are computed and excluded rather than
assumed absent: a forward-modelled shared satellite-error kernel is flat
through the signal band and decays only on the constellation-footprint scale,
and the weather-organized partitions contradict the atmospheric signature —
cross-hemisphere pairs retain the cleanest decay (λ = 2101 km) while
same-latitude-band pairs beyond 8,000 km return −0.118 rather than the
positive zonal coherence an atmospheric driver requires. A supplementary satellite-clock
analysis does not recover the exponential structure, consistent with the shared
reference-time architecture of the combined MGEX product: because satellite
clocks are not independently estimated in that product, it tests only the
receiver-clock channel, and the satellite-clock null is a product limitation
rather than a prediction test.

*Keywords:* Temporal Equivalence Principle, MGEX, multi-GNSS, receiver clock, correlation length, anisotropy, orbital velocity, ionosphere, CMB frame, temporal shear, temporal topology

## 1. Introduction

Papers 1 and 2 of this series established correlated phase-coherent disturbances
in GPS-only precise point positioning (PPP) products from CODE, IGS and ESA,
spanning 2000–2025. The disturbances exhibit a spatial correlation length of
order thousands of kilometres, an east–west anisotropy exceeding the north–south
counterpart, coupling to Earth orbital velocity, and a preferred axis near the
CMB dipole. Paper 3 moved to raw RINEX single-point positioning (SPP),
demonstrating that the same signatures persist without reliance on
analysis-centre orbit or clock products. Between them, those papers already
provide multi-centre, multi-decade, raw-observation, and GPS-only product-family
independence. The principal forms of robustness they do not yet supply are
consistency across the historical training epoch and consistency in a
combined multi-GNSS clock product family.

This paper (Paper 14) conducts a stress test across those two
dimensions. The analysis uses the public MGEX combined multi-GNSS
receiver-clock solution—the CODE COD0MGXFIN CLK product distributed by NASA
CDDIS—for the window 2025-01-01 to 2026-05-01. The present MGEX analysis uses a
combined multi-GNSS clock product distinct from the GPS-only PPP and raw RINEX
products analysed earlier, and a largely held-out 2025–2026 window, though it is
not fully temporally disjoint from all earlier analyses. The MGEX clock product is
generated by combining GPS, GLONASS, Galileo and BeiDou observations into a
single receiver-clock estimate per station, using software, orbit and combination
strategies that differ from the GPS-only PPP of Papers 1–2. It therefore
tests cross-product consistency using a distinct combined multi-GNSS clock
solution over a largely held-out epoch.
Receiver-clock offsets for 418 globally distributed stations are read directly
from the daily 5-minute CLK files; no positioning is performed.

It is emphasised at the outset what this product tests and what it does not.
Because the MGEX clock files contain one combined multi-GNSS solution per station,
they do not permit an independent per-constellation comparison; that
form of independence is the province of raw per-system processing, which lies
outside the scope of this combined-product analysis. The independence offered
here is temporal (a largely held-out epoch) and product-type (a combined multi-GNSS
clock product line from the same analysis centre, CODE, whose IGS-final products anchor Papers 1–2), not per-constellation.

The empirically derived spatial correlation length $\lambda \approx 1862$ km is consistent with the abstract environmental operator $\mathcal{S}_\Sigma(\mathcal{E})$ defined in the foundational TEP framework. In the near-Earth environment, this correlation length traces the macroscopic geometric signature of the continuous saturation of Temporal Topology, anchoring the differential clock correlations independently of the specific multi-GNSS processing strategy. The precise transfer map from the microscopic screening operator to the empirically fitted correlation length remains to be derived.

## 2. Data and Product Scope

### 2.1 MGEX Combined Clock Product

All data are public and freely available from NASA CDDIS
([cddis.nasa.gov](https://cddis.nasa.gov)) and the International
GNSS Service (IGS). The analysis uses the Multi-GNSS Experiment (MGEX) combined
receiver-clock product COD0MGXFIN produced by the Center for Orbit Determination
in Europe (CODE), with the Wuhan (WUM) and CNES (GRG) MGEX products as fallbacks
on the few days CODE is unavailable (WUM 20 days, GRG 9 days out of 486 calendar
days in the window, i.e. ≈6 % fallback usage). These are daily clock (CLK) files at a
5-minute sampling interval. Each file reports a single combined multi-GNSS
receiver-clock estimate per station, derived by CODE from GPS, GLONASS, Galileo
and BeiDou observations; the combined solution is the quantity analysed. No
positioning is performed and no raw RINEX observation files are used.

### 2.2 Analysis Window and Station Set

The analysis period is 2025-01-01 to 2026-05-01 (462 days with usable data out of 486 calendar days in the window),
yielding approximately 16 months of daily files. The window is chosen to
provide a largely held-out 2025–2026 epoch, though it is not fully temporally
disjoint from Papers 1–2, which span 2000–2025. Receiver-clock offsets are extracted from the AR (receiver-clock)
records of each CLK file for 418 globally distributed stations, converted to
nanoseconds, and assembled into per-station daily time series.

### 2.3 Product Limitation: No Per-Constellation Separation

Unlike Papers 1–3, which can isolate individual systems, the MGEX CLK product
contains one combined receiver clock per station rather than per-constellation
clocks. The four constellations therefore cannot be separated in this product,
and the analysis is performed once on the combined solution. The
independence offered here is temporal (a largely held-out epoch) and product-type (a
combined multi-GNSS clock product line), not
per-constellation and not analysis-centre: COD0MGXFIN is produced by CODE,
the same centre whose IGS-final products anchor Papers 1 and 2.

### 2.4 Data Provenance and Shared Error Budgets

The provenance matrix across the GNSS papers makes explicit which analyses
share an error budget:

| Paper | Data layer | Product / processing chain | Epoch | Stations | Error budget shared with |
| --- | --- | --- | --- | --- | --- |
| 1 | AC precise clocks | CODE final, IGS combined, ESA final (Bernese / IGS combination / NAPEOS) | 2023-01 to 2025-06 | 364 | CODE family shared with Papers 2, 14; AC ecosystem shared with Paper 33's centre arms |
| 2 | AC precise clocks | CODE IGS-final 30 s clocks (Bernese) | 2000-03 to 2025-06 | 474 | Same CODE product family as Papers 1 and 14; epoch superset of Paper 1 |
| 3 | Raw RINEX | Broadcast ephemerides + own RTKLIB SPP | 2022 to 2024 | 539 | Estimation-chain independent; shares the physical tracking network |
| **14 (this study)** | AC combined clocks | CODE COD0MGXFIN (WUM/GRG fallback &asymp;6% of days) | 2025-01 to 2026-05 | 418 | Same CODE centre as Papers 1&ndash;2; epoch largely held out |
| 33 | Mixed | Own RTKLib SPP arms + AC receiver clocks (COD/GFZ/GRG/JPL) | 2022 to 2025 | 5&ndash;160 per arm | AC arms share the IGS/MGEX ecosystem; SPP arms are estimation-independent |

Agreement with Papers 1&ndash;2 is therefore corroboration within a shared
product ecosystem &mdash; a different product line and a largely held-out
epoch over substantially the same tracking network &mdash; rather than
replication on independent data. The GRG and WUM products used for fallback
and matched-day comparison run different processing software on that same
network, so they bound software- and implementation-specific artifacts but
not network-level or datum-level common modes. The estimation-chain-independent
axis in the corpus is Paper 3 (raw RINEX, broadcast ephemerides, own
processing); Paper 33's broadcast-SPP arms provide a second such axis, while
its analysis-centre arms share this ecosystem.

### 2.5 Evaluated Signatures

The analysis evaluates the following signatures, consistent with the approach
in earlier papers: correlation length, azimuthal anisotropy, orbital-velocity
coupling, CMB-frame alignment, ionospheric independence, and geometric robustness.

## 3. Methods

### 3.1 Receiver-Clock Extraction

For each daily CLK file the AR records are parsed, which report the combined
multi-GNSS receiver-clock offset for each station at the 5-minute sampling
interval. The offsets are converted from seconds to nanoseconds and ordered in
time, giving one clock time series per station per day. No positioning is
performed and no coordinate residuals are computed; the receiver clock is the
sole observable.

Screening in TEP is represented at the theory level by the environmental operator
*S*<sub>&Sigma;</sub>(*&Epsilon;*).
Quantities such as
&rho;<sub>T</sub>,
*R*<sub>T</sub>(*M*),
*S*<sub>&oplus;</sub>(*r*),
compactness &Phi;/*c*<sup>2</sup>,
local stellar density,
geometric coherence length,
and channel-specific response coefficients
are domain-specific projections of *&Epsilon;*,
not independent screening mechanisms
and not interchangeable universal thresholds.
Each is an observational transfer model
that parameterizes the same underlying operator
in a regime-appropriate form.

### 3.2 Spectral Pre-processing

Each daily series is linearly detrended to remove the dominant clock drift.
Cross- and auto-spectra are estimated with Welch's method (Hann window,
segment length up to 96 samples). Spectral estimates are restricted to the
nominal TEP-sensitive band [10 µHz, 500 µHz], but with the daily Welch
segmentation used here the effective lowest non-zero resolved frequency is set
by the segment length, approximately 35 µHz for 96 five-minute samples. The
practical analysed band is therefore approximately 35–500 µHz, corresponding to
periods from about 8 hours down to 33 minutes, which lies safely within the
5-minute Nyquist limit.

### 3.3 Spectral Phase-Alignment Metric

For every station pair, the magnitude-weighted circular mean of the
cross-spectral phase is computed across the selected frequency band. The
phase-alignment index is the cosine of the weighted mean phase, ranging from
−1 (anti-phase) to +1 (in-phase). It measures signed mean phase-offset
alignment, not phase concentration; the latter requires the magnitude of the
complex resultant. The index is bounded and saturating — it reports the
directional dominance of the in-phase component of the band-integrated
cross-spectrum rather than the correlated-power fraction — so its fitted
exponential amplitude is a phase-direction measure, not a calibrated
coupling, and is decomposed against a linear coherence estimator in §4.5b.
The established phase-alignment statistic is retained
consistently with the earlier GNSS analyses, and is the quantity carried
through all subsequent tests. Pairs are
binned by great-circle distance into 40 logarithmically spaced bins from 50 km
to 13,000 km, with a minimum of 10 pairs per bin.

### 3.4 Correlation-Length Model

The distance-binned phase alignment is modelled as

$C(r) = A \exp(-r / \lambda) + C_0 \tag{1}$

where A is the amplitude, λ is the correlation length, and C₀ is an incoherent
offset. Fitting uses bounded non-linear least squares (Trust Region Reflective),
weighting each bin by the inverse of its standard error of the mean
(1/SEM, where SEM = σ_bin / sqrt(count)). Fit quality is
quantified by R², and the uncertainty on λ is propagated from the covariance
matrix. Because pair-day records share stations within each day and
orbit/product states across consecutive days, the formal covariance error is
supplemented by cluster-robust resampling of whole epochs: a day-block
bootstrap, a calendar-month-block bootstrap, and a leave-one-month-out
jackknife, each refitting Eq. (1) on the resampled record pool. The choice of
morphology is likewise tested rather than assumed: on identical bins and
weights the exponential is compared by AIC against a stretched exponential
with free exponent, the algebraic family A/[1+(r/r_t)^k]+C₀ with free k
and at the derived Paper-0 pairwise exponents k = 4 and k = 4/3, and a
Gaussian. Sensitivity rows at raised minimum-baseline cutoffs are flagged
when the fitted λ contacts the lower parameter bound (which equals the
cutoff), since a boundary-pinned value is not an interior measurement.

### 3.5 Anisotropy, Orbital, CMB, Ionospheric, and Null Tests

The EW > NS test fits λ separately in eight 45°-wide azimuth sectors and compares
the mean of the east and west sectors with the mean of the north and south
sectors. Because ionospheric decorrelation inverts the anisotropy at long
baselines, a longitude-matched subset is also reported (station-pair longitude
difference below 30°), which, following Paper 3, isolates pairs at similar local
solar time. Significance and confidence intervals for the ratio are obtained from
500 bootstrap resamples. For anisotropy, pair-bootstrap results are reported as sensitivity tests, while station/spatially-clustered resampling is treated as the more conservative inference because station pairs are not independent. The ratio is additionally reported in three distance
strata (0–500, 500–1,000, >1,000 km).

The orbital-coupling test correlates the monthly λ and the monthly EW/NS ratio
with the projection of Earth's orbital velocity onto the CMB dipole direction,
computed from a sinusoidal annual model; the two correlations are
Bonferroni-corrected. Because exponential-fit noise on short baselines can
degrade these metrics, an exploratory PA-difference metric is also computed:
the monthly mean EW minus NS phase alignment, which avoids fitting entirely
and directly measures anisotropy strength. The physical hypothesis and sign
were motivated by earlier TEP-GNSS analyses, but this operationalisation was
selected after the primary MGEX estimators proved unstable; it is therefore
exploratory within MGEX and excluded from confirmatory family-wise inference.
Because the 16-month series is short and its monthly values are plausibly
autocorrelated, the nominal Gaussian p-values for the exploratory
PA-difference associations are accompanied by four robustness diagnostics
computed in the same step: an effective-degrees-of-freedom correction from
the residual lag-1 autocorrelation (n_eff = n(1−ρ̂)/(1+ρ̂), t-test on
n_eff−2), a Monte-Carlo ensemble of AR(1) surrogates matched to the observed
residual autocorrelation, a phase-preserving circular-shift permutation of
the monthly series against the fixed periodic predictor (resolution floor
1/n), and a leave-one-month-out correlation range testing whether a single
phase alignment dominates. Two operational-asymmetry controls further price
the vulnerability of the raw EW − NS contrast to unequal east/west station
populations: a composition-matched difference computed as the equal-weight
mean of cell-level EW − NS contrasts over (distance × |midpoint-latitude|)
cells, and a uniformly-tracked-station subset restricted to pairs whose two
endpoints both sit above the median pair-day participation count; each
variant is correlated with the same two predictors under the same
diagnostics.
Months are retained chronologically (YYYY-MM) and a partial trailing month is
excluded. The CMB-frame test performs a full-sky
grid search for the anisotropy axis on a 37 × 19 grid in right ascension and
declination (10° steps), correlating the monthly EW−NS phase-alignment
difference with the projection of Earth's orbital velocity onto each candidate
axis. The difference metric is used in place of the EW/NS ratio because
phase alignment ranges from −1 to +1, making the ratio unstable when the NS
mean is near zero; the difference is bounded and requires no sign-based
selection filter. Monthly aggregation reduces daily noise, while the daily
analysis is retained as a supplementary diagnostic. The grid search tracks both
the most positive and most negative correlation extrema. Because orbital
velocity samples essentially the ecliptic plane, the scan is interpreted as a
test of a plane-projected annual predictor, not as a unique three-dimensional
axis localization; a nominal sky-angle threshold is reported descriptively
only. A permutation null performs a complete grid search over all 703
directions for each of 200 permutations, accounting for the look-elsewhere
effect directly. Ionospheric controls stratify the fit by
geomagnetic Kp index and exclude storm days; geometry controls repeat the fit on
hemisphere-balanced and distance-matched subsets; a forward-modelled
product-error kernel tests the shared satellite-error projection channel
directly; hemisphere-class and latitude-band within-product partitions test
the weather-organized atmospheric channel; and four null tests (temporal
shuffle, spatial shuffle, phase randomisation, and solar-rotation reassignment)
test whether the recovered structure vanishes under label permutation.

### 3.6 Sensitivity and Robustness Analysis

Because the 418-station network contains dense regional clusters that
disproportionately populate short baselines, a sensitivity analysis is
performed to verify that the reported results are not artifacts of
short-baseline fitting noise. For each primary result the analysis is
repeated at five minimum-baseline cutoffs (50, 200, 500, 1,000, and
2,000 km), tabulating the correlation length, anisotropy ratio, and
CMB correlation as functions of the cutoff. A robust supplementary
analysis uses a 500-km minimum baseline and a 50-pair minimum bin
count, eliminating bins with fewer than 50 pairs where small-sample
fluctuations can produce extreme coherence means. For the anisotropy
ratio a station-declustered point estimate is also computed, weighting
each pair by the inverse square root of the product of its two
endpoints' pair counts so that stations in dense clusters do not
dominate the exponential fit. The robust null tests apply the same
elevated thresholds to all four shuffle controls. These sensitivity
analyses do not replace the pre-specified primary analysis but provide
a conservative supplementary assessment of each result's stability.

## 4. Results

The analysis comprises 1,827,555 station-pair-day observations drawn from a network
of 418 stations over 462 usable days. Because the MGEX product is a single combined solution, the
results are reported once rather than per constellation. The analysis evaluates the core correlation-length signature, directional anisotropy, orbital-velocity coupling, CMB-frame alignment, and ionospheric/geometric robustness controls. The correlation length and ionospheric persistence are observed; the longitude-matched anisotropy is not robust under spatially clustered resampling. The primary pre-specified monthly λ and EW/NS orbital metrics are null (Bonferroni p > 0.5). A model-free PA-difference metric (mean EW minus NS phase alignment) gives a suggestive orbital association across 16 chronological full months that survives autocorrelation-aware corrections for the CMB-projection channel (corrected p = 0.018; AR(1) p = 0.003; circular-shift floor), but it is exploratory within this MGEX study because it was adopted after the primary fits proved unstable. The direct CMB-projected association is consistent with that signal. The full-sky scan detects directional modulation, but orbital velocity constrains only a plane-projected direction, so its nominal sky axis is not uniquely localizable.

### 4.1 Correlation Length

The isotropic fit gives a well-constrained correlation length.

| $\lambda$ (km) | $\sigma_\lambda$ (km, formal) | A | C₀ | R² | Pair-day observations |
| --- | --- | --- | --- | --- | --- |
| 1862 | 112 | 0.827 | −0.065 | 0.699 | 1,827,555 |

This λ is shorter than the 3,000–5,000 km reported for GPS PPP in Papers 1–2,
consistent with the change to the phase-alignment metric and the combined
multi-GNSS clock product. The fitted amplitude A = 0.827 is likewise a
bounded phase-direction index rather than a calibrated correlation
amplitude; its ~7× offset from the IGS-final family is decomposed
estimator-axis by estimator-axis in §4.5b (Step 3.6). Until the transfer map from
$\mathcal{S}_\Sigma(\mathcal{E})$ to $\lambda$ is derived, the fitted
correlation length is therefore a product- and metric-dependent
phenomenological observable: the physical TEP signature carried across
products is the order-of-magnitude (thousands-of-kilometres) spatial
coherence, not the precise numerical scale, which responds to the
observable definition and clock-product construction. The ±112 km uncertainty
is the formal covariance-based fit error, which treats the 1.83 million
pair-day records as independent within-bin samples. Because pair-days share
stations within each day and share orbit/product states across consecutive
days, cluster-robust resampling of whole epochs was used to measure the true
statistical scatter: a day-block bootstrap over the 462 days gives
$\sigma_\lambda = \pm 86$ km (95% CI 1,639–1,932 km), a month-block
bootstrap over the 17 calendar months gives $\sigma_\lambda = \pm 155$ km
(95% CI 1,363–1,971 km), and a leave-one-month-out jackknife gives
$\sigma_\lambda = \pm 94$ km with leave-out values confined to
1,794–1,896 km. The within-product statistical precision is therefore about
±8%, only modestly above the formal error; the dominant scale uncertainties
are the epoch and product systematics quantified below. A sensitivity analysis
demonstrates that the correlation length is stable across minimum-baseline
cutoffs from 50 to 2,000 km, with fit quality improving substantially as the
shortest baselines are excluded (R² = 0.699 at 50 km to approximately 0.92).

| Min. baseline (km) | $\lambda$ (km) | R² |
| --- | --- | --- |
| 50 (primary) | 1862 | 0.699 |
| 200 | 1893 | 0.791 |
| 500 | 1819 | 0.848 |
| 1000 | 1850 | 0.917 |
| 2000 | 2000&dagger; | 0.912 |

&dagger; The 2,000 km cutoff row is a boundary contact: the fitted value
(2000.000 km) sits on the lower parameter bound ($\lambda \geq$ minimum
baseline) rather than interior to it, so it indicates only that the data
prefer $\lambda \geq 2{,}000$ km at that cutoff — it is flagged as pinned
in the archived output and is not an independent measurement of the
scale.

The correlation length varies by less than 8% across the full range of
cutoffs ($\lambda$ = 1,819–2,000 km), while R² improves from 0.70 to
0.92. The stability of the fitted correlation length across these cutoffs
demonstrates that the recovered distance-decay scale is not driven solely by
the shortest-baseline bins.

Two further stability checks bound the epoch and processing dependence of
λ (Step 3.1). First, splitting the 462-day pair record into quarterly
windows returns λ = 1,050–1,935 km (R² from −0.18 to 0.80), with
λ = 1,320 km for 2025 alone and 1,935 km for the 2026 subset: the fitted
scale is epoch-dependent at the ~30–40% level, and any single ~1-year
estimate — including the 1,396 km figure reported for the original
~1-year span — is an unstable realisation rather than a product-fixed
value. Second, on strictly matched days the software-independent CNES (GRG)
processing chain — a different code base on the same tracking network —
returns pooled λ = 1,933 km (R² = 0.58) against 2,744 km
for the CODE product on verified-independent matched days (~44,000 pairs
each; same-product fallback days excluded by epoch-aligned
detrended-residual comparison). The GRG estimate lands within 4% of the
pooled MGEX scale (1,862 km), while the ~30% residual offset to the CODE
chain on identical days is the direct product-conditioning signature
decomposed in §4.5b: the recovered scale is a property of the data rather
than of one center's constellation weighting or software. These numbers
fold in the cadence correction of §4.5b — GRG products are delivered at
30-s sampling and per-pair series are now built by intersecting station
epochs on the shared 300-s grid; under the earlier
length-truncated handling the same diagnostic returned 3,257 km against
2,514 km, so the correction tightens the GRG–pooled agreement from ~75%
to ~4%.
The dominant systematic on the λ scale is therefore epoch/window
composition (~30–40%), with analysis-center processing contributing a
further ~±600 km, both well above the ±155 km cluster-robust
statistical error of the pooled estimate.

The choice of decay morphology is tested rather than assumed. On identical
log-spaced bins and identical SEM weights, six functional forms are compared
by AIC: the exponential (AIC = 1,664.7, $\lambda = 1{,}862$ km), a
stretched exponential with free exponent (AIC = 1,665.2; the fit selects
$b = 0.98$, i.e. the pure exponential is recovered rather than assumed),
the algebraic suppressor family $A/[1+(r/r_t)^k]+C_0$ with free
exponent (AIC = 1,893.8 at $k = 1.78$), the same family at the derived
Paper-0 pairwise exponents $k=4$ (AIC = 4,209.0) and $k=4/3$
(AIC = 2,125.1), and a Gaussian (AIC = 3,355.2). The exponential is the
AIC winner and the nested stretched-exponential family independently
returns $b \approx 1$: the measured profile is a finite-range exponential
correlator — the two-point morphology of a screened spatial field — rather
than an algebraic tail. The derived suppressor operators
$\mathcal{S}_\Sigma$ and $y(s)$ describe the response as a function of
the ambient environment, a different observable; the correlation kernel
measured here is the field's spatial covariance, and it is cleanly
exponential.

### 4.2 EW / NS Anisotropy

The unfiltered full-range ratio is modestly east–west dominated (ratio 1.29)
but does not reach significance under the station-clustered bootstrap
(p = 0.49). A suggestive east–west excess emerges once pairs are
matched in longitude (fitted ratio 1.32, bootstrap mean 1.71, pair-bootstrap p = 0.058), though the
significance dissolves under spatial-clustered resampling (95% CI [0.17, 13.28],
p = 0.214). The 500–1,000 km stratum also shows the expected excess.

| Subset | $\lambda_{EW}$ (km) | $\lambda_{NS}$ (km) | Ratio | 95% CI | p |
| --- | --- | --- | --- | --- | --- |
| Full range | 1810 | 1402 | 1.29 | [0.18, 2.60] (station-boot) | 0.49 |
| Longitude-matched (<30°) | 1753 | 1327 | 1.32 | [0.91, 2.65] (pair-boot, mean 1.71) | 0.058 |

| Distance stratum | EW/NS ratio |
| --- | --- |
| 0–500 km | 0.73 |
| 500–1,000 km | 1.38 |
| >1,000 km | 0.82 |

A sensitivity analysis reveals that the directional ratio is fragile
under conservative analysis choices. The ratio declines from 1.29 at
the primary 50-km cutoff to 1.08 at 500 km and 0.89 at 1,000 km, with
station de-clustering yielding a ratio of 0.99. The anisotropy is
therefore classified as not robust; the directional signal is
interpreted as null in the combined MGEX product, consistent with the
non-significant station-clustered bootstrap (p = 0.49).

| Analysis | $\lambda_{EW}$ (km) | $\lambda_{NS}$ (km) | Ratio |
| --- | --- | --- | --- |
| Primary (min 50 km) | 1810 | 1402 | 1.29 |
| Min 200 km | 2098 | 1384 | 1.52 |
| Min 500 km | 1880 | 1734 | 1.08 |
| Min 1000 km | 2211 | 2479 | 0.89 |
| Min 2000 km | 2496 | 2671 | 0.93 |
| Station de-clustered | — | — | 0.99 |
| Robust (min 500 km, bin 50) | 1878 | 1733 | 1.08 |

### 4.3 Orbital-Velocity Coupling

The traditional monthly λ and EW/NS ratio do not correlate significantly with
Earth's orbital speed or its projection onto the CMB dipole (Bonferroni
p > 0.5 for both). These primary pre-specified metrics are degraded by
exponential-fit parameter degeneracy on the short 16-month baseline, where
bounded non-linear least squares fitting of $A \exp(-r/\lambda) + C_0$
suffers from severe variance inflation with too few monthly points.
The monthly means retain January 2025 through April 2026 as 16 separate
chronological observations; the one-day May fragment is excluded. A model-free
PA-difference metric (monthly mean EW − NS phase alignment) avoids exponential
fitting and has a negative association with orbital speed (r = −0.716,
p = 0.0018) and with the velocity projection onto the CMB dipole (r = −0.767,
p = 0.0005). The nominal p-values assume independent monthly points; under
autocorrelation-aware diagnostics the associations are tempered but not
erased. With effective degrees of freedom from the residual lag-1
autocorrelation (ρ̂ = 0.40 and 0.29), the corrected values are p = 0.072
(speed) and p = 0.018 (v_proj); matched AR(1) surrogates give p = 0.019 and
p = 0.003; the phase-preserving circular-shift null gives p = 0.125 and
p = 0.063 — the 1/16 resolution floor for the CMB-projection channel.
Leave-one-month-out correlations stay within r = −0.66 to −0.86 for both
predictors, so the association is not carried by a single favourable phase
alignment. The physically motivated velocity-projection channel is thus the
one that survives every correction, while the scalar-speed control is
marginal. The physical hypothesis and expected sign were motivated by the
earlier TEP-GNSS papers, which reported comparable orbital associations; this
is therefore corroborative across products. It remains exploratory within the
MGEX analysis because this particular operationalisation was selected after
the primary λ and EW/NS estimators proved unstable. It is not included in the
family-wise correction applied to confirmatory tests and should be tested
prospectively in a longer MGEX series. The raw EW − NS contrast is the
estimator most sensitive to any operational asymmetry between the east and
west station populations (network density, tracking schedules); two controls
price that channel directly. A composition-matched difference — EW − NS of
cell means averaged with equal weight over (distance × |latitude|) cells —
retains the association (r = −0.586, p = 0.017 vs speed; r = −0.650,
p = 0.006 vs v_proj, n = 16), and a uniformly-tracked-station subset (both
endpoints above the median pair-day participation) retains it more strongly
(r = −0.685, p = 0.003; r = −0.739, p = 0.001; corrected p = 0.029 and
AR(1) p = 0.007 for the v_proj channel). The association is therefore not
carried by population or scheduling asymmetry between the sectors.
Two further diagnostics bound the interpretation. First, the fitted-scale
nulls are consistent with estimator dilution rather than signal absence:
the monthly λ sequence scatters with CV ≈ 21% about its mean while
remaining uncorrelated with the PA-diff sequence (r = −0.01), and the
per-month per-sector fits that enter the EW/NS λ ratio carry median
propagated ratio errors of ≈ 27% — and one month (2025-08) has a sector
fit pinned at the λ upper bound, rendering that month's ratio formally
uninformative — while the resulting ratio
sequence is uncorrelated even with the fit-free anisotropy it is meant
to measure (ratio vs PA-diff r = −0.20, p = 0.45) and remains null after
dropping the degenerate month (r = −0.10 vs speed, r = −0.17 vs v_proj).
A velocity-driven modulation at the observed contrast amplitude cannot
propagate through the fitted decay length on a 16-month baseline — the
unfitted contrast is the sensitive channel by construction. Second, a
hemisphere decomposition of the same contrast (Step 3.4) returns the
signed orbital-speed anticorrelation in all three pair classes —
intra-northern (r = −0.45), intra-southern (r = −0.88), and
cross-hemisphere (r = −0.49) — with the northern and southern annual
sinusoid phases consistent with a common value (Δφ = 22° ± 44°). A
hemisphere-local seasonal mechanism would appear in anti-phase between
the hemispheres; the observed common phase instead indicates a global
heliocentric driver.

| Quantity | r | p (raw) | p (corr.) | Interpretation |
| --- | --- | --- | --- | --- |
| λ vs speed | −0.16 | 0.554 | 1.000 | Primary (pre-specified): not significant |
| EW/NS ratio vs speed | 0.14 | 0.610 | 1.000 | Primary (pre-specified): not significant |
| λ vs v_proj | −0.05 | 0.850 | 1.000 | Primary (pre-specified): not significant |
| EW/NS ratio vs v_proj | −0.03 | 0.911 | 1.000 | Primary (pre-specified): not significant |
| PA_diff vs speed | −0.716 | 0.0018 | 0.072 | Exploratory; n_eff = 6.9; AR(1) p = 0.019; shift p = 0.125; LOO r ∈ [−0.84, −0.66] |
| PA_diff vs v_proj | −0.767 | 0.0005 | 0.018 | Exploratory; n_eff = 8.8; AR(1) p = 0.003; shift p = 0.063 (floor); LOO r ∈ [−0.86, −0.72] |
| PA_diff (cell-matched) vs speed | −0.586 | 0.0170 | — | Operational-asymmetry control; composition-matched over distance × |latitude| cells |
| PA_diff (cell-matched) vs v_proj | −0.650 | 0.0064 | — | Operational-asymmetry control; n_eff p = 0.065; AR(1) p = 0.028; shift p = 0.063 (floor) |
| PA_diff (uniform-tracking) vs speed | −0.685 | 0.0034 | — | Operational-asymmetry control; both endpoints above median participation |
| PA_diff (uniform-tracking) vs v_proj | −0.739 | 0.0011 | — | Operational-asymmetry control; n_eff p = 0.029; AR(1) p = 0.007; shift p = 0.063 (floor) |

### 4.4 CMB-Frame Alignment

The CMB analysis uses the same PA-difference observable as Section 4.3 and is
therefore a directional characterization, not an independent test. With 16
chronological full-month observations, the direct CMB-projection correlation is
r = −0.733, computed from daily EW−NS differences averaged equally within
each calendar month. This differs from the −0.767 result in Section 4.3,
which averages all pair observations within each month (so months are
pair-count weighted) before correlating PA difference with the CMB-projected
velocity. The latter is the stronger weighted PA-difference result (p = 0.0005).
A 37 × 19 scan gives a nominal maximum
at RA = 230°, Dec = +60° (r = −0.745), 82.7° from the CMB dipole. This apparent
angular discrepancy cannot reject CMB consistency: Earth’s orbital velocity
samples essentially one plane, so the scan identifies only a direction’s
projection onto that plane, not its full three-dimensional vector. The CMB and
nominal maximum projections correlate at r = 0.984, and 61 of 703 grid
directions have |r| at least as large as the direct CMB value. The 30°
three-dimensional threshold is consequently not an identifiable test in this
geometry. The grid detects directional annual modulation (plus-one permutation
p = 0.005 from 200 complete-grid shuffles), while the pre-specified CMB
projection remains compatible with it. Widening sectors to 60° gives r =
−0.791 at another member of the same projection-degenerate family. The short
series further limits precision, so this is corroborative CMB consistency, not
a unique sky-axis localization.

| Min. baseline (km) | Direct r(CMB) | Months |
| --- | --- | --- |
| 50 (primary) | −0.733 | 16 |
| 200 | −0.731 | 16 |
| 500 | −0.730 | 16 |
| 1000 | −0.723 | 16 |
| 2000 | −0.697 | 16 |

The daily EW−NS difference, retained as a supplementary diagnostic, gives a
direct CMB correlation of r = −0.495 (n = 462). Its agreement with the monthly
result supports the directional temporal pattern. Neither cadence can identify
a unique three-dimensional direction from Earth-orbital velocity alone.

| Cadence | Sector width | Nominal RA (°) | Nominal Dec (°) | r | Δθ_CMB (°) | p (perm.) | Interpretation |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Monthly | 45° | 230 | 60 | −0.745 | 82.7 | 0.005 | Axis not uniquely identifiable; CMB projection consistent |
| Monthly | 60° | 230 | 60 | −0.791 | 82.7 | 0.005 | Axis not uniquely identifiable; CMB projection consistent |
| Daily | 45° | 230 | 60 | −0.505 | 82.7 | — | Axis not uniquely identifiable |

### 4.5 Ionospheric and Geometry Controls

The correlation length persists on geomagnetically quiet days ($\lambda = 1708$ km,
Kp ≤ 2, 712,186 pairs, R² = 0.743). Active days (Kp ≥ 5, 39,093 pairs) give a longer
$\lambda = 1931$ km (R² = 0.295); the smaller sample increases the variance and the value awaits confirmation with additional active-period data. Storm-excluded days (Kp < 5) give
$\lambda = 1858$ km (R² = 0.705, 1,788,462 pairs), consistent with the base value. The distance-matched subset
yields a comparable $\lambda = 1889$ km (R² = 0.599), indicating the signal is not a
product of network geometry. Hemispheric $\lambda$ differs north to south, as expected
from the uneven station distribution; the distance-matched control is the
geometry-robust comparison.

| Control | $\lambda$ (km) | R² |
| --- | --- | --- |
| Base (all days) | 1862 | 0.699 |
| Quiet days (Kp ≤ 2) | 1708 | 0.743 |
| Active days (Kp ≥ 5) | 1931 | 0.295 |
| Storm-excluded (Kp < 5) | 1858 | 0.705 |
| Distance-matched | 1889 | 0.599 |
| Northern hemisphere | 1365 | 0.678 |
| Southern hemisphere | 2353 | 0.038 |

### 4.5a Product-Error and Atmospheric-Channel Controls

Two mundane channels could in principle mimic a distance-ordered clock
correlation; each is evaluated against a computed forward model or a
within-product partition rather than by label-permutation nulls alone.
The leading channel is shared satellite orbit and clock error projecting
into each station's receiver-clock estimate. The forward model (Step 3.2)
constructs the expected kernel directly: the pairwise cosine similarity
of the stations' satellite-weight vectors, built from real CODE SP3
orbits (122 satellites, 300 s cadence, 97 epochs spanning a full day)
under standard elevation weighting and the 513 MGEX station coordinates.
The kernel is flat above 0.93 through the 0–2,000 km band where the
measured coherence carries most of its structure, decays only on the
constellation-footprint scale (0.70 at 3–5,000 km, 0.45 at 5–8,000 km,
0.18 at 8–12,000 km — an effective scale an order of magnitude longer
than the observed $\lambda = 1862$ km), and is everywhere non-negative
whereas the measured profile crosses zero. A shared satellite-error
projection can weight this kernel but cannot steepen it, so the
geometric product-error channel cannot supply the measured decay.

The remaining mundane candidate is the spatially-correlated atmospheric
component — tropospheric delay terms absorbed into the per-station clock
estimates. An atmospheric driver is weather-organized rather than
distance-organized: synoptic fields are hemisphere-confined, and the
large-scale modes (annular modes, jet-stream teleconnections) are
coherent along latitude bands over $\sim 10^4$ km, not across them.
Three within-product controls on true station latitudes (Step 3.3) test
this directly. First, refitting the profile within hemisphere classes,
the cross-hemisphere subset — pairs that cannot share any synoptic
weather — returns the cleanest structure in the archive
($\lambda = 2101$ km, R² = 0.769, 655,068 pair-days), against intra-NH
$\lambda = 1308$ km (R² = 0.653). The intra-SH binned profile likewise
decays monotonically (+0.63 to −0.14); its weak exponential-fit quality
(R² = 0.032) reflects the long-baseline-dominated coverage of the
southern network rather than absent structure. Second, at
$d > 4000$ km, same-latitude-band pairs ($|\Delta\mathrm{lat}|  30°$) return −0.026 — the opposite ordering of
the zonal-teleconnection prediction, and not a distance artifact because
the same-band subset is the closer one (median 7,963 vs 9,129 km).
Third, the decisive cell: same-band pairs at $d > 8000$ km — where a
zonal atmospheric kernel must retain positive coherence — return mean
phase alignment −0.118 on 107,013 pair-days, opposite in sign to the
atmospheric requirement. The weather-organized component is therefore
excluded at the level of its required structure; what remains of the
tropospheric channel is per-station incoherent residual error, which
cannot generate a distance-ordered pair correlation by construction.

| Control subset | $\lambda$ (km) | R² | Pair-days | Mundane prediction |
| --- | --- | --- | --- | --- |
| Intra-NH pairs | 1308 | 0.653 | 850,240 | Structure expected |
| Intra-SH pairs (profile decays +0.63 → −0.14) | 2437 | 0.032 | 322,247 | Structure expected |
| Cross-hemisphere pairs (no shared weather) | 2101 | 0.769 | 655,068 | No structure expected — cleanest fit observed |

| Latitude-band control | Mean PA | Pair-days | Median d (km) | Atmospheric requirement |
| --- | --- | --- | --- | --- |
| Same band, $|\Delta\mathrm{lat}| < 10°$, $d > 4000$ km | −0.096 | 217,430 | 7,963 | Positive zonal coherence — contradicted |
| Different band, $|\Delta\mathrm{lat}| > 30°$, $d > 4000$ km | −0.026 | 899,956 | 9,129 | Lower coherence than same band — contradicted |
| Same band, $|\Delta\mathrm{lat}| < 10°$, $d > 8000$ km | −0.118 | 107,013 | — | Positive coherence required — contradicted |

### 4.5b Estimator Decomposition of the Fitted Scale

The pooled MGEX scale ($\lambda = 1862$ km) sits below the IGS-final
family reported in Papers 1–2 (3,330–4,549 km). Both analyses measure
the same statistic — the cosine of the magnitude-weighted circular mean
cross-spectral phase — so the offset can be decomposed
estimator-axis by estimator-axis on identical real station-days
(Step 3.5). Three documented settings differ between the pipelines: the
Welch segment length (96 samples here, $\Delta f \approx 34.7$ µHz,
versus a full-day segment in Papers 1–2, $\Delta f \approx 11.6$ µHz),
the per-epoch common-mode removal (applied here; not applied in Papers
1–2), and the clock products themselves.

On a 58-day sample spanning the archive (211,012 pair-days), the
factorial decomposition returns: canonical settings
$\lambda = 1285 \pm 233$ km (day-block bootstrap; the full-archive
estimate on all 486 days is 1,862 km); common-mode removal switched off,
$\lambda = 1368$ km; full-day spectral resolution with common-mode,
$\lambda = 1363$ km; and the Papers-1–2 estimator combination
(full-day resolution, no common-mode removal), $\lambda = 1598$ km.
Restricting the full-day spectrum to the resolved MGEX floor
($\geq 34.7$ µHz) returns $\lambda = 1045$–$1049$ km, confirming
that the two lowest resolved bins carry the longest-range content — as
required if the correlation is band-ordered. The estimator axes thus
move the fitted scale over $1.0$–$1.6 \times 10^3$ km on identical
data (≲±25% around the canonical value) and cannot supply the ~2×
offset to the IGS-final family; the residual resides in the clock
products themselves — their multi-GNSS combined construction, station
complement, and epoch coverage — rather than in the correlation
estimator.

The same ordering is recovered per sub-band: on the full-day matrix the
fitted scale declines from $\lambda \approx 1.8$–$2.0 \times 10^3$ km
below 30 µHz to $\approx 0.5$–$1.2 \times 10^3$ km above 50 µHz,
mirroring the IGS-final profile, which declines from
≈3.9–5.9×10<sup>3</sup> km in the tidal bands to ≈1.0–1.4×10<sup>3</sup> km
at intermediate frequencies (Paper 1 §3.5.1). Both products therefore
return the same band-ordered, distance-structured correlation; the MGEX
profile is a scaled-down version of the IGS-final one. The defensible
cross-product statement is a thousands-of-km coherence scale with a
shared band ordering — consistency with the geometric saturation anchor
at the ~1.9–4.5×10<sup>3</sup> km precise-product family level (raw-SPP
channels bound the scale from below at ~0.7–1.1×10<sup>3</sup> km and
the raw ionofree channel returns 4,767 ± 835 km; Paper 3) — rather
than a point-scale identification.

Two provenance details are bounded by the same diagnostic. First, GRG
product-days are natively 30-s sampled, so on the configured
5-min-labelled frequency grid those station-days resolve a
~0.35–5 mHz band rather than the nominal 10–500 µHz band; substituting
their correctly decimated contribution moves the pooled scale by ~1%
(1,271 → 1,285 km on the sampled day set; 9 of 486 days carry GRG data).
The per-center consequence is larger — the GRG-only matched-day scale
moves from 3,257 km to 1,933 km under the corrected pairwise epoch
intersection (Step 3.1) — because the per-center path additionally
truncated unequal-length series without timestamp alignment, a fault the
shared-grid construction removes. The pairing step now thins
finer-cadence products to the 300 s grid
before spectral estimation. Second, single-window realisations remain
epoch-dependent at the ~30–40% level, as reported in §4.1.

The fitted amplitude A decomposes along the same axes (Step 3.6;
46 COD days, 190,867 pair-days, all statistics evaluated on identical
pair sets). The phase-alignment index is bounded in [−1, +1] and
saturates once the in-phase component dominates the band-integrated
cross-spectrum, so A is not a linear estimator of the shared coherent
power: on the same bins the canonical configuration returns
A = 0.789 while the band-averaged real coherency — the linear
in-phase coherence measure — returns A = 0.250. The same
statistic-level inflation is visible inside the corpus on Paper 3's
raw-SPP analysis, which reports both statistics on identical pair-days:
normalised coherence amplitude 0.182 against phase-alignment 1.056,
a ~5.8× spread on 61.8 million pairs. Estimator conditioning
contributes the next factor: under the Papers-1–2 settings on the same
MGEX station-days (single full-day Welch segment, no common-mode
removal, no flat-line rejection) the fitted amplitude falls to
A = 0.466, and the flat-line rejection alone moves
A = 0.436 → 0.789 — held-constant clocks constitute 11.5% of
COD and 18.2% of GRG station-days and inject near-random pair phases.
Per-product ordering inside MGEX confirms the remaining conditioning:
on identical canonical settings GRG returns A = 0.462 against
COD 0.789. A labelled transfer calibration (Step 3.6; synthetic
red-noise pairs, estimator calibration only) inverts the published
amplitudes to underlying band coherences of ρ ≈ 0.30 (MGEX) and
ρ ≈ 0.05 (IGS finals) under each family's own spectral window —
a product-level difference in the shared in-phase budget, expected of a
multi-GNSS combined solution relative to GPS-only finals, rather than a
contradiction in the measured structure. The fitted offset shows the
same bookkeeping: all three product families return negative C₀, and the
deeper MGEX asymptote is the zero-sum signature of its per-epoch median
datum — toggling the common-mode subtraction on identical pair-days
moves C₀ from −0.015 to −0.066, the same distant-pair anticorrelation
produced by the conditioned datum kernel (Paper 1, §4.6). The
cross-product comparison is therefore carried by the signed,
distance-ordered in-phase decay and its scale ordering, not by amplitude
equality, and the ~7× headline ratio is a property of the statistic's
saturation plus product and estimator conditioning.

### 4.6 Null Tests

Null tests were run using the same MGEX combined-clock phase-alignment pipeline.
Three of the four nulls (temporal shuffle, spatial shuffle, and phase
randomisation) collapsed to R² ≤ 0, confirming that the recovered correlation
length depends on the true temporal, spatial and phase information. The
solar-rotation reassignment null, however, retains R² = 0.50 under the
primary parameters. A sensitivity analysis attributes this to a small-sample
bin artifact: short-baseline bins at 100–115 km (40 pairs) and 132–152 km
(23 pairs) receive extreme coherence values (0.29 and 0.21) under the
specific shuffle sequence, producing a spurious short-baseline spike that
the exponential fit interprets as structure. Under a robust analysis (minimum baseline 500 km, minimum bin
count 50) all four nulls give R² ≤ 0.122, confirming that the
distance-decay structure depends on the true temporal, spatial, and phase
information and is not a solar-period artifact.

| Null | R² (primary) | Passes? | R² (robust) | Passes? |
| --- | --- | --- | --- | --- |
| Temporal shuffle | −0.343 | Yes | −0.011 | Yes |
| Spatial shuffle | −0.000 | Yes | 0.028 | Yes |
| Phase randomisation | −0.001 | Yes | 0.122 | Yes |
| Solar-rotation reassignment | 0.499 | No | 0.010 | Yes |

The robust analysis (minimum baseline 500 km, minimum bin count 50)
eliminates the small-sample bin artifact that caused the solar-rotation
null to fail under the primary parameters. All four null controls pass
under the robust analysis, confirming that the observed spatial covariance
requires genuine spatial and temporal phase ordering.

## 5. Discussion

The combined MGEX receiver-clock product shows a distance-dependent structure
with $\lambda = 1862$ km. The estimate remains within 1,819–2,000 km
when the minimum baseline is raised from 50 to 2,000 km, while the binned-fit
R² improves from 0.70 to about 0.92. It persists on quiet days and after
distance matching, and it survives the two computed mundane channels rather
than only label-permutation nulls: the forward-modelled shared
satellite-error kernel is flat through the signal band and decays an order of
magnitude more slowly than the measured profile; conditioning the same
channel on the product's own ephemerides and satellite-clock content decays
no faster than λ = 6,650 km (~3.6× the measured scale; Paper 1, Step 5.2);
and the weather-organized
partitions return the opposite of the atmospheric signature — cross-hemisphere
pairs, which share no synoptic weather, carry the cleanest decay, and
same-latitude-band pairs beyond 8,000 km return negative coherence where a
zonal atmospheric kernel requires positive coherence. The ±112 km formal fit
uncertainty treats 1.83 million pair-days as independent; cluster-robust
resampling of whole epochs inflates the statistical error only modestly —
±86 km under day-block bootstrap, ±155 km under month-block bootstrap, and
±94 km under leave-one-month-out jackknife — so the within-product
statistical precision is about ±8%. The dominant scale uncertainties are
systematic rather than statistical: single-window realisations are
epoch-dependent at the ~30–40% level, and the matched-day cross-center
comparison contributes ~±600 km. A factorial decomposition on identical
station-days (Step 3.5, §4.5b) shows that the estimator axes — Welch
segment length and per-epoch common-mode removal — move the fitted scale
by ≲±25%, so the residual offset to the IGS-final family
(3,330–4,549 km, Papers 1–2) resides in the clock products rather than
in the correlation estimator, and the same band-ordered profile (longest
scale in the lowest resolved sub-band) is recovered in both product
families. The amplitude is conditioned the same way and more strongly:
the fitted A is a bounded, saturating phase-direction index (§4.5b, Step
3.6) whose ~7× offset from the IGS-final family decomposes into a ~3×
statistic nonlinearity, the documented estimator axes, and a
product-level difference in the shared in-phase budget — the linear
coherence measured on identical pair-days is A = 0.250, and the
deeper negative MGEX asymptote is the signature of its own zero-sum
datum, not a sign discrepancy. The result is nonetheless a useful
cross-product corroboration of a correlation scale in the thousands-of-
kilometres range, not a direct measurement-equivalence claim with the earlier
GPS PPP products.

The conventional directional measures do not replicate robustly. The full-range
EW/NS ratio is 1.29 but is non-significant under station-clustered bootstrap
(p = 0.49). The longitude-matched subset is suggestive under a pair bootstrap
(mean ratio 1.71, p = 0.058) yet non-significant under spatially clustered
resampling (p = 0.214). Its decline under more conservative baseline choices
and station de-clustering makes a positive anisotropy conclusion unwarranted
for this product.

The orbital finding is more informative but has a narrower evidential status.
The primary monthly $\lambda$ and EW/NS estimators are null after Bonferroni
correction, unsurprising with only 16 complete chronological months and
unstable monthly exponential fits. The model-free PA-difference statistic,
monthly mean EW minus NS phase alignment, has r = −0.716 with orbital speed
(p = 0.0018) and r = −0.767 with velocity projected onto the CMB dipole
(p = 0.0005). Because the nominal p-values assume independent monthly points,
they are reported alongside autocorrelation-aware diagnostics: effective-
degrees-of-freedom values p = 0.072 and p = 0.018, matched AR(1)-surrogate
values p = 0.019 and p = 0.003, and circular-shift values p = 0.125 and
p = 0.063 (the 1/16 floor for the projection channel). The association is
therefore tempered but not erased, and it is the physically motivated
velocity-projection channel that survives; leave-one-month-out correlations
remain within r = −0.66 to −0.86, excluding single-month dominance. These
signs and magnitudes are consistent with the preceding
TEP-GNSS analyses, which motivated the physical prediction before this MGEX
analysis. That history makes the result corroborative across products. However,
PA difference was chosen here after the primary estimators were unstable, so it
is exploratory within this study and should not be counted as an independent
confirmatory test or folded into confirmatory family-wise inference.

The corrected CMB analysis keeps January 2025 through April 2026 as separate
months. Its direct CMB-projection correlation is r = −0.733, consistent with
the stronger weighted orbital result (r = −0.767, p = 0.0005). The full-sky
scan has a nominal maximum at RA = 230°, Dec = +60° (r = −0.745), 82.7° from
the CMB dipole. That angular separation must not be interpreted as a failed
CMB test. Earth’s orbital velocity remains in the ecliptic plane, so dot-product
regression identifies only an axis’s in-plane projection; directions separated
substantially in three dimensions can give nearly identical annual predictors.
Here the CMB and nominal-best projected predictors correlate at r = 0.984, and
61 of 703 grid directions equal or exceed the direct-CMB |r|. The 30°
three-dimensional rule is thus not an identifiable criterion for this design.
The permutation result (plus-one p = 0.005 over 200 complete-grid shuffles)
supports directional annual structure, while the direct CMB projection remains
consistent with that structure. The 60° sensitivity result and daily analysis
are corroborative diagnostics, not independent confirmations. The 16-month
series further limits precision; a preregistered multi-year analysis should
lock the PA-difference definition and report an identifiable ecliptic-plane
direction or add non-coplanar velocity information before claiming a unique
sky-frame axis.

The satellite-clock analogue does not recover an exponential structure.
GLONASS, BeiDou, and Galileo have negative fitted R² values, while GPS has
R² = +0.13 with boundary convergence. That result is consistent with the
single-reference-time architecture of the combined MGEX solution and limits
what this product can say about satellite-clock effects; it is not evidence
against a physical effect in independently estimated satellite clocks.

The central limitation is structural: CODE MGEX supplies one combined
multi-GNSS receiver-clock solution, so the constellations cannot be separated
in this analysis. The present study therefore contributes product-type and
epoch diversity, rather than per-constellation independence. Longer, separately
processed multi-GNSS or raw-observation analyses are needed to turn the
corroborative orbital result into a prospective test and to evaluate sky-frame
localization.

## 6. Conclusions

This product-type independence stress test uses the CODE COD0MGXFIN combined
multi-GNSS receiver-clock product from 2025-01-01 to 2026-05-01. Its principal
result is a distance-dependent correlation scale of $\lambda = 1862$
km (R² = 0.699; formal within-fit σ = ±112 km, cluster-robust
σ = ±155 km from month-block epoch resampling), stable across the tested
baseline cutoffs and present on quiet
geomagnetic days, with epoch- and product-level scale systematics of ~30–40%
reported alongside the statistical error. Under the primary short-baseline parameters the
solar-rotation reassignment null retains R² = 0.50, traced to a small-sample
bin artifact at 100–152 km; the robust analysis (minimum baseline 500 km,
minimum bin 50) passes all four null controls and is reported as the clean
result, with the fitted scale essentially unchanged (1,819 km). This is cross-product corroboration of a spatial structure at
the thousands-of-kilometres scale reported in the earlier TEP-GNSS papers,
while remaining subject to the different metric and combined-clock product
construction. The fitted scale is a product-dependent phenomenological
observable — the transfer map from the microscopic screening operator
$\mathcal{S}_\Sigma(\mathcal{E})$ to $\lambda$ remains to be derived — so
the result is reported as corroboration of the TEP scale, not as a
parameter-level prediction test.

The secondary results are deliberately qualified. Standard anisotropy and
primary orbital metrics are not robust in this short MGEX interval. A
model-free PA-difference statistic has strong associations with orbital speed
(r = −0.716, p = 0.0018) and CMB-projected orbital velocity (r = −0.767,
p = 0.0005) across 16 chronological months. Because those signs were
anticipated by previous TEP-GNSS work, the result is corroborative across
products; because this operationalisation was selected after primary MGEX fits
were unstable, it remains exploratory within this study and needs a
prospective, longer-baseline replication.

The CMB-projected orbital association is likewise consistent: the direct
monthly CMB projection has r = −0.733, while the weighted PA-difference test
has r = −0.767 (p = 0.0005). The grid’s nominal RA = 230°, Dec = +60° maximum
is 82.7° from the CMB dipole, but this is not contrary evidence because orbital
velocity constrains only an axis’s ecliptic-plane projection. The grid cannot
uniquely localize a three-dimensional sky direction; its within-grid
permutation statistic (plus-one p = 0.005 from 200 shuffles) measures annual
directional structure. The satellite-clock analysis is limited by the combined
product’s shared reference time and does not recover an exponential signal; the
MGEX product therefore cannot test the satellite-clock channel at all, and a
per-constellation solution with independently estimated satellite clocks would
be required for that test. The
next decisive test is a preregistered chronological analysis over multiple
annual cycles, reporting an identifiable plane-projected direction or using
non-coplanar velocity information.

## 7. Data and Code Availability

All data are publicly available from NASA CDDIS
([cddis.nasa.gov](https://cddis.nasa.gov)) and the IGS
([igs.org](https://igs.org)) under their standard open-data
policies. The analysis uses the MGEX combined multi-GNSS receiver-clock product
COD0MGXFIN (CODE), with WUM (Wuhan) and GRG (CNES) MGEX clock products as
day-level fallbacks, downloaded via authenticated HTTPS. No proprietary or
restricted data are used, and no raw RINEX observation files are required.

Analysis code and pipeline outputs are available at
[github.com/matthewsmawfield/TEP-GNSS-MGEX](https://github.com/matthewsmawfield/TEP-GNSS-MGEX).

## 8. References

### 8.1 TEP Series

- Smawfield, M. L. (2025). *Temporal Equivalence Principle: Dynamic Time & Emergent Light Speed*. v0.14 (Jakarta). DOI: 10.5281/zenodo.16921911.

- Smawfield, M. L. (2025). *Global Time Echoes: Distance-Structured Correlations in GNSS Clocks* (Paper 1, TEP-GNSS). DOI: 10.5281/zenodo.17127229.

- Smawfield, M. L. (2025). *Global Time Echoes: 25-Year Temporal Evolution* (Paper 2, TEP-GNSS-II). DOI: 10.5281/zenodo.17517141.

- Smawfield, M. L. (2025). *Global Time Echoes: Raw RINEX / SPP Consistency Test* (Paper 3, TEP-GNSS-RINEX).

### 8.2 GNSS Data and Products

- Montenbruck, O., Steigenberger, P., & Hauschild, A. 2014, *GPS Solutions*, 18, 49, "Multi-signal GNSS data from the IGS MGEX experiment — status and outlook".

- International GNSS Service (IGS). MGEX Multi-GNSS Experiment. [https://igs.org/mgex/](https://igs.org/mgex/).

- NASA CDDIS. Crustal Dynamics Data Information System. [https://cddis.nasa.gov](https://cddis.nasa.gov).

- Center for Orbit Determination in Europe (CODE). MGEX combined products. [https://www.aiub.unibe.ch](https://www.aiub.unibe.ch).

### 8.3 Cosmology and CMB

- Planck Collaboration, Aghanim, N., Akrami, Y., et al. 2020, *A&A*, 641, A6, "Planck 2018 results. VI. Cosmological parameters", doi:10.1051/0004-6361/201833910.

### 8.4 Methodological References

- Welch, P. D. 1967, "The use of fast Fourier transform for the estimation of power spectra: a method based on time averaging over short, modified periodograms", *IEEE Trans. Audio Electroacoust.*, AU-15, 70–73.

Smawfield, M. L. 2026. Temporal Equivalence Principle series, Papers 0–25. Zenodo preprints and associated repositories.

## 9. Reproducibility Statement

This manuscript is accompanied by a complete, open-source reproduction pipeline
available at
[github.com/matthewsmawfield/TEP-GNSS-MGEX](https://github.com/matthewsmawfield/TEP-GNSS-MGEX).
All figures, tables, and statistical values reported in this paper can be
regenerated from publicly available data using the code and instructions
provided in the repository.

**Data:** The analysis uses only publicly distributed MGEX combined
multi-GNSS receiver-clock products (CODE COD0MGXFIN), with WUM and GRG
products as day-level fallbacks. Data were downloaded from NASA CDDIS
([cddis.nasa.gov](https://cddis.nasa.gov)) and the IGS
([igs.org](https://igs.org)) via standard authenticated HTTPS
endpoints. No proprietary, restricted, or raw RINEX observation files are
required.

**Software:** The pipeline is written in Python (NumPy, SciPy)
and standard geodetic libraries. All dependencies are listed in
`requirements.txt` in the repository root. No commercial or
closed-source software is required.

**Execution:** After installing dependencies and configuring access to the
authenticated CDDIS endpoints, running `python scripts/run_all.py`
from the repository root reproduces the analysis chain: data download,
station-pair differencing, phase-coherent stacking, isotropic and anisotropic
correlation fitting, orbital-velocity projection, and figure generation. Run
`python scripts/generate_figures.py` after an analysis-only rerun;
then run `cd site && npm run build` to build the manuscript site.
Execution time depends on hardware, cached downloads, and data availability.

**Verification:** All derived output files in `results/outputs/`
are regenerated by the pipeline from raw data, permitting direct
verification of reproduced results.

**Environment:** The pipeline has been tested on macOS 14 and
Ubuntu 22.04 with Python 3.11. All dependencies are listed in
`requirements.txt` for reproducible installation.

---

*This document was automatically generated from the TEP-GNSS-MGEX research site. For the interactive version with figures and enhanced formatting, visit: https://mlsmawfield.com/tep/gnss-mgex/*

*Source code and data available at: https://github.com/matthewsmawfield/TEP-GNSS-MGEX*
