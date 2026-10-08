# Pre-registration: ENSO flavor and where Pacific storms deepen fastest (RA-10, ERA5 proxy)

Written and committed **before any deepening position was joined to any ENSO index**. Anything changed afterwards is logged under "Deviations (post hoc)" at the bottom and both versions are reported.

Question (agenda RA-10, following PR 38): PR 38 found no effect of the Modoki index (EMI) on where North Pacific cyclones *form* over the Kuroshio (RR 0.995 per SD, 0.917-1.079). Does EMI move where storms *deepen fastest* (the jet exit sets that), or where they reach lowest pressure, even though genesis does not move?

Plain-words answer I will give: **yes** (a shift in longitude or latitude that survives FDR and passes the consistency checks), **no, well-powered** (the whole 95% interval inside the shift the idea predicts), or **can't tell**, in each case with the size in degrees per SD.

## Data and pipeline

- **Tracks and fixes:** PR 38's committed tracker output, `research/era5/enso_kuroshio/data/{tracks,fixes}.csv.gz`. It is a **new MSLP-only tracker run that reuses pipeline A's detector and linker**; fixes are 6-hourly, Dec 1 to Mar 31, winters 1979-80..2025-26 (n = 47). The HF label (S11 only) comes from **pipeline A** (`pipelineA_match.csv`, gust index at 71.7 kt). Pipeline B is not used. Everything is an **ERA5 proxy**, not the archive. Tropical cyclones are in (as in pipeline A); the 25N floor on the deepening point removes most of them.
- **Indices:** `enso_kuroshio/data/{winter_indices.csv,sst_boxes_raw.csv}` (ERA5 SST, 1991-2020 monthly anomalies; EMI of Ashok et al. 2007; Nino3.4; Nino4 minus Nino3; DJF PNA from CPC).
- No pull. All inputs are committed.

## What had been looked at before writing this

PR 38's results (genesis counts, null). From the fix table only: the marginal distribution of the Bergeron-normalized 24 h deepening over all tracks globally (percent above 6, 12, 18, 24 hPa per winter), and that every track has contiguous 6 h fixes, to choose the population threshold. **No deepening position, longitude, latitude, or any join of deepening with an ENSO index, a winter, or a basin had been computed.**

## Definitions

- **Deepening window of a track.** For each fix t with a fix at t+24 h: ΔP24 = MSL(t) − MSL(t+24 h). Bergeron-normalized B = ΔP24 · sin 60° / sin(φ̄), φ̄ = mean latitude of the two fixes. A track's **max-deepening window** is the one with the highest B (first if tied). **Deepening point** = midpoint (lat, lon) of that 24 h window; its time is the window midpoint. A track needs at least 5 fixes (24 h) to have one.
- **Population.** Pacific box: deepening point in 25-67N, 120E-240E. Primary population **bombs: B ≥ 24 hPa/24 h** (the standard Sanders-Gyakum bomb; expected about 35 a winter). Secondary population **rapid: B ≥ 12**.
- **Window.** The deepening-window midpoint falls from 03 Dec 00Z to the last 18Z of February (DJF). Secondary JFM: 01 Jan to 28 Mar 18Z. (The first two days are dropped because tracks already present on 1 Dec start truncated; the pull ends 31 Mar 18Z.)
- **Position of lowest pressure** (S5): the fix with the minimum MSL along the track, for the same bomb population.
- **Outcome per winter:** the mean longitude (deg E) and mean latitude of deepening points of the population. Unit of inference is the **winter** (n = 47); storms in a winter are weights, not sample size.

## Predictor and lag

- **Primary: DJF-mean EMI**, with DJF Nino3.4 and a linear year term as controls (the PR 38 setting). ENSO is a seasonal boundary condition: the project default of an index averaged over days -10..-4 before a storm is meaningless for a three-month SST mean, and a few storms cannot change a three-month tropical SST mean, so the primary is concurrent DJF, as in PR 38.
- **Lagged contrast, registered (S3, S4): SON-mean EMI** (Sep-Nov preceding the winter, with SON Nino3.4), which strictly precedes every storm (lag one to five months). If DJF and SON disagree, I report both and rely on neither alone.
- Only indices from distinct regions are combined; PNA enters only in the mechanism check M1.

## Model and inference

- **P1 (longitude)**: weighted least squares, winter-mean longitude of bombs ~ EMI + Nino3.4 + year, weights = number of bombs that winter, predictors standardised over the 47 winters. Inference: Freedman-Lane permutation (10,000) of reduced-model residuals across winters, two-sided, plus a winter pairs-bootstrap 95% interval (5,000). The effect is **degrees per SD of EMI**.
- Same model for every test below unless stated.

## Tests (BH-FDR over all 18)

| id | outcome | change from P1 |
|---|---|---|
| P1 | longitude of bomb deepening point | none |
| P2 | latitude of bomb deepening point | none |
| S1 | longitude, rapid population (B ≥ 12) | population |
| S2 | latitude, rapid population | population |
| S3 | longitude, bombs | SON-lag EMI and SON Nino3.4 |
| S4 | latitude, bombs | SON-lag EMI and SON Nino3.4 |
| S5 | longitude of the position of lowest pressure, bombs | outcome |
| S6 | longitude, bombs | N4-N3 replaces EMI |
| S7 | longitude, bombs | no year term |
| S8 | longitude, bombs | no Nino3.4 term |
| S9 | B-weighted centroid longitude of rapid storms (weights = B) | outcome |
| S10 | longitude, bombs, JFM window | window |
| S11 | longitude of the deepening point of HF-reaching storms (pipeline A, 71.7 kt, Pacific), **2004-05 on, 22 winters**, no year term | population, era (gust label, Decision 1) |
| S12 | El Nino winters only (Nino3.4 ≥ +0.5): CP minus EP mean bomb longitude, permutation of class labels (Kug rule from the same indices) | class |
| S13 | winter 75th percentile of bomb longitude (eastern extent) | outcome |
| S14 | winter 25th percentile of bomb longitude (western extent) | outcome |
| S15 | longitude, bombs, 1979-80..2003-04 only (25 winters, pressure-only, no gust label) | era |
| S16 | longitude, bombs, 2004-05..2025-26 only (22 winters) | era |

Mechanism, outside the FDR family: M1 adds DJF PNA to P1 (report EMI effect with and without it, the EMI-PNA correlation, PNA's own coefficient). No jet-level field is pulled.

## Decision rules (fixed now)

- **Yes** (longitude): P1 permutation p < 0.05 **and** q < 0.10 over the 18 tests, the same sign in S7 and S8, and the same sign in S15 and S16. Same rule for latitude with P2.
- **No, well-powered:** the 95% interval of P1 lies inside ±3.0 degrees longitude per SD (the shift the idea predicts: at least 3 degrees per SD) **and** the minimum detectable shift at 80% power is ≤ 3.0. For latitude the stated size is ±2.0 degrees per SD (registered now; the agenda gives none).
- **Can't tell:** anything else, and always for S11, S12 if their intervals are wider than the stated size.
- A secondary result "holds up" only at q < 0.10; every test is reported with p, q, interval, and the minimum detectable shift.
- **Power:** minimum detectable shift by planted-effect simulation: add δ · EMI_z(winter) degrees to the longitude of every bomb, refit with the permutation test, 300 plants per δ, δ from 0.5 to 6 degrees per SD. Reported for P1, P2 and S3.

## What this cannot show

The depth of a 24 h window at 0.25 degrees from a tracker that smooths on about 0.5 degrees; a proxy for observed deepening. No land mask on the deepening point, so a bomb whose window is over a coastal low counts. Atlantic not tested. Not a test of jet position itself (no wind fields). ENSO flavor has few independent events (13 El Nino winters, 4 CP), so S12 is weak; the continuous index is the better-powered test. The 18 tests share outcomes, so BH is a mild correction.

## Looks

No model is fitted on one block and scored on another; all 47 winters, including 2015-25, are used once in the regression, so this is logged as one further look at the 2015-2025 block (the agenda's fifteenth) and, because P1 and its pressure-only variants use 1979-2000, as a further look at pre-2001 seasons. Counts are added to `research/era5/hemispheric/results/heldout_looks.log` at close-out.

## Deviations (post hoc)

1. Longitude midpoint was first an arithmetic mean, wrong across the 0/360 meridian (40 Atlantic-Europe tracks landed in the Pacific box); corrected to a circular midpoint. 2. S5's first run required the minimum-pressure fix to lie in the box; the registered definition does not, and is restored. Both versions are in the README; no decision-rule outcome changed. A 200-permutation smoke run preceded the reported run.
