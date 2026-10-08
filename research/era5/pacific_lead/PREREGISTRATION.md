# Does a Pacific hemispheric state lead the Atlantic HF state by about a week? Pre-registration (RA-5)

Written and committed **before any Pacific predictor was joined to any Atlantic outcome or index in this thread**.
What had been looked at: the file structure of `hemispheric/results/{oos_index_*.csv,weekly_table.csv.gz,frozen_primary.json}`
(column names, row counts), and nothing else. No cross-basin correlation, lagged or not, has been computed by this thread.
Background that predates this plan: PR 41 (each basin's pattern index against its own basin's counts) and PR 64 (channels).

ERA5 is a **proxy** for the atmosphere. Pipeline A (`research/era5/hf_history`, 800 km ocean gust index, 71.7 kt) appears only
in secondary checks S3 and S4 and is named when it does. The archive is the primary outcome.

## Question (agenda RA-5)

After removing month and each basin's own pattern, does the Pacific state of one or two weeks earlier predict the Atlantic HF
state, as a Rossby wave packet crossing the hemisphere (about a week) would imply?

Answer given in plain words: **yes** (a positive lead, replicated), **no** (well powered), or **can't tell**, and by how much
(rate ratio per SD of the Pacific index, and partial correlation, with intervals).

Prediction (written before data): a positive lag-1 association, so the Pacific index leads the Atlantic one by one week. A
flat association means the Pacific ridge in the Atlantic pattern (PR 41) is coincident with, not a precursor of, the Atlantic
state. I have no evidence for either; a null is plausible.

## Data and timing

- **Indices.** The out-of-sample weekly pattern index of each basin from PR 41 (`oos_index_atl.csv`, `oos_index_pac.csv`):
  each season's index comes from a fit on the other 21 seasons (leave-one-season-out). 22 seasons (2004-05 to 2025-26, label =
  starting year), 30 weeks of 7 days from 1 October, 660 weeks per basin. The index of week `w` is the 7-day mean of the daily
  anomaly fields over days S-7..S-1 (S = start of week w). Both indices are standardised over the 660 weeks as in
  `hem_channels/chanlib.index`.
- **Lag.** The Pacific index of week `w-k` covers days S-7(k+1)..S-7k-1 before week w: k=1 is days -14..-8, k=2 is days
  -21..-15. This extends the project's default lag range (-10..-4) on purpose, because the question is about a one- to
  two-week lead of one basin over the other; the Atlantic week's own index stays at the default days -7..-1. The lag range is
  not extended further.
- **Outcome, primary: the archive.** Weekly counts of HF lows per basin (the `y` column of the index files). Pacific archive
  count of week `w-k` is the count-predictor. The Atlantic count of the 7 days before the week (`prev_atl` of
  `weekly_table.csv.gz`) is the autoregressive term, entered as log(1+count).
- **Common sample.** Weeks 2..29 of each season (28 weeks, 616 rows), so that lag 2 never reaches into the previous season. The
  same rows are used for every primary test.
- **Season split.** All 22 seasons; the Pacific and Atlantic indices are already out of sample by leave-one-season-out. No hold-out
  is created here.
- **Replication (secondary S3): 1979-80 to 2000-01**, within-era only (decision 1). Frozen primary-split patterns of PR 41
  (fitted 2004-14) applied to the weekly PC scores in `weekly_table.csv.gz`, pipeline A depth counts as outcome (cuts 966.2 hPa
  Atlantic, 965.0 hPa Pacific, as fixed in `freq_split`, not re-chosen). The pre-2001 archive is short-counted and is not used.
  Pressure/field quantities may cross eras; no gust level or trend is claimed.

## Models

Let `A(w)`, `P(w)` be the Atlantic and Pacific pattern indices, `yA(w)` the Atlantic archive count, `yP(w)` the Pacific archive
count, `M(w)` month-of-week dummies (as in `chanlib.month_dummies`), `prevA(w)` the Atlantic count in the 7 days before the week.

| test | outcome | model | tested term (one-sided, positive) |
|---|---|---|---|
| P1 | `yA(w)` | Poisson: M + log(1+prevA) + A(w) + P(w-1) | slope of P(w-1), reported as rate ratio per SD |
| P2 | `yA(w)` | Poisson: M + log(1+prevA) + A(w) + P(w-2) | slope of P(w-2) |
| P3 | `A(w)` | OLS: M + A(w-1) + P(w-1) | slope of P(w-1), reported with partial r |
| P4 | `A(w)` | OLS: M + A(w-1) + P(w-2) | slope of P(w-2) |
| P5 | `yA(w)` | Poisson: M + log(1+prevA) + A(w) + z[log(1+yP(w-1))] | slope of the Pacific count, per SD |
| P6 | `yA(w)` | Poisson: M + log(1+prevA) + A(w) + z[log(1+yP(w-2))] | slope of the Pacific count, per SD |

P1, P2, P5, P6 remove "the pattern in each basin" by holding the Atlantic pattern of the same week and the Atlantic's own
recent count fixed; P3 and P4 remove month and the Atlantic index's own persistence. One-sided because the prediction is
positive. A significantly negative slope is reported and counts as "no lead of the predicted sign".

**Inference.** Season-block permutation: whole 30-week Pacific seasons (index and archive count together) are shuffled among
seasons before the lag is formed within season, 10,000 shuffles, one-sided p = (1 + #{null >= observed}) / (1 + 10,000). Because
whole seasons keep their week-of-season position, the seasonal cycle shared by the basins stays in the null and only
interannual and subseasonal association is tested. Interval: 95% season-block bootstrap (2,000 draws). Cluster unit =
season, so effective n is the number of seasons (22), not weeks.

**Multiplicity.** Family 1 = P1..P6 (BH q over six). Family 2 = every pre-registered test below (P1..P6, S1..S4 rows; BH q over
all). Both are reported for every test, and passing is counted as "k of N".

## Secondary checks (labelled; all included in Family 2 where a p is computed)

- **S1 shared drivers.** P1, P2, P3, P4 refitted with the named indices NAO, PNA, ONI, MJO1 and MJO2 (as in `weekly_table.csv.gz`,
  same window convention as the pattern) added as controls at both the Atlantic week and the Pacific lag window. 4 tests.
- **S2 direction.** The reverse of P1 and P3: Atlantic index/count of week w-1 against the Pacific index/count of week w, with
  the symmetric controls. 2 tests. If the Pacific is a precursor of the Atlantic, P1/P3 exceed these; the difference of the two
  slopes (P to A minus A to P) is given with a bootstrap interval, not tested.
- **S3 within-era replication, 1979-2000.** P1, P2, P3, P4 on pipeline A depth counts with the frozen patterns. 4 tests. This is
  the third look at pre-2001 seasons (PR 41 S5 and PR 64 S7 were the first two).
- **S4 proxy outcome 2004-2025.** P1 and P2 with the pipeline A HF-equivalent weekly count (gust >= 71.7 kt, tracks counted by start date,
  as in `weekly_table.csv.gz`'s `pAhf_atl`) in place of the archive. 2 tests.
- **S5 lag profile (descriptive, not tested).** Partial correlation of Pacific index at lag k = -3..+3 (negative = Atlantic
  leads) with the Atlantic index, controls as in P3, with the 95% envelope of the permutation null. No test of any single lag.
- **S6 same-week contrast (k = 0), labelled circularity contrast**, not used for inference.
- **S7 stability.** P1 and P3 in 2004-14 and 2015-25 separately (descriptive).
- **Not tested here: the daily structure of the lag.** The literature's wave-breaking lag of about six days cannot be separated
  from the weekly design: the lag-1 Pacific window (days -14..-8) and the Atlantic window (days -7..-1) are centred seven days
  apart, which brackets it. A daily-resolution version needs the daily fields (about 17 GB, the `hemispheric/fields.py` pull);
  it is not part of this answer and is offered separately.

## Decision rules (thresholds fixed now)

1. **Yes, replicated.** A positive-signed P1 or P3 has Family 1 q < 0.05, and the matching S3 analogue (P1 or P3) has the same
   sign with one-sided p < 0.05. State the lag (1 or 2) and the size.
2. **Yes, not replicated.** P1..P6 contains a q < 0.05 positive slope, but S3 does not meet the rule above. Reported as a lead
   in 2004-2025 only.
3. **No (well powered).** No primary test has q < 0.05 in the positive direction, the planted-effect power of P1 at rate
   ratio 1.10 per SD is at least 80% and of P3 at partial r = 0.10 is at least 80%, and the 95% upper bounds of the P1 and P3
   estimates are below 1.10 and 0.10.
4. **Can't tell.** Everything else. The minimum detectable effects are stated with the answer.

Power is estimated by planted-effect simulation: P1 by adding the planted log-RR times the standardised P(w-1) to the fitted
null linear predictor and drawing Poisson counts (the training-period dispersion in PR 41 was 1.0; the dispersion of the
fitted null is checked and reported); P3 by adding the planted slope to the null fit and resampling null residuals in season
blocks. 400 simulations per grid point (RR 1.03, 1.05, 1.08, 1.10, 1.15, 1.20; r 0.03, 0.05, 0.10, 0.15), 1,000 permutations
inside each simulation.

## Looks spent

This plan spends: one look at the 2015-16 to 2025-26 seasons (the seventh; they are inside the 22 seasons of every primary
test, with the index computed leave-one-season-out), one look at the 2004-05 to 2014-15 seasons, and the third look at
1979-80 to 2000-01. They are appended to `hemispheric/results/heldout_looks.log` when the results are written.

## What this study cannot show

- A cause or a physical path. A positive lead is compatible with the Pacific state shifting the Atlantic pattern, and with a
  shared slow driver; S1 removes the named indices only.
- Anything about the daily wave-breaking structure (see above).
- That the leading index is the right Pacific variable: it is the one fitted to Pacific counts, so a lead of a different
  Pacific feature would not show here.
- 22 seasons give about 22 independent seasonal anomalies, but the week-to-week autocorrelation of an index is high; the
  permutation unit is the season for that reason.

## Deviations (post hoc)

Logged when the results were written. None changes a test, a threshold or a sample chosen in advance.

1. **Two MJO weeks are missing** in `weekly_table.csv.gz` (2021-22 week 13 and 2022-23 week 13, both MJO components). They are set
   to 0 (the mean) after standardising. Affects S1 only (2 of 660 weeks).
2. **Power level.** The plan said "power of at least 80%" without a significance level. Both are reported: one-sided p < 0.05
   (the level a single pre-specified test would use) and p < 0.05/6 (Bonferroni over the six primary tests). Rule 3 is met at the
   first and not at the second; see the README.
3. **Power computation.** The 400 simulations and 1,000 shuffles are as planned. For the OLS tests the shuffled slopes are
   computed by residualising on the controls (Frisch-Waugh), algebraically the same slope as a refit; Poisson tests refit.
4. **S5 sample.** The lag profile uses weeks 3..26 of every season so that all lags -3..+3 share one sample.
5. **Not in the plan, run after the primary results were in (post hoc):** the month-only lagged correlation of the two indices
   (`posthoc.py`), shown beside S5, and the half-by-half S7 numbers were computed with the primary ones. Both are labelled.
