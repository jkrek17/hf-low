# Do archive HF fixes per event drift, and is it a step at a scatterometer change? (RA-11, pre-registration)

Written 2026-10-08 before any change-point, step, covariate or proxy-contrast outcome was computed. The archive is the OPC HF
archive (`data/hf_lows/HF_Data_-_{Atl,Pac}.csv`). ERA5 appears only as a **proxy** contrast (pipeline A,
`research/era5/hf_history`). Changes made after the first outcome go under "Deviations (post hoc)" with both versions shown.

## Question (RA-11)

The archive's HF fixes per event fall about 0.18 per decade over 2004-05 to 2025-26 (t about -2.5; `drift_check-result.txt`) while
events and total HF fixes per season are flat. Is that a gradual drift or a step at a known satellite-wind change? Does it matter for
Decision 1 (gust drift) and the sustained-HF results (PR 79)?

Plain-words answer to give: "step at a sensor date", "gradual drift", "decline is real but its form can't be told", or "no
detectable decline", each with the size (fixes per event) and the minimum detectable effect.

## What was looked at before writing this

- Reported only: fixes per event -0.181 per decade (t -2.54) on the built payload, -0.180 (t -2.51) from the CSVs, season means for
  2004-2010 (3.05, 3.39, 3.35, 3.30, 3.26, 3.21, 2.95) in `drift_check-result.txt`; the single-fix share (~26%) and counts at >=2 and >=3
  fixes in `count_reconcile` and PR 79. No step, change-point, covariate, basin split, latitude split or proxy contrast has been computed.
- Data definitions only (columns, ID format, fix counts). Event = (basin, ID); season = year in the ID's first four digits (June-May);
  HF fix = row with Category `HF` (the rule `drift_check.py` uses; `DHF` is not counted).

## Sensor dates (verified before use; the agenda recalled them from memory)

| change | date used | status |
|---|---|---|
| QuikSCAT stops (antenna rotation lost) | 2009-11-23 | verified: CIMSS satellite blog, "QuikSCAT satellite ceases operations" |
| Metop-A (first ASCAT) launch | 2006-10-01 (month level) | launch month verified (NOAA OSPO ASCAT overview); the operational-use date is NOT verified |
| Metop-B (second ASCAT) launch | 2012-09-01 (month level) | launch month verified (NOAA OSPO); operational-use date NOT verified |
| Metop-C (third ASCAT) launch | 2018-11-01 (month level) | from memory, NOT verified |

Month-level dates and the unverified operational lags mean the ASCAT steps are placed on launch month, not on the date forecasters
saw the data; a step that lags by months to a year is not excluded and is reported as a limit. OSCAT (2009-2014, India) is not in
the plan; no source checked says whether OPC used it.

## Outcome and sample

Archive events (all HF, both basins pooled, Atlantic and Pacific separately as secondary), seasons 2004-05..2025-26 (22 seasons).
Outcome n_hf = HF fixes in the event. The step is applied by the date of the event's first HF fix. Seasons before 2004-05 are not used
(Decision 1 and the archive's incomplete 2001-04 seasons). No gust-based pre-2001 data enter. No held-out model is scored; there are no
looks at 2015-2025 or pre-2001 to log (a descriptive record, nothing fitted on one block and scored on another).

## Duplicate archive IDs

`count_reconcile/results/linked_ids_for_review.csv` lists 21 ID pairs that look like one low split in two (11) or recorded twice
(10). Primary: raw IDs, as `drift_check.py` does, so the -0.18 reproduces. Sensitivity (not in the FDR family): the 21 pairs merged
into one event (earlier ID, union of fixes with identical timestamps counted once); P1 and P2-P5 repeated, and the conclusion is
flagged if it changes (q crossing 0.05 or the preferred form flipping). Merging moves seasonal event counts by 0.6 (Atlantic) and
0.3 (Pacific) per season (count_reconcile), so a change is not expected.

## Primary tests (family P, 5 tests, Benjamini-Hochberg q)

- **P1 linear drift.** OLS of the 22 season means of n_hf on season; slope per decade, t on n-2 df, plus an event-level bootstrap interval resampling seasons (2000 draws; seed 11). Durbin-Watson reported.
- **P2-P5 steps** at 2006-10-01, 2009-11-23, 2012-09-01, 2018-11-01. Statistic: mean n_hf after minus before, event level. p from the season-block bootstrap (2000 draws; seed 11, resampling seasons whole; two-sided percentile p). Number of seasons each side is reported. P2 has two seasons before it, so it is a weak test and said to be.

## Secondary (family S, BH within S and across P+S)

- **S1 unknown-date scan.** Split points at season boundaries 2007-08..2021-22 (at least 3 seasons each side); maximum |t| of the
  step in season means; p by permuting season means (exchangeability under "no structure"; approximate because it ignores a trend).
- **S2 form comparison.** AIC on season means (Gaussian): null, linear, each P step, the S1 best step. Rule: a step is preferred to
  linear if its AIC is lower by at least 2; linear preferred to a step if lower by at least 2; otherwise "form can't be told".
- **S3 basins.** P1 and the QuikSCAT step (P3) separately for Atlantic and Pacific.
- **S4 shape.** Season-share of one-fix events and of events with at least 3 fixes: trend per decade (links to PR 79).
- **S5 Atlantic north of 60N.** Atlantic events whose first HF fix is north of 60N versus south: the QuikSCAT step in each and the
  difference (the agenda's QuikSCAT-label-bias item). Exploratory; Atlantic north of 60N is the lower-confidence group.
- **S6 covariates.** Event-level OLS of n_hf on season (per decade) with event minimum archive pressure, basin and calendar month of
  the first HF fix: does the slope change? Also the trend of event minimum pressure itself. If the adjusted slope drops, part of the decline is storm character, not practice.
- **S7 proxy contrast.** The same P1 on ERA5 pipeline A HF-equivalent tracks (`lifecycle_events.csv`, `n_hf`, seasons 2004-2025,
  basins pooled) and the archive-minus-proxy slope difference (season-block bootstrap). A flat proxy next to a falling archive is
  evidence for recording practice, with the caveat that the proxy has its own cut (71.7 kt gust) and never sees the archive's labels.

## Decision rules (thresholds)

- **Gradual drift**: P1 q < 0.05, no P step q < 0.05, and S2 does not prefer any step.
- **Step at a sensor date**: that P step q < 0.05 and S2 prefers it to linear.
- **Decline real, form can't be told**: P1 q < 0.05 and the above two not met (steps and linear fit within 2 AIC).
- **No detectable decline**: P1 q >= 0.05, and a "no" only if power >= 80% at the stated effect; otherwise "can't tell".
- A step with q >= 0.05 is not "no step" unless its MDE is stated and below the effect predicted by the decline (-0.4 fixes over the record).

## Power plan

Planted-effect simulation: take the real events, remove the fitted linear trend, resample whole seasons with replacement, add (a) a
linear decline of 0.181 per decade and (b) a step of 0.40 (the total linear change over 22 seasons) at each sensor date, 1000 runs
each. Report detection rates for P1 and each step, the MDE (smallest step detected 80% of the time at alpha 0.05) per date, and a
confusion table for S2 (does AIC pick linear when truth is linear, step when truth is step).

## What this cannot show

Why the practice changed (forecaster behaviour, warning criteria, the spreadsheet's own conventions) and when forecasters actually
gained a sensor. A decline and a null are both consistent with a slow change in how forecasters pick the last HF fix. The archive
is the reference for ERA5 evaluation, so a falling fixes-per-event matters for fix-hour and duration measures, not for event counts.

## Deviations (post hoc)

None yet.
