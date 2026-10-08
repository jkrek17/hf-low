# Why is high-latitude Atlantic HF wind late in the life cycle? Pre-registration (RA-26)

Written and committed **before any outcome (LATE) was tabulated against any location predictor in this thread**.

What had been looked at: (1) the RA-7 write-up (PR 86): late share in the Atlantic at or north of 60N 52.7% against 26.9% south, OR 3.10, 2.85 without terrain-flagged fixes;
(2) **counts only**, which are not outcome tests: Atlantic events 1,039, of which 283 have their onset gust maximum at or north of 60N; of those 151 are LATE (h_on_minp >= 0), 66 LATE6, 30 LATE12;
145 of the 283 are terrain-flagged (maximum within 100 km of Greenland); (3) column names and missing-value counts of `hf_structure/results/fixes.csv` and `highlat/gustloc_fixes.csv`.
Not looked at: the share of the barrier pattern, coast distance, centre distance or wind direction by LATE status in the high-latitude group. (RA-7's T1, gust maximum more than 400 km from the centre, was
looked at for all events pooled, not for this group.)

All results are an ERA5 **proxy**, **pipeline A** (`research/era5/hf_history`, 800 km ocean gust index, threshold 71.7 kt), seasons 2004-05 to 2025-26 (decision 1), non-TC events (as RA-7). No pull.

## Question (agenda RA-26) and the answer I will give

Are late high-latitude Atlantic HF storms the ones whose gust maximum sits in northerly barrier flow near the coast and far behind the centre, or slow deepeners with the maximum near the centre?
Answer in plain words: yes, no (well powered), or can't tell, and by how much (late share against non-late share, with interval).

## Data

- Events and LATE: exactly RA-7's table and filters (`late_wind.load()`; events dropped for onset mismatch are counted).
- Group: **Atlantic events whose onset-fix gust maximum is at or north of 60N** (RA-7's T2, `gmax_lat >= 60`). Expected n = 283: late (h_on_minp >= 0) 151, non-late 132.
- Location of the maximum, from `hf_structure/results/fixes.csv` at the onset fix: `gmax_r` (km from centre), `gmax_coast_km` (km to nearest land, any coast), `terrain` flag (within 100 km of Greenland, as RA-7).
- Wind direction at the maximum: `highlat/gustloc_fixes.csv` (`max_wdir`, direction the wind blows FROM), joined on track and onset time; coverage reported; a missing value is excluded from that test, not imputed.

## Definitions

- LATE (primary) = `h_on_minp >= 0` (agenda wording "after the pressure minimum"; **63% of these are same-fix ties**, so LATE6 = `h_on_minp >= 6` is a registered sensitivity).
- BARRIER pattern = `gmax_coast_km <= 300` and `gmax_r > 400`.
- NORTHERLY = wind from 315 to 45 degrees at the maximum.
- Terrain-flagged fixes: **kept in the primary run** (RA-7's T2 included them), removed in a registered sensitivity run (S5). Removing them also removes most near-Greenland maxima, so S5 is expected to lower the BARRIER share in both groups; it is read as a difference between groups, not as a level.

## Tests (one family, BH-FDR over all 9; q reported within and across)

Unit: events. Intervals: season-block bootstrap (2,000 draws, seed 7), 22 seasons. Difference = late share minus non-late share.

- **P1 (primary)**: BARRIER share, late against non-late, in the high-latitude group.
- S1: coast <= 300 km alone. S2: centre distance > 400 km alone. S3: NORTHERLY share (among fixes with a wind direction).
- S4: P1 with LATE6. S5: P1 with terrain-flagged fixes removed. S6: P1 at the peak-gust HF fix instead of the onset fix.
- S7: near-centre alternative: maximum within 400 km of the centre (share), late against non-late.
- S8: the same P1 in the Atlantic south of 60N (specificity: barrier pattern should not distinguish late events there; reported as a contrast, not as support).
- S9: season-clustered logistic, LATE on BARRIER adjusted for marginal storm (peak gust < 76.7 kt), translation speed and month group; odds ratio with CR1 interval.

## Decision rules (thresholds fixed now)

- **YES** if late BARRIER share >= 0.60, non-late <= 0.35 (the agenda's prediction, both halves), and the lower 95% bound of the difference > 0.
- **Partly** if the difference is positive with its interval above 0 but a threshold is missed: report "directionally yes, smaller than predicted" with the shares.
- **No (well powered)** if the upper 95% bound of the difference is below 0.15 and power to detect the predicted difference of 0.25 is at least 80%.
- Otherwise **can't tell**.
- The slow-deepener alternative is supported only if S7 is higher in late events with its interval above 0 and P1 is not.

## Power

Planted-effect simulation from the observed group sizes (binomial draws at 0.60 against 0.35, and the season bootstrap): report power at difference 0.25 and the minimum detectable difference at 80%, for P1 and each secondary.
Count rule from the agenda: 151 late cases is above ~100, so the test is run.

## Looks

Uses all 22 seasons including 2015-25 once, descriptively; no model trained on one block and scored on another. By the agenda's count it is look #16 at 2015-25; no pre-2001 season is used (zero looks).

## What this cannot show

ERA5 under-resolves the gust peak; the proxy near Greenland has the highest false-alarm ratio (0.44); no terrain mask. The "coast" is the nearest land of any kind. A pattern in the proxy is not a demonstration of physical barrier flow without a wind-profile check.

## Deviations (post hoc)

(none yet)
