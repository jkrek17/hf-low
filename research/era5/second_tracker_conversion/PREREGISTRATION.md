# Pre-registration: does the pattern's effect on conversion survive a second tracker's denominator? (RA-30)

Written and committed **before any tracker track was joined to a deepening flag, an HF label or a pattern index**, and before the
full M and V track sets were regenerated or reopened. ERA5 **proxy** throughout. Two pipelines are involved and are named with every number:
**pipeline A** (`research/era5/hf_history`, 0.25 degrees, 800 km ocean gust, HF at 71.7 kt) supplies the HF label and the PR 85
reference numbers; **trackers M and V** (`research/era5/second_tracker`, PR 116, 1.5 degrees, 6-hourly WeatherBench2) supply the cyclone
denominator. Pipeline B is not used.

## Question

PR 85 found that the large-scale pattern (PR 41 hemispheric index) acts on the **conversion** of deepening cyclones to HF lows,
P(HF | deepening), about x1.25 per SD in the Atlantic and x1.14 in the Pacific, and not on whether cyclones deepen. Both the
denominator (the deepening cyclones) and the numerator came from pipeline A's own tracks. Does the conversion effect keep its sign
and size when the deepening cyclones are found by a different tracker? Answer in plain words per basin: holds, does not hold, or can't tell, and by how much.

## Already looked at

- PR 116 inventory only: event-level recall of A's HF events by M and V (88.9-97.6%), track counts (M 83,038, V 50,040), the
  columns of the committed `events_{M,V}.csv` (track, matched, tr_track, n_hit). Nothing about deepening, HF share per tracker track, or any index.
- The PR 85 reference results (A2 Atlantic x1.252, 1.189-1.331; Pacific x1.137, 1.074-1.223, 22 seasons, pipeline A).
- **Not** looked at: any tracker track's MSLP change, any tracker-based weekly count of deepening or HF-and-deepening cyclones.

## Data and units

- **Seasons.** 2004-05 to 2021-22 (label 2004..2021), October-April, 18 seasons. RA-15 coverage ends January 2023 (WeatherBench2),
  so seasons 2022-23 to 2025-26 are not covered. Decision 1: fit and test 2004-05 onward only; no pre-2004 season is used.
- **Cyclones (denominator).** Every tracker track of 24 h or more (the PR 116 minimum), placed in the week of its first fix; weeks as
  PR 85 (30 weeks of 7 days from 1 October, effective n 18 seasons). Basin by pipeline A's rule at the track's intensity-peak fix
  (as `match.py` S1/S2). Transitioning tropical cyclones are in. Basins are never pooled.
- **HF label.** A tracker track is HF if it is the matched track (`tr_track`, rule of PR 116: at least half of the event's HF fixes
  within 500 km (M) or 700 km (V) at the same time) of a pipeline A HF event; otherwise not HF. The label is A's gust label carried over
  through the PR 116 match; the trackers carry no gust. The ~3-11% of A events that no tracker track matches are lost from the numerator
  (recall is the PR 116 figure per tracker and basin) and are reported, as is the chance-match risk (V's looser rule).
- **Deepening cyclone (D).** As PR 85: the track has at least one 00/12 UTC fix whose predecessor 12 h earlier is a fix of the same
  track and whose 12 h central-pressure change `dp12 <= -3.6` hPa. Whole-track flag (PR 85's stated limit is kept).
  - **M (primary for this rule):** central pressure is the track's own MSLP minimum value (1.5 degree).
  - **V:** the V track has no pressure, so each V fix gets the minimum of the 1.5 degree MSLP within 500 km of the fix (WeatherBench2
    `mean_sea_level_pressure`, October-April 2004-2021 only; a small additional pull, stated before it is made: under 2 GB, a subset of the 4.5 GB M pull).
  - 1.5 degree pressure is shallower than A's 0.25 degree pressure, so the same hPa rule selects a different set. The registered
    primary keeps -3.6 hPa per 12 h (comparability with PR 85); the shift is reported, and -2.4 hPa per 12 h is a registered secondary
    (S2) so the reader can see the sensitivity to the coarser grid. The primary is not re-chosen after results.
- **Pattern index.** `idx` of `hemispheric/results/oos_index_{atl,pac}.csv` (PR 41 leave-one-season-out), standardised over the
  weeks of the 18 seasons in each basin (so per +1 SD of this sample). Not refitted. Lag as built into that index (days -10 to -4 before the week). No new pattern is searched.

## Primary tests

Weekly Poisson models exactly as PR 85 Part 1 (month effects Oct..Apr, log1p of the same outcome in the previous 7 days, `idx`; season-clustered
sandwich; `share_mechanism`'s `chanlib` reused, `joint_se` for contrasts). Outcomes per basin and tracker: N (all tracks), D, HD (HF and deepening).
With b_X the `idx` coefficient: **A2 = b_HD - b_D**, the log rate ratio of P(HF | deepening) per SD.

- **P1-P4:** A2 for {V, M} x {Atlantic, Pacific}. V is primary (agenda), M a registered secondary reading of the same question, and both are reported.
  Family of the four tests: p by season-block permutation of `idx` (2,000 permutations, floor 1/2001), Benjamini-Hochberg q within the 4.
- **Reference.** Pipeline A on the same 18 seasons and the same weekly design (re-run, not quoted from PR 85's 22 seasons), so the comparison is like for like.
  The paired difference A2(tracker) - A2(pipeline A) has a 95% season-block bootstrap interval (2,000 draws, seed 20261011).

## Decision rule (per basin and tracker, fixed now)

1. **Holds** if A2 > 0 and q < 0.05 and the point estimate lies inside PR 85's 95% interval (Atlantic 1.189-1.331, Pacific 1.074-1.223 as rate ratios)
   (the agenda prediction). If the point lies outside but the paired difference to pipeline A on the same seasons includes 0, report "holds in sign, size differs from PR 85's but not from pipeline A on these 18 seasons".
2. **Does not hold** if the 95% interval of the rate ratio lies wholly at or below 1.0, or the sign reverses with q < 0.05.
3. **Can't tell** otherwise. Power: the minimum detectable effect is 2.8 x the clustered SE of A2 (rate-ratio scale), reported for each test; a "does not hold" or null needs the interval to exclude the PR 85 point, else "can't tell".
4. A well-powered null is a full result.

## Secondary tests (own BH families, none changes a primary)

- **S1 Atlantic north of 60N.** A2 for Atlantic tracks with peak latitude at or north of 60N and, separately, south of it (V and M). M recall there was 74% near Greenland and 85% north of 60N, so a M-based label loses more events there; reported separately and not pooled in the primary.
- **S2 Deepening threshold -2.4 hPa per 12 h** (1.5 degree grid sensitivity), both trackers.
- **S3 P(deepening)** (A1 = b_D - b_N), the PR 85 well-powered-null check (interval inside 0.95-1.05), both trackers.
- **S4 Match rule x2 and x0.5 distance** (PR 116 S4), both trackers, to bound HF-label noise.
- **S5 Tracker-agnostic HF label**: HF if matched by either tracker (a track is HF if matched by M or V) for the denominator of the other (a cross-check that the label, not the denominator, drives any difference).

## Power and limits

- Report the clustered SE and MDE of A2 per test and basin, and a planted-effect check: scale A2 in the HF label of V/M tracks to the PR 85 value and to 1.0 and report how often the decision rule returns holds / does not hold.
- The trackers carry no gust: recall below 100% and chance matches both move the HF label; this noise biases A2 toward 1 and is reported, not corrected.
- Only 18 of 22 gust-era seasons are covered; the 2022-23 to 2025-26 seasons cannot be added without a new pull.
- Not a new pattern search and not an independent sample: it reuses the 22 gust-era seasons' pattern index and A's HF label. It is a sensitivity re-run, counted as a look at 2015..2025 (conservatively, number assigned by the agenda thread), in this study's own log `research/era5/looks/second_tracker_conversion.log`. Zero looks at pre-2004 seasons.

## Deviations (post hoc)

Logged after the plan was committed. None yet.
