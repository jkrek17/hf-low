# Does the NAO x PNA shift of Atlantic storm longitude replicate on 1979-2000? (pre-registration, RA-8)

Written 2026-10-08, before any 1979-2000 position was related to any index. ERA5 **proxy**. Not observations.

## The lead being retested

`research/era5/tele_combos` (hf-low PR 40, `results/results.csv`, tier S4, test T3): for **all cyclones** in the **Atlantic**, the interaction NAO x PNA
(both averaged over days -10..-4 before the track's first fix, standardised) moves the longitude of the **gust-index peak** by **+1.11 deg per SD x SD** east
(95% season-pairs interval +0.51 to +1.86, joint lon+lat p 0.003, BH q 0.18 across 68 tests, 25 seasons 2001-02..2025-26, Oct-Apr, 10,988 tracks). It was found by
searching 68 tests, so it is a post hoc lead and its size is probably inflated (winner's curse). For HF-equivalent lows the same coefficient was +0.17 (no effect).

## Why the replication cannot be identical (stated up front)

1. **Position.** Gust-based quantities before 2004 are barred (Decision 1), and pipeline A's `all_tracks.csv.gz` stores only the gust-peak position. The replication
   uses the **position of the minimum MSLP** fix of each track. A failure is therefore partly ambiguous: the lead may not exist, or it may exist for the gust peak and not for the pressure minimum.
2. **Tracks.** `all_tracks.csv.gz` has no fixes. The committed fix-level table that covers 1979-2000 for all lows is the MSLP-only tracker run of PR 38
   (`research/era5/enso_kuroshio/data/{tracks,fixes}.csv.gz`): pipeline A's detector and linker (`hf_history/extract.py::lows_at`, `track.py::link`; PR 38 recovered 1.0000 of
   pipeline A's Pacific peak fixes) on ERA5 MSLP only. No pull is made. The same linker means the same population rule (low below 1010 hPa, track of at least `MIN_LEN` fixes, 20-75N).
3. **Season.** That table holds **1 Dec - 31 Mar** only. The lead used 1 Oct + 210 d (Oct-Apr). Genesis days 3 Dec - 31 Mar are used (first two days dropped, as PR 38 did, because a
   low already present on 1 Dec is not new). Winter s is Dec s to Mar s+1.

## What was looked at before this was written

Column names, row counts, month coverage and the longitude split of genesis points of the PR 38 tables (42,923 tracks, 398,611 fixes, Dec-Mar, 47 winters). The PR 40 results
table, row for this lead. No position was related to any index, and no 1979-2000 position was summarised.

## Question and answer rule

Question: in the Atlantic, in winters 1979-80..2000-01 (22 winters, never used for this lead), does the product of lagged NAO and lagged PNA shift the longitude of
the minimum-pressure position of all cyclones east, as in 2001-26?

**Prediction (fixed now): the same sign, positive.** Size: the nominal +1.11 deg per SD x SD, with a shrunken expectation of about +0.5 to +0.6 because of the winner's curse
(PR 40 itself gave +0.60 over 1979-2025 with an era term, which includes the source seasons and so is not independent). A null is an equally full answer.

## Data and variables (fixed)

- Tracks: PR 38 `tracks.csv.gz` and `fixes.csv.gz`. Unit = one track, dated at its first fix. Outcome = (lon, lat) of the fix of **minimum MSLP** (the earliest if tied).
  Basin = pipeline A's `track.basin(lat, lon)` applied **at that minimum-pressure fix** (Atlantic: 30-67N, lon >= 262 or <= 10). Atlantic longitude unwrapped (+360 if < 180). Tropical-cyclone
  transitions are not removed (as in the lead).
- Winters: 1979..2000 (primary). Genesis days 3 Dec - 31 Mar.
- Predictors: CPC daily NAO and PNA (`/mnt/project-files/teleconnection-test/cpc_indices/`), mean of days -10..-4 before the genesis day, standardised over the analysis days of the sample.
  Same-time contrast: days -3..+3.
- Model (the lead's, from `tele_combos/core.py::location_test`): OLS of (lon, lat) on A + B + A.B + month fixed effects + linear winter trend. Inference: sign-flip score test with winters as clusters
  (50,000 draws, seed 20261008), 95% interval by winter-pairs bootstrap (2,000 draws), CR1 SE alongside.

## Tests

**Primary (the only one that decides the answer): P1.** Atlantic, lagged NAO x PNA, longitude of minimum pressure, winters 1979..2000. Gamma, interval, one-sided p (the sign is predicted:
H1 gamma > 0; one-sided p = half the sign-flip two-sided p when gamma > 0, else 1 minus that half).

Secondary, labelled so, all reported; BH over the whole list including P1:
- S1 latitude of minimum pressure, same model. S2 joint lon+lat (2 df) p.
- S3 basin assigned at the **first fix** instead of the minimum-pressure fix (robustness of the basin rule).
- S4 same-time predictors (NAO/PNA days -3..+3): the contrast; circularity is the point, not a test of the lead.
- S5 **bridge (positive control, not independent):** the same model on the source era 2001..2013 winters (Dec-Mar, minimum-pressure position). It shows whether the lead exists with this
  position definition and this season where it was found. 2014 onward is left out on purpose so that no 2015-25 look is spent. If S5 shows nothing, a failure of P1 cannot be read as a failure to replicate.
- S6 **Pacific** NAO x PNA, longitude of minimum pressure, 1979..2000 (a control: no effect is predicted; Pacific lon 135-240).

Not run: gust-based position (barred pre-2004), Oct, Nov, Apr (not in the committed table), HF-equivalent lows (PR 40 found nothing, and the gust label is barred before 2004).

## Decision rules (fixed now)

SESOI = **0.55 deg** per SD x SD (half the lead's point estimate).

- **Replicates:** gamma > 0 in P1 AND one-sided p < 0.05 AND BH q < 0.10 over the 7 tests. (The prediction is directional, a two-sided test would punish the pre-registered sign.)
- **Does not replicate at face value:** the 95% interval's upper bound is below +1.11.
- **Well-powered null:** the whole 95% interval lies inside +/-0.55. This needs a null SE of about 0.28 or less.
- **Can't tell:** anything else, said in those words, with the minimum detectable shift (2.8 x null SE) and the interval stated.
- If S5 (bridge) is itself null, the report says the lead does not show up under this position definition even where it was found, and P1 is then reported as uninformative about the gust-peak lead.

## Power plan

Before P1 is fitted: season-block permutation of the index series (A and B moved together between winters, day of season kept), 400 permutations, as `tele_combos/power.py`. The spread of gamma is the
null SE; detectable shift at 80% power = 2.8 x SE. Reported for P1, S5 and S6 first, and written to `results/power.txt`.

## Looks

This is the 6th look at pre-2001 seasons (RA-10 was the 5th). No 2015-25 season is scored (S5 stops at 2013). Appended to `research/era5/hemispheric/results/heldout_looks.log` by timestamp when the result is closed.

## Deviations (post hoc)

(none yet)
