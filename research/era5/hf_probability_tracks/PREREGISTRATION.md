# P(HF) along ERA5 tracks back to 1979, and candidate missed observations: pre-registration

ERA5 PROXY, pipeline A (`research/era5/hf_history`). Written before any P(HF) value was computed or joined to
an archive match. What had been looked at: the file formats; the track counts per era; the chunk sizes of the
WeatherBench2 store (to size the pull). Nothing about P(HF) by track, nothing about match rates.

## Question (Jason, 2026-10-08)

"Plot the storm tracks with their probability of HF. There might be something to find there even if it's
missing observation." And, going back to 1979, "list the strongest storms in the ERA5."

In plain words the answer will be one of **yes / no / can't tell**, with how many storms, to the question:
*are there storms the near-storm model rates likely-HF where the archive has nothing?*

## Models (the PR 12 code, `research/era5/intensity/model.py`, refitted unchanged: logistic, L2, C = 1)

Target: pipeline A gust index reaches 71.7 kt within 24 h (`hf24`). Fit and score on seasons 2004-05 to
2025-26 only (decision 1), leave-one-season-out (LOSO) for every 2004+ number.

| Name | Predictors | Role |
|---|---|---|
| F | PR 12 `full` (state + Hart + environment, **includes the current gust index g800**) | 2004+ headline, as PR 12 |
| N | `full_nogust` (F without g800) | primary for every test below, 2004+ |
| S | `state_nogust` (storm state without g800: central pressure, 12 h tendency, age, latitude, speed, basin, day of year) | the only model usable before 2004 without a ~111 GB environment pull; also run on 2004+ so S and N can be compared |

Why a no-gust model before 2004: the gust index drifts upward at fixed depth before 2001 (decision 1), so a
predictor or target built on it would carry the drift into P. N and S use pressure-based information only.
The *target* still means "ERA5's gust index would reach 71.7 kt as it does in 2004+ at that state": P is
read as "this state is HF-like by the 2004+ standard", which is era-stable, and it is not a count.

Check before anything else: the refit of F must reproduce PR 12's LOSO BSS 0.423 (hf24) to +/-0.002 and HSS
0.59 to +/-0.01. If not, stop and explain.

## Sample

- 2004-05 to 2025-26: all 159,430 in-domain 00/12 UTC fixes (`intensity/results/fixes_2004.csv.gz`, env table).
- 1979-80 to 2003-04: the 00/12 UTC in-domain fixes of the 8,311 catalog tracks (events and matched null
  cases; `hf_history/results/era5_hf_catalog_tracks.csv`). **Not the full low population**: the other ~31k
  pre-2004 tracks have no committed fixes (re-extraction ~200 GB, not run). A null case is by construction a
  track that never reached HF, so "high P on a pre-2004 null" is a model false alarm by definition and the
  pre-2004 sample cannot be used for rates, only for ranking and mapping.
- Predictors for the pre-2004 fixes are built exactly as `fixes.py` builds them (motion from the previous 6 h).

## Track-level quantity

P_pre(track) = the largest LOSO P(hf24) over the track's in-domain 00/12 fixes that are not already HF
(`hf_now` false). Tracks with no such fix are listed apart and get no P_pre. Reported for F, N and S.
Candidate threshold: **P_pre >= 0.5** (fixed here, round, not tuned). Sensitivity at 0.3 and 0.7 is reported
and labelled as sensitivity.

## Archive matches

`research/era5/gust_depth/results/matches.csv` (PR 52): each archive event of seasons 2001-02 on is matched
to one pipeline A track (peak position within 500 km, track span within 12 h of the archive span, same
basin). A track is "archive-listed" if it appears there. I use it as is; I do not re-match. Its limit: the
match is on the *gust peak* position, so a track whose gust peak sits far from the archive position is
counted unlisted. Tests use seasons 2004-05 on (the archive is incomplete before).

## Tests (2004-05 to 2025-26, n = 22 seasons; season-block bootstrap 1,000 draws for intervals)

**T1. Does P see archive HF storms that the gust threshold misses?** Among tracks with gust index < 71.7 kt
and a P_pre, is P_pre higher for archive-listed tracks? Statistic: AUC of N's P_pre for archive-listed vs
not. Rule: **yes** if the 90% interval of AUC excludes 0.5 and the Pmax >= 0.5 listed rate is at least twice
the rate below 0.5; **no (well-powered)** if the upper end of the AUC interval is at most 0.55; otherwise
**can't tell**. Same test for S.

**T2. What are the unlisted high-P tracks?** Take tracks with N P_pre >= 0.5 and gust >= 71.7 kt (proxy
events), 2004+. Split into archive-listed and unlisted. Remove from the unlisted group, in this order, the
known false-alarm sources: (a) tropical-cyclone-linked (within 400 km of an IBTrACS point at the same time at
any fix; `tc_candidates.csv`), (b) Atlantic with the gust peak at or north of 60N (the lower-confidence
Greenland-barrier group), (c) tracks with fewer than 3 in-domain fixes. What remains is the **residual**.
Report the residual count, its share of the high-P events, and by basin and season.

**T3. Does the residual look like the listed events or like a proxy artefact?** Compare residual with listed
events on four measures: gust margin over 71.7 kt, P_pre, number of in-domain fixes, minimum pressure.
Mann-Whitney with season-block permutation, four tests, Benjamini-Hochberg q. Reading, fixed now:
- no measure with q < 0.05 and the residual has at least 30 tracks: **consistent with missed observations**
  (indistinguishable from listed storms; it still cannot be proven without a third source);
- two or more measures differ at q < 0.05 in the direction of weaker/shorter: **looks like a proxy artefact**;
- anything else, or fewer than 30 residual tracks: **can't tell**.

**The answer to the headline question** (yes / no / can't tell) is built only from T1-T3:
- *yes* needs T1 yes **or** (T2 residual of at least 30 tracks and T3 consistent with missed observations);
- *no* needs T1 no and T2 residual below 10 tracks or T3 artefact;
- otherwise *can't tell*. Even a *yes* means "storms the archive does not list that look like listed storms",
  not "confirmed missing observations": ERA5 alone has no independent wind observation, and this document
  says so wherever the answer is quoted.

**Pre-2004 (no archive, no test).** The 1979-2003 catalog tracks with S P_pre >= 0.5 are listed and mapped
as candidates, with the 2004+ listed rate at the same S P_pre as the calibration of how likely the archive
would have listed them. Counts per season before 2001 are not compared with later seasons (decision 1).

## Strongest storms (descriptive, no test)

From `all_tracks.csv.gz` (75,087 tracks, 1979-80 to 2025-26), per basin, top 25 by:
1. minimum ERA5 central pressure (`minp`);
2. background-adjusted depth, `minp` minus the calendar-month MSLP climatology at the deepest fix
   (`tele_intensity/results/event_table.csv`, `anom_clim`); available for catalog events only, so I first
   check that all top-25 by `minp` are catalog events and say so if not.
Peak gust shown, flagged for seasons before 2001-02. Storm names only where date and position match a
known event, written as inferred. ERA5 central pressure in the 1980s over open ocean was constrained by
fewer observations, so depths may be too shallow rather than too deep.

## Not done / cannot show

- Whether any unlisted storm really had hurricane-force wind. That needs scatterometer, ship or buoy data
  or the OPC text, none reachable here.
- P for the other ~31k pre-2004 lows (needs a re-extraction) and the environment-model P before 2004 (needs
  the 111 GB pull; Jason has been asked).
- Any HF level or trend across 2001.

## Deviations (post hoc)

(none yet)
