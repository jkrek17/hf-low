# P(HF within 24 h) along ERA5 tracks back to 1979, and the strongest storms (ERA5 proxy, pipeline A)

Pre-registered in `PREREGISTRATION.md` (commit `6fcddf1`) before any P(HF) was computed. Code: `score.py` (per-fix
probabilities), `analyse.py` (tests, strongest storms), `figs.py` (figures). Outputs in `results/`. Reads committed files
only; nothing here feeds `docs/` or the site.

**This is a proxy.** "HF" is pipeline A's gust index reaching 71.7 kt in ERA5 (`research/era5/hf_history`), not the
archive. P is the PR 12 near-storm model's probability of that within 24 h. Tropical cyclones are in.

## Answers in plain words

**Are there storms the model rates likely-HF where the archive has nothing? Mostly no, and the few that might be
cannot be confirmed from ERA5.**

- *Where the archive exists* (2004-05 to 2025-26): 848 proxy events reach P >= 0.5 before they turn HF (model N, no gust
  predictor). 632 are archive-listed. Of the 216 unlisted, removing tropical-cyclone-linked (23) and Atlantic-north-of-60N
  (38) leaves **155 (18%)**. They are not like the listed storms: median 3.3 kt over the threshold against 8.2 kt,
  shorter (16 against 18 fixes) and shallower (963.9 against 957.6 hPa), all four at q < 0.001. That is what threshold-edge
  proxy storms look like, not what missed strong storms would look like. **Only 9 are 10 kt or more over the threshold
  under both ways of matching to the archive**, and one of those reaches 950 hPa (`results/residual_unlisted_both_rules_2004.csv`).
  At least one of the nine is a tropical cyclone the flags miss (track 249925, 29 Aug 2005, Gulf of Mexico, inferred to be Katrina).
- *Before 2004* there is no archive, so nothing can be tested. Jason approved the 111 GB environment pull, so the 1979-2003
  catalog tracks are scored with model N too. 962 of them reach P_pre >= 0.5 (864 proxy events, 98 matched null cases;
  `results/tracks_pre2004_summary.csv`; storm-state-only model S gives 706). The 2004+ calibration says how to read them: of
  catalog-type tracks at the same P, **74.5% of events and 15.7% of null cases were archive-listed**. So of the 98 pre-2004
  null cases rated >= 0.5, about 15 would have been archive storms had the archive existed. That is the extent of "something
  to find" before 2004 in this sample; the 40 tropical-cyclone-linked tracks among the 962 are in the file with a flag.
- **What the pre-registered rule returned, and why I do not lean on it.** T1 passed, so the written combination rule
  returns "yes". T1 shows P sees archive-listed storms that ERA5's gust field misses: among 33,383 sub-threshold tracks,
  P_pre >= 0.5 picks out 20% archive-listed against 1.1% below (AUC 0.911, 90% interval 0.899-0.922). That is the reverse
  of a missed observation (the archive has it, ERA5's gust does not), so I consider the rule badly specified. Logged as a
  post hoc deviation. Example, inferred from date and position: the Christmas Eve 2013 Atlantic storm (track 330264, ERA5
  minimum 928.9 hPa, gust index 70.7 kt, one knot under the threshold, P_pre 0.86, archive-listed).

Nothing in ERA5 alone shows a hurricane-force wind that an observation missed. Confirming any of the nine, or any
pre-2004 candidate, needs scatterometer, ship or buoy data or the OPC text, none reachable from here.

## Skill of the models (leave-one-season-out, 22 seasons 2004-05 to 2025-26, 159,430 fixes, hf24)

| | BSS vs basin-month climatology | HSS | note |
|---|---|---|---|
| F, PR 12 full | 0.4228 | 0.59 | reproduces PR 12 (0.423, 0.59) |
| N, full without the gust index | 0.3804 | 0.56 | used for 2004+ tests |
| S, storm state without the gust index | 0.2624 | 0.46 | before-2004 fallback, kept for comparison |

Dropping the gust index costs 0.04 BSS with the environment and 0.11 without it (PR 12's `state` with gust is 0.370).
S is weak (its P correlates 0.85 with N on 2004+ catalog fixes, 0.82 at track level before 2004). Before 2001 the gust
index drifts upward at fixed depth (decision 1), so no gust predictor is used there; P before 2004 reads as "this
state is HF-like by the 2004+ standard", not as a count and not comparable across 2001 as a level.

## What was run

- **Plots** (`results/figs/`, also in the project folder): Fig 1 1979-2003 catalog tracks coloured by N (Fig 1b: by S); Fig 2 all 2004-2025
  tracks coloured by N; Fig 3 where the gust peaks of listed and unlisted high-P events sit; Fig 4 archive-listed rate by P
  and how the unlisted events differ; Fig 5 the 25 deepest lows per basin.
- **Sample before 2004 is the catalog only** (events and one matched null per event: 4,311 of about 39,600 tracks, 30,388 fixes).
  The other lows have no committed fixes; re-extracting them is about 200 GB. Null cases are non-events by construction, so
  the pre-2004 sample can be ranked and mapped but not used for rates.
- **Per-fix tables for the map thread:** `results/fix_probs_2004.csv.gz` (columns: track, time, basin, season, lat, lon, msl,
  g800, hf_now, hf24, P_F, P_N, P_S, cut_*; P is leave-one-season-out; 00/12 UTC in-domain fixes only) and
  `results/fix_probs_pre2004.csv.gz` (P_N, P_S, and P_Sg with the gust index as a labelled within-era contrast; role =
  event or null_case).
- **Tests.** T1 (above). T2 residual 155 of 848 (atl 91, pac 64), 3-10 per season, every season has some. Sensitivity: P_pre >= 0.3
  gives 253 residual of 1,230; >= 0.7 gives 67 of 464; model S gives 97 of 536. T3 four tests, all q < 0.001 (permutation
  p at its floor, 1/2001). Post hoc, looser matching (any archive HF fix within 400 km of any 00/12 fix): residual 145, and
  106 unlisted under both rules.
- The share of proxy events that are archive-listed over 2004-2025 is 61% (peak-position rule) or 62% (loose); the
  calibration season window gave POD 0.77, so a share of the unlisted events is the proxy's known surplus (ERA5 counts 6-19%
  more events than the archive).

## Strongest storms (descriptive; `results/strongest_storms.csv`, per basin top 25 by minimum ERA5 central pressure and by depth against the monthly climatology)

| | Atlantic | Pacific |
|---|---|---|
| Deepest (min ERA5 MSLP) | 912.3 hPa, 15 Dec 1986 (59.8N 34.8W) | 922.6 hPa, 31 Dec 2020 (52.0N 173.3E) |
| 2nd | 914.9 hPa, 10 Jan 1993, Braer storm (inferred) | 923.5 hPa, 30 Oct 1989 |
| 3rd | 920.3 hPa, 15 Feb 2020, Storm Dennis (inferred) | 923.9 hPa, 3 Oct 1981 |
| 25th | 932.1 hPa | 938.2 hPa |

- Of the Atlantic 25 deepest by pressure, 11 are before 2001-02 and 12 from 2004-05 (Pacific 12 and 11). None of that is a
  trend claim. ERA5 pressure in the 1980s was constrained by fewer observations, so older depths may be too shallow.
- Names are inferred from date and position only: Braer (Jan 1993, ERA5 914.9 hPa), Storm Dennis (Feb 2020, 920.3), ex-Hurricane
  Fiona (24 Sep 2022, 934.2) and the Bering Sea bomb of Nov 2014 (ex-Typhoon Nuri, 925.5). Not checked against any record in
  this repository.
- Ranking by depth against the monthly climatology favours autumn storms and tropical transitions: 4 of the Atlantic and 7 of
  the Pacific top 25 are tropical-cyclone-linked. One of the Atlantic 25 deepest by pressure (track 330264, Dec 2013) is not a
  catalog event, so it has no climatology-adjusted depth.
- Peak gust is in the file with a flag for seasons before 2001-02, where the gust index drifts upward at fixed depth.

## Not done

(The environment pull for the catalog fixes was run after Jason's go-ahead: 2,848 chunks, about 111 GB; derived table `results/env_pre2004.csv.gz`, 1.4 MB.)

- P for the other ~31,000 pre-2004 lows (a re-extraction of about 200 GB) and so any complete 1979-2003 near-miss list.
- Any independent wind observation of a candidate.

## Independent check

See the verifier comment on the pull request. Numbers in this README that were not recomputed independently are listed
there.
