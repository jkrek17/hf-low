# Does the basin pattern score tell a forecaster which weeks have more bombs?

**ERA5 proxy, pipeline A**, 22 seasons 2004-05 to 2025-26, Oct-Apr weeks, tropical transitions in, 00/12 UTC fixes only. Plan: `PREREGISTRATION.md` (`699bb6a`), committed before any bomb count met a score. Bomb = a fix at latitude 30N or more with 24 h deepening of at least 1 Bergeron (24 hPa in 24 h scaled to 60N, `intensity/fixes.py`); one count per track, in the week of its first bomb fix. Totals inside the 30 Oct-Apr weeks: 921 bombs Atlantic, 1,136 Pacific (a further 111 and 108 first-bomb fixes fall after week 29 and are not counted).

## Answer in plain words
**Yes, but weakly, and much less than for hurricane-force lows.** The same basin pattern score that picks out HF-heavy weeks (PR 41, trained on HF lows) also picks out bomb-heavy weeks, mostly in the Atlantic.
- Atlantic: bottom fifth of weeks 1.1 bombs per week, top fifth 1.8; month-adjusted ratio **1.6** [1.40, 1.85]. Weeks with at least one bomb: 69% (bottom) to 86% (top).
- Pacific: 1.6 to 2.0 bombs per week; ratio **1.2** [1.09, 1.39]. Weeks with at least one bomb: 81% to 91%.
- For comparison HF lows rose about 2 times in the top fifth (`pattern_quintiles`). A bomb week is common whatever the pattern, so the score moves the odds a little.
- Out-of-season-sample skill over month and last week's count: Atlantic +2.0% [+0.3, +3.6], Pacific +1.1% [-0.1, +2.3] (permutation p 0.001 and 0.002, BH q 0.0013 and 0.0020). On the same weeks and the same method the score gains +10.6% [7.7, 12.9] and +4.5% [0.6, 7.7] for HF counts (proxy counts, not comparable to PR 41's archive numbers).
- Rule verdict: "scores bombs" in both basins by the registered rule; the Pacific skill interval touches zero, so treat it as marginal.
- **A score refit for bombs adds nothing.** S1 (30 pattern components, ridge Poisson, nested leave-one-season-out) minus the HF-trained score: Atlantic +0.6 points [-2.1, +3.1], Pacific -0.5 [-2.3, +1.4]. The ridge strength chose the largest value on the grid (300) in every fold, so a heavier penalty was not tried (a limit, not a finding).

Verified by a fresh Sonnet agent from the committed fixes and index files: totals, bombs per week and weeks with a bomb by fifth, the month-adjusted ratios and the two skills all matched. Not independently checked: bootstrap and permutation intervals, the ridge refit (S1), the HF comparison (T4), the power run.

Power (planted effect, before any association was read, `results/power.json`): at a rate ratio of 1.15 per SD the T2 test finds it 98% (Atlantic) and 99% (Pacific) of runs; false positives 5% and 3%. Upper bound, since draws ignore season-to-season variation in the effect.

## Limits
- Another look at the PR 41 index on all 22 seasons (`looks/basin_bomb_score.log`); not a fresh hold-out.
- Weekly, not daily: it does not say which day or where. 00/12 UTC fixes miss bombs that peak between; ERA5 under-resolves the strongest deepening. The 30N floor and the bomb-per-track rule were chosen in the plan before any count.
- Proxy storms only; no GFS-fed version yet (`gfs_test/results/gfs_tracks_*.csv.gz` hold GFS pressure tracks, so one could be built).
- Month and last week's count are the only controls; no position control.

Files: `bomblib.py`, `power.py`, `analysis.py`, `results/tests.csv`, `results/extra.json`, `results/power.json`.

## Forecast version (GEFS-fed score; Jason, 2026-10-19 "try forecast")
**ERA5 proxy truth, pipeline A**, 15 seasons 2004-05 to 2018-19, 450 weeks per basin. Plan: `PREREGISTRATION_FORECAST.md` (`6e676f7`, one deviation logged before scoring). No new data: the index is `gfs_test/results/index_zu.csv`, the PR 41 z500 + u250 pattern computed from GEFS v12 control forecasts for days 1-7 issued on day S-8, so a week is scored 8 to 14 days before it starts. Reference: the same index computed from ERA5.

### Answer in plain words
**The score works just as well from a forecast as from reanalysis, but the size is inflated and the clean test is too small to confirm it.**
- Top vs bottom fifth of weeks, month-adjusted, GEFS index: bombs Atlantic **1.8** [1.43, 2.12], Pacific **1.35** [1.13, 1.60]; HF lows Atlantic **2.3** [1.76, 2.79], Pacific **1.9** [1.47, 2.53].
- Leave-one-season-out skill: bombs +4.8% [+2.5, +6.9] Atlantic and +2.5% [-0.2, +4.8] Pacific; HF lows +10.7% [6.3, 14.5] and +6.2% [2.5, 9.8]. Permutation p 0.001 for all four (BH q 0.001).
- Forecast minus reanalysis (same weeks, same weights): bombs -0.2 points [-1.1, +0.8] Atlantic and -0.2 [-1.1, +0.7] Pacific, so nothing lost; HF lows -2.0 [-3.2, -0.9] Atlantic (lost about 1 in 6 of the skill) and -0.6 [-1.8, +0.7] Pacific.
- Rule verdict: "works" for all four basin-outcome pairs; "loses most of it": no.
- **Caution:** the score's weights were fitted on ERA5 HF counts in 2004-05 to 2014-15, so 11 of the 15 seasons are in-sample for the weights, which is why these numbers are larger than the ERA5 bomb result above (+2.0%, +1.1%). The four seasons the weights never saw (2015-16 to 2018-19) give, for GEFS: bombs +1.4% [-3.1, +5.2] Atlantic and +4.6% [-0.2, +11.4] Pacific; HF lows +5.1% [+0.5, +11.4] and -1.2% [-8.9, +6.1]. Four seasons is too few to say; only the Atlantic HF interval excludes zero. Same pattern with the ERA5 input, so the forecast is not the weak link.

### Limits
- One GEFS control member, 00 UTC initialisation, weekly, 8 to 14 days ahead; it says nothing about which day or where. No GFS-fed version (2021 onward needs its own freeze).
- Seasons 2019-20 onward not covered; the 2015-18 look is the only clean one and it is small.
- Another look at the same outcomes: the bomb and HF counts for 2004-2018 were already scored with the ERA5 index (PR 275 and `pattern_quintiles`). Log: `looks/basin_bomb_score.log`.

Files: `forecast.py`, `results/forecast_tests.csv`.
