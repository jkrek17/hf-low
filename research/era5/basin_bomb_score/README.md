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
