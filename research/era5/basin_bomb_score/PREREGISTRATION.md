# Basin pattern score against bombs (24 hPa in 24 h): pre-registration

Asked by Jason 2026-10-19 ("what about a pattern recognition type engine that scores the basin wide pattern for bomb potential?", then "you can explore the bomb thing" after the quintile table, `pattern_quintiles/`, PR 273). Committed before any bomb count was joined to any score. Existing inventory: bomb fixes already exist in `intensity/results/fixes_2004.csv.gz` (`ndr24`, Bergerons); no weekly bomb table exists.

Question: does a basin-wide weekly pattern score tell a forecaster which weeks will have more bombs, and is a score refit for bombs any better than the HF-trained PR 41 index?

**ERA5 PROXY, pipeline A**, 22 seasons 2004-05 to 2025-26, Oct-Apr weeks (Oct 1 + 7k days, k 0-29), Atlantic and Pacific separately, tropical transitions in, 00/12 UTC fixes only (some bombs between those times are missed; a 6-hourly pull is not planned).

**Outcome.** Bomb = a pipeline A fix at 00 or 12 UTC with `ndr24` >= 1 Bergeron (24 hPa in 24 h scaled by sin(lat)/sin 60; `intensity/fixes.py`) at latitude >= 30N. One count per track, in the week of its first bomb fix. Weekly bomb count per basin.

**Scores.**
- S0 (primary): the PR 41 leave-one-season-out index (`hemispheric/results/oos_index_{atl,pac}.csv`), trained on HF counts, no refit.
- S1 (secondary, bomb-refit): Poisson ridge regression of weekly bomb counts on the 30 PR 41 pattern components (`z1-10,u1-10,s1-10`, `hemispheric/results/weekly_table.csv.gz`, lag-1 window), ridge strength chosen by inner season-block cross-validation inside each leave-one-season-out fold (nested), out-of-fold score for every season.

**Tests (per basin, 2 basins).**
- T1: month-adjusted top / bottom quintile ratio of weekly bomb counts for S0 (as in `pattern_quintiles`), season-block bootstrap (2,000).
- T2: leave-one-season-out Poisson deviance skill of S0 added to a baseline of month + previous-week bomb count, with a season-block permutation p (2,000).
- T3: same skill for S1 and the difference S1 minus S0 (bootstrap by season).
- T4: HF versus bomb on the same weeks: skill of S0 for the weekly HF count (reproduces PR 41) next to T2, so the two are on the same scale.
Multiplicity: BH over T1 and T2 in both basins (4 tests); T3 and T4 reported unadjusted as secondary.

**Decision rules.** "Scores bombs" if T2 skill > 0 with q < 0.05 and the T1 month-adjusted lower bound > 1. "Weak" if only one holds. "Can't tell" otherwise; "no" only if a planted-effect run says 80% power at a per-SD rate ratio of 1.15 and the interval excludes it. S1 "beats S0" if the T3 difference interval excludes 0.

**Allowed before the association is read:** marginal bomb counts per week and a planted-effect power run (`power.py`) on the real week structure. No index-outcome association is computed before this file is committed.

**Not done in this pass:** a GFS-fed version (needs GFS pressure tracks, `gfs_test/results/gfs_tracks_*.csv.gz` exist; would be a follow-up with its own freeze), a day-level score, position control, pre-2004 seasons (decision 1; pressure-only bomb counts could go back to 1979 later).

Known conflict from earlier work: the HF index moves with month and latitude (`hf_pattern_phase`); the baseline carries month. The Bergeron scaling inflates rates at low latitude, hence the 30N floor (chosen here, before any count). This is another look at PR 41's index on all 22 seasons (log `looks/basin_bomb_score.log`).

## Deviations (post hoc)
(none yet)
