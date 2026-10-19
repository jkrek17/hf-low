# Gale / storm / HF shares by pattern-index quintile (pre-registration)

Asked by Jason 2026-10-19 (basin-wide "bomb potential" score; option A). Committed before any class count was joined to the index.

Question: does the PR 41 basin pattern index (leave-one-season-out, `hemispheric/results/oos_index_{atl,pac}.csv`) separate quiet weeks from loaded ones for gale, storm and HF lows, beyond what the calendar month already says?

**ERA5 PROXY, pipeline A** (`hf_history/results/all_tracks.csv.gz`), seasons 2004-05 to 2025-26, Oct-Apr weeks (Oct 1 + 7k days, k 0-29), Atlantic and Pacific separately, tropical transitions in. Track assigned to the week of its peak time. Class by track peak gust index: HF >= 71.7 kt, storm 54 to <71.7 kt (as in `hf_vs_storm`), gale 40 to <54 kt (**40 kt is a new, uncalibrated cut chosen here before any count; no archive gale series exists**).

Index: PR 41 out-of-sample index, quintiles per basin from the 22-season weekly index. Predictors read days -7..-1 before the week (lag already built in). No refit, no new model. This is the 26th-plus look at the PR 41 index across all 22 seasons (logged in `looks/pattern_quintiles.log`); not a fresh hold-out.

Primary outputs per basin and class (3 x 2 = 6 tests): (1) mean tracks per week by quintile; (2) ratio top / bottom quintile, and month-adjusted ratio (observed / expected from the month mean, pooled over the quintile); (3) the same as share of all tracks in the quintile ("N in 10"). Season-block bootstrap (2,000 draws) intervals. BH across the 6 month-adjusted top/bottom ratios.

Decision rule: "separates" if the month-adjusted top/bottom ratio is >= 1.5 with lower interval bound > 1 and q < 0.05; "weak" if the lower bound > 1 but ratio < 1.5; otherwise "can't tell". Power is not simulated; effective n is 22 seasons, so intervals are the guide.

Known conflict from earlier work: index and latitude/season move together (`hf_pattern_phase`), hence the month adjustment. Not done: position control, GFS-fed index.

## Deviations (post hoc)
(none yet)
