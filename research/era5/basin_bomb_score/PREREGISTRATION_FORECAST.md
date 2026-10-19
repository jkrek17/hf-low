# Basin pattern score from GEFS forecasts against bombs and HF lows (pre-registration)

Asked by Jason 2026-10-19 ("try forecast") after the ERA5 bomb result (PR 275). Committed before any bomb or HF count was joined to the GEFS index. The bomb counts and the ERA5-index results are already known (PR 275); the GEFS index has been scored only for archive HF counts (`gfs_test` stage B).

Question: if the pattern score is computed from forecast fields instead of reanalysis, does it still tell which weeks will have more bombs and more HF lows?

**Forecast input.** `gfs_test/results/index_zu.csv` (committed, no new pull): the z500 + u250 part of the PR 41 index (frozen weights fitted on ERA5 2004-05 to 2014-15, SST block dropped) computed from the GEFS v12 control reforecast, 00 UTC initialisation on day S-8, forecast days 1-7 covering the PR 41 window S-7..S-1, so it is available 8 to 14 days before the weeks it scores. Seasons 2004-05 to 2018-19 (15 seasons, 450 weeks per basin). Column `atl_gefs`, `pac_gefs`; reference `atl_era5`, `pac_era5` is the z+u-only index on ERA5 (not the full PR 41 index).

**Outcomes.** Weekly bomb count exactly as in `PREREGISTRATION.md` of this folder (pipeline A proxy, fix at 30N or more, 24 hPa/24 h Bergeron-scaled, one per track, first bomb fix week). Weekly HF count from pipeline A (gust index >= 71.7 kt, peak-time week), same weeks. ERA5 PROXY truth, not the archive.

**Tests (per basin, per outcome; 2 x 2 families).**
- F1: month-adjusted top / bottom fifth ratio with the GEFS index, season-block bootstrap (2,000).
- F2: leave-one-season-out Poisson deviance skill of the GEFS index over month + previous-week count (as T2), season-block permutation p (1,000).
- F3 (reference, not a test): the same F1 and F2 with the ERA5 z+u index on the same weeks, and the GEFS minus ERA5 skill difference with a season bootstrap.
- F4 (guard): skill on 2015-16 to 2018-19 only (4 seasons, the part the frozen weights never saw), reported with its interval, no decision rule attached; power is poor.
BH over F1 and F2 for the four basin-outcome pairs (8 tests).

**Decision rules.** Forecast score "works" for an outcome if F2 skill > 0 with q < 0.05 and the F1 lower bound > 1 in that basin. "Loses most of it" if the GEFS-minus-ERA5 skill difference interval excludes 0 on the negative side. Otherwise "can't tell". No power claim for a "no": 15 seasons.

**Limits stated now:** the index weights were fitted on 2004-14 ERA5 discovery seasons, so 2004-14 skill is in-sample for the weights (as in stage B deviation 7); a GEFS control member only, not an ensemble; no GFS (2021 onward) because the frozen weights and baseline would need a separate freeze; 00 UTC init only; weekly, not daily; seasons 2019-20 onward not covered.

## Deviations (post hoc)
1. (before any association was read) BH is applied to the four F2 permutation p-values only; F1 has no p-value and is judged by its interval, as in the main bomb plan's decision rule. The "8 tests" in the Tests section is therefore 4.
