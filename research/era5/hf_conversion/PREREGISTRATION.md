# Pre-registration, test 1: conversion composite (deepening storms that become HF against deepening storms that do not)

Committed before any converter and non-converter were compared on any environment variable or field. **ERA5 proxy, pipeline A**, seasons 2004-05 to 2025-26, Atlantic and Pacific separately, transitioning tropical cyclones in.
Plan from `hf_vs_storm/PLAN_THREE_TESTS.md` (test 1). Context: PR 85 found the hemispheric pattern acts on conversion of deepening storms to HF; this asks what the converting storms' environment looks like.

Question (Jason's test 1): among storms already deepening fast, which environment features separate those that reach HF within 24 h from those that do not? Answer forms: yes (q < 0.05, size), no (well-powered), or can't tell.

What had been looked at: counts only (3,316 deepening fixes on 2,274 storms; 1,287 converting fixes on 966 storms; about 450 matched pairs per basin expected) and the column list of `intensity/results/env_2004.csv.gz`. No environment value was compared between converters and non-converters.

## Population, groups, matching
- Population: 00/12 UTC fixes in `intensity/results/fixes_2004.csv.gz` with `g800` < 71.7 (not yet HF), 24 h normalised deepening rate `ndr24` >= 1.0 (the deepening class; Bergeron units) and `hf24` defined. **Converter**: `hf24` true (gust index reaches 71.7 kt within the next 24 h). **Non-converter**: `hf24` false.
- One fix per storm, drawn at random among its eligible fixes in the group (seed 20261012). A storm that has both converting and non-converting eligible fixes can appear in both groups (counted).
- Matching: 1:1 nearest neighbour on `ndr24` within basin x calendar month (pairs closer than 0.25 Bergeron), without replacement, converters matched in random order. Unmatched converters are dropped and counted. Pairs: all that match (stage A); stage B uses a random subset of the pairs sized below.
- Anchor time t0 = the fix. Lags: t0, t0 - 12 h, t0 - 24 h (same track; a missing lag fix is missing).
- Note that `ndr24` is the realised deepening over the next 24 h, so the groups are matched on how fast they deepen, not on how they looked beforehand. The question is what else differs.

## Stage A: environment table, no pull
- Variables from `env_2004.csv.gz` (1.5 degree, within 500 to 1000 km of the fix) joined at each lag: `jet250` (250 hPa jet maximum), `eady` (850-500 hPa Eady growth rate), `sstgrad` (SST gradient), `sst_t500` (SST minus 500 hPa temperature, instability), `tcwv` (column water vapour), `flux` (surface heat flux), `vadv500` (500 hPa vorticity advection), `div300` (300 hPa divergence). Missing (sst, sstgrad, sst_t500 over ice or land, 11-20%) drops the pair from that variable; counted.
- **Primary family (48 tests):** matched converter minus non-converter difference (stratified by month; season-block bootstrap of the 22 seasons, 2,000 resamples; paired structure kept by resampling seasons) for 8 variables x 2 basins x 3 lags, BH over 48. MDE = 2.8 x SE, in SD of the variable.
- Decision: a variable "separates" if q < 0.05 at the lag; "does not" only if MDE <= 0.2 SD, else can't tell. Also reported: the change over the 24 h (value at t0 minus value at t0 - 24 h), converter minus non-converter (secondary family, 8 x 2 = 16 tests, own BH).

## Stage B: 0.25 degree storm-relative fields (pull)
- Subset of pairs: pairs per basin = the largest N <= 250 whose pull totals at most 40 GB by HEAD requests (sizing script run before any field is pulled, and N written in the README before the pull starts); random draw seed 20261012. Times: t0, t0 - 12 h, t0 - 24 h for both members, that is 3 lags x 2 members x N pairs x 2 basins distinct times at most.
- Fields per time: 0.25 degree MSLP, 2 m dewpoint, 2 m temperature, SST (surface, 4 x about 2.4 MB); 250 hPa wind (u, v) and 850 hPa temperature and specific humidity at 0.25 degree by blosc-block range requests of single levels (about 2.5 MB each) from ARCO-ERA5.
- Low re-detected with pipeline A's rules, matched to the track position within 25 km (else dropped, counted); position at a lag by interpolation of the 00/12 track and snapping to a detected low within 150 km. Rotation by heading at t0.
- **Primary family (30 tests):** five scalars x 2 basins x 3 lags, stratified difference as above, BH over 30: (1) distance from the low to the 250 hPa wind maximum within 1,500 km; (2) cosine of that maximum's bearing relative to the low's motion (+1 ahead); (3) mean 850 hPa equivalent potential temperature within 500 km (from T850, q850 and 850 hPa); (4) near-surface stability, mean (SST minus 2 m temperature) within 500 km over ocean; (5) 850 hPa baroclinicity, mean magnitude of the horizontal 850 hPa temperature gradient within 500 km.
- Maps (descriptive): converter and non-converter composites rotated to motion, their difference with pixelwise BH stipple, at lags 0 and -24 h.
- The 800 km owned gust is not recomputed here (gust is not pulled); the 25 km position check stands in. Stated as a limit.

## Power, confounds and limits
- Matching is on deepening rate, basin and month, not on latitude, storm age or central pressure; those are reported per group with numbers, and a latitude-matched sensitivity (strata basin x 5 degree latitude band) is run for stage A.
- Hindsight caveat: converters are defined by their next 24 h, so a difference at t0 is an association usable as a forecast signal only because the variables are known at t0 or earlier; lags -12 and -24 h are before the label window.
- 22 seasons are the bootstrap unit; storms in the same week share an environment.
- Held-out looks: none scored, no model fitted. Transitioning tropical cyclones are in.

## Stage B size and N (written before the pull)

- Sized by HEAD requests (`results/sizes_surface.txt`: MSLP 2.17, dewpoint 2.37, t2m 2.34, SST 1.37 MB per time) and by three real single-level fetches (u250 3.6, v250 3.7, T850 2.6, q850 3.4 MB; `arco_level`), 21.5 MB per time. The plan assumed about 2.5 MB per level; the measured cost is 3.3 MB, a 30% higher level cost.
- **N = 250 pairs per basin** (the cap): 500 pairs, 1,000 members, 1,756 lag rows, **1,576 distinct times**, 33.9 GB (under the 40 GB cap; the lag rows are fewer than 3 x 1,000 because a storm's t0 - 12 h or t0 - 24 h fix often does not exist). Pairs drawn at random with seed 20261012 from the 788 stage A pairs (`results/conv_sample.csv`, `times_conv.csv`; `select_conv.py`).
- Whole-set pulls for the three follow-up tests so far: Test 2 15.3 GB (done), Test 3 stage B about 16.7 GB (sized, running), this test about 33.9 GB; about 66 GB against the 94 GB Jason approved.
- Dry run on 8 times passed (all 8 lows snapped, no all-NaN scalar). Wind speeds at 250 hPa are in m/s, theta-e is Bolton (1980), stability is SST minus 2 m temperature over ocean points (at least 50 in the 500 km disc), baroclinicity is |grad T850| in K per 100 km.

## Deviations (post hoc)

None yet.
