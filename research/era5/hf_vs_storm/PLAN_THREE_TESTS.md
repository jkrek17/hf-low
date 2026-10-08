# Plan: three follow-up composite tests (run on Jason's go; results in hf_lifecycle, hf_pattern_phase, hf_conversion, PR 133)

ERA5 proxy, pipeline A, seasons 2004-05 to 2025-26, Atlantic and Pacific separately, transitioning tropical cyclones in. Each test gets its own PREREGISTRATION.md committed before any field is read.
Common rules: season-block bootstrap (2,000 resamples of 22 seasons) for every difference; Benjamini-Hochberg q < 0.05 within each test's family and across grid points within each map;
well-powered null reported as a result; lagged environment at days -10 to -4 where a pattern index is used; ERA5 reads low in extreme storms. Sample sizes below are counted from committed tables (labels only); detectable differences are estimates to be confirmed by bootstrap SE (MDE = 2.8 x SE).
Pull sizes are estimates from the measured bytes per time of the running composites (0.25 degree gust 3.2 MB, MSLP 2.2, u10 3.3, v10 3.3, 2 m dewpoint 2.4 per time); each will be sized by HEAD requests before pulling. Every pull stays under 50 GB or stops for a card.

## Test 1: Conversion composite (deepening storms that become HF vs those that do not)
- **Hypothesis:** among storms already deepening fast, the ones that reach HF within 24 h sit in a different environment (jet-streak position, low-level baroclinicity, near-surface stability, 850 hPa moisture), which is the mechanism behind the main result that the pattern acts on conversion.
- **Groups and matching:** fixes with 24 h Bergeron deepening rate (ndr24) >= 1.0, not yet HF (g800 < 71.7). Converters (hf24 true) vs non-converters. Counted in `intensity/results/fixes_2004.csv.gz`: 3,316 fixes on 2,274 storms (Atlantic 1,475, Pacific 1,841); 1,287 fixes (966 storms) convert. One fix per storm, matched on ndr24 (+-0.25), basin and calendar month, 1:1 without replacement.
- **Fields and times:** at the fix and 12 h and 24 h earlier. Storm-relative (rotated to motion) 0.25 degree MSLP, 2 m dewpoint, 2 m temperature, SST; 1.5 degree 250 hPa wind (jet-streak position relative to the low), 850 hPa theta-e and the 850-1000 hPa baroclinicity; stability = SST minus 2 m temperature.
- **Statistic and field significance:** stratified difference converters minus non-converters, season-block bootstrap per pixel, BH within each map. Primary family: 6 scalar summaries (jet-streak distance and bearing, baroclinicity, stability, theta-e, 850 hPa moisture flux) x 2 basins x 3 lags, BH over 36.
- **Sample and detectable difference:** about 450 matched pairs per basin available (limited by non-converters at high ndr24); pulling 250 pairs per basin. Expected MDE about 0.2 to 0.3 storm-to-storm SD per scalar; a smaller difference is "can't tell".
- **Pull size:** 1,000 storms x 3 times x (4 surface fields ~10.5 MB) at 0.25 degree is about 31 GB, plus the 1.5 degree WeatherBench2 pressure-level chunks for the jet and theta-e (about 8 GB, from the earlier RA-20 pull rate). About 39 GB.

## Test 2: Pattern-phase composite (HF lows forming under high vs low hemispheric-pattern index)
- **Hypothesis:** the PR 41/85 pattern changes the structure of HF lows (size, intensity, asymmetry, moisture), not only how many storms convert.
- **Groups and matching:** all HF onsets (2,000: Atlantic 1,104, Pacific 896) split into top and bottom tercile of the out-of-sample pattern index (lagged days -10 to -4 before the track's first fix, the PR 85 index), within basin and calendar month. About 370 per tercile in the Atlantic and 300 in the Pacific. Middle tercile shown, not tested.
- **Fields and times:** storm-relative gust, wind speed, MSLP, 2 m dewpoint at HF onset (the same measures as the running composites), plus structure scalars (central MSLP, gradient, radius of maximum gust, 48 kt outer radius, right-of-motion share, dewpoint, gust factor).
- **Statistic and field significance:** high minus low tercile, season-block bootstrap, BH over 7 scalars x 2 basins = 14 (primary) and per pixel within maps. Because the index has about one independent value per week, seasons are the resampling unit.
- **Sample and detectable difference:** 300 to 370 per tercile gives MDE about 0.3 to 0.4 SD per scalar after the season bootstrap; "no difference" only below 0.25 SD, otherwise "can't tell".
- **Pull size:** the running pull already holds 800 HF onsets (400 per basin); the remaining 1,200 onsets need 5 fields x 14.4 MB: about 17 GB. If the running onsets suffice (about 130 per tercile per basin, MDE roughly 0.5 SD) the test needs no new pull.

## Test 3: Life-cycle time-lapse (when do HF and storm-force lows first differ?)
- **Hypothesis:** HF and storm-force-only lows are indistinguishable early and separate at an identifiable lead time before peak, which is the lead a forecaster can use.
- **Groups and matching:** HF vs storm-force-only (54 to 71.7 kt), anchored at peak gust time, matched on basin and month as in the running composites.
- **Fields and times:** every 12 h from -48 h to +24 h around peak (7 times). Stage A uses committed track tables only (central pressure, 12 h deepening, g800, size from `fixes_2004`): no pull. Stage B: storm-relative gust and MSLP at 0.25 degree.
- **Statistic and field significance:** per lag, stratified difference with season-block bootstrap and BH within each map; the headline is the first lag at which a pre-defined scalar family (central MSLP, 12 h deepening, gradient, gust-area 48 kt) differs at BH q < 0.05, with the bootstrap interval on that lag (resampling storms' lags jointly). The "first lag" is defined and fixed before looking. Lead times are reported as differences in lag, not as forecast skill.
- **Sample and detectable difference:** stage A, all 2,000 HF and 8,136 SF (MDE below 0.1 SD); stage B, 250 per group per basin (MDE about 0.25 SD per scalar).
- **Pull size:** stage B is 7 lags x 1,000 storms x (gust + MSLP 5.4 MB) = about 38 GB. Stage A is free.

## Order and total
Test 3 stage A needs no pull and can start first; test 2 may need none. If all three run, the pulls total about 39 + 17 + 38 = 94 GB, so they will run one at a time, each under 50 GB, and the combined plan needs your go-ahead only if any single pull passes 50 GB (none does).
