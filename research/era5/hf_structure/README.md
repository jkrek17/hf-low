# Where the hurricane-force winds sit, and how big the area is (ERA5 proxy, pipeline A)

Pipeline A = `research/era5/hf_history`. Everything here is an ERA5 **proxy**: ERA5
instantaneous 10 m gust and MSLP around the lows pipeline A tracked, at the times
pipeline A's gust index reaches its calibrated threshold of 71.7 kt. It is not a
measurement of hurricane-force wind. The plan, choices and decision rules are in
`results/PLAN.md`, committed before the full run (`2076a8d`, `cf4615b`). The
deviations log there is empty: nothing was changed after the results were seen.

## Files

| File | What |
|---|---|
| `extract.py` | Fetches gust and MSLP at each time, re-detects the lows with pipeline A's own rules, computes the per-fix statistics and motion-relative / north-up boxes. Raw boxes stay in `$ERA5_WORK` (not committed). |
| `composite.py` | Life-cycle stage, Hart phase, composites, storm-resampled intervals, figures. |
| `results/fixes.csv` | One row per HF-strength fix (5,983 rows, 1.9 MB). |
| `results/structure.txt` | Every number quoted below, plus every comparison in the plan. |
| `results/composite.npz`, `composite_motion.png`, `composite_northup.png` | Composite fields and figures. |
| `results/pull_bytes.csv` | Bytes streamed per time. |

## Sample and pull

Every in-domain point of a pipeline A event track with g800 >= 71.7 kt, seasons
2004-05 to 2025-26: **5,983 fixes at 5,311 distinct times on 2,000 storms**
(Atlantic 3,303 fixes / 1,103 storms; Pacific 2,680 / 897). The agenda's 6,075
times counted points outside the basin boxes, where pipeline A's index is
undefined; those are dropped. **Pull: 33.3 GB** for 5,311 times (gust + MSLP from
ARCO-ERA5, plus WeatherBench2 10 m wind speed at 1,500 of them). The agenda's
22 GB used 3.7 MB per time, measured 6.3 MB here. Under the 50 GB gate. No existing
cache covered it: the near-storm thread's fields are 1.5 degree pressure-level
data at 338,511 fixes, with no 10 m wind or gust.

**Reproduction check.** Re-detecting lows and recomputing the 800 km owned gust
reproduces the catalog's g800 at 99.98% of fixes (max difference 0.1 kt) with the
same centre at every fix.

## Results (intervals: 90%, resampling storms)

Pipeline A, ERA5 proxy, median unless stated. Area = owned ocean area with gust
>= 71.7 kt within 1,200 km. RMG = radius of the maximum gust.

| | Atlantic | Pacific |
|---|---|---|
| RMG, km | 251 [238, 263] | 222 [205, 231] |
| HF-equivalent area, 1e3 km2 | 13.1 [11.9, 14.4] (mean 32.6) | 10.9 [10.0, 12.0] (mean 26.6) |
| Core area (gust >= 90% of storm maximum), 1e3 km2 | 23.8 [22.8, 25.3] | 26.1 [24.6, 27.6] |
| HF-equivalent outer radius, km | 362 [341, 379] | 310 [294, 324] |
| Share of HF-eq area south of the centre | 77.5% [74.7, 80.0] | 82.0% [79.7, 84.2] |
| Share right of motion | 76.7% [74.5, 78.9] | 81.3% [79.4, 83.1] |
| Share rear-right quadrant | 54.1% [51.6, 56.5] | 61.1% [58.8, 63.6] |
| Share front-right / rear-left / front-left | 22.6 / 13.8 / 9.5% | 20.2 / 14.0 / 4.6% |
| Maximum gust in rear-right quadrant | 51.5% [49.6, 53.5] | 60.4% [58.3, 62.4] |

By stage (hours from the track's minimum MSLP: deepening <= -12 h, mature within
6 h, filling >= 12 h), Atlantic / Pacific:

| | deepening | mature | filling |
|---|---|---|---|
| fixes | 1,048 / 1,140 | 1,584 / 1,215 | 671 / 325 |
| RMG, km | 175 [158, 195] / 192 [176, 209] | 262 [251, 275] / 235 [223, 253] | 337 [306, 373] / 251 [208, 283] |
| HF-eq area, 1e3 km2 | 10.0 [8.4, 11.9] / 9.3 [7.9, 10.3] | 14.9 [13.3, 17.1] / 12.7 [11.2, 14.9] | 13.6 [11.0, 16.2] / 12.0 [8.5, 17.1] |
| Share of area south of centre | 81% / 83% | 76% / 82% | 77% / 79% |

## Outcome of each pre-registered hypothesis

- **H1, HF winds south of the centre: supported in both basins.** South share 77.5% and 82.0%, intervals exclude 50%; maximum gust south of centre at 65.5% and 75.6%. This is the crescent in the QuikSCAT composite (Von Ahn et al. 2006), reproduced in a sample about 120 times larger. The paper gives no radius, area or quadrant figure, so the comparison is of shape only. Because most storms move eastward or northeastward, "south" and "right of motion" overlap heavily; the two are not independent confirmations.
- **H2, rear-right quadrant: supported in both basins.** Rear-right holds 54% and 61% of the area against 25% for equal quadrants. The maximum gust sits there in 51.5% (Atlantic, interval just above 50%) and 60.4% (Pacific).
- **H3, life cycle: partly supported.** The area is smaller while deepening than at maturity in both basins (non-overlapping intervals). Mature versus filling is **not** separated in either basin (Atlantic 14.9 vs 13.6, Pacific 12.7 vs 12.0, intervals overlap), so the rule "mature largest" fails. RMG grows from deepening to filling: Atlantic 175 to 337 km, clearly; Pacific 192 to 251 km, intervals touching ([176, 209] and [208, 283]), so by the rule it is **not** supported in the Pacific. The concentration south of the centre eases from deepening to filling in the Atlantic (81% to 77% by area; 77% to 53% by position of the maximum).
- **H4, warm-core compact: contrary.** Thermally symmetric warm-core fixes have a **larger** RMG (256 km [249, 278]) and a larger HF-equivalent area (16.5 vs 6.4 x 1e3 km2) than asymmetric cold-core fixes (RMG 196 km [178, 209]). The opposite of the expectation of a compact seclusion core. The split is mainly warm versus cold core: asymmetric warm-core fixes are similar to symmetric warm-core ones (RMG 233 km, area 17.9e3 km2). Hart terms exist at 00/12 UTC only (2,921 of 5,983 fixes); the symmetric cold class is small (101 fixes).
- **H5, basin difference (no prior): Atlantic broader.** Atlantic RMG 29 km larger and area 2.2 x 1e3 km2 larger; Pacific more concentrated rear-right (61% vs 54%). Intervals on the RMG difference barely overlap; treat as suggestive. Part of the Atlantic excess is terrain (below).

## Caveats that bind the numbers

1. **ERA5 under-resolves the peak.** At 31 km ERA5 sustained 10 m wind reaches 64 kt at only 1.3% (Atlantic) and 0.6% (Pacific) of these HF-strength fixes; the median maximum sustained wind is about 50 kt, a gust factor of 1.5 (1,833 fixes with WeatherBench2 wind, one per storm). Only 45% / 38% reach 51.2 kt (64/1.25, the factor Jelenak et al. found). So an HF **area at 64 kt sustained is essentially zero in ERA5**, and every HF-equivalent area above is an area of gust >= 71.7 kt, calibrated to the storm maximum. Absolute areas are biased low relative to scatterometer or buoy winds; shapes and relative differences are more trustworthy than the km2 values. The HF-equivalent area is about 13e3 km2 (a disc of about 65 km radius), but the HF-equivalent *outer radius* is 360 km because the area is patchy; the core-area statistic (gust >= 90% of the storm's own maximum) does not depend on the absolute bias.
2. **Stage is defined here, provisionally.** Minimum MSLP over the in-domain track, not the life-cycle thread's definitions (not committed when this was run). Tracks truncated at the domain edge can mislabel stage.
3. **Fixes within a storm are correlated.** Intervals resample storms; the sample averages 3.0 fixes per storm and some storms contribute many.
4. **"Episode position" rows in `structure.txt` are a selection artefact.** The first and last HF fix of an episode are by construction the ones barely over the threshold, so they have small areas and low gust. They are shown because the plan said every comparison would be, not as life-cycle results.
5. **Terrain.** 468 Atlantic fixes (14% of Atlantic, 7.8% of all) have their maximum gust within 100 km of Greenland or Iceland (flag, not studied; owned by the high-latitude thread). Without them Atlantic RMG falls from 251 to 231 km and the south share rises from 77.5% to 81.7%, i.e. the Atlantic broad, less-southerly picture is partly orographic.
6. **Pre-2004 and drift.** Seasons 2004-05 on only, where the gust index does not drift. Nothing here bears on the pre-2001 drift question.
7. **Transitioning tropical cyclones are in** (pipeline A includes them, 5-7% of events).
8. **Pipeline A fixes, not archive fixes.** A fix here is an ERA5 track point whose gust index passed the calibrated threshold; the archive's own HF fix is OPC's judgement.

## Reproduce

`ERA5_WORK=<dir> python3 extract.py 8` (about 10 minutes at 8 processes, 33 GB) then
`ERA5_WORK=<dir> python3 composite.py <env_2004.csv.gz from research/era5/intensity/results>`.
Needs numpy, scipy, pandas, numcodecs, matplotlib and network access to the public
ARCO-ERA5 and WeatherBench2 buckets. `fixes.csv` alone reproduces the tables except
the composites (which need the boxes).
