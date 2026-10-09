# Jet-streak quadrant and upstream-trough terms: do they add skill to the HF model?

ERA5 **proxy**, **pipeline A** tracks (`../hf_history`), outcome = the 800 km ocean gust index reaching 71.7 kt. WeatherBench2 1.5 degree fields, October to April, seasons 2004-05 to 2021-22 (18; 1 June to 31 May), 74,009 fixes on 16,874 tracks, 4,012 reaching HF within 24 h. Leave-one-season-out. The plan was committed before any field was read: `PREREGISTRATION.md` (commit 78c7a3b). Numbers below are from `results/skill.txt`, `extras.txt`, `qc.txt`, `descriptive.txt`.

## Answer in plain words
1. **Jet-streak quadrant (is the low in the right entrance or left exit?): no, for the chance of hurricane force within 24 hours.** Adding the jet-streak terms (12 at the fix and 12 at 24 h earlier, 26 columns with the missing-lag flag and the 24 h change) to the PR 12 model changes the Brier skill score by +0.0016 [-0.0007, +0.0036] (base 0.413). The upper end is below the 0.005 that the plan called worth having, so this is a well-powered null: effects of half a point of skill or more are excluded.
2. **Upstream trough: no useful gain for HF within 24 hours.** A small positive +0.0015 [+0.0003, +0.0028], below the 0.005 worth having and not significant after the multiple-test correction (q 0.19); same well-powered-null verdict. Both together: +0.0029 [-0.0002, +0.0056]; the upper end is just above the 0.005 threshold, so "can't tell"; any gain is at most about half a point or slightly more.
3. **The trough does help predict rapid deepening, by a small amount.** Trough terms add +0.0080 [+0.0056, +0.0104] to the rapid-deepening score (0.300 to 0.308), in 17 of 18 seasons and in both basins (Atlantic +0.0071, Pacific +0.0086). As a yes/no forecast the gain is tiny (HSS 0.474 to 0.481). Jet-streak terms add nothing there (+0.0007).
4. **HF within 48 hours: borderline.** Jet plus trough add +0.0051 [+0.0023, +0.0078] (q 0.016 in its family, 0.044 across all 19 tests); it clears the 0.005 threshold only marginally, and the yes/no change is nil (HSS +0.002).
5. **Descriptively, HF lows sit more often in the ascent quadrants.** At HF onset the low is in the right-entrance or left-exit quadrant in 41% of Atlantic cases against 30% of matched storm-force-only peaks, and 42% against 26% in the Pacific (sum of the RE and LX rows of `descriptive.csv`). The streak is about 4 kt stronger and the trough 14 m (Atlantic, q 0.054) and 26 m (Pacific, q 0.002) deeper. That is not a skill result: HF storms are stronger by definition, and the model already knows the storm's strength.

Likely reason the HF result is flat: most of the model's skill comes from the storm's own current gust and pressure, and the quadrant and trough terms mostly repeat what those already say. This explanation was not tested here.

## What was done
- Pull: 1,946 WeatherBench2 chunks, 55.50 GB (limit 60 GB; approved about 55 GB), plus 12 chunks (about 0.34 GB, estimated) for the QC drawings. Nothing raw is kept; `results/features.csv.gz` (7.7 MB) holds the extracted features (75,000 fixes, including 28-30 September lag sources).
- Features (`features.py`): jet streak in the vortex-removed 250 hPa wind (strength, half-length, quadrant of the low, along-flow speed gradient, 300 hPa divergence) and the upstream trough in Z500 anomaly (depth, distance, bearing, amplitude, tilt, phase). Definitions are in the plan; lags 24 h (primary), 0, and 0/12/24/48.
- Tests: base + group, logistic, 18 folds, season-block intervals, sign-flip p, Benjamini-Hochberg within each pre-registered family. Primary family (3 tests): J2, T2, J2+T2 on HF within 24 h with lags 0 and 24 h. 16 secondary tests in 5 further families. `results/gains.png` shows all.
- **Primary family result: 0 of 3 meet the rule** (gain >= 0.005, interval above 0, q < 0.05). Across the 19 tests, 3 meet the rule: trough and jet+trough for rapid deepening, jet+trough for HF within 48 h. Across all 19 with one BH, q < 0.05 in 4 and q < 0.10 in 8.

## Detection checks (made before any outcome was read)
- Q1 passes: a streak is found for 99.97% of fixes (21 of 74,008 miss) with a 60 kt wind maximum within 2,500 km (this threshold is easy to meet in winter; the check is weak), and a trough for 76% of all fixes (84% of those with a heading; 9% have no heading).
- Q2 passes: mean 300 hPa divergence is higher in the right-entrance and left-exit quadrants than in the other two, by +0.18 (Atlantic) and +0.20 (Pacific) x 1e-5 /s, q 0.0005 in both basins. So the quadrants are physically real in these data.
- Q3 (12 drawn cases, `results/qc_cases.png`; Storm Dennis, ex-Typhoon Nuri and 10 random): the streak outline, axis and quadrant look right in all 12. **The trough finder often picks the western side of the storm's own upper low rather than a separate upstream trough** (Dennis, and random cases 5 and 7 have the minimum within about 900 km of the low); Nuri has no trough in the sector. Logged, not tuned. This is the main weakness of the trough terms.

## Sensitivity (pre-registered; 15 seasons 2004-05..2018-19 only, see deviations)
HF within 24 h: no variant changes the null (jet with vortex removal 1000 km, +0.0009; none, +0.0010; trough with exclusion 500 km, +0.0008; 1000 km, +0.0005). Rapid deepening: trough gain +0.0078 with the default 750 km exclusion, +0.0088 at 500 km and +0.0050 [+0.0023, +0.0078] at 1000 km, so the gain falls as the exclusion grows (the difference was not tested), consistent with part of it coming from the storm's own upper low; a positive gain remains with a trough at least 1000 km away.

## Limits
Perfect-prognosis test: ERA5 at the fix and earlier, not a forecast test. 1.5 degree grid smooths streak cores. The storm's own circulation is removed from the wind only as an azimuthal mean and not from Z500. Oct-Apr only. Nothing here uses pre-2004 seasons or 2022-23 on. Pipeline A only; the ERA5 record is a proxy.

## Deviations (post hoc)
- The sensitivity runs cover 15 of 18 seasons. The first 304 chunks of the pull (newest seasons, 8.7 GB) were processed before the variants were added; repeating them would have passed the 60 GB limit. The base was refitted on the same 15 seasons.
- Q3: only two of the four named storms (Dennis, Nuri) fall inside the sample (Braer 1993 and Fiona 2022 do not); ten random cases were drawn instead of eight (twelve in total).
- Q1 trough criterion was applied to all fixes (fixes with no heading count as not found); both versions are shown.
- The test-count line of `skill.txt` was corrected after the run (a script bug counted 0 q(all) passes); the table values were unaffected.

## Verification
Fresh-agent recomputation and an independent re-implementation of the features: see `VERIFICATION.md`.

## Reproduce
`ERA5_WORK=DIR python3 extract.py 4` (55 GB, resumable), `ERA5_WORK=DIR python3 qc.py`, `python3 qc_cases.py`, `SCRATCH=DIR python3 skill.py results`, `python3 desc.py results`, `SCRATCH=DIR python3 extras.py results`.
