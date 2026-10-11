# Pre-registration: ERA5 vertical motion (omega) and ascent coupling as predictors of P(HF)

Status: written before any omega field was read. ERA5 proxy, pipeline A (`research/era5/hf_history`); outcome = pipeline A gust index reaching 71.7 kt. Follows PR 149 (`hf_jet_trough`, same fixes, same model, same rule). Requested by Jason 2026-10-09 18:38Z ("is the area big enough, could an HF storm be coupled to jets, is there a better measure of vertical lift"); go-ahead by card tap, 2026-10-11 00:10Z, for "omega + coupling".

## Question
Does ERA5 vertical motion around the low, and whether it has more than one ascent centre, add skill to the PR 12 near-storm model, with and without the jet and trough features of PR 149?

## What has been looked at already
- PR 12 base (`fixes_2004`, `env_2004`) has only 500 km means of 300 hPa divergence and 500 hPa vorticity advection at the fix. No vertical motion is in the base.
- PR 149: jet-streak (J2) and trough (T2) features, box +-4000 km, one streak within 2,500 km, one trough in a 750-3,000 km rear sector. Base BSS 0.4133 for hf24 on 74,009 fixes, 18 seasons; J2 +0.0016, T2 +0.0015, J2+T2 +0.0029, none passing the rule. Rapid deepening gained +0.0080 from T2.
- Nothing about omega has been computed on these fixes. A cost sizing was done (HEAD requests only).

## Deviation from the card text (stated before anything is pulled)
The card said "omega plus a wider box and a flag for a second streak, about 22 GB". The wind-based second-streak flag needs u and v at 250 hPa again (20.9 MB per chunk, +40.6 GB), so omega plus winds is 62.5 GB, over the gate. Only the omega part is within what was approved. Here coupling is defined from the ascent field itself. A jet-wind second-streak flag is not part of this test and would need its own go-ahead.

## Sample
Same fixes as PR 149: Oct-Apr (plus 28-30 Sep as lag sources), seasons 2004-05 to 2021-22 (18), 00/12 UTC, transitioning tropical cyclones in. Leave-one-season-out. Fit and test 2004-05 on only (decision 1).

## Field and frame
WeatherBench2 1.5 deg, 6-hourly, `vertical_velocity` (Pa/s, negative = ascent), levels 500 and 700 hPa. Define ascent A = -omega (positive = rising). Box +-5000 km about the low at 100 km spacing (101 x 101), bilinear, great-circle placement as in `hf_jet_trough/features.py`; 3x3 binomial smoothing. Motion frame: +x along the heading (previous 6 h); fixes without a heading use the north frame and the motion-relative features are missing with an indicator, as in PR 149.

## Predictors
Group W (strength and position of lift), at 500 and 700 hPa:
1. `a500_1000`, `a700_1000`: maximum A within 1,000 km.
2. `a500_500`, `a700_500`: mean A within 500 km.
3. `aarea500`: area (1e6 km2) with A500 >= 0.3 Pa/s within 3,000 km.
4. `adx_km`, `ady_km`: position of the A500 maximum (within 2,500 km) relative to the low, along and across the heading (km/1000).
5. `afwd500`: mean A500 in the forward half (+-90 deg of heading) minus the rear half, within 2,000 km.
Group C (coupling), from A500 smoothed:
6. `nasc`: number of distinct ascent centres within 4,000 km: local maxima (8-neighbour) with A500 >= 0.4 Pa/s, at least 1,000 km apart (keep the stronger of two closer ones).
7. `asc2`: A500 of the second strongest centre (0 if none); `asc2dist` its distance (km/1000, median-filled if none); `asc2cos`, `asc2sin` its bearing relative to heading.
8. `couple_w` = `nasc` >= 2 and `asc2` >= 0.4 Pa/s and `asc2dist` <= 2.5 (one pre-registered indicator).
Lags: all W and C values at 0 and 24 h by lookup of the track's earlier fix (as PR 149); 24 h change in `a500_1000`, `asc2`. Continuous values standardised per basin and calendar month from predictors only.
Storm-separation sensitivity (logged, not used to select): W recomputed with the ring 750-2,500 km only (`a500_x750`), because the storm's own latent-heating ascent is in the 0-750 km core.

## Detection QC (before any outcome is joined)
Run with outcome columns not loaded.
- Q1 an A500 maximum >= 0.2 Pa/s exists within 2,500 km for >= 80% of fixes; `nasc` >= 2 for between 10% and 70% of fixes (a flag that is nearly always or never on cannot be used).
- Q2 physical sign: using PR 149's quadrants (features.csv.gz), mean `a500_1000` is higher in RE and LX than in LE and RX, in both basins, season-block bootstrap q < 0.05. A check that the lift field and the jet quadrants agree, not a skill test.
- Q3 12 cases (the 4 named storms of the probability README plus 8 random) drawn with the ascent centres marked, reviewed by eye; failures logged, not tuned away.
If Q1 or Q2 fails the features are not used for the skill test and the failure is the result. No threshold is changed after a QC result; any change is a logged deviation with Q1-Q3 rerun.

## Skill tests
Logistic, L2, C = 1, median impute, clip 0.5/99.5, standardised, the PR 149 code (`skill.py`, `../intensity_extra/evaluate.py`). Base = PR 12 full on the same fixes. Gain = LOSO BSS minus base BSS, 90% season-block interval, exact sign-flip p over 18 seasons, BH per family. Rule: gain >= 0.005, interval above 0, q < 0.05 = "adds skill"; upper interval end < 0.005 = well-powered null.
- Primary family (BH over 3): W, C, W+C added to the base; target hf24; lags 0 and 24 h.
- Secondary families (BH within each): the same three for rapid deepening (ndr24 >= 1) and hf48; the same three added on top of base + J2 + T2 (incremental to the jet and trough features); lag sets {0} and {0,12,24,48}; `couple_w` alone; W on the x750 ring.
- Season-level gains for each of 18 seasons; by basin descriptive.
All tests are reported and none is promoted afterwards. Total count and the count passing the rule out of the total are stated.

## Descriptive companion
Mean A500 within 1,000 km and `nasc` by PR 149 jet quadrant, and at HF onset against the storm-force-only peak (same anchors as `hf_vs_storm`, matched on basin and month; season-block bootstrap, BH). The ascent area and share with two centres at HF onset, by basin.

## Pull size (measured)
Compressed omega per two-day chunk, mean of 12 random chunks: 11.28 MB (all 13 levels in a chunk). The PR 149 set is 1,946 chunks, so 1,946 x 11.28 MB = 21.9 GB. Under the 50 GB gate. Resumable, per chunk, only features kept, nothing raw committed, stop if the running total passes 24 GB. Every sensitivity variant is computed in the same pass.

## Not done and limits
Not a forecast test (perfect prognosis from ERA5). 1.5 deg omega is smooth: narrow frontal ascent is under-resolved. Omega contains the storm's own ascent and its latent heating, so it is partly an outcome-side diagnostic; the x750 ring and the lags are the guards. No Q-vector, no second-streak detection from wind (needs +40.6 GB), no pull after 2023-01-09, pipeline A only. Pre-registered power: the PR 149 base left an interval of about +-0.002 BSS around zero gains; the same data are used, so an effect of 0.5 BSS point or more is detectable.

## Deviations (post hoc)
None yet.
