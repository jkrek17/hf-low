# Gulf Stream and marine wind, region-wide: satellite winds and GFS

Plan: `PREREGISTRATION.md` (committed first, `402bec7`, before any satellite wind, ERA5 temperature or GFS field was read; nine deviations logged there: items 1-4, 6 and 7 before any outcome, item 5 after one indirect look at the discarded first power run, items 8 and 9 after the primary results). Follows the 7-buoy study in `../gulf_stream/` (PR 263), which Jason judged too few buoys to be a test.

Data: ASCAT satellite winds (Metop-A, -B, -C), averaged over 1-degree boxes for each satellite pass, western Atlantic 29-41N 80-55W, October to April. ERA5 is a control only (pressure-gradient wind, air and sea temperature) and is a proxy, never the outcome. Part B compares with 24 h GFS forecasts. No ERA5 gust, no hurricane-force definition, neither gust pipeline. ASCAT has no data above about 45 kt, so **nothing here applies at storm or hurricane force.**

## What a forecaster can use

**In the satellite winds, when the water is 6 C or more warmer than the air, the wind is about 3 kt higher than when it is only 0 to 3 C warmer, at the same pressure gradient. GFS at 24 h already shows most of that. What it misses is under 1 kt in the satellite's own wind, and it disappears when the satellite wind is converted to the true 10 m wind. These data do not support adding wind to GFS.**

Raw medians (satellite wind, pressure-gradient wind 25 to 35 kt, nothing fitted):

| Water minus air | All boxes | Within 1 degree of latitude of the Stream front | 2 degrees or more from the front |
|---|---|---|---|
| colder than air (below 0 C) | 19.3 kt | 18.9 | 19.3 |
| 0 to 3 C (reference) | 21.6 | 22.0 | 21.5 |
| 3 to 6 C | 23.6 | 22.6 | 23.8 |
| 6 to 9 C | 25.0 | 23.4 | 25.7 |
| 9 C or more | 25.6 | 24.0 | 26.8 |

The fitted comparison below compares each 1-degree box with itself, holds the pressure-gradient wind fixed, and contrasts 6 C or more (average of the 6-9 and 9+ classes) with 0 to 3 C. 12 winters, 2013-14 to 2024-25. "q" is the probability value after allowing for running five tests at once; below 0.05 counts as a real difference.

| Test | Wind change | Knots (reference 21.6 kt, gradient wind 30 kt) | Plan verdict |
|---|---|---|---|
| A1 all boxes | +15.9% (+13.9 to +17.9) | +3.4 (+3.0 to +3.9) | yes (q 0.0002) |
| A1 after the shuffle check, see below | **+13.1%** | **+2.8** | yes |
| A1 after converting to true 10 m wind (S5, approximate) | +11.3% (+9.3 to +13.4) | +2.4 | secondary |
| A2 within 1 degree of the front | +8.2% (+6.2 to +10.1) | +1.8 (+1.3 to +2.2) | yes |
| A3 2 or more degrees from the front | +19.1% (+17.3 to +21.1) | +4.1 (+3.7 to +4.6) | yes |
| A4 near the front minus away | -9.2% (-11.4 to -7.0) | -2.0 (-2.5 to -1.5) | **can't tell** (significant, but the plan expected the opposite sign, from the buoy study) |
| A5 near the front, wind from 270-360 degrees | +9.8% (+8.4 to +11.3) | +2.1 (+1.8 to +2.5) | yes (positive, as registered) |

- **Shuffle check.** If the temperature differences are shuffled among boxes of the same winter, month and pressure-gradient third, the method still reports +2.5% (range +2.3 to +2.6, 20 shuffles), so about 2.5 points of A1 come from the method, not the temperature. The plan's rule is to give the estimate with that offset removed: +13.1%. It was run for A1 and B1 only; for A2, A3, A4 and A5 the offset is unknown. For B1 the shuffle gives +0.5%, under the plan's 1 point, so no correction (corrected would be +1.7%). The shuffle that keeps the whole winter together gives about 0, so the offset comes from differences between months and gradient classes.
- **Where the bump is bigger.** In this fit it is larger away from the front (+19%) than near it (+8%), and the plan reads a difference in that direction as can't tell. Almost all of the "away" boxes lie on the warm, south side of the Stream (97.5% of them, and 95% of those with water 6 C or more warmer), so "away" means the warm side, not both sides. The front position comes from a daily 0.25-degree sea-surface-temperature map, and a box 1 degree wide can be given the wrong temperature difference where the water changes sharply, which would shrink the near-front estimate. The two groups also differ in air mass and in how well the pressure gradient predicts the wind. The data cannot say why they differ. A5 (northwest wind) is +9.8%, close to the +8.2% for all wind directions near the front, so northwest flow does not stand out here. The buoy study's northwest lead (+19.7%, post hoc) is not reproduced at that size. The satellite is an independent instrument, but its winters overlap the buoy study's.
- **The satellite reads high in unstable air.** ASCAT reports the "equivalent-neutral" wind, the wind that would give the same drag on the sea surface in neutral air. In unstable air that is above the true 10 m wind. Converting (S5, a standard approximate stability correction, `region_lib.convert_en`) lowers A1 from +15.9% to +11.3%. The first version of this conversion was wrong by roughly a factor of ten (it moved A1 by only 0.8 points); a reviewer caught it and it was corrected after the primary results (deviation 9, `test_convert_en.py`).
- **Checks on A1 (secondary, labelled as such):** Dec-Mar only +14.9%; each satellite on its own +14.3 to +17.4%; gradient wind 30 kt or more +15.2%; 2021-25 alone +14.0% (post hoc, not in the plan). Leaving out one winter at a time: A1 +15.1 to +16.2%, A2 +7.4 to +8.8%, A3 +18.4 to +19.5%. The plan's continuous version (S1, effect per +5 C) was not run.

### How much to add to GFS (24 h forecast)

In the 0-3 C class ASCAT reads 3.6% above GFS (ratio 1.036). Some or all of that may be the satellite: on 348 buoy matches ASCAT read 6% above the buoys (deviation 1). It is not a reason to add wind to GFS. The extra that the satellite shows over GFS when the model's own water-minus-air difference is 6 C or more:

| Test (5 winters, 2021-22 to 2025-26) | Satellite-minus-GFS gap grows by | Knots at model wind 25 kt (35 kt) | After converting to true 10 m wind (S5) | Plan verdict |
|---|---|---|---|---|
| B1 all boxes | +2.2% (+0.9 to +3.5) | +0.6 (+0.2 to +0.9) (+0.8) | -1.9% (-3.2 to -0.8) | yes (q 0.0008; interval above 0 in 5 of 5 refits that leave one winter out) |
| B2 within 1 degree of the front | +3.6% (+1.6 to +5.2) | +0.9 (+0.4 to +1.3) (+1.3) | -0.2% (-2.2 to +1.5) | yes (q 0.0008; 5 of 5 refits) |
| B3 2 or more degrees away | +1.4% (0.0 to +2.5) | +0.3 (0.0 to +0.6) | -2.4% (-3.7 to -1.2) | can't tell (q 0.060; only 2 of 5 refits) |
| B4 near the front minus away | +2.2% (-0.1 to +4.4) | +0.6 (-0.0 to +1.1) | not run | can't tell (q 0.064) |

How to read it: GFS's own winds rise with the water-minus-air difference nearly as much as the satellite's do (an inference from comparing two separate fits, not a direct test). In the satellite's own terms the remaining gap is 0.5 kt over all boxes to about 1 kt near the front, and once the satellite wind is converted to true 10 m wind the gap is zero or slightly negative. Whether the front needs more than elsewhere: can't tell. When the water is colder than the air the satellite reads about 3.5% below GFS (7% lower than in the reference class; descriptive). A stable-air satellite wind is expected to read low, so this is not shown to be a model error. With wind from 270-360 degrees near the front the gap is +6.3% (B5, secondary). The 48 h forecast test (B6) was eligible because B1 passed, and was not run (left unfinished).

Part B has 5 winters, so its intervals are wide. The plan requires the interval to exclude 0 in 5 of 5 refits that leave one winter out; B1 and B2 meet it (`lowo_B.py`). The conversion (S5) versions were not checked that way.

## Power (`power-result.txt`)

Planted-effect simulation (add a known +3% where the water is 6 C or more warmer, then see how often the test finds it): found 97% of the time for A1, 92% for A2, 100% for B1 and 90% for B2. A planted +3% only near the front is found 80% (A4) and 86% (B4) of the time, and the difference tests stay quiet (1 to 2%) when the same +3% is planted everywhere, as they should. This keeps day-to-day weather dependence only within winters, so it flatters the real tests somewhat. B3 and B4 did not pass. The plan's rule for a well-powered "no" (a +4% plant, interval inside +-1 kt) was not run, so they stay can't tell. The simulation was run at 100 by 200 resamples instead of 200 by 400 (deviation 7); a first run with a wrong design is kept in `power-result-first-run-discarded.txt`.

## Limits

- ASCAT: two passes a day, footprints about 25 km, rain-flagged cells dropped (rain probably differs with the temperature difference), no data above about 45 kt, and equivalent-neutral (above). The north wall of the Stream is a few tens of km wide.
- The pressure-gradient wind from pressure alone is an imperfect yardstick. In the buoy study the wind-to-gradient ratio fell from 0.84 to 0.46 as gradient wind rose (post hoc check D5). The fit is a straight line in the logarithm of the gradient wind, so part of the temperature effect can be the yardstick missing something else; the shuffle check measures part of that.
- On some days the front cannot be located, mostly east of 65W. Those boxes drop out of the front tests, but about a third of the near-front boxes (30%) and two fifths of the away boxes (43%) are still east of 65W.
- Waves were not studied here (Jason chose A and B only). The buoy study's wave result stays can't tell. Wind against current is still untested with real currents.
- GFS is one model at one lead. The marine wind (GFE blend) project's own results are separate and not repeated here.

## Verification

Two fresh agents recomputed numbers from the box tables with their own code. **Matched** (within 0.2 points, medians within 0.1 kt): both sample sizes and band counts; all 15 raw medians; A1 to A5 and B1 to B4 estimates; the knots for A1, B1 and B2; the leave-one-winter-out ranges for A1-A3 and B1-B3; A1 and B1 winter-resample intervals (+13.95 to +17.83 and +0.92 to +3.47 against +13.9 to +17.9 and +0.9 to +3.5); the S5 refits for A1, B1, B2, B3 (+11.28, -1.88, -0.21, -2.37, exact); the placebo means (A1 +2.42, B1 +0.53, winter-only shuffle +0.13); the front-position shares (all 6 exact). The conversion formula was judged a sensible approximation (single fixed surface roughness, dry-air Richardson number). The first verifier's placebo disagreed (+0.06) because its shuffle wrote permuted values to the wrong rows; the second verifier, told how to write them back, matched.
**Not independently recomputed:** the power simulation; the A4 and B4 interval ends and the q values (the reviewer checked they follow from the listed p); the secondary S2, S3, S4, S6 and B5 figures; the leave-one-winter-out intervals for B (`lowo_B.py`, own code, one run); the 20-shuffle placebo ranges; the S5 intervals; the raw A1-B1 knots at 35 kt. **Not independent:** the choice of bins, the 1-degree front band, the 25 and 35 kt reference winds, and the 6 C cut were set in the plan. A reviewer checked every claim in the first draft against the result files; the fixes are in this version.

## Reproduce

```
python3 research/era5/gulf_stream_region/pull_ascat.py /mnt/project-files/gulf-stream/ascat
python3 research/era5/gulf_stream_region/pull_era5_temp.py /mnt/project-files/gulf-stream/era5temp
python3 research/era5/gulf_stream_region/pull_gfs.py /mnt/project-files/gulf-stream/gfs
python3 research/era5/gulf_stream_region/build_boxes.py   # takes SHARED_DIR OUT_DIR; writes boxes_A.csv.gz, boxes_B.csv.gz (59 MB, 28 MB; not in git; copies in /mnt/project-files/gulf-stream/tables)
python3 research/era5/gulf_stream_region/power_region.py boxes_A.csv.gz boxes_B.csv.gz power-result.txt
python3 research/era5/gulf_stream_region/analyse_region.py boxes_A.csv.gz boxes_B.csv.gz results
python3 research/era5/gulf_stream_region/lowo_B.py boxes_B.csv.gz results/lowo-B-result.txt
python3 research/era5/gulf_stream_region/s5_conversion.py boxes_B.csv.gz results/s5-result.txt boxes_A.csv.gz
python3 research/era5/gulf_stream_region/test_convert_en.py
```

(`build_boxes.py` reads the land-sea mask from `/mnt/project-files/era5-cache/lsm_nh_0p25.npy`.)
