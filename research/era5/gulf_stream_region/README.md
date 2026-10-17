# Gulf Stream and marine wind, region-wide: satellite winds and GFS

Plan: `PREREGISTRATION.md` (committed before any ASCAT wind, ERA5 temperature or GFS field was read; eight deviations logged there; items 1-7 before any outcome, item 8 (the registered 5-of-5 interval rule and the B version of S5, added as scripts) after the primary results). Follows the 7-buoy study in `../gulf_stream/` (PR 263), which Jason judged too few buoys to be a test.

Data: ASCAT scatterometer winds (Metop-A, -B, -C), 1 degree box means per satellite pass, western Atlantic 29-41N 80-55W, October to April. ERA5 enters only as a control (pressure-gradient wind, air and sea temperature); it is a proxy and is never the outcome. Part B compares against 24 h GFS forecasts. No ERA5 gust, no hurricane-force definition, neither gust pipeline.

## What a forecaster can use

**Over the open Gulf Stream region, wind is much higher than the pressure gradient suggests when the water is 6 C or more warmer than the air: about 3 to 4 kt more in the raw numbers. GFS at 24 h already carries most of that; the part it misses is small, about 0.5 to 1 kt, and a little larger right at the Stream's north wall.**

Raw medians (ASCAT wind, pressure-gradient wind 25 to 35 kt, nothing fitted):

| Water minus air | All boxes | At the Stream's north wall | 2 degrees or more away |
|---|---|---|---|
| colder than air (below 0 C) | 19.3 kt | 18.9 | 19.3 |
| 0 to 3 C (reference) | 21.6 | 22.0 | 21.5 |
| 3 to 6 C | 23.6 | 22.6 | 23.8 |
| 6 to 9 C | 25.0 | 23.4 | 25.7 |
| 9 C or more | 25.6 | 24.0 | 26.8 |

Fitted (12 winters, pressure-gradient wind held fixed, box effects removed), wind when water is 6 C or more warmer than air versus 0 to 3 C warmer:

| Test | Change | Knots (reference wind 21.6 kt at gradient 30 kt) | Plan verdict |
|---|---|---|---|
| A1 all boxes | +15.9% (+13.9 to +17.9) | +3.4 (+3.0 to +3.9) | yes (q 0.0002) |
| A2 at the north wall (within 1 degree latitude of the front) | +8.2% (+6.2 to +10.1) | +1.8 (+1.3 to +2.2) | yes |
| A3 2 or more degrees away from the front | +19.1% (+17.3 to +21.1) | +4.1 (+3.7 to +4.6) | yes |
| A4 north wall minus away | -9.2% (-11.4 to -7.0) | -2.0 (-2.5 to -1.5) | **can't tell: significant, opposite sign** |
| A5 north wall, wind from 270-360 degrees (the buoy lead) | +9.8% (+8.4 to +11.3) | +2.1 (+1.8 to +2.5) | yes (confirms the buoy lead on new data, at about half the buoy size) |

- **The bump is biggest away from the Stream front, not at it.** That is the opposite of the "stronger over the Stream" idea, so by the plan's rule A4 is "can't tell". Away from the front, 6 C or more of air-sea difference mostly means cold air from the continent over warm water, and the pressure-gradient wind (geostrophic, from pressure alone) is a poorer yardstick in that air than elsewhere. The data cannot separate instability from the other differences between those two groups.
- **A shuffled-temperature check (placebo) gives +2.5% for A1** (range +2.3 to +2.6, 20 shuffles), more than the plan's 1 point, so about 2.5 points of the A1 figure comes from the method. The plan's real-minus-offset reading is about +13% (about +2.8 kt). The placebo was run for A1 only.
- **Equivalent-neutral caveat.** ASCAT reports the wind that would give the same stress in neutral air, so in unstable air it reads high against the true 10 m wind. The approximate correction (S5) takes A1 from +15.9% to +15.1%. For part B it takes B1 from +2.2% to +1.4% and B2 from +3.6% to +2.7% (B3 +0.6%, interval spans 0). The B numbers below are therefore upper limits for the true 10 m wind; for waves, which respond to stress, they are closer to right.
- Robustness (secondary, raw p): A1 holds in Dec-Mar only (+14.9%), in each satellite separately (+14.3 to +17.4%), at gradient wind 30 kt or more (+15.2%), and over 2021-25 alone (+14.0%). Leaving out one winter at a time moves A1 from +15.1 to +16.2%, A2 from +7.4 to +8.8%, A3 from +18.4 to +19.5%.

### How much to add to GFS (24 h forecast)

GFS runs about 3.6% below ASCAT even in the reference air-sea class (ASCAT/GFS ratio 1.036). That is a general model-versus-satellite difference, and the same in every class, so it is not a Stream effect. The extra when the model's own water-minus-air difference is 6 C or more:

| Test (5 winters, 2021-22 to 2025-26) | Extra wind | Knots at model wind 25 kt (35 kt) | Plan verdict |
|---|---|---|---|
| B1 all boxes | +2.2% (+0.9 to +3.5) | +0.6 (+0.2 to +0.9) (+0.8) | yes (q 0.0008; interval above 0 in 5 of 5 leave-one-winter-out refits) |
| B2 north wall | +3.6% (+1.6 to +5.2) | +0.9 (+0.4 to +1.3) (+1.3) | yes (q 0.0008; 5 of 5 refits) |
| B3 away from the front | +1.4% (0.0 to +2.5) | +0.3 (0.0 to +0.6) | can't tell (q 0.060; only 2 of 5 refits exclude 0) |
| B4 north wall minus away | +2.2% (-0.1 to +4.4) | +0.6 (0.0 to +1.1) | can't tell (q 0.064) |

How to read it: at 24 h, GFS already includes most of the Stream's extra wind. Adding about 1 kt near the north wall, less elsewhere, is what these five winters support. There is no support here for adding several knots to GFS. Water colder than air is where the model errs the other way: GFS runs about 7% high against ASCAT there (ratio 0.966 versus 1.036, descriptive, not tested). With wind from 270-360 degrees at the north wall the extra is +6.3% (B5, secondary, raw p). B6 (a 48 h forecast) was not run because it needed B1 to pass first and the extra pull was not made.

Part B uses 5 winters, so its intervals are wide by construction. The 5-of-5 rule from the plan is met for B1 and B2 (`lowo_B.py`, `results/lowo-B-result.txt`).

## Power (`power-result.txt`)

Planted-effect simulation, 100 outer by 200 inner winter resamples (the plan said 200 by 400; reduced for run time, logged below). A planted +3% where the water is 6 C or more warmer is found 97% (A1), 92% (A2), 100% (B1) and 90% (B2) of the time (unadjusted p < 0.05, correct sign); the difference tests A4 and B4 find a planted +3% only in the band 80% and 86% of the time, and correctly stay quiet (1 to 2%) when the same +3% is planted everywhere. The plan's "well-powered no" rule needed a 95% interval inside +-1 kt at 25 kt; no test here is a null, so it was not needed. The simulation keeps day-to-day dependence only within winters, so it is an upper bound.

## Limits

- ASCAT: two passes a day, footprints about 25 km, rain-flagged cells dropped, no data above about 45 kt (saturates near 50), and equivalent-neutral (above). The north wall is a few tens of km wide; its position here comes from daily 00 UTC ERA5 SST at 0.25 degree.
- The pressure-gradient yardstick in part A is imperfect. The buoy study found the wind-to-gradient ratio falls from 0.84 to 0.46 as gradient wind rises. The model is log-linear, so some stability can leak into the temperature effect, and the shuffle check shows about 2.5 points of the A1 estimate comes from the method.
- Longitudes east of 65W have no defined front (the north wall is outside the 32-42N search window), so A2 to A5 and B2 to B4 use the western part of the region only.
- Waves were not studied here (Jason chose A and B only). The buoy study's wave result stays "can't tell". Wind against current is still untested with real currents.
- GFS is one model at one lead. The marine wind (GFE blend) project's own results are separate and not repeated here.

## Verification

Filled in after the fresh-agent recompute; see the section below.

## Reproduce

```
python3 research/era5/gulf_stream_region/pull_ascat.py /mnt/project-files/gulf-stream/ascat
python3 research/era5/gulf_stream_region/pull_era5_temp.py /mnt/project-files/gulf-stream/era5temp
python3 research/era5/gulf_stream_region/pull_gfs.py /mnt/project-files/gulf-stream/gfs
python3 research/era5/gulf_stream_region/build_boxes.py        # writes boxes_A.csv.gz, boxes_B.csv.gz (59 MB, 28 MB; not in git, in /mnt/project-files/gulf-stream/tables)
python3 research/era5/gulf_stream_region/power_region.py work/boxes_A.csv.gz work/boxes_B.csv.gz power-result.txt
python3 research/era5/gulf_stream_region/analyse_region.py work/boxes_A.csv.gz work/boxes_B.csv.gz results
python3 research/era5/gulf_stream_region/lowo_B.py work/boxes_B.csv.gz results/lowo-B-result.txt
python3 research/era5/gulf_stream_region/s5_B.py work/boxes_B.csv.gz results/s5-B-result.txt
```
