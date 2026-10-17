# Gulf Stream effects on marine wind and waves

Plan: `PREREGISTRATION.md` (committed first, `d648513`; two deviations and one count correction logged there, all before any outcome). Data: NDBC buoy reports (observations) at 00 and 12 UTC, October to April, 21 winters 2004-05 to 2024-25. ERA5 enters only as a pressure-gradient control (a proxy, not an outcome; ERA5 assimilates buoy pressures). No ERA5 gust, no hurricane-force definition, neither gust pipeline.

## What a forecaster can use

**Wind: there is no single "add N knots for instability" number in these buoys. At the three Stream buoys the data show a bump of about 1 kt; at the shelf buoys they do not, and the split may reflect where the buoys sit, not the Stream.**

Compared with the same pressure-gradient wind and water 0 to 3 °C warmer than the air, wind when the water is 6 °C or more warmer than the air was:

| Buoys | Change in wind | Knots at pressure-gradient wind 30 kt (reference wind 20.6 kt) | Pre-registered verdict |
|---|---|---|---|
| All 7 pooled (W1) | -0.2% (-2.3 to +2.2) | 0.0 (-0.5 to +0.5) | no pooled effect, but it averages opposite-signed groups |
| 3 Stream buoys (W2: 41002, 41001, 41025) | +4.7% (+1.2 to +8.7) | +1.0 (+0.2 to +1.8) | yes |
| 4 shelf buoys (W3: 41013, 44014, 41004, 44009) | -2.8% (-4.5 to -0.9) | -0.6 (-0.9 to -0.2) | can't tell: significant, opposite sign |
| Stream minus shelf (W4) | +7.8% (+3.6 to +12.2) | +1.6 (+0.7 to +2.5) | yes, but see the buoy-by-buoy spread below |

- **Gusts** (secondary): buoy gust rises 7.1% (+4.5 to +10.7), more than mean wind, pooled over all 7 buoys.
- **Very stable air** (water colder than air) cuts wind at shelf buoys by 11% at the same gradient, as expected, but at the Stream buoys the same class is +4.2% (-0.5 to +9.3), not lower.
- The bump does not keep growing: the 9 °C+ class is 2.6% below the 6-9 class (W5, significant, opposite sign to the expectation).
- **Cold air off the continent is where it shows (post hoc lead, recomputed for the Stream figure only).** With wind from 270-360 degrees, the same comparison is +19.7% at the Stream buoys (roughly +4 kt at a 20 kt wind) and +4.2% at the shelf buoys. Pooled over all buoys, other wind directions show a 12% decrease (S1b, secondary, no group split, not recomputed). High air-sea temperature difference and northwest flow occur together, so any direction bias in the pressure-gradient control loads onto the temperature effect. Not a rule to apply to all flow.
- **The size depends on which buoy.** One buoy at a time, the effect ranges from -8.8% to +19.2% (41002 +14.9, 41001 +19.2, 41025 -7.0, 41013 -5.6, 44014 +7.1, 41004 -8.8, 44009 +9.7; post hoc). 41025 sits at Diamond Shoals and 41013 and 41004 near the coast, so wind direction and land probably contribute. The winter bootstrap does not allow for this buoy-to-buoy spread, so W4 is in effect three buoys against four, and the Stream-versus-shelf split may reflect buoy position (41025 is a Stream buoy at -7.0%, while shelf buoys 44014 and 44009 are +7.1% and +9.7%).

How to read it: with water well above the air temperature, and a cold continental flow crossing the Stream, the data support adding a few knots over the Stream above what the pressure gradient suggests. They do not support a uniform bump across the shelf, and they do not cleanly separate instability from the other differences between buoys.

**Waves: can't tell.** Wind against the current (steady wind of 20 kt or more, Stream water of 22 °C or warmer, three Stream buoys, 614 hours, 14 winters):

| Test | Result |
|---|---|
| V1 opposed minus following wave height | +0.5 ft (-0.8 to +1.7); raw means at 20-30 kt winds: 8.3 ft opposed, 7.7 ft following |
| V2 same, minus the same difference in non-Stream water | -2.9 ft (-6.4 to +1.5) |
| V3 placebo: same wind sectors at shelf buoys with no strong current | **+2.9 ft (+0.4 to +6.1), fails** |

The placebo difference is as large as anything we are trying to find, so wind-direction sectors differ in more than the current (fetch, swell from other storms). The comparison cannot isolate a current effect. The power run for waves is not reliable (its detection rates do not fit its own standard errors, and it assumes the sector comparison is unbiased, which the placebo shows it is not), so no statement is made about what size of effect would have been seen. The V1 interval (-0.8 to +1.7 ft) says only that this data cannot separate the current from other effects. At 30-40 kt there are 73 opposed or following Stream-water hours and above 40 kt only 5, so there is almost no information at the wave heights that matter most.

## Tests and corrections

Family W (5 tests): W1 -0.2 (q 0.79); W2 +4.7 (q 0.0085); W3 -2.8 (q 0.0077); W4 +7.8 (q 0.0050); W5 -2.6 (q 0.0065). By the plan's rule, only W2 and W4 are "yes" with the expected sign. W3 and W5 are significant in the opposite direction, so the family as a whole is mixed and does **not** support "yes" as a general rule. Family V (3 tests): V1 q 0.72, V2 q 0.72, V3 q 0.0006 (the placebo failure). Winter-block bootstrap, 10,000 draws; the 21 winters, not the 25,139 hours, are the sample. Secondary and post hoc results (S1-S4, D1-D5 in `results/`) are leads, not findings. A placebo that shuffles temperature within each buoy, winter and pressure-gradient tercile gives +1.4% at the Stream buoys (range -0.3 to +3.1) and +0.2% at the shelf buoys, so about 1 point of the Stream result can come from the method.

## Power and limits

- Power (planted-effect simulation, `power-result.txt`, unadjusted p < 0.05, correct sign): the pooled wind test finds +3% 86% of the time and +6% 100%. The Stream-only test finds +3% 32% of the time and +6% 95%. Wave power (V1 27% at +1 ft, 99% at +2 ft; V2 never above 40%) is shown in `power-result.txt` but not used: see the waves section. Wave verdict: can't tell.
- Pressure-gradient control (post hoc check D5): the wind-to-gradient ratio falls from 0.84 (Vg 15-25 kt) to 0.46 (50+ kt), so the log-linear control is imperfect and stability can leak into the estimate.
- The Stream is classed by buoy water temperature of 22 °C and a fixed current direction, which will misclassify meanders and rings. Probably dilutes the wave effect.
- Does not test whether GFS, the GFE blend or ECMWF already include the effect. The bump needed in a given model is a separate hindcast (the marine wind project).
- The ERA5 wave model, as far as we know, does not include ocean currents, so ERA5 waves are not used. No ocean current analysis was reachable from this container.
- Anemometer heights 3.2-4.1 m, adjusted to 10 m with a neutral log law; stability changes that adjustment by a few per cent.
- Buoy mean wind is a sustained measure; the gust comparison is a lead only.

## Verification

A fresh agent recomputed the counts, the median-wind table, the group wind effects, the per-buoy effects and the wave counts from the raw buoy and ERA5 files with its own code. Everything matched within 0.3 percentage points. Exceptions: the wave hour count (611 against 614, depends on how "within 30%" is read), and two buoys (44014, 41013) that move by 1-2 points depending on how the pressure grid is sampled at the buoy (nearest point, as used here, against bilinear). **Not recomputed by the verifier:** the power simulation, the secondary tests S1-S4 (S3 gust +7.1%, S2 Dec-Mar 0%), D2 beyond the Stream figure, D3 leave-one-winter-out (Stream +3.7 to +5.4%, shelf -3.3 to -2.5%), the D4 placebo, and the family q values (the reviewer checked that the q values follow from the listed p values). **Not independent:** the choice of the 22 °C and sector definitions was made in the plan.

## Reproduce

```
python3 research/era5/gulf_stream/pull_ndbc.py /mnt/project-files/gulf-stream/ndbc
python3 research/era5/gulf_stream/pull_msl.py /mnt/project-files/gulf-stream/msl 12      # 19.4 GB streamed
python3 research/era5/gulf_stream/build_dataset.py /mnt/project-files/gulf-stream/ndbc /mnt/project-files/gulf-stream/msl data/stationtimes.csv.gz
python3 research/era5/gulf_stream/power.py data/stationtimes.csv.gz power-result.txt
python3 research/era5/gulf_stream/analyse.py data/stationtimes.csv.gz results
python3 research/era5/gulf_stream/posthoc.py data/stationtimes.csv.gz results
```

(`data/stationtimes.csv.gz`, 2.1 MB, is committed. The NDBC and ERA5 inputs are in the shared folder, not in git.)
