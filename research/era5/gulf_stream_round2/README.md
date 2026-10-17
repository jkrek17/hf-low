# Gulf Stream marine wind, round 2: sharper tests of the satellite and GFS result

Plan: `PREREGISTRATION.md` (committed first, `fc4f007`, before any dewpoint, satellite sea-surface temperature, current or extra GFS field was read; four deviations logged there, items 1 to 3 before any outcome, item 4 after the primary results). Follows PR 269 (`../gulf_stream_region/`), after Jason asked how a PhD meteorologist would review that work and said "Go" to the three fixes proposed.

Data as in round 1 (ASCAT satellite winds in 1-degree boxes per pass, western Atlantic 29-41N 80-55W, October to April, 12 winters 2013-14..2024-25; GFS 24 h forecasts, 5 winters 2021-22..2025-26), plus ERA5 dewpoint, JPL MUR satellite sea-surface temperature, NOAA surface currents, and four more GFS fields. ERA5 is a control only and a proxy. No ERA5 gust, no hurricane-force definition, neither gust pipeline. ASCAT has no data above about 45 kt, so **nothing here applies at storm or hurricane force.** 24 h forecasts only.

## What a forecaster can use

**GFS at 24 h does not need a Gulf Stream "instability bump" added to its wind. When the satellite and the model are compared in the same terms, the gap does not grow with the water-minus-air difference: -0.9% over all boxes, -0.2 kt at a 25 kt model wind (-0.5 to +0.1). The earlier "add about 1 kt near the front" came from comparing a satellite wind that assumes neutral air with a model wind that does not.**

Satellite winds against the pressure gradient (what the atmosphere does, not what GFS does):

**The wind excess over what the pressure gradient suggests, when the water is 6 °C or more warmer than the air, is a cold-air-outbreak effect. With cold air advecting (about 9 in 10 of those times) the wind is about 15% higher (+3 kt at a 30 kt gradient wind). Without cold advection it is about 4% lower. It is not stronger over the Stream front than elsewhere.**

### GFS against the satellite, like with like (family F1, 5 winters)

Both winds are turned into the "neutral-air" wind that gives the same drag on the sea surface (COARE 3.6 bulk flux scheme). Wind change when the GFS water-minus-air difference is 6 °C or more versus 0 to 3 °C:

| Test | Change | Knots at model wind 25 kt | Plan verdict |
|---|---|---|---|
| B1e all boxes | -0.9% (-2.2 to +0.2) | -0.2 (-0.5 to +0.06) | no (well powered): within 1 kt and the test finds a planted +3% every time |
| B2e within 1 degree of latitude of the front | +0.3% (-1.7 to +1.8) | +0.1 (-0.4 to +0.5) | no, with the power taken from round 1's B2 run (90% at +3%), not rerun for this outcome |
| B3e 2 degrees or more from the front | -1.4% (-2.8 to -0.2) | -0.3 (-0.7 to -0.05) | can't tell (q 0.08, and the sign is opposite to a missing bump) |
| B4e near minus away | +1.7% (-0.7 to +3.8) | +0.4 (-0.2 to +1.0) | can't tell (q 0.21) |

- This reverses round 1's B1 and B2 (+2.2% and +3.6%). Those compared an equivalent-neutral satellite wind with a model wind that includes the stability effect. Like with like, the difference is zero.
- The same comparison with the satellite converted to true 10 m wind against the GFS 10 m wind gives the same answer (-0.9%, secondary).
- **Model-native check (secondary, `stress-check-result.txt`).** Using GFS's own surface stress instead (an 18-24 h average, so not the same time as the satellite pass, and including the model's gustiness) the GFS neutral wind rises 6.7% *more* than the satellite's with the temperature difference (-1.7 kt). It points the same way: no bump to add, if anything the model is too strong in unstable air. The time-average mismatch makes this a check, not a result.
- The registered near-minus-away test with the satellite front (B4m, +2.9%, +0.8 to +4.9) is "can't tell": only 4 of 5 leave-one-winter-out refits keep the interval above 0. Its near-front companion B2m (+4.0%, +2.3 to +5.3, 5 of 5 refits) is a "yes" under the plan, but it compares unlike winds; the like-with-like version with the same front (B2f) is +0.7% (-1.1 to +2.3).

### The excess wind in the satellite data (families F2 and F3, 12 winters)

Wind change when the water is 6 °C or more warmer than the air, versus 0 to 3 °C, at the same pressure-gradient wind. Reference wind 21.6 kt at gradient 30 kt. "Controls" are thermal advection and isobar curvature.

| Test | Change | Knots | Plan verdict |
|---|---|---|---|
| Round 1 baseline, all boxes (A1) | +15.9% (+14.0 to +17.9) | +3.4 | yes |
| A1a satellite converted to true 10 m wind | +12.2% (+10.2 to +14.4) | +3.1 | yes |
| A1c with air-mass controls | +10.3% (+8.6 to +12.2) | +2.2 | yes; +8.0% after the shuffle offset |
| A1-CAA cold advection (about 89% of the warm-water times) | +14.6% (+12.9 to +16.3) | +3.2 | yes |
| A1-nonCAA no cold advection | -4.4% (-6.2 to -2.3) | -0.9 | can't tell: significant, opposite sign |
| A1f all three fixes together (true wind, controls, satellite SST) | +6.9% (+5.2 to +8.7) | +1.5 | secondary |

Near the front versus away from it (plan expected near > away, from the buoy study):

| Specification | Near front | Away | Near minus away | Plan verdict |
|---|---|---|---|---|
| Round 1 (ERA5 SST front) | +8.2% | +19.1% | -9.2% (-11.5 to -7.0) | can't tell: significant, opposite sign |
| With air-mass controls (A2c to A4c) | +2.1% (+0.2 to +3.9) | +13.6% | -10.2% (-12.3 to -8.1) | same |
| Satellite SST front and temperature (A2m to A4m) | +11.1% | +18.4% | -6.2% (-8.5 to -3.9) | same |
| All three fixes (A2f to A4f, secondary) | +0.7% (-1.4 to +2.8) | +9.7% | -8.2% (-10.5 to -5.8) | same |

- **The air-mass controls do not explain the gap.** Cold advection explains why the effect exists (it is there with cold advection and absent without it), but the near-front boxes are lower than the away boxes in every version.
- **The satellite SST front shrinks the gap by about a third** (-9.2 to -6.2), which fits the old front having blurred the near-front group. It does not remove it.
- **Ocean current does not explain it either (secondary, 7 winters from 2018-19).** Adding the surface current (NOAA geostrophic, which leaves out the wind-driven part) to the satellite's water-relative wind moves the near-front estimate from +8.3% to +5.8% and the gap from -8.5 to -10.5. The current near the front averages 0.64 m/s against 0.27 m/s away.
- Almost all "away" boxes (98%) lie on the warm south side of the front, so "away" means the warm Sargasso side. Post hoc, not in the plan (`posthoc-latitude-result.txt`): inside 33-37N alone, near-front is -1.7% against +21.9% away, and inside 37-41N +11.4% against +15.6%, so the gap is not a latitude artefact, and is largest where "near the front" is the Gulf Stream off Cape Hatteras. This is a lead, not a finding.
- **Shuffle check.** If the temperature differences are shuffled among boxes of the same winter, month and pressure-gradient third, the method still reports +2.1% (A1c), +2.5% (A1m) and +0.4% (B1e). For A1c and A1m the plan's real-minus-offset rule applies, giving A1c +8.0% (+1.7 kt) and A1m +13.2% (+2.9 kt). For B1e the offset is under 1 point. The shuffle gives near-minus-away offsets of -0.2, +0.6 and -0.7 points, so the gap is not a method artefact.

### Reading the results together

- The pressure-gradient comparison measures something about the real atmosphere (cold air over warm water produces more wind than the gradient suggests), not something GFS misses. GFS already carries it in the quantity that matters here.
- "Stronger instability bump over the Gulf Stream" is not supported anywhere in this data: near the front the excess is the smallest group in every specification, and about zero in the all-fixes version.
- The answer to Jason's original question for the forecaster, "should I bump wind over the Stream in unstable air?": for GFS at 24 h, no, on these data; the model's wind already rises with the temperature difference as much as the satellite's. For waves the question is still untested.

## Power

`power-result.txt` (planted-effect simulation on the per-winter sufficient statistics, 200 by 400): a planted +3% where the temperature difference is 6 °C or more is found 100% of the time for the all-box tests (A1c, A1m, B1e), and a +3% planted only near the front is found 80% (A4c), 82% (A4m) and 81% (B4e) of the time; the difference tests correctly stay quiet (1 to 2%) when the same +3% is planted everywhere. This keeps day-to-day weather dependence only within winters, so it flatters the real tests somewhat. B2e was not simulated in this round.

## Limits

- The bulk-flux conversion uses fixed mean winter radiation and boundary-layer height, no sea-state effect on roughness, and 2 m humidity from ERA5 or GFS. The inversion from neutral to actual wind fails by more than 0.2 kt for 0.1% of rows (stable air, light wind), all in the colder-than-air class that no contrast uses.
- No boundary-layer depth or lapse-rate field was pulled (no pressure-level fields). Thermal advection and curvature are crude air-mass controls.
- The thermal advection is the pressure-gradient wind carrying the 2 m temperature gradient, computed at 0.5-degree spacing, and can be spiky next to the coast.
- MUR is a gap-filled analysis (true resolution 4-10 km). The box-mean SST and the front come from it.
- Geostrophic currents miss the Ekman part. Only 7 winters.
- B has 5 winters. ASCAT: two passes a day, rain-flagged cells dropped, no data above about 45 kt. GFS: one model, one lead. Waves are untested.

## Verification

To be filled in after the fresh-agent recompute and the claim review.

## Reproduce

```
python3 research/era5/gulf_stream_round2/pull_era5_d2m.py /mnt/project-files/gulf-stream/era5d2m
python3 research/era5/gulf_stream_round2/pull_mur.py /mnt/project-files/gulf-stream/mur
python3 research/era5/gulf_stream_round2/pull_currents.py /mnt/project-files/gulf-stream/currents
python3 research/era5/gulf_stream_round2/pull_gfs2.py /mnt/project-files/gulf-stream/gfs2
python3 research/era5/gulf_stream_round2/build_round2.py /mnt/project-files/gulf-stream /mnt/project-files/gulf-stream/tables work      # boxes2_A/B (not in git; copies in /mnt/project-files/gulf-stream/tables)
CONV_ONLY=1 python3 research/era5/gulf_stream_round2/analyse_round2.py work/boxes2_A.csv.gz work/boxes2_B.csv.gz work/cache work/results   # conversion cache only
python3 research/era5/gulf_stream_round2/power_round2.py work/boxes2_A.csv.gz work/boxes2_B.csv.gz work/cache power-result.txt
python3 research/era5/gulf_stream_round2/analyse_round2.py work/boxes2_A.csv.gz work/boxes2_B.csv.gz work/cache results
python3 research/era5/gulf_stream_round2/lowo_round2.py work/boxes2_B.csv.gz work/cache results/lowo-B-round2-result.txt
python3 research/era5/gulf_stream_round2/stress_check.py work/boxes2_B.csv.gz results/stress-check-result.txt
python3 research/era5/gulf_stream_round2/posthoc_round2.py work/boxes2_A.csv.gz results/posthoc-latitude-result.txt
python3 research/era5/gulf_stream_round2/test_bulk.py
```
