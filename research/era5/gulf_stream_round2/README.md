# Gulf Stream marine wind, round 2: sharper tests of the satellite and GFS result

Plan: `PREREGISTRATION.md` (committed first, `fc4f007`, before any dewpoint, satellite sea-surface temperature, current or extra GFS field was read; four deviations logged there, items 1 to 3 before any outcome, item 4 after the primary results). Follows PR 269 (`../gulf_stream_region/`), after Jason asked how a PhD meteorologist would review that work and said "Go" to the three fixes proposed.

Data as in round 1 (ASCAT satellite winds in 1-degree boxes per pass, western Atlantic 29-41N 80-55W, October to April, 12 winters 2013-14..2024-25; GFS 24 h forecasts, 5 winters 2021-22..2025-26), plus ERA5 dewpoint, JPL MUR satellite sea-surface temperature, NOAA surface currents, and four more GFS fields. ERA5 is a control only and a proxy. No ERA5 gust, no hurricane-force definition, neither gust pipeline. ASCAT has no data above about 45 kt, so **nothing here applies at storm or hurricane force.** 24 h forecasts only.

## What a forecaster can use

**Over all boxes, GFS at 24 h needs no added Gulf Stream bump. Once the satellite wind and the model wind are put in the same terms (COARE 3.6 bulk flux scheme), the satellite-minus-GFS gap does not grow with the water-minus-air difference: -0.9%, -0.2 kt at a 25 kt model wind (-0.5 to +0.1 kt), a well-powered "no". Close to the front this is not settled: the like-with-like test is within about half a knot but its own power was not simulated, and the registered test with the satellite front on unconverted winds (B2m) is a "yes" at +4.0%. The earlier "add about 1 kt near the front" disappears when both winds are compared in the same terms. Winds below about 45 kt, 5 winters; nothing here applies at storm or hurricane force.**

Satellite winds against the pressure gradient (what the atmosphere does, not what GFS does):

**When the water is 6 °C or more warmer than the air, the satellite wind runs above what the pressure gradient suggests: +10% after the air-mass controls (+8%, about +1.7 kt, after the shuffle offset; +7%, about +1.5 kt, with all three fixes). The excess is found when cold air is being carried over the water (air temperature advection of -1 K/day or colder; about 9 in 10 of these cases): +15% in the satellite's neutral-equivalent wind, about +3 kt on a 22 kt wind. Without cold advection it is -4% (can't tell: significant, opposite sign). Near the front the excess is smaller than on the warm side 2 or more degrees south in every version (can't tell: significant, opposite sign), so there is no sign of a stronger effect at the front.**

### GFS against the satellite, like with like (family F1, 5 winters)

Both winds are turned into the "neutral-air" wind that gives the same drag on the sea surface (COARE 3.6 bulk flux scheme). Wind change when the GFS water-minus-air difference is 6 °C or more versus 0 to 3 °C:

| Test | Change | Knots at model wind 25 kt | Plan verdict |
|---|---|---|---|
| B1e all boxes | -0.9% (-2.2 to +0.2) | -0.2 (-0.5 to +0.06) | no (well powered): within 1 kt and the test finds a planted +3% every time |
| B2e within 1 degree of latitude of the front | +0.3% (-1.7 to +1.8) | +0.1 (-0.4 to +0.5) | can't tell (interval -0.4 to +0.5 kt, but power was not simulated for this test; round 1's B2, a different outcome, had 90% at +3%) |
| B3e 2 degrees or more from the front | -1.4% (-2.8 to -0.2) | -0.3 (-0.7 to -0.05) | can't tell (q 0.08 within the family, 0.039 across all 17; the sign is opposite to a missing bump) |
| B4e near minus away | +1.7% (-0.7 to +3.8) | +0.4 (-0.2 to +1.0) | no (well powered), marginal: upper end +0.96 kt against the 1 kt bound, power 0.81 from a simulation that is an upper bound; q 0.21 |

- This replaces round 1's B1 and B2 (+2.2% and +3.6%). Those compared an equivalent-neutral satellite wind with a model wind that includes the stability effect. Like with like, no change larger than about 0.5 kt (B1e -0.9%, -2.2 to +0.2). The leave-one-winter-out rule (5 of 5) is applied to "yes" verdicts only, as in round 1; B1e has 1 of 5 refits below 0.
- The same comparison with the satellite converted to true 10 m wind against the GFS 10 m wind gives the same answer (-0.9%, secondary).
- **Model-native check (secondary, `stress-check-result.txt`).** Using GFS's own surface stress instead (an 18-24 h average, so not the same time as the satellite pass, and including the model's gustiness) the GFS neutral wind rises 6.7% *more* than the satellite's with the temperature difference (-1.7 kt). Run after seeing B1e (deviation 4). Its sign agrees but its size does not (-6.7% against -0.9%), so this route neither confirms nor overturns the conversion; the time-average mismatch makes it a check, not a result.
- The registered near-minus-away test with the satellite front (B4m, +2.9%, +0.8 to +4.9) is "can't tell": only 4 of 5 leave-one-winter-out refits keep the interval above 0. Its near-front companion B2m (+4.0%, +2.3 to +5.3, 5 of 5 refits) is a "yes" under the plan, but it compares unlike winds; the like-with-like version with the same front (B2f) is +0.7% (-1.1 to +2.3).

### The excess wind in the satellite data (families F2 and F3, 12 winters)

Wind change when the water is 6 °C or more warmer than the air, versus 0 to 3 °C, at the same pressure-gradient wind. Reference wind 21.6 kt (the median satellite wind at 25-35 kt gradient wind). "Controls" are thermal advection and isobar curvature.

| Test | Change | Knots | Plan verdict |
|---|---|---|---|
| Round 1 baseline, all boxes (A1) | +15.9% (+14.0 to +17.9) | +3.4 | yes |
| A1a satellite converted to true 10 m wind | +12.2% (+10.2 to +14.4) | +3.1 | yes |
| A1m satellite SST for temperature difference | +16.0% (+14.2 to +18.1) | +3.5 | yes; +13.2% after the shuffle offset |
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

- **These crude air-mass controls do not explain the gap.** The effect is found with cold advection (+14.6%) and not without it (-4.4%, can't tell: opposite sign); the controls shrink it (A1c +10.3%) but do not remove it, and the near-front boxes are lower than the away boxes in every version. Near the front on its own, A2c is a registered "yes" at +2.1% (q 0.029); A2m and A3m are "yes" at +11.1% and +18.4%.
- With the satellite SST front and temperature difference the gap is -6.2% instead of -9.2% (both the front and the temperature changed; the change itself was not tested, and the intervals overlap). It does not go away.
- **The geostrophic current does not explain it either (a lower bound on the current effect) (secondary, 7 winters from 2018-19).** Adding the surface current (NOAA geostrophic, which leaves out the wind-driven part) to the satellite's water-relative wind moves the near-front estimate from +8.3% to +5.8% and the gap from -8.5 to -10.5. The current near the front averages 0.64 m/s against 0.27 m/s away.
- Almost all "away" boxes (98%) lie on the warm south side of the front, so "away" means the warm Sargasso side. Post hoc, not in the plan (`posthoc-latitude-result.txt`): inside 33-37N alone, near-front is -1.7% against +21.9% away, and inside 37-41N +11.4% against +15.6%, so the gap persists inside 4-degree latitude bands and is largest at 33-37N, off Cape Hatteras. Away boxes within a band are still mostly the warm south side. This is a lead, not a finding.
- **Shuffle check.** If the temperature differences are shuffled among boxes of the same winter, month and pressure-gradient third, the method still reports +2.1% (A1c), +2.5% (A1m) and +0.4% (B1e). For A1c and A1m the plan's real-minus-offset rule applies, giving A1c +8.0% (+1.7 kt) and A1m +13.2% (+2.9 kt). For B1e the offset is under 1 point. The shuffle gives near-minus-away offsets of -0.2, +0.6 and -0.7 points, so the gap is not produced by the fitting design. The shuffle does not test sampling that differs between near and away boxes: rain-flagged cells dropped over the Stream, the ASCAT ceiling near 45 kt, coastal boxes off Hatteras, and MUR gap-filling under cloud; those remain untested.

### Reading the results together

- The pressure-gradient comparison measures something about the real atmosphere (cold air over warm water produces more wind than the gradient suggests), not something GFS misses. GFS already carries it in the quantity that matters here.
- In the satellite-versus-pressure-gradient data the excess is smaller near the front than away from it in every specification (to about zero with all three fixes). In the GFS comparison the near-minus-away estimates are positive in every version (B4 +2.2, B4e +1.7, B4m +2.9, B4f +2.4%) but none meets the plan's rules, and B2m, on unconverted winds, is a "yes".
- For the forecaster's question, "should I bump GFS wind over the Stream in unstable air?": over all boxes, no, on these data; the model's wind rises with the temperature difference to within about 0.5 kt of the satellite's. Close to the front it is unresolved. For waves the question is untested.

## Power

`power-result.txt` (planted-effect simulation on the per-winter sufficient statistics, 200 by 400): a planted +3% where the temperature difference is 6 °C or more is found 100% of the time for the all-box tests (A1c, A1m, B1e), and a +3% planted only near the front is found 80% (A4c), 82% (A4m) and 81% (B4e) of the time; the difference tests correctly stay quiet (1 to 2%) when the same +3% is planted everywhere. The simulation is an upper bound (it keeps weather dependence only within winters and ignores the 5-of-5 rule). The 80% for A4c sits at the plan's threshold and is printed to two decimals, so it could be just under. B2e was not simulated in this round.

## Limits

- The bulk-flux conversion uses fixed mean winter radiation and boundary-layer height, no sea-state effect on roughness, and 2 m humidity from ERA5 or GFS. The inversion from neutral to actual wind fails by more than 0.2 kt for 0.1% of rows (stable air, light wind), all in the colder-than-air class that no contrast uses.
- No boundary-layer depth or lapse-rate field was pulled (no pressure-level fields). Thermal advection and curvature are crude air-mass controls.
- The thermal advection is the pressure-gradient wind carrying the 2 m temperature gradient, computed at 0.5-degree spacing, and can be spiky next to the coast.
- MUR is a gap-filled analysis (true resolution 4-10 km). The box-mean SST and the front come from it.
- Geostrophic currents miss the Ekman part. Only 7 winters.
- B has 5 GFS seasons across model upgrades; A mixes ASCAT-A/B/C. B has 5 winters. ASCAT: two passes a day, rain-flagged cells dropped, no data above about 45 kt. GFS: one model, one lead. Waves are untested.

## Verification

A fresh agent (no access to the results or this README) recomputed from `boxes2_A/B.csv.gz` with its own regression: A1, A1a, A1c, A1-CAA, A1-nonCAA, the cold-advection share, A1m, A2m, A3m, A4m, B1e, B2e, B2m, Cabs1 and the A1c shuffle offset (5 draws). All 15 matched within 0.04 points. Notes from that check: A4m is the percent of the difference of log contrasts, not 11.07 - 18.44 (which would be -7.4); the B1e/B2e conversion was checked against the same cached function, so that part is not independent; Cabs1 reads ASCAT direction as "from" (+16.2 if read as "to"), uses all 724,984 rows that have a current and plain ASCAT plus the geostrophic current (on the 7-winter sample plain A1 is +15.0).

An independent claim reviewer found overstated verdict words and headlines; all were reworded here (B2e and B4e verdicts, the two forecaster headlines, "explains", "not a latitude/method artefact", the stress-check reading, missing A1m row).

Not independently recomputed: bootstrap intervals, the power simulation, the leave-one-winter-out intervals, the stress check, the post hoc latitude bands, the placebo ranges, the q values, the A2c/A3c/A2f/A3f rows and the current-sample rows other than Cabs1. B2m and B4m knots in `round2-result.txt` use the A reference wind (21.6 kt), not B's 25 kt; this README quotes percentages only.

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
