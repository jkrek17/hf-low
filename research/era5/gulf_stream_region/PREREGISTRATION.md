# Pre-registration: Gulf Stream and marine wind, region-wide (ASCAT) and against GFS

Written and committed before any ASCAT wind, any ERA5 air or sea temperature, or any GFS field was read for this
question. It extends `research/era5/gulf_stream/` (7 buoys, PR 263), which is cited and not repeated. No ERA5 gust, no
hurricane-force definition, neither gust pipeline.

## What has been seen

- Buoy study (PR 263): pooled water-minus-air effect about 0; Stream buoys about +5%, shelf buoys about -3%, one buoy at
  a time -9% to +19% (post hoc). Post hoc lead: northwest flow at the Stream buoys +20%. Test **A5** below registers
  that lead as a confirmatory test on new data, so a positive A5 is not a fresh discovery of the same data.
- Sizing only: ASCAT inventory file sizes, the GFS and ERA5 field sizes, and that the sources are reachable. Not seen:
  any ASCAT wind in this region and period, any ERA5 temperature field for this question, any GFS wind.

## Question (Jason, 2026-10-17: "a few buoys isn't much of a test; what about high-res model data?")

A) Across the whole western Atlantic Stream, how much stronger is the satellite-observed wind when the water is much
warmer than the air, at the same pressure-gradient wind, and does it differ close to the Stream's north wall?
B) How much does GFS (0.25 degree) forecast wind miss, against the same satellite winds, as a function of the
model's own air-sea temperature difference? That is the bump a forecaster would add to the model.
Plain answers to give: how many knots to add, yes / no (well-powered) / can't tell.

## Data

- **ASCAT**: RSS ASCAT v02.1 daily 0.25 degree bytemaps, Metop-A, -B, -C (`/mnt/project-files/ascat-check/rss_inventory.csv`),
  October to April. Wind is the scatterometer's **equivalent-neutral** 10 m wind (see "What this cannot show"). A cell is
  used if wind code is valid, scatterometer rain flag is 0 and radiometer rain code is 0.
- **Region**: 29-41N, 80-55W, ocean only. Cells within 0.5 degree of land (ERA5 land-sea mask, `era5-cache/lsm_nh_0p25.npy`,
  any land fraction >= 0.5 within two cells) are dropped.
- **Unit**: 1 x 1 degree box mean per satellite pass, requiring at least 8 valid cells of 16. Cells in one box and pass are
  not independent, and neither are neighbouring boxes.
- **A seasons**: winters 2013-14 to 2024-25 (12), when Metop-B is operating and ERA5 MSLP crops exist. **B seasons**: 2021-22
  to 2025-26 (5): the GFS bucket starts October 2021.
- **Controls**: ERA5 MSLP (crops in `/mnt/project-files/gulf-stream/msl/`) at 00 and 12 UTC gives the pressure-gradient wind
  `Vg` at the box centre, central differences +-0.5 degree, **bilinearly interpolated to the box centre** (the buoy study
  showed nearest-grid sampling moves some single-station results by 1-2 points), rho = 1.25 kg/m^3. ERA5 2 m
  temperature at 00 and 12 UTC and ERA5 SST at 00 UTC (daily) give `dT = SST - T2m`, clipped to [-3, 12] degrees C. The
  control fields are matched to the pass by the nearest 00 or 12 UTC within 4 h; others are dropped. ERA5 assimilates ASCAT,
  but wind speed, direction and ASCAT enter neither `Vg` nor `dT`; SST and T2m are used only to classify, not as outcomes.
- **Stream band**: from the daily 00 UTC ERA5 SST, smoothed with a 1 degree box, the latitude of the strongest northward
  SST decrease between 32N and 42N at each 1 degree of longitude, if that gradient is at least 1 C per 100 km (otherwise the
  longitude is undefined and dropped). **In band** = box centre within 1 degree of latitude of that front; **outside** = at
  least 2 degrees away, either side. In between is not used in the band tests.
- **B model**: GFS 0.25 degree, **24 h forecast**, valid 00 and 12 UTC, fields `UGRD`, `VGRD` 10 m, `TMP` 2 m, `TMP:surface`
  (skin temperature, equal to SST over open ocean), by HTTP range request from the NOAA GFS AWS bucket. ASCAT boxes within
  4 h of the valid time are matched to the model value at the box centre (bilinear). Model `dT` = surface minus 2 m
  temperature, clipped as above.
- **Pull plan**: ASCAT 9.2 GB (inventory sizes, A+B+C, October-April 2013-10 to 2026-04); ERA5 SST 3.7 GB (13 winters x 212
  days x 1.36 MB, 00 UTC) and T2m 13.0 GB (13 x 424 x 2.35 MB); GFS 7.1 GB (2,120 valid times x 3.35 MB measured from the
  index); total about 33 GB, under the 50 GB gate and below the 35 GB announced to Jason. Stop and ask above 40 GB.

## Model A and tests (12 winters)

For boxes with `Vg >= 15 kt` and ASCAT wind `>= 5 kt`: `log(V_ascat) = a_box + b log Vg + g(dT) + e`, with five `dT` bins
(<0; 0-3 reference; 3-6; 6-9; >=9 C) and a fixed effect for each 1 degree box. **W-contrast** = mean of the 6-9 and >=9
bins minus the reference bin, in log wind. Same contrast, three groups:

Family **A** (5 tests):
- A1: all boxes.
- A2: in band.
- A3: outside the band.
- A4: A2 minus A3.
- A5: in band and flow from 270-360 degrees (ASCAT direction), the buoy lead; expected positive.

Secondary (labelled, raw p): S1 effect per +5 C of `dT` as a continuous term; S2 Dec-Mar only; S3 each satellite
separately; S4 only winds >= 30 kt (ASCAT saturates near 50 kt, so no data above 45 kt); S5 an approximate conversion from
equivalent-neutral to actual 10 m wind with a bulk-Richardson / Businger-Dyer stability function (the expected effect is
to shave a few per cent off unstable-air winds).

## Model B and tests (5 winters)

`log(V_ascat / V_GFS) = a_box + g(dT_GFS) + e`, same bins and boxes, no `Vg` term (the model already contains the
pressure field). **B-contrast** = the high-`dT` bins minus the reference, which is the extra per cent a forecaster adds to
the model. In knots at model wind 25 and 35 kt. Family **B** (4 tests): B1 all boxes; B2 in band; B3 outside the band;
B4 B2 minus B3. Secondary: B5 NW flow in band; B6 24 h versus a 48 h forecast (a second, smaller pull of the same size, only
if B1 passes). A bias that is the same in every `dT` bin is a model wind bias, not a Stream effect, and is reported as such
(the reference-bin ratio).

## Inference and multiplicity

Effective n is winters (12 for A, 5 for B), not boxes or passes. Resample winters with replacement, 10,000 draws,
percentile intervals; with 5 winters the B intervals are reported as wide by construction and the result is "can't tell"
unless the interval excludes 0 in 5 of 5 leave-one-winter-out refits. Benjamini-Hochberg within family A and family B; q also
across all 9. Placebo (must be run and reported): shuffle `dT` among boxes of the same winter, month and `Vg` tercile, 20
times; if the shuffled mean differs from 0 by more than 1 point, report the real result minus that offset beside it. A
result stands as "yes" only if its sign is the expected one and q < 0.05; a significant result in the opposite direction is
reported as "can't tell: opposite sign" (lesson from PR 263).

## Power

Planted-effect simulation before any coefficient is read: add +3% and +6% to wind where `dT >= 6`, noise from residuals
taken after the full fit, with the 12 (A) or 5 (B) winters resampled, 200 outer x 400 inner. The decision rule for a
"well-powered no": the 95% interval inside +-1 kt at 25 kt wind and power at least 80% at a planted +4%. If not, the
answer is "can't tell". The sim is an upper bound on power because it keeps the real day-to-day dependence only within
winters.

## What this cannot show

- ASCAT reports an **equivalent-neutral** wind, which is the wind that would give the same surface stress in neutral air. In unstable air the actual 10 m wind is
  lower by roughly 2-8% (more at light wind and large `dT`). An ASCAT-derived bump therefore overstates the 10 m wind bump
  and is closer to the stress bump that drives waves; S5 gives the approximate conversion. B compares to a model 10 m wind,
  which includes stability, so B has the same issue and S5 applies to it too.
- ASCAT footprints are about 25 km; the north wall is a few tens of km wide. Meander position in the band tests is
  approximate (daily 00 UTC ERA5 SST at 0.25 degree).
- ASCAT is two passes a day; a bump present only at night or in weather ASCAT does not sample is not seen. Rain-flagged cells
  are dropped, so heavy-rain convective winds are missing.
- Neither A nor B says anything about waves. A wave test with the real current needs the CoastWatch currents and buoy waves
  and was left out by Jason's choice (A and B only).
- GFS is one model at one lead; the GFE blend project's own result is not repeated here.

## Deviations (post hoc)

Logged in order, all *before any outcome* (only decoding checks and descriptive counts of the control fields had been
looked at):

1. ASCAT direction byte is oceanographic (direction the wind blows towards). Checked on 348 buoy matches (41013, 44014,
   41002, 41004) in the 2013-14 winter: median direction difference 174 deg as stored, 6 deg after +180. The builder adds 180
   degrees before the NW-flow split. (ASCAT / buoy wind ratio 1.06, correlation 0.89 on the same matches.)
2. The Stream-band front is set undefined when the strongest gradient sits on the edge of the 32-42N search window (the
   first descriptive look showed it pinned at 42N for every longitude east of 65W, where the window does not contain the
   north wall). This removes those longitudes from the band tests (A2-A5, B2-B4); A1 and B1 still use every box.
3. Secondary S4 as written ("only winds >= 30 kt") would select on the outcome. It is run instead as the pressure-gradient
   wind `Vg >= 30 kt` (a control, not the outcome); no ASCAT-wind cut is used anywhere.
4. Implementation: the 1 degree box fixed effects are absorbed by within-box demeaning (exact for the full-sample OLS fit);
   winter resamples reuse the full-sample box means (an approximation to refitting the dummies, chosen because the table has
   1.3 million box-passes). The planted-effect simulation removes the fitted dT-bin effect from the outcome and adds the planted
   effect to the remaining residuals, which is the same as base fit plus residual.
5. The first power run removed the fitted dT effect of one pooled fit before planting; it returned zero detection for the
   in-band test at a planted +6% (impossible if the in-band effect equalled the pooled one), so the groups differ and the pooled
   removal was the wrong design (this first run is also one indirect look at the direction of the A4/B4 difference: in-band A below the pooled effect, in-band B above outside) for the group tests. Power is re-run removing each band group's own fitted dT effect (in, out,
   mid, undef). No coefficient was read or printed; the first run's output is kept in `power-result-first-run-discarded.txt`.
6. The ASCAT pull is slower than sized (about 3.7 min per satellite-winter); the sizes and the pull plan are unchanged.
7. *Before any outcome.* Power simulation run at 100 outer x 200 inner (the plan said 200 x 400), for run time on 1.3M box-passes.
8. *After the primary results.* The analysis script printed only the point range for leave-one-winter-out, so `lowo_B.py` was added to compute the registered 5-of-5 interval rule for B, and `s5_B.py` to apply the registered S5 conversion to B (the plan says S5 applies to B; the first script ran it for A1 only). Both are registered quantities, run after seeing the family results; they could not change a primary estimate.
9. *After the primary results.* The first S5 (equivalent-neutral to actual wind) used the bulk Richardson number itself as the stability parameter and so understated the correction roughly tenfold (a reviewer caught it: at 12 m/s the correction for water 7.5 versus 1.5 C warmer than air was 0.5 points, should be about 3). It is replaced by `region_lib.convert_en` (zeta = 10.8 x Rb), run by `s5_conversion.py` and checked by `test_convert_en.py`. The first version's output is kept in `results/s5-first-version-result.txt`. S5 is a secondary; no primary estimate changed. Also corrected in the analysis script: a printed string called the ASCAT/GFS ratio "same in every bin" as if found, which it was not.
