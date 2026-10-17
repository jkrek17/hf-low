# Pre-registration: Gulf Stream marine wind, round 2 (sharper tests of the satellite and GFS result)

Written and committed before any ERA5 dewpoint, satellite sea-surface temperature, ocean current or extra GFS field was read for this question. It follows `research/era5/gulf_stream_region/` (PR 269), whose boxes tables (`/mnt/project-files/gulf-stream/tables/boxes_A.csv.gz`, `boxes_B.csv.gz`) are reused unchanged for ASCAT winds, ERA5 pressure-gradient wind and temperatures, and GFS winds. No ERA5 gust, no hurricane-force definition, neither gust pipeline.

## Why (Jason, 2026-10-17: "How would a PhD meteorologist review your work", then "Go")

The reviewer-style critique of round 1 named three things that could change its conclusion. Each gets one fix here, as its own test family:

1. **The satellite-to-true-wind step decided the GFS answer.** Round 1 got +2.2% before and -1.9% after a crude conversion. Fix: use the COARE 3.6 bulk flux algorithm (`pycoare`), and compare like with like (stress-equivalent wind on both sides).
2. **The pressure-gradient yardstick ignores air mass.** High water-minus-air difference mostly means cold-air outbreaks. Fix: add thermal advection and isobar curvature as controls, and split by cold advection.
3. **The Stream front and the temperature difference came from a smooth 0.25-degree ERA5 SST map.** Fix: use daily satellite SST (JPL MUR, 0.05 degree) for both.

A fourth issue is added by the same logic: ASCAT measures wind relative to the moving water, and the Gulf Stream moves at up to 2 m/s. A secondary test adds the geostrophic surface current.

## What has been seen

- Round 1 results (PR 269), including the numbers that motivated this round: A1 +15.9%, A2 +8.2%, A3 +19.1%, A4 -9.2% (opposite sign), B1 +2.2%, B2 +3.6%, S5 after the corrected conversion B1 -1.9%; the shuffle check gives +2.5 points on A1; 97.5% of "away" boxes lie south of the ERA5 front.
- Sizing only for the new data: the GFS index lists DPT, PRES, UFLX, VFLX; ERDDAP lists MUR (`jplMURSST41`, 2002 onward) and NOAA geostrophic currents (`nesdisSSH1day`, `ugos`/`vgos`, from February 2017); `pycoare` installs and runs.
- Not seen: any dewpoint, any MUR value, any current, any COARE-converted wind in this region.

## Data to pull (about 25 GB streamed, under the 50 GB gate, nothing cached covers it)

- ERA5 2 m dewpoint, 00 and 12 UTC, 1 Oct-30 Apr, winters 2013-14..2024-25 (12), cropped as in round 1: 12 x 424 x 2.35 MB = 12 GB. Surface pressure is taken as ERA5 MSLP already pulled (error under 1% of the stability terms).
- MUR daily SST, 25-45N, 85-50W, stride 5 (0.05 degree), the same 12 winters, ERDDAP: about 2,500 days x 1.7 MB = 4.3 GB (probe: 0.57 MB per day at stride 10).
- GFS 24 h forecasts (valid 00 and 12 UTC, 2021-22..2025-26): add `DPT:2 m`, `PRES:surface`, `UFLX:surface` and `VFLX:surface` (18-24 h average, a secondary), 2,120 times x 4 fields x about 0.84 MB = 7 GB.
- NOAA geostrophic currents (`ugos`, `vgos`), daily, winters 2017-18..2024-25 (8): about 0.5 GB.
- Stop and ask if the measured total passes 35 GB.

## Method

**Bulk conversion.** `pycoare.coare_36`, wind at 10 m, air temperature and humidity at 2 m (relative humidity from 2 m temperature and dewpoint, Magnus formula), pressure from MSLP, latitude at the box centre, fixed mean winter radiation (shortwave 100, longwave 310 W/m2) and boundary-layer height 600 m (these only enter the cool-skin and gustiness terms), ocean temperature as sea surface temperature with the cool skin switched on for ERA5 and MUR (bulk or subskin values) and off for the GFS skin temperature. The neutral 10 m wind (`u_n_rf`) is the stress-equivalent wind. Forward direction (actual to neutral) is used for GFS. For ASCAT (neutral to actual) the actual wind is found by three secant steps so that its neutral wind equals ASCAT's; a test checks that the conversion returns the input at zero air-sea difference and lowers the wind in unstable air.

**Round 2 boxes.** Same 1-degree box-pass units, bins (<0, 0-3 reference, 3-6, 6-9, >=9 C), box fixed effects, contrast (mean of the 6-9 and >=9 bins minus reference), winter-block bootstrap (12 winters; 5 for B, 10,000 draws), p-values and Benjamini-Hochberg q within family and across all tests. Verdict words as in round 1: yes (q < 0.05 and the expected sign), can't tell, can't tell: significant but opposite sign, no (well-powered) only if the 95% interval lies within +-1 kt at 25 kt and power is at least 80% at a planted +3%.

**Family F1, like-with-like wind (5 tests).** B1e to B4e: `log(U_EN_ascat / U_EN_gfs)` against the GFS dT bins (GFS skin minus 2 m temperature), all boxes, in band, away, in minus away (bands as round 1, ERA5 front). A1a: round 1's A model with ASCAT converted to actual 10 m wind by COARE, all boxes. Secondary: B on the actual-wind side (ASCAT converted to actual 10 m against GFS 10 m); direct GFS stress (UFLX, VFLX, 18-24 h average) turned into a neutral wind and compared with ASCAT, as a check of the conversion route.

**Family F2, air mass (6 tests).** Round 1's A model plus controls: thermal advection of 2 m temperature by the pressure-gradient wind (`-Vg . grad T2m`, in K per day, central differences of +-0.5 degree, five bins: <= -4, -4 to -1, -1 to +1 (reference), +1 to +4, >= +4) and curvature (Laplacian of MSLP, three equal-count bins). Tests: A1c all boxes, A2c in band, A3c away, A4c in minus away; A1-CAA (advection <= -1 K/day, cold air) and A1-nonCAA (> -1), each with the same controls. Reading rule written now: if the away-minus-near gap (A4) shrinks to a 95% interval that includes 0 after the controls, the round 1 gap is read as air mass; if it stays about the same, it is not explained by cold advection.

**Family F3, front and temperature from satellite SST (6 tests).** Front per 1-degree longitude strip each day from MUR SST smoothed with a 0.5-degree box: the latitude of the strongest warm-to-south gradient between 32N and 42N, at least 0.01 K/km, undefined on the window edge (round 1's rule). dT = box-mean MUR SST minus ERA5 2 m temperature (A) or GFS 2 m temperature with MUR SST replacing the GFS skin value only in a secondary. In band: within 1 degree of latitude of the MUR front; away: 2 or more. Tests: A1m, A2m, A3m, A4m (round 1's A model with MUR dT and MUR front, no air-mass controls); B2m, B4m (round 1's B with GFS dT and the MUR front bands, GFS winds in actual terms). Secondary: the same with the controls of F2 and the conversion of F1 together, labelled "all three fixes".

**Currents (secondary, 2017-18 onward only).** ASCAT wind vector relative to the surface plus the NOAA geostrophic current at the box centre gives an absolute wind; round 1's A1/A2/A3 repeated on that wind. Reports whether the Stream's motion explains any of the near-front gap. Geostrophic currents omit the Ekman part, so this is a lower bound on the current effect.

**Placebo.** Shuffle dT among boxes of the same winter, month and Vg tercile (A) or winter, month (B), 20 draws, for A1c, A1m and B1e. If the shuffled mean differs from 0 by more than 1 point, the real result minus the offset is reported beside it, as in round 1. The offset for the near/away difference tests is also reported (shuffled A4 variants), because round 1 left it unknown.

**Multiplicity.** Families F1 (5), F2 (6), F3 (6): 17 tests. BH within family and across all 17.

## Power

Before any coefficient: the planted-effect simulation of round 1 (remove each band group's own fitted dT effect, plant +3% and +6% where dT >= 6, 100 outer x 200 inner winter resamples) for A1c, A4c, A4m, B1e and B4e. If power for A4c or A4m at a planted difference between groups of +3% is under 80%, the answer for "does the front matter" is can't tell regardless of the estimate.

## What this cannot show

- Boundary-layer depth and lapse rate (no pressure-level fields pulled); advection and curvature are only crude air-mass controls.
- COARE assumptions: fixed radiation and boundary-layer depth, no sea-state or wave-age effect on roughness, no rain, MUR is a subskin analysis interpolated under cloud (true resolution 4-10 km).
- Geostrophic currents without the Ekman part; absolute-current test covers 8 winters only.
- ASCAT sampling (two passes a day, rain-flagged cells dropped, no data above about 45 kt) and GFS as one model at one lead are unchanged. Nothing here applies at storm force. Waves are still untested.

## Deviations (post hoc)


1. *Before any outcome.* Power simulation done on the per-winter sufficient statistics (a winter resample is a count vector, as in the analysis bootstraps), 200 outer x 400 inner, instead of refitting resampled rows at 100 x 200. Same planted design (each band group's own fitted dT effect removed, +3% and +6% planted). Output: `power-result.txt`.
2. *Before any outcome.* 0.1% of A rows (stable air, light wind, dT near -3) do not invert within 0.2 kt (worst 1.5 kt in a 100,000-row sample); they stay in the sample and sit in the dT < 0 bin, outside every contrast. MUR for winter 2025-26 was pulled after the first table build (the B sample covers it).
3. *Before any outcome.* The NOAA geostrophic-current dataset has data only from February 2018, so the current secondary covers 7 winters (2018-19 to 2024-25), not 8. Noticed while building the tables.
4. *After the primary results.* Added, all labelled: `lowo_round2.py` (the registered 5-of-5 interval rule for the B tests, also for B1e, B2e, B4e); `stress_check.py` (the GFS native-stress route, listed in the plan as a secondary; run after seeing the primary B1e); `posthoc_round2.py` (post hoc, not in the plan: near-minus-away inside 4-degree latitude bands, run because the gap persisted after every registered fix).

