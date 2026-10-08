# Is the unexplained part of the HF conversion step noise at the threshold? Pre-registration (agenda RA-27)

Written and committed **before any ERA5 or HRES wind field for this sample was read**. ERA5 proxy throughout; the outcome is
**pipeline A** (`research/era5/hf_history`: 800 km ocean gust index, HF-equivalent at 71.7 kt), a proxy for the archive, not the
archive. Transitioning tropical cyclones are in. Atlantic fixes north of 60N are lower confidence (no terrain mask).

**What had been looked at.** Only labels, positions and sizes: the committed fix table (`hf_probability_tracks/results/fix_probs_2004.csv.gz`:
position, `g800`), the track table, the variable list of the second system's store, and HEAD-request sizes (`results/sizes.txt`).
No wind field of either system at these times.

## Feasibility finding that changes the test (reported first)

The question needs a second analysis run "through the same crossing rule". The rule is a threshold on the **ERA5 instantaneous 10 m
gust**. WeatherBench2's IFS HRES t0 analyses (`hres_t0`, 2016-01-01 to 2023-01-08, 0.25 degree, 6-hourly, both the 0.25 degree
archive and the regridded sets) carry **no gust variable** (variables: 10 m u, v, wind speed; 2 m temperature; MSLP; surface
pressure; total precipitation; pressure-level fields). So the rule cannot be transferred unchanged. What can be done, and what the
agenda named as the fallback, is to run the **same construction on sustained 10 m wind** in both systems. That tests how reproducible
a threshold crossing of a resolved-wind index is between two analyses; it does not test the gust index itself, whose gust
parametrisation is a model diagnostic that HRES's analysis does not provide here. I will say so wherever the result is quoted.

## Question and the answer I will give

If a second analysis is run through the same crossing construction, do the same cyclones cross? Answer, one of: **mostly noise**
(low agreement), **not noise** (high agreement: crossing is a reproducible storm property), or **can't tell / in between**, with the
agreement number. It is compared with the skill the forecast models reach (PR 12 HSS 0.59 for P(HF within 24 h)) to say how much of
the unexplained part the noise ceiling can account for.

## Sample (fixed now)

- Pipeline A tracks, seasons 2016-17 to 2021-22 (June to May; the HRES store ends 2023-01-08, so 2022-23 is incomplete and left out).
- A track is in the sample when its maximum 00/12 UTC in-domain gust index `g800` (as in the fix table) is at least the **floor**.
  Floor rule, decided now: **55 kt** (the agenda's near-threshold band starts there) if the pull at 55 kt is at most 40 GB,
  otherwise 60 kt (but see deviation 1: 65 kt was used). The sample is all 00/12 UTC in-domain fixes of those tracks (00/12 only because those are the fixes in the table;
  06/18 UTC fixes are not used in either system, so timing is identical in both).
- Positives for matching: tracks in the sample with gust index at least 71.7 kt (count per `results/sample_fixes_<floor>.csv`).

## Indices (frozen)

Fix-centred, over ocean cells within 800 km that the fix owns (pipeline A's nearest-low ownership, from ERA5 lows re-detected with
`hf_structure`'s code, which reproduced the catalog `g800` for 99.98% of fixes). **Centres and ownership are ERA5's for both
systems**, so only the wind field differs; this isolates field noise and ignores any difference in where HRES would place the low.

- `E` = ERA5 (ARCO 0.25 degree) 10 m wind speed maximum over those cells, in kt, from u and v.
- `H` = IFS HRES t0 analysis (0.25 degree) 10 m wind speed maximum over the same cells, in kt, from u and v.
- Track index = maximum over the track's sample fixes.
- Validation before use: on 150 random sample times the ERA5 gust re-derived with this ownership must return `g800` within 0.5 kt
  for at least 99% of fixes; if not, stop and report.

## Crossing rule (frozen)

A track "crosses" in a system when its index is at least that system's cut. The cut is **count-matched**: in each system it is set
so that the number of crossing tracks equals the number of sample tracks with ERA5 gust index at least 71.7 kt (the way pipeline A's
cut was set against the archive). Matching removes the mean offset between systems (HRES is a finer model) and keeps the base rate the
same, so the comparison is about ranking near the cut, not about a bias. Pooled across basins for the primary; by basin as secondary.

## Primary tests (family of 4, BH-FDR)

Statistic: Heidke skill score of the 2x2 table crossing(ERA5) against crossing(HRES) (with equal marginals this equals Cohen's kappa).
Interval: bootstrap over tracks (5,000 resamples, seed 20261008) with the cuts re-matched inside every resample; a season-block
bootstrap over the 6 seasons is reported beside it (few blocks, so wide).

- T1 pooled: H0 HSS <= 0.6 (the agenda's prediction), one-sided bootstrap p.
- T2 Atlantic and T3 Pacific: same null.
- T4 basin difference in HSS (two-sided).

Decision rule on the **pooled HSS and its track-bootstrap 95% interval**:

- **Mostly noise**: upper bound < 0.70 (consistent with the agenda's "no higher than about 0.6").
- **Not noise**: lower bound >= 0.75.
- Anything else: **in between / can't tell**, stated with the number.

Reading: if reproducibility between two analyses is HSS r, no field-based forecast of the crossing in one system can be expected to
beat about r, so r is compared with the forecast skill already reached (0.59). It bounds, but does not apportion, the unexplained
two thirds of PR 85's effect.

## Secondary (descriptive, not part of the FDR family)

- S2 P(HRES crosses | ERA5 index relative to its cut) in bins <0.85, 0.85-0.95, 0.95-1.05, 1.05-1.15, >=1.15.
- S3a ERA5 sustained vs ERA5 gust crossing (same system, different index); S3b HRES sustained vs ERA5 gust crossing.
- S4 HSS among tracks with ERA5 gust index 65-80 kt (a band around the gust cut chosen on the gust index, not on either wind index).
- S5 fix-level agreement with the same cuts.
- S6 linear relation of the track indices (slope, correlation, residual SD, mean HRES minus ERA5) and the ratio of 99th percentiles.
- S7 the 99th percentile of the owned cells instead of the maximum.
- S8 sample restricted to gust index >= 70 kt.

## Power and what this cannot show

Power by normal approximation from the bootstrap SE of the pooled HSS (minimum detectable HSS above 0.6 at 80% power, one-sided 5%),
reported with the result. With roughly 430 positives the SE should be a few hundredths, so a conclusive answer on T1 is expected;
the weak points are not statistical:

1. **The gust is not tested**, only sustained wind (above). If the gust index's own noise (parametrised, model-dependent) is larger,
   this result is a lower bound on the noise.
2. **The systems are not independent.** HRES t0 and ERA5 assimilate largely the same observations and ERA5 uses an older IFS, so
   agreement overstates independence and reproducibility; the noise found here is a lower bound on representation noise.
3. **Resolution**: ERA5 is native 31 km, HRES 9 km regridded to 0.25 degree; count-matching removes the level difference, not
   the spatial one.
4. Centres are ERA5's, so a centre displacement in HRES is not counted. Timing noise (00/12 vs the true peak) is the same in both.
5. Six seasons, 2016-2021; no climate index, so storms are the unit and the independent sample is storms, not seasons.
6. Season 2021-22 was in the gust-threshold calibration window (2021-22 on); nothing is fitted here, but the gust count is matched to it.

## Held-out look

This design fits nothing and scores no model on a separate block. It does use seasons 2016-17 to 2021-22, inside the 2015-2025
block, so it is counted conservatively as a **look at 2015-2025** (the agenda thread numbers it at merge; the entry is in
`research/era5/looks/second_analysis.log`). Zero looks at pre-2001 seasons (decision 1).

## Pull plan

See `results/sizes.txt`. No go-ahead needed if total streaming is under 50 GB. Reduced on the fly; no raw field is written.

## Deviations (post hoc)

1. **Floor set to 65 kt, not 55 (made before any field was read).** The floor rule above was written against the u and v sizes
   alone (`results/sizes.txt`: 46.7 GB at 55 kt, 41.6 GB at 60, 35.6 GB at 65). It left out ERA5 mean sea-level pressure, which the
   extraction needs to re-detect lows and assign ownership (HEAD sizes 7.6 GB at 60 kt, 6.5 GB at 65 kt), and the gust objects for
   the 150 validation times (about 0.5 GB). With those, 60 kt would stream about 49.8 GB, too close to the 50 GB line to rely on an
   estimate, so I take 65 kt: 35.6 + 6.5 + 0.5 = about 42.6 GB. Consequence: tracks with a gust index of 55 to 65 kt are out, so
   the agenda's "55 to 71.7 kt" band is only covered from 65 kt; 944 tracks, 429 with a gust index at or above 71.7 kt. S8 now
   restricts to gust index at or above 70 kt instead of 60.
