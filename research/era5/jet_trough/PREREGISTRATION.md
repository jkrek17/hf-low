# Pre-registration: jet speed, upstream trough depth, and bomb / hurricane-force odds

Written and committed before any outcome (bomb, HF) was joined to the jet or trough
predictors. ERA5 proxy, pipeline A tracks (`research/era5/hf_history`), the
near-storm framework's fix table (`research/era5/intensity/results/fixes_2004.csv.gz`,
`env_2004.csv.gz`). "HF" is pipeline A's gust index reaching 71.7 kt. It is not the archive.

## What has been seen

- The framework's base rates and skill table (README): 3.4% of fixes reach HF within 24 h,
  3.4% of class-eligible fixes are rapid deepeners, and the full model's skill.
- Not read before this file: `coefficients.csv`, any cross-tabulation of an outcome
  against `jet250`, and any trough variable (not yet extracted).
- The predictors defined here are computed first, from ERA5 only, and checked for coverage
  and distribution without outcomes (step 1 below).

## Question (Jason, 2026-10-08)

Is there a threshold jet speed and upper-level trough depth that tilts a cyclone toward being a
bomb (rapid deepener) or reaching hurricane force? Is there a probability of HF?

The second part is answered by the framework (hf-low PR 12). Here we report what it gives
conditional on jet speed and trough depth, and whether trough depth adds skill to it.

## Hypotheses

- **H1 (jet).** The odds of bomb and of HF onset increase with the jet speed near the low,
  at fixed storm state. The question is whether the curve has a threshold (a break) or is a smooth ramp.
- **H2 (trough).** Same for upstream trough depth.
- **H3 (combination).** Jet speed and trough depth amplify each other beyond their separate
  effects (positive interaction on the logit scale).
- **H4 (added skill).** Adding trough depth to the framework's full predictor set improves its
  held-out skill.
- **Null for H1, H2.** A smooth ramp: a straight line or smooth curve in the focal variable
  fits as well as a piecewise linear fit with one knot, on held-out seasons.

Direction is pre-specified as positive for H1-H3. A negative sign is reported, not reinterpreted.

## Definitions

**Time.** Predictors are measured at the fix time t (00 or 12 UTC). Outcomes cover the next 24 h,
(t, t+24 h]. Lead is 0-24 h, mean about 12 h before the outcome is reached. The brief's "12 to 24 h
before" is therefore approximated, not met exactly; a lead-24 sensitivity below shifts the outcome
window to (t+24, t+48 h].

**Population.** The framework's 159,430 fixes (every pipeline A low fix below 1010 hPa, in domain,
at 00/12 UTC, seasons 2004-05 to 2025-26). Per basin (Atlantic, Pacific); never pooled.

**Outcomes.**
- `BOMB`: `ndr24` >= +1.0 Bergeron (24 hPa in 24 h scaled by sin(lat)/sin(60 deg), the
  Sanders-Gyakum criterion), the framework's rapid-deepening class. Defined on fixes whose track
  survives 24 h with |ndr24| <= 3 (the framework's class-eligible fixes).
- `HFON`: the gust index reaches 71.7 kt in (t, t+24 h], among fixes that are not HF at t
  (`hf_now` false). A track that ends inside the window counts as not reaching it.
- `BOTH/neither` table: on fixes that are class-eligible and not HF at t, the four cells
  bomb only, HF only, both, neither. Reported as a table, not tested.

**Jet speed (primary).** `jet250`, already in `env_2004.csv.gz`: maximum 250 hPa wind speed within
1000 km of the fix (kt), 1.5 deg ERA5. Reused as is.

**Upstream trough depth (primary), `trough_up`.** On the 5.625 deg ERA5 grid (WeatherBench2's
64 x 32 conservative regrid, and the same regrid of ARCO-ERA5 after 2023-01-09), 500 hPa height
anomaly (m) = Z500 minus the calendar-month mean of Z500 at that cell over the fit seasons
2004-05 to 2014-15 (all 00 and 12 UTC times of those months). `trough_up` = minus the lowest anomaly
among grid cells 500-2500 km from the fix at bearings 225 deg to 315 deg (west, upstream in
westerlies). Larger = deeper trough. Cells beyond the grid (south of 14 N) are ignored; a fix with
no cell in the sector gets a missing value and is dropped from that analysis (count reported).

**Secondary variables (separate FDR family, labelled).**
- `trough_loc`: minus the mean Z500 anomaly within 1000 km of the fix. Partly the storm's own
  signature, so it is not a clean predictor of the storm's later intensity.
- `jet_up`: maximum 250 hPa wind speed (kt) in the same upstream sector, 5.625 deg grid.

**Storm-state covariates `C`** (always in the adjusted models): `msl`, `dp12` (0 if the track is
under 12 h old, plus the `young` flag), latitude, cos and sin of day of year. For HFON add the
current gust index `g800`. Standardised with the fit-season median and 0.5-99.5 percentile clip,
as the framework does.

## Design

**Split.** Fit seasons 2004-05 to 2014-15 (11). Held-out seasons 2015-16 to 2025-26 (11). Every
number that is called a result comes from the held-out seasons or from a model fitted on the fit
seasons only. The held-out seasons are scored once per model; every model is declared here.

**Raw curves (descriptive).** Observed frequency of each outcome by decile of the focal variable,
fit and held-out seasons separately, per basin, with 90% intervals from resampling seasons.

**Adjusted models.** For each basin x outcome x focal variable `x` (standardised on fit seasons),
unpenalised logistic regressions with `C` plus:
- `M_C`: nothing (covariates only);
- `M_lin`: linear in `x`;
- `M_spl`: natural cubic spline in `x`, 4 degrees of freedom, boundary knots at the 1st and 99th fit percentiles
  (the smooth alternative);
- `M_hinge`: continuous piecewise linear, `x` plus max(x - k, 0), with one knot `k` chosen by maximum
  likelihood on the fit seasons over a grid from the 10th to the 90th fit percentile in steps of 2.5 percentiles.

**Tests (held-out log-loss, summed by season, 2,000 resamples of the 11 held-out seasons).**
- T1: `M_lin` beats `M_C` (does `x` carry anything beyond storm state?).
- T2: `M_hinge` beats `M_lin`.
- T3: `M_hinge` beats `M_spl`.
Each p is the share of resamples in which the held-out log-loss gain is <= 0. A gain is the difference in
total log-loss; a positive value favours the first-named model.

**Rule for the plain-words verdict, per basin x outcome x variable (decided now):**
- **Threshold supported**: T1, T2 and T3 all pass FDR (q < 0.05), the hinge slope after the knot is
  larger than before (same sign as H1), and the knot lies inside 15th-85th fit percentile.
- **Smooth ramp (no threshold)**: T1 passes and T2 or T3 does not pass. If `M_spl` beats
  `M_hinge` with q < 0.05, curvature is real but smooth.
- **No relation beyond storm state**: T1 does not pass and the power analysis says an effect of the
  size below was detectable.
- **Cannot tell**: otherwise.

**Sharpness (reported for any supported or smooth case).** From the fit-season model at median
covariates: predicted probability at the 10th, 25th, 50th, 75th, 90th percentile of `x`; the odds ratio
for `x` going from 0.5 SD below the knot to 0.5 SD above it, and the same two points in the
spline model. The knot location has a 90% interval from 500 resamples of the fit seasons.

**H3 interaction.** On the primary variables, `M_add` (spline in jet, spline in trough, plus `C`) against
`M_int` (adds the product of the two standardised variables): held-out log-loss gain, same bootstrap.
Also a 3 x 3 grid of fit-season terciles of jet x trough giving observed frequencies of BOMB, HFON
and the four-cell table, with counts and held-out-season intervals.

**H4 added skill.** The framework's `full` predictor set (`model.SETS["full"]`, L2 C = 1, unchanged)
fitted on the fit seasons, versus the same plus `trough_up` (and `trough_loc` and `jet_up` in a
second variant), scored on held-out seasons: Brier skill score against basin-month
climatology and the yes/no HSS at the count-matched cut for HFON-like `hf24` and for BOMB, pooled across
basins as the framework does with a basin term, plus per basin. Intervals from resampling held-out seasons.
The framework's own mean predicted P(HF24) is also tabulated in the jet x trough terciles next to the
observed frequency.

**Multiplicity.** One Benjamini-Hochberg family over the primary tests: 2 basins x 2 outcomes x 2
primary variables x 3 tests (24), 4 interaction tests, and 4 added-skill tests (basin x outcome): 32
tests. Secondary variables are a second family (2 variables x 2 x 2 x 3 = 24) with their own q. Every
test is reported, passed or not.

**Sensitivities (reported, not in the family).**
1. Lead 24: outcome window (t+24, t+48 h].
2. Oct-Apr fixes only (the HF season), against all months.
3. Unadjusted hinge, `x` only.
4. Zonal-eddy height (Z500 minus its zonal mean at the same latitude and time) in place of the
   climatological anomaly, since the warming trend of Z500 over 2004-2025 shifts the anomaly baseline.

**Power (stated now, run before reading the T1-T3 results).** Simulation on the real fit and held-out fixes:
outcomes are redrawn from each basin's fitted `M_lin`, plus a hinge with the post-knot slope raised
so that the odds ratio across 0.5 SD either side of a knot at the median is 2.0 on top of the linear
trend. The share of 200 simulations in which the T2 and T3 rule above declares a threshold is the power;
the share under no added hinge is the false-positive rate. The "no relation" verdict needs T1 power of
at least 80% for a per-SD odds ratio of 1.15.

## Rules for changes

Any change after the held-out seasons have been scored is logged below as post hoc, with both versions
reported. The held-out look count is kept in `results/looks.log`.

## Deviations log

(empty at commit)
