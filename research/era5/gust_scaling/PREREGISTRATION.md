# Does the gust follow gradient-wind scaling in HF-strength fixes? (RA-9 first pass, pre-registration)

Written 2026-10-08, before any regression, correlation or truncated fit of gust on depth, size, latitude or an index was
computed. ERA5 **proxy**, **pipeline A** (`research/era5/hf_history`; transitioning tropical cyclones are in, 5-7% of
events). No pull; committed tables only.

**What was looked at before this file was committed.** Column names and one-variable summaries (`describe()`) of
`hf_structure/results/fixes.csv`; that all 2,000 storms of that file are rows of `tele_intensity/results/event_table.csv`
(the join works); the layout of `docs/data/teleconnections.json`. No gust was related to any other variable.

## Question (agenda RA-9)

At fixed depth, how much of the gust is set by storm size, latitude and the index state (ONI, PNA), and does ENSO shift
storm size? Plain-words answer to be given: yes / no (well-powered) / can't tell, and by how much.

## Sample and what it can and cannot show

The 5,983 HF-strength fixes of PR 27 (2,000 storms; Atlantic 3,303 / Pacific 2,680; seasons 2004-05 to 2025-26). Every
fix has pipeline A's 800 km gust index g800 >= 71.7 kt **by construction**. So:

- The outcome is **left-truncated at 71.7 kt**. An ordinary regression of g800 on its drivers in this sample is biased
  towards zero (selection on the outcome). The primary estimator is therefore a **truncated-normal regression on
  ln g800** (truncation point ln 71.7 known, mean linear in the predictors, constant residual SD per basin); OLS is the
  contrast, not the answer.
- Truncation on the outcome alone is handled by that estimator **if** the residual is roughly normal. Selection that works
  through anything else is not. The test is run **inside the HF tail**: it says whether gust scales with depth, size and
  latitude among storms that reach HF. It cannot say whether the scaling holds for the whole cyclone population (weaker
  storms, or storms that never reach 71.7 kt), nor whether an index changes the *probability* of reaching HF. A
  predictor that matters mostly for getting past the threshold is attenuated here. Size is available only at these fixes.
- Size and gust come from the same ERA5 field: the radius of maximum gust (RMG) is read off the gust field whose maximum
  is the outcome. That can build in a correlation. ERA5 at 31 km under-resolves the peak (PR 27 caveat 1).

## Variables, all fixed here

- **y** = ln g800 (kt). Truncation point ln 71.7.
- **Depth** Dp = `ring_storm - msl_fix` in hPa, where `msl_fix` is the fix's central MSLP (`fixes.csv: msl`) and
  `ring_storm` is the 900-1100 km ring-mean MSLP at the storm's deepest fix (`event_table.csv: ring`; one value per storm,
  applied to all of its fixes, an approximation because the surroundings change between fixes). Fixes with Dp < 5 hPa are
  dropped (count reported). x1 = ln Dp. Alternative depths (secondary): Dp_clim = `clim - msl` (month climatology at the
  deepest fix) and Dp_raw = 1013 - msl.
- **Size** R = max(`gmax_r`, 25) km, the radius of the maximum gust from the centre (25 km = grid spacing). x2 = ln R.
  Alternative size (secondary): R_eq = sqrt(`a_g48` / (pi * `own_ocean_frac`)) km, the radius of the equivalent disc of
  ocean with gust >= 48 kt; truncated at the 1,200 km box in the largest storms.
- **Latitude** x3 = ln sin(lat), so that f enters as in the gradient-wind relation.
- **ONI** of the calendar month before the fix month (CPC, `docs/data/teleconnections.json`), **PNA** the mean of daily
  values over days -10..-4 before the fix time (same file). Both standardised over the basin's fixes. The two come from
  distinct regions. NAO is not used (the question names ONI and PNA).
- **Basins** are analysed separately; pooled results are secondary.
- **Unit and inference.** Fixes within a storm and storms within a season are correlated, and ONI has about one value per
  season. Every interval and p value is from a **season-block bootstrap (22 blocks, 2,000 draws)**; p is two-sided from the
  bootstrap normal approximation. Effective n for the index tests is 22 seasons, not 5,983 fixes.

## Predictions (falsifiable)

Gradient-wind balance, V^2/r + fV = (1/rho) dp/dr, with a pressure profile of depth Dp and scale R, gives
V = -fR/2 + sqrt((fR/2)^2 + k Dp/rho). Evaluated over Dp 25-60 hPa, R 150-350 km, latitude 45-60 N (k = 1,
rho = 1.25 kg/m3), the log-log slopes are **d ln V / d ln Dp = +0.55 to +0.72**, **d ln V / d ln R = -0.11 to -0.44**,
and the same for **d ln V / d ln sin(lat)** (R and f enter as the product fR). ERA5 gust is not gradient wind, RMG is not
the pressure-profile scale, and the Dp here is approximate, so the sign is the test and the magnitude is descriptive.

- **H1 (depth).** Slope on x1 is positive, in the order of +0.5 to +0.75. Falsified if the estimate is not positive or
  q >= 0.05.
- **H2 (size).** Slope on x2 is negative, in the order of -0.1 to -0.45: at fixed depth a larger radius gives less gust.
  Falsified if not negative or q >= 0.05.
- **H3 (latitude).** Slope on x3 is negative (same range). Note PR 52 found 1.0 kt less gust per degree at fixed *raw*
  depth; with ring-relative depth in the model that background effect should be gone and only the Coriolis effect should
  remain. If the slope is still negative and of that size, the Coriolis reading holds; if it vanishes, PR 52's slope was
  background pressure.
- **H4 (ENSO and size).** ONI correlates with size: ln R rises with ONI (El Nino, larger storms), Pacific primary,
  Atlantic secondary. This was the agenda's stated prediction; "no ONI size effect" means the ONI gust-versus-depth gap of
  PR 52 / PR 34 is not about size.
- **H5 (index state at fixed depth and size).** ONI and PNA coefficients in the gust model are zero once Dp, R and
  latitude are in. Not directional; reported both ways. If ONI moves ln R (H4) and ln R moves gust (H2), the ONI
  coefficient in the gust model should shrink when x2 is added; reported as the ratio of the ONI slope with and without x2.

## Models and tests

- **M1** (scaling alone): y ~ x1 + x2 + x3. Contrast, not tested.
- **M2** (gust): y ~ x1 + x2 + x3 + ONI + PNA. Truncated-normal MLE, per basin.
- **M3** (size): ln R ~ x1 + x3 + ONI + PNA, OLS, per basin. (Selection on gust acts on this too; noted in the report.)
- **Primary family, 14 tests:** per basin, the five M2 slopes (x1, x2, x3, ONI, PNA) and the two M3 slopes (ONI, PNA).
  BH-FDR across the 14; counts passing reported as "n of 14".
- **Secondary, labelled as such and all reported with their own q:** (S1) OLS instead of truncated; (S2) one fix per storm
  (the deepest, fewest-correlation sample); (S3) Atlantic fixes north of 60 N dropped (no terrain mask; lower confidence),
  and fixes with the maximum gust within 100 km of coast or flagged terrain dropped; (S4) alternative depths; (S5) R_eq as
  the size; (S6) pooled both basins with a basin term; (S7) month-of-year effects added; (S8) stage (deepening / mature /
  filling, PR 27 `stage`) interaction on x2 and x1; (S9) same-time PNA (day 0) as a labelled contrast.
- Sensitivity columns are compared to the primary on sign and size; a sensitivity that flips a primary sign is reported
  in the headline, not in a footnote.

## Decision rules

- "Follows the scaling" for an exponent: sign as predicted, q < 0.05 in the primary family, and the 95% interval
  overlaps the predicted range; "does not follow" if sign is wrong with q < 0.05, or the interval excludes the range with
  the right sign; otherwise "cannot tell".
- A null on H4/H5 is "no" only if the bootstrap minimum detectable effect (2.8 x SE, an approximation for 80% power at
  alpha 0.05, stated as such) is below 0.05 in ln R per SD of ONI (a 5% change in radius) or below 0.01 in ln g per SD
  for the gust model; otherwise "cannot tell".
- No tuning after the fact; any change to variables or the model is entered below as post hoc and both versions are shown.

## Held-out seasons

No hold-out split is used: this is an estimation on all 22 seasons, with no model frozen or scored. It spends **no look**
at the 2015-25 block in the sense of `heldout_looks.log` (the 2015-25 seasons enter the sample, as they entered PR 27).

## Power plan

Minimum detectable effects from the bootstrap standard errors (2.8 x SE). Not a planted-effect simulation; the
approximation is labelled where quoted.

## Deviations (post hoc)

Logged 2026-10-08, after the committed plan and before the full 2,000-draw run. Estimates from smoke runs (20-40 draws)
had been seen when these were written, so none of these is blind to the signs and rough sizes of the primary slopes.

1. **Optimiser, not estimator (no change to the plan).** The first implementation of the pre-registered truncated-normal
   MLE (joint BFGS on uncentred covariates) stopped early on a flat ridge and gave start-dependent slopes. Fixed by
   centring the covariates and profiling sigma (log-concave location problem solved by damped Newton at each sigma,
   bounded 1-D search over sigma); a profile-likelihood table over sigma 0.07-0.5 and six Nelder-Mead starts agree on
   one maximum per basin (Atlantic sigma 0.128, Pacific 0.124). The estimator is the one in the plan.
2. **Added secondary S10 (post hoc): generalised Pareto regression of the exceedance g800 - 71.7 kt**, log scale linear
   in the same covariates, common shape. Reason: the truncated-normal slopes depend on the normal tail being right
   (the fitted slope scales with sigma along the profile), so a second tail model is a check. Its slopes are on the
   log-scale of the exceedance; they are converted to the elasticity of the mean gust by x mean(exceedance)/mean(gust)
   (about 0.09). Reported next to the primary, not instead of it.
3. **Added diagnostic (post hoc): probability-integral-transform check** of the truncated-normal and GPD fits against the
   observed gusts (KS distance), to say which tail describes the data.
4. **Implementation notes.** Alternative depths use the same Dp >= 5 hPa rule. S9 drops fixes without a day-0 PNA value.
   49 of 5,983 fixes have Dp < 5 hPa and are dropped from every model.
