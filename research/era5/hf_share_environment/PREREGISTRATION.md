# Pre-registration: which environmental ingredients track the midwinter rise in the hurricane-force share

Written and committed before any environment value was compared with any outcome. ERA5 proxy throughout.
Cyclone tracks, events and the share are **pipeline A** (`research/era5/hf_history`). The environment values are the near-storm
fields of the intensity framework (`research/era5/intensity/results/env_2004.csv.gz`: ERA5 1.5 degree, WeatherBench2 regrid up to
2023-01-09 and coarsened ARCO-ERA5 after; no new ERA5 pull).

## Question

The seasonal-cycle thread (hf-low PR 45) found that the HF peak is set by the **share** of cyclones that reach HF, not by their number
(Pacific about 7-8% in October-November to 11.4% in December-January; Atlantic about 7% to 15-16% in December-February), and that
jet speed, baroclinicity and the like were not examined. Which ingredients of a cyclone's environment, measured at the time and place
it first deepens, rise from autumn to midwinter, predict which cyclones reach HF within a period, and so account for the rise in the share?

## What had been seen before this file

Column names, row counts and overall missing fractions of the intensity tables (`sst` and `sst_t500` 11%, `sstgrad` 20%, others 0).
The PR 45 monthly shares and latitude result. Nothing about any ingredient against month or against HF.

## Data and unit

- Unit: one pipeline A track in `all_tracks.csv.gz` (a low below 1010 hPa for at least 24 h with at least two in-domain fixes),
  seasons 2004-05 to 2025-26 (22), genesis month (month of `start`) October to April, analysed separately for each basin.
- Outcome HF: the track's `gust800_kt` >= 71.7 kt (the catalog event rule; a proxy). Transitioning tropical cyclones are left in, as in PR 45.
- Periods: ON = October-November, DJF = December-February, MA = March-April. The primary contrast is **DJF versus ON** (the rise);
  DJF versus MA (the fall) is secondary.
- Real sample size: 22 seasons. All intervals and p values are 2,000 season-block bootstrap draws (seed 20261010), resampling whole seasons.

## Reference fix (primary): the first deepening fix

The environment is read at one fix per cyclone, chosen using only the past of the track so that it is the same rule for every cyclone:
the first 00/12 UTC in-domain fix whose 12 h pressure change `dp12` <= -3.6 hPa (the lower edge of the deepening class, 0.3 Bergeron,
at 60N). A cyclone that never reaches it ("non-deepening") uses its first 00/12 in-domain fix and stays in every denominator, so the share
is the same HF/cyclones as in PR 45. The -3.6 hPa rule was set before its frequency was looked at.

## Ingredients (primary set, 8) and expected direction

All are columns of the env table, z-scored within basin over the included cyclones, so effects are per SD.

| code | ingredient | sign expected to help HF | expected seasonal change |
|---|---|---|---|
| jet250 | 250 hPa maximum wind within 1000 km | + | up in DJF |
| eady | 850-500 hPa Eady growth rate, 500 km mean | + | up in DJF |
| sstgrad | maximum SST gradient within 500 km | + | up in DJF |
| sst_t500 | SST minus 500 hPa temperature (air-sea instability) | + | up in DJF |
| flux | upward sensible plus latent heat flux | + | up in DJF |
| tcwv | total column water vapour | + | down in DJF |
| div300 | 300 hPa divergence | + | not stated |
| vadv500 | 500 hPa absolute vorticity advection | + | not stated |

H1: jet250 and eady are the leading candidates (they are the baroclinic supply). H2: tcwv rises into autumn and falls into winter, so it
cannot carry the rise even if it matters within a period. Hypotheses are directional only to define "helps"; tests are two-sided.
Secondary ingredients, not in the primary family: latitude of the reference fix, Hart B, -V_T lower, -V_T upper, SST.
Jet position and storm position relative to the jet are not in the table; the reference fix latitude is the only position measure.
Not examined for that reason.

## Quantities, per basin and ingredient X

1. **Alignment** Delta_X: mean of X (SD) in DJF minus ON. Tests: Delta_X != 0.
2. **Within-period association** beta_X: log-odds of HF per SD of X from `logit P(HF) = a + d_DJF + d_MA + beta X` (period fixed effects, so only
   storm-to-storm variation inside a period). Tests: beta_X != 0.
3. **Explained fraction** F_X = 1 - c1/c0, where c0 is the DJF coefficient without X and c1 with X (period effects, ON baseline). Joint F with all 8
   ingredients, and the drop in F when each is removed from the joint model (unique part). Reported with 95% percentile intervals, not tested.
4. **Two stages** (descriptive): the share = P(deepening fix exists) x P(HF given a deepening fix). Log-odds change DJF minus ON of each factor, with intervals.

Models are unpenalised logistic regressions on complete cases for the ingredient(s) in the model. Missing rate of each ingredient by period and its HF rate
are reported; a model on the common complete-case sample is a sensitivity.

## Family and verdict rules

Primary BH family: Delta_X and beta_X tests, 8 ingredients x 2 basins x 2 = 32 two-sided bootstrap p values, one q across all 32.
An ingredient **lines up with the rise** in a basin when Delta_X and beta_X both have q < 0.05, the product beta_X * Delta_X is positive, and the 95%
interval of F_X is above 0. It **does not line up** when the F_X interval lies entirely below 0.15, or when either q >= 0.05 with a well-powered interval
(beta_X interval within +-0.10 log-odds per SD); otherwise **cannot tell**. A joint F whose interval includes 1 means the ingredients may account for the whole rise; below 0.5
upper bound means a large part of the rise is not in these ingredients. A well-powered null is a full result.

## Secondary analyses (each its own BH family; not pooled with the primary)

- S1: entry environment, the first 00/12 in-domain fix of every track (earlier than the deepening).
- S2: the fix of greatest 24 h deepening (uses the future of the track, so it is partly a response of the environment to the storm; labelled).
- S3: fix-level onset model on `fixes_2004`: every 00/12 fix not already HF, outcome HF within 24 h, month of the fix, period fixed effects, same F.
- S4: the fall contrast (DJF versus MA) for the primary reference fix.
- S5: secondary ingredients.
- S6: common complete-case sample; excluding non-deepening cyclones.
- S7: does beta_X differ between ON and DJF (X by period interaction)? Reported as a descriptive check on the tail argument in PR 45.

## Limits stated in advance

- The environment is read at the time of deepening, close to the low, so part of it (moisture, flux, divergence) may be the storm's own doing, not a
  cause. S1 is the earlier reading.
- Period contrasts are three coarse bins; the correlation of two months' environments inside a season is handled by the season block.
- ERA5 proxy: nothing here is the archive. The Pacific and Atlantic HF share are pipeline A's.
- Departures from this plan are logged below as post hoc, with both versions reported.

## Departures log

(empty at commit)
