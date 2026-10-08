# Greenland barrier and tip jets: which ingredients give more hurricane-force gusts? (pre-registration)

Written 2026-10-08, before any gust, pressure or station value for this question was
read. Everything ERA5 here is a **proxy**: ERA5 instantaneous 10 m gust and surface
fields (ARCO-ERA5, 0.25 degree), not an observation of hurricane-force wind. Pipeline A
(`research/era5/hf_history`) supplies the low-pressure tracks and the 71.7 kt
threshold; the *outcome* below is a regional gust metric defined here, not pipeline A's
index, and is called "G_T" throughout. Anything changed after the first outcome is read
goes in `README.md` under "Departures from this plan", with both versions reported.

## What was looked at before writing this

- `research/era5/highlat/FINDINGS.md` and `PLAN.md` (pipeline A, Atlantic HF-strength
  fixes north of 60N: 31% terrain type, 79% of those from N/NE, 19% in the Cape Farewell
  box). This is why the outcome is built the way it is. `gustloc_fixes.csv` was **not**
  opened for this question.
- Counts of pipeline A low fixes by box and month in `intensity/results/fixes_2004.csv.gz`
  (to size the pull): 4,559 distinct 00/12 UTC times Nov-Mar with a low at 55-67N,
  50-15W. The `g800` column of that file was counted once for sizing in a wider box and
  Oct-Apr (778 fixes at or above 71.7 kt); it was not related to any ingredient.
- Nothing about ingredients versus outcome.

## Question

Is there a set of ingredients that makes hurricane-force gusts from Greenland's
barrier jet or tip jet more likely, for example a strong high over Greenland? Do
ingredients amplify or damp each other beyond their separate effects?

## Sample

- Times: 00 and 12 UTC, **November to March**, seasons 2004-05 to 2025-26 (Oct and Apr
  dropped to keep the pull small; they hold few Greenland events). A season runs 1 June
  to 31 May, so each winter is one season. Decision 1: gust-based, 2004-05 on only.
- Candidate time: at least one pipeline A track fix (the intensity-framework fix table,
  basin atl) with its centre at 55-67N, 50W-15W. This is "a low in the region". The
  reference low at a time is the deepest such fix.
- Pull: gust and MSLP for every candidate time (about 4,560 times x 5.4 MB = about 25
  GB), then SST and 2 m temperature (about 4.4 MB) for the cases and a control sample
  (about 7-10 GB). Total about 32-35 GB, under the 50 GB gate. No existing cache holds
  these fields (the 27 GB of the high-latitude thread were not kept). Only small
  reduced windows are stored, in the ignored `$ERA5_WORK`.

## Outcome

Computed on the 0.25 degree grid at each candidate time.

- **Eligible points**: ocean (land-sea mask < 0.5, as pipeline A), within 300 km of
  Greenland's coast (Greenland only, not Iceland), latitude 58-72N, longitude 52W-10W,
  and at least 400 km from the centre of **every** low detected in that hour's MSLP
  field by pipeline A's own rule (smoothed minima below 1010 hPa, 400 km de-duplication,
  `highlat/gustloc.py: lows()`). The 400 km rule is the high-latitude thread's
  terrain-type test applied point by point: it keeps gust that the storm's own wind
  field does not explain.
- **G_T** = the largest ERA5 gust (kt) over eligible points; 0 if none.
- **Y_T = 1 if G_T >= 71.7 kt** (the pipeline A threshold), else 0. Cases are Y_T = 1
  times; controls are Y_T = 0 times, from the same candidate set, so every control also
  has a low in the region.
- **Class** of a case, by where the largest eligible gust sits: **tip jet** if inside
  59-61.5N, 50W-35W (the Cape Farewell box of the high-latitude thread), **barrier**
  otherwise (east coast and Denmark Strait). This is by position only; wind direction
  at the maximum is not fetched.
- Limits stated now: ERA5 at 0.25 degree smooths steep coastal terrain; the proxy's
  false-alarm ratio is highest here (0.44 north of 60N against 0.26 south, pipeline A);
  the classes are positional proxies for the physical jets.

## Ingredients (fixed, computed at the same hour)

Standardised (z-score) with the discovery half's mean and SD, then applied unchanged to
the held-out half. Month indicators are in every model.

| Code | Ingredient | Definition |
|---|---|---|
| GH | Greenland high | Mean MSLP over the ice sheet, land points with 64-78N, 55W-30W (hPa). MSLP over the ice sheet is a reduction to sea level and is an index, not a pressure that exists; used because it is available in every season. |
| GRAD | Coastal pressure gradient | Mean MSLP over ocean points 0-200 km from Greenland between 62-70N, 40W-10W minus mean MSLP over ocean points 400-600 km from Greenland in the same sector (hPa). Positive = higher pressure at the coast, which drives southward barrier flow. |
| STAB | Cold-air supply | Mean of (SST - 2 m air temperature) over ocean points 60-70N, 40W-10W within 600 km of Greenland (K), the standard cold-air-outbreak measure. Sea-ice points are included at ERA5's frozen SST and are noted as a limit. |
| NAO | NAO phase | CPC daily NAO index on the calendar day of the time. |
| MOT | Low motion | Eastward and northward components (km per h) of the reference low's velocity from the fix table `speed` and `heading`. Two terms, tested jointly. |
| BASE | Low depth and position | Reference low's MSLP (hPa) and its position as km east and north of Cape Farewell (60N, 44W), with the squares and cross product. This is the baseline "comparable lows" adjustment, not an ingredient. A season trend term (years since 2004-05) is also in the baseline. |

Secondary, in its own family: **GBI(Z500)**, the mean 500 hPa height over 60-80N,
80W-20W, from WeatherBench2's 1.5 degree store, available only for times before
2023-01-10, so about 18.3 of the 22 seasons. It is the conventional blocking index and
checks the MSLP-based GH.

## Design

- **Split.** Discovery = seasons with an even start year (2004-05, 2006-07, ... 2024-25).
  Held-out = odd start years (2005-06, ... 2025-26). Primary decisions use this
  assignment only. The swap (fit on odd, test on even) is run afterwards as a
  replication and labelled as one. Pooled all-season estimates are reported and are not
  out of sample.
- **Case-control sampling.** After the gust pass, all cases are kept and controls are a
  seeded (20261008) random sample of 3 per case, drawn within each half, among Y_T = 0
  times. The ingredient fields are fetched for those times only. Logistic regression is
  unchanged by this sampling apart from its intercept; the continuous analysis below
  uses inverse-sampling weights (cases 1, controls 1/f).
- **Models (logistic, Y_T).**
  - M0 = month + BASE.
  - M1 = M0 + GH + GRAD + STAB + NAO + MOT (the additive ingredient model).
  - M1 + one product term for each interaction X1-X4 (below).
- **Interactions (the "amplify or damp" tests)**, each added to M1 on its own as a
  product of standardised terms: X1 GH x depth (MSLP of the reference low), X2 GH x STAB,
  X3 GH x distance of the reference low from Cape Farewell, X4 depth x STAB. Chosen
  for mechanism: a high is more useful if there is a deep low to build a gradient with,
  if the air is cold, and if the low is in the right place.
- **Primary family F1 (9 tests, Benjamini-Hochberg q within the discovery half):** the
  effect of GH, GRAD, STAB, NAO (z, per SD) and MOT (joint Wald, 2 df) in M1, and X1-X4.
  p-values come from a season-block bootstrap (2,000 resamples of the 11 discovery
  seasons, percentile, two-sided), because consecutive 12-hour times are not
  independent. The effective sample is **seasons (11 per half)**; the number of case
  times is reported, and so is the number of case episodes (cases joined when 24 h or
  less apart).
- **Held-out rule, fixed now.** An ingredient is called **supported** if its discovery
  q < 0.10 **and** in the held-out half it has the same sign with a season-block
  bootstrap p < 0.05 (one-sided in the discovery direction). An interaction is called
  supported by the same rule. Anything else is "not supported": reported with its
  interval and, where the interval is narrow enough, as a null with the effect size it
  excludes.
- **Out-of-sample skill.** Fit M0 and M1 on the discovery half; on the held-out half
  report AUC, Brier score, Brier skill of M1 against M0, and against a month-only
  climatology, each with a season-block bootstrap 95% interval (weights for the control
  sampling). Plain-words effect: the probability of Y_T at the 10th and 90th percentile
  of each ingredient with the others at their means.
- **Classes (secondary family F2).** The same M1 for tip-jet cases against controls and
  barrier cases against controls, 5 main-effect tests each (10 tests, own BH). If fewer
  than 40 cases of a class sit in the discovery half, that class is described only and
  not modelled.
- **Continuous outcome (secondary family F3).** G_T (kt) with the M1 ingredients, weighted
  least squares, the same split and bootstrap, 5 tests. A more powerful view of the same
  question; it is not used to rescue a null in F1, and is read only alongside it.
- **Matched sensitivity.** For each case in the held-out half, the nearest 2 controls in
  the sample by standardised distance on reference-low MSLP, latitude and longitude,
  within the same month; conditional logistic M1. Reported beside the adjusted result.
- **GBI(Z500) (family F4, 2 tests):** M1 with GBI in place of GH, then GBI added to M1,
  both on the available times.
- **Station check (family F5).** ISD hourly reports from Prins Christian Sund (04390099999,
  60.0N 43.1W) and Tasiilaq (04360099999, 65.6N 37.6W), held-out seasons only, nearest
  report within +-1 h, quality codes 2, 3, 6, 7 dropped.
  - Contrast A: station mean wind (kt) at case times against control times (does the case
    definition see stronger observed wind?). Tip cases are compared at Prins Christian
    Sund, barrier cases at Tasiilaq, all cases at both.
  - Contrast B: among the sampled held-out times, the top quintile of the ingredient
    score (M1 linear predictor minus its M0 part) against the rest.
  - Each contrast gives the difference in mean wind and in the share of reports at or
    above 34 kt, with a season-block bootstrap interval. 2 stations x 2 metrics x 2
    contrasts = 8 tests, own BH. A coastal station in a fjord does not represent the
    open-water maximum; the check is whether the sign and rough size follow ERA5.
- **Power.** After the gust pass, case counts per half and per class are reported, and
  the smallest odds ratio per SD detectable at 80% power for the discovery half is
  computed from them (normal approximation, case-control with 3 controls per case,
  the design effect from the case episodes). The ingredient fields are fetched only
  after that is written down, and the plan above does not change with the case count
  except the 40-case rule.

## Decision summary (plain language the final answer will use)

For each ingredient: yes, no or cannot tell; how much (odds ratio per SD and the change
in probability from low to high); whether it held out of sample; whether pairs amplify
beyond additive. A well-powered null is a result. Every test run is listed with its q.

## Departures from this plan

None yet.
