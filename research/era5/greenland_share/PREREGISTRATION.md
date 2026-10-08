# Why does a strong Greenland high lower the Atlantic HF share? (written plan, EXPLORATORY)

Written 2026-10-08 before any storm characteristic or HF flag was related to the Greenland-high index in this
analysis. Everything ERA5 is a **proxy**: **pipeline A** (`research/era5/hf_history`, 800 km ocean gust index,
HF-equivalent at 71.7 kt). Transitioning tropical cyclones are in, as in PR 14 and PR 63.

**Status: exploratory.** All 22 gust-era seasons (2004-05 to 2025-26, Decision 1) were already used by hf-low
PR 58 and PR 63, and PR 63 reported the share effect (0.870 per SD of the high) that this study explains. There
is no held-out sample. No result here is confirmatory; it is written down in advance so that the comparisons are
fixed, all of them are reported, and anything changed afterwards is logged as post hoc.

## What was looked at before writing this

- Published results: PR 63 (`nao_share_barrier/README.md`: share ratio 0.870 (0.82-0.92) per SD of the lagged
  Greenland high GH; NAO share ratio 1.125 alone and 1.039 with GH in the model; r(NAO, GH) = -0.69), PR 58
  (`greenland_jets/README.md`: GH raises coastal-band gust odds about 4 times per SD at fixed low depth and
  position), PR 14, PR 11.
- Column names, row counts and join coverage of `hf_history/results/all_tracks.csv.gz`,
  `intensity/results/fixes_2004.csv.gz` and `env_2004.csv.gz` (16,992 Atlantic tracks 2004-05 on; 99.7% have fix
  rows). No outcome was related to GH, NAO or any other predictor, and nothing was fitted.
- GH, NAO and the sample are those of PR 63 (same code path: `freq_split/split.py`, `nao_share_barrier/gh_daily.csv`),
  so the PR 63 numbers are reproduced first (RR(share) NAO 1.125, GH 0.870, NAO|GH 1.039). If they do not
  reproduce, stop and explain before going on.

## Question and the three explanations

PR 63: a strong lagged Greenland high **lowers** the share of Atlantic cyclones that reach HF (0.870 per SD), also
with barrier-type fixes removed. PR 58: the same high **raises** the odds of a barrier-jet gust about 4 times at
fixed low depth and position. Which explains the share effect?

| Code | Explanation | Signature the data must show |
|---|---|---|
| **STEER** | The high steers cyclones south or east, away from the places where they become HF | cyclone positions shift with GH (NAO held fixed); adding position to the share model removes most of the GH share effect |
| **WEAKEN** | The high makes cyclones shallower (blocking, less deepening, weaker jet) | lower peak intensity or deepening rate with GH; adding depth/deepening/environment removes the effect, and at fixed depth the GH coefficient on HF is **zero or positive** (the PR 58 sign, which reconciles the two papers) |
| **STAND-IN** | The high is a proxy for the NAO phase | GH loses its share effect when NAO is in the model, and NAO keeps its own |

These are not exclusive. The plain-words answer to give: which of the three the data favour, or that they cannot
be told apart, with the interval that supports it.

## Sample, indices, model (PR 63 / PR 14 primary)

Atlantic pipeline A tracks (lows below 1010 hPa, at least 24 h), genesis day = first fix, window 1 October plus 210
days, seasons 2004-05 to 2025-26 (22). NAO and GH are daily series averaged over days -10..-4 before the genesis
day, standardised over the analysis days (PR 63's `NAO_lag`, `GH_lag`); the effects below are **per +1 SD**. HF =
track peak 800 km gust >= 71.7 kt (1,003 of 9,636 tracks in PR 63's window). Calendar-month effects and a linear
season trend are in every model. Resampling is by **season block** (22 blocks, 2,000 draws, seed 20261008); the
effective n is 22 seasons times about 5 independent index values per season, never the 9,636 cyclones.
Storm-level models use the same month/trend terms; fix-level quantities are collapsed to one value per storm
first, so fixes are never treated as independent.

## Block A: how much can the two indices be separated (A1-A4)

- **A1** correlation of the two lagged indices over analysis days, over storm genesis days, and over the 22
  season means; variance inflation factor and the partial SD of GH given NAO (sqrt(1 - r^2)).
- **A2** log share coefficient, by Poisson daily-count models as in PR 63 (share = HF coefficient minus all-cyclone
  coefficient): GH alone, GH given NAO, NAO alone, NAO given GH. Tests A2a-A2d. Also the same for the HF count and
  the all-cyclone count separately (descriptive, so that "fewer HF" and "more cyclones" can be told apart).
- **A3** joint season-block bootstrap of (GH|NAO, NAO|GH) share coefficients: the 95% region, and whether it
  includes the three corners "GH only", "NAO only", "neither".
- **A4** minimum detectable effect (2.8 x bootstrap se) of each coefficient alone and with the other in the model:
  how large a GH effect independent of NAO this sample could have shown.

Rule for "can be separated": the GH|NAO share coefficient excludes 0 **and** the NAO|GH coefficient interval
excludes the PR 63 NAO-alone value; "cannot be told apart" if the GH|NAO interval includes both 0 and the GH-alone
value.

## Block B: what the high does to the storms (B1-B12), NAO held fixed

Storm-level outcomes (one value per Atlantic storm, from `fixes_2004.csv.gz` and `all_tracks.csv.gz`; the first
row of a track in the fix table is "first fix below 1010 hPa at 00/12 UTC", called **first position**). Each is
regressed on GH and NAO jointly (OLS or logistic as stated), with month and trend terms; the test is the GH coefficient
**given NAO**. GH alone and NAO given GH are reported as contrasts.

| Code | Outcome | Reads as |
|---|---|---|
| B1, B2 | first-position latitude, longitude | where storms start |
| B3, B4 | latitude, longitude of the minimum-pressure fix | where storms mature |
| B5, B6 | latitude, longitude of the peak-gust position (`peak_lat`, `peak_lon`) | where the wind is strongest |
| B7 | track minimum pressure (hPa) | intensity |
| B8 | track peak 800 km gust (kt) | wind intensity |
| B9 | maximum 24 h normalised deepening rate (`ndr24`), tracks with at least one value | deepening |
| B10 | storm passes through the deep-low box 55-67N, 50-15W (PR 58's region) at its minimum-pressure fix (logistic) | steering into/away from the deep-low region |
| B11 | `jet250` at the first position (`env_2004`) | upstream jet |
| B12 | `eady` at the first position | baroclinicity |

Composites (descriptive, no tests): cell means and counts of the storm quantities above by GH tercile within NAO
tercile (3x3, with the sparse corners flagged), and a map of fix density for the top minus bottom GH tercile.

## Block C: how much of the GH share effect do the candidate mediators carry (C1-C5)

Track-level logistic HF ~ GH + NAO + month + trend (+ mediators), season-block bootstrap. Quantity: the GH
coefficient (log odds per SD) unadjusted, then with mediators added, and the **fraction retained** (adjusted /
unadjusted) with its interval. Mediator sets:

- **C1 STEER:** first-position latitude and longitude (linear plus squares), nothing post-genesis.
- **C2 WEAKEN-depth:** C1 plus track minimum pressure and peak-pressure-fall rate (max `ndr24`).
- **C3 WEAKEN-environment:** C1 plus first-position `jet250`, `eady`, `sstgrad`, `tcwv`.
- **C4 all:** C1 + C2 + C3.
- **C5 PR 58 sign:** in C2, is the GH coefficient on HF at or above 0 (one-sided p, bootstrap)? Prediction under
  WEAKEN: yes. Minimum pressure and HF are close (gust is set largely by pressure gradient), so C2 may over-explain by
  construction; the sign of the residual GH effect is the informative part, not the retained fraction.

The same models with NAO removed from the right-hand side are a contrast, not a test.

## Tests, families, rules

- Family A: A2a-A2d (4). Family B: B1-B12 (GH given NAO; 12). Family C: C1-C5 (the paired difference of the GH
  coefficient with versus without the mediators, and C5; 5). BH q within family and across all 21. Count passes out of 21.
- No test is "confirmed". With q < 0.05 across all 21 an item is "supported in this sample"; q between 0.05 and
  0.20, a lead. A null needs the minimum detectable effect stated; otherwise "inconclusive".
- Verdict mapping (stated now, applied after): STEER favoured if C1 retains under 50% (upper interval bound under 75%)
  and at least one of B1-B6, B10 passes q < 0.05; WEAKEN favoured if C2 or C3 retains under 50% and at least one of
  B7-B9, B11, B12 passes and C5 holds; STAND-IN favoured if the GH|NAO share coefficient retains under 50% of GH alone
  with an interval including 0 while NAO|GH keeps at least 75% of NAO alone; "cannot be told apart" when A3's region
  includes both corners, or two verdict conditions hold together. More than one can hold: say so.
- What this cannot show: causation (the high is partly storm-made, which the 4-10 day lag lowers, not removes);
  whether proxy error near Greenland (no terrain mask, Atlantic fixes north of 60N lower confidence) matters; anything
  about seasons outside 2004-05 to 2025-26. Mediator timing is post-genesis, so "mediation" is statistical
  association, not a causal chain.

## Deviations (post hoc)

None yet.
