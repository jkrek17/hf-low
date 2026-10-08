# Greenland barrier and tip jets: which ingredients give more hurricane-force gusts?

ERA5 **proxy** throughout: ARCO-ERA5 0.25 degree gust, MSLP, SST and 2 m temperature; WeatherBench2 1.5 degree
Z500 for one check; pipeline A's (`research/era5/hf_history`) low tracks. The outcome **G_T** is a regional gust
metric defined here, not pipeline A's index. Nov-Mar, seasons 2004-05 to 2025-26 (decision 1: gust-based, 2004-05
on). The plan is `PREREGISTRATION.md`, committed first (`438172f`); stage 1 and the case counts were committed
(`0c5b458`) before any ingredient was related to the outcome.

## Answer in plain words

- **The low itself is most of it.** A deep low (reference-low MSLP of 972 hPa or lower, the deepest third) over the Irminger Sea and Denmark Strait
  (the composite pressure minimum at cases is near 63N 32W) holds 214 of the 269 cases. Depth and position alone separate the times with a
  hurricane-force Greenland-coast gust from the rest well (held-out AUC 0.81).
- **A strong high over Greenland adds a lot on top: yes.** Per SD (11 hPa) of ice-sheet MSLP the odds of a
  gust of 71.7 kt or more in the coastal band rise about 4 times (odds ratio 3.7 with the even seasons, 4.3 with
  the odd seasons; 4.0 and 3.4 in the swap). From the 10th to the 90th percentile of the index the chance goes from
  3.7% to 19.9% of low-in-region times, with the low's depth and position held to the same distribution. It
  held out of sample, replicated with the seasons swapped, and was confirmed by an independent index (Z500
  blocking index, OR 2.7 per SD) and by tight matching on all 4,559 times (post hoc, OR 4.6-8.9).
  The pre-registered matched check did not show it and is explained below.
- **Cold air over the sea: probably a modest help, not confirmed.** (SST minus 2 m temperature upstream.)
  Odds ratio 1.4-1.9 per SD in all four fits, but it passed the held-out rule in neither direction of the swap.
  Pooled 1.54 [1.10, 2.18]; q 0.065. Cannot tell more than that with 11 seasons per half.
- **The pressure gradient across the strait: no.** Once the high is in, the coast-versus-offshore gradient adds
  nothing (OR 0.97 [0.64, 1.39] pooled).
- **NAO phase: no in the pre-registered test** (0.98 [0.79, 1.18] per SD pooled, excluding more than about 20%
  per SD). A post hoc fit on all times without the cold-air term gave 1.23 [1.05, 1.47]; the two can both be right
  because NAO and the high are correlated (-0.59), so read NAO as small and unresolved.
- **The low's motion: cannot tell.** Direction changed between the adjusted and the matched fits.
- **Do they amplify each other?** Only one pair shows anything. The high matters more when the low is not at its
  deepest (high x depth interaction +0.19 per SD^2 on the logit scale, same sign in all four fits; q 0.044 in the
  fit half, one-sided held-out p 0.067, so **not supported** by the rule fixed in advance; pooled OR 1.19 [1.02,
  1.43], q 0.065). In words: a strong high and a very deep low partly substitute for each other, so together
  they are a little less than the product of their separate odds. The other three pairs are null or flip sign.
- **Tip jet against barrier jet.** 69 of the 269 cases peak in the Cape Farewell box (tip), 200 elsewhere
  (barrier). The barrier class carries the high effect (OR 4.1 even, 6.3 odd). The tip class has 29 cases in the
  even seasons, under the 40-case rule, so it is described only: tip cases have a lower ice-sheet pressure
  (-0.69 SD) and a more positive NAO (+0.57 SD) than controls. Do not read the tip jet as an effect of the
  Greenland high.
- **Observations agree with the case definition, but not clearly with the ingredient ranking.** Prins
  Christian Sund (Cape Farewell) reports 25.1 kt at case times against 14.6 kt at control times, and 38.1 kt at tip-case
  times (47.8% of reports at or above 34 kt against 4.3%). The ingredient-favoured top quintile does **not** show
  stronger wind at either station (-0.5 kt at Prins Christian Sund, -1.8 at Tasiilaq). A post hoc regression gives
  +2.9 kt per SD of the high at Prins Christian Sund [-0.3, +5.5] and +1.2 at Tasiilaq [+0.1, +2.6].
  The ERA5 high effect is therefore **not** independently confirmed by station wind.

Everything is a statement about ERA5 gusts at 0.25 degree near steep terrain, where the proxy's false-alarm ratio
is highest (pipeline A 0.44 north of 60N against 0.26 south). No station sees a hurricane-force wind here, so the
validation is of sign and rough size only.

## Numbers (ERA5 proxy; 22 seasons; 11 per half; season-block bootstraps, 2,000 resamples)

| Quantity | Value |
|---|---|
| Candidate times (a pipeline A low at 55-67N, 50-15W, 00/12 UTC, Nov-Mar) | 4,559 |
| Cases, G_T >= 71.7 kt (5.9%): tip / barrier | 269: 69 / 200 |
| Cases, even seasons (fit) / odd seasons (held out) | 127 (87 episodes) / 142 (85 episodes) |
| Smallest OR per SD detectable at 80% power, fit half | 1.33 (case times), 1.41 (episodes) |
| GH odds ratio per SD, fit / test / swap fit / swap test | 3.69 [2.21, 11.06] / 4.34 [2.18, 11.53] / 3.99 [2.04, 10.32] / 3.42 [2.09, 9.92] |
| P(G_T >= 71.7 kt), GH at 10th to 90th percentile | 3.7% to 19.9% (x5.35) |
| STAB odds ratio per SD, same four | 1.40 [0.72, 2.90] / 1.86 [1.27, 2.86] / 1.76 [1.27, 2.62] / 1.36 [0.75, 2.61] |
| GRAD, NAO pooled OR per SD | 0.97 [0.64, 1.39]; 0.98 [0.79, 1.18] |
| Held-out AUC, M1 / M0 (depth and position only) | 0.864 [0.828, 0.896] / 0.811 [0.758, 0.849] |
| Held-out Brier skill, M1 vs M0 / vs month climatology | +0.085 [-0.015, +0.185] / +0.197 [+0.068, +0.295] |
| Z500 blocking index replaces GH, OR per SD, fit / test | 2.79 [1.33, 9.60] / 2.66 [1.25, 7.48] |
| Continuous G_T, GH, kt per SD, fit / test | +4.52 [+3.23, +5.84] / +4.46 [+2.13, +6.69] |
| Pull | 35.5 GB (24.4 stage 1, 11.1 stage 2), under the 50 GB gate |

Full tables, every test and its q: `results/results.txt` (primary), `results/results_swap.txt` (replication),
`results/tests.csv`. Composite pressure map: `results/composite.png` (deep low near 63N 32W; ice-sheet box mean
-3.2 hPa at all cases against same-month non-case times, which reflects the low's own field, so the high effect
is the effect at fixed depth and position, not a marginal difference).

## What the pre-registered rule says, test by test

F1 (9 tests): **GH supported** (q 0.005 in the fit half; same sign and one-sided p 0.001 held out; same in the
swap). STAB, GRAD, NAO, MOT and the four interactions: not supported. Counting the swap as a replication, GH is
the only ingredient that held in both directions.

## Matched sensitivity: pre-registered, and it failed to match

The plan's matched check used the 3:1 control sample, which is too thin to match on. The cases are 13.5 hPa deeper
than their matched controls (held-out) and 15.7 hPa (fit), so it left residual confounding, and its null for GH
(OR 0.90 and 1.02) cannot be read as a result. `posthoc2.py` repeats the match on all 4,559 times (GH, GRAD and NAO
are known at every time): with up to 5 controls and no caliper the OR is 2.09 [1.38, 3.11] (cases still 6.3 hPa
deeper); with a 0.5 SD caliper 4.63 [2.90, 7.21] (-1.0 hPa); with a 0.25 SD caliper 8.86 [4.41, 27.96]
(-0.2 hPa; 127 sets). Both versions are reported. Tighter matching moves the estimate up, not down.

## Clarifications fixed before the first model was fitted

- The fix table's `speed` is in kt over the previous 6 h (the plan said km per h) and is missing at a track's first fix;
  those 57 rows keep zero motion and get an indicator in the baseline.
- z-scores use the fit half's mean and SD weighted by the inverse sampling weight, so "per SD" is per SD of all
  low-in-region times.
- Logistic fits use a ridge of 1e-6 on all terms for numerical stability only.
- Station family: the plan said 8 tests; with the tip and barrier variants of contrast A it is 12.
- The CPC NAO file has two rows with the missing value `-99.000` glued to the day; they are read as missing (none fall in the sample).
- Two times (2004-11-01 00 and 12 UTC) were printed to test `stage1.py` before the full run, and the analysis was
  run once with 40 bootstraps to test the code before the 2,000-bootstrap run; the final numbers are the 2,000 run.

## Added after the first results (post hoc, not in any family)

`matched.py` diagnostic lines (match quality), `posthoc.py` (does the GH effect depend on how depth is adjusted? No:
3.6-5.0 with depth squared and cubed, depth x position, depth terciles; 3.8 inside the deepest tercile), `posthoc2.py`
(all-times fits and matching), `tables.py` (descriptive risk tables by tercile), and the station-wind regression in
`station.py`. In `tables.py`, within the deepest depth tercile the share of times with G_T >= 71.7 kt is 13.2%, 17.0%
and 12.1% from low to high GH (marginal, unadjusted), so the high effect is carried by the adjusted comparison and
by the middle depth tercile (0.6% to 4.9%), not by the deepest tercile's raw rates. Read the tables as
description, and the adjusted odds ratios as the estimate. They disagree in that one cell; the cause was not
established (position and NAO are correlated with GH within the deepest tercile).

## Limits

- 11 seasons per half. The held-out test for an effect near OR 1.4 per SD is underpowered.
- G_T is a gust over ocean within 300 km of Greenland, 58-72N, not pipeline A's 800 km index. It excludes
  points within 400 km of any detected low, so a storm's own wind field is not counted.
- Class is by position of the maximum only; wind direction at the maximum was not fetched. "Tip" and "barrier" are
  positional proxies for the jets.
- GH is MSLP over the ice sheet, a reduction to sea level, an index not a pressure that exists. The Z500 check
  (times before 2023-01-10) agrees in sign and size.
- Stations: ISD files for PCS (04390099999) and Tasiilaq (04360099999); only 308 and 335 of the 568 held-out rows have a
  report within +-1 h (ISD has gaps and some years are missing). Tasiilaq's
  reports never reach 34 kt in the sample, so its share-of-reports tests are degenerate (p set to 1).
- Pipeline A's domain stops at 67N and the lows come from it; a gust driver north of 67N is not a candidate.

## Reproduce

```
ERA5_WORK=/path/to/scratch python3 research/era5/greenland_jets/stage1.py 10        # 24.4 GB, resumable
ERA5_WORK=... python3 research/era5/greenland_jets/build_cases.py                    # counts, control sample, power
ERA5_WORK=... python3 research/era5/greenland_jets/stage2.py 8                       # 11.1 GB
ERA5_WORK=... python3 research/era5/greenland_jets/analyse.py 2000 [--swap]          # reads results/stage2_ingredients.csv
ERA5_WORK=... python3 research/era5/greenland_jets/{matched,station,composite,tables,posthoc,posthoc2}.py
```
`analyse.py` needs the CPC daily NAO file (`CPC_DIR`). `station.py` needs the ISD files in `$ERA5_WORK/isd/`
(`https://noaa-global-hourly-pds.s3.amazonaws.com/<year>/<station>.csv`, years 2004-2026).

## Departures from the plan

None to the pre-registered tests. Additions are listed above and labelled.
