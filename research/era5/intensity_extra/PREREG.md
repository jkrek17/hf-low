# Pre-registration: do more predictors improve the near-storm P(HF) model?

Written 2026-10-08 ~10:10 UTC, committed before any candidate feature was
extracted and before any skill was computed with a candidate in the model.

**What had been looked at when this was written.** The PR 12 results
(`../intensity/results/skill.txt`, `coefficients.csv`) were seen, so the
baseline's skill and which of its own predictors carry weight are known. No
candidate below has been computed on these fixes. WeatherBench2 chunk sizes
and one chunk's value ranges were inspected to cost the pull (this is how the
potential-vorticity candidate was found to be empty, see "Deviations before
the run").

**This is a proxy study.** Every outcome is pipeline A's (`../hf_history`)
ERA5 gust index reaching 71.7 kt. Every predictor is ERA5 at the fix time or
earlier. It is not the archive and not a forecast-model test.

## Question

The model of PR 12 gives P(gust index reaches 71.7 kt within 24 h) at every
pipeline A low fix from storm state, Hart phase space and nine environment
terms. Which further predictors add out-of-sample skill, and by how much?

## Sample

- Pipeline A fixes at 00/12 UTC from `../intensity/results/fixes_2004.csv.gz`,
  seasons 2004-05 to 2021-22 (18 seasons). Fixes before 2023-01-10 only,
  because WeatherBench2's 1.5 degree ERA5 ends there; 2022-23 is dropped as an
  incomplete season rather than patched from the 0.25 degree store (a
  pressure-level pull from that store is far above the 50 GB gate).
  Expected: 130,271 fixes on 29,820 tracks; 4,377 reach HF within 24 h.
- The baseline is refitted on this sample, so its skill will differ from PR 12's
  (22 seasons). Both are reported.
- Independent unit: the season (18).

## Baseline and candidate models

- Baseline `base` = PR 12's `full` set: STATE + HART + ENV columns of
  `../intensity/model.py`, logistic regression, L2, C = 1, median-impute,
  clip 0.5-99.5 %, standardise on the training seasons. Unchanged and not tuned.
- Candidate model = `base` plus one candidate group's columns, same fit.
- Fitting: leave-one-season-out (LOSO), 18 folds. Reference for skill:
  basin x calendar-month climatology fitted on the training seasons, as in PR 12.

## Candidate groups (fixed list; no additions after results are seen)

All fields are the WeatherBench2 1.5 degree ERA5 at the fix time (00 or 12 UTC).
Discs use great-circle distance from the fix, cos(latitude) weighting, as in
`../intensity/env.py`. The storm's heading is `heading` from the fixes file
(previous 6 h). Nothing after the fix time is used. Total precipitation is the
6 h accumulation ending at the fix time (WeatherBench2 convention, stated here
as an assumption; if its stamp turned out to mark the start of the window the
group is dropped and logged).

Tier 1 (WeatherBench2 single-level fields plus 250/850 hPa wind speed;
estimated 42 GB, measured below; under the 50 GB gate):

| Group | Columns | Definition | Why |
|---|---|---|---|
| H history | `dg12`, `dpacc` | `dg12` = g800(t) - g800(t-12 h); `dpacc` = [msl(t)-msl(t-12 h)] - [msl(t-12 h)-msl(t-24 h)]; both from the same track, 0 where the track is younger | rate and acceleration of the storm's own development; free, no pull |
| A air-sea | `airsea` | mean over ocean cells (land-sea mask < 0.5) within 500 km of SST - 2 m temperature (K) | low-level instability and heat-flux driver; the base has SST minus 500 hPa T and the flux but not this |
| M moisture transport | `ivt500`, `ivtmax` | WeatherBench2 integrated vapour transport, mean within 500 km and maximum within 1000 km (WB2 units; standardised) | atmospheric-river / warm-sector moisture feeding the low; the base has column water vapour but no transport |
| D diabatic | `precip6` | mean 6 h precipitation within 500 km (m) | latent heating already under way |
| J jet geometry | `jet_along`, `jet_cross` | location of the maximum 250 hPa wind speed within 1500 km, as displacement from the fix in the heading frame (100 km): along = toward the direction of motion, cross = to the right of motion | position relative to the jet (entrance / exit, left / right of the jet), which `jet250` (maximum speed) does not carry |
| L low-level wind | `ws850` | maximum 850 hPa wind speed within 500 km (kt) | low-level jet, the wind that mixes down; expected to be largely redundant with the current gust index |
| P pressure environment | `dphigh`, `ringgrad` | `dphigh` = maximum MSLP within 2000 km minus minimum MSLP within 300 km of the fix (hPa); `ringgrad` = mean MSLP gradient magnitude on the 300-800 km annulus (hPa per 100 km) | the neighbouring high and the background gradient |
| K land | `landfrac` | area fraction of land within 500 km | coastal and orographic influence; the gust index is ocean-only |

Tier 2 (pre-registered now, **not run until Jason gives the go-ahead**: needs
about 88 GB more, above the 50 GB gate):

| Group | Columns | Definition |
|---|---|---|
| N stability | `stab` | mean over ocean cells within 500 km of (theta 700 hPa - theta 925 hPa) / 225 hPa (K per hPa) |
| R trough | `trough` | mean 500 hPa height on the 1500-2500 km annulus minus the minimum 500 hPa height within 1000 km (m) |
| W warm conveyor | `omega700` | minimum (strongest ascent) 700 hPa vertical velocity within 500 km (Pa s-1) |

Not in the list, and why:

- Translation speed is already in the base (`speed`).
- Surface heat flux, SST, SST gradient, column water vapour, Eady growth,
  300 hPa divergence, 500 hPa vorticity advection, 250 hPa maximum wind,
  Hart B / -V_T are already in the base.
- Upper-level potential vorticity: the WeatherBench2 PV array is empty
  (all NaN) in the chunk inspected; deriving it needs u, v and T at several
  levels (about 70 GB more). Trough (R) and jet geometry (J) stand in. It is
  not a Tier 2 member.
- Threshold behaviour of jet speed and trough depth belongs to the thread
  "Jet and trough thresholds for bombs and HF".

## Evaluation

**Primary target.** `hf24`: gust index reaches 71.7 kt in (t, t+24 h], all fixes.
**Primary metric.** Gain = BSS(candidate) - BSS(base), both against the same
climatology, from the pooled LOSO Brier scores. The 90 % interval resamples
whole seasons (2,000 draws).
**Test.** One p-value per group from an exact two-sided sign-flip test on the
18 per-season Brier-score differences (2^18 flips of their sum). Benjamini-
Hochberg across groups for the FDR.

**Secondary targets, all reported.** (a) `hf24` restricted to fixes not already
HF at t (onset). (b) Rapid deepening (normalised deepening rate >= 1 Bergeron,
binary logistic with the same columns), Brier gain. (c) Yes/no skill at the
count-matched cut: HSS change for hf24 and for rapid deepening, with season-
bootstrap intervals (descriptive, outside the FDR). (d) Gain by basin
(descriptive).

**FDR family.** Every group x {hf24, hf24 onset, rapid deepening}. q is
reported both within the primary target and over the whole family. If Tier 2
is run later, q is recomputed over all groups and the Tier 1 values are marked
provisional.

**Decision rule (fixed now).**

- *Helps*: primary gain >= +0.005 BSS, 90 % interval above zero, and q < 0.05
  within the primary family.
- *Does not help*: the 90 % upper bound of the primary gain is below +0.005
  (a well-powered null: gains of half a point of BSS or more are excluded).
- Otherwise *inconclusive*.

0.005 is about 1 % of the base's BSS (0.42). It was chosen as a size a
forecaster would notice, before any candidate was scored.

**Combinations.**

- `all1`: base plus every Tier 1 group (no selection).
- `sel1`: base plus the groups judged *Helps* above. Selected on the same
  seasons it is scored on, so its gain is optimistic; this is stated wherever
  it appears.
- `sel1_nested`: the honest version of `sel1`. In each outer fold the groups
  are chosen from the 17 training seasons by six-block inner cross-validation
  (blocks of three seasons) with the rule "inner gain >= +0.005", then fitted
  and scored on the held-out season.
- If no group is judged *Helps*, `sel1` and `sel1_nested` are not built and
  that is the result.

**Cross-checks planned.** The base's LOSO on the 18-season sample against
PR 12's 22-season value (0.423). A fresh Sonnet verifier recomputes the
committed numbers from the committed feature table.

## Deviations before the run

- PV candidate removed (empty array), see above. Recorded before any skill was
  computed.

## Post hoc

None yet. Anything added after the first candidate skill is computed will be
logged here with the date and reported separately.
