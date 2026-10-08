# What turns sub-55 kt storms into HF lows within a day? Pre-registration (agenda RA-20, Arms A and B)

Written and committed **before any ERA5 field for these times was read**. ERA5 proxy throughout; outcomes are **pipeline A**
(`research/era5/hf_history`: 800 km ocean gust index, HF-equivalent at 71.7 kt), a proxy for the archive, not the archive.
Transitioning tropical cyclones are in. Atlantic fixes north of 60N are lower confidence (no terrain mask).

**What had been looked at.** Only labels and sizes: the counts below, the pull sizes (`results/sizes.txt`, HEAD requests, no
field downloaded), and two power simulations that use labels and the leave-one-season-out (LOSO) baseline probabilities
(`power.py`, `power_bss.py`). Earlier work on the same fixes: PR 12 (logistic P(HF within 24 h)), PR 56 (eight extra groups),
PR 68 (error analysis: about 41% of misses are storms below 55 kt, POD 0.06 there), PR 76 (boosted trees, BSS +0.032), PR 85
(the pattern acts on conversion; about two thirds of that effect sits outside the 1.5 degree near-storm fields). None of them saw
a 0.25 degree field, a vertical velocity or a mesoscale feature at these times.

## Question and the answer I will give

Among storms still below 55 kt, do (B) mesoscale surface structure at native 0.25 degrees, or (A) fields aloft (stability, trough
depth, warm-conveyor ascent), separate those that reach HF within 24 h from those that do not, beyond the PR 76 boosted-tree
model? The answer will be one of **yes** (name the group and the size), **no (well-powered)** (state the minimum detectable
effect), or **can't tell**, per group.

## Population, cases, controls (all fixed now from `intensity/results/fixes_2004.csv.gz`)

- **Stratum S.** Pipeline A fixes at 00/12 UTC, seasons 2004-05 to 2025-26 (decision 1), current gust index `g800` < 55 kt and not
  already HF (`hf_now` false). 135,663 fixes. Fixes whose track leaves the domain within 24 h are kept, as in PR 12, 68 and 76
  (32 of the cases). Stratum base rate 0.0068.
- **Cases.** Fixes in S with `hf24` true (the gust index reaches 71.7 kt in the next 24 h): **927** (Atlantic 531, Pacific 396), on
  803 tracks, 31 to 54 per season, all 22 seasons. This is above the 400 the agenda set as the floor for a conclusive test.
- **Pull design.** Times with at least one case: **895** (census). Plus **600** other times drawn at random (seed 20261008,
  `select_times.py`; list in `results/times_armB.csv`), each carrying weight 25.29 (15,173 other times / 600). The weighted
  sample is an unbiased stratified sample of S, so a stratum Brier score is estimable. Times were chosen by the outcome (cases) on
  purpose; nothing is tuned afterwards.
- **Controls.** Every other fix in S at the same time as a case (6,369 controls, 6.9 per case): the synoptic background is shared.
  **Design refinement against the agenda wording** ("controls matched on current gust, deepening rate and month, deepening fixes"):
  a 1:k match (same time, gust within 3 kt, 12 h pressure change within 2 hPa, or both missing) exists for only 198 of 927 cases
  (393 cases have no 12 h change because the storm is new), so exact matching would discard most cases. The adjustment instead
  goes through the baseline model's logit, which already holds gust, 12 h change, age, latitude and month. The 198-case matched
  subset is reported as a sensitivity result. Decided before any field was read.

## Baselines (existing, not refitted)

`hf_boosted/results/loso_probs.csv.gz`: column **M1** (boosted trees, PR 76's settings) is the baseline that must be beaten;
**M0** (PR 12 logistic) is the secondary baseline. Both are LOSO by season, so a fix is never scored by a model that saw its
season. Those probabilities enter every model below as a fixed offset.

## Features (frozen; every statistic is computed from the field at the fix time t only, never later)

Fix-centred, over **ocean cells within 800 km of the fix that the fix owns** (nearest of the lows present at that time, pipeline A's
ownership rule). The extraction must reproduce the table's `g800` to within 0.5 kt for at least 99% of fixes before any feature is
used; if it does not, stop and report.

**Arm B, ARCO-ERA5 0.25 degree hourly** (instantaneous 10 m gust, 10 m u and v, mean sea-level pressure). Four groups, 11 features:

| Group | Feature | Definition |
|---|---|---|
| B1 gust-field concentration | `A45` | log(1 + area in km2 of owned ocean cells with gust >= 45 kt) |
| | `PK` | gust maximum minus the 99th percentile of owned gust (kt) |
| | `RG` | distance from the fix to the gust maximum (km) |
| B2 pressure structure | `SHARP` | MSLP at the fix (0.25 deg) minus the Gaussian-smoothed (sigma 1.5 deg) MSLP at the fix (hPa) |
| | `CURV` | Laplacian of MSLP at the fix from a quadratic fit to the 0.25 deg points within 150 km (hPa per (100 km)^2) |
| | `GRAD` | 90th percentile of the MSLP gradient magnitude within 400 km (hPa per 100 km) |
| B3 sustained wind vs gust | `WS99` | 99th percentile of owned 10 m wind speed (kt) |
| | `GF` | gust at the gust-maximum cell divided by 10 m wind speed there |
| | `A34` | log(1 + area in km2 of owned cells with 10 m wind >= 34 kt) |
| B4 shear and convergence | `SHR` | maximum horizontal gradient of 10 m wind speed over owned cells (kt per 100 km) |
| | `CONV` | minus the minimum 10 m horizontal divergence within 400 km (1e-5 per second; larger is stronger convergence) |

**Arm A, WeatherBench2 1.5 degree 6-hourly** (temperature, geopotential, vertical velocity; the Tier 2 groups of PR 56's plan,
`intensity_extra/PREREG.md`, unchanged): `stab` (mean over ocean cells within 500 km of (theta700 - theta925) / 225 hPa),
`trough` (mean 500 hPa height on the 1500-2500 km annulus minus the minimum within 1000 km), `omega700` (minimum 700 hPa vertical
velocity within 500 km). Three groups of one feature each. WeatherBench2 ends 2023-01-09, so Arm A covers the times up to then
(1,250 of 1,495 times, 740 of 895 case times, about 18.6 seasons) and is scored on that subset against the same baseline there.

## Models and metrics

For each group g with standardised features z (mean and SD from the training seasons): **logit p = logit(p_M1) + a + b'z**, with
a ridge penalty of 1 on b (fixed, never tuned), fitted by weighted maximum likelihood on the sampled fixes of the training seasons,
LOSO by season. Also an "all" model per arm (all that arm's features). Same models over M0 as the secondary baseline. One nonlinear
check, fixed now: `HistGradientBoostingClassifier(max_iter=100, learning_rate=0.06, max_leaf_nodes=8, l2_regularization=1.0,
early_stopping=False, random_state=0)` on [logit p_M1, the arm's features], sample-weighted, LOSO.

1. **Primary test, within-time (case-control).** For each group, a likelihood-ratio test (df = number of features in the group) of
   the group's coefficients in a Breslow-type conditional logit on the case-time risk sets (each case against all S fixes at its
   time, baseline logit as an offset), fitted on all 22 seasons. This asks whether the group separates cases from same-time
   controls beyond M1. Fixes within a storm are not independent, so the p-value is the larger of the model-based LR p and a
   season-cluster bootstrap p (2,000 resamples of the 22 seasons, LR recomputed); the decision uses the larger. Out-of-sample
   evidence goes beside it: the held-out log-likelihood gain per case, LOSO, with a sign-flip p over the 22 season differences
   (100,000 flips), reported for every group but not used for the decision because it has far less power (below).
2. **Primary effect size, forecast.** Stratum BSS gain DeltaBSS_S = (Brier of M1 - Brier of augmented) / Brier of the stratum
   climatology, on the weighted sample, pooled LOSO, with a 90% season-bootstrap interval. The agenda's bar is +0.005. **This is
   a weak test**: see the power section; its sign-flip p is reported but is not the decision test.
3. Secondary (all reported, same family accounting): the same two against M0; the matched 198-case subset; the tree check; per
   basin; cases split by `dp12` available or missing (new storms).

**Families.** Primary tests are test 1 for the four B groups and the B "all" model (5 tests) and, if Arm A runs, the three A groups
and the A "all" model (4 tests). Benjamini-Hochberg q within an arm, and over all primary tests together; both reported. Secondary
tests form a second family with its own q, and the total count "passing FDR out of all tests run" is reported.

## Decision rules (thresholds on estimates)

- **Yes, the group separates and meets the agenda bar:** test 1 q < 0.05 and DeltaBSS_S >= 0.005 (point estimate; the interval is
  reported and a lower bound at or below 0 is said in the same sentence). The agenda prediction for Arm A (warm-conveyor ascent adds
  at least 0.005) is read this way for `omega700`.
- **Separates, forecast gain below the bar:** test 1 q < 0.05 and DeltaBSS_S < 0.005.
- **No (well-powered):** test 1 q >= 0.05 and the group's 80% minimum detectable effect (below) is at or below 0.15 SD per feature.
  Otherwise **can't tell**.
- If no group of either arm passes, the answer is that explosive onsets from sub-55 kt are **not predictable from the resolved
  fields tested**, which closes the Tier 2 question with evidence. A null in Arm B means ERA5's own 0.25 degree structure does not
  separate them either, and the limit is sub-grid or random.
- Hypotheses stated in advance: **H_B** at least one B group passes; **H_A** `omega700` adds at least +0.005 BSS. If both are
  false, "no" for both.

## Power (labels and baseline probabilities only; `power.py`, `power_bss.py`, `results/power*.txt`)

Planted effects: an independent standard-normal feature changes the log-odds of a case by beta per SD. All numbers below are
from the committed result files; none used a field or an outcome other than the labels.

- **Within-time likelihood-ratio style test** (`power.py`, `results/power.txt`; score test on the real 895 risk sets, 927 cases,
  6,369 controls, feature independent of the baseline, 400 simulations): standardised case-minus-control shift of 0.05 SD has
  power 0.24, **0.10 SD 0.79**, 0.15 SD 1.00 at alpha 0.05; at alpha 0.05/6 the values are 0.08, 0.55 and 0.97. **The 80% minimum
  detectable effect is about 0.10 SD per feature (about 0.12 SD at the multiplicity-corrected level).** A feature correlated with the
  baseline offset has less independent information and needs a larger shift.
- **Stratum BSS gain** (`power_bss.py`, `results/power_bss.txt`; full re-draw of outcomes for the 135,663 stratum fixes, the pull
  design applied, LOSO fit, 60 simulations per beta, sign-flip over 22 seasons): beta 0.2 gives DeltaBSS +0.002 (detected 3%),
  **beta 0.3 gives +0.005 (detected 7%)**, beta 0.45 gives +0.014 (22%). The agenda's +0.005 bar is therefore reachable by an effect
  that this metric has about 7% power to confirm; **DeltaBSS_S is an effect-size estimate, not a decision test.**
- **Held-out within-time log-likelihood gain** (same simulations): beta 0.2 detected 18%, beta 0.3 62%, beta 0.45 98%.
- Consequence: a single feature that shifts a case by 0.10 SD or more is found if present; one that meets the +0.005 bar (beta about
  0.3, about 0.3 SD of shift) is found by the within-time test with near certainty and by the out-of-sample versions with moderate
  power. The multi-feature groups have different MDEs; `power.py` is re-run with the real feature correlation matrix (features
  joined to the risk sets, not to the outcome) before any test is run, and its output is committed.


## Looks and what this cannot show

- Spends **one look at the 2015-25 block per arm run** (LOSO scores all 22 seasons, including 2015-25), logged in
  `hemispheric/results/heldout_looks.log` by timestamp, using that log's true numbering. No pre-2001 season is used (gust-based).
- Seasons are not independent of PR 12, 56, 68 and 76, which used the same fixes; the baseline is not fresh data and this is
  not a replication.
- ERA5 under-resolves the strongest winds and the melting layer (so `A45`, `PK`, `omega700` may be damped); a null does not rule
  out structure real storms have. Arm A is 1.5 degree and cannot see mesoscale ascent.
- Features are measured at t, in the same fields that define the gust, so Arm B is partly a different view of the same ERA5 gust
  field. It is a test of what that field's structure adds to the forecast, not an independent observation of the storm.
- The `omega700` minimum and `trough` may be correlated with the baseline's jet and Eady terms; the model handles it, but a gain
  is then a partial correlation, not proof of mechanism.
- Effective n is seasons (22; 18.6 for Arm A) for the forecast metric, and case-time risk sets (895) for the within-time one.
  Fixes within a storm are not independent; the season-block sign-flip and bootstrap respect that, the conditional likelihood
  itself does not.

## Deviations (post hoc)

None yet.
