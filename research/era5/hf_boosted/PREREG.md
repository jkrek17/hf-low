# Is the boosted-tree gain on P(HF within 24 h) real? Pre-registration (RA-21)

Written before any boosted-tree or flexible-logistic model was scored under leave-one-season-out. **What had been looked at:** PR 12's logistic headline (hf24 BSS 0.423, onset-only 0.291, HSS 0.59), and PR 68's post hoc boosted-tree line (one default run, 11-fold season-grouped CV: BSS 0.456, HSS 0.61; with the tier-1 fields 0.455). No LOSO score of any tree or spline model, no partial dependence, and no per-stratum tree score existed. The code in this directory was run only for timing and shape checks with scores suppressed.

**Scope.** Pipeline A (`research/era5/hf_history`), ERA5 proxy. Target hf24 = the gust index reaches 71.7 kt in (t, t+24 h]. Sample, loader and design are PR 12's: `intensity/model.py` `load()`, 159,430 in-domain 00/12 UTC fixes, 22 seasons 2004-05 to 2025-26 (Decision 1: nothing earlier), the `full` predictor set (23 columns), leave-one-season-out (LOSO). Skill = Brier skill against basin-month climatology fitted on the training seasons. Effective n is 22 seasons; every interval resamples whole seasons (1,000 draws, seed 20261010, 90% intervals as in PR 12).

**Looks spent.** The 22 seasons have been used by PR 12, 56, 58, 64, 68 and others (the 2015-25 block is on its sixth logged look in `hemispheric/results/heldout_looks.log`). This note spends **two looks**: (L1) LOSO over all 22 seasons, (L2) one temporal split, train 2004-14, test 2015-25, 2025-26 included. They are not independent of PR 68's run, which used the same seasons. They test whether the gain survives the stricter fold design and shows up season by season, not whether it holds on seasons nobody has seen. Fresh seasons need RA-4 or time.

## Models (settings fixed now, none tuned)

| Name | Model |
|---|---|
| M0 | PR 12 logistic: L2, C = 1, median-impute, clip 0.5-99.5%, standardise (`model.Prep`, `model.fit`), the `full` set |
| M1 | `HistGradientBoostingClassifier(max_iter=300, learning_rate=0.06, max_leaf_nodes=15, l2_regularization=1.0, early_stopping=False, random_state=0)` on the same 23 raw columns (NaN passed through). These are PR 68's settings; nothing was tuned on any score |
| M2 | M1 with `interaction_cst` = each feature alone: additive trees, any single-feature nonlinearity, no interaction |
| M3 | M0 plus, for each of the 18 continuous predictors, a cubic spline basis (`SplineTransformer`, 5 knots at quantiles, degree 3, constant extrapolation) fitted on the training rows after `Prep`; same L2, C = 1. Smooth curvature, no interaction. Continuous predictors: msl, dp12, lat, speed, g800, logage, B, VTL, VTU, jet250, div300, vadv500, eady, sst, sstgrad, sst_t500, flux, tcwv. The rest (young, pac, doy_c, doy_s, nosst) stay linear |

The logistic M0 is refit inside this note's LOSO so the comparison is on identical fixes and folds; its BSS must reproduce PR 12's 0.423 to within 0.002, otherwise the loader or fold design is wrong and no tree score is reported until that is fixed.

## Primary test (one confirmatory test)

**P1.** Pooled LOSO BSS gain ΔBSS = BSS(M1) - BSS(M0), season-bootstrap 90% interval, and the number of seasons (of 22) in which M1 has the lower Brier sum. Added: a sign-flip permutation p on the 22 season-level Brier differences (100,000 flips).

Decision rule, written as thresholds on the point estimate and the interval:
- **Real, as predicted:** ΔBSS ≥ +0.020 and the lower 5% bound > 0.
- **Real, smaller than predicted:** +0.010 ≤ ΔBSS < +0.020 and the lower bound > 0.
- **Not real (the 0.03 was selection or a fold artefact):** ΔBSS < +0.010.
- **Can't tell:** ΔBSS ≥ +0.010 with the lower bound ≤ 0.

Agenda prediction (RA-21): at least +0.02, concentrated in fixes whose current g800 is 55-71.7 kt. A gain under +0.01 means the post hoc 0.03 was selection.

**P1b (the stratum prediction, part of the same prediction).** Strata by current g800: below 55 kt, 55-71.7 kt, at or above 71.7 kt. The prediction is met if the 55-71.7 kt stratum carries at least 60% of the total Brier reduction (it holds about 13% of fixes and about 54% of the Brier loss in PR 68, so 60% means more than proportional) and its per-fix gain exceeds both other strata's. Season-bootstrap interval on the share. Reported as met / not met; not entered into an FDR family.

**Power and what it can show.** The interval comes from 22 seasons, so its width is the season-to-season SD of the difference, not the fix count. After the run the standard error of ΔBSS is reported and the 80% minimum detectable gain is stated as 2.8 x SE. A result under +0.010 is "not real" only if the interval's upper bound is also under +0.020; otherwise it is "inconclusive".

## Secondary tests (one Benjamini-Hochberg family, all reported; q across all tests with P1 included also reported)

All are LOSO differences in Brier sum as a share of the climatology Brier sum, season-bootstrap 90% interval, sign-flip p.
- **S1** onset only (not HF at t): M1 - M0.
- **S2a, S2b** Atlantic only, Pacific only: M1 - M0.
- **S3** temporal split L2: M1 - M0, trained 2004-14, tested 2015-25. M0 here is PR 12's early-trained logistic refit (its reported 0.426 is the contrast).
- **S4** HSS at the pooled count-matched cut (forecast count = observed count over all held-out predictions, the same cut rule for every model, as PR 68 used): HSS(M1) - HSS(M0), season bootstrap. The training-quantile cut (PR 12's rule) is shown as a contrast.
- **S5** M2 - M0: how much is single-ingredient nonlinearity. **S6** M1 - M2: how much is interaction. **S7** M3 - M0: how much of S5 smooth curvature explains. **S8** M1 - M3: what trees add beyond smooth additive curves, that is sharp thresholds plus interactions.
- **S9** (18 tests) M0 with a spline for one continuous predictor at a time (same basis as M3, one feature) minus M0: the single ingredient whose nonlinearity is worth the most.
- **S10** (5 tests) M0 plus one product term of two standardised predictors, minus M0. Pairs fixed now, chosen for mechanism not for results: g800 x dp12 (strong and already deepening), msl x dp12, jet250 x eady (jet and baroclinicity), sstgrad x flux (surface heating), lat x g800 (high latitude).

Count: P1 (1) + S1 (1) + S2a/b (2) + S3 (1) + S4 (1) + S5-S8 (4) + S9 (18) + S10 (5) = 33 tests, 32 of them secondary. (The first draft of this line said 32 in total; corrected before any model was scored.)

## Partial dependence and interactions (run only if P1 gives "real, as predicted" or "real, smaller"; descriptive, no p values)

M0 and M1 refit on all 22 seasons. For each of the 18 continuous predictors: partial dependence on 5,000 randomly drawn fixes (seed 20261010), 25 grid points from the 2nd to the 98th percentile, shown on the logit scale. Report where M1's curve departs from M0's: the grid value of the steepest rise (a candidate threshold), the range of the curve, and the departure from a straight line (R² of a linear fit to the M1 logit curve). Friedman's H² for the five S10 pairs under M1. Statements about which thresholds the logistic form misses are tied to the S9/S10 and S5-S8 results above, not to the plots alone. The "Jet and trough thresholds for bombs and HF" thread (pre-registered, `research/era5/jet_trough`) tests threshold against ramp for jet speed and trough depth using Z500 features I do not have here; this note does not test those thresholds and does not use its features. jet250 here is the PR 12 column, and any statement about it is reported as a contrast to that thread's result, not a duplicate.

## What this study cannot show
Seasons are the same ones PR 68 saw; no fresh data. Tree hyperparameters are PR 68's untuned defaults, so a tuned model could gain more and this is a floor for tree skill, not a ceiling. The target is pipeline A's own gust index, so a "threshold" is a property of the proxy's relation to the 1.5 degree fields, not of the atmosphere. Fixes within a storm are correlated; intervals resample seasons, which absorbs that.

## Deviations (post hoc)
Written after the scores were seen. Nothing registered was changed; these are additions and clarifications.
1. P-value display: the sign-flip test uses 100,000 flips, so a p of 0 is displayed as < 1e-5.
2. P1b (the 55-71.7 kt prediction) came out NOT MET as registered (share 40%, rule 60%). Not redefined. Post hoc PH1 shows where the beyond-smooth part of the gain sits.
3. S6 (M1 - M2) mixes interaction with under-fitted additive trees (S5 is less than half of S7); interpret S8 and PH3 instead. Reported as registered.
4. The partial-dependence grid (2nd-98th percentile) ends at 71.1 kt for g800, before the already-HF fixes; PH2 extends it to the 99.9th percentile.
5. Post hoc additions (`posthoc.py`, `results/posthoc.txt`): PH1 strata shares of M1 - M3 and M0 - M3; PH2 extended g800 partial dependence; PH3 M3b = M3 + the five pre-specified products, LOSO; PH4 spline blocks for B, VTL, VTU (with and without g800) and for seven other predictors picked from the S9 table. Both the registered and the post hoc numbers are reported; none changed a verdict.
6. The run_loso.py `clim()` helper computes the basin-month climatology with the same +0.5 smoothing as `model.clim_probs`; the M0 reproduction of PR 12 (0.4228 against 0.423) passed the pre-set band before any tree score was read.
