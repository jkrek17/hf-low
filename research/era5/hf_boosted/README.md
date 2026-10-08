# Is the boosted-tree gain on P(HF within 24 h) real? (RA-21)

Plan: [PREREG.md](PREREG.md), committed before any tree or spline model was scored (`242a764`; test count corrected in the next commit, still before scoring). Its deviations log was written afterwards and says so.
Code: [run_loso.py](run_loso.py) (all LOSO fits, about 3 min), [analyse.py](analyse.py), [pdp.py](pdp.py) (about 2 min), [posthoc.py](posthoc.py) (post hoc). They read committed tables only, no ERA5 pull:

    W=$(mktemp -d)
    python3 -I research/era5/hf_boosted/run_loso.py . $W
    python3 -I research/era5/hf_boosted/analyse.py . $W research/era5/hf_boosted/results
    python3 -I research/era5/hf_boosted/pdp.py . research/era5/hf_boosted/results
    python3 -I research/era5/hf_boosted/posthoc.py . $W research/era5/hf_boosted/results

**Pipeline A, ERA5 proxy.** PR 12's `hf24` sample: 159,430 fixes, 22 seasons 2004-05 to 2025-26, 23 predictors, leave-one-season-out (LOSO); Brier skill against basin-month climatology; intervals resample whole seasons (1,000 draws, 90%). The logistic refit reproduces PR 12 (BSS 0.4228 against 0.423). Results: [results/loso_results.txt](results/loso_results.txt), `tests.csv`, `strata.csv`, `pd.txt`, `pd.png`, `posthoc.txt`, per-fix held-out probabilities in `loso_probs.csv.gz`.

## Answer

**Yes, the gain is real, and it is about the size PR 68 found (+0.032 BSS), not a fluke of one run. It is not new information: most of it is curvature the straight-line logistic cannot draw, concentrated in the storm-structure variables.**

| | BSS | gain over PR 12 logistic |
|---|---|---|
| M0 logistic (PR 12 refit) | 0.4228 | |
| **M1 boosted trees** (PR 68 settings, untuned) | **0.4549** | **+0.0321 [+0.0261, +0.0375]**, better in 21 of 22 seasons, sign-flip p < 1e-5 |
| M3 smooth additive curves (splines, logistic) | 0.4455 | +0.0226 [+0.0202, +0.0253], 22 of 22 seasons |
| M2 additive trees | 0.4334 | +0.0106 |

Pre-registered verdict (P1): **real, as predicted** (at least +0.020 and lower bound above 0). The standard error is 0.0034, so a gain of 0.0095 would have been detectable at 80% power. It holds in both basins (Atlantic +0.036, Pacific +0.028), for onset fixes alone (+0.026), and out of time: trained 2004-14, tested 2015-25, +0.026 (BSS 0.4523 against 0.4264). HSS at the count-matched cut rises from 0.588 to 0.615. 18 of the 32 secondary tests have q < 0.05 (19 of 33 with P1); all are in `tests.csv`.

**The stratum half of the prediction was not met.** The agenda said the gain would sit in the 55-71.7 kt stratum. It holds 40% of the gain [28, 52] (13% of fixes); fixes already at 71.7 kt or more hold 47% [36, 59] (2% of fixes), and fixes below 55 kt hold 13%. Rule was at least 60% in the middle stratum: not met. Per fix, the gain is largest for storms already above the cut (they either stay HF or decay within 24 h).

## What the trees capture

Decomposition by nested LOSO fits (S5-S8, and post hoc PH3-PH4; not exactly additive, but the pieces add to 0.032):

| Piece | Gain in BSS | Share of +0.032 |
|---|---|---|
| Smooth single-ingredient curvature (M3 over M0) | +0.0226 | about 70% |
| ... of which the three Hart phase-space parameters B, VTL, VTU alone | +0.0112 (22 of 22 seasons) | about 35% |
| ... Hart parameters plus the current gust g800 | +0.0141 | about 44% |
| The five pre-specified pairs as product terms (PH3) | +0.0041 [+0.0029, +0.0055] | about 13% |
| Left for trees beyond both (sharp edges, other interactions; PH3) | +0.0053 [+0.0006, +0.0098] (13 of 22 seasons) | about 17%, borderline |

- **Which single ingredient (S9, one spline at a time, BH q across 18):** B +0.0059, VTL +0.0057, VTU +0.0054 (q < 0.001 each), g800 +0.0026 (q 0.001), sea-surface temperature +0.0021 (q 0.04), surface heat flux +0.0011 (q 0.04). Jet speed at 250 hPa +0.0007 (q 0.13), column vapour -0.0002, sea-surface-temperature gradient +0.0006 and the Eady growth rate +0.0011 (q 0.10) do not clear. Nothing about the jet speed column suggests a threshold here; the "Jet and trough thresholds" thread tests that on its own Z500 features and this note does not duplicate or test it.
- **The shapes (`pd.png`, descriptive).** The logistic draws a straight line; the trees show (i) a floor: P(HF within 24 h) is flat below about 35 kt current gust and only then rises, steepening toward the cut; (ii) tails in the phase space: the risk rises in the strongly cold-core tail (VTL below about -250, VTU below about -200; positive means warm core) and for very asymmetric storms (B above about 75), and is flat otherwise; (iii) saturation: heat flux adds little above about 100 W m-2 and the pressure effect flattens below about 980 hPa. The latitude and 12 h pressure-change curves are close to the logistic's.
- **Interactions are small.** Friedman H-squared under the trees is at most 0.003 for the five pairs (gust x deepening rate 0.003, heat flux x SST gradient 0.0025, the rest under 0.0015). As product terms they are worth +0.0041 together; gust x deepening rate +0.0016 (q 0.001), latitude x gust +0.0014 (q 0.010), pressure x deepening rate +0.0007 (q 0.007). Jet x Eady -0.0001 (q 0.19) adds nothing.
- **Latent heat.** Surface heat flux has a small saturating nonlinearity (+0.0011); column water vapour has none (-0.0002, q 0.71). Consistent with PR 68: the surface proxies do not carry a hidden threshold. Warm-conveyor ascent aloft remains untested (Tier 2).
- Where the beyond-smooth part sits (PH1, post hoc): 83% [60, 127] of M1 over M3 is in fixes already at 71.7 kt or more, i.e. whether a storm that is already HF stays HF, which the smooth curves do not fit sharply. The smooth part (M0 to M3) is spread over the near-threshold and below: 49% in 55-71.7 kt, 32% above, 19% below.

Plain reading: the forecast form leaves about +0.03 BSS on the table, mostly because the phase-space parameters act in the tails and the gust acts from a floor, not as straight lines. A model that adds spline terms for those four variables recovers about 44% of it with no new data, a full boosted model recovers it all, and nothing here points to a missing physical ingredient: the gain is in how the existing ones enter.

## Looks, deviations, limits

- **Looks spent: 2** on the 22 gust-era seasons (LOSO over all 22; one temporal split, 2015-25 test). The same seasons were used by PR 68's 11-fold run, so this confirms the gain is robust to the fold design and visible season by season, not that it holds on unseen seasons. Appended to `hemispheric/results/heldout_looks.log`.
- **Post hoc, labelled:** PH1-PH4 (`posthoc.txt`); the p-value display floor (100,000 flips, so p = 0 is shown as < 1e-5); the stratum prediction was reported as not met rather than redefined; the partial-dependence grid (2nd-98th percentile) stops at 71.1 kt for g800, before the HF fixes, so `posthoc_g800_pd.csv` extends it to the 99.9th percentile (M1 rises from about -4.5 at 67 kt to about -1.0 at 91 kt on the logit scale, steeper than the logistic's line).
- **S6 (M1 minus M2, +0.0215) overstates interaction.** The additive trees (300 trees spread over 23 features) fit single-feature curves worse than splines do (S5 +0.0106 against S7 +0.0226), so the difference of M1 and M2 mixes interaction with under-fitted curves. S8 and PH3 are the cleaner reading. Reported as it was registered.
- Partial dependence on collinear predictors (SST, column vapour, SST minus T500, Eady) is not a causal curve; the logistic's large opposite-signed weights on them compensate each other. The H-squared and PD are in-sample and descriptive.
- Tree settings are PR 68's untuned defaults, so +0.032 is a floor for tree skill. The target is pipeline A's gust index: a "floor" or "threshold" describes the proxy's relation to 1.5 degree fields, not the atmosphere. Atlantic fixes north of 60N are included as in PR 12 (lower confidence). Transitioning tropical cyclones are in.
- Hart sign convention: positive VTL/VTU is warm core (`intensity/hart.py`).

## Verification

VERIFY_PLACEHOLDER
