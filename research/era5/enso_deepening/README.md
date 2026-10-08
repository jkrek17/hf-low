# Does ENSO flavor move where North Pacific storms deepen fastest? (ERA5 proxy, pipeline A detector)

Plan, committed before any deepening position was joined to an index: [PREREGISTRATION.md](PREREGISTRATION.md) (commit `c6fcd34`).
Code: `deepening_points.py` (one max-deepening point per track), `analysis.py` (all 18 tests, power). Results: [results/](results/) (`summary.txt` is the readable report, `results.csv` has every test, `power.csv` the planted-effect power, `winter_table.csv` the winter means).

All numbers are an **ERA5 proxy**. Tracks are PR 38's MSLP-only tracker output, which reuses **pipeline A's** low detector and linker (`research/era5/enso_kuroshio/data/`); the HF label (S11 only) is pipeline A's gust index at 71.7 kt. Pipeline B is not used. Tropical cyclones are left in (the 25N floor on the deepening point removes most). No data was pulled.

## Answer

**No, and the test is strong enough to say so for the size the idea predicts.** Over 47 winters (1979-80 to 2025-26), one SD more of the Modoki index (EMI, higher = warm anomaly farther west), with Nino3.4 and a linear trend held, moves the longitude where Pacific bombs (24 h Bergeron-normalized deepening of at least 24 hPa, about 33 a winter) deepen fastest by **-1.09 degrees (95% CI -2.44 to +0.34; permutation p = 0.17; q = 0.47)**, and the latitude by **-0.07 degrees (-0.56 to +0.37)**. The registered rule for a well-powered null (interval inside the shift the idea predicts, 3 degrees longitude and 2 degrees latitude per SD, and 80% power at no more than that) is met for both: planted-effect power is 80% at about 2.2 degrees longitude per SD and under 1 degree latitude per SD. **None of the 18 registered tests passes FDR** (smallest q = 0.35).

What the data leave open: a westward lean of **up to about 2.4 degrees of longitude per SD** (the interval's lower end) cannot be excluded, and 14 of the 15 longitude tests have a negative sign (higher EMI, deepening point farther west). Those tests share the same outcome, so this is one weak lean seen many times, not many pieces of evidence. The winter-to-winter spread of mean bomb longitude is 3.7 degrees SD, so a 1 degree lean per SD would explain under 10 percent of its variance.

This matches PR 38 (genesis over the Kuroshio, RR 0.995 per SD, well-powered null): the flavor of El Nino does not move where storms form and does not move, at the 3-degree scale, where they explode. It closes the ENSO-flavor idea for the Pacific at that scale. It does not test jet position itself (no wind fields were used).

## Numbers (degrees per +1 SD of the index over the 47 winters; 10,000 permutations; 5,000 winter bootstrap resamples; BH over 18 tests)

| id | what | n winters | effect | 95% CI | detectable at 80% | p (perm) | q (BH) |
|---|---|---|---|---|---|---|---|
| **P1** | bomb deepening longitude ~ EMI + N34 + year | 47 | -1.09 | -2.44 to +0.34 | 1.97 (planted: 2.2) | 0.171 | 0.47 |
| **P2** | bomb deepening latitude, same model | 47 | -0.07 | -0.56 to +0.37 | 0.66 | 0.753 | 0.90 |
| S1 | longitude, rapid storms (B >= 12) | 47 | -0.12 | -1.41 to +0.86 | 1.56 | 0.852 | 0.94 |
| S2 | latitude, rapid storms | 47 | -0.03 | -0.41 to +0.35 | 0.54 | 0.892 | 0.94 |
| S3 | longitude, **SON-lagged** EMI and N34 (strictly before the storms) | 47 | -0.06 | -1.37 to +1.32 | 1.87 | 0.943 | 0.94 |
| S4 | latitude, SON-lagged | 47 | -0.07 | -0.63 to +0.28 | 0.65 | 0.748 | 0.90 |
| S5 | longitude of the lowest-pressure fix of bomb tracks | 47 | -1.63 | -3.20 to +0.09 | 2.31 | 0.069 | 0.42 |
| S6 | longitude, N4-N3 in place of EMI | 47 | -1.08 | -2.18 to -0.01 | 1.54 | 0.063 | 0.42 |
| S7 | longitude, no year term | 47 | -1.06 | -2.36 to +0.36 | 1.90 | 0.167 | 0.47 |
| S8 | longitude, no N34 term | 47 | +0.40 | -0.60 to +1.38 | 1.39 | 0.481 | 0.67 |
| S9 | B-weighted centroid longitude, rapid storms | 47 | -0.46 | -1.63 to +0.38 | 1.39 | 0.420 | 0.63 |
| S10 | longitude, JFM window | 47 | -0.80 | -2.23 to +1.27 | 2.43 | 0.330 | 0.59 |
| S11 | longitude, HF-reaching storms (pipeline A), 2004-05 on | 22 | -4.03 | -8.18 to +0.91 | 6.42 | 0.102 | 0.46 |
| S12 | El Nino winters: CP (4) minus EP (9), degrees (normal interval) | 13 | -5.25 | -9.88 to -0.63 | 6.60 | 0.019 | 0.35 |
| S13 | winter 75th percentile of bomb longitude | 47 | -1.60 | -4.14 to +1.79 | 4.13 | 0.250 | 0.56 |
| S14 | winter 25th percentile of bomb longitude | 47 | -0.79 | -2.16 to +0.94 | 2.24 | 0.308 | 0.59 |
| S15 | longitude, 1979-80..2003-04 (fresh, pressure-only) | 25 | -0.69 | -1.96 to +1.54 | 2.48 | 0.394 | 0.63 |
| S16 | longitude, 2004-05..2025-26 | 22 | -2.37 | -6.04 to +1.21 | 4.98 | 0.183 | 0.47 |

"Detectable at 80%" is 2.8 times the bootstrap SE (S12: the difference's SE). P1's planted-effect power (`power.csv`, 300 plants per size, null surrogate by permuting reduced-model residuals): 0.22 at 1 degree, 0.48 at 1.5, 0.70 at 2, 0.89 at 2.5, 0.96 at 3, 1.00 at 3.5. P2 is at 0.65 at 0.5 degree and 1.00 at 1; S3 tracks P1.

Mechanism (outside the family): adding DJF PNA leaves the EMI effect at -0.97 (-2.31 to +0.61); PNA's own coefficient is +0.42 (-0.45 to +1.38, p = 0.51). EMI correlates 0.27 with PNA, 0.73 with Nino3.4, 0.88 with SON EMI.

Decision-rule check: P1 and P2 are not "yes" (p > 0.05, q > 0.10); both meet "no, well-powered" (interval inside +-3 and +-2 degrees; 80% power at 2.2 and under 1 degree). S11 and S12 have intervals wider than the stated size, so they are "can't tell" by rule. S12 (p = 0.019, q = 0.35) says CP-El Nino storms deepen about 5 degrees west of EP ones, from 4 CP winters; a lead only.

## Power and what limits it

- The sample is 47 winters with one flavor value each; storms (about 33 bombs a winter) only sharpen each winter's mean. EMI and Nino3.4 correlate 0.73, so "flavor at fixed strength" costs precision (S8 drops the control: +0.40, still nothing).
- The HF-only outcome (S11, 22 winters, Decision 1) and the CP-versus-EP contrast (S12, 4 CP winters) are underpowered and say little on their own.
- 1979-2003 winters are fresh data for this pressure-only outcome; S15 uses them alone and agrees in sign and size with the full sample.

## What this does not show

- A proxy: deepening is from a tracker that smooths on about 0.5 degrees, 6-hourly, no land mask on the deepening point, tropical cyclones in. Atlantic not tested.
- Not a jet test: the idea is that EMI shifts the jet exit; only the storms' own deepening position is observed here.
- Pressure-only outcomes use 1979 onward; only S11 uses the gust label and is restricted to 2004-05 on (Decision 1).

## Looks

No model was fitted on one block and scored on another. All 47 winters, including 2015-2025, enter the regressions once: one further look at the 2015-2025 block (the 15th by the agenda's count) and one at pre-2001 seasons (the fifth; P1 and S15 use 1979-2000 pressure-only data). Logged in `research/era5/hemispheric/results/heldout_looks.log`.

## Verification

A fresh Sonnet agent that had not seen `analysis.py`, `deepening_points.py`, `results/` or this README recomputed from the committed `enso_kuroshio/data/` files with its own implementation (`verify/verify_numbers.py`, `verify/VERIFICATION.md`). After the fix under "Deviations" below, it matches: track and deepening-window counts (35,773), bombs per winter (33.15; min 21, max 49), rapid (70.3), HF-reaching (21.95); the EMI coefficient of P1, P2, S1-S4, S6-S11 and S13-S16 and the CP-minus-EP difference of S12, to three decimals; P1 permutation p (0.175 against 0.171) and bootstrap interval (-2.47 to +0.37 against -2.44 to +0.34); P2 p (0.752 against 0.753); the P1 planted-effect power at 2, 2.5 and 3 degrees (0.72, 0.88, 0.97 against 0.70, 0.89, 0.96); SD of winter-mean longitude (3.72 weighted) and latitude (1.01); the correlations of EMI with Nino3.4, PNA and SON EMI. S5 matched (-1.630) once the registered definition was restored (see below).

**Not independently checked:** permutation p-values and intervals of every test but P1 and P2, the BH q-values, the power of P2 and S3 and the P1 values at 0.5-1.5 and 3.5-6 degrees, the mechanism check M1, the S12 permutation p, and the detectable-effect figures (2.8 times the bootstrap SE).

## Post-hoc deviations (two code errors found by the verifier; both versions kept)

1. **Longitude midpoint across 0/360.** `deepening_points.py` first took the arithmetic mean of the two window longitudes, which puts Atlantic-to-Europe tracks crossing the Greenwich meridian at about 180E, inside the Pacific box (40 spurious bomb tracks, 1,598 against 1,558). Fixed with a circular midpoint. The first-run results are kept in `results/superseded_wrapbug/`. Effect: P1 -1.05 became -1.09, P2 +0.03 became -0.07; no conclusion or decision-rule outcome changed.
2. **S5 box.** The first run also required the minimum-pressure fix itself to lie in the Pacific box, which the plan does not say. Restored the registered definition (the bomb population is defined by the deepening point). S5 -1.71 (p 0.051) became -1.63 (p 0.070); both are above q = 0.10.

The planned analysis was otherwise unchanged. A 200-permutation smoke run preceded the reported run.

## Reproduce

    python3 research/era5/enso_deepening/deepening_points.py research/era5/enso_kuroshio/data/fixes.csv.gz research/era5/enso_deepening/data/points.csv.gz
    python3 research/era5/enso_deepening/analysis.py research/era5/enso_deepening/data/points.csv.gz research/era5/enso_kuroshio/data/pipelineA_match.csv research/era5/enso_kuroshio/data/winter_indices.csv research/era5/enso_kuroshio/data/sst_boxes_raw.csv research/era5/enso_deepening/results 10000 5000

Runs in about 30 seconds from committed files.
