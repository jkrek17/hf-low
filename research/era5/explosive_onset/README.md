# What turns sub-55 kt storms into HF lows within a day? (RA-20, Arms A and B)

Plan: [PREREGISTRATION.md](PREREGISTRATION.md), committed before any field was read (`fe1f63d`); features and per-group power committed
before any test (`9f7aae3`). **ERA5 proxy, pipeline A** (`research/era5/hf_history`, 800 km ocean gust index, HF at 71.7 kt). Seasons
2004-05 on. Transitioning tropical cyclones are in. Baseline: the PR 76 boosted-tree LOSO probabilities (M1; M0 = PR 12 logistic).

## Answer in plain words

**No usable gain. Neither the 0.25 degree surface structure (Arm B) nor the fields aloft (Arm A) improves the 24 h forecast for storms
below 55 kt beyond the tree model, and the +0.005 Brier-skill bar is not met by any group.** Arm B has a small, borderline signal that
separates cases from same-time controls; Arm A has none that survives multiplicity. Neither is large enough to explain the missed
explosive onsets.

| | Result (stratum S: 135,663 fixes below 55 kt and not yet HF; 927 cases) |
|---|---|
| Arm B, B1 gust-field concentration | separates cases from same-time controls: LR 12.6 (3 df), p 0.0056, q 0.028 within the arm, **q 0.050 over all nine primary tests**. Per SD, odds of conversion: gust maximum farther from the centre +0.155 log-odds (se 0.053), smaller area of 45 kt gusts -0.109 (0.046). Forecast gain DeltaBSS_S -0.00004 (90% interval -0.0039 to +0.0036) |
| Arm B, B4 shear and convergence | LR 8.4 (2 df), p 0.015, q 0.029 within, 0.052 overall. Weaker 10 m wind shear -0.193 (0.065), stronger low-level convergence +0.125 (0.066). DeltaBSS_S -0.0020 (-0.0062 to +0.0025) |
| Arm B, all 11 features | LR 23.1 (11 df), p 0.017, q 0.029 within, 0.052 overall. DeltaBSS_S -0.0021 (-0.0094 to +0.0048) |
| Arm B, B2 pressure structure, B3 wind vs gust | no (p 0.17 and 0.24; MDE for the near-collinear features is 0.2 to 0.7, so "can't tell" for these) |
| Arm A, omega700 (warm-conveyor ascent) | p 0.033, q 0.13 within the arm, 0.075 overall. Stronger ascent raises the odds, as predicted: -0.131 log-odds per SD of omega (se 0.059), about 1.14 per SD, below the 0.166 MDE. DeltaBSS_S -0.0057 (-0.0142 to +0.0018): **the +0.005 prediction is excluded** (upper 95% bound +0.0018) |
| Arm A, trough depth, stability | p 0.21 and 0.87; DeltaBSS_S -0.0012 and -0.0029; MDE 0.165 and 0.227 per SD, so "can't tell" by the plan's rule (0.15 or less is needed for "no") |
| Against the logistic M0 instead | 8 of 9 groups separate cases (q < 0.05 in the secondary family), DeltaBSS_S +0.002 to +0.010 with intervals that cross 0. So the features carry the kind of curvature PR 76's trees already absorb; that reading is a hypothesis, not tested here |

Decision-rule outcome: **H_B** (at least one B group passes): passes within arm, borderline overall, no forecast gain. **H_A**
(omega700 adds +0.005 BSS): **not supported**. Reading: the one-day explosive onsets from weak storms are not recoverable from the
resolved surface structure or the Tier 2 fields tested here; what separation exists is about 1.1 to 1.2 per SD and does not help a
Brier score. The Tier 2 pull for the rest of the record is not justified by this result.

## Tests and multiplicity

- 64 tests were run (9 primary within-time vs M1, 9 vs M0, 18 held-out and 18 BSS sign-flip tests, 10 subset tests). BH over all 64:
  12 pass q < 0.05; of them only **B1 vs M1 (q 0.030)** is a primary test against the tree baseline. The others are within-time tests
  against the logistic M0 (7), the held-out gain for B1 vs M0 and B_all vs M0, and `stab` vs M1, whose significant held-out gain is
  **negative** (the feature hurts out of sample). `results/summary.txt`, `results/all_tests_fdr.csv`.
- Primary family: 3 of 9 pass within arm (all Arm B), 0 of 9 over all nine.
- Subsets (descriptive, the "all" model): Atlantic p 0.13, Pacific 0.18, storms with a 12 h pressure change p 0.018 (Arm B) and 0.040
  (Arm A), new storms without one 0.74 and 0.87, the 198 matched-control cases p 0.12.
- Tree check (fixed settings): ΔBSS -0.19 (Arm B) and -0.18 (Arm A). It overfits at this sample size and weighting; uninformative.

## Power

Minimum detectable effect at 80% power, per feature, adjusted for the others in its group (log-odds per SD, `power_groups_*.txt`):
B1 0.16 to 0.17, B4 0.19 and 0.18, B3 `GF` 0.14 but `WS99` 0.39; B2 `SHARP`, `CURV` 0.55 (near-collinear); Arm A 0.17 (trough,
omega700), 0.23 (stab). The within-time test has 79% power at a 0.10 SD shift for a single independent feature. The BSS sign-flip has
only about 7% power at the +0.005 bar (`power_bss.txt`), which is why DeltaBSS_S is an effect size, not the decision test.

## What was run

Arm B: 1,495 times, 17.93 GB streamed (sized 17.94); Arm A: 1,007 WeatherBench2 chunks, 26.96 GB (sized 26.96). 895 case times plus 600
random other times (weight 25.29). Arm A covers times to 2023-01-09 only: 740 of 895 case times, 766 of 927 cases, 19 of 22 seasons;
the 155 later case times would need about 72 GB from the 0.25 degree store and were not pulled. The re-extracted gust matched the table's
`g800` within 0.5 kt for 100% of 9,266 fixes with an in-domain value (3,127 fixes have a table value of 0 and are outside the check).

## Limits

- Features come from the same ERA5 gust field that defines the label, at the same time, so Arm B asks what that field's own structure adds.
  ERA5 under-resolves the strongest winds and the melting layer; a null is not proof the real storm has no such structure.
- Arm A is 1.5 degree and cannot see mesoscale ascent; it covers 19 seasons.
- The seasons were used by PR 12, 56, 68 and 76. This spends two looks at 2015-25 (logged as #15 and #16 by the agenda count).
- Fixes within a storm are not independent: the p-value is the larger of the model-based LR and a season bootstrap, and still treats risk
  sets as independent. The 12 h pressure change is missing for 42% of the cases (new storms); there the features show nothing.
- Which q decides was fixed after the results (see Deviations in the plan); both are above.

## Verification

A fresh Sonnet agent, without reading the plan, README or result files, recomputed from the committed inputs and agreed within rounding on: the
counts (135,663 / 927 / 531 / 396 / 803 tracks / 895 case times / 15,173 other times / 7,296 fixes), the LR statistics and p for B1, B2, B3, B4,
the 11-feature model, omega700, trough and stab, the B1 coefficients, and DeltaBSS_S for the 11-feature model (-0.0021), B1 (-0.00004) and
omega700 (-0.0057). **Not independently checked:** the season-bootstrap p-values and the larger-of rule, every q value and the 64-test
accounting, the held-out log-likelihood gains and sign-flip p, the DeltaBSS intervals, the subset tests, the tree check, all power runs and
MDEs, the pull sizes (HEAD-based; the streamed bytes matched them, 17.93 and 26.96 GB).

## Reproduce

    python3 -I research/era5/explosive_onset/select_times.py .          # times, labels only
    python3 -I research/era5/explosive_onset/size_pull.py .             # HEAD requests, exact sizes
    ERA5_WORK=/path python3 research/era5/explosive_onset/extract_armB.py . 8    # 17.9 GB
    ERA5_WORK=/path python3 research/era5/explosive_onset/extract_armA.py . 4    # 27.0 GB
    ERA5_WORK=/path python3 -I research/era5/explosive_onset/analyse.py . --power BA
    ERA5_WORK=/path python3 -I research/era5/explosive_onset/analyse.py . --run B   # then --run A
    python3 -I research/era5/explosive_onset/summarise.py .

`results/features_armB.csv.gz` and `features_armA.csv.gz` hold the extracted feature tables (0.5 MB and 0.3 MB); raw fields are not committed.
