# Do the two basins trade HF activity week to week? (RA-25 stage 1; archive and ERA5 proxy)

Plan, committed before any count correlation was computed: [PREREGISTRATION.md](PREREGISTRATION.md) (`3b7a49e`). Code: `run.py`
(16 pre-registered tests, lag profile, halves, coverage, power), `posthoc.py` (labelled post hoc). Results: [results/](results/).
Inputs: `hemispheric/results/weekly_table.csv.gz` and the PR 41 out-of-sample indices. ERA5 is a **proxy**. Pipeline A
(`research/era5/hf_history`, 800 km ocean gust index at 71.7 kt, depth cuts 966.2 hPa Atlantic and 965.0 hPa Pacific) supplies the two
proxy samples. Pipeline B is not used. Stage 2 (a sector index from Z500) needs the RA-2 pull and is not started.

## Answer

**Yes, but weakly: the weekly HF counts do trade off between the basins, at r of about -0.10, which is a sixth of the -0.60 seen between the
pattern indices.** The pre-registered prediction (r at most -0.15 in all three samples) was **not met**.

| sample (same-week r after removing the seasonal cycle, 660 weeks) | r | 95% season-block bootstrap | one-sided p | q (3 primary) |
|---|---|---|---|---|
| ARCH: archive counts, 2004-05 to 2025-26 | **-0.093** | -0.186, 0.000 | 0.018 | 0.018 |
| PRX: pipeline A gust-proxy counts, 2004-05 to 2025-26 | **-0.106** | -0.165, -0.044 | 0.007 | 0.018 |
| DEP: pipeline A depth-proxy counts, 1979-80 to 2000-01 (within-era only) | **-0.097** | -0.192, -0.010 | 0.015 | 0.018 |

- **By the pre-registered rules: rule 2, "partly / not at the predicted size".** All three primary tests pass FDR (3 of 3, q 0.018), but
  none reaches -0.15, and none meets rule 3 (a well-powered null), so neither "a count seesaw as predicted" nor "the seesaw is construction"
  is the right summary.
- **Relation to the index seesaw.** Over the same 660 weeks the two PR 41 indices correlate at -0.600; the archive count correlation is 0.155 of
  that (bootstrap 0.004 to 0.317). Most of the index seesaw is therefore not visible in the counts. Whether it is "construction" (both
  indices fitted to counts in the same fields) or only too noisy to see in counts is not decided by these numbers; see the post hoc block.
- **Power.** Planted-effect simulation (500 runs per point): a true count correlation of -0.15 is detected with 95% power in ARCH, 96% in
  PRX and 94% in DEP at one-sided p < 0.05 (90%, 91%, 87% at the Bonferroni-3 level). Minimum detectable r at 80% power: **-0.119 ARCH, -0.113 PRX,
  -0.118 DEP**. The observed r (about -0.10) is just under that, which is why the intervals touch zero in ARCH and sit near -0.15 at their lower end.
- **Pipeline note.** ARCH and PRX describe the same 22 winters, so they are not independent confirmations of each other; DEP, in a different
  era and from fields only, is the independent replication. The three agree to within 0.013.
- **Prediction scorecard.** r at most -0.15 in all three: failed in all three. Sign and significance as expected, size not.

## All 16 pre-registered tests

13 of 16 pass FDR across the whole family (q < 0.05); the three that do not are the interannual archive and proxy correlations and the
DEP subseasonal one. Full table: `results/tests.csv`.

| test | ARCH r (p) | PRX r (p) | DEP r (p) |
|---|---|---|---|
| primary: week-of-season means removed | -0.093 (0.018) | -0.106 (0.007) | -0.097 (0.015) |
| S1: month dummies removed instead | -0.083 (0.018) | -0.093 (0.006) | -0.083 (0.014) |
| S2a: season totals (n = 22 seasons) | +0.136 (0.73) | -0.278 (0.10) | **-0.629 (0.0004)** |
| S2b: each season's own mean also removed (within-season only) | -0.107 (0.008) | -0.098 (0.011) | -0.072 (0.055, q 0.063) |
| S3: NAO and PNA same-week removed | -0.087 (0.025) | -0.102 (0.008) | -0.096 (0.020) |
| S4: depth proxy in the fit era 2004-2025 | | | -0.127 (0.003) |

- **The weekly signal is subseasonal, not interannual.** Within seasons the correlation is similar in all samples (S2b -0.07 to -0.11). Season
  totals disagree: archive +0.14, gust proxy -0.28, depth proxy 1979-2000 -0.63 (q 0.007). With 22 seasons per sample and no replication of
  the DEP value in ARCH or PRX, I do not treat the interannual value as a finding.
- **It is not NAO against PNA.** Removing both indices leaves r unchanged to the second decimal (S3).
- **Depth versus gust is not the explanation.** The depth proxy in the fit era (S4, -0.127) is as negative as the gust proxy (-0.106).
- **Lag profile (S5, descriptive, `results/lag_profile.csv`).** Strongest at lag 0 for ARCH (-0.093, with -0.080 at one week);
  PRX and DEP have their largest values at +1 week (-0.129, -0.131: the Pacific week after the Atlantic one), against a null
  envelope of about +-0.09. Nothing at +-2 weeks beyond the envelope in ARCH. No test of any single lag was planned.
- **Halves (S6, descriptive).** ARCH -0.136 in 2004-2014 and -0.050 in 2015-2025; PRX -0.151 and -0.055. The count coupling is weaker in the
  later half in both; with 11 seasons per half I cannot say whether that is real.
- **Coverage (S7).** Weeks with both counts zero are 12-13% in every sample; dropping them makes r more negative (-0.194, -0.187, -0.155),
  not less, so shared quiet spells or coverage gaps do not create the correlation. The archive counts are not overdispersed after the
  seasonal cycle is removed (variance to mean 0.95 Atlantic, 0.87 Pacific).

## Post hoc (not in the plan; both views shown)

Weekly counts average 1.2 to 1.5 storms per week and their deseasonalised variance equals their mean (var/mean 0.87 to 1.13), so they are
close to pure Poisson noise around the seasonal mean. That limits any correlation of counts. To see whether the correlation grows when
the noise is averaged down, `posthoc.py` recomputed r on non-overlapping blocks of 1, 2, 3, 5, 6 and 10 weeks (all divisors of 30 are
shown, none picked; `results/posthoc.json`):

| block (weeks) | 1 | 2 | 3 | 5 | 6 | 10 |
|---|---|---|---|---|---|---|
| ARCH r | -0.093 | -0.103 | -0.134 | -0.187 | -0.117 | -0.265 |
| PRX r | -0.106 | -0.214 | -0.192 | -0.316 | -0.247 | -0.389 |
| DEP r | -0.097 | -0.187 | -0.190 | -0.308 | -0.227 | -0.469 |

The correlation grows with block length in all three, which is what a real but noisy coupling would do, and what a pure sampling artefact
would not. It is also what slower variability (interannual and seasonal-scale) would add at long blocks, so the longer blocks mix the
subseasonal coupling with the interannual one. The 5-week values (-0.19 to -0.32) are the nearest to the index scale. Block intervals are
wide (`results/posthoc.json`; e.g. ARCH 5 weeks -0.37 to -0.01) and the p values were not corrected. This is a lead on how noisy counts hide
a coupling, not a result.

## What this means for the basin results

- **The Atlantic and Pacific pattern results are not independent stories, but they are not one story either.** There is a real, small trade-off in
  the storms themselves (r about -0.10 per week, replicated in an independent era), and most of the -0.60 seen between the indices is not
  present in the counts. How much of the rest is a feature of how the indices were fitted, and how much is hidden by the noise of weekly
  counts, cannot be separated with counts alone.
- Stage 2 (a sector index that reproduces the seesaw) would separate the two, and is waiting on the RA-2 pull.

## What this does not show

A mechanism; that the true underlying rates are not much more strongly coupled than the counts (post hoc block table suggests they may be);
anything about fields (stage 2 not run); any pre-2001 gust-based level or trend (DEP is pressure-depth, within-era). Power is for correlations of
deseasonalised counts; small-count discreteness is not simulated.

## Looks spent

One further look at 2015-16 to 2025-26 (agenda numbering: the tenth; it is the eleventh line of `hemispheric/results/heldout_looks.log`) and the
fourth look at pre-2001 seasons. No model was fitted, so no hold-out was consumed in the fitting sense.

## Reproduce

`python3 research/era5/basin_counts/run.py` (about 2 minutes) then `python3 research/era5/basin_counts/posthoc.py`.

## Verification

See `VERIFICATION.md`.
