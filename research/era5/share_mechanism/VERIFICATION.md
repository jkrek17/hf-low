# Independent verification of the share-mechanism numbers

Own code (not importing part1.py/part2.py), run with `python3 -I`; `chanlib.py` used only for loading tracks, the OOS index, month dummies,
the weekly counts and the Poisson IRLS. Inputs: committed `all_tracks.csv.gz`, `fixes_2004.csv.gz`, `env_2004.csv.gz`, `oos_index_{atl,pac}.csv`.
Intervals: 300-draw season-block bootstrap, percentile (the README's are 2,000 draws, so small differences at the interval ends are expected).
Scripts: `/tmp/claude-0/verify/v1.py` (Part 1), `v2.py` (Part 2), `v3.py` (window check). Point estimates match to the quoted precision everywhere.

## Part 1 (weekly design, rate ratio per SD)

| Quantity | Quoted | Mine | Match |
|---|---|---|---|
| Atl A4 share (PR 64 T-values) | 1.271 (1.19-1.35) | 1.271 (1.210-1.353) | match |
| Pac A4 share | 1.129 (1.06-1.21) | 1.129 (1.063-1.215) | match |
| Atl A1 P(deepening) | 1.009 (0.988-1.029) | 1.009 (0.987-1.029) | match |
| Pac A1 | 0.982 (0.964-1.002) | 0.982 (0.964-1.000) | match |
| Atl A2 P(HF given deepening) | 1.252 (1.189-1.331) | 1.252 (1.197-1.339) | match |
| Pac A2 | 1.137 (1.074-1.223) | 1.137 (1.073-1.230) | match |
| Atl A3 residual | 1.006 (0.98-1.03) | 1.006 (0.981-1.034) | match |
| Pac A3 | 1.011 (0.995-1.028) | 1.011 (0.996-1.027) | match |
| Atl f_conv | 0.94 (0.83-1.06) | 0.936 (0.841-1.066) | match |
| Pac f_conv | 1.06 (0.87-1.31) | 1.063 (0.883-1.375) | match (upper end differs with 300 vs 2,000 draws; the ratio is unstable when A4 is small) |
| Atl f_deep* | 0.037 (-0.056 to 0.113) | 0.037 (-0.061 to 0.117) | match |
| Pac f_deep* | -0.153 (-0.451 to 0.013) | -0.153 (-0.519 to 0.000) | match |
| Atl f_res (sum check) | 0.027 | 0.027 | match |
| Pac f_res | 0.090 | 0.090 | match |

Weekly totals I get (window Oct 1 + 210 d): Atlantic N 9,656, D 4,521, HD 901, H 1,007; Pacific N 10,040, D 4,947, HD 788, H 825.
A4 reproducing PR 64 (1.271, 1.129): **confirmed**, and A1+A2+A3 = A4 exactly in log terms.

**143 HF cyclones without a deepening fix: confirmed**, 106 Atlantic + 37 Pacific (HF tracks in the weekly window with no 00/12 fix having dp12 <= -3.6).
Note: the PREREGISTRATION quotes Pacific HF-with-deepening 789; the weekly window contains 788. The one missing track starts 2022-04-29 06 UTC (day 210, outside
the 210-day window), so the difference is the window edge, not an error.

## Part 2 (deepening cyclones, reference fix 12 h before the first deepening fix)

First deepening fixes 14,043 (all seasons, both basins); 12,791 have an earlier fix in the track, 12,641 of them exactly 12 h earlier (150 have a longer gap and are
excluded by the "00/12 fix immediately before" reading; using the previous row regardless of gap was not needed, since the sample sizes below match exactly).

| Quantity | Quoted | Mine | Match |
|---|---|---|---|
| Atl cyclones / HF (complete case) | 3,285 / 754 | 3,285 / 754 | match |
| Pac cyclones / HF | 4,247 / 691 | 4,247 / 691 | match |
| Atl beta_tot* | 0.343 (0.26-0.43) | 0.343 (0.255-0.420) | match |
| Pac beta_tot* | 0.200 (0.12-0.30) | 0.200 (0.131-0.284) | match |
| Atl M* (g-formula, 8 ingredients) | 0.31 (0.19-0.47) | 0.312 (0.191-0.461) | match |
| Pac M* | 0.27 (0.02-0.72) | 0.273 (0.025-0.692) | match |
| Atl difference method | 0.23 (0.08-0.42) | 0.233 (0.069-0.401) | match |
| Pac difference method | 0.28 (-0.07 to 0.77) | 0.275 (-0.068 to 0.750) | match |
| Atl OR jet250 (b-path) | 2.05 (1.88-2.26) | 2.050 (1.870-2.250) | match |
| Atl OR eady | 1.96 (1.81-2.12) | 1.957 (1.797-2.129) | match |
| Pac OR jet250 | 2.25 (2.06-2.50) | 2.245 (2.069-2.509) | match |
| Pac OR eady | 1.63 (1.50-1.80) | 1.634 (1.512-1.785) | match |

My implementation: z standardised within basin on the analysis sample; a_j by OLS of z_j on month dummies + x; outcome model logit(month + 8 z + x) with ridge 1e-6;
NIE = mean[p(z+a, x) - p(z, x)], NDE = mean[p(z, x+1) - p(z, x)], M = NIE/(NIE+NDE); difference method 1 - beta_dir/beta_tot on the same sample.
Bootstrap draws with beta_tot < 0.01: 0 of 300 in both basins (the "unstable" flag is not triggered).

## Discrepancies

None beyond rounding and the 300-vs-2,000-draw interval noise. The Pacific f_conv upper end (1.375 vs 1.31) is the only interval end more than 0.05 apart; it is
a ratio with a small denominator and the interval is compatible.

## Not checked

- Permutation p-values and BH q-values (q 0.0015, etc.), a-path and b-path tests other than jet250/eady, the tests count (28 of 40).
- All secondary variants (S1-S5 of both parts, S3 other indices, halves), the single-ingredient and unique M_j, the secondary-ingredient and no-month variants.
- Minimum detectable effects, `power_m.py` results, the "well-powered null" classification (only the intervals of A1 were compared; both lie within 0.95-1.05 in my run
  apart from the Pacific lower-end 0.964 and upper 1.000, which are inside).
- The verdict wording against the decision rules (I confirm the numbers they rest on, not the rule application), and the heldout_looks.log entry.
- The pre-registration timing claim (commit order) and the statement that no 1979-2000 replication was possible.
