# Do the two basins trade HF activity week to week? Pre-registration (RA-25, stage 1)

Written and committed **before any cross-basin count correlation was computed by this thread**.
What had been looked at: the column names and row count of `hemispheric/results/weekly_table.csv.gz` (1,410 rows = 47 seasons x 30 weeks),
the RA-5 plan (`pacific_lead/PREREGISTRATION.md`), and the RA-25 agenda entry. Nothing else. RA-5 (PR 81) used the Pacific count only as a
control inside regressions (tests P5, P6) and never reported the raw weekly count-to-count correlation; it reported the index-to-index
correlation (same week r = -0.60, 2004-2025). That is the number this study is trying to interpret.

ERA5 is a **proxy**. The archive (OPC warnings) is the observation. Pipeline A (`research/era5/hf_history`, 800 km ocean gust index, 71.7 kt;
depth cuts 966.2 hPa Atlantic and 965.0 hPa Pacific as fixed in `freq_split`) is used for the two proxy outcomes and is named when it appears.
Pipeline B is not used.

## Question (agenda RA-25, stage 1 only)

The two basin pattern indices of PR 41 anticorrelate at -0.60 in the same week. Both were fitted to counts in the same fields, so some of
that could be built in. Are the **weekly HF counts themselves** anticorrelated between the basins once the seasonal cycle is removed?

Answer given in plain words: **yes** (a count seesaw, by how much), **no** (well powered: the index seesaw is construction and the two
basin results should not be read as one story), or **can't tell**.

Prediction (written before data): r <= -0.15 in all three samples below. I have no evidence for it; a near-zero r is plausible because
counts are small integers and the two basins' storm tracks are driven by partly independent upstream states.

Stage 2 (a sector index from Z500 fields) needs the RA-2 pull and is **not started**.

## Samples (all from `weekly_table.csv.gz`; 30 weeks of 7 days from 1 October, season label = starting year)

| id | outcome | seasons | pipeline | note |
|---|---|---|---|---|
| ARCH | archive weekly HF-low counts (`y_atl`, `y_pac`) | 2004-05 to 2025-26 (22) | observation | 2004-05 on only, because the archive was short-counted before |
| PRX | pipeline A gust-equivalent weekly counts (`pAhf_*`) | 2004-05 to 2025-26 (22) | A, proxy | fit/test era (decision 1) |
| DEP | pipeline A pressure-depth weekly counts (`pAdepth_*`) | 1979-80 to 2000-01 (22) | A, proxy | within-era only (decision 1); no level or trend is claimed |

The pre-2001 gust counts are **not** used. DEP is the pressure-depth count, permitted back to 1979. Tracks are counted by start date in
the week (as the table was built).

## Primary test (one per sample, three in all)

1. Deseasonalise: subtract, per basin and per week-of-season (0..29), the mean over the 22 seasons of the sample. Counts are not transformed.
2. Statistic: Pearson r between the Atlantic and Pacific deseasonalised weekly counts over all 660 weeks, same week (no lag; counts are
   outcomes of the same period, so the question of a precursor does not arise).
3. Inference: one-sided (negative) **season-block permutation**: whole 30-week Pacific seasons are shuffled among seasons (20,000
   shuffles; p = (1 + #{r_perm <= r_obs}) / (1 + 20,000)). The unit is the season, so weekly autocorrelation is respected and the effective
   n is the 22 seasons. 95% interval: season-block bootstrap (resample seasons with replacement, 5,000 draws, percentile).
4. Report r, interval, p, BH q (family 1: the three primary tests).

## Secondary tests (labelled; all enter family 2)

- **S1 alternative deseasonalisation.** Month dummies (7 months, as in `chanlib`) removed by OLS instead of week means. 3 tests (ARCH, PRX, DEP).
- **S2a interannual.** Correlation of season totals (30-week sums) between the basins, n = 22 per sample, exact season permutation. 3 tests.
- **S2b subseasonal.** Week-of-season means **and** each season's own mean removed, so only within-season weekly co-variation remains. 3 tests.
- **S3 named indices removed.** Both deseasonalised counts residualised on same-week NAO and PNA (columns `nao`, `pna` of the table) before the
  correlation. If the count seesaw is just NAO against PNA it vanishes here. 3 tests.
- **S4 depth in the fit era.** DEP's definition on 2004-05 to 2025-26 (`pAdepth_*`), separating "depth versus gust" from "era". 1 test.
- **S5 lag profile (descriptive, untested).** r between Atlantic week w and Pacific week w+k, k = -2..+2, within season, with the permutation 95% envelope.
- **S6 halves (descriptive).** The primary r for ARCH and PRX in 2004-2014 and 2015-2025 separately.
- **S7 coverage.** Share of weeks with both counts zero, per sample, and the primary r after dropping those weeks (descriptive; a shared quiet
  spell can create correlation, and a shared coverage gap would too).
- **S8 relation to the index seesaw.** The same-week r of the two out-of-sample pattern indices of PR 41 (`hemispheric/results/oos_index_{atl,pac}.csv`)
  over the same 660 weeks, and the ratio r_counts(ARCH) / r_index with a season-block bootstrap interval (descriptive).
- A **positive** primary r is reported with a two-sided p from the same permutation distribution (|r|), labelled as a non-predicted sign.

Family 1 = the three primary tests. Family 2 = primary + S1 + S2a + S2b + S3 + S4 = 16 tests. Both q values are reported for each, and the
number passing is stated as k of N.

## Decision rules (thresholds fixed now)

1. **Yes, a count seesaw.** In ARCH, PRX and DEP the primary r is <= -0.15 and family 1 q < 0.05. The index seesaw is then at least partly in the
   storms; its size relative to the index r is r_counts / r_index (S8).
2. **Partly / not replicated.** At least one sample has family 1 q < 0.05 with r < 0, but rule 1 is not met. Report which samples, and the answer is
   "weaker or less consistent than the index seesaw".
3. **No (well powered).** No primary test has q < 0.05 in the negative direction, **and** the planted-effect power at rho = -0.15 (one-sided alpha 0.05)
   is at least 80% in the sample, **and** the 95% bootstrap interval of r excludes -0.15 (lower bound above -0.15). Required in the archive and
   at least one proxy sample for the overall answer "no, the index seesaw is construction".
4. **Can't tell.** Everything else. The minimum detectable correlation (80% power) is stated for each sample with the answer.

## Power

Planted-effect simulation per sample: P* = rho * Z_A + sqrt(1 - rho^2) * Z_P(seasons permuted), where Z are the standardised deseasonalised
counts, so the planted correlation is rho and the Pacific keeps its seasonal-block autocorrelation; the primary permutation test (1,000
shuffles inside each simulation) is then applied. 500 simulations per grid point, rho in {-0.05, -0.10, -0.15, -0.20, -0.25}; power at
one-sided p < 0.05 and at p < 0.05/3. The minimum detectable r is read off the grid by linear interpolation. The planted correlation is a
correlation of deseasonalised counts, not a Poisson rate ratio; small-count discreteness is not simulated (a limitation, stated in the README).

## Looks spent

No model is fitted, but the 2015-16 to 2025-26 weeks are inside ARCH and PRX and the 1979-2000 weeks inside DEP. So this plan counts **one further
look at 2015..2025** (agenda numbering: the tenth; the log has ten lines, this will be the eleventh) and **one further look at pre-2001
seasons** (the fourth). The entry is appended to `hemispheric/results/heldout_looks.log` when the results are written.

## What this study cannot show

- A mechanism. A count correlation compatible with a hemispheric state feeding one basin at the other's expense, with a common modulator of
  the archive's recording, or with one season's overall level.
- Anything about fields: stage 2 (sector index) is not run.
- Counts of a few per week have little dynamic range; a real but weak coupling (|r| < the MDE stated with the answer) cannot be excluded.
- DEP says nothing about the archive before 2004 (short-counted) and the proxy before 2001 is used within-era only.

## Deviations (post hoc)

None yet. Logged here, with both versions shown, if any are needed.
