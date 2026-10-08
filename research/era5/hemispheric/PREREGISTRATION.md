# Is there a hemispheric state that leads to more hurricane-force lows? Pre-registration

Written and committed **before any outcome was joined to any predictor**. At that point only the field extraction
(`fields.py`, predictors only) had been run, and its overlap check against WeatherBench2 (a predictor-only check).
No archive count, no pipeline A count and no index had been looked at in this thread. What the project already knows from
other threads (NAO/PNA main effects, the ONI x MJO null, the Q1 and Q2 results) is background and is listed in "What was
already known".

Question (Jason, 2026-10-08): *is there a state of the atmosphere across the hemisphere that leads to more HF lows,
whatever the combination?* Data-driven: the pattern is found from the fields, not from named indices. Then: is it just a
known index combination, or something else?

Everything below is fixed in advance. Anything changed later is logged under "Deviations" as post hoc, and both versions
are reported. "No pattern beats the named indices out of sample" is a full answer.

## Data

**Outcome, primary: the archive** (`docs/data/hf-lows.json`, rows with `cls == "low"`, the same file and the same rule the
additive teleconnection test used). Weekly counts per basin, Atlantic and Pacific separately, weeks of 7 days from 1 October,
30 weeks per season (Oct-Apr), by the date in the `start` field. Seasons 2004-05 to 2025-26 (season label = starting year
2004..2025), 22 seasons, 660 weeks per basin. Seasons before 2004-05 are never used to fit or to test (decision 1).

**Outcome, proxy (secondary checks): ERA5 pipeline A** (`research/era5/hf_history/results/all_tracks.csv.gz`, 800 km ocean
gust index, HF-equivalent at 71.7 kt; a proxy). Counted by the track's `start` date into the same weeks. For the pre-2004
within-era check only, the depth version: tracks with `minp` at or below 966.2 hPa (Atlantic) or 965.0 hPa (Pacific), the
count-matched cuts fixed in `freq_split` (hf-low PR 14), not re-chosen here.

**Predictor fields** (ERA5 proxy of the atmosphere; reanalysis, not observations): 12 UTC daily Z500 (m), U250 (m/s, the
jet), MSLP (hPa) and SST (K) on WeatherBench2's 64 x 32 conservative grid (5.625 degrees); 12 rows of cell centres
25.3 to 87.2 N (the cells that overlap 20-90 N), 64 longitudes. Source: WeatherBench2 to 2023-01-09; ARCO-ERA5 hourly
0.25 degrees regridded onto the same cells by area-weighted overlap from 2023-01-10 (checked against WeatherBench2 on 12
overlap days: rms difference 0.001 m in Z500, `overlap_check.txt`). Volume about 10 GB, below the 50 GB gate. Seasons
1979-2025, days 5 Sep to 10 May.

**Anomalies.** Per grid point and field: minus the calendar-day mean over 1991-2020, smoothed with a 15-day running mean
(Feb 29 uses the mean of Feb 28 and Mar 1), then minus a linear trend in the Oct-Apr seasonal mean fitted over seasons
1979-2021 and extrapolated to 2022-2025 (removes warming drift; the trend is fitted on predictors only). SST uses ocean cells
that are valid on every day.

**Named indices** (the comparison model): NAO and PNA (CPC daily, `/mnt/project-files/teleconnection-test/cpc_indices/`),
ONI (monthly, at the month of the window midpoint), and the MJO as the two-coordinate longitude projection used by the
additive test (MJO1, MJO2 at the pentad of the window midpoint). AO is used only in one attribution sensitivity. The loader
is the additive test's (`extend_from_cpc`), copied with its conventions. Each index is standardised with the discovery-period
mean and SD.

## Units and windows

Week `w` of season `s` starts on `S`. The predictor for a week is the mean of the daily anomaly fields over **days S-7 to
S-1** (the seven days before the week starts, primary lag; so the field leads the week's onsets by 1 to 13 days, 4 days
on average). This is the guard against the circularity found in Q1 and Q2, where storms feed the indices. Secondary lag:
days S-14 to S-8. The concurrent window (S to S+6) is run once as a labelled contrast to show how much circularity there is
and is not tested or used for inference.

## Split

- **Discovery**: seasons 2004-05 to 2014-15 (11 seasons, 330 weeks per basin).
- **Held-out test**: seasons 2015-16 to 2025-26 (11 seasons, 330 weeks per basin). Looked at **once** for the primary
  analysis, after everything below is frozen, and the number of looks is stated with the result.
- Within discovery, leave-one-season-out cross-validation chooses the penalty. That is allowed; the held-out seasons never
  enter it.

## Primary method: a ridge-penalised Poisson pattern model on field EOFs

1. EOFs: for each of Z500, U250, SST, the leading **K = 10** principal components of the discovery weekly anomaly maps (cos
   latitude square-root weighted), each score scaled to unit variance in discovery. Held-out weeks are projected onto the
   discovery EOFs. 30 predictors. EOFs never see the outcome.
2. Model per basin (log link): month fixed effects (Oct..Apr, by the month of the week's midpoint) + `log(1 + count in the
   previous 7 days, same basin)` + the 30 PC scores with a ridge penalty `lambda * sum(beta^2) / 2`. The month effects and
   the previous-week term are not penalised. `lambda` is picked from {1, 3, 10, 30, 100, 300, 1000, 3000, no PCs} by
   leave-one-season-out deviance within discovery. "No PCs" is allowed to win, which is the null result.
3. The previous-week count is in every model on purpose: storms in the previous week imprint on the previous week's fields,
   so without it the pattern could be only storm clustering. The model that omits it (month only) is reported as a
   secondary contrast.
4. The **pattern** shown on maps is the model's linear predictor mapped back to the grid, `sum_j beta_j EOF_j`, per field.

**Comparison models**, all fitted on discovery only:

- `B0` month effects only; `B1` month + previous-week count (the baseline).
- `N` the named-index model: `B1` + NAO + PNA + ONI + MJO1 + MJO2, main effects, unpenalised. `N+I`: `N` + NAO x PNA, NAO x
  ONI, PNA x ONI (ONI x PNA was untested before).
- `P` the pattern model (`B1` + penalised PCs). `P+N`: `N` + penalised PCs, with its own `lambda` chosen the same way.

**Skill** on the held-out weeks: Poisson deviance `D`, skill score `SS(X|Y) = 1 - D_X / D_Y`.

## Secondary method: self-organising-map regimes

A 3 x 3 SOM (batch updates, 200 epochs, neighbourhood radius from 1.5 to 0.5, initialised on the first two PCs of the
training data, seed 20261008) trained on the 10 Z500 PC scores (unit variance) of the discovery weeks. Each held-out week
takes the best-matching node. Model: `B1` + node effects (unpenalised). Report each regime's HF rate ratio in discovery and
held-out, and `SS(SOM|B1)`. The regime is a hemispheric Z500 state with no named index in it.

## Tests

Per basin, on held-out weeks. One-sided permutation p-values, 10,000 permutations, **season-block**: the PC rows of held-out
season j are swapped with those of a randomly chosen season (week index kept, so the seasonal cycle and the previous-week
count stay aligned with the outcomes). Seed 20261008.

- **P1** `SS(P|B1) > 0`: does the field pattern predict HF counts beyond the month and the previous week?
- **P3** `SS(P|B1) - SS(N|B1)`, 95% interval from a 10,000-draw bootstrap over held-out seasons; one-sided p = share of
  draws at or below 0: does the pattern do better than the named indices?
- **P4** `SS(P+N|N) > 0`: does the pattern carry information the named indices do not (PC rows permuted, indices kept)?
- Informational: `SS(N|B1)` with its own permutation p (index rows permuted), `SS(P|B0)`.

**Family 1** = P1, P3, P4 in both basins (6 tests): Benjamini-Hochberg q across the six. **Family 2** = everything
pre-registered below (Family 1 + secondary tests): BH across all of them. Both are reported for every test.

Pre-registered secondary tests (each in both basins):

- **S1** swap split: discovery 2015-16 to 2025-26, test 2004-05 to 2014-15 (P1, P3, P4). Uses the same seasons again;
  labelled as a replication, not an independent test.
- **S2** lag 2 (days S-14..S-8), primary split (P1, P3, P4).
- **S3** SOM regimes, primary split (`SS(SOM|B1)`).
- **S4** proxy outcome: the primary-fitted model applied to pipeline A HF-equivalent counts in the held-out seasons (P1
  analogue; the model is the archive-fitted one, not refitted).
- **S5** within-era 1979-80 to 2000-01 with pipeline A *depth* counts (decision 1: variation within the era only): month +
  previous-week + one scalar `gamma` times the frozen primary pattern index, refitted on that era as nuisance; test is the
  one-sided deviance gain for `gamma > 0` with season-block permutation of the index. Only per-season relative variation is
  used; no level or trend is read from it.

## Decision rules (fixed now)

Per basin:

- **Does a hemispheric state exist that leads to more HF lows?** *Yes* if P1 has Family 1 q < 0.05. *No (well-powered
  null)* if q >= 0.05 and the held-out data had 80% power at a pattern effect of RR 1.15 per SD of the pattern index or
  smaller. Otherwise *can't tell*.
- **Known combination or something new?** Only if *yes*: *known* if P4 q >= 0.05 and the attribution R-squared below is
  at least 0.5; *partly new* if P4 q < 0.05; otherwise *unclear*. P3 states whether it beats the named model.
- The replications (S1, S2, S4, S5) must agree in sign with the primary for a *yes* to be called robust; a *yes* whose swap
  replication (S1) fails is reported as "not replicated".

## Power

Planted-effect simulation, 2,000 draws per level: held-out counts are replaced by Poisson draws from the `B1` mean times
`RR^z` with `z` the standardised held-out pattern index (the penalised-PC part of the discovery-fitted linear predictor),
RR in {1.05, 1.10, 1.15, 1.20, 1.30, 1.50}. Power = share of draws with P1 statistic above the 95th percentile of the
permutation null. The smallest RR with 80% power is the detectable effect. The effective sample is 11 seasons; weekly
values within a season are autocorrelated, so the permutation is by season and not by week.

## Attribution

On all 22 seasons (660 weeks per basin), the *out-of-sample* pattern index: for each left-out season, the pattern model is
refitted on the other 21 (EOFs and `lambda` as in the primary split, rechosen by inner leave-one-season-out), and the left-out
season's penalised-PC linear predictor is its index. The index is regressed (OLS, month effects) on NAO, PNA, ONI, MJO1 and
MJO2 (A1), plus the three products (A2), plus AO (A3). Report the adjusted R-squared with a 95% season-block bootstrap
interval, standardised coefficients, and the share of held-out *skill* that the named model carries (`SS(N|B1)` against
`SS(P|B1)`). "Known" = R-squared of A2; "residual" = 1 minus that. Also reported: the map correlation of the discovery and
held-out regression maps (Z500, month-adjusted, B1-residual), as a replication diagnostic.

## Maps (descriptive, after all tests)

Z500 / U250 / SST pattern maps per basin from the model refitted on all 22 seasons, and the Z500 composite difference between
the highest and lowest quintile of the index, hatched where the field-significance FDR (BH across grid points, season-block
permutation p) holds. Maps are descriptive; they are not a test.

## What this does not do

- It does not claim a cause. A pattern that leads HF counts could be a state that favours cyclogenesis, one that steers
  storms into the basin, or one the archive's own recording favours.
- 22 seasons is a small sample for a hemispheric search: 11 per side. A null here is a null at the stated power, not a
  proof of absence.
- Pipeline A is a proxy and pre-2004 gust counts are not used for levels or trends (decision 1); S5 uses the depth version
  and variation within the era only.
- Coarse 5.6 degree fields resolve planetary and synoptic-scale patterns; they cannot resolve jet-exit details or fronts.

## What was already known (before this thread)

Weekly Poisson tests (REPORT.md in `teleconnection-test`): Atlantic NAO RR 1.128 per SD, Pacific PNA 1.217, ONI and MJO no
main effect, no robust interaction, a vanished ONI x MJO hit. Q2 (hf-low PR 14) halved the published effects once the index
was lagged. Q1: NAO/PNA shift where storms are and the background pressure, not intensity. These shaped the choice of lag and
the previous-week covariate; they do not enter the pattern search.

## Deviations (post hoc; filled in as they happen)

(none yet)
