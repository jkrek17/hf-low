# Do hurricane-force lows cluster in time and along the same track? (pre-registration)

Written 2026-10-08, before any statistic in this directory was computed. Climatology agenda
item 8 (`/mnt/project-files/climatology-agenda/hf-climatology-agenda.md`); ledger Future-task 9.
Anything that changes after the first computed result is logged in `README.md` under
"Departures from this plan", with the reason, and both versions are reported.

**What has been seen already.** Only schemas and counts of the data below, and the agenda's raw archive
variance-to-mean ratios (Oct-Mar, seasonal cycle *not* removed, so upper bounds): weekly 1.14 (Atlantic) and
1.04 (Pacific), monthly 1.76 and 1.64. No statistic with the nulls below has been computed. Nothing in this
plan was tuned to a result.

Archive numbers are observations of OPC warnings. Pipeline A numbers are an **ERA5 proxy** (800 km ocean
gust index, HF-equivalent at 71.7 kt; `research/era5/hf_history`). Pipeline B is not used.

## Questions

- **Q-time.** Do HF lows arrive in bursts, more than a Poisson process with the known seasonal cycle allows? Per basin.
- **Q-space.** Is the bunching also spatial, i.e. do lows that arrive close in time also sit close together, or follow the same track?
- **Q-cond.** How much of any excess dispersion is accounted for by large-scale state (NAO for the Atlantic, PNA for the Pacific, lagged)? How much remains?

## Data (fixed now)

| Tier | Source | Events | Event time | Seasons |
|---|---|---|---|---|
| 1 ARCH (primary) | `docs/data/hf-lows.json`, class `low`, at least one fix with category `HF` (Atlantic 969, Pacific 856 in 2004-05 on, all months) | one per archive ID | first `HF` fix | 2004-05 to 2025-26 |
| 2 PIPE-A (replication) | `research/era5/hf_history/results/all_tracks.csv.gz`, `gust800_kt >= 71.7` | one per track | `peak_time` (time of the track's maximum gust index) | 2004-05 to 2025-26 (R1); 1979-80 to 2003-04 (R2, fitted separately, within-era only, decision 1) |
| 2b PIPE-A deep (pressure based) | same file, `minp <= 1000` and `n_fix >= 8` | one per track | `peak_time` | 1979-80 to 2025-26, fitted with an era term (R3) |

- Window: 1 October to 28 April (210 days, day 0 to 209), as in the additive test and PR 14. Events outside it are dropped from the main analysis and their number is reported. Season = start year.
- Position of an event (space tests): mean position (unit-vector mean) of its HF-category fixes (archive) or its track points at or above 71.7 kt (pipeline A, from `era5_hf_catalog_tracks.csv` for events; for deep cyclones, `peak_lat`, `peak_lon`).
- Period 1979-2003 gust-based counts are used only inside that era (R2); no level or trend is carried across 2001 (decision 1).
- Tropical-cyclone-origin events are in (pipeline A does not mask them); not separated here.

## Measures and nulls (fixed now)

Daily event counts per basin and season form the base table (210 days x seasons).

**M1. Dispersion of block counts.** Blocks: W7 (30 non-overlapping 7-day blocks per season from 1 Oct) and W30 (7 blocks of 30 days).
Statistic: Pearson dispersion phi = sum over blocks (y - mu)^2 / mu, divided by (blocks - parameters); the **excess dispersion** E = phi - E_null[phi]
(about variance/mean - 1). Only W7 and W30 enter the decision. A curve of E at block lengths 1, 2, 3, 7, 14, 30 days is drawn for description and carries no test.
- **N1 (primary null).** Poisson regression of the *daily* counts on calendar-month dummies and a centred linear season trend (the design of PR 14 and the additive test, without any index); mu = fitted daily means summed over each block. N1 treats interannual rate variation as part of the dispersion.
- **N2 (within-season null).** As N1 with a dummy for every season instead of the trend. E under N2 is clustering that remains after each season's rate is allowed to be its own.
- Reference distribution: 2,000 parametric draws from the fitted Poisson, each refitted and rescored; one-sided p for clustering (E larger than the null). Interval: 2,000 season-block bootstrap resamples (percentile 95%).
- Sensitivity (reported, no decision weight): two-harmonic seasonal cycle instead of month dummies; blocks shifted by 3 days; the Oct-Mar window; archive with ID-quality flags removed (`idOk` false) and the second member of each pair in `data/hf_lows/collision_pairs.csv` dropped.

**M2. Inter-arrival times.** Per basin and season, order the events and take gaps between successive events inside the window.
Primary statistic: the number of gaps of 48 h or less, compared with its expectation under N1 (events drawn from the fitted daily intensity, uniform within the day, 2,000 draws; one-sided p).
Secondary (reported): coefficient of variation of the gaps after time rescaling by the N1 intensity (1 under Poisson), a KS distance of the rescaled gaps from Exp(1) with the same reference distribution, and the share of events followed by another within 72 h.

**M3. Space-time interaction (Knox).** Pairs of events in the same basin and season. K = number of pairs with time separation at most 3 days and great-circle separation at most 1,000 km. Null: event times permuted among the events of the same season-month block (positions kept, 4,000 permutations), which holds the rate by season and month fixed; one-sided p; report K / E_null[K].
The grid of (2, 3, 5 days) x (500, 1000, 1500, 2000 km) is reported as a map of K / E[K] and has no decision weight; (3 d, 1000 km) is the single primary.

**M4. Local clustering.** Grid of 5 degrees latitude x 10 degrees longitude cell centres inside the basin domain. An event passes a cell if any of its HF fixes (archive; pipeline A track points at or above 71.7 kt) lies within 700 km of the centre; its local time is the first such fix. Weekly (W7) passage counts per cell; cells with fewer than 40 passages in the 22 seasons are not tested. Per-cell statistic E (as M1) against N1 refitted per cell (500 parametric draws per cell, fitted mean reused, no refit: stated approximation). Field significance: Benjamini-Hochberg FDR across the tested cells at q = 0.05, and the number of cells with p <= 0.05 compared with its distribution under a joint position-shuffle null (positions permuted among events of the same season-month, 1,000 permutations; the per-cell p-values for each permutation are ranked among the others). The position-shuffle null also gives a second map: local excess beyond the basin-wide timing of events.

**M5. Conditioning.** Daily models for mu, all with month dummies and the season trend:
- M0: none (= N1).
- M-idx (primary): basin index, NAO for the Atlantic and PNA for the Pacific, the CPC daily value averaged over days -10 to -4 before each day (the PR 14 lag), standardised over the analysis days.
- M-idx2 (secondary): both NAO and PNA, plus monthly ONI (lagged one month) in both basins.
- M-split (secondary): the seasonal mean of the index and the within-season anomaly as separate terms, to see whether any explained dispersion is interannual or sub-seasonal.
- M-season: season dummies (upper bound for what any interannual state can explain).
For each, E at W7 and W30 with the same reference distribution as M1. **Key number: share explained = 1 - E(M-idx) / E(N1)**, with a season-block bootstrap interval (both models refitted in each resample). It is interpreted only if E(N1) has a bootstrap interval above 0; otherwise reported as "no dispersion to explain".
A pattern index from the "Hemispheric patterns behind HF-active periods" thread would be added as a secondary conditioning term only after that thread has committed it. No jet-state index is available now; it is not substituted.

## Hypotheses, stated before results

- **H1 (time, basin scale).** Weekly excess dispersion is small: E(W7) under N1 at most 0.15 in the Atlantic and 0.10 in the Pacific. Monthly excess is larger (W30 greater than W7) but mostly interannual: E(W30) under N2 is under half E(W30) under N1.
- **H2 (inter-arrival).** The number of gaps of 48 h or less exceeds the N1 expectation by no more than 15%.
- **H3 (space).** Atlantic lows that arrive within 3 days of each other sit closer than chance (K / E[K] above 1.2); the Pacific shows less. This is the most uncertain prediction; a ratio near 1 in both basins is a full result.
- **H4 (conditioning).** NAO and PNA account for between 10% and 50% of E(N1) at W30 and for less at W7; they do not account for most of the excess.
Each can fail. A well-powered null is a full result.

## Decision rules (fixed now)

- **Clustered in time at scale W** in a basin: M1 N1 one-sided q < 0.05 (Benjamini-Hochberg over the tier-1 family below) and the bootstrap interval of E excludes 0. **Clustering within season** additionally needs the N2 interval of E to exclude 0. If N1 is positive and N2 is not, the dispersion is interannual rate variation, not bursts.
- **Bursts at the storm scale:** M2 primary q < 0.05.
- **Spatial:** M3 primary q < 0.05 and K / E[K] above 1.
- **Local map:** a cell is called out only if its BH q < 0.05 and the field count of p <= 0.05 cells exceeds the 95th percentile of the shuffle distribution.
- **Conditioning:** "mostly explained" if the whole interval of the share is above 0.5; "minor" if below 0.5; otherwise "not resolved".
- **Multiplicity.** Tier-1 family (BH, q = 0.05): M1 (W7 and W30, 2 basins, N1) = 4, M2 primary x 2 = 2, M3 primary x 2 = 2, the index coefficient of M-idx x 2 = 2: 10 tests. Tier 2 (pipeline A R1) gets the same 10 as its own family and counts as replication, which is weaker than independent: the atmosphere is the same, only event detection differs. R2 and R3 are reported as further replications, each its own family. M4 is its own family per basin. Sensitivities are not tested.
- Effective sample size is seasons (22 in tier 1), not events or weeks; intervals are season-block bootstraps and p-values are simulations that keep seasons whole.

## Power (computed in the analysis run, from fitted seasonal-cycle rates only)

A Neyman-Scott burst process is simulated on the fitted N1 intensity: a fraction e of events are secondary lows arriving after a parent with an exponential delay of mean 1.5 days and a position displaced by a normal 500 km. For e in {0.05, 0.10, 0.20, 0.30} the power of the W7 M1 test, the M2 test and the M3 test is reported at the real sample size, with the smallest e detected at 80% power. A null is read against that number.

## Seeds, code, checks

`cluster.py` (analysis), `power.py` (power), outputs in `results/`; seed 20261008. A fresh agent recomputes every quoted number from the committed files with its own code before merge; numbers it did not check are listed in the README.
