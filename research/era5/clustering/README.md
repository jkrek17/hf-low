# Do hurricane-force lows cluster in time, and along the same track?

Climatology agenda item 8 / ledger Future-task 9. Plan: `PREREGISTRATION.md` (commit `efc74a5`, before any statistic was computed). Code: `cluster.py`, `power.py`, `sensitivity.py`, `posthoc.py`, `figs.py`, `summarize.py`. Outputs: `results/`. Seeds fixed (`SEED = 20261008`). All **archive** numbers are OPC warnings; all **pipeline A** numbers are an **ERA5 proxy** (800 km ocean gust index, HF-equivalent at 71.7 kt), not observations. Pipeline B is not used.

## Answer in plain words

1. **In time: no.** HF lows do not arrive in bursts. At weekly scale the counts are, if anything, slightly *less* variable than a Poisson process with the seasonal cycle (excess dispersion E = -0.06 Atlantic, -0.03 Pacific, archive, 2004-05 to 2025-26), and neither basin has more gaps of 48 h or less than chance (0.91 x expected in both). At monthly scale there is a positive point estimate (E = +0.11 Atlantic, +0.08 Pacific) that is not distinguishable from zero and mostly vanishes once each season is given its own rate. The pipeline A proxy agrees (2004-05 on, and 1979-80 to 2003-04 within the era). The test would have caught a burst process in which about 10% of HF lows are followed by a second one within roughly 1.5 days (80% power at the weekly scale); the archive's upper 95% bounds on the weekly excess are E of about 0 (Atlantic) and 0.06 (Pacific), which corresponds to roughly 4% secondary lows or fewer. The proxy is looser in the Pacific (upper bound 0.19).
2. **In space: yes, modestly, and mainly in the Atlantic.** Two HF lows that occur within 3 days of each other are about 1.2 times as likely to lie within 1000 km of each other as chance allows (archive: Atlantic 1.20, Pacific 1.22; permutation p = 0.0002 and 0.0025, Benjamini-Hochberg q = 0.0025 and 0.010). The Atlantic value replicates in the proxy twice over (1.20 in 2004-05 on, 1.20 in 1979-80 to 2003-04). The Pacific one does not replicate in the proxy (1.13, q = 0.09; 1.07, p = 0.21). The excess sits at 2-4 days separation and is gone by about a week (post hoc, figure `results/knox_by_lag.png`). This is not "bursts in time": the arrival rate is flat, but successive lows tend to pass through the same region, so the pattern is persistence of where the storm track lies plus lows following one another along it. This design cannot tell those two apart.
3. **Where.** Beyond the basin-wide timing of events, the places that see HF lows in runs are in the central and eastern North Atlantic (about 45-60N, 10-35W), in both the archive and the proxy (field significant: archive 6 cells against 1.9 expected, p = 0.045; proxy 10 against 2.0, p = 0.001). No Pacific field is significant. By the pre-registered rule (a cell needs a Benjamini-Hochberg q below 0.05 against the per-cell Poisson null) **no cell is called out**: the per-cell Poisson test is weak at these counts (q at least 0.43). The map is the second, position-shuffle one (`results/local_cells.png`).
4. **How much do NAO and PNA explain?** The pre-registered quantity (share of overdispersion explained) **is undefined**: there is no overdispersion to explain at the weekly scale and the monthly one has an interval spanning zero. NAO and PNA do shift how many HF lows occur (archive per SD, lagged days -10 to -4: Atlantic NAO RR 1.064, permutation p = 0.12; Pacific PNA RR 1.119, p = 0.003; proxy 1.124 and 1.105), but that is a change in the mean rate, not in how bunched the arrivals are. For the spatial excess (post hoc): holding the lagged index tercile fixed removes about 17% of the Atlantic excess (1.20 to 1.17 archive, 1.20 to 1.16 proxy). Holding the same-time index fixed (circular, an upper bound) removes about 32% (archive) to 54% (proxy). The conditioned ratios stay above 1 with intervals that exclude 1 in the Atlantic (archive 1.02-1.24 same-time). Most of it remains. Pacific: 17% lagged and 38% same-time in the archive, intervals include 1 in the proxy.

Caution for item 1: deep cyclones of all kinds (pipeline A, central pressure at most 1000 hPa and at least 8 fixes, 47 seasons) are strongly *more regular* than Poisson (E = -0.37 at weekly scale in both basins). That is too large to be weather alone; the tracker's linking and 24 h minimum probably impose a minimum spacing. It does not affect the HF subset's near-zero result but means a Poisson null is a poor description of tracked lows in general.

## Design, as committed

Daily event counts per basin, 1 October to 28 April, seasons 2004-05 to 2025-26 (22). Archive: class `low` with at least one `HF` fix, time of the first HF fix (Atlantic 902 events in the window, Pacific 810; 66 and 46 outside it). Pipeline A (tier 2): `gust800_kt >= 71.7`, time of the peak gust index (Atlantic 1,013, Pacific 830; 89 and 62 outside). R2 1979-80 to 2003-04 fitted on its own (decision 1: within-era only). R3 deep cyclones, pressure based, 1979-80 to 2025-26 with an era term. Measures M1 to M5 and the nulls are in `PREREGISTRATION.md`; in short: M1 Pearson dispersion of 7-day and 30-day counts against a Poisson fit with month dummies and a season trend (N1) and with season dummies (N2), reference distribution from 2,000 refitted simulations, intervals from 2,000 season-block bootstrap resamples; M2 gaps of 48 h or less against N1; M3 Knox pair counts (3 days, 1000 km) against times permuted within season-month; M4 700 km cells on a 5 x 10 degree grid; M5 conditioning on the lagged NAO or PNA. Effective n is 22 seasons.

## Results

Archive (tier 1; `results/summary.txt` has every tier). Excess dispersion E is Pearson dispersion minus its null mean (about variance/mean - 1).

| | Atlantic | Pacific |
|---|---|---|
| Events in window | 902 | 810 |
| E, 7-day blocks, N1 (95% bootstrap) | -0.056 (-0.142, -0.001) | -0.028 (-0.148, 0.060) |
| E, 30-day blocks, N1 | 0.113 (-0.218, 0.341) | 0.077 (-0.155, 0.216) |
| E, 30-day blocks, N2 (season dummies) | 0.009 (-0.267, 0.215) | 0.029 (-0.204, 0.182) |
| Gaps of 48 h or less, observed vs N1 expectation | 311 vs 340.2 | 255 vs 280.6 |
| Knox pairs (3 d, 1000 km), observed vs permutation expectation | 251 vs 209.0 (1.20) | 135 vs 111.0 (1.22) |
| Knox permutation p; BH q (tier-1 family of 10) | 0.0002; 0.0025 | 0.0025; 0.010 |
| Index per SD (lagged), RR; permutation p | NAO 1.064; 0.12 | PNA 1.119; 0.003 |
| Field test, cells with shuffle p <= 0.05 vs null | 6 vs 1.9 (p = 0.045) | 3 vs 2.2 (p = 0.36) |

Tier-1 family (10 tests, BH q = 0.05): significant are the Knox tests (both basins) and the Pacific PNA coefficient. No temporal-clustering test (M1 7 d and 30 d, M2) is anywhere near significant (smallest q 0.32).

Replication (pipeline A proxy), same family: Atlantic Knox 1.20 (q = 0.0025) in 2004-05 on and 1.20 (q = 0.0025) in 1979-80 to 2003-04; Pacific Knox 1.13 (q = 0.09) and 1.07 (p = 0.21). 7-day E: Atlantic -0.023 and -0.054; Pacific 0.063 and -0.118 (the last significantly negative). Index coefficients 1.124 (Atl NAO, 2004-05 on) and 1.105 (Pac PNA). Atlantic field test 10 cells vs 2.0 (p = 0.001); Pacific 2 vs 2.05.

Power (`results/arch_power.json`, `r1_power.json`; 150 replicates per cell; burst process with a mean delay of 1.5 days and 500 km displacement): the 7-day test reaches 80% power when about 10% of events are secondary members (0.79 Atlantic, 0.84 Pacific); the 30-day test and the gap test need about 20%; the Knox test reaches about 60% at 5% and 95% at 10%. Size at e = 0 is 3-8% (nominal 5%, 150 replicates).

Sensitivity (pre-registered variants, archive; no decision weight): Knox ratio 1.18-1.22 in all variants including dropping the 41 Atlantic and 5 Pacific events with ID-quality flags or the second member of a collision pair (clean: 1.18, 1.20). The 7-day E stays between -0.095 and -0.005. The 30-day E ranges 0.03 to 0.18 across variants, raw p 0.08 to 0.40, never FDR-significant (`results/arch_sensitivity.json`).

### Post hoc (written after seeing the results above; exploratory)

`posthoc.py`, `results/{arch,r1}_posthoc.json`, `results/knox_by_lag.png`.
- Knox ratio by time separation (pairs within 1000 km, archive Atlantic / Pacific): up to 24 h 1.11 / 0.94; 24-48 h 1.18 / 1.21; 48-72 h 1.30 / 1.45; 72-96 h 1.12 / 1.25; 4-7 d 1.01 / 1.00; 7-14 d 0.94 / 1.01. The proxy has the same shape in the Atlantic (1.23, 1.08, 1.27, 1.21, 1.07, 0.92). The 2-4 day peak is what one expects from successive lows along a track; an excess at under 24 h (duplicate records of one storm) is not what drives it.
- Per season: the pooled 72 h ratio is above 1 in 16 of 22 Atlantic and 15 of 22 Pacific seasons (archive); 19 of 22 and 14 of 22 in the proxy. Season-block bootstrap intervals for the pooled ratio: archive Atlantic 1.08-1.31, Pacific 1.09-1.36; proxy Atlantic 1.07-1.31, Pacific 0.97-1.31.
- Lower tail (regularity): archive Atlantic p = 0.17 (7-day), 0.10 (gaps of 48 h or less); Pacific 0.31, 0.11. Not established for the HF subset.
- Conditioning of the Knox excess on NAO/PNA terciles: shares given above.
- None of these decides a pre-registered question.

## Hypotheses against results

- **H1** (weekly E at most 0.15 Atlantic and 0.10 Pacific; E at 30 days larger but mostly interannual): weekly part met (-0.06, -0.03). Monthly E is larger than weekly, and E under N2 is below half the N1 value (0.009 vs 0.113; 0.029 vs 0.077), but the monthly E is itself not distinguishable from zero, so "mostly interannual" is consistent but not resolved.
- **H2** (gaps of 48 h or less at most 15% above expectation): met (0.91 x).
- **H3** (Atlantic Knox ratio above 1.2, Pacific less): Atlantic 1.202 (1.08-1.31), nominally at the boundary and not distinguishable from 1.2; Pacific 1.22 in the archive, so "less" **failed** there (it is less in the proxy, 1.13).
- **H4** (NAO/PNA explain 10-50% of the 30-day excess): not assessable, there is no resolvable excess.

## Departures from the plan, and what was not done

- Post hoc analyses above were added after the pre-registered results. The pre-registered statistics are unchanged.
- The index coefficient's p-value method was not stated in the plan; the season-block permutation p (2,000 permutations, as PR 14) is used, with the cluster-robust SE reported alongside.
- R3 was run with lighter settings (500 simulations, 500 permutations, 500 bootstrap resamples) because of its 10,000 events per basin; M4 was not run for R2 or R3.
- Calendar months are taken on a fixed non-leap calendar (day index 0 to 209), so in leap years the March and April month boundaries shift by one day.
- Bootstrap intervals for E use the null mean from the point-estimate simulations; the cell tests (M4) reuse the fitted mean without refitting, as stated in the plan.
- M4's pre-registered call-out rule produced no cell. The shuffle-null map is reported as a second map, as planned, not as a call-out.
- Not done: the plan's "by subregion" dispersion beyond the M4 grid; a jet-state conditioning term (no index available; the hemispheric-patterns thread's index was not used, being uncommitted); the within-season clustering of tropical-cyclone-origin lows separately.
- Interpretation limits: a Knox excess says nothing about mechanism; storm-track persistence, repeated lows along one track and a shared remote driver are indistinguishable here. The deep-cyclone tier suggests the tracker enforces some spacing; the HF subset does not show it.

## Verification

See `VERIFICATION.md` (a fresh agent recomputing the quoted numbers from the committed inputs).

## Reproduce

```
cd research/era5/clustering
python3 cluster.py arch results; python3 cluster.py r1 results; python3 cluster.py r2 results
python3 cluster.py r3 results --nperm 500 --nboot 500 --nsim 500
python3 power.py arch results; python3 power.py r1 results; python3 sensitivity.py results
python3 posthoc.py arch results; python3 posthoc.py r1 results
python3 figs.py results results; python3 summarize.py results > results/summary.txt
```
Inputs: `docs/data/hf-lows.json`, `data/hf_lows/collision_pairs.csv`, `research/era5/hf_history/results/{all_tracks.csv.gz,era5_hf_catalog_tracks.csv}`, CPC NAO and PNA and ONI from `/mnt/project-files/teleconnection-test/cpc_indices/`. About 10 minutes per tier on 4 cores.
