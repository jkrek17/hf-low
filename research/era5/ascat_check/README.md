# ASCAT check of the ERA5 gust proxy (status: RUN; numbers recomputed by a second agent, see "Verification")

Question (Jason, 2026-10-13, item 2): does pipeline A's gust index (**ERA5 proxy, pipeline A**; transitioning tropical cyclones are in)
agree with scatterometer (ASCAT) winds? **ERA5 assimilated scatterometer winds** (ASCAT on Metop-A and B; ERA5 documentation, Hersbach et al.
2020; the exact dates and whether Metop-C enters were not confirmed here) inside the same 12 h windows used here, so the ASCAT winds are **not independent
of ERA5**: this is a consistency check, probably flattering, not independent validation. Plan: `PREREGISTRATION.md`, committed first
(`5440a7b`); deviations 1-7 logged in it (1-6 before any wind comparison, 7 post hoc). It tests the proxy against ASCAT winds, not whether
OPC's hurricane-force call was right (that call is an analyst's judgement).

## Answer (the registered verdict rule gives "yes")

**Yes, as a consistency check:** ERA5's maximum gust on the cells ASCAT saw ranks storms by ASCAT's maximum wind very well (rank correlation
0.92 to 0.95, AUC 0.97), and storms with a gust index at or above 71.7 kt have a clearly higher ASCAT footprint maximum than storms of the same
depth with a lower index. At matched quantiles ERA5 gusts read about 15 to 25 % above the ASCAT footprint maximum (54 to 71.7 kt range). That
ERA5 itself assimilated ASCAT probably inflates these agreements; agreement may not carry over to years when no scatterometer was assimilated.

Sample: pipeline A storms from 2006-07 to 2025-26 (20 seasons; ASCAT-A starts 2007-03, so 2004-05 and 2005-06 have no pass) whose peak-index fix is at
00 or 12 UTC with a stored index above 0, and an ASCAT pass (+-1.5 h, cell by cell) that sees at least 40 % of the ocean cells within 800 km.
Candidates: 13,851 from 2004-05, **12,568 from 2006-07**; 12,327 storm-pass footprints, 5,791 with >= 40 % coverage. The analysis keeps the
best-covered pass per storm: **4,252 storms** (34 % of the 12,568; Atlantic 2,352, Pacific 1,900), **430 with index >= 71.7 kt** (265 and 165).
Metop-A/B/C chosen for 1,886 / 1,632 / 734.

| Registered test (BH-FDR over these 5) | Result | Threshold | q |
|---|---|---|---|
| P1 rank correlation of ERA5 vs ASCAT footprint maxima, Atlantic | **0.946** [0.940, 0.952] | >= 0.6, lower >= 0.5 | < 0.001 |
| P1 same, Pacific | **0.924** [0.913, 0.932] | same | < 0.001 |
| P2 AUC, ERA5 footprint max separating ASCAT max >= 45 kt (pooled) | **0.971** [0.967, 0.975] | >= 0.85, lower >= 0.80 | < 0.001 |
| P3 ASCAT footprint max (25 km mean wind over the ~52 % of the disc ASCAT saw), index >= 71.7 minus depth-matched lower-index storms | **+12.0 kt** [11.1, 13.1] (n 423); [10.4, 13.8] if the matching is redone inside each bootstrap draw | > 0 and q < 0.05 | 0.0012 |
| P4 ASCAT / ERA5 ratio differs among Metop-A, B, C | F 5.33, p 0.037 | none (stability check) | 0.037 |

P1-P3 meet the "yes" rule. P4 is also significant (q 0.037), which means the sensors differ; that is the adverse direction for stability, not a pass
(see below). Intervals are season-block bootstrap (resampling the 20 seasons); bootstrap p of 0.000 means below 1/1,000. P1 uses the ERA5 maximum on
the cells ASCAT saw, not the stored full-disc index (the stored index against the ASCAT maximum gives a lower correlation, partial coverage
and neighbouring lows; not reported as a result here). The rule asks for P2 in both basins; P2 was computed pooled. Post hoc per basin:
Atlantic 0.971 [0.965, 0.976], Pacific 0.972 [0.966, 0.977].
Power: the P3 interval gives a standard error of about 0.5 kt (0.9 kt with the matching redone in each draw), so the smallest difference
detectable at 80 % power is about 1.4 to 2.4 kt (2.8 x SE), far below the +5 kt in the "no" rule.

## What the numbers mean in wind terms (descriptive, no verdict)

- Where the index is at or above 71.7 kt the ASCAT maximum on the same cells has median 62.8 kt; 91 % reach 45 kt and 84 % reach 50 kt. In the 54-64 kt
  class the median ASCAT max is 44.3 kt, in the < 54 class 31.9 kt.
- Matched quantiles (a distribution mapping, not storm by storm; both basins): an ERA5 footprint maximum of 71.7 kt sits at the same quantile as an ASCAT
  maximum of **62.2 kt** (ratio 0.87); 64 kt with 54.8 kt (0.86); 54 kt with 43.5 kt (0.81). Atlantic 63.0 / Pacific 61.4 kt at 71.7. ASCAT loses
  sensitivity and tends to read low above about 50 kt (values to 97 kt are reported); winds above 60 kt, including 62.2 and 62.8 kt, are in the range the
  plan said would be reported but never used as truth. The 71.7 kt level cannot be tested directly.
- Pointwise (8.2 million cell pairs on the analysis set; **not independently recomputed**, no intervals): Pearson 0.94, mean wind/gust ratio 0.78 (slope 0.77),
  bias -6.5 kt south of 40N to -5.2 kt north of 60N, and -5.2 kt within 200 km of the fix, -6.1 at 200-400 km, -5.8 at 600-800 km (not monotonic).

## Stability and sensitivity

- **P4** (registered): the ASCAT/ERA5 ratio for ERA5 maxima >= 54 kt differs among sensors within seasons: A 0.831, B 0.838, C 0.824 (1.4 point spread,
  1.7 % relative; q 0.037). It pools all 5,791 passes from the 4,252 storms, so a storm can count more than once, and it tests a sensor difference, not
  an era drift. The plan said a drift here "would reopen Decision 1's question for 2004 on" and set no size threshold: **whether this does is for Jason**.
  Descriptive, no interval or test: ratios by season run 0.81 to 0.86, slope +0.006 log units per decade.
- Post hoc, from a run that enters no verdict and was **not independently recomputed** (10.3 % of cells dropped instead of 0.9 %; 3,450 storms, P3 n 314):
  with radiometer-flagged rain also dropped, P1 0.949 / 0.931, P2 0.976, P3 +10.9 kt [9.8, 12.3], and P4 p 0.32.
- Post hoc, not registered: the rank correlation is not an artefact of footprint size (0.94-0.96 Atlantic, 0.90-0.94 Pacific within thirds of cell
  count; partial rho 0.945 and 0.921).
- Post hoc, Atlantic split at 60N (Decision 6 registered it; first reporting omitted it; `results/post_hoc.txt`): rho 0.939 [0.929, 0.946] at or south of 60N
  (n 1,566), 0.964 [0.955, 0.970] north of 60N (n 786); P3 with the northern Atlantic dropped +12.5 kt [11.1, 13.8]. No change to the verdict.

## Limits (read before quoting)

- ASCAT here is the RSS v02.1 0.25 degree daily product. It can confirm ordering and a storm-force-equivalent (45 kt) level but not the 71.7 kt level;
  QuikSCAT (2004-2009) was not run.
- The 4,252 analysed storms are 34 % of the 12,568 candidates from 2006-07 on; coverage depends on orbit, land, ice and storm size, and the accepted
  sample is enriched in strong storms (index >= 71.7 is 10.1 % of analysed storms against 6.9 % of candidates). It is not the archive's mix.
- The index is a gust over the ocean within 800 km, the ASCAT value a 25 km area-mean equivalent-neutral wind: the ratio mixes gust factor, area
  averaging, ASCAT's high-wind response and ERA5's assimilation of ASCAT; partial coverage cancels because the comparison is on the same cells.
  Converting an ERA5 value to an ASCAT-scale wind needs this mapping; converting to the sustained wind OPC uses is not tested here.
- Seasons: 20 (2006-07 to 2025-26). This supports the gust index as a wind-class ordering from 2006-07 on and says nothing about pre-2001 levels.

## Verification

A fresh Sonnet agent, given the footprint table, the raw cells, the ERA5 cache and the definitions but not the analysis code or this README,
recomputed (`verify/VERIFICATION.md`, scripts in `verify/`): storm counts, both Spearman correlations, AUC, P3 (+12.003 kt, n 423), the class
table, the matched quantiles and the per-sensor ratios: all match; bootstrap intervals for both correlations and AUC match to 0.001. The
**P3 interval differs by method** ([11.1, 13.1] with the matching fixed once; [10.4, 13.8] redone per draw); both are shown. It also rebuilt 40
random footprints (10 with index >= 71.7) from the raw cells and ERA5 cache: 40 of 40 match. The statistics were recomputed from the analysis's own footprint
table; only the 40 were rebuilt from raw cells. **Not independently recomputed:** pointwise cell-pair table, rain-flag loss shares (0.9 % of cells), season-ratio
series and slope, tercile and partial correlations, the radiometer-rain run, the synthetic test, the post hoc items above (an Opus reviewer's own numbers for
per-basin AUC and the 60N split agree), the coverage fix's effect on counts beyond the 40 rebuilt rows, and the file counts and GB.
The commit subject of `a5b8d9d` ("ERA5 gust proxy tracks ASCAT winds") overstates this; do not quote it.

## Files

- Code: `blosc_nolib.py` (cache decoder via system libzstd; PyPI is blocked in this environment), `rss_adapter.py` (RSS bytemaps to cell rows),
  `collocate.py`, `analyse.py` (P1-P4), `descriptives.py`, `post_hoc.py`, `run_all.py`, `recover_index.py` (R0), `test_synthetic.py`, `verify/`.
  Results: `results/p1_p4.txt`, `results/descriptives.txt`, `results/post_hoc.txt`, `results/scat_footprints.csv.gz`,
  `results/sensitivity_scat_plus_radiometer.txt`.
- Shared folder `/mnt/project-files/ascat-check/` (README there): `scat_footprints.csv`, `cell_pairs_ascat_vs_era5.csv.gz` (15.3 million matched
  cells; 8.2 million of them in the analysis set), `cells/` (15.3 million ASCAT cell rows within 800 km and +-1.5 h of the fixes), `rss_inventory.csv`.
  Re-run without downloading: `python3 run_all.py <cells dir> <out prefix>` then `python3 analyse.py <out>_footprints.csv`.
- Download: 13,059 RSS daily files listed (20.0 GB), 11,974 opened (18.3 GB), nothing kept but the cell rows. The plan said to stop and ask Jason over
  20 GB; the listing was 20.0 GB, the opened files 18.3 GB, and no ask was made (the standing gate is 50 GB).
