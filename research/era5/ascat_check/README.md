# ASCAT check of the ERA5 gust proxy (status: RUN; verification noted below)

Question (Jason, 2026-10-13, item 2): does pipeline A's ERA5 gust index agree with independent scatterometer winds?
**ERA5 proxy, pipeline A** (transitioning tropical cyclones are in). Plan: `PREREGISTRATION.md`, committed first (`5440a7b`),
deviations 1-6 logged in it, all before any wind value was compared. This tests whether the proxy tracks observed winds. It does
not test whether OPC's hurricane-force call was right (that call is an analyst's judgement).

## Answer (verdict rule from the plan: "yes")

**Yes: the ERA5 gust index sorts lows by observed wind speed very well, and the hurricane-force class carries clearly stronger
observed winds.** It does not read ASCAT wind in knots: ERA5 gusts are about 15-25 % higher than the ASCAT footprint maximum.

Sample: pipeline A storms from 2006-07 to 2025-26 (20 seasons; ASCAT-A starts 2007-03, so 2004-05 to 2005-06 have no pass) whose
peak-index fix is at 00 or 12 UTC with a stored index above 0, and an ASCAT pass (+-1.5 h, cell by cell) that sees at least
40 % of the ocean cells within 800 km. 13,851 storms in the candidate set, 12,327 storm-pass footprints, 5,791 with >= 40 %
coverage; the analysis keeps the best-covered pass per storm: **4,252 storms** (Atlantic 2,352, Pacific 1,900), of which
**430 have an index at or above 71.7 kt** (265 and 165). Metop-A/B/C chosen for 1,886 / 1,632 / 734.

| Registered test (BH-FDR over these 5) | Result | Threshold | q |
|---|---|---|---|
| P1 rank correlation of footprint maxima, Atlantic | **0.946** [0.940, 0.952] | >= 0.6, lower >= 0.5 | < 0.001 |
| P1 same, Pacific | **0.924** [0.913, 0.932] | same | < 0.001 |
| P2 AUC, ERA5 footprint max separating ASCAT max >= 45 kt | **0.971** [0.967, 0.975] | >= 0.85, lower >= 0.80 | < 0.001 |
| P3 ASCAT max, index >= 71.7 minus depth-matched lower-index storms | **+12.0 kt** [11.1, 13.1] (n 423) | > 0 and q < 0.05 | 0.0012 |
| P4 ASCAT / ERA5 ratio differs among Metop-A, B, C | F 5.33, p 0.037 | none (stability check) | 0.037 |

Intervals are season-block bootstrap (resampling the 20 seasons). 5 of 5 tests pass q < 0.05 (the bootstrap p of 0.000 means
below 1/1,000). Power: the P3 interval gives a standard error of about 0.5 kt, so the smallest difference detectable at 80 %
power is about 1.4 kt (2.8 x SE), far below the +5 kt in the "no" rule; the "no" branch could not have been missed for lack of power.

## What the numbers mean in wind terms (descriptive, no verdict)

- Where the index is at or above 71.7 kt (hurricane force proxy) the ASCAT maximum on the same cells has median 62.8 kt;
  91 % reach 45 kt and 84 % reach 50 kt. In the 54-64 kt class the median ASCAT max is 44.3 kt, in the < 54 class 31.9 kt.
- Matched quantiles (both basins): an ERA5 footprint maximum of 71.7 kt corresponds to an ASCAT maximum of **62.2 kt** (ratio 0.87),
  of 64 kt to 54.8 kt (0.86), of 54 kt to 43.5 kt (0.81). Atlantic 63.0 / Pacific 61.4 kt at 71.7. ASCAT winds saturate and read low
  above about 50 kt and the footprint covers about half the disc, so these ASCAT numbers are lower bounds on the true maximum:
  winds above 60 kt were reported but never used as truth.
- Pointwise (8.2 million cell pairs on the analysis set): Pearson 0.94, wind = 0.77 x gust on average, bias -6.5 kt south of
  40N to -5.2 kt north of 60N, and -5.2 kt within 200 km of the fix to -5.8 kt at 600-800 km.

## Stability and sensitivity

- **P4 is narrowly significant but the size is small:** ratio ASCAT/ERA5 for ERA5 maxima >= 54 kt is 0.831 (A, 2006-2021),
  0.838 (B) and 0.824 (C), a spread of 1.4 %; by season it stays between 0.81 and 0.86 with a slope of +0.006 log units per decade
  (descriptive). The test turns not significant (p 0.32) when radiometer-flagged rain cells are also dropped
  (`results/sensitivity_scat_plus_radiometer.txt`), so it is a fragile detection of a 1 % difference, not evidence of a sensor
  or era drift. It does not reopen Decision 1 on this evidence.
- Radiometer-rain sensitivity (10.3 % of cells dropped instead of 0.9 %): P1 0.949 / 0.931, P2 0.976, P3 +10.9 kt [9.8, 12.3].
  Same verdict.
- Post hoc, not registered: the rank correlation is not an artefact of footprint size. Within thirds of the number of ASCAT cells it
  is 0.94-0.96 (Atlantic) and 0.90-0.94 (Pacific); partial rho given cell count 0.945 and 0.921.

## Limits (read before quoting)

- ASCAT here is the RSS v02.1 0.25 degree daily product (C-2015 model function). It saturates near 50 kt, so it can confirm
  ordering and a storm-force-equivalent level but cannot test the 71.7 kt level directly; QuikSCAT (2004-2009) was not run.
- Coverage rule keeps 31 % of the candidate storms (4,252 of 13,851); the pass geometry is set by the orbit, not by the storm, but
  the sample is not the full archive. Atlantic fixes north of 60N are included (lower confidence; not split out).
- The index is a gust over the ocean within 800 km, the ASCAT value a 25 km area-mean wind: the ratio 0.87 mixes gust factor,
  saturation and partial coverage. The proxy tracks observed winds; calling a given ERA5 value "X kt sustained" needs this
  mapping, not 1:1.
- Seasons: 20 (2006-07 to 2025-26). 2004-05 to 2005-06 have no scatterometer pass here. Decision 1 holds; this check supports the
  gust index as a wind-class ordering on 2006-07 on and says nothing about pre-2001 levels.
- Not independently recomputed until the verification section below says so.

## Files

- Code: `blosc_nolib.py` (cache decoder via system libzstd; PyPI is blocked in this environment), `rss_adapter.py` (RSS
  bytemaps to cell rows), `collocate.py`, `analyse.py` (P1-P4), `descriptives.py`, `run_all.py`, `recover_index.py` (R0),
  `test_synthetic.py`. Results: `results/p1_p4.txt`, `results/descriptives.txt`, `results/scat_footprints.csv.gz`.
- Shared folder `/mnt/project-files/ascat-check/` (README there): `scat_footprints.csv`, `cell_pairs_ascat_vs_era5.csv.gz` (15.3
  million matched cells; 8.2 million of them in the analysis set), `cells/` (15.3 million ASCAT cell rows within 800 km and +-1.5 h of the fixes),
  `rss_inventory.csv`. Re-run without downloading: `python3 run_all.py <cells dir> <out prefix>` then `python3 analyse.py <out>_footprints.csv`.
- Download: 11,974 of 13,059 daily files opened (18.3 GB of the 20.0 GB listed), nothing kept but the cell rows.
