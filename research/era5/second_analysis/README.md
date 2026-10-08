# RA-27: is the unexplained part of the HF conversion step noise at the threshold?

Pipeline A proxy (ERA5), not the archive. Plan: `PREREGISTRATION.md` (committed first, `c33dd9f`; floor change committed before any field read).
Looks entry: `../looks/second_analysis.log`.

## Answer in plain words

**Mostly noise, for sustained 10 m wind: two analyses agree on which storms cross only at HSS 0.61 (95% CI 0.55 to 0.66).** The gust
rule itself could not be run, because the IFS HRES t0 analyses in WeatherBench2 carry no gust variable. So this bounds the noise of a
sustained-wind crossing, not of the gust index. Treat it as a lower bound on representation noise (the two systems share observations),
and as an upper-bound-on-skill argument: no field-based forecast of ERA5's crossing is expected to reach much past 0.6 to 0.7 when a
second analysis of the same storms reproduces it only at 0.61. Forecast skill so far is 0.59 (PR 12). Big caveat, below: the sample is storms
whose ERA5 gust is already at least 65 kt, which compresses the range.

## What was run

944 pipeline A tracks (Atlantic 499, Pacific 445), seasons 2016-17 to 2021-22, ERA5 gust index (00/12 UTC) at least 65 kt, 429 of them at
or above 71.7 kt. Per fix: maximum 10 m wind speed over owned ocean cells within 800 km, in ERA5 (ARCO 0.25 degree) and the IFS HRES t0
analysis (WeatherBench2, 0.25 degree). Centres and ownership are ERA5's for both. Cuts count-matched to 429 crossings in each system
(ERA5 48.6 kt, HRES 53.1 kt). Pull 42.6 GB (u, v, MSLP, plus gust for 150 validation times). Code: `select.py`, `size_pull.py`, `extract.py`, `analyse.py`;
outputs `results/analysis_65.txt`, `fix_index_65.csv.gz`, `track_index_65.csv.gz`.

## Numbers (`results/analysis_65.txt`)

| Test | Result |
|---|---|
| T1 pooled HSS | 0.611, CI [0.549, 0.660], p(H0 <= 0.6) 0.42, q 0.74 |
| T2 Atlantic | 0.599 [0.518, 0.671], q 0.74 |
| T3 Pacific | 0.597 [0.522, 0.691], q 0.74 |
| T4 Atlantic minus Pacific | +0.002 [-0.123, +0.102], q 0.83 |

Pre-registered rule: upper bound 0.66 < 0.70, so **mostly noise**. 0 of 4 tests reject "HSS <= 0.6"; the prediction (about 0.6) is not contradicted.
Power: bootstrap SE 0.028, so 80% power against a true HSS of 0.67; a well-powered statement that the agreement is not above about 0.67.
Season-block bootstrap (6 seasons) [0.566, 0.646].

Secondary (descriptive):
- S2: P(HRES crosses | ERA5 index / cut): <0.85: 0.00 (n 8); 0.85-0.95: 0.10 (276); 0.95-1.05: 0.39 (399); 1.05-1.15: 0.92 (171); >=1.15: 1.00 (90). The transition zone is about +-10% of the cut.
- S6: HRES index = 1.17 x ERA5 index - 3.4 kt, correlation 0.79, residual SD 4.4 kt against an ERA5 spread of 3.3 kt among tracks within 15% of the cut. HRES runs 5.0 kt higher on average (finer model); 99th percentiles ratio 1.28.
- S3a ERA5 sustained vs ERA5 gust crossing: HSS 0.72 [0.67, 0.76]. S3b HRES sustained vs ERA5 gust: 0.58 [0.52, 0.63].
- S4 (gust index 65-80 kt, the band around the cut): HSS 0.51. S8 (gust index >= 70): 0.53. S5 fix level: 0.69.
- **S7, 99th percentile of the owned cells instead of the maximum: HSS 0.80 [0.75, 0.83].** The maximum is a single-cell statistic and noisier than the 99th percentile.

## How far to trust it (read before quoting)

1. **Range restriction, a design weakness I did not foresee.** The sample keeps only tracks with gust index at least 65 kt, and 429 of 944 cross, so the count-matched wind cut sits near the middle of the sample (about a median split), not at a rare-event tail. Agreement among storms that are already similar is lower than agreement across all cyclones. The 0.61 is therefore a statement about discriminating among strong storms near the cut, which is the question, but it is not comparable with the label-versus-archive HSS of 0.66 (PR 68) that covers a wider range.
2. **The index choice matters a lot**: with the 99th percentile instead of the maximum, agreement is 0.80. The "noise" is partly the single-cell maximum. A gust index is also a maximum, so its noise is not bounded by the 0.80.
3. Sustained wind, not gust: shown above. Systems share observations; HRES regridded from 9 km.
4. Centres are ERA5's; timing noise is the same in both systems.
5. The floor was 65 kt, not 55 (pre-run deviation in the plan), so tracks of 55 to 65 kt are not covered.

## Deviations and checks

- Pre-run: floor 65 kt (sizes omitted MSLP); logged in the plan before any field was read.
- Post hoc, labelled: the plan's validation ("ERA5 gust re-derived with this ownership must return `g800` within 0.5 kt for >= 99%") read 98.8% raw (407 of 412 usable). The 5 misses are all fixes whose catalog `g800` is 0.0 because pipeline A only assigns cells with gust at least 17 m/s (33 kt) to a low, while the re-derived value was 18 to 32 kt. Excluding those 5 (re-derived value below 33 kt and catalog 0), all 407 match to within 0.05 kt. The definition is clarified, not the threshold loosened; the raw figure is shown above.
- `analyse.py` was fixed for a column collision before the first successful run (no result changed by it).
- Independent recomputation: see the section below (filled in by the fresh verifier).
