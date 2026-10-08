# Does a farther-west El Nino mean more cyclone genesis over the Kuroshio? (ERA5 proxy)

Plan, committed before any genesis count or SST index was computed: [PREREGISTRATION.md](PREREGISTRATION.md) (commit `e317203`).
Code: `common.py`, `sst_indices.py`, `lows_msl.py`, `tracks.py`, `indices.py`, `validate.py`, `analysis.py`, `maps.py`.
Results: [results/](results/) (`summary.txt` is the readable report, `results.csv` has every test, `genesis_maps.png` the descriptive maps).
Small inputs behind them: [data/](data/) (SST box means, winter indices, tracks, fixes, pipeline A match table).

All numbers are an **ERA5 proxy**. Cyclones come from a **new MSLP-only tracker run that reuses pipeline A's detector and linker** (`hf_history/extract.py::lows_at`, `track.py::link`). Hurricane-force (HF) labels come from **pipeline A** (800 km ocean gust index at 71.7 kt). Pipeline B is not used.

## Answer

**No detectable effect, and the test can exclude anything much bigger than +8% per SD.** Over 47 winters (1979-80 to 2025-26), with Nino3.4 and a linear trend held, one SD more of the Modoki index (EMI) changes DJF genesis in the Kuroshio box by **0.995x (95% CI 0.917-1.079)**, permutation p = 0.89. The registered rule for a well-powered null (upper bound below +10% per SD) is met. The data do not support "more genesis when the warm anomaly sits farther west" at an effect of that size, and every one of the 15 registered tests has q (BH) of 0.96.

For hurricane-force lows the answer is **cannot tell**: 22 usable winters, and the effect detectable at 80% power is about +41% per SD for the count of K-genesis storms that later reach HF and +62% (odds) for the HF share. The point estimates (1.055, 1.047) are far inside the noise.

## Numbers (per +1 SD of the flavor index; 10,000 permutations; BH over 15 tests)

Box K = ocean first fixes in 25-35N 120-142E (East China Sea and south of Japan) or 30-42N 142-165E (Kuroshio Extension). Mean K genesis per DJF winter is 28.1 (SD 5.0) against 164.8 in the whole North Pacific. Mean of K-genesis storms that reach HF: 10.8 a winter 1979-2025, 11.5 over 2004-05 on (about 40% of K-genesis storms).

| id | what | n winters | effect | 95% CI | detectable at 80% | p (perm) | q (BH) |
|---|---|---|---|---|---|---|---|
| **P1** | DJF count in K ~ EMI + N34 + year | 47 | 0.995 | 0.917-1.079 | 1.12 | 0.895 | 0.96 |
| S1 | share of North Pacific genesis in K | 47 | 0.992 | 0.908-1.085 | 1.14 | 0.850 | 0.96 |
| S2 | K-genesis storms reaching HF, 2004-05 on | 22 | 1.055 | 0.829-1.343 | 1.41 | 0.656 | 0.96 |
| S3 | HF share of K-genesis storms (odds), 2004-05 on | 22 | 1.047 | 0.749-1.464 | 1.62 | 0.792 | 0.96 |
| S4 | P1 for Jan-Mar | 47 | 0.958 | 0.884-1.038 | 1.12 | 0.235 | 0.96 |
| S5 | P1 for Dec-Mar | 47 | 0.965 | 0.899-1.036 | 1.11 | 0.321 | 0.96 |
| S6 | P1 with N4-N3 instead of EMI | 47 | 0.998 | 0.939-1.062 | 1.09 | 0.962 | 0.96 |
| S7 | box shifted 5 N | 47 | 0.994 | 0.919-1.074 | 1.12 | 0.864 | 0.96 |
| S8 | box shifted 10 E | 47 | 0.997 | 0.905-1.098 | 1.15 | 0.949 | 0.96 |
| S9 | East China Sea and south of Japan only | 47 | 0.920 | 0.788-1.075 | 1.25 | 0.298 | 0.96 |
| S10 | Kuroshio Extension only | 47 | 1.030 | 0.936-1.133 | 1.15 | 0.519 | 0.96 |
| S11 | HF-reaching count 1979-2025 with an era term | 47 | 0.985 | 0.865-1.122 | 1.21 | 0.819 | 0.96 |
| S12 | CP (n=4) minus EP (n=9) El Nino winters, ratio of means | 13 | 0.978 | | | 0.877 | 0.96 |
| S13 | P1 without the year term | 47 | 0.998 | 0.920-1.081 | 1.12 | 0.951 | 0.96 |
| S14 | P1 without the Nino3.4 term | 47 | 0.983 | 0.930-1.038 | 1.08 | 0.513 | 0.96 |

S12: CP winters 1987, 2002, 2014, 2018; EP winters 1982, 1986, 1991, 1994, 1997, 2006, 2009, 2015, 2023 (Kug rule on ERA5 SST indices). Mean DJF K genesis per 88 days: CP 27.2, EP 27.8. With four CP winters this says almost nothing on its own.

**Mechanism (secondary, outside the FDR family).** Adding DJF PNA to P1 leaves the EMI effect at 0.999 (0.919-1.085); PNA itself has 1.015 per SD (0.951-1.084, p = 0.64) for K genesis. EMI and PNA correlate at 0.28, EMI and Nino3.4 at 0.73. There is nothing to mediate. A jet-level pull was not made.

## Power and what limits it

- The sample is 47 winters, one flavor value each. A winter's genesis count (about 28) is nearly Poisson (variance to mean 0.90, so the dispersion is held at 1). On that, the registered 80%-power effect for P1 is **+12% per SD**; the interval's upper end is 8%.
- EMI and Nino3.4 correlate at 0.73 (variance inflation 2.15), so "flavor at fixed strength" costs precision. S14 drops the control and gets the same answer (0.983).
- Only 13 winters have Nino3.4 >= +0.5 and just 4 of those are CP by the Kug rule, so any CP-vs-EP statement is weak. The continuous index is the better-powered test.
- The HF outcomes use 22 winters (2004-05 on), as Decision 1 requires, and are underpowered. S11 adds 1979-2000 with an era term (within-era information only) and does not change the picture.

## Validation of the tracker (before any outcome test, as registered)

All 12,355 pipeline A Pacific tracks whose peak fix falls in the pulled windows are recovered as a fix of a new track (match on peak time and position, rate 1.0000 in every decade). The new run therefore reproduces pipeline A's low detection and linking. ERA5 SST box indices correlate 0.99 with the official ONI (Nino3.4 vs ONI, DJF, 47 winters).

## What this does not show

- **A proxy.** Genesis here is the first fix of an ERA5 MSLP low below 1010 hPa that lives at least 24 h and starts over ocean. It is not the archive and not a forecaster's genesis point. The detector smooths on about 0.5 degrees and drops lows within 400 km of a deeper one, which can hide a weak low next to a strong one.
- **Pressure-tracked genesis is used for all 47 winters.** The HF label is gust-based and carries the pre-2001 drift (Decision 1); only S11 touches pre-2004 HF and only with an era term.
- **Box choice.** The box was fixed from the literature before looking. The shifted and split boxes (S7-S10) agree; S9 (East China Sea and south of Japan) has the widest interval and the lowest point estimate (0.92, p = 0.30).
- **Not tested:** genesis at other lags than concurrent DJF, ENSO phase interactions (the thread "El Nino and PNA together" owns ONI x PNA), a Kuroshio SST-front-relative definition of the box, and Atlantic. Nov is not in the pull, because it would have taken the MSLP pull from 49.5 GB to 61.9 GB.
- **ERA5 SST** switches product around 2007 (HadISST2, then OSTIA); monthly means here use 4 samples a month. The 0.99 correlation with the official ONI shows the amplitude is right; the flavor contrast (EMI) is not checked against another SST product.
- **The maps** (`results/genesis_maps.png`) are descriptive, with 4 CP and 9 EP winters and 5 degree bins; they show noise-level differences and no coherent pattern near the Kuroshio.

## Verification

A fresh Sonnet agent that had not seen the code or my numbers recomputed from the committed `data/` files with its own implementation, and its values agree with the ones above: winter indices (to 1e-15), first fixes of all 42,923 tracks against `fixes.csv.gz`, K genesis per winter (mean 28.106, SD 5.018; 23, 32, 25 in 1982, 1997, 2015), basin mean 164.77, HF-reaching mean 10.830 (1979-2025) and 11.455 (2004-05 on), P1 0.995 (0.917-1.079) with permutation p 0.888, S1 0.992 (0.908-1.085), S2 1.055 (0.829-1.343), S3 1.047 (0.749-1.464), the CP and EP winter lists and means (27.17, 27.78), and that all 12,355 pipeline A matches exist in `tracks.csv.gz`. Box edges are inclusive on all sides in both implementations; its HF label joined on track id, which is equivalent to the registered peak-fix match because every match is a track id.

**Not independently checked:** S4-S11, S13, S14, the detectable-effect figures, the BH q-values, the mechanism numbers (M1, M2), the maps, and the SST-to-ONI correlation (0.99).

## Post-hoc deviations

1. Permutation scheme: the registered Freedman-Lane scheme is used for the count outcomes (P1, S2, S4-S11, S13, S14). For the binomial outcomes (S1, S3) it would put pseudo-observations outside [0, 1], so the permutation is of the part of the predictor not explained by the controls. P1 under that scheme gives p = 0.86 against 0.89 for Freedman-Lane.
2. Dispersion is floored at 1 (counts are under-dispersed, variance to mean 0.90); this keeps the intervals from shrinking below Poisson.
3. The first run of `analysis.py` used 200 permutations as a smoke test and gave the same effects; the numbers above come from the 10,000-permutation run.

## Reproduce

    python3 research/era5/enso_kuroshio/sst_indices.py work/sst_boxes_raw.csv                 # about 3 GB
    python3 research/era5/enso_kuroshio/lows_msl.py 1979 2025 4                               # about 49.5 GB, resumable, writes work/lows_msl/
    python3 research/era5/enso_kuroshio/tracks.py work work/out
    python3 research/era5/enso_kuroshio/indices.py work/sst_boxes_raw.csv <cpc_indices> <repo root> work/winter_indices.csv
    python3 research/era5/enso_kuroshio/validate.py research/era5/hf_history/results/all_tracks.csv.gz work/out
    python3 research/era5/enso_kuroshio/analysis.py work/out/tracks.csv.gz work/out/pipelineA_match.csv work/winter_indices.csv results 10000

Run from `research/era5/enso_kuroshio/` (the scripts import `common.py` and `hf_history/track.py`). `lsm.npy` is the ERA5 land-sea mask chunk 876576. Both pulls were measured at 2.17 MB per MSLP step and 1.37 MB per SST step; no cache covered them.
