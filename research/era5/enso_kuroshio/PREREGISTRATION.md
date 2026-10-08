# Pre-registration: ENSO flavor and cyclogenesis over the Kuroshio (ERA5 proxy)

Written and committed **before any genesis count, SST index, or outcome was computed**. Anything changed afterwards is logged under "Post-hoc deviations" at the bottom, and both versions are reported.

Question (Jason, 2026-10-08): does an El Nino whose warm anomaly sits farther west (central-Pacific, "Modoki") lead to more cyclone genesis over the Kuroshio and Kuroshio Extension, and does that carry through to more Pacific hurricane-force (HF) lows forming there?

All cyclone numbers come from a **new MSLP-only low tracker that reuses pipeline A's detector and linker** (`research/era5/hf_history/extract.py::lows_at`, `track.py::link`), run 6-hourly Dec-Mar 1979-80..2025-26 on the ERA5 0.25 degree hourly store. HF labels come from **pipeline A** (`all_tracks.csv.gz`, gust index at 71.7 kt). ERA5 is a **proxy**, not the archive.

## What I had looked at before writing this

Only inventories (what files exist, their columns, chunk sizes). I read the existing Pacific frequency-split result (`freq_split/README.md`), which reports ONI against basin-wide Pacific cyclone and HF counts as null (RR 1.00-1.02 per SD, October-April 2004-05..2025-26). That is prior information that basin-wide counts do not respond to ONI; this study asks about one region and about flavor, which that test did not. No genesis position of any cyclone, no SST, and no flavor index had been computed.

## Why a new tracker run is needed

`all_tracks.csv.gz` keeps only each track's peak position, start time and summary; the first-fix position (genesis) is not in it, and `era5_hf_catalog_tracks.csv` holds full fixes only for the 8,311 catalog events and their time-matched null cases, not a population. Re-running the tracker is the only way to get genesis for every low.

Pull estimate: MSLP only, one global 2.17 MB chunk per 6-hourly step (measured mean of 48 chunks, SD 0.012 MB). Dec 1-Mar 31, 47 winters = about 22,800 steps = **49.5 GB**. SST for the indices: 4 samples a month, all months, 1979-2025 = about 2,260 steps x 1.37 MB = **3.1 GB**, a separate extraction. No cache covers either. Nov is left out because ENSO teleconnections are strongest in DJF and late winter and Nov would take the MSLP pull to 61.9 GB, over the gate in `CLAUDE.md`.

## Winters, windows, genesis

- Winter s runs Dec of year s to Mar of s+1 (matches pipeline A's June-May seasons). s = 1979..2025, **n = 47 winters**.
- **Primary window DJF**: genesis from 03 Dec 00Z to the last 18Z of February. The first two days are dropped because a low already present on 1 Dec would be mis-labelled as new. Counts carry an offset for window length (leap years).
- Secondary windows: **JFM** (1 Jan to 28 Mar 18Z) and **Dec-Mar** (03 Dec to 28 Mar 18Z). The last 3 days are dropped because a track needs 24 h (4 fixes) to exist and the pull ends 31 Mar 18Z.
- **Genesis** = the first fix of a track (pipeline A detector: smoothed MSLP minimum below 1010 hPa; link rules as in `track.py`; at least 4 fixes). It must lie over the ocean (ERA5 land-sea mask < 0.5), so a lee low that forms over land and crosses the coast is not a Kuroshio genesis. Tropical cyclones are left in, as in pipeline A.

## Kuroshio box (set from the literature, not from the data)

The Kuroshio runs from east of Taiwan through the East China Sea shelf break and south of Japan to the Kuroshio Extension east of 142E. The documented genesis sites on it are the East China Sea and south-of-Japan lows (about 25-35N, 120-142E) and Kuroshio Extension lows (about 33-42N, 142-165E).

- **Primary box K** = ocean points in A or B: A = 25-35N, 120-142E; B = 30-42N, 142-165E. (The Sea of Japan and the Yellow Sea are excluded by the latitude and longitude limits; land by the mask.)
- Sensitivity boxes: K shifted 5 degrees north (K_N), K shifted 10 degrees east (K_E), A alone (K_ECS), B alone (K_KE).
- Basin denominator for shares: all tracks with first fix over ocean in 25-67N, 120E-240E, same window.

## Flavor indices (ERA5 SST, 4 samples per month, anomalies from a 1991-2020 monthly climatology, DJF mean of Dec, Jan, Feb of winter s)

- Nino3 (5S-5N, 150W-90W), Nino4 (5S-5N, 160E-150W), Nino3.4 (5S-5N, 170W-120W); N4-N3 = Nino4 minus Nino3.
- **EMI** (Ashok et al. 2007) = SSTA[165E-140W, 10S-10N] - 0.5 SSTA[110W-70W, 15S-5N] - 0.5 SSTA[125E-145E, 10S-20N].
- **Primary flavor predictor: EMI.** Amplitude control: Nino3.4 (same SST source). I also report the correlation with the official ONI in the repo, as a check on the ERA5 SST. All predictors are standardised over the 47 winters before fitting.
- Winter classes for S12 (Kug et al. 2009 rule on the same indices): El Nino if Nino3.4 DJF >= +0.5; **CP** if Nino4 anomaly > Nino3 anomaly, **EP** otherwise. Neutral/La Nina are shown for reference only.

## Primary hypothesis and test

**H1 (directional):** with ENSO amplitude fixed, a more westward warm anomaly (higher EMI) gives more DJF genesis in box K: coefficient b_EMI > 0.

**P1.** Quasi-Poisson GLM, log E[N_s] = offset(log days) + a + b_EMI x EMI_s + b_34 x N34_s + b_yr x year_s, n = 47. Inference on b_EMI by **Freedman-Lane permutation** (10,000 permutations of the reduced-model residuals across winters), two-sided; the quasi-Poisson Wald interval is reported beside it. The winter is the unit; cyclone count is not the sample size.

## Pre-specified secondary tests (all will be reported; BH-FDR over the 15 tests of P1 and S1-S14)

| id | outcome | predictor | note |
|---|---|---|---|
| S1 | DJF share of basin genesis in K (binomial, quasi) | EMI | model as P1 |
| S2 | DJF K-genesis storms that later reach HF (pipeline A, 71.7 kt, in-domain) | EMI | **2004-05..2025-26 only** (22 winters), model without year |
| S3 | HF share among K-genesis storms | EMI | same 22 winters |
| S4 | P1 with window JFM | EMI | |
| S5 | P1 with window Dec-Mar | EMI | |
| S6 | P1 | N4-N3 replaces EMI | |
| S7 | P1 | EMI | box K_N |
| S8 | P1 | EMI | box K_E |
| S9 | P1 | EMI | box K_ECS |
| S10 | P1 | EMI | box K_KE |
| S11 | HF-reaching count, 1979-2025, with an era dummy (before/after 2001-02), EMI | within-era variation only, Decision 1 |
| S12 | El Nino winters only: CP vs EP mean DJF K genesis (permutation of class labels) | class | n of each class stated |
| S13 | P1 without the year term | EMI | |
| S14 | P1 without the Nino3.4 term | EMI | |

HF matching: a new track is HF if its fix at pipeline A's `peak_time` and `peak_lat/peak_lon` matches a pipeline A track with `gust800_kt` >= 71.7 and basin `pac`. **Validation of the tracker first**: at least 95% of pipeline A Pacific tracks whose peak_time falls in the pulled windows must be recovered this way; if not, I report the rate and stop before running the outcome tests.

**Mechanism (secondary, labelled, outside the FDR family).** M1: refit P1 adding DJF-mean PNA (CPC daily index in `cpc_indices/`); report b_EMI with and without it and the correlation of EMI with PNA. M2: the PNA coefficient. No jet-level pull is planned. The ONI x PNA question belongs to the thread "El Nino and PNA together" and is not repeated here.

## Decision rules (fixed now)

- **Yes** for the primary question: b_EMI > 0 with permutation p < 0.05 in P1, and the sign unchanged in S13 and S14.
- **No (well-powered null):** 95% interval for the effect per SD of EMI has an upper bound below +10% in count.
- **Cannot tell:** anything else. With 47 winters and a handful of flavor-defining events this is the likely outcome for S2, S3, S12; an underpowered null is reported as inconclusive.
- A result counts as holding up only if its BH q < 0.10 over the 15 tests. P1's own p is reported beside its q.
- **Power** is reported for every test as the effect per SD detectable at 80% (2.8 x SE from the fitted model); it uses no estimated effect.
- Maps (genesis density, 5 degree bins, per EMI SD and CP minus EP) are descriptive, with no significance claim.

## Not done on purpose

No re-choosing of box, windows, indices or controls after seeing counts. No Atlantic. No new pulls beyond the two above. Pipeline B is not used.

## Post-hoc deviations

(none yet)
