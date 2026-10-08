# Teleconnections and HF-low intensity (ERA5 proxy, pipeline A)

Question: do NAO and PNA change how intense hurricane-force lows get, once the index's own pressure signature and the shift in where storms sit are taken out?

Everything here is an **ERA5 proxy** built on **pipeline A** (`research/era5/hf_history`, 800 km gust index, HF-equivalent at 71.7 kt). The sample is the 4,157 catalog events, so this is intensity *given* a storm reached the HF-equivalent threshold, not intensity of all cyclones.

## Answer

**The teleconnections move HF lows and change the pressure background they sit in. They do not make them deepen faster, and only the Pacific shows a hint of stronger wind.**

Per SD of the index taken over days -10 to -4 before the track's first fix (so the storm cannot feed it). NAO and PNA are fitted together, with month fixed effects and a linear season term. p is from 1,000 season-block permutations, so the floor is 0.001. q is Benjamini-Hochberg over the 28 primary tests.

| | Atlantic, NAO | Pacific, PNA |
|---|---|---|
| Central pressure (hPa) | **-0.92** (p 0.003, q 0.017) | **-0.90** (p 0.001, q 0.009) |
| against climatology at that place and month | -0.07 (p 0.82) | **-0.91** (p 0.002, q 0.014) |
| against the 900-1100 km ring at that time | **-0.72** (p ≤ 0.002, q 0.009) | +0.12 (p 0.59) |
| background: ring minus climatology | +0.65 (p 0.011, q 0.051) | **-1.03** (p 0.001, q 0.009) |
| Fastest 24 h deepening (Bergerons) | -0.012 (p 0.43) | +0.010 (p 0.48) |
| Peak gust index, 2004-05 on (kt) | -0.01 (p 0.98) | +0.52 (p 0.037, q 0.13) |
| Hours at HF-equivalent, 2004-05 on | +0.02 (p 0.98) | +0.87 (p 0.092, q 0.29) |
| Deepest fix moves (deg) | +1.46 lat, +1.90 lon (p ≤ 0.002) | +1.15 lon (p 0.050), lat -0.13 |

**Atlantic.** NAO+ storms are about 0.9 hPa deeper per SD, but they are as deep as the climatology of where they end up: NAO+ moves them north and east toward the Icelandic low. Against their immediate surroundings they are 0.7 hPa deeper, or 1.0 hPa when position is also held fixed, but the surroundings are 0.6 hPa *higher* than climatology. One reading is a sharper pressure contrast across the storm. Yet the gust index, the HF-equivalent duration, and the deepening rate show nothing. A wind response, if there is one, is under the detection floor of about 0.65 kt/SD (2.8 × SE).

**Pacific.** PNA+ storms are about 0.9 hPa deeper per SD, and that holds against climatology at the same place. It is entirely the background: the storm-relative depth does not change (+0.1 hPa). The storm sits in a deeper Aleutian low. The gust index and HF duration rise (+0.5 kt, +0.9 h per SD) but do not survive the false discovery correction; treat that as a lead, not a finding.

**Neither index changes how fast storms deepen.** The deepening rate is near zero in both basins, with |t| < 1.

The contrast run with the index taken over the five days ending at the deepest fix gives the circular answer: Atlantic -1.93 hPa and Pacific -1.76 hPa raw, about twice the lagged values. Half of the same-time signal is the storm, or its parent pattern, written into the index.

## Effective sample size

Season-to-season variance of event intensity is tiny: the ICC is -0.01 to 0.05 for every outcome (`meta.json`). So the information is in the index's variation *within* seasons, about 5.3 (NAO) and 5.7 (PNA) independent values per season. That gives about 250 to 270 effective values over 47 seasons for pressure outcomes, and about 115 to 125 over 22 seasons for gust outcomes. The bootstrap and permutation resample whole seasons, so they respect this.

## Method

- `background.py` takes each event's in-domain fix of lowest pressure (it equals the catalog's `minp`) and adds two ERA5 MSLP values on the 1.5° grid:
  - `clim`: the calendar-month climatology at that position, from every 4th WeatherBench2 time chunk over 1979-01 to 2023-01 (1,244 to 1,372 fields per month).
  - `ring`: the cos-lat weighted mean over the 900 to 1,100 km annulus at the same time.
- Events from 2023-01-10 on (320 of them) use hourly ARCO coarsened with `research/era5/intensity/env.py`'s `coarsen()`. On 24 WeatherBench2 times in 2021 to 2022, ARCO and WeatherBench2 ring means agree to within 0.008 hPa (mean -0.001).
- About 6 GB was streamed in total. Only the monthly climatology was kept, in `work/`, which is ignored.
- `analyse.py` fits OLS per basin and outcome on NAO + PNA (z-scored), calendar-month fixed effects (month of the deepest fix), and a linear season term. The `lagged+position` variant adds lat, lat², lon, lon² of the deepest fix. Inference:
  - season-block bootstrap SE (1,000)
  - season-block permutation p (1,000; each event takes the index for the same day of season in a permuted season)
  - leave-one-season-out range
  - BH q over the primary rows: lagged index, no position, pressure outcomes on 1979+, gust outcomes on 2004+.
- Pressure outcomes use all 47 seasons, because depth did not drift before 2001 (`STATUS.md`). Gust outcomes use 2004-05 on (gust-drift gate and `RECORD_START`).
- Deepening rate is `(p(t) - p(t+24h)) / 24 × sin 60° / |sin φ|` over the whole track, positive for deepening.
- Transitioning tropical cyclones are **in**. Pipeline A does not mask them.

## Limits

- Ring-relative depth is not a pure storm measure. A storm on the flank of a large-scale low or high has the curvature of that pattern in its ring mean. That is a plausible reason the Atlantic ring depth and background move in opposite directions.
- Intensity given HF: if a teleconnection pushes marginal storms over the threshold, the conditional distribution can shift even without any change in the storms. Question 2 in `/mnt/project-files/science-questions/` (all cyclones, not just HF events) is the remedy and is queued separately.
- The CPC daily NAO and PNA files are read from a directory argument. They are the files the teleconnection test used, which match the repository's copies wherever they overlap. They are not committed here.

## Files

- `background.py`, `analyse.py`: the code.
- `results/summary.txt`: every fit (lagged, lagged+position, at-deepest), plus position shifts and ICC.
- `results/results.csv`, `results/position_shift.csv`, `results/meta.json`: the same in machine form.
- `results/event_table.csv`: the per-event table the fits use, including `clim` and `ring`.

Reproduce:

```
python3 background.py ../hf_history/results/era5_hf_catalog.csv ../hf_history/results/era5_hf_catalog_tracks.csv work/event_background.csv 32
python3 analyse.py ../hf_history/results/era5_hf_catalog.csv ../hf_history/results/era5_hf_catalog_tracks.csv work/event_background.csv CPC_DIR results 1000
```

## Verification and order of work

A fresh agent recomputed every number quoted above from `results/event_table.csv`, the catalog and the CPC files, with its own code. All of them agree. Its one caveat was that the permutation p for Atlantic ring depth came out at 2/1001 in its run, so it is quoted as ≤ 0.002. `results/HYPOTHESIS_AND_ORDER.md` records the hypothesis, the planned test, and the order of the exploratory look, the plan, and the run. It was written after the run.
