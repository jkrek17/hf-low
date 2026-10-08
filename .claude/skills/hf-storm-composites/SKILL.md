---
name: hf-storm-composites
description: Use when a question needs ERA5 fields around pipeline A cyclone tracks at a set lead time (composites, near-storm environment, storm-relative structure), to reuse committed tables and the proven pull methods before streaming anything.
---

# Storm-relative ERA5 composites on pipeline A tracks

Load `hf-harness` first. For pull mechanics also see the account skill `era5-data-acquisition` (its repo pointer is the old awips-tools path; use `research/era5/` in jkrek17/hf-low). Everything here is an ERA5 **proxy**, pipeline A.

## Step 1: check what already exists (do not re-pull)

Committed on `claude/exciting-fermat-8vcvgq`, all under `research/era5/`:

| File | Holds |
|---|---|
| `hf_history/results/all_tracks.csv.gz` | 75,087 tracks 1979-2025 (track, basin, gust800_kt, minp, peak time/position, start, end, n_fix). No genesis positions. |
| `hf_history/results/era5_hf_catalog_tracks.csv` | Every 6 h point on event and null-case tracks (6 MB). |
| `intensity/results/fixes_2004.csv.gz` | 159,430 fixes at 00/12 UTC, 2004-05 on: msl, dp12, age, speed, heading, g800, 24 h outcomes (dp24, ndr24, hf24, hf48). |
| `intensity/results/env_2004.csv.gz` | Same 159,430 fixes (join on track, time): Hart B, VTL, VTU, jet250, div300, vadv500, eady, sst, sstgrad, flux, tcwv, sst_t500 within 500-1000 km, 1.5 degree. |
| `hf_structure/results/fixes.csv` | 5,983 HF-strength fixes with gust-field statistics, quadrants, RMG, stage. |
| `hemispheric/results/weekly_table.csv.gz` | Weekly counts with the lagged hemispheric predictors. |

Raw fields, boxes and per-season caches live in `$ERA5_WORK` (default `<module>/work/`), which is ignored and **does not survive the session**. Nothing raw is shared between sessions; if a table above covers the quantity, join to it. A lead time is a join: take the fix at t minus lead on the same track.

## Step 2: if you must pull

Stores (anonymous GCS, `gcsfs`/`urllib`):

- ARCO-ERA5 0.25 degree hourly, `full_37-1h-0p25deg-chunk-1.zarr-v3` (`era5lib.py`, `hf_structure/extract.py`). Gust and MSLP cost about 6.3 MB per time; the HF-structure pull was 33.3 GB for 5,311 times.
- WeatherBench2 1.5 degree 6-hourly (`intensity/env.py`, `seasonal_cycle/extract_eddy.py`) up to 2023-01-09; ARCO coarsened with `env.py coarsen()` after. 5.625 degree daily fields: `hemispheric/fields.py` (about 15 GB for 1979-2025; HTTP range requests fetch one level).
- Outputs are per-time or per-season files, resumable, reduced on the fly (boxes, ring means), never committed.

Rules: estimate GB before starting (measured bytes per time x times) and state whether a table above covers part of it. **Over 50 GB needs Jason's go-ahead**: prepare, state the estimate and how you got it, stop. Stream newest seasons first.

## Step 3: method conventions

- Seasons 2004-05 on for levels and fits (decision 1).
- Re-detect lows with pipeline A's own rules (smoothed MSLP minima below 1010 hPa, 20-75N, 400 km de-duplication; ocean points owned by the nearest low within the radius) and **reproduce the catalog g800** before using any field (HF structure matched 99.98% of fixes, max 0.1 kt).
- Motion from the track point 6 h before to 6 h after (previous 6 h only if the predictor must not see the future).
- Composite statistics: resample storms, not fixes; report the composite of null cases or the climatology beside it; say transitioning tropical cyclones are in.
- Terrain: no mask on pipeline A; flag fixes with the gust maximum within 100 km of Greenland or Iceland, report with and without.
- ERA5 under-resolves peak wind; areas in km2 are biased low, shapes and differences are more trustworthy.
- Test claims follow `hf-preregistered-test`; verify and close with `hf-result-closeout`.
