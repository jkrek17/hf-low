---
name: hf-storm-composites
description: Use when a question needs ERA5 fields around pipeline A cyclone tracks at a set lead time (composites, near-storm environment, storm-relative structure), to reuse committed tables and the proven pull methods before streaming anything.
---

# Storm-relative ERA5 composites on pipeline A tracks

Load `hf-harness` first. For pull mechanics also see the account skill `era5-data-acquisition` (its repo pointer is the old awips-tools path; use `research/era5/`). ERA5 **proxy**, pipeline A.

## Step 1: check what already exists (do not re-pull)

Committed on `claude/exciting-fermat-8vcvgq`, under `research/era5/`:

| File | Holds |
|---|---|
| `hf_history/results/all_tracks.csv.gz` | 75,087 tracks 1979-2025 (basin, gust800_kt, minp, peak time/position, start, end). No genesis positions. |
| `hf_history/results/era5_hf_catalog_tracks.csv` | Every 6 h point on event and null tracks. |
| `intensity/results/fixes_2004.csv.gz` | 159,430 fixes at 00/12 UTC, 2004-05 on: msl, dp12, age, speed, heading, g800, 24 h outcomes. |
| `intensity/results/env_2004.csv.gz` | Same fixes (join on track, time): Hart B, VTL/VTU, jet250, div300, vadv500, eady, sst, flux, tcwv. |
| `hf_structure/results/fixes.csv` | 5,983 HF-strength fixes with gust-field statistics, quadrants, RMG. |
| `hemispheric/results/weekly_table.csv.gz` | Weekly counts with the lagged hemispheric predictors. |
| `hf_vs_storm/results/` | HF (2,000) and storm-force (8,136) storms and anchors (`storms.csv`); 2,400 storm-scale anchors with structure scalars, lat, age (`storm/anchors.csv`); box-mean series (`large/series.csv`); `strat.py` (stratified difference, season bootstrap, BH). |

Raw fields in `$ERA5_WORK` are ignored and **gone after the session**. If a table covers the quantity, join to it; a lead time is a join to the fix at t minus lead on the same track.

## Step 2: if you must pull

Stores (anonymous GCS; code in `era5lib.py`, `hf_structure/extract.py`, `intensity/env.py`): ARCO-ERA5 0.25 degree hourly `full_37-1h-0p25deg-chunk-1.zarr-v3`; WeatherBench2 1.5 degree 6-hourly to 2023-01-09, ARCO coarsened after; 5.625 degree daily (`hemispheric/fields.py`, or `hf_vs_storm/fields_large.py` with water vapour).

Measured MB per time at 0.25 degree: gust 3.2, MSLP 2.2, u10 3.3, v10 3.3, dewpoint 2.4, t2m 2.3. Pulls are per-time files, resumable, reduced on the fly, never committed.

- Size by HEAD requests on a sample (`hf_vs_storm/size_pull.py`), not from memory. Size **every** part: the large-scale pull came in at 16.7 GB against a 12-14 GB guess.
- **The 50 GB gate counts all of a thread's pulls for the question together** (large + storm scale = 43.9 GB in the composites job; t2m was dropped to keep margin). Over 50 GB needs Jason's go-ahead: prepare, state the estimate and how you got it, stop. Do not split a multi-test plan into sub-50 GB pulls.
- Stream newest seasons first.

## Step 3: method conventions

- Seasons 2004-05 on for levels and fits (decision 1). Re-detect lows with pipeline A's rules (smoothed MSLP minima below 1010 hPa, 20-75N, 400 km de-duplication; ocean points owned by the nearest low) and **reproduce the catalog g800** before using any field (within 0.5 kt for 99% of anchors; max 0.1 kt last time).
- Motion from the track point 6 h before to 6 h after (previous 6 h only if the predictor must not see the future). Report how many short-lived tracks have no heading and drop out of rotated composites.
- Statistics: resample storms, not fixes; match groups on basin x month, stratified difference, season-block bootstrap, BH per family, MDE (2.8 x SE). Show the comparison group beside each composite; say transitioning tropical cyclones are in.
- Pre-register exclusion and missing-data rules (season window, missing lags, no heading) as well as groups; the composites plan omitted them (deviation 2) and the verifier's n differed by 2 to 14.
- Dry-run every test on a few times and assert no all-NaN output before the full run. The first run had every SST box mean NaN (land cells); the fix was an ocean-only, weight-renormalised mean. Log such fixes as post hoc deviations.
- **HF versus a weaker group:** gust, 10 m wind and HF-equivalent area differ by construction (HF is defined by the gust index): descriptive, outside the test family. Pre-register a **pressure-matched** variant (5 hPa bins of central MSLP) for size and structure; it cut the 48-kt radius gap from +158 km to about +40 (Atlantic), 0 (Pacific).
- Confounders to match or report with numbers: anchor latitude (HF onsets 2.6 degrees south in the Atlantic, 5.7 in the Pacific; dewpoint differences are raw, not latitude-adjusted), storm age at the anchor (+4.8 h at onset, +14 h at HF peak versus SF peak), owned-ocean fraction (+0.13).
- Terrain: no mask on pipeline A; flag gust maxima within 100 km of Greenland or Iceland, report with and without.
- ERA5 under-resolves peak wind; km2 areas are biased low, differences more trustworthy. Pixel maps are spatially correlated: descriptive.
- Test claims follow `hf-preregistered-test`; verify and close with `hf-result-closeout`. The verifier gets the pre-registration and committed tables and recomputes counts, means, differences **and** bootstrap p, q, MDE (skipped last time), for every comparison and secondary family. Extraction, heading, maps and figures stay unchecked: list them in `verify/VERIFICATION.md` and the README.
