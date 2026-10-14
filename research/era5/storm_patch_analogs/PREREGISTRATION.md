# Storm-centred pattern score and analogs: pre-registration

Status: **plan only. No field has been pulled and no patch has been read.** Written 2026-10-14 for Jason's question
"take a 500 nm square around the storm, match it against previous HF lows, and score the comparison". ERA5 proxy,
pipeline A, WeatherBench2 1.5 degree fields. Every score here is against the proxy's own gust index, not the archive.

## Question

Does the shape of the weather pattern around a low, taken as a whole, tell us anything about its next 24 h (gale,
storm or hurricane force) beyond what the existing storm-state and near-storm summary numbers already tell us? And
can that pattern be used to retrieve similar past storms ("analogs") whose outcomes are a usable forecast?

Prior (written before any data): small. The upper-air and trough groups added -0.006 to +0.002 BSS over the trees
(RA-20, `../explosive_onset`), jet times vorticity advection and jet times Eady added +0.016 and +0.015 on the GFS
outlook (PRs 178/179), and perfect knowledge of the next 24 h of deepening adds only +0.07 BSS (`../hf_skill_limits`).
My expectation is +0.00 to +0.015 outlook RPSS. A null is a full result.

## Population (outcome counts were read; no pattern was)

- Pipeline A fixes at 00/12 UTC, joined from `../intensity/results/fixes_2004.csv.gz` and `env_2004.csv.gz`, the
  same table and loader as `../bomb_env/score.py` (`model.load`).
- Months October to April, seasons 2014-15 to 2021-22 (8 seasons, the latest complete ones in WeatherBench2, which
  ends 2023-01-09). May to September is left out: it holds almost no HF and would cost about 40% more pull.
- **32,288 fixes on 7,379 tracks.** Outcome classes (next 24 h peak gust, the PR 170 definition: below gale <38 kt,
  gale 38 to <54, storm 54 to <71.7, hurricane force >=71.7): 0.381 / 0.352 / 0.225 / 0.041; 1,330 hurricane-force
  outcomes on 513 tracks (291 Atlantic, 222 Pacific). Fixes per season 3,878 to 4,159.
- The independent sample is **8 seasons**. Tracks and fixes within a track are not independent.

## Data and pull

- WeatherBench2 `era5/1959-2023_01_10-6h-240x121_equiangular_with_poles_conservative.zarr`: u and v at 250 hPa,
  geopotential at 500 hPa, mean sea-level pressure, at the fix time. One chunk file per variable holds 8 six-hourly times
  at all 13 levels.
- Climatology: WeatherBench2 1990-2019 6-hourly by day of year (`../hf_env_composites/clim.py`, `../hf_basin_composites/clim_msl.py`).
- **Size estimate: 853 two-day chunks (counted from the fix times) x 29.08 MB (the measured u, v, z, msl cost per chunk
  in `../hf_basin_composites`) = about 24.8 GB, plus 0.5 GB climatology: about 25.3 GB, under the 50 GB gate.** Estimate,
  not measured on this run; the extractor stops at 30 GB.
- No existing cache covers it: the earlier storm-centred pulls kept only sums and means (`hf_env_composites/results/comp_north.npz`).
  The shared-folder `near-storm-environment/` holds figures only.
- Output: one compact file per season in `/mnt/project-files/storm-patches/` (float16, about 90 MB in all) with a README,
  not committed to git. Scripts and small result tables are committed.

## Patch

Centre = the fix; north-up (earth-relative, not rotated to motion); 13 x 13 points at 154.3 km spacing, so the
edge is 926 km = 500 nautical miles from the centre; bilinear sampling from the 1.5 degree grid
(`../hf_env_composites/common.py` `sample`, `destination`). Channels, as anomalies from the 1990-2019 climatology:
500 hPa height (m), sea-level pressure (hPa), 250 hPa u and v (m/s). Raw 500 hPa height and sea-level pressure are
stored too, for display only. 4 x 169 = 676 numbers per fix. Heading-rotated patches are **not** in this plan.

## Models (all leave-one-season-out within the 8 seasons)

Same logistic form as `../bomb_env/score.py`: multinomial, L2, C = 1, `model.Prep` standardising, the four-class
outlook, RPSS against basin-month climatology.

- **BASE = "full"**: the storm-model predictors (`model.SETS["full"]`), refitted on the 8 seasons. Reference only:
  the same fixes under the 22-season fit from `../bomb_env/results/oof_outlook_full.csv.gz`.
- **PC15 features**: the first 15 principal components of the patch (each channel scaled by its training standard
  deviation; scaler and PCA fitted on the training seasons only), added to BASE.
- **Primary test T1:** outlook RPSS of BASE + PC15 minus BASE, pooled over fixes, 90% interval by resampling the
  8 seasons (2,000 draws, seed 20261014). One primary test.
  - "Pattern adds skill": gain >= +0.010 and interval lower bound > 0.
  - "No useful gain": interval upper bound < +0.005.
  - Otherwise "cannot tell". I do not claim a null as well-powered unless the upper bound is below +0.005.

## Secondary tests (BH-FDR over S1 to S3; each is a gain over BASE, same interval method)

- **S1 environment only:** PC15 from the patch with the inner 300 km removed (points at r >= 300 km), so the core
  circulation, which the state numbers already hold, is not compared.
- **S2 composite projection:** one feature, the dot product of the standardised patch with the training composite
  difference (mean patch of hurricane-force-outcome fixes minus the mean over all fixes), recomputed within each fold.
- **S3 analogs:** k = 30 nearest neighbours in the whitened PC15 space; the library excludes the test season, the
  fix's own track and (for training rows) the row's own season; features = the neighbours' outcome-class fractions
  (3 free). Training rows get the same leave-own-season-out neighbours, nested inside the fold.
- **S4 placebo (reported, not in the FDR family):** PC15 from patches permuted among fixes of the same basin, calendar
  month and hour (seed 20261014). Its gain shows the over-fit cost of 15 extra columns; the T1 gain must exceed it.

## Descriptive, not tested

- State-only plus PC15 against state-only (how much of the pattern is already in the summary numbers).
- Gain by class (hurricane-force Brier skill), by basin, and for fixes whose track continues 24 h.
- **Analog calibration:** the stand-alone k = 30 analog probability of hurricane force: Brier skill against climatology,
  reliability slope and a decile table, to say whether "6 of 20 closest storms reached HF" means about 30%.
- **Analog display (no skill claim):** five nearest other-season analogs (track, time, outcome class) for every fix
  with a hurricane-force outcome, and storm-by-storm case pages for the 6 storms with the lowest minimum pressure in
  the 8 seasons, chosen before looking at their analogs.

## Limits stated now

- ERA5 proxy, perfect-prognosis: ERA5 patches in, ERA5 gust outcomes out. A live product would use GFS patches and a
  tracker; the earlier full-model test kept about 85% of the ERA5 gain on GFS, but that is a later step.
- 8 seasons; the 8-season subset is a convenience sample of the latest complete seasons. If T1 shows signal I will
  ask before extending to 2004-05 on (about 58 GB by the same estimate, over the gate).
- All 22 seasons were already used in earlier leave-one-season-out looks, so there is no fresh held-out data. This is
  one more look at them, logged in `../hemispheric/results/heldout_looks.log` when run.
- 1.5 degree grid: synoptic pattern, not mesoscale structure. Tropical cyclones are in. Units are nautical miles for
  the display and kilometres in the code.
- Fix selection requires `track` and `time` to match `env_2004.csv.gz`; fixes with missing environment are dropped as in `model.load`.

## Deviations log

(empty; any change after the first field is read goes here with the reason)

## Status at commit

Blocked before the pull: this session's network policy denies `storage.googleapis.com` (WeatherBench2) and `pypi.org`,
and scipy, scikit-learn, joblib, gcsfs, numcodecs and matplotlib are not installed. Nothing has been run.
