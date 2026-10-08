---
name: hf-nwp-forecast-test
description: Use when scoring an ERA5-trained hf-low product (state model, pattern index) on GFS or GEFS forecasts, to freeze, pull, match tracks and judge degradation versus ERA5.
---

# Scoring an ERA5-trained product on operational or reforecast NWP

Load `hf-harness` first; the statistics layer is `hf-preregistered-test`, closing is `hf-result-closeout`. This skill adds only what is specific to feeding forecast fields to a model that never saw them. Worked example: `research/era5/gfs_test/` (README, PREREGISTRATION.md, `results/stageA.txt`); PR 128, prereg 95035fe, frozen model 27082e9. Pipeline A (ERA5) is a proxy: name it with every number.

## 1. Plan, freeze, then pull (order matters)

- Commit `PREREGISTRATION.md` before any field is pulled (95035fe). Draft plan `/mnt/project-files/gfs-test/PLAN.md` changed twice after the pipeline code was read (tracks need PRMSL+GUST every 6 h): read the tracker and model code before sizing.
- **Freeze before scoring.** `a_fit_model.py` fits on ERA5 pipeline A fixes of seasons 2004-2020 only and writes `results/gfs_model.json` (coefficients, `prep` median/clip/mu/sd, climatology). Commit it alone (27082e9) before any GFS fix is scored. Never refit, never tune. Drop predictors a forecast track cannot have (`young`, `logage`); use a no-`dp12` fallback model for lead 0 and for tracks younger than 12 h and report its share (13.7%).
- **Freeze the whole object.** The PR 41 pattern weights were committed as betas only; the EOFs and PC scalers lived in an uncommitted workdir, which cost a 7.7 GB ERA5 rebuild (deviation 5-6). Check that a frozen artifact can be applied from the repo alone before registering.
- **Check which seasons fitted the frozen thing.** B.2 scored 2004-2018 but the index weights were fitted on 2004-2014, so 11 of 15 seasons were in-sample (deviation 7, found after scoring). Register the clean subset (2015-18) as primary up front.
- Fix ambiguous scoring rules as dated prereg deviations before any score (deviation 9: model per lead, 00/12 UTC valid times, truth match, both cycles count as different leads).

## 2. Cheap pulls (anonymous S3, byte ranges)

- GFS 0.25 degree: bucket `noaa-gfs-bdp-pds`, `gfs.YYYYMMDD/CC/atmos/gfs.tCCz.pgrb2.0p25.fFFF`, `.idx` beside it. Read the `.idx`, take the byte range of one record (`PRMSL:mean sea level`, `GUST:surface`, `LAND:surface` once), decode with eccodes and assert shortName, 721x1440, scan flags. Code: `a_fetch.py` (`get`, `idx_rows`, `message`). Archive starts about May 2021, so 5 seasons (2021-22 to 2025-26). Process on the fly, store only lows and basin gust; work dir outside the repo (`GFS_WORK`, default `/home/claude/gfs_work`).
- GEFS v12 reforecast: bucket `noaa-gefs-retrospective`, `GEFSv12/reforecast/YYYY/YYYYMMDD00/c00/Days:1-10/hgt_pres_abv700mb_YYYYMMDD00_c00.grib2` (also `ugrd_`), 00 UTC only, ends 2019. Code: `b_fetch.py --probe` prints bytes first.
- **Sizes: probe real `.idx` ranges, then log actual bytes per file** (`bytes.log`). Plan said GEFS 11-13 GB; actual 1.345 GB (0.5 degree, single-level records, 7,567 responses). GFS PRMSL+GUST f000-f048 every 6 h, 2,122 cycles: 31.2 GB unique (34.9 GB logged with re-pulls) against a 31 GB estimate.
- The 50 GB gate counts the whole thread: 31.2 + 1.3 + 7.7 (ERA5 rebuild) = 40.3 GB.
- **Pull pitfalls.** Transient proxy 403/404 on files that exist: retry with backoff (6x for 404) before declaring Missing; count and report missing files (more than 2% per season must be reported; got 0). Stale `.idx` offsets in late Nov 2022 files: walk GRIB section-0 lengths (`message_walk`).

## 3. Forecast tracks and matching

- Run the pipeline A detector (`hf_history/extract.py lows_at` constants) and linker (`hf_history/track.py link`, MAXD 900 km per 6 h, `MIN_LEN` 4) on each forecast cycle separately, f000-f048 (`a_tracks.py`; `tests/test_a_tracks.py` plants a Gaussian low and gust patch, no network). GFS 0.25 degree is pipeline A's grid, so constants carry over. Tracks born after lead 36 drop out at MIN_LEN 4.
- Truth match (`a_score.py match_era5`): ERA5 pipeline A fix at the same valid time, nearest, within 500 km; no match counts as 0 and the count is reported (11,091 of 90,131). Secondary truth: any OPC archive HF fix within 600 km in (t, t+24 h]. Match on valid time, not on forecast cycle.
- Report counts at every filter step and have a fresh agent recompute them: the build report said 395k rows, truth was 318,612 (unfiltered 470,136; of these 151,401 have an empty basin).
- Cheap no-tracker check (A1): max ocean gust in the basin box vs "archive HF fix in basin within 3 h of valid time", logistic plus day-of-year harmonics, leave-one-season-out.

## 4. Degradation versus ERA5

- Feed the same frozen model the matched ERA5 state on the same fixes; the criterion is GFS-fed minus ERA5-fed BSS at lead 24 (A2.2: within 0.05). Result: -0.041 Atlantic, -0.061 Pacific.
- **Judge it by the interval, not the point.** Registered as a point rule, Atlantic passed while no interval sat inside +-0.05. Register an equivalence interval (both bounds inside the margin) or call it "cannot tell". State what ERA5-fed means when ERA5 `dp12` is missing (+0.377 vs +0.381).
- A "skill falls with lead" criterion (A1.2) failed: no measurable decay over 0-48 h (diff +0.005 [-0.004, 0.016]). Atlantic point values were ordered but the paired 0-vs-48 difference spanned 0, Pacific was unordered. Expect a flat curve to 48 h and register the criterion so a flat result reads as a finding.
- Recalibration (intercept and slope on the first season, applied later; A2.4) is a labelled secondary; it gained nothing here.

## 5. Seasons used for calibration, looks, reporting

- GFS seasons 2021-22 to 2025-26 are the gust-index calibration seasons: say so with every result; calibrate single thresholds on the first season only (A1.3, 2021-22), score the rest.
- The Stage B window 2015-18 is PR 41's held-out block: count each scoring as a look. Append one line per scoring to `research/era5/looks/gfs_test.log` (ISO time, tag, topic, seasons, overlap with calibration/held-out, number of code runs, "nothing changed between"). Timing or test runs of scoring code count.
- Season-block bootstrap (2,000 resamples, fixed seed), BH-FDR over the registered family (14 tests; 11 met criterion, 10 q < 0.05), planted-effect power.
- Malformed archive rows (e.g. date `20241101018`) are dropped; report the effect on labels.
- GFS model upgrades inside 2021-26 are not a reforecast; list them and split if one falls inside a season.
- Close with `hf-result-closeout`; the verifier did not check bootstrap intervals, power runs or reliability tables, so list them.
