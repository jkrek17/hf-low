# Do the forecaster products work on operational model forecasts? Pre-registration

Written and committed **before any GFS or GEFS field was pulled or read**. Draft plan and size: `/mnt/project-files/gfs-test/PLAN.md`
(Jason chose "Run both", 2026-10-08 18:44Z). Two things changed from that draft after reading the pipeline code, both
recorded here before any pull: (1) track-level scoring needs 6-hourly forecast tracks, so Stage A pulls only PRMSL and GUST
(not 10 m wind or upper air) at every 6 h; wind-structure product 3 is therefore not tested here; (2) Stage B asks the question the
product actually poses (a forecast pattern for days 1-7 predicting archive counts in days 8-14), see below.
Revised size: Stage A about 31 GB, Stage B about 11 GB, total about 42 GB, under the 50 GB gate. No ERA5 is pulled.

Pipeline A = `research/era5/hf_history` (ERA5, a proxy). All ERA5-side numbers here are that proxy. Name the pipeline with every number.
The truth in this test is the OPC archive (primary where stated) or the pipeline A label (where stated). Nothing was scored yet.

## What is already known (background, not tested here)
Intensity state-only model, PR 12 (leave-one-season-out BSS 0.221 onset, full 0.423 HF24); pattern effect PR 41/85. This test
asks whether skill survives when the inputs are **GFS or GEFS fields**, which the models never saw.

## Data
**GFS** (NOAA open-data bucket `noaa-gfs-bdp-pds`, 0.25 degree, byte-range reads by `.idx`, same method as `jkrek17/wind-particles/server.py`):
`PRMSL:mean sea level` and `GUST:surface`, cycles 00 and 12 UTC, forecast hours 0, 6, ..., 48, October to April, seasons
2021-22 to 2025-26 (season label 2021..2025; five seasons; the archive starts about May 2021). `LAND:surface` once for the ocean mask.
If a cycle or hour is missing it is skipped and counted; more than 2% missing in a season is reported, not hidden.

**GEFS v12 reforecast** (`noaa-gefs-retrospective`): control member `c00`, `Days:1-10` files, variables `hgt_pres_abv700mb` (500 hPa
only) and `ugrd_pres_abv700mb` (250 hPa only), 00 UTC, 12-hourly steps, to day 7. Initialisations at 00 UTC on day S-8 for every week start S
(S = 1 October + 7k, k = 0..29) of seasons 2004-05 to 2018-19 (2004..2018; reforecast ends 2019): 450 initialisations.

**Archive**: `data/hf_lows/HF_Data_-_Atl.csv`, `HF_Data_-_Pac.csv` (fixes with category HF and DHF; HF = hurricane force; six-hourly).
**Pipeline A tables**: `intensity/results/fixes_2004.csv.gz` (ERA5 fixes with hf24), `hf_history/results/all_tracks.csv.gz`.
**PR 41 frozen pattern**: `hemispheric/results/frozen_primary.json` (index weights, discovery seasons 2004-14).

## Basin boxes (pipeline A, `track.py`)
Atlantic: 30-67N and (lon >= 262E or lon <= 10E). Pacific: 27-67N and 135-240E. Gusts are taken over ocean points (GFS `LAND == 0`) only.

## Stage A1: basin verification, no tracker
Unit: cycle x basin x lead tau in {0, 24, 48} h; valid time = cycle + tau. Label: any archive HF fix (category HF) in the basin box within +-3 h
of the valid time. Predictor G: maximum GFS GUST (kt) over ocean points in the basin box at the valid time. Model: logistic on G plus month
harmonics (cos, sin of day of year); comparison model month harmonics only. Leave-one-season-out over the five seasons. Scores: Brier skill
score (BSS) of the G model against the month-only model, and AUC; season-block bootstrap (2,000 resamples of 5 seasons; with five seasons
the intervals are wide and stated as such).
- A1.1 At tau = 24, BSS > 0, bootstrap 95% interval above 0, in each basin.
- A1.2 AUC is ordered tau 0 >= 24 >= 48 in each basin (tested as the paired bootstrap difference 0 vs 48).
- A1.3 A single G threshold chosen on 2021-22 only (the value that makes the share of cycles above it equal the label frequency that season, tau = 0)
  and then applied to 2022-23 to 2025-26: HSS at tau = 24 has a lower bootstrap bound of 0.30 or more (half of PR 12's 0.59).
Decision rule: "works" only if the A1.1 interval is above 0 in both basins; otherwise "cannot tell" if the interval spans 0 and power to see BSS 0.05 is
under 80% (planted-effect simulation), else "no".

## Stage A2: track-level HF onset probability on GFS forecast tracks
Tracks: pipeline A detector (`hf_history/extract.py` `lows_at`: 2-point Gaussian smoothing, 21x21 minimum filter, MSLP < 1010 hPa, 400 km
separation, gust ownership within 800 km over ocean) and linker (`hf_history/track.py` `link`) applied to each forecast cycle separately, f000-f048,
GFS grid 0.25 degree, 20-75N. A gust index `g800` is the ocean GUST maximum owned within 800 km.
Model: the PR 12 `state` model **without the history terms that a forecast track cannot have**: predictors `msl, dp12 (lead >= 12 h only),
lat, speed, g800, pac, doy_c, doy_s` (no `young`, `logage`); logistic, L2 C = 1, never tuned, standardised, fitted on **ERA5 fixes of seasons 2004-2020**
(fixes_2004.csv.gz, 00 and 12 UTC), target `hf24` (pipeline A gust index reaches 71.7 kt in (t, t+24 h]). Frozen and written to
`gfs_model.json` and committed before any GFS fix is scored. Applied unchanged to GFS forecast fixes at 00/12 UTC valid times, lead tau in {0, 12, 24, 36, 48} h.
Truth (primary): the pipeline A `hf24` of the ERA5 fix nearest the GFS fix at the same valid time within 500 km; no ERA5 fix within 500 km counts as 0
(the count of such fixes is reported). Truth (secondary): any archive HF fix within 600 km during (t, t+24 h].
Climatology reference: basin x calendar-month hf24 frequencies from the ERA5 training seasons.
- A2.1 Primary: BSS against climatology > 0 at tau = 24 (interval above 0, season bootstrap), in each basin.
- A2.2 GFS-fed vs ERA5-fed: on the same fixes, the same frozen model fed the matched ERA5 state; the BSS difference (GFS minus ERA5) is within 0.05
  in magnitude. This is the cost of using a forecast instead of reanalysis at that valid time.
- A2.3 BSS is non-increasing in tau across 0, 12, 24, 36, 48 (reported with intervals; no test beyond A2.1 for the slope).
- A2.4 (secondary) Recalibration: intercept and slope fitted on 2021-22 GFS fixes against the primary truth, applied to 2022-23 to 2025-26. Reported alongside A2.1,
  not as a replacement. Archive-truth BSS and HSS reported as secondary.
Notes fixed now: GFS GUST is the model's own gust parameterisation, not ERA5's; a level difference is expected and is why A2.4 exists. Hart
phase space and the environment predictors are not used (they need fields not pulled). The full PR 12 model is therefore not tested; the state-only model is
the cheapest operational form and is the lower bound.

## Stage B: pattern product (PR 41) with a forecast pattern
Index: PR 41's frozen weights and anomaly convention applied to GEFS control Z500 and U250 averaged over forecast days 1-7 from the 00 UTC initialisation
on day S-8 (so the forecast days cover S-7 to S-1, exactly the primary lag window), regridded to the 5.625 degree grid by area-weighted overlap; anomalies by
removing the GEFS-reforecast calendar-day mean of the same fields over the 15 seasons (a model-climatology anomaly, fixed in advance, so GEFS bias is removed
but not fitted to ERA5). The ERA5 index for the same windows comes from `oos_index_{atl,pac}.csv`.
Outcome: archive weekly counts (PR 41 weeks, 30 per season, Oct-Apr), seasons 2004-2018 (15 seasons; includes 2015-18 of PR 41's held-out seasons, so
the number of looks at those is stated).
- B.1 GEFS-forecast index vs ERA5 index for the same windows: Spearman correlation above 0.5 in each basin (ERA5 index is the leave-one-season-out index).
- B.2 Held-out deviance skill of count ~ month + previous-week count + GEFS index (Poisson, same machinery as PR 41, leave-one-season-out over 15 seasons) versus
  month + previous-week count: positive, interval above 0, in each basin. Reference: PR 41 on ERA5 index, +5.6% Atlantic and +2.3% Pacific (22 seasons).
- B.3 The conversion effect cannot be tested with GEFS control alone (no tracks); not attempted.
Decision rule as A1; power by planted effect at the PR 41 size.

## Multiplicity and reporting
Family: A1.1-A1.3 x 2 basins (6), A2.1 x 2 (2), A2.2 x 2 (2), B.1-B.2 x 2 (4): 14 primary tests; Benjamini-Hochberg q, and the count passing out of 14.
Secondary and recalibration results are labelled as such. A well-powered null is a full result; where the power is low it is "cannot tell". Look count at
2015-2025 (PR 41 held-out seasons) for Stage B: one (2015-2018 only); recorded in `research/era5/looks/gfs_test.log`. The GFS seasons 2021-22 to 2025-26 are
inside the block that calibrated the gust threshold; this limit is stated with every GFS result.
Known limits: five GFS seasons; GFS model upgrades inside the period are listed and, if one falls in a season, results are also reported without that season;
operational GFS is not a reforecast, so forecast model version changes over 2021-26 are not controlled; GEFS v12 is a different, coarser model from GFS.
Verification: a fresh Sonnet agent recomputes sample counts, the frozen model's BSS and AUC values and the Stage B correlations from committed tables;
numbers it did not recompute are listed.

## Deviations (post hoc)
Logged 2026-10-08 ~19:10Z, before any Stage B score or any Stage A score was computed.
1. **Stage B size.** GEFS reforecast records are 0.5 degree and one level only: actual pull 1.345 GB (7,567 responses), not the 11 GB estimated.
2. **Stage B fields.** 12 UTC step of each forecast day only (forecast hours 12, 36, ..., 156), the plain mean of seven days; no 00 UTC steps.
3. **Regridding.** GEFS 0.5 degree fields are linearly interpolated onto the 0.25 degree grid and then passed to `fields.regrid` unchanged (area-weighted onto 5.625 degrees).
4. **No detrend.** PR 41's per-cell linear trend over 1979-2021 ERA5 has no GEFS equivalent on 15 seasons; anomalies are model-climatology anomalies only (as written above).
5. **Frozen weights cannot be applied as committed.** `frozen_primary.json` holds only the ridge betas; the EOF vectors, mean map, PC standard deviations and latitude weights were pickled in an ERA5 workdir that is not in the repository, and the committed PC scores cannot be inverted. Stage B's index is therefore not computed yet. Planned fix (needs the thread's pulls to stay under 50 GB; decided after Stage A's real byte count): rebuild the primary-split EOFs by re-running `hemispheric/run.py freeze` (about 15-17 GB of ERA5, deterministic), then apply only the Z500 and U250 blocks of the frozen betas (the SST block is dropped because GEFS has no SST here; stated, and reported with a note that the ERA5 index with SST is a different quantity). B.1 compares the GEFS-based index with the same z+u-only index computed from ERA5 for the same windows (discovery-period fit), not with `oos_index`, because `oos_index` is a leave-one-season-out index from 22 separate fits and cannot be reproduced from one frozen fit.
