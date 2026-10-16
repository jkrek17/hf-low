# Gulf of Alaska jets, follow-up: mean wind instead of gust (lead 1) and GFS forecast pressure (lead 2)

Written and committed BEFORE any mean-wind number and any GFS field was looked at. Asked by Jason 2026-10-16 19:25Z ("1 and 2. Worth seeing how
ascat looks?") after the study in `PREREGISTRATION.md` (PR 255). Everything about ERA5 is a PROXY. Gust is the shared-cache ERA5 gust (pipeline: this study's
own zone maxima, not pipeline A or B). The first study's Cook Inlet result (sea warmer than air, 13 strong episodes) was the largest of 15 tests and is a
lead; this follow-up is labelled POST HOC for that reason.

## Lead 1: is the Cook Inlet sea-minus-air result a gust-formula artefact? (no new pull)

Data: the 1,800 fields already pulled (`/mnt/project-files/goa-jets/fields_1800_times.tar`). Sample: the same Z3 GAP episodes as `results/episodes.csv`
(regime assigned at the peak-gust time, as before). For each pulled time: **WMAX** = largest ERA5 10 m mean wind speed (kt, sqrt(u^2+v^2)) over Z3 ocean
cells; **GF** = Z3 zone-maximum gust / WMAX at that time (ratio of two maxima, which need not be the same cell).

Tests (family of 5, BH; ridge-penalised 0.5 as in `analyse.py`; season-block bootstrap 2,000, seed 16; standardised factors; each adjusted for DEPTH):
1. STAB -> **wind-strong** (top 13 of the 109 episodes by WMAX at the peak time, ties by gust; 13 = the strong-by-gust count, so the share matches), logistic.
2. INTERIOR_DP -> wind-strong, logistic.
3. STAB -> WMAX (kt), linear (ordinary least squares, bootstrap).
4. INTERIOR_DP -> WMAX, linear.
5. STAB -> GF, linear (does the gust factor itself move with sea-air difference).

Reading rule (fixed now): STAB is called **not a gust artefact** if test 1 has OR >= 1.5 with an interval excluding 1 AND test 3 has a positive slope with
an interval excluding 0. If only tests 1/3 fail and test 5 has an interval excluding 0, it is **partly the gust formula**. Otherwise **cannot tell**
(13 strong episodes). The overlap between the wind-strong set and the gust-strong set is reported.

## Lead 2: do GFS forecast pressure patterns at 24 and 48 h separate event times from non-events?

Purpose: a table a forecaster can use from the model, without waiting for the analysis. Outcome (truth, a PROXY): ERA5 zone gust >= 54 kt (event) and >= 64 kt
(strong, descriptive only) at the valid time, from `results/stage0_zone_gust.csv` (known for every 00/12 UTC time). Valid times: every 00 and 12 UTC,
October to April, seasons 2021-22 to 2025-26 (GFS archive starts about May 2021): 5 seasons x 212 days x 2 = about 2,120 times. Forecasts: GFS 0.25
degree, the cycle 24 h and 48 h before the valid time (12 UTC cycle for a 12 UTC valid time, 00 for 00), records PRMSL (mean sea level), UGRD and VGRD (10 m),
GUST (surface), by `.idx` byte ranges from `noaa-gfs-bdp-pds` with the repository's reader (`grib2_t53.py`). Window cropped to 44-70N, 175-125W, the same
rows and columns as the ERA5 study. **Pull estimate: 4 records (0.97+0.94+1.0+0.61 = 3.5 MB per file, probed from real `.idx` files) x 2 leads x about 2,120
times = about 14.8 GB, under the 50 GB gate; no existing cache holds these fields** (`gfs-raw-gust` and `gfs-daily-max-gust` hold gust only, coarser or
daily-max). Missing files are counted and reported.

Predictors, computed from the GFS forecast exactly as `features.py` did for ERA5 (same boxes, same land mask, same definitions): forecast regime from the
zone-mean 10 m wind direction; DEPTH; LOW_DIST (to the forecast zone-maximum gust location); COASTPILE (Z1, Z2) or INTERIOR_DP (Z3); GFSGUST = forecast zone-maximum gust
(kt; GFS gust is a different index from ERA5 gust and is not corrected). STAB is NOT used (needs SST; GFS surface skin temperature is not pulled).

Analysis (unit: valid time; independent sample: **5 seasons**, so every interval is wide and is reported as such):
- **T1 table.** For each zone, lead and forecast regime: chance of an ERA5-proxy event (and strong) given the forecast band of each ingredient. Bands: bottom half /
  50th-80th / 80th-95th / top 5% of the forecast-ingredient distribution within the forecast regime over all 2,120 times (cuts fixed from forecasts, no outcomes). Plain
  rates. Season-block bootstrap intervals (5 blocks) plus the range across the five seasons.
- **T2 skill (registered family of 6: Z1, Z2, Z3 x 24 h, 48 h).** Within the forecast barrier regime (Z1, Z2) or forecast gap regime (Z3), leave-one-season-out
  ridge logistic of event on (a) GFSGUST alone versus (b) GFSGUST + DEPTH + LOW_DIST + COASTPILE (or INTERIOR_DP). Statistic: out-of-sample log-loss (b) minus (a), per
  season-block bootstrap (5 seasons, 2,000 resamples), BH over 6. The pressure ingredients are said to **add** only if the interval excludes 0 with (b) better. A
  well-powered null is a result; with 5 seasons a tie is more likely "cannot tell" and will be called that. Also reported: AUC of (a), (b) and of the climatology-by-regime.
- **T3 same-ordering check (descriptive).** Do the forecast bands order the same way as the ERA5 same-time table of the study (`q1_table.csv`)?

Caveats stated now: (1) the ingredients were chosen from the ERA5 study, whose 22 seasons include these 5, so this is not independent of the ingredient choice, only of
the field source; (2) GFS model changes inside 2021-26 are not reconstructed here (listed by their release notes only if checked); (3) truth is the ERA5 proxy, not
observation; (4) GFS 0.25 degree smooths narrow straits like ERA5, so Cook Inlet and Shelikof remain under-resolved; (5) the 5 seasons are also the gust-calibration
seasons of the outlook work. One scoring run; looks are logged in `research/era5/looks/goa_jets_followup.log`.

## ASCAT: sizing only, not run

Count first what the existing ASCAT cache (`/mnt/project-files/ascat-check/cells`, cells within 800 km and +-1.5 h of a pipeline A fix) already holds in the Gulf of
Alaska zones at episode times, and what a fresh RSS pull for those days would cost. Nothing is pulled without Jason's word beyond that count.

## Deviations

- 2026-10-16, fixed in the analysis script before the scoring run (POST HOC wording, found by the claim review): GFSGUST was added to the T1 band table (it is a registered predictor but T1 listed only ingredients);
  the T2 verdict label "cannot tell / no gain" differs from the registered wording, and a significantly worse result (Cook Inlet 48 h, q 0.003) had no registered label and is reported as "worse"; T3 was not run;
  the climatology AUC was not reported (0.5 within a regime). Lead 1: p values are floored at 1/2,000 (no bootstrap draw on the other side).
