# Does the 24 h gale / storm / HF outlook survive GFS fields? Pre-registration

Written and committed before any new GFS field is pulled and before any outlook is scored on GFS. Pipeline A (ERA5, a proxy) is the training and
truth side everywhere; name it with every number. Follows `../bomb_env` (PREREGISTRATION.md 7c0669a; models, classes and outcome are defined there)
and `../gfs_test` (GFS tracks, matching, bootstrap machinery). Jason asked 2026-10-13 ("start the GFS").

## Question
The 24 h outlook (below gale <38, gale 38-54, storm 54-71.7, HF >= 71.7 kt on the pipeline A gust index g800, peak over the next two 00/12 UTC
fixes) was fitted and scored on ERA5 only (perfect prognosis). Fed GFS forecast fields instead, how much skill is left, and does the upper-level
environment still add anything to the storm's own state?

## Models (frozen before scoring; written to `results/outlook_models.json`, committed alone before any GFS fix is scored)
Multinomial logistic, L2 C = 1, never tuned, `intensity/model.py` Prep (median impute, clip 0.5-99.5 pct, standardise), fitted on ERA5 pipeline A fixes of
seasons 2004-2020 only (00/12 UTC), target the four-class outcome (gcls). Predictors without the history terms a forecast track cannot have
(no `young`, `logage`): S = msl, dp12, lat, speed, g800, pac, doy_c, doy_s; A = S + jet250 + vadv500; B = S + jet250 + eady. A fallback without dp12 is
frozen for each (lead 0 and tracks younger than 12 h; share reported). Not testable from the fields pulled: the Hart terms, div300, SST, SST gradient,
flux, tcwv, so the full model is not tested; S, A and B are.
Climatology reference: basin x month class frequencies from the ERA5 training seasons.

## Data
- GFS forecast tracks already committed: `../gfs_test/results/gfs_tracks_2021..2025.csv.gz` (2,122 cycles 00/12 UTC, Oct-Apr, seasons 2021-22 to 2025-26, 6-hourly leads
  0-48; pipeline A detector and linker on GFS PRMSL and GUST). No new pull is needed for S.
- New pull for A and B only: GFS 0.25 degree, same bucket and byte-range method, leads 0 and 24 h of every cycle, records UGRD and VGRD at 250, 500 and 850 mb,
  TMP at 500 and 850 mb, HGT at 500 and 850 mb (10 records). Probed `.idx` sizes: 7.5 to 8.3 MB per cycle-lead, so about 2,122 x 2 x 8 MB = 34 GB (below the
  50 GB gate; no existing cache holds these fields). Fields are decimated to the 1.5 degree grid and the three terms computed exactly as `intensity/env.py`
  (jet250 within 1000 km, vadv500 and eady mean within 500 km) at each GFS fix; only a per-fix table (cycle, lead, track, jet250, vadv500, eady) is kept.
  The compact table is saved to `/mnt/project-files/gfs-test/outlook-env/` with a README. Leads 12, 36 and 48 are not pulled; the environment models are scored at leads 0 and 24 only.
- Truth: the pipeline A four-class outcome (gcls) of the ERA5 fix nearest the GFS fix at the same valid time within 500 km (no match within 500 km: dropped
  from the outlook score and counted; the A2 convention counted no match as 0 for hf24, here a missing storm has no class, so it is dropped and reported).

## Tests (family of 5, Benjamini-Hochberg q within the 5; season-block bootstrap 2,000 draws over the 5 GFS seasons, fixed seed, 90 % intervals)
Score: ranked probability score skill (RPSS) against the basin-month climatology.
- G1. S fed GFS, lead 24, pooled basins: RPSS > 0, interval above 0.
- G2. S, lead 24: (GFS-fed minus ERA5-fed, same matched fixes) RPSS, equivalence: interval inside +-0.05 (both bounds). Otherwise "cannot tell" or "degrades".
- G3. A over S, GFS-fed, lead 24: RPSS gain > 0, interval above 0. G4. B over S, same. 
- G5. S fed GFS: RPSS ordered lead 0 >= 24 >= 48, tested as the paired bootstrap difference lead 0 minus lead 48 (a flat curve is a finding).
Secondary, labelled: S at leads 12 and 36; each basin; recalibrated outlook (class intercepts fitted on 2021-22, applied 2022-23 on); reliability of the GFS-fed
HF probability; GFS-fed A and B gain at lead 0. Decision rules as `hf-preregistered-test`: a well-powered null is a result; with 5 seasons the intervals are wide
and the page of results says "cannot tell" where they span the margin.
Expectation written down: S survives with a lower score than ERA5 (GFS gusts differ in level), the A/B gains survive in sign but shrink; a flat skill-versus-lead curve to 48 h as in `../gfs_test` Stage A2.

## Limits stated in advance
Five GFS seasons, all inside the block that calibrated the 71.7 kt gust cut. Operational GFS is not a reforecast; model upgrades over 2021-26 are listed and not controlled.
GFS GUST is the model's own parameterisation, not ERA5's. Truth is ERA5 (a proxy), not observations. 38 and 54 kt are nominal. Fixes within a storm are not independent.
Looks: GFS seasons are not PR 41's held-out block; the ERA5 side (2004-2020 training) is not scored here. One scoring run of the family; later runs are logged in `research/era5/looks/gfs_test.log`.
Verification: a fresh Sonnet agent recomputes sample counts, S, A and B RPSS values and the G1 to G4 point estimates from the committed tables; intervals and the pull itself it is not asked to repeat.
