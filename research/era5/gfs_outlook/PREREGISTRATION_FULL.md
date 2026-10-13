# GFS test of the full 24 h outlook model, and do the environment gains survive? Pre-registration (part 2)

Written and committed before any additional GFS field is pulled and before any full-model or ERA5-fed A/B comparison is scored. Extends
`PREREGISTRATION.md` (same dir; models S, A, B and their frozen coefficients `results/outlook_models.json`, commit 7870142, are not changed).
Pipeline A (ERA5) is a proxy: the training record and the truth. Jason's go: 2026-10-13 19:55Z ("1, 2, 3", item 1 of the whats-next list).

## Question
Part 1 left two things open: (1) the full model (state + Hart + all ten environment terms) was not run on GFS; (2) whether the A and B gains over S on GFS
(+0.0157, +0.0146 at lead 24) are as large as the same gains when the same frozen models are fed ERA5 on the same fixes was not measured.

## Models (frozen before scoring; added to `results/outlook_models_full.json`, committed alone)
Multinomial logistic, C = 1, ERA5 pipeline A fixes of seasons 2004-2020, `intensity/model.py` Prep, four-class gcls target, no `young`, `logage`:
- FULL = S + Hart (B, VTL, VTU) + jet250, div300, vadv500, eady, sst, sstgrad, sst_t500, flux, tcwv, nosst; fallback without dp12 as in part 1.
- FULL_NOSF = FULL without sst, sstgrad, sst_t500, flux, nosst (the terms whose GFS counterpart is least like ERA5's: GFS skin temperature for SST, a 6-hour average flux).
Climatology and S, A, B: part 1, unchanged.

## GFS fields (pull; same bucket and byte-range method; leads 0 and 24 of the same 2,122 cycles)
New records only, existing ones are not re-pulled (the 32.8 GB part-1 pull kept only a per-fix table, so these fields are not cached anywhere):
HGT 925, 700, 600, 400, 300 mb (Hart; 850, 500 already used), UGRD and VGRD 300 mb (div300), PWAT entire atmosphere (tcwv), TMP surface (SST proxy over ocean points: LAND = 0
and ICEC < 0.15), LHTFL and SHTFL surface (flux, lead 24 only: the 18-24 h average; no flux exists at lead 0), LAND and ICEC once. Probed `.idx` sizes: about 7.0 to 7.2 MB per
cycle at lead 0 and 8.1 to 8.7 MB at lead 24, so 2,122 x about 15.5 MB = about 33 GB (below the 50 GB gate; HARD STOP at 40 GB).
Predictors are computed exactly as `intensity/env.py` and `hart.py` on the 1.5 degree grid (as part 1); GFS HGT in gpm is used as ERA5 z/g; SST = TMP surface minus 273.15 over ocean;
flux = LHTFL + SHTFL (GFS fluxes are positive upward; ERA5's are positive downward and env.py negates, so the sign is already upward-positive: keep upward-positive). A compact per-fix table
goes to `/mnt/project-files/gfs-test/outlook-env/` (appended as new columns or a second file) with its README.
Lead 0 has no flux: the full-model score at lead 0 sets flux to the ERA5 training median (documented, secondary only). The primary full-model tests are at lead 24.

## Tests (family of 6, Benjamini-Hochberg q within the 6; season-block bootstrap, 2,000 draws, 5 GFS seasons, fixed seed; 90 % intervals; RPSS vs basin-month climatology; lead 24, pooled basins, matched fixes as in part 1)
- F1. FULL fed GFS: RPSS > 0, interval above 0.
- F2. FULL minus S, GFS-fed: gain > 0, interval above 0.
- F3. FULL, GFS-fed minus ERA5-fed (same fixes): equivalence, interval inside +-0.05.
- F4A, F4B, F4F. Gain retention for A, B, FULL: (gain over S when fed GFS) / (gain over S when fed ERA5, same fixes), season-block bootstrap of the ratio; "survives" if the 90 % interval lower bound is at least 0.5; "shrinks" if the upper bound is below 0.5; otherwise "cannot tell".
Secondary, labelled: the gains themselves (GFS-fed and ERA5-fed) with intervals; FULL_NOSF in the same comparisons (does dropping the SST and flux terms cost anything on GFS); each basin; lead 0; the gust-class reliability of FULL's HF probability on GFS.
Expectation written down: FULL survives with a lower score than the ERA5-fed FULL (probably about -0.03 to -0.05); A and B retain most of their gain (ratio 0.7 to 1.0); the SST and flux terms lose the most on GFS.
Truth, matching, dropped counts and dp12 fallback rules: part 1 (Deviations there apply).

## Limits stated in advance
Five GFS seasons (126 distinct season resamples, so intervals are narrow and q floors at 0.0006: the verdict is stability of sign across five seasons, not coverage of a future season), all inside the gust-calibration block; GFS GUST, skin temperature and averaged flux are the
model's own, not ERA5's; ERA5 is truth and a proxy; 38 and 54 kt are nominal; flux differs in definition at lead 24 (18-24 h average vs ERA5 hourly-mean); looks: one scoring run, logged in `research/era5/looks/gfs_test.log`.
Verification: a fresh Sonnet agent recomputes counts, FULL, FULL_NOSF, A, B and S RPSS (GFS-fed and ERA5-fed), and F1 to F4 point estimates from committed tables; intervals and the pull are not re-run.

## Deviations
Dated notes, each written before the step it affects. At the time of the first note nothing new had been pulled or scored; `results/outlook_models_full.json` was committed 8d2e1bf (before any pull).

- 2026-10-13, pull implementation and budget (before the pull):
  1. Part 1 kept only the per-fix jet250, vadv500, eady table, so the Hart terms and div300 / sst_t500 need HGT 850, HGT 500 and TMP 500 again. These three records are re-pulled at leads 0 and 24 with the rest. The plan's 33 GB estimate did not include them. Re-estimated from the `.idx` of 12 random cycles (both leads): new records 7.17 MB (lead 0) and 8.69 MB (lead 24) per file, re-pulled records 2.43 MB per file, so about 33.6 GB + 10.3 GB = about 43.9 GB over at most 4,244 files (plus the `.idx` files). That is above the registered 40 GB HARD STOP and below the 50 GB gate in CLAUDE.md. The stop for this pull is therefore raised to 46 GB (cumulative `bytes.log`, including every `.idx`); it is still a hard stop and the pull halts and reports if it is exceeded. No record is dropped to save bytes: all of them enter FULL.
  2. ICEC is read from every file (each lead) instead of once per cycle, so the sea-ice mask is the one valid at the field's own time; LAND once globally (`a_fetch.land_mask`, cached as `land.npy` in GFS_WORK). Ocean point: LAND == 0 and ICEC < 0.15. SST = TMP surface (K) at ocean points, NaN elsewhere, set to NaN BEFORE `env.coarsen` (which averages over the finite points and returns NaN for cells with under half the weight finite), so sst, sstgrad and sst_t500 behave as ERA5's land-masked SST does; `env.derived` then gives sst (C) and sstgrad (K per 100 km) unchanged.
  3. Sign of flux: GFS LHTFL and SHTFL are positive upward; they are coarsened and passed to `env.derived` as lhf = -LHTFL, shf = -SHTFL so its own negation returns flux = LHTFL + SHTFL, upward-positive, the same sign as the ERA5 training column. Lead 0 has none (flux = NaN). 18-24 h average, GFS fixed parameter units W m-2.
  4. Heading for Hart's B: as `intensity/fixes.py`, direction of motion over the previous 6 h of the GFS track, arctan2(dlon cos(lat), dlat) mod 360, from the fix at lead - 6 of the same cycle and track; NaN at lead 0, where `env.features` itself uses 90 degrees as for a first ERA5 fix. Nothing about the future of the track is used.
  5. `env.derived` / `env.features` are called on the fields as in part 1 (gcsfs and numcodecs stubbed). The inputs of the terms already tabulated in part 1 (u250, v250, u500, v500, u850, v850, t850) are given zero stand-ins and their outputs (jet250, vadv500, eady) are discarded; only B, VTL, VTU, div300, sst, sstgrad, sst_t500, flux, tcwv are kept. GFS HGT (gpm) is used as ERA5 z / g, PWAT (kg m-2) as tcwv.
  6. Missing files are counted per season and left NaN (a fix with a missing file has no term and drops out of the FULL tests, counted).
