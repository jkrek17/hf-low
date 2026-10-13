# GFS_outlook: does the 24 h gale / storm / HF outlook survive GFS fields?

Pre-registration `PREREGISTRATION.md` (plan 2ac631b, written before any pull or score; dated Deviations appended before scoring). ERA5 pipeline A is a **proxy**:
it is both the training record of the models and the truth here; nothing below is validated against observations. Pipeline A numbers throughout.

## What was done
1. `fit_models.py` (commit 7870142, before any GFS field was pulled): multinomial logistic, L2, C = 1, `intensity/model.py` Prep, on ERA5 pipeline A 00/12 UTC fixes of seasons 2004-2020 (123,105 rows; 100,386 with dp12), four-class target gcls (peak g800 at +12 and +24 h, cuts 38/54/71.7 kt, `bomb_env/score.load`). S = msl, dp12, lat, speed, g800, pac, doy_c, doy_s; A = S + jet250 + vadv500; B = S + jet250 + eady; each with a no-dp12 fallback. `results/outlook_models.json`.
2. `fetch_env.py`: GFS 0.25 degree, leads 0 and 24 of the 2,122 cycles of `../gfs_test/results/gfs_tracks_2021..2025.csv.gz`; 10 byte-range records per file; decimated with `intensity/env.py` `coarsen` to 1.5 degrees; jet250, vadv500, eady at all 69,174 GFS track fixes of those leads (27,943 at lead 0, 41,231 at lead 24). Pulled 32.84 GB over 4,244 files (below the 40 GB stop and the 50 GB gate); **0 missing files** in every season; 0 NaN terms. Table: `results/gfs_outlook_env.csv.gz` (1.1 MB) and `/mnt/project-files/gfs-test/outlook-env/`. As a plausibility check, jet250 mean over these GFS fixes is 108.6 kt against 108.4 kt in the ERA5 env table for 2021-25 fixes (different fix sets, not a validation); vadv500 mean 0.32 vs 0.27 and eady mean 0.73 vs 0.68, i.e. GFS 7-16 % higher on the two derived terms. The GFS fixes differ from the ERA5 ones (forecast leads), so this is not a like-for-like comparison and the difference is not explained here.
3. `score_outlook.py`: one scoring run of G1-G5 and the labelled secondaries. RPSS against basin x month climatology of the training seasons; season-block bootstrap over the 5 GFS seasons (2,000 draws, seed 20261013, 90 % intervals), BH over the 5. A plumbing dry run on a partial environment table preceded it; its outputs were not read except the count lines.

## Counts (all from `results/outlook.txt`)
GFS track rows 318,612; Atlantic/Pacific fixes at valid 00/12 UTC, leads 0/12/24/36/48: 90,131 (lead 0: 15,202; 12: 19,291; 24: 21,149; 36: 19,255; 48: 15,234). No pipeline A fix within 500 km at the same valid time, dropped: 11,091 (lead 0: 1,211; 12: 2,075; 24: 2,982; 36: 2,667; 48: 2,156). Matched 79,040 (median distance 82 km; none matched across basins). Matched truth classes at lead 24 (below gale, gale, storm, HF): 6,220 / 6,745 / 4,366 / 836 (n 18,167); lead 0: 3,119 / 5,966 / 4,087 / 819. Of the matched ERA5 fixes, 24,839 have a track ending within 24 h (class 0 by the `bomb_env` convention). Environment terms present for all matched lead 0 (13,991) and lead 24 (18,167) fixes. No-dp12 fallback share at leads 12/24/36/48: 0.178 / 0.157 / 0.098 / 0.000 (lead 0 is always fallback; speed is missing at lead 0 and imputed). ERA5-fed S used the dp12 model for 54,357 of 65,049 lead >= 12 fixes.

## Results (RPSS vs basin-month climatology; point [90 % interval])
| Test | Estimate | n | q | Verdict (rule in Deviations 5) |
|---|---|---|---|---|
| G1 S fed GFS, lead 24, RPSS | +0.3147 [+0.3105, +0.3194] | 18,167 | 0.0006 | yes |
| G2 S lead 24, GFS-fed minus ERA5-fed RPSS | -0.0286 [-0.0311, -0.0264] | 18,167 | 0.0006 | equivalent (inside +-0.05) |
| G3 A minus S, GFS-fed, lead 24 | +0.0157 [+0.0142, +0.0174] | 18,167 | 0.0006 | yes |
| G4 B minus S, GFS-fed, lead 24 | +0.0146 [+0.0129, +0.0163] | 18,167 | 0.0006 | yes |
| G5 S lead 0 minus lead 48 | -0.0291 [-0.0424, -0.0152] | 27,069 | 1.0 | reversed (skill is higher, not lower, at 48 h) |

S RPSS by lead, GFS-fed: 0 +0.3024; 12 +0.3168; 24 +0.3147; 36 +0.3250; 48 +0.3315. ERA5-fed on the same fixes: 0.3362, 0.3371, 0.3433, 0.3586, 0.3751. ERA5-fed lead 24 RPSS +0.3433 [+0.3400, +0.3465]. The GFS-fed deficit is 0.020 to 0.044 across leads (0.034, 0.020, 0.029, 0.034, 0.044), largest at 48 h.

Secondary (labelled, no multiplicity control): lead 24 S RPSS Atlantic +0.2832 [+0.2638, +0.3025], Pacific +0.3441 [+0.3278, +0.3618]; A minus S Atlantic +0.0117, Pacific +0.0193; B minus S Atlantic +0.0128, Pacific +0.0163 (all intervals above 0). At lead 0: A minus S +0.0385 [+0.0363, +0.0413], B minus S +0.0429 [+0.0409, +0.0452]. Recalibration (class offsets fitted on 2021-22, applied 2022-2025 at lead 24): RPSS +0.3128 to +0.3209, gain +0.0081 [+0.0062, +0.0102] (4 seasons). Reliability of the GFS-fed S hurricane-force probability at lead 24: slope 0.96, mean forecast 0.0427 against base rate 0.0460; the top decile forecasts 0.325 against observed 0.313, the middle deciles under-forecast (e.g. forecast 0.0035 against observed 0.0105). Mean forecast class probabilities at lead 24 under-forecast gale (0.343 vs 0.371) and over-forecast below gale (0.389 vs 0.342).

Reading: S keeps most of its ERA5-fed skill when fed GFS (RPSS 0.31 against 0.34, a loss of 0.03 that sits inside the registered +-0.05 margin), and the upper-level terms add a small gain over the storm's own state (about +0.015 RPSS at lead 24, +0.04 at lead 0). The A/B gains were not computed ERA5-fed on this sample in part 1, so whether they shrank against ERA5 was not tested there (part 2, section below, does it). The expected flat-or-falling skill-versus-lead curve is not what was found: RPSS rises with lead. A likely contributor (not tested) is that lead 0 uses the fallback model with imputed speed, and that fixes at later leads are a different, longer-lived subset of storms with a different climatology reference; the comparison across leads is between different samples.

## Deviations
See the dated notes in `PREREGISTRATION.md`: scoring rules the plan left open (model choice per fix, the ERA5-fed rule, G3/G4 sample, p-value and verdict definitions, recalibration), and the pull implementation (eccodes installed with pip; `env.py` imported with its zarr readers stubbed). No deviation from the registered tests, data or seeds otherwise.

## Limits
- Five GFS seasons, all inside the block that calibrated the 71.7 kt gust cut; not independent of that calibration. GFS seasons are not PR 41's held-out block.
- With 5 season blocks the bootstrap has at most 3,125 ordered resamples (126 distinct season multisets), so the intervals are narrow and the p values bottom out at 1/2001 (q 0.0006 for all four tests that pass). They describe season-to-season sampling in this set, not uncertainty about a future season, and fixes within a storm are not independent. A verdict of "yes" here means the sign is stable across these five seasons.
- GFS GUST is the model's own parameterisation, not ERA5's. Operational GFS is not a reforecast; model upgrades over 2021-26 are not controlled.
- Truth is ERA5 pipeline A, a proxy, not observation. 38 and 54 kt are nominal classes; only 71.7 kt is tied to the archive. Fixes with no match within 500 km are dropped (11,091 of 90,131), so the score covers storms the ERA5 proxy also tracked.
- A and B use only jet250, vadv500 and eady; the Hart terms, div300, SST, flux and tcwv were not pulled, so the full model is not tested in part 1 (it is in part 2, section below).
- Leads 0, 12, 36, 48 for A and B were not pulled (A/B scored at leads 0 and 24 only).

## Not independently checked
Every number above comes from `results/outlook.txt` and `results/outlook.json` produced by one run of `score_outlook.py` and has not been recomputed by an agent that did not produce it (the plan asks for a fresh agent to recompute counts, S/A/B RPSS and the G1-G4 point estimates; intervals and the pull are not part of that check). The environment terms were not compared against an independent implementation.
Reproduce: `python3 fit_models.py` (needs `results/` ERA5 tables), `GFS_WORK=... python3 fetch_env.py`, `python3 score_outlook.py`.

## Independent check (added after the first scoring)
A fresh Sonnet agent, writing its own code from the committed tables and the frozen coefficients, recomputed the counts (318,612 track rows; 90,131 fixes in basins; 11,091 dropped; 79,040 matched; lead-24 classes 6,220 / 6,745 / 4,366 / 836), G1 to G4 point estimates (+0.3147, -0.0286, +0.0157, +0.0146), the GFS-fed S RPSS at leads 0, 12, 24, 36 and 48, and the ERA5-fed S values by lead. All match to the digits shown. Not independently checked: bootstrap intervals, q values, per-basin, recalibration and reliability results, the G5 interval, and the GFS environment terms against a second implementation.

# Part 2: the FULL model on GFS, and do the environment gains survive? (`PREREGISTRATION_FULL.md`)

ERA5 pipeline A is a **proxy**: it is the training record, the truth, and (for the "ERA5-fed" side) the source of the environment terms. Nothing below is validated
against observations. Plan `PREREGISTRATION_FULL.md` (commit 0fae1f8, with dated Deviations written before the pull and before scoring). Numbers come from
`results/outlook_full.txt` / `.json`, one scoring run of `score_outlook_full.py`.

## What was done
1. `fit_models_full.py` (commit 8d2e1bf, pushed before any new GFS field was pulled or scored): FULL = S + Hart (B, VTL, VTU) + jet250, div300, vadv500, eady, sst, sstgrad, sst_t500, flux, tcwv, nosst, and FULL_NOSF (without sst, sstgrad, sst_t500, flux, nosst), each with a no-dp12 fallback; multinomial logistic, C = 1, ERA5 pipeline A seasons 2004-2020, four-class gcls, `intensity/model.py` Prep, no `young` / `logage`. dp12 models fitted on 100,386 rows (classes 53,103 / 30,428 / 14,457 / 2,398), fallbacks on 123,105. `results/outlook_models_full.json`; `outlook_models.json` was not modified.
2. `fetch_env_full.py`: HGT 925/850/700/600/500/400/300 mb, UGRD/VGRD 300 mb, TMP 500 mb, PWAT, TMP surface, ICEC (every file), LHTFL and SHTFL (lead 24 only), LAND once; leads 0 and 24 of the 2,122 cycles; byte-range method of `gfs_test/a_fetch.py`; `intensity/env.py` coarsen / derived / features at all 69,174 GFS fixes (27,943 lead 0; 41,231 lead 24). **44.19 GB pulled** over 4,244 files (mean 9.65 MB at lead 0, 11.17 MB at lead 24), above the plan's 33 GB because HGT 850, HGT 500 and TMP 500 had to be pulled again (part 1 kept only a per-fix table); the plan's 40 GB hard stop was raised to 46 GB in a dated Deviation before the pull (still under the 50 GB gate). **0 missing files in every season**, no NaN in B, VTL, VTU, div300, tcwv; sst / sst_t500 NaN at 20,532 fixes (no ocean cell within 500 km), flux NaN at all lead 0 fixes. Table `results/gfs_outlook_env_full.csv.gz` (2.7 MB) and `/mnt/project-files/gfs-test/outlook-env/` (README there lists the columns).
3. `score_outlook_full.py`: F1-F4F and the labelled secondaries, part 1 matching and truth (nearest pipeline A fix at the same valid time within 500 km, gcls), dp12/fallback rules and season-block bootstrap (5 GFS seasons, 2,000 draws, seed 20261013, 90 % intervals, BH over the 6). ERA5-fed side: every environment term and the state terms come from the matched pipeline A fix (`env_2004.csv.gz`, `fixes_2004.csv.gz` via `bomb_env/score.load`); no ERA5 was pulled. A plumbing run on a partial GFS table (804 rows) was made before the full table existed to catch crashes; it was not scored, its outputs were not used and it was not committed.

## Counts (from `results/outlook_full.txt`)
GFS track rows 318,612; atl/pac fixes at valid 00/12 UTC, leads 0-48: 90,131; no pipeline A fix within 500 km (dropped): 11,091; matched 79,040 (lead 0: 13,991; lead 24: 18,167). At both leads every matched fix has all GFS terms and all ERA5 terms finite (13,991 and 18,167), so the scored samples are the full matched sets, no extra filter dropped anything. Lead 24 scored sample by season: 2021: 3,488; 2022: 3,685; 2023: 3,849; 2024: 3,562; 2025: 3,583. Truth classes at lead 24 (below gale, gale, storm, HF): 6,220 / 6,745 / 4,366 / 836. dp12 model used for 15,322 of 18,167 GFS-fed fixes and 14,453 of 18,167 ERA5-fed fixes (the ERA5 fix may be `young`). No sst (no ocean within 500 km): 7.1 % of GFS-fed fixes against 4.2 % of the matched ERA5 fixes.

## Results (lead 24, pooled basins, RPSS vs basin-month climatology; point [90 % season-block interval], n = 18,167)
| Test | Estimate | q | Verdict (registered rule) |
|---|---|---|---|
| F1 FULL fed GFS, RPSS | +0.3509 [+0.3450, +0.3576] | 0.0005 | yes |
| F2 FULL minus S, GFS-fed | +0.0362 [+0.0332, +0.0396] | 0.0005 | yes |
| F3 FULL, GFS-fed minus ERA5-fed RPSS | -0.0352 [-0.0381, -0.0319] | 0.0005 | equivalent (inside +-0.05) |
| F4A retention, A | 0.765 [0.723, 0.805] | 0.0005 | survives (lower bound >= 0.5) |
| F4B retention, B | 0.934 [0.866, 1.018] | 0.0005 | survives |
| F4F retention, FULL | 0.845 [0.802, 0.890] | 0.0005 | survives |

Gains over S and the ratios (same fixes, same bootstrap draws): A +0.0157 [+0.0142, +0.0174] GFS-fed, +0.0205 [+0.0192, +0.0218] ERA5-fed; B +0.0146 [+0.0129, +0.0163] and +0.0157 [+0.0132, +0.0179]; FULL +0.0362 [+0.0332, +0.0396] and +0.0428 [+0.0402, +0.0458]. RPSS GFS-fed | ERA5-fed: S 0.3147 | 0.3433; A 0.3304 | 0.3638; B 0.3293 | 0.3590; FULL 0.3509 | 0.3861; FULL_NOSF 0.3391 | 0.3770. The S, A, B GFS-fed numbers reproduce part 1 (G1, G3, G4) exactly because the sample is the same 18,167 fixes.

Secondary (labelled; no multiplicity control):
- FULL_NOSF: gain over S +0.0244 [+0.0200, +0.0291] GFS-fed, +0.0337 [+0.0295, +0.0378] ERA5-fed, retention 0.725 [0.672, 0.776]. FULL minus FULL_NOSF: +0.0118 [+0.0099, +0.0138] GFS-fed, +0.0092 [+0.0071, +0.0111] ERA5-fed, i.e. the sst / sstgrad / sst_t500 / flux / nosst terms add at least as much on GFS as on ERA5 here, not less.
- Basins: Atlantic (n 8,607) FULL RPSS +0.3146, gain over S +0.0314 GFS-fed / +0.0387 ERA5-fed, retention 0.810 [0.732, 0.887]; Pacific (n 9,560) FULL RPSS +0.3848, gain +0.0406 / +0.0467, retention 0.871 [0.813, 0.935].
- Lead 0 (n 13,991; fallback model; GFS-fed flux imputed to the ERA5 training median, ERA5-fed uses its own flux): gain over S, GFS-fed | ERA5-fed: A +0.0385 | +0.0313 (ratio 1.229 [1.138, 1.311]); B +0.0429 | +0.0341 (1.257 [1.206, 1.323]); FULL +0.0528 | +0.0683 (0.774 [0.732, 0.820]); FULL_NOSF +0.0542 | +0.0544 (0.997 [0.948, 1.061]). The lead 0 ratios above 1 and the FULL / FULL_NOSF gap come in part from the imputed flux, so they are not comparable with lead 24.
- Reliability of FULL's hurricane-force probability at lead 24: slope 0.99 GFS-fed (mean forecast 0.0449, base rate 0.0460; top decile forecast 0.3546, observed 0.3511) and 1.09 ERA5-fed (mean 0.0418; top decile 0.3343 vs 0.3665). Mean class probabilities GFS-fed against observed: below gale 0.382 vs 0.342, gale 0.336 vs 0.371, storm 0.237 vs 0.240, HF 0.045 vs 0.046.

Reading against the expectation written in advance: FULL is lower fed GFS than fed ERA5 by 0.035 (inside the written range 0.03-0.05); A and B keep 0.77 and 0.93 of their ERA5-fed gain (written range 0.7-1.0); FULL keeps 0.85. The expectation that the SST and flux terms would lose most on GFS is **not** supported at lead 24: dropping them lowers the GFS-fed RPSS by 0.012 and the ERA5-fed by 0.009, and FULL_NOSF keeps a smaller share of its ERA5-fed gain (0.73) than FULL (0.85); the Hart, div300, tcwv part loses more of its ERA5-fed gain than the SST / flux part. This is a descriptive comparison of two models, not a test of a single term.

## Deviations
Dated notes in `PREREGISTRATION_FULL.md`: (1) the re-pull of HGT 850, HGT 500, TMP 500 (not in the plan's 33 GB) and the raised stop (46 GB), ICEC per file, SST masking before coarsening, flux sign, track heading, zero stand-ins for the part 1 terms' inputs; (2) scoring rules the plan left open (sample, model choice per fix, ratio definition, p values, verdicts, secondaries). The registered tests, seeds, models and data are otherwise as written.

## Limits
- Five GFS seasons, all inside the block that calibrated the 71.7 kt gust cut. With 5 season blocks the bootstrap has only 126 distinct resamples, so intervals are narrow and p floors at 1/2001 (q 0.0005): "yes", "survives" and "equivalent" mean the sign and size are stable across these five seasons, not that a future season is covered. Fixes within a storm are not independent.
- GFS skin temperature (TMP surface) is not ERA5 SST; the GFS flux is an 18-24 h average, ERA5's is the hourly mean at the fix time; GFS GUST is the model's own. These are the differences that the FULL_NOSF comparison was meant to expose; it did not show a larger loss for them at lead 24.
- Truth is ERA5 pipeline A (a proxy), not observation; 38 and 54 kt are nominal; fixes without a pipeline A match within 500 km are dropped (11,091 of 90,131), so the score covers storms the ERA5 proxy also tracked.
- The ERA5-fed side feeds the models the environment of the matched ERA5 fix (a fix up to 500 km away, median 82 km), so it measures "same storm, ERA5 analysis" and not a different track.
- F4 ratios are ratios of small gains (+0.015 to +0.04 RPSS); the bootstrap interval reflects season sampling only. The F4B upper bound exceeds 1.
- Descriptive means of the terms on GFS-fed vs matched ERA5 fixes are in `outlook_full.txt` (e.g. sstgrad 1.387 vs 1.547, vadv500 0.357 vs 0.362); they are not a validation of the GFS terms.

## Not independently checked
Every number in this section comes from `results/outlook_full.txt` / `.json`, one run of `score_outlook_full.py`; none has been recomputed by an agent that did not produce it (the plan asks a fresh agent to recompute counts, the RPSS of FULL, FULL_NOSF, A, B, S in both feeds, and the F1-F4 point estimates from the committed tables; intervals and the pull are not part of that check). The GFS environment terms were not compared against an independent implementation. S, A, B GFS-fed numbers match part 1's independently checked values. Only one scoring run was made.
Reproduce: `python3 fit_models_full.py`; `GFS_WORK=... python3 fetch_env_full.py` then `--assemble results/gfs_outlook_env_full.csv.gz`; `python3 score_outlook_full.py`.

## Independent check of part 2 (added after the scoring)
A fresh Sonnet agent, writing its own code from the committed tables and the frozen coefficients, recomputed the counts (318,612 track rows; 90,131 fixes in basins; 11,091 dropped; 79,040 matched; lead 24: 18,167 fixes, classes 6,220 / 6,745 / 4,366 / 836, by season 3,488 / 3,685 / 3,849 / 3,562 / 3,583; dp12 model used for 15,322 GFS-fed and 14,453 ERA5-fed fixes; no-SST share 7.1 % vs 4.2 %), the lead-24 RPSS of S, A, B, FULL and FULL_NOSF fed GFS and fed ERA5, the gains over S, the retention ratios (0.765, 0.934, 0.845, 0.725) and F3 (-0.0352). All match to the digits shown. Not independently checked: bootstrap intervals, q values and the "survives" and "equivalent" verdicts, the secondaries (per basin, lead 0 and the flux imputation, reliability), the pull, and the GFS environment terms against a second implementation.
Note: jet250, vadv500 and eady for the scored fixes are in the part 1 table `results/gfs_outlook_env.csv.gz`, not in `gfs_outlook_env_full.csv.gz`; join the two on cycle, lead and track.
