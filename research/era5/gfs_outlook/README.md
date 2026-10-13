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

Reading: S keeps most of its ERA5-fed skill when fed GFS (RPSS 0.31 against 0.34, a loss of 0.03 that sits inside the registered +-0.05 margin), and the upper-level terms add a small gain over the storm's own state (about +0.015 RPSS at lead 24, +0.04 at lead 0). The A/B gains were not computed ERA5-fed on this sample, so whether they shrank against ERA5 is not tested here (the expectation written in advance cannot be judged). The expected flat-or-falling skill-versus-lead curve is not what was found: RPSS rises with lead. A likely contributor (not tested) is that lead 0 uses the fallback model with imputed speed, and that fixes at later leads are a different, longer-lived subset of storms with a different climatology reference; the comparison across leads is between different samples.

## Deviations
See the dated notes in `PREREGISTRATION.md`: scoring rules the plan left open (model choice per fix, the ERA5-fed rule, G3/G4 sample, p-value and verdict definitions, recalibration), and the pull implementation (eccodes installed with pip; `env.py` imported with its zarr readers stubbed). No deviation from the registered tests, data or seeds otherwise.

## Limits
- Five GFS seasons, all inside the block that calibrated the 71.7 kt gust cut; not independent of that calibration. GFS seasons are not PR 41's held-out block.
- With 5 season blocks the bootstrap has at most 3,125 ordered resamples (126 distinct season multisets), so the intervals are narrow and the p values bottom out at 1/2001 (q 0.0006 for all four tests that pass). They describe season-to-season sampling in this set, not uncertainty about a future season, and fixes within a storm are not independent. A verdict of "yes" here means the sign is stable across these five seasons.
- GFS GUST is the model's own parameterisation, not ERA5's. Operational GFS is not a reforecast; model upgrades over 2021-26 are not controlled.
- Truth is ERA5 pipeline A, a proxy, not observation. 38 and 54 kt are nominal classes; only 71.7 kt is tied to the archive. Fixes with no match within 500 km are dropped (11,091 of 90,131), so the score covers storms the ERA5 proxy also tracked.
- A and B use only jet250, vadv500 and eady; the Hart terms, div300, SST, flux and tcwv were not pulled, so the full model is not tested.
- Leads 0, 12, 36, 48 for A and B were not pulled (A/B scored at leads 0 and 24 only).

## Not independently checked
Every number above comes from `results/outlook.txt` and `results/outlook.json` produced by one run of `score_outlook.py` and has not been recomputed by an agent that did not produce it (the plan asks for a fresh agent to recompute counts, S/A/B RPSS and the G1-G4 point estimates; intervals and the pull are not part of that check). The environment terms were not compared against an independent implementation.
Reproduce: `python3 fit_models.py` (needs `results/` ERA5 tables), `GFS_WORK=... python3 fetch_env.py`, `python3 score_outlook.py`.
