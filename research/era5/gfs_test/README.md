# gfs_test: Stage A data side (phase 1)

Pre-registration: `PREREGISTRATION.md`. Phase 1 builds, freezes and pulls; nothing here is scored.

- `a_fetch.py` GFS pull (PRMSL, GUST:surface by .idx byte ranges; LAND once), processed on the fly into per-cycle lows and per-file basin gust. Work area `/home/claude/gfs_work` (not committed).
- `a_tracks.py` pipeline A detector (`hf_history/extract.py` `lows_at`, same constants, on the GFS 721x1440 0.25 degree grid, which is the grid pipeline A uses) and linker (`track.link`, MAXD 900 km per 6 h) applied per cycle over f000-f048; writes `results/gfs_tracks_<season>.csv.gz` (tracks of at least 4 fixes, `track.MIN_LEN`, as pipeline A) and `results/gfs_basin_gust.csv.gz`.
- `a_fit_model.py` fits and freezes `results/gfs_model.json` (ERA5 pipeline A proxy fixes, seasons 2004-2020, target hf24, never scored). Lead 0 uses `model_nodp12`; lead >= 12 h uses `model_dp12`.
- `tests/test_a_tracks.py` synthetic Gaussian low plus planted gust patch (no network).

Training rows (fixes_2004.csv.gz joined to env_2004.csv.gz as `intensity/model.py:load`, which gives 159,430 rows for 2004-2025, the count in the model's skill.txt; 123,105 rows for 2004-2020):

| model | rows | Atlantic | Pacific | hf24 = 1 |
|---|---|---|---|---|
| model_dp12 (dp12 present) | 100,386 | 47,461 | 52,925 | 3,405 (atl 1,888, pac 1,517) |
| model_nodp12 (all rows) | 123,105 | 57,878 | 65,227 | 4,085 (atl 2,270, pac 1,815) |

22,719 rows dropped from the dp12 fit because dp12 is missing (tracks younger than 12 h).


## Results (written 2026-10-08; pipeline A is the ERA5 proxy; archive truth where stated)

Pre-registration `PREREGISTRATION.md` (committed 95035fe before any pull; frozen model 27082e9 before any score). 14 registered tests; 11 meet their registered criterion, 10 have BH q < 0.05 (A2.2 is an equivalence criterion judged on the point difference). Five GFS seasons 2021-22 to 2025-26 (the gust-calibration seasons, so not independent of the calibration); GEFS v12 reforecast 2004-2018 for the pattern product. Pulled: GFS 31.2 GB, GEFS 1.3 GB, ERA5 7.7 GB (40.3 GB, under the 50 GB gate).

| Test | Atlantic | Pacific | Criterion |
|---|---|---|---|
| A1.1 basin BSS at +24 h (GFS max gust vs month base rate) | +0.241 [0.211, 0.278] | +0.275 [0.237, 0.316] | pass |
| A1.2 AUC falls with lead (0 vs 48) | 0.882/0.881/0.877, diff +0.005 [-0.004, 0.016] | 0.900/0.904/0.899, diff +0.001 [-0.017, 0.013] | fail: no measurable decay |
| A1.3 HSS at +24 h, threshold set on 2021 only | 0.502 [0.474, 0.543] | 0.489 [0.452, 0.522] | pass (POD 0.69 / 0.51, FAR 0.50 / 0.41) |
| A2.1 track onset-probability BSS vs climatology at +24 h | +0.327 [0.266, 0.376] | +0.316 [0.264, 0.359] | pass |
| A2.2 GFS-fed minus ERA5-fed BSS at +24 h | -0.041 [-0.066, -0.024] | -0.061 [-0.077, -0.036] | Atlantic pass on the point rule, Pacific fail; no interval sits inside +-0.05 |
| B.1 GEFS vs ERA5 pattern index (Spearman) | 0.960 [0.948, 0.968] | 0.942 [0.926, 0.952] | pass |
| B.2 week-ahead deviance skill (GEFS index) | 9.14% [5.39, 12.87] | 6.78% [2.63, 10.76] | pass, but in-sample for the weights (deviation 7); 2015-18 only 4.5% / 0.7% |

Other: A2.3 skill does not decay across leads 0-48 (BSS 0.29-0.34 Atlantic, 0.30-0.35 Pacific). A2.4 recalibration changes lead-24 BSS by -0.001 and +0.005 (no gain). Against the archive (secondary truth, any HF fix within 600 km in the next 24 h, base rate 4.9% / 3.6%): BSS +0.17 / +0.25 against the hf24 climatology (indicative only), HSS 0.43 / 0.49. 11,091 of 90,131 scored fixes have no ERA5 fix within 500 km and count as 0. GFS minus ERA5 g800 averages -0.5 kt; 13.7% of lead >= 12 fixes used the no-dp12 fallback.

What this shows and does not show: the state-only onset-probability model, trained only on ERA5 2004-2020 and never refitted, keeps most of its skill when fed GFS forecast tracks out to +48 h (loses 0.04-0.06 BSS against feeding it ERA5 for the same storms). The pattern index computed from GEFS forecasts tracks the ERA5 one closely, but the registered skill test was partly in-sample, and the clean 2015-18 subset is small and its Pacific interval spans zero. Not tested: the Hart phase-space and environment terms of the full PR 12 model, the wind-structure checklist, GFS fields at sub-6 h, or any season after 2025-26.

Verification: a fresh Sonnet agent recomputed the A1/A2/B point estimates, counts, the 2021 thresholds and the frozen-model hash from the committed files; all matched (ERA5-fed A2.2 value depends on which model is used when ERA5 dp12 is missing: +0.377 by the registered rule, +0.381 by the alternative). Not independently checked: bootstrap intervals and p values, the A1 power simulation, B.2 power and minimum detectable effect, the reliability tables.
