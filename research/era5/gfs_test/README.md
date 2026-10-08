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
