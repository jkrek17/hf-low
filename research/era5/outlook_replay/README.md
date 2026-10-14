# Calibrated outlook shading from GFS gust: how often OPC agreed

**GFS forecast gust (a model forecast, PROXY for sustained wind) against OPC archive fixes; 5 seasons (Oct-Apr 2021-22 to 2025-26), 00 UTC cycles, days 1 to 4.** "OPC agreed" = the archive has a Category HF (or S or HF) low fix within 324 nm of the point during the forecast day; the call is an analyst judgement, so this is agreement, not verification. Pre-registered in `PREREGISTRATION.md` (deviations at its end).

Files: `grib2_t53.py` (GRIB2 reader), `fetch_dmax.py` (pull; the 129 MB product is in `/mnt/project-files/gfs-daily-max-gust/`), `tabulate.py` -> `results/point_days.csv` (counts per cycle, day, basin, gust band), `analyse.py` -> `results/calibration.txt`, `bands_pooled.csv`, `loso_brier.csv`, `select_cases.py` -> `results/replay_cases.csv`, `build_payload.py` -> `results/replay_payload.json` (Replay page input).

## Result (both basins; share of 1-degree ocean point-days with an OPC HF low within 324 nm)
Max GFS gust within 60 / 120 / 270 / 270 nm (days 1 to 4): below 40 kt 0.2 / 0.1 / 0.0 / 0.0%; 56-64 kt 11.9 / 7.6 / 2.6 / 3.0%; 64-72 kt 29.0 / 22.1 / 12.4 / 12.0%; 72-80 kt 49.5 / 44.6 / 31.3 / 27.8%; 80 kt and above 57.6 / 57.0 / 46.7 / 38.9%. Each band is monotone in every day and basin. Held-out Brier skill of the banded table against the held-out season's own base rate: HF 0.141 / 0.169 / 0.184 / 0.143 (both basins; seasons range 0.123 to 0.196). Full table and per-season range in `results/calibration.txt`.

Limits: five seasons, one model with version changes; point-days within a cycle are not independent, the sample is cycles and seasons; 6-hourly gust samples can miss a peak; OPC's HF call is sometimes made as a service, and a model can also under-forecast; 324 nm truth radius is the one used in all earlier GFS tables, wider than the 60 and 120 nm shading radii; nothing beyond day 4 or outside October to April.

Check: a fresh agent recomputed four cycles with great-circle distances (the counts here use a local grid approximation for the radius): same evaluation points, band counts within a median 1.4% (mean 2.8%) per row, HF and storm totals per day identical, and the pooled shares, the held-out Brier skill (0.169, seasons 0.149 to 0.183, day 2, both basins) and the monotone band shares agree. The day window includes the fix at +6 h (lower end inclusive).
