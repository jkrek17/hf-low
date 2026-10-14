# How wide should an outlook neighbourhood be? (position error only)

**ERA5 PROXY gust, pipeline A tracks; not direct observation.** Pre-registered in `PREREGISTRATION.md` (deviations listed at its end). All distances in nautical miles.

## What was done
1. `step1_position_error.py`: GFS-low versus ERA5-low position difference by lead (5 GFS seasons 2021-22 to 2025-26), `results/position_error.txt`.
2. `step2_coverage.py`: take the ERA5 gust field (22 seasons, Oct-Apr, 00/12 UTC, ocean only), move it by a random Step 1 position error, shade everything within R of the moved hazard, and count how much of the real hazard falls inside. Needs `blosc_decode.js` and the cache in `/mnt/project-files/era5-cache/` (not in git).
3. `step3_summarise.py`: season-block bootstrap (22 seasons, 2,000 resamples), `results/radius_results.csv` and `.txt`.

## Result (HF = gust >= 71.7 kt, share of real HF area inside the shading)
Decision rule (smallest R with coverage >= 0.80 in both basins): 0 h 30 nm, 24 h 60 nm, 48 h 120 nm, 72 h 270 nm, 96 h 270 nm (the lower end of the 95% interval gives the same answer). Shaded area is about 5x (60 nm), 11-12x (120 nm), 33-42x (270 nm) the forecast HF area.

This is an **upper bound**: only the low-position part of the error is simulated. Intensity errors, missed storms and invented storms are not in it, so a real GFS map covers less at the same radius.
