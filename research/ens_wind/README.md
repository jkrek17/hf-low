# Ensemble 10 m wind chances: is "6 of 31 members" really about 2 in 10?

Asked by Jason 2026-10-28 ("10 m winds from the GFS ... do we calibrate ensemble probabilities?"). Plan, choices and every deviation are in `PREREGISTRATION.md`, committed before any field was pulled. Sustained 10 m wind from the models, in knots; no ERA5 is involved. Pulled: GEFS 31 members, six winters (2020-21 to 2025-26), 00 UTC runs every 5th day from 15 Oct to 15 Apr, leads 24, 48, 72, 96 h (222 runs); ECMWF ENS and AIFS ENS (51 members each) for 2025-26 at 48 and 96 h. About 58 GB of byte-range reads (logged in `bytes.log` during the run), plus about 6 GB of ASCAT daily files.

## What to take from it
- **48 kt: the GEFS member fraction is already about right.** Where ASCAT saw the cell, Atlantic day 1: no member at 48 kt, it happened 0.2 % of the time (about 2 in 1,000); up to 1 in 10 members, 11 %; 1 to 3 in 10, 35 %; 3 to 6 in 10, 64 %; over 6 in 10, 91 %. Pacific is the same shape, slightly lower. The calibration slope is 0.92 to 1.11 in all eight cells and none is flagged overconfident. A fitted correction adds a little skill in the Atlantic (Brier skill +0.02 to +0.035, intervals above zero) and about nothing in the Pacific (+0.000 to +0.017, intervals cross zero).
- **64 kt: members under-call it when only a few show it, and over-call it a bit at the top.** Truth here is "OPC had a hurricane-force low within 324 nm", so these are neighbourhood chances, not chances at a point. Up to 1 in 10 members reaching 64 kt: it still happened 13 to 18 % of the time (days 1 to 3, both basins); 1 to 3 in 10: 33 to 54 %; over 6 in 10: 69 to 84 % (Atlantic day 3: 84 %). The fitted slope is below 1 in 5 of 8 cells, but after the multiplicity correction only Atlantic day 3 stays flagged (the other overconfident-looking cells could be noise, 6 winters). A fitted correction helps in 7 of 8 cells (Brier skill +0.02 to +0.09; Pacific 48 h crosses zero). The correction table is `results/gefs_forecaster_tables.txt`; at high member counts it is thin, read neighbouring entries.
- **ECMWF ENS is about as good as GEFS; AIFS ENS is not usable raw.** On 2025-26 at 48 and 96 h, ENS minus GEFS Brier skill is between -0.07 and +0.005 and every interval crosses zero. AIFS ENS is 0.17 to 0.33 worse, with intervals below zero. It never put a member at 64 kt in any cell of the winter, and at 48 kt it says "none" for almost everything, yet 1.2 % of its "none" cells had 48 kt (GEFS and ENS: 0.1 to 0.2 %). One winter, wide intervals.
- **Check on leakage.** Taking the truth from a date 23 days away, calibrated skill is 0.000 (-0.004 at worst) in all six cells that ran.

## Numbers
`results/gefs_primary.txt` (all 16 registered cells: cases, events, rate by bin, slope with interval, Brier skill raw and calibrated, BH decision), `results/compare_2025-26.txt`, `results/gefs_forecaster_tables.txt`. `results/gefs_primary.json` has the per-season pieces.

## Not usable: 48 kt against OPC
Reported in `results/gefs_primary.txt` (secondary block) but not interpreted. OPC's archive holds storm-force only on lows that earned a category, so even with every member at 48 kt it records an event only 28 to 37 % of the time and raw skill is strongly negative. Calibrating to that would calibrate to the archive's completeness, not to the wind.

## Limits
- ASCAT sees about half of any area, saturates near 50 kt, and compares stress-equivalent wind with the models' actual 10 m wind (no COARE here, see deviation 1). 48 kt rows are "cells ASCAT observed". 64 kt cannot be scored against ASCAT.
- Forecast cell = maximum of the ocean points in a 1 degree cell; ASCAT's cell maximum covers at least half of the ocean sub-cells. This favours the model seeing more.
- Six winters, 37 run dates per winter, every 5th day. Intervals resample whole winters, so they are wide for the 64 kt rows. Climatology for Brier skill is one number for the whole basin, so skill includes knowing where storms go.
- GEFS model version changes inside 2020-26 were not checked.
- OPC "agreed" is an analyst call. Numbers from the marine wind project's GFE blend (about 15 % of ASCAT 48 kt+ peaks at about 4 % false alarm; nothing skilful at 64 kt) are reported there, not rechecked here, and are not the same quantity as these cell chances.

## Checked and not checked
A fresh Sonnet agent recomputed from the raw member maps, ASCAT files and OPC csv, with its own code and the written spec only: cases, events, bin rates and raw Brier skill for GEFS Atlantic 48 kt day 1, Pacific 64 kt day 2, Atlantic 64 kt day 3, and the GEFS/ENS/AIFS ENS Atlantic 48 kt 48 h comparison. All matched. **Not rechecked:** slopes, bootstrap intervals, BH decisions, isotonic tables, calibration gains, the placebo run, and the other cells.

## Reproduce
`ens_fetch.py` (pull; resumable), `ascat_truth.py` (ASCAT cells), `build_counts.py` (count tables), `run_primary.py`, `compare.py`, `stats.py`, `test_stats.py` (planted-effect test, no network). Set `ENS_WORK` to a scratch directory (default `work/ens_wind`, ignored by git). Pulled fields (about 200 MB of cell maxima) are not committed.
