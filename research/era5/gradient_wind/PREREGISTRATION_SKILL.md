# Pre-registration: does knowing the pressure gradient add skill to the GFS wind forecast?

Written and committed before any model was fitted to these data. Follows `PREREGISTRATION.md` and PR 282 (Jason, 2026-10-28: "see if there is any added skill by knowing the pressure gradient"). Reuses cached tables only, no download.

## What has been seen

Sizing only: the columns and row counts of `/mnt/project-files/gulf-stream/tables/boxes_A.csv.gz` and `boxes_B.csv.gz`. From earlier work on the same tables (PR 269, 271): GFS 24 h wind against ASCAT was compared by stability, not with a gradient term. No fit of ASCAT on GFS wind plus gradient has been run.

## Data

- Outcome: ASCAT 1-degree box-pass mean wind (kt), equivalent-neutral, rain-free, land-masked, western Atlantic 29-41N, 80-55W, October to April (`v_ascat`).
- Baseline predictor: GFS 24 h forecast 10 m wind speed valid at the 00/12 UTC nearest the pass (`v_gfs`, table B).
- Gradient: ERA5 geostrophic wind from MSLP at the valid time, in kt (`vg`, table A, built by `gulf_stream_region/build_boxes.py::vg_field`). Joined to table B on day, segment, box, satellite. Winters 2021-22 to 2024-25 (4 winters; table A stops at 2024).
- Stability for the secondary model: GFS skin minus 2 m temperature, clipped -3 to 12 C (`dT`, table B).
- ERA5 is a proxy. `vg` is the **analysis** gradient at the valid time, not GFS's own 24 h forecast gradient: it tells what knowing the true gradient adds, an upper bound on what a forecaster who correctly reads the same chart adds. GFS's own forecast pressure is not cached; a variant with it needs a new pull (not run, estimated at about 2 GB, not probed).
- ASCAT saturates near 45 to 50 kt, so everything here describes winds below about 45 kt; nothing applies at storm or hurricane force. GFS wind is the real 10 m wind and ASCAT is equivalent-neutral, so a fitted linear scale absorbs part of the stability difference; the stability model M2 is secondary.

## Models (ordinary least squares, ASCAT wind in kt)

- M0: 1 + `v_gfs` (the calibrated GFS baseline).
- M1: M0 + `vg`.
- M1q (secondary): M1 + `v_gfs`^2 + `vg`^2 + `v_gfs` x `vg`.
- M2 (secondary): M1 + `dT`.
- M3 (secondary): 1 + `vg` alone.

## Scoring

Leave-one-winter-out: fit on three winters, score on the fourth, pool the four. Skill = 1 - MSE(model) / MSE(M0), reported with RMSE in kt. Interval: bootstrap over **days** (box-passes on one day share a weather situation), 2000 draws. With only 4 winters the interval understates the uncertainty; the rule below also needs all four winters to agree.

## Planned tests (two primary, Benjamini-Hochberg across them)

- T1: M1 against M0, all box-passes. 
- T2: M1 against M0, box-passes with GFS wind of 25 kt or more (where forecasters care).
- A meaningful gain is **2% of MSE or more** (about 1% in RMSE). Verdicts: "yes" if the corrected interval excludes 0, the gain is at least 2%, and all four winters gain; "no (well powered)" if the upper end of the interval is below 2%; otherwise "cannot tell".

## Changes after seeing results

Logged below as post hoc, with the date.
