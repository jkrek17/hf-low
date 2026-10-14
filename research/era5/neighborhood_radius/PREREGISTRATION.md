# Pre-registration: how wide should an outlook neighbourhood be? (position error only)

Written 2026-10-14 before any score by radius or any position error was computed. Jason's ask (thread "What we learned summary page", "lets think about what the radius should be", then "do the radius test"). The question: if an outlook map shades the area within R nautical miles of the model's predicted hazard, what R is needed so that the shading covers where the hazard really occurs, given how far GFS puts lows from where ERA5 has them? **ERA5 PROXY throughout; no direct observation. Pipeline A for tracks.**

## What this can and cannot show
- It measures the **displacement (position) part of forecast error only**. The "forecast" is the ERA5 gust field itself, moved by a position error drawn from GFS-versus-ERA5 low-position differences. Intensity error, shape error, missed storms and storms that do not exist in the forecast are NOT in it. So coverage here is an **upper bound** on what a real GFS-driven map would achieve at that radius.
- Truth is the ERA5 gust proxy, not OPC (OPC records low positions, not wind areas) and not measured wind. The ASCAT check could later replace it.
- No hypothesis test. This is a design measurement, reported as curves with season-block intervals. Nothing below is a multiplicity family.

## Step 1: position error by lead (existing tables, no pull)
- GFS fixes: `research/era5/gfs_raw_gust_96h/results/raw_gust_fixes_96h.csv.gz` (5 seasons 2021-22 to 2025-26, leads 0 to 96 h, every 12 h; `lat`, `lon`, `msl`, `valid_time`, `G`).
- ERA5 fixes at 00/12 UTC: `research/era5/intensity/results/fixes_2004.csv.gz` (`time`, `lat`, `lon`, `msl`).
- Match each GFS fix to the nearest ERA5 fix at the same valid time within **1,200 km** (great-circle). No pressure test (it would shrink the errors it measures). Report the matched share by lead, so unmatched cases are visible (they have larger errors or no ERA5 low).
- Primary group: GFS fixes with `G` >= 64 kt (near-hazard storms). Secondary: all fixes. Report n, median, mean, 75th and 90th percentile of distance in nm, by lead, per basin, with the signed (dlat, dlon) means to show any systematic shift.
- The matched (dlat, dlon) vectors of the **primary group** are the displacement sample for Step 2 (pooled over both basins if a basin has fewer than 100 matches at a lead, stated).

## Step 2: shading coverage by radius (existing ERA5 gust cache)
- Truth fields: `/mnt/project-files/era5-cache/gust_nh_12h/` (0.25 degree, 00/12 UTC, uint8 in 0.5 kt steps). **October to April only**, seasons 2004-05 to 2025-26 (22 seasons, about 9,300 times). Domains: North Atlantic 100W to 10E, North Pacific 120E to 120W, both 20N to 75N; **ocean cells only** (`lsm_nh_0p25.npy`, ocean where < 0.5). Pacific and Atlantic scored separately.
- Hazard classes (gust index): **gale >= 38 kt, storm >= 54 kt, hurricane force >= 71.7 kt** (the cache's 0.5 kt steps put the HF cut at stored value 144, i.e. 72.0 kt; stated). Only HF is tied to OPC; gale and storm are nominal.
- Forecast field F = truth field T shifted by one displacement vector drawn at random (seeded) from the Step 1 sample for that lead. Leads 0, 24, 48, 72, 96 h. **Two independent draws per time, per lead** (variance reduction only). A shift moves the whole basin field together; fill outside is empty.
- For radii R = **30, 60, 120, 270, 432 nm** (270 and 432 are the stored search radii of the GFS tables), shading S_R = all ocean cells within R of any F cell. Distances use grid-latitude bands of 5 degrees with the band-centre longitude spacing (stated approximation, a few per cent).
- Per basin, class, lead and R, over all times with at least one T cell: 
  - **Coverage** = (T cells inside S_R) / (T cells). The share of the real hazard the shading would have caught.
  - **Shading cost** = (cells in S_R) / (F cells). How much larger the shaded area is than the forecast hazard area it was built from.
  - **Neighbourhood precision** = (F cells with a T cell within R) / (F cells). How much of the forecast hazard has real hazard nearby.
- Pooled as ratios of sums over seasons; **90% intervals by resampling the 22 seasons** (2,000 resamples, seeded). Fixes and times within a season are not independent, the season is the unit.

## Decision rule (fixed now)
For each lead, the **recommended radius** is the smallest R on the grid with HF **coverage >= 0.80 in both basins** (point estimates). 0.80 is a design target chosen before looking, not an optimum. If none of the grid reaches 0.80 at a lead, report "above 432 nm". Gale and storm curves are reported next to it with no separate rule. The cost curve is always shown beside the radius.

## Not claimed
No claim about skill of GFS gust, no calibration of probabilities, no truth beyond the proxy. A recommended radius says how wide to draw a footprint to absorb **position error**, not to absorb intensity error or missed storms. Real GFS shading would need to be checked separately (scored GFS gust fields beyond 96 h do not exist here).

## Deviations (post hoc)
(none yet)
