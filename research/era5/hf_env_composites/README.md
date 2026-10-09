# Storm-centred composites of every HF low: Z500, 250/500 hPa wind, divergence (ERA5 proxy, pipeline A)

Plan: [PREREGISTRATION.md](PREREGISTRATION.md) (committed `a4deb71` before any field was read). **ERA5 proxy**: an
"HF low" is a pipeline A low whose 800 km ocean gust index reaches 71.7 kt at a 00/12 UTC fix; not the archive.
Tropical cyclones are in. Descriptive study: averages with no comparison group, so nothing here says what
*distinguishes* HF lows (the Composites page does that against storm-force-only lows).

Population: 1,223 HF lows (664 Atlantic, 559 Pacific), onsets October to April, seasons 2004-05 to 2022-23
(WeatherBench2 ends 2023-01-09; the 2023 onward 0.25 degree gap is not covered). Times: onset, peak, and the same
track 12, 24 and 48 h before onset. Fields: WeatherBench2 ERA5 at 1.5 degrees. Anomaly = ERA5 minus the WeatherBench2
1990-2019 climatology (6-hourly by day of year).

## Answer in plain words (what the average HF low looks like)
1. **It forms on the equatorward side of a deep upper low, under the jet, and the two stack by onset.** Two days before
   onset the mean 500 hPa low sits about 2,000 km to the left of the storm's track (Atlantic: 5,175 m at x = +800,
   y = +2,000 km; Pacific 5,223 m at +1,000, +2,100) with the surface storm on its southern flank. By onset the
   minimum is over the storm in the Atlantic (5,125 m at onset, 5,092 m at peak) and in the Pacific is still about 1,400 km left at onset (5,193 m) and over the centre at peak (5,188 m).
   The Z500 anomaly at the composite minimum deepens from -100 m (48 h before) to -261 m (onset) and -280 m (peak) in
   the Atlantic, and from -74 to -197 and -218 m in the Pacific. A positive Z500 anomaly (a ridge, tens of metres) sits to
   the right of the storm's track throughout.
2. **A strong jet runs over the storm, and by onset it is behind and to the right of it.** The mean 250 hPa wind over
   the storm is 39 to 45 m/s (Atlantic) and 53 to 58 m/s (Pacific) 48 to 24 h before onset, then drops to 26 and 35 m/s
   at onset and 23 and 30 m/s at peak. The jet maximum in the mean field moves from about 600 to 900 km ahead of the
   storm (48 and 24 h before, Atlantic; 900 km ahead at 48 h in the Pacific) to about 800 to 900 km behind and 200 to
   600 km to the right from 12 h before onset (Atlantic), and 1,000 to 1,100 km behind from 24 h before (Pacific). Part of that shift is the
   storm's own circulation, which is not removed here: its cyclonic wind adds to the jet on its south side and gives a
   wind minimum at the centre, so these maps cannot say how much is jet and how much is storm.
3. **Upper-level divergence sits just ahead of the surface low.** The 250 hPa divergence maximum is +1.8 to +2.5
   (Atlantic) and +1.9 to +3.2 (Pacific) x 1e-5 /s, 200 to 500 km ahead of the centre at every time, with convergence
   behind it. At 500 hPa the divergence is an order of magnitude smaller and noisy at 1.5 degrees; do not read detail
   into it.
4. **What this cannot show.** An average of many storms smears jet streaks and troughs, so entrance and exit quadrants,
   trough tilt and trough position cannot be read from it. The per-storm check (S3) says why: the 250 hPa maximum
   within 2,000 km is ahead of the storm in 54% (Atlantic, 48 h before) falling to 33% at peak, with a median distance of
   about 1,340 to 1,550 km, so individual storms differ a lot from the mean. The pre-registered detection test of jet
   quadrant and trough features is a separate piece of work.

## Numbers
n = storms with a valid centre value (rotated frame needs a heading, so it is a little smaller than the north-up n).
Rotated frame, `results/scalars.csv` (S1 to S4 as defined in the plan); all S2/S4 positions are x ahead, y left, km.

| | Atl 48 h | Atl 24 h | Atl 12 h | Atl onset | Atl peak | Pac 48 h | Pac 24 h | Pac 12 h | Pac onset | Pac peak |
|---|---|---|---|---|---|---|---|---|---|---|
| n (rotated) | 151 | 370 | 514 | 645 | 650 | 92 | 325 | 469 | 551 | 555 |
| S1 mean 250 hPa speed at centre (m/s) | 38.9 | 45.3 | 39.8 | 26.1 | 22.5 | 53.2 | 57.7 | 49.2 | 34.6 | 30.2 |
| S1 mean 500 hPa speed at centre (m/s) | 25.9 | 29.7 | 25.9 | 17.6 | 15.5 | 30.3 | 33.0 | 28.4 | 20.5 | 18.3 |
| S2 max of mean 250 hPa speed (m/s) | 44.8 | 47.7 | 48.6 | 46.3 | 45.5 | 56.0 | 62.4 | 63.3 | 59.7 | 59.1 |
| S2 its position (x, y) | +900, +200 | +600, +100 | -900, -200 | -900, -400 | -800, -600 | +900, +100 | -1000, 0 | -1100, 0 | -1000, -400 | -1100, -400 |
| S4 min of mean Z500 anomaly (m) | -100 | -149 | -197 | -261 | -280 | -74 | -100 | -138 | -197 | -218 |
| S4 its position (x, y) | -400, +400 | -400, +400 | -300, +200 | 0, 0 | 0, 0 | +1600, +1600 | -500, +500 | -400, +200 | -100, 0 | 0, 0 |
| S4 max of mean 250 hPa divergence (1e-5/s) | 1.8 | 2.0 | 2.5 | 2.3 | 2.1 | 1.9 | 2.4 | 3.0 | 3.2 | 2.9 |
| S3 share of storms with the 250 hPa max ahead (90% interval) | .54 (.48-.59) | .45 (.40-.51) | .34 (.30-.38) | .34 (.31-.37) | .33 (.30-.35) | .59 (.50-.66) | .42 (.38-.45) | .31 (.27-.35) | .26 (.22-.29) | .24 (.21-.27) |

Post hoc (not in the plan, labelled so): the position of the raw mean Z500 minimum in the rotated frame, quoted in
point 1, was read from the same composite grids after the figures were seen.
Season-block bootstrap (2,000 resamples of seasons, seed 20261009), pixelwise BH q < 0.05: nearly every pixel of the
speed anomalies passes (the mean wind in HF-low neighbourhoods is higher than the climatology almost everywhere), so
stippling on the speed figures is uninformative; the Z500 and divergence stippling is selective. Pass counts are in
`results/summary.txt`. Pixels are correlated, so stippling is descriptive.

## Figures (`results/figs/`; rotated frame is primary, north-up for orientation)
Names: `{rot,north}_{atl,pac}_{z500,wind250,wind500,div}.png`. Columns 48, 24, 12 h before onset, onset, peak. Top
row = mean field (contours of Z500 every 100 m, or white arrows for the wind relative to the storm's motion); bottom row
= mean anomaly from climatology, stippled where the season-block bootstrap passes BH q < 0.05.
- `rot_*_z500`: the 500 hPa low starts to the storm's left and moves onto it; the ridge anomaly ahead and right.
- `rot_*_wind250`, `wind500`: the jet band across the storm and the weaker centre at onset and peak.
- `rot_*_div`: 250 hPa (rows 1 and 2) and 500 hPa (rows 3 and 4) divergence; read the 250 hPa rows.

## Deviations and checks
- Chunks pulled: 1,193, not the 1,356 estimated in the plan (the estimate counted lag fixes that do not exist for
  young tracks). WeatherBench2 34.0 GB + climatology 1.96 GB = 36.0 GB in all, under the 50 GB gate; one test chunk
  was read twice (0.03 GB).
- Rotated-frame wind components are rotated by the heading at the storm and the angle is not adjusted across the 8,000 km
  box (meridian convergence ignored). Points beyond 85 degrees latitude are masked.
- 19 seasons (2004-05 to 2022-23). The onset heading is missing for 27 lows and the peak heading for 18 (they are left
  out of the rotated frame only).
- Survivorship: a lag exists only if the track existed then; at 48 h before onset n is 151 (Atlantic) and 92 (Pacific),
  against 664 and 559 at onset, and they are older storms.
- Verification: see `VERIFICATION.md`.

## Reproduce (about 36 GB streamed; raw fields are not committed)
    python3 -I select.py ../intensity/results/fixes_2004.csv.gz results/storm_times.csv
    ERA5_WORK=$W python3 -I clim.py            # 2.0 GB climatology
    ERA5_WORK=$W python3 -I extract.py 8       # 34 GB, resumable, stops at 48 GB
    ERA5_WORK=$W python3 -I analyse.py 2000    # needs about 6 GB of memory
    python3 -I figs.py
