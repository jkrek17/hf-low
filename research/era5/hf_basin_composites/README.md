# Basin-wide typical weather pattern for HF lows (ERA5 proxy, pipeline A, Oct-Apr)

Descriptive averages, no test and no comparison group (nothing pre-registered; the choices below are defaults set for a request from Jason on 2026-10-11 and can be changed). **ERA5 proxy**: an "HF low" is a pipeline A low whose 800 km ocean gust index reaches 71.7 kt at a 00/12 UTC fix; not the archive. Tropical cyclones in. WeatherBench2 ERA5 at 1.5 degrees.

Population: the 1,223 HF lows of `../hf_env_composites/results/storm_times.csv` (first HF fix, October to April onsets, 2004-05 to 2022-23): 664 North Atlantic, 559 North Pacific. Two times per low: HF onset, and 48 hours before onset (the calendar time, whether or not the track existed then). Maps are Earth-relative (fixed longitude and latitude), not storm-centred. Shading = anomaly from the WeatherBench2 1990-2019 climatology at the same day of year and hour; contours = the mean field. MSLP, 500 hPa height, 250 hPa wind speed.

Pull: 1,332 chunks, 38.72 GB (sized 38.7 GB by HEAD request; WeatherBench2 u, v, z, msl per two-day chunk, 29.08 MB), plus 0.5 GB of climatology. Only sums were kept. Reproduce: `ERA5_WORK=DIR clim_msl.py` (and `../hf_env_composites/clim.py`, `clim_ws.py`), `extract.py 4`, `figs.py`.

## What the average looks like (`summary.txt`, `basin_atl.png`, `basin_pac.png`)
- Atlantic, onset: mean low of about 992 hPa near 63N 33W (6.9 hPa below climate), a 500 hPa trough 80 m below climate near 60N 34W, a jet of up to 48 m/s from the US east coast across the ocean (+6 m/s). Two days before, the same Icelandic-low and trough pattern is already there but weaker (pressure 3.3 hPa and height 45 m below climate) and the jet is 5 m/s stronger than climate over the western ocean. Onsets sit along the southern flank of the trough and low (latitude 39.8 to 64.2N for the middle 80%, mean 53.5N).
- Pacific, onset: mean low of about 1,001 hPa near 42N 150E (2.8 hPa below climate), trough 29 m below climate near 43.5N 145E, jet up to 65 m/s off Japan (+5 m/s). The signal is weaker than in the Atlantic because HF onsets are spread across 130E to 140W, so the average smears (onset latitude 35.2 to 52.2N for the middle 80%, mean 43.4N).
- These are averages of many storms in different places; they show the typical background, not a single storm.

## Correction found on the way
The wind-speed anomaly in the storm-centred composites (PR 146, `hf_env_composites/extract.py`) is |V| minus |V_clim| with V_clim the climatological mean vector. That is biased high everywhere (at 45N 30W on 1 January 12 UTC the climatological mean speed is 37.4 m/s but the speed of the mean vector is 28.8 m/s). The first version of these maps used the same definition and showed +7 to +15 m/s over whole basins. Here the anomaly is the mean |V| minus the climatological mean |V| (WeatherBench2 `wind_speed`, `clim_ws.py`) at each event's own day and hour. The PR 146 speed-anomaly panels (and the two 250 hPa wind figures on the site Composites page built from them) carry the old, biased shading; the mean wind-speed numbers in their captions are not affected. Not yet redone.

## Limits
Population and times are those of PR 146 (WeatherBench2 ends 2023-01-09; later seasons not covered). 1.5 degree grid. The 48 h maps use all 1,223 lows, unlike PR 146's 48 h column (storms that already had a track). Anomaly is vs 1990-2019, so a warming trend adds a little (the global mean 500 hPa height anomaly over 2004-23 is about +13 m).

## Verification
See `VERIFICATION.md`.
