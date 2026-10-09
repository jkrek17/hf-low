# Storm-centred composites of every HF low: upper-level environment (ERA5 proxy, pipeline A)

Committed before any field is read. Descriptive study: the average environment around each HF low, no
hypothesis test against another group. Jason asked for it directly (2026-10-09: "a storm centered plot with
different meteorology variables composited for every HF low"; Oct-Apr is fine).

**Proxy.** An HF low is a pipeline A (`../hf_history`) low whose ERA5 800 km ocean gust index reaches 71.7 kt
at a 00/12 UTC fix. Not the archive. Tropical cyclones are in.

## What had been looked at
The existing storm-scale HF composites (`../hf_vs_storm`) are surface fields only. The only storm-relative
250 hPa composite is for deepening storms that went on to HF (`../hf_conversion`, 250 pairs per basin).
PR 12 coefficients, PR 78 jet and trough results and the fixed-box Z500/U250 composites were seen. No
storm-relative Z500 or 500 hPa composite has been computed.

## Population (`select.py`, from `../intensity/results/fixes_2004.csv.gz`)
- HF low = a track with at least one 00/12 UTC fix with g800 >= 71.7 kt. Onset = its first such fix. Peak =
  its fix with the largest g800.
- October to April onsets only; onset and peak before 2023-01-10 (WeatherBench2 ends 2023-01-09); seasons
  2004-05 on. The post-2023 gap is left for later (0.25 deg store needed).
- Times: onset, peak, and the same track's fix 12, 24 and 48 h before onset (missing for young tracks; n is
  reported at each time). Atlantic and Pacific separately.

## Fields (WeatherBench2 1.5 deg, 6-hourly) and climatology
Z500 (m); u, v at 250 and 500 hPa, speed (m/s); divergence of the wind at 250 and 500 hPa (1e-5 /s).
Anomaly = field minus the WeatherBench2 ERA5 climatology (1990-2019, 6-hourly by day of year and hour, as
published): u, v and Z500 anomalies directly; speed anomaly = |V| - |V_clim|; divergence anomaly =
divergence of the anomaly wind.

## Frames
Storm-centred at the fix. Box +-4000 km, 100 km spacing, 81 x 81, bilinear on the 1.5 deg grid. Two frames:
(1) rotated, +x along the motion over the previous 6 h (the fix table's `heading`), +y to the left; vector
components rotated into along and cross; storms with no heading at that time are left out of this frame
only. (2) north-up, +x east, +y north. Points beyond 85 degrees latitude are masked. The rotated frame is the
primary one.

## Statistics
Mean over storms of each field. Anomaly composites: season-block bootstrap (2,000 resamples of seasons with
replacement, seed 20261009), two-sided pixelwise p against zero, Benjamini-Hochberg q < 0.05 within each
field x basin x time; stippling is descriptive because pixels are correlated. The independent sample is
seasons, not storms.

Descriptive scalars, defined now:
- S1 value at the centre of the mean 250 hPa speed and 500 hPa speed (rotated frame).
- S2 location (x, y) and value of the maximum of the mean raw 250 hPa speed within 2500 km (rotated frame),
  at each time.
- S3 per-storm location of the maximum raw 250 hPa speed within 2000 km (the storm's own upper circulation is
  not removed, so a storm-centred maximum can be the storm itself; reported as is): share ahead (x > 0) and
  behind, median distance, with season-block intervals.
- S4 location and value of the minimum of the mean Z500 anomaly within 3000 km, and of the maximum 250 hPa
  divergence.
No threshold, group or window is changed after the composites are seen; a change is logged as a post hoc
deviation.

## Limits stated in advance
1.5 deg smooths jet cores. Averaging smears streaks and troughs, so quadrant or tilt statements cannot be
made from these composites. Frames are tied to the track heading, not to the jet or the steering flow. Gust
index values drift before 2001 (decision 1); the window starts 2004-05.

## Pull
WeatherBench2 u, v, z whole chunks for the 1,356 chunks that hold the needed times (about 38.7 GB at 28.5 MB
per chunk, measured on 30 chunks) plus about 1.5 GB of the climatology. Total about 40 GB, under the 50 GB
gate; the run stops if the running total passes 48 GB. Raw boxes stay in `$ERA5_WORK`, not committed.
