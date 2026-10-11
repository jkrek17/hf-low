# Verification (basin composites, pipeline A ERA5 proxy, descriptive)

A fresh agent recomputed from the WeatherBench2 stores and the committed files only (no notes read, nothing changed in the repo).

Recomputed and matching: onset counts 664 Atlantic / 559 Pacific; onset latitude mean 53.5N / 43.4N and 10th-90th percentiles 39.8-64.2N / 35.2-52.2N; onset MSL mean field and anomaly, both basins (max difference 6e-5 hPa); Atlantic onset Z500 field and anomaly (max difference 2e-4 m); deepest Atlantic Z500 anomaly (-80.2 m, 60.0N); deepest MSL anomalies (-6.9 hPa at 63.0N Atlantic, -2.8 hPa at 42.0N Pacific); climatological 250 hPa mean speed at 45N 30W, 1 Jan 12 UTC (37.38 m/s) against the speed of the mean vector (28.76 m/s), which is the bias that made the first wind-anomaly shading too high.

Found and fixed: the plotted longitude axis was not on the 1.5-degree grid, so the maps and the quoted longitudes were up to 0.5 degrees off. `figs.py` now snaps the axis to grid points; values are unchanged (positions moved 0.5 degrees; Pacific Z500 field minimum 5110.0 m instead of 5109.3 m because the window edge moved).

Not independently recomputed: 250 hPa wind-speed mean and anomaly fields (about 19 GB of u, v) including the +6.4 / +5.4 m/s maxima and the 47.5 / 65.1 m/s field maxima; all 48-hour-before fields; Pacific Z500 and the Z500 field minima; the claim that the first wind-anomaly version showed +7 to +15 m/s; the 1,332 chunks / 38.72 GB; whether the PNGs show these fields (checked by eye only).

Known issue: the wind250 anomaly rows of the storm-centred composites (PR 146, and the site figures built from it) use the biased reference (speed of the climatological mean vector); they are not redone.
