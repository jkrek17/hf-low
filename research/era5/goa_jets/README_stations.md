# Do coastal weather stations back HRRR's extra shoreline gusts? (Gulf of Alaska)

HRRR values are a 3 hour MODEL forecast; ERA5 is a PROXY; the stations are observations from exposed capes and islands (anemometer heights differ and are not adjusted). Plan: `PREREGISTRATION_stations.md` (committed `cd8bfde` before any station value was compared; deviations at its end). Every 00 and 12 UTC time, October to April, 2018-19 to 2025-26 (3,394 of 3,396 times), NDBC records within 30 minutes of the hour. Pipeline: this study's own zone maxima, neither pipeline A nor B. Distances in nautical miles where quoted. Code: `stations_pull.py`, `stations_analyse.py`; tables `results/stations_pairs.csv`, `stations_h1.csv`, `stations_h2.csv`, `stations_analysis.txt`.

## Short answer

Yes, partly: when HRRR reaches the 54 kt level and ERA5 does not, the stations really do see strong gusts far more often than in ordinary weather, but less often than when both agree. And at the stations HRRR does not read higher than the observations. What this does not show: that HRRR's shoreline cell values are accurate, or that the strongest gusts it puts a mile off the coast happen.
- **Stations vs the models (H1, gust 35 kt or more at the station).** HRRR's gust at the station cell runs 0.75 to 0.99 of the station gust (median), ERA5's 0.67 to 1.03; both read at or below the stations, neither above. HRRR is closer to 1 than ERA5 at 4 of 9 stations. Registered rule: closer at 5 or more means "HRRR is the better coastal reference"; 4 means **no better than ERA5** at these sites. Rank agreement is similar for both (Spearman 0.3 to 0.6 for HRRR, 0.3 to 0.75 for ERA5).
- **Are HRRR-only events real (H2)?** Share of times with a station gust of 45 kt or more, by what the two models said about the zone:

| Zone (stations) | HRRR 54 kt, ERA5 not | both 54 kt | neither | HRRR-only vs neither |
|---|---|---|---|---|
| North Gulf (Pilot Rock, Middleton) | 9 in 100 (494 times) | 14 in 100 | 0.9 in 100 | 9.9 times [6.7, 15.4] |
| Panhandle (Cape Decision, Lincoln, Guard) | 12 in 100 (190) | 28 in 100 | 0.5 in 100 | 24 times [15, 52] |
| Cook Inlet (Pilot Rock, Amatuli, Flat, Augustine) | 21 in 100 (467) | 51 in 100 | 1.9 in 100 | 11 times [7.5, 20] |

Registered reading, all three zones: **supported** (3 times or more; the counts could separate ratios of 2.9, 5.5 and 2.3 times). But "neither" is mostly calm weather, a low bar; the better yardstick is that HRRR-only times show about 0.4 to 0.65 of the "both" rate. So HRRR-only events are real stormy times, weaker on average than the ones both models call.

## What it does and does not say

- The stations are a handful of points. A station that sees nothing does not rule out a jet elsewhere in the zone; a station that sees 45 kt does not prove HRRR's peak cell (which sat within about 2 nautical miles of land) was right.
- Both models read lower than the stations at the strongest observed gusts, which fits stations on exposed capes and gust definitions that differ (a short peak gust at 5 to 10 m against a model diagnostic).
- The result supports a reading: ERA5 under-counts storm-level gusts at the coast, and HRRR's higher counts are not just a gust-diagnostic bias at these sites. Event counts between the two models remain not directly comparable (cell counts, definitions).
- Cape Spencer is outside the zone boxes and used in H1 only. Middleton Island has data for two winters only.
- Intervals are 8-season block bootstraps (2,000 draws) and have not been independently recomputed unless a later line here says so.
