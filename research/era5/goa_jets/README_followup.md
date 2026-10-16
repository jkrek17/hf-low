# Gulf of Alaska jets, follow-up: mean wind, and what the GFS shows 24 and 48 hours ahead

ERA5 is a PROXY, not observation, and GFS values are model forecasts. Pipeline: this study's own zone maxima (`results/stage0_zone_gust.csv`, shared ERA5 gust cache), not pipeline A or B.
Plans and rules: `PREREGISTRATION_followup.md` (committed 1279f20 before any number was seen; deviations logged at its end). Zone, regime and ingredient definitions are those of
the first study (`README.md`). Rates are the share of 00 and 12 UTC times, October to April. A fresh agent recomputed the point estimates and counts quoted below; it did **not** recompute
the bootstrap intervals, BH q values or season ranges.

## 1. Is the Cook Inlet "sea warmer than the air" result real wind, or the gust formula? (post hoc follow-up of one of 15 tests)

Pre-registered reading: **partly the gust formula**. Among the 109 Cook Inlet northerly (gap) episodes, the first study found 3.7 times the odds of a 64 kt-gust episode per standard
deviation (about 2.2 K) of sea-minus-air temperature difference. Redone with the 10 m mean wind, the strongest 13 episodes by mean wind (43 kt or more) show no link: odds 0.82 per SD
(interval 0.41 to 1.71), and the change in the zone's peak mean wind is -0.2 kt per SD (interval -0.75 to 0.35), so an effect larger than about +0.35 kt per SD is outside the interval.
Only 5 of the 13 gust-strong episodes are among the 13 wind-strong ones. What does move with sea-minus-air difference is the ratio of the zone's peak gust to its peak mean wind: up 0.043 per SD
(interval 0.036 to 0.050) on a typical ratio of 1.49 (the two peaks need not be in the same cell; adjusted for depth only). So where the sea is warmer than the air, ERA5's gust runs larger relative to its
mean wind. Whether real gusts do the same is not tested (no gust observations in the straits). Use: do not read "warm sea, cold air" as a sign of stronger sustained wind in Cook Inlet from this proxy.
The interior-pressure term fared no better (odds 1.49, interval 0.82 to 3.25).
Limit: the mean wind is read at the episode's peak-gust time, not its peak-wind time, which favours the gust. 21 season blocks (one season has no episode) for the intervals.
Files: `results/followup_wind_tests.csv`, `followup_wind.txt`, `followup_wind_episodes.csv`.

## 2. Do GFS forecasts at 24 and 48 hours separate event times (ERA5-proxy zone gust 54 kt or more) from non-events?

5 seasons (2021-22 to 2025-26), 2,122 valid times, 15.01 GB pulled, 4,244 files, none missing (`results/gfs_pull_summary.txt`). The GFS forecast regime and ingredients were computed as for ERA5.
The forecast gust and the proxy truth are both model gusts, so this measures how well the forecast matches the reanalysis proxy, not what a ship felt. Five seasons is small: bootstrap
intervals from 5 season blocks under-cover and some collapse to a point; read them as rough. Ranges across seasons are in `followup_gfs_t1.csv`.

**What tracks the ERA5 proxy best.** In the forecast barrier (North Gulf, Panhandle) or gap (Cook Inlet) regime, the GFS gust forecast alone separates events well at 24 h and about as well at 48 h.
North Gulf barrier regime, 24 h, chance an ERA5-proxy 54 kt event happens: forecast gust in the lower half of forecasts 0 in 311 times; next 30 percent 5 in 100; next 15 percent about 3 in 10;
top 5 percent about 8 in 10 (26 of 32). Panhandle top 5 percent 33 of 33; Cook Inlet 30 of 38 (48 h: 26/32, 31/34, 30/38).

**If you are looking at the pressure map instead of the model gust** (forecast bands, chance of an event; 24 h unless said; bootstrap interval in brackets):
- Coast pile-up (pressure near the coast minus 300 to 400 km out), top 5 percent of forecasts: North Gulf 11 of 32 [22, 79], Panhandle 12 of 33 [27, 48], about 1 in 10 below the top 5 percent.
  Cook Inlet interior-minus-exit pressure, top 5 percent: 16 of 38 [32, 47], against about 4 in 100 in the lower half.
- Low in the nearer half of forecasts (within about 620 to 730 nm of the zone's peak gust; cut 1,145 to 1,347 km by zone): about 15 to 16 in 100 (48 h up to 18), against 0 to 10 in 100 by band farther out.
- Deeper forecast low: about 4 to 5 in 100 in the shallower half, 10 to 14 in 100 in the next band, 15 to 30 in 100 in the deepest 20 percent. At 24 h the top 5 percent is no higher than the 80th to 95th band;
  at 48 h it is higher in North Gulf (25 vs 18) and Cook Inlet (21 vs 15), 32 to 38 times.
- Base rates, forecast barrier regime: North Gulf about 10 in 100, Panhandle 12 in 100; Cook Inlet gap 9 in 100; other regimes 3 to 6 in 100. Strong (64 kt or more): 1 to 2 in 100 in barrier and gap regimes, 0 to 1 elsewhere.

**Do the pressure ingredients add to the model gust?** Cannot tell, and no gain was seen. Leave-one-season-out mean log-loss, gust alone versus gust plus depth, low distance and the coast or interior
pressure term: North Gulf 24 h 0.154 vs 0.160, 48 h 0.168 vs 0.169; Panhandle 0.141 vs 0.140, 0.172 vs 0.173; Cook Inlet 0.128 vs 0.134, 0.138 vs 0.149. In no cell did the extra terms improve
out of sample; the intervals rule out a gain larger than about 0.006 in log-loss (about 4 percent of the gust-only value). At Cook Inlet 48 h the extra terms were **worse** (+0.011, interval 0.002 to 0.020, all 5 seasons agree).
AUC for gust alone is 0.94 to 0.97. The pressure bands are associated with higher event rates (above); using them as a cross-check on the gust is a suggestion, not something tested here.
Files: `results/followup_gfs_t1.csv`, `followup_gfs_t2.csv`, `followup_gfs.txt`, `gfs_features.csv`, `gfs_pull_summary.txt`.

## Not done / limits
- T3 (do the forecast bands order like the ERA5 same-time table) was not run; the climatology AUC was not reported (it is 0.5 within a regime by construction).
- The ingredients came from the 22-season ERA5 study that includes these 5 seasons, so the field source is independent, the choice of ingredients is not. These 5 seasons are also the gust-calibration
  seasons of the outlook work.
- GFS model changes inside 2021-26 were not reconstructed. GFS gust is not corrected to ERA5 gust. STAB (needs SST) was not used. GFS 0.25 degree smooths Shelikof and Cook Inlet like ERA5.
- An ERA5-proxy 54 kt reads as a buoy median gust of 44 to 48 kt (first study).
- ASCAT: the existing cache (`/mnt/project-files/ascat-check/cells`) holds 0 of its 15,255,551 cells in the Gulf of Alaska zones; a fresh RSS pull for the 1,380 episode days is 2,139 files, 3.28 GB
  (`ascat_sizing.py`, `results/ascat_sizing.txt`). Not run.
