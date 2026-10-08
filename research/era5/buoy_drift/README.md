# ERA5 surface wind against moored buoys, 1979–2004

Does ERA5's own surface wind, or its gust, strengthen relative to in-situ wind over 1979–2000? `research/era5/drift_check.py` found that the ERA5 gust index rises at fixed storm depth over those seasons in both pipelines (A, `research/era5/hf_history`: +0.665 kt/decade on a 64.38 kt mean; B, `event_fields`/`criterion`/`series`: +0.645 kt/decade on 58.40 kt; both about +1.0 to +1.1 % per decade). If ERA5's winds drift against the real wind, that would explain it and the gust index could not be carried back. If they do not, the drift lies elsewhere: in the storms, or in something the buoys cannot see.

ERA5 is a **proxy** record here as everywhere in this project. The buoys test it; nothing here validates the HF proxy itself before the archive.

Results: [`results/buoy_drift-result.txt`](results/buoy_drift-result.txt). Station-season values behind every trend: [`results/station_seasons.csv`](results/station_seasons.csv).

## Findings (2026-10-08)

Section numbers refer to `results/buoy_drift-result.txt`. Percent per decade is 100 × the slope of a log ratio. U_e and G_e are ERA5 wind and gust, U_b the buoy wind at 10 m.

1. **At the buoys, ERA5's surface wind and gust do not strengthen against the in-situ wind over 1979–2000.** With one fixed effect per station segment and a season-clustered se (section 4), ERA5 gust against buoy wind at pair-mean winds of 15 m/s and above moves −1.80 %/decade (se 1.92, t = −0.94; leave-one-season-out −3.19 to −1.21). ERA5 wind against buoy wind in the same bin moves −1.13 %/decade (se 1.73, t = −0.66). Letting each station have its own slope, the gust ratio averages −3.80 %/decade (se 2.27, t = −1.67), negative at 8 of 12 stations.
2. **The test cannot detect a drift as small as the gust-index drift.** Pipeline A's drift is +1.03 %/decade and B's +1.10 %/decade (section 6). The buoy standard errors are 1.7–2.3 %/decade, so power against +1.07 %/decade is under 10 %. The buoy record rules out an upward drift of ERA5 against the buoys larger than about +2.2 %/decade (upper 95 % limit, section 6), not one of +1 %/decade.
3. **Where the buoy record is cleanest, the movement is the other way.** Over 1985–2000 (a window chosen after seeing station steps of 10–20 % in the early 1980s, so read it as a check), ERA5 against the buoys falls in the 15 m/s bin at all 11 stations with eight or more seasons in that window: wind −6.83 %/decade (se 1.04, t = −6.59), gust −8.10 %/decade (se 1.28, t = −6.31) (section 9). Pipeline B's gust at fixed depth rises +1.16 %/decade over the same seasons (t = +1.46).
4. **The same holds next to deep lows.** For pairs with a pipeline B centre of 955–975 hPa within 500 km (2,124 pairs, 1979–2000), ERA5 gust against buoy wind moves −3.49 %/decade (se 1.75, t = −2.00); over 1985–2000, −6.01 %/decade (t = −3.32) (section 10).
5. **ERA5's own gust factor falls.** G_e/U_e, which involves no buoy at all, declines −1.53 %/decade over 1979–2000 (segment fixed effects, t = −2.86) and is negative at 9 of 12 stations. So the gust parameterisation is not inflating ordinary ocean gusts over time.
6. **Season-to-season, the buoy comparison does not track the gust-index drift.** The correlation between the buoys' season effects (gust ratio, 15 m/s bin) and B's fixed-depth gust is +0.13 over 1979–2000 and −0.44 over 1985–2000 (section 8).
7. **The one positive result is not credible as an ERA5 drift.** In the gale bin (pair mean of 17.5 m/s and above) the gust ratio rises +4.55 %/decade (t = +2.31, segment fixed effects). Before 1985 that bin holds one to five station-seasons per season, mostly from stations with the early-1980s steps, and dropping one station brings it to +0.23 %/decade.

**What this means for the gust drift.** Nothing at the buoys says that ERA5's surface winds or gusts rose against reality before 2001. If anything ERA5 fell against the buoys over 1985–2000, by three to seven times the size of the index drift. That decline is uniform across stations, which points to ERA5 (the 1990s brought scatterometer and microwave winds into the assimilation), but uniform buoy hardware or reporting changes could also produce it, and the NDBC hull and sensor history needed to tell the two apart was not reachable. Either way the gust-index drift at fixed depth is not explained by a general upward drift in ERA5's surface wind. It may sit in the storms themselves, or in ERA5's storm cores away from the buoys. The buoys cannot say which.

## Run

```
python3 research/era5/buoy_drift/fetch_buoys.py     # ISD buoy reports -> work/isd/   (~1.5 GB streamed, 30 MB kept)
python3 research/era5/buoy_drift/extract_era5.py    # ERA5 at the buoys -> work/era5/ (~30 GB streamed, 20 MB kept)
python3 research/era5/buoy_drift/analyse.py         # -> results/
```

`work/` is ignored. Both fetchers are resumable. Needs `numpy`, `numcodecs` and `lz4`.

## Data

**Buoys.** NOAA's Integrated Surface Database, global-hourly CSV on AWS (`s3://noaa-global-hourly-pds`), which carries the GTS reports (FM-18 BUOY) of NDBC and Canadian moored buoys under USAF ids 99xxxx. The NDBC, MEDS, ICOADS and Copernicus hosts were all refused by the work container's network policy on 2026-10-08, so ISD is the in-situ source. 26 buoys in the two archive basins were fetched; the ocean screen in `analyse.py` keeps those whose nearest ERA5 point and its eight neighbours are all sea.

**ISD holds no moored-buoy reports after 2004.** The comparison period after the drift therefore covers only seasons 2001–2004 (2004 is October to December only), not 2004 onward. A longer later period needs NDBC's own archive, which was not reachable.

**ERA5.** Google's ARCO-ERA5 stores: 10 m wind speed from the 6-hourly store, instantaneous 10 m gust (the gust variable both pipelines use) from the hourly store, at 00/06/12/18 UTC, October to March, 1979–2004. Each Zarr chunk is one global field; `extract_era5.py` fetches only the two blosc blocks that hold 22–67 N, which cuts the pull from about 120 GB to about 30 GB. The block decoder was checked bit-for-bit against a full-chunk decode for both variables.

## Buoy metadata and how breaks are handled

- **Averaging period.** NDBC moored buoys report an 8-minute mean wind. ERA5 10 m wind is an instantaneous model value representing a grid-box mean. The difference is a constant offset, not a trend.
- **Resolution.** For the NDBC buoys with records from the 1970s and 1980s, every FM-18 wind in ISD is a whole number of m/s in every era. The Canadian buoys, South Nomad, 44014 and 44025 mix whole m/s with values converted from knots in the 1990s and 2000s. Rounding adds noise; a change of reporting unit can add a small step, which the break search below is there to catch.
- **Anemometer height.** ISD does not carry it, and NDBC's per-deployment hull history was not reachable. Every buoy is adjusted from an assumed 5 m (3 m discus, 6 m NOMAD) to 10 m with a neutral log profile, z0 = 2e-4 m, a factor of 1.068. A 10 m and 12 m discus hull carries its anemometer near 10 m; for such a station the absolute bias here is about 6 % too high. A constant height factor multiplies a station's whole record and is absorbed by its fixed effect, so it cannot create or remove a trend.
- **Hull or sensor changes** change the height factor partway through a record, which is a step. Steps are found in each station's monthly ln(ERA5/buoy) after subtracting the cross-station median of the same month, so a drift common to all stations (which is what an ERA5 drift would be) is not mistaken for a station break. Trends are then reported twice: with one fixed effect per station, and with one per station segment between detected breaks.
- **Gust definition.** The buoy gust groups in ISD change code between eras (type 3 in the 1980s, type 5 from 1990, type 4 added around 2000), so the buoy gust is not a consistent quantity and is not used. ERA5 gust is compared with the buoy's mean wind instead, and split into ERA5's wind bias and ERA5's own gust factor: ln(G_e/U_b) = ln(U_e/U_b) + ln(G_e/U_e).

## Limits

- **The buoys are not independent of ERA5.** As ERA5's documentation describes it (Hersbach et al. 2020; not checked from here), ERA5 assimilates surface wind and pressure from buoy reports. At a buoy, ERA5 is pulled toward that buoy. A drift in ERA5 away from observations, such as in storm cores far from any buoy, could be larger than the drift visible at the buoys. A null result here means ERA5's surface wind did not drift where it was observed. It does not rule out drift elsewhere.
- **Point winds, not storm maxima.** The gust index is a storm maximum within 500 or 800 km. The buoys sample the whole cold-season distribution at fixed points; the highest bins hold few pairs per station-season.
- **Sample size.** The independent sample for the trend is seasons (22 for 1979–2000), and every t uses the season count. Stations are not independent of each other within a season.
- **Era coverage.** Fewer stations report before 1985; the Canadian buoys begin in 1988–1991.
