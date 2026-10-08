# ERA5 wind-gust climatology and within-era gust rankings

**ERA5 PROXY (reanalysis instantaneous 10 m gust), NOT DIRECT OBSERVATION.** Reanalysis gusts near mountains, coasts and sea ice are not reliable; ERA5 smooths small
intense features. Pipeline A is the gust index used elsewhere in this repository; the grids here are plain ERA5 gust, ocean points only.

Code: [extract.py](extract.py) (ERA5 pull), [analyse.py](analyse.py) (maps, tables, rankings). Results: `results/`.

Reproduce (about 10 minutes, then about 1 minute; the pull is 29.7 GB, under the 50 GB gate):

    ERA5_WORK=<dir> python3 -I research/era5/gust_climo/extract.py 6
    ERA5_WORK=<dir> python3 -I research/era5/gust_climo/analyse.py

## What was done

- **Grids** (fig15, `gust_climo_grid.npz`): ARCO-ERA5 instantaneous 10 m gust at 00 and 12 UTC, 1 October to 30 April, seasons 2004-05 to 2025-26 (9,338 times; 29.7 GB streamed; no earlier cache
  covered it). Per ocean cell: the mean annual maximum gust, the share of sampled times at or above 34, 50, 64, 71.7, 80 and 90 kt, and the highest sampled gust. Sampling is 12-hourly and October to April only, so
  true maxima are higher and the late-summer tropical-cyclone season is **not** in the climatology.
- **Where the strongest gust sits in an HF-strength storm** (fig16): hf_structure's gust-maximum position at 5,983 HF-strength fixes (all months, 2004-05 on, 3,303 Atlantic and 2,680 Pacific; PR 27).
- **Rankings** (`top_gust_storms_by_era.csv`, `gust_climo.txt`): the 25 highest pipeline A gust indices, ranked **within** 1979-2000 and within 2004-2025 only, never against each other (decision 1: pre-2001 gust at fixed depth drifts upward).
  The strongest-storm list by depth and background-adjusted depth belongs to another thread; use it for one list across the whole period.

## What it shows (all in `results/gust_climo.txt`)

- **Atlantic:** the hot spots are the Greenland coast, not the open ocean. The share of sampled times with gust at or above 64 kt peaks at 2.84% near 69N, 21.8W (Jan Mayen and the east Greenland coast); the highest mean annual maximum is 83.9 kt at 65.3N, 36.5W; the highest single sampled gust is
  112.9 kt at 64.8N, 36.3W. The Cape Farewell, Denmark Strait and east Greenland coastal strip stands out in fig15. These are consistent with Greenland tip-jet and barrier winds (PR 29), but ERA5 gusts there are unvalidated (no terrain mask). Excluding cells within about 200 km of land, the Atlantic
  maximum of the mean annual maximum is 73.7 kt at 46.3N, 39.8W, and the highest frequency of gust at or above 64 kt is 1.22% at 59.0N, 40.8W (south of Cape Farewell). 20.1% of Atlantic ocean cells (20-75N) have a mean annual maximum at or above 64 kt, and 1.0% at or above 71.7 kt.
- **Pacific:** weaker and broad. The open-ocean maximum of the mean annual maximum is 71.3 kt at 40.8N, 167.3E (east of Japan, the same region as the busiest HF centre box), and the highest share at or above 64 kt is 0.58% at 40.8N, 166.8E. The only coastal maximum is the Gulf of Anadyr
  area (60.3N, 171.3E; 71.8 kt, 0.69%). 10.8% of Pacific ocean cells (20-75N) have a mean annual maximum at or above 64 kt, and none reaches 71.7 kt. The highest single sampled gust, 116.7 kt at 20.0N, 124.5E, is a tropical cyclone at the southern edge of the domain.
- **In HF-strength storms** the maximum gust sits a median 251 km (Atlantic) and 222 km (Pacific) from the centre, with a median maximum of 76.8 and 76.2 kt. 17.7% of Atlantic and 6.2% of Pacific maxima are within 100 km of a coast. The busiest 5 x 10 degree box for the position of the maximum gust is
  65-70N, 40-30W (70.9 h per season; coastal Greenland, see above) in the Atlantic and 40-45N, 160-170E (46.9 h per season) in the Pacific.
- **Rankings.** Within 1979-2000 the highest gust indices (1-3) are Pacific typhoons transitioning near Japan (7 Nov 1983, 119.6 kt; 2 Oct 1981, 116.1 kt) and Atlantic hurricane-linked 9 Sep 1995 (115.5 kt); the first non-tropical Atlantic storm is 8 Jan 1990 (109.8 kt, 934 hPa) and the 10 Jan 1993 storm (915 hPa) is eleventh. Within 2004-2025 the highest are 23 Sep 2022 (122.9 kt, tropical-cyclone linked), 16 Sep 2004 (115.0 kt, tropical-cyclone linked), then non-tropical
  Atlantic storms on 11 Feb 2011 (114.4 kt) and 7 Feb 2022 (112.9 kt) and the Pacific 14 Dec 2024 storm (112.3 kt). Tropical-cyclone-linked events (flagged) are 15 of the top 25 in 1979-2000 and 11 of the top 25 in 2004-2025; the list is not a list of the deepest extratropical storms.

## Caveats

- The gust climatology is **2004-05 on only** and October to April only; a 1979-2000 gust climatology would be another 30 GB and is not comparable in level (decision 1), so it has not been pulled.
- Instantaneous gust, 12-hourly: a lower bound on true maxima. Not sustained wind. ERA5 gusts at coasts, over terrain and over sea ice (Greenland, Chukotka) are least reliable; fig15's brightest cells are there.
- The within-era ranks use pipeline A's 800 km gust index, which is a storm-level maximum, not a grid-point gust.

## Verification

A fresh Sonnet agent, without access to the code, recomputed: (A) a 100-time slice of the pull from ERA5 itself (the first 50 days of 2004-05, 00 and 12 UTC): the sums and maxima agree to float32 rounding
(8.5e-5 and 7.6e-6 kt) and the six exceedance-count grids are identical in all 6 x 361 x 1,440 cells; (B) from the committed grid, the number of times, each domain's highest mean annual maximum, highest 64-kt frequency and highest sampled gust with their
locations, and the shares of ocean cells at or above 64 and 71.7 kt; (C) the ranks 1-3, 8 and 11 (1979-2000) and 1-5 (2004-2025), and the HF-fix counts, median distance and gust of the maximum, coastal share and busiest box. All matched.

**Not independently checked:** the full 9,338-time aggregation (only a 100-time slice was re-streamed, so the combining of the 94 chunks into the grid is checked only through the committed-grid results above); the open-ocean (no land within about 200 km) figures; the ocean mask (ERA5 land-sea mask below 0.5);
the 15 and 11 tropical-cyclone counts were counted by the verifier without an expected value, and agree with the text; the figures were inspected by eye.
