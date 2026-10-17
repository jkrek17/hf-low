# Gulf of Alaska jets seen by HRRR-Alaska (3 km model): does the 0.25 degree ERA5 proxy miss them?

HRRR values are a 3 hour MODEL forecast, not observation; ERA5 is a PROXY. Pipeline: this study's own zone maxima (zones, regimes and definitions of the first study, `README.md`), neither pipeline A nor B.
Plan: `PREREGISTRATION_hrrr.md` (committed 5b43cbc before the full pull; a few records had been decoded earlier to check the grid, packing and wind rotation; deviations are logged at its end). Every 00 and 12 UTC time, October to April, 2018-19 to 2025-26: 3,394 of 3,396 times
(2 files missing, both in March 2019), 12.64 GB streamed, kept at `/mnt/project-files/goa-jets/hrrr/`. A fresh agent recomputed the point estimates and counts quoted here; the bootstrap intervals and season ranges were **not** independently recomputed.
"HRRR 54 kt" means the first study's 54 kt applied to HRRR's own gust, a different diagnostic; it is not ERA5's 54 kt and is not calibrated. Distances are in kilometres in the result files (1 nautical mile = 1.85 km).

## Short answer

HRRR reads higher than ERA5 for the same zone, mostly in a thin band along the shoreline in the North Gulf and Cook Inlet. These data cannot say whether the shoreline gusts are real; the Panhandle difference is not a shoreline effect.
- **At the seven buoys** (HRRR cell 0.7 to 1.7 km from the buoy; but the buoys are 9 to 286 km from land, so they do not test the shoreline band), on times chosen by the buoy gust (35 kt or more, 123 to 244 times each), the median HRRR gust is 0.99 to 1.15 times the buoy gust and the median ERA5 gust 0.96 to 1.14 times. The two models read about the same against the buoys. The buoy gust is a short peak at about 5 m, the model gust a 10 m diagnostic, so a ratio near 1 is not proof that either is right.
- **HRRR reaches the 54 kt level more often than ERA5** (North Gulf 741 vs 226 times in 3,394; Panhandle 376 vs 205; Cook Inlet 719 vs 222), and nearly every ERA5 event is also an HRRR event (92, 88 and 98 percent). HRRR's zone holds about 44 to 48 times as many ocean cells (`results/hrrr_support.txt`); a maximum over more cells is larger, but that was not tested as the cause.
- **In the North Gulf and Cook Inlet most of the extra events sit at the coast.** Of the times HRRR reaches 54 kt and ERA5 does not, the peak-gust cell is within 10 km of land at 81 percent (North Gulf) and 78 percent (Cook Inlet), median 3 km. In the Panhandle it is 33 percent, median 49 km. For events in barrier or gap flow, North Gulf southeasterly flow: HRRR peak median 3 km from land, 74 percent within 10 km (ERA5: 39 km, 27 percent); Cook Inlet northerly flow: 3 km, 81 percent (ERA5 34 km, 6 percent); Panhandle peaks are far out (median 100 km).
- **Leave out HRRR cells within 25 km of land** (post hoc; ERA5 left unrestricted, so this is not a like-for-like comparison) and the gap shrinks: HRRR zone gust is 1.04 to 1.09 times ERA5's, and 54 kt happens at 15.6 percent of North Gulf barrier times, 16.2 Panhandle, 13.2 Cook Inlet gap, against 10.1, 11.3 and 9.2 for ERA5.
A coastal station (Middleton Island, Cape Hinchinbrook and similar) or SAR winds would settle whether the shoreline gusts are real; none is in this study.

## H1: HRRR against ERA5 zone-maximum gust, same times

| Zone, flow | Spearman | HRRR/ERA5 median [8-season interval] | 54 kt: HRRR / ERA5 / both | 64 kt: HRRR / ERA5 |
|---|---|---|---|---|
| North Gulf, all | 0.82 | 1.19 [1.18, 1.21] | 741 / 226 / 207 | 265 / 28 |
| Panhandle, all | 0.87 | 1.12 [1.11, 1.13] | 376 / 205 / 181 | 81 / 30 |
| Cook Inlet, all | 0.90 | 1.17 [1.16, 1.18] | 719 / 222 / 217 | 241 / 20 |

By flow the ratios are 1.15 to 1.27 (North Gulf), 1.10 to 1.17 (Panhandle), 1.14 to 1.22 (Cook Inlet) (`results/hrrr_h1.csv`). HRRR reaches 64 kt 2.7 to 12 times as often as ERA5 (20 times for Cook Inlet inflow; `results/hrrr_support.txt`), again mostly a shoreline-band effect in the North Gulf and Cook Inlet.

## H3: buoys (results/hrrr_h3.csv)

Per buoy, buoy gust 35 kt or more: HRRR gust / buoy gust median 0.99 (46077) to 1.15 (46082), quartiles 0.90 to 1.27; HRRR wind / buoy wind 0.99 to 1.11; ERA5 gust / buoy gust 0.96 to 1.14 on the same times.
At buoy gust 45 kt or more (9 to 39 times per buoy) HRRR's zone gust was 54 kt or more at 79 to 100 percent of those times (35 kt or more: 61 to 84 percent). Each buoy's HRRR cell is water (land flag false).

## H4: chance of an HRRR zone gust of 54 kt or more (every 00/12 UTC time; exact rates; 8 season blocks)

| Zone, flow | times | HRRR 54 kt | HRRR 64 kt | ERA5 54 kt, same times |
|---|---|---|---|---|
| North Gulf, southeasterly (barrier) | 989 | 28.7 in 100 [26.3, 30.7] | 10.7 | 10.1 |
| Panhandle, southeasterly | 842 | 17.0 [14.4, 20.7] | 4.6 | 11.3 |
| Cook Inlet, northerly (gap) | 1,055 | 24.9 [19.7, 28.9] | 8.9 | 9.2 |

These include the shoreline band; the 25 km figures are in the short answer. By ingredient, top 5 percent against the lower half of the same flow (same-time values, not a forecast test):
- Panhandle coast pile-up: 86 in 100 [74, 96] (43 times) against 4.5 in 100 [2.6, 6.5]; the 80th to 95th band is 44 in 100.
- North Gulf coast pile-up: 60 in 100 [47, 78] (50 times) against 28 in 100 [25, 31]; only the top 5 percent stands out.
- Cook Inlet interior-minus-exit pressure: 64 in 100 [50, 77] (53 times) against 22 in 100 [17, 28].
Season ranges are wide: each season holds only 1 to 20 times in the top-5-percent cell (many hold 2 to 5), and the share reaching 54 kt runs from 0 to 100 percent across seasons, so read the bands as ordering, not precise rates (`results/hrrr_support.txt`). The intervals in this table were not independently recomputed. Files: `results/hrrr_h4.csv`, `hrrr_analysis.txt`.

## H2: where the peak sits (results/hrrr_h2.csv)
Distance from the peak-gust cell to the nearest land for events in barrier or gap flow, same land mask for both: North Gulf HRRR 284 events median 3.0 km (79 percent within 25 km), ERA5 100 events 38.6 km (30 percent); Panhandle HRRR 143 events 100 km (29 percent), ERA5 95 events 50 km (28 percent); Cook Inlet HRRR 263 events 3.0 km (90 percent), ERA5 97 events 34 km (28 percent).

## Not done / limits
- Strong-versus-storm odds ratios, depth and low distance (the HRRR domain covers only 66 percent of the reference-low box), and stability were not run; this is a 3 hour forecast, not a forecast skill test.
- Eight seasons; HRRR-Alaska versions changed in 2018-2026 and were not reconstructed; HRRR assimilates little over the open Gulf; the maximum over many more cells is larger for the same field, so event counts are not comparable between grids without the near-shore and cell-count caveats above.
- The near-coast sensitivity (`hrrr_sensitivity.py`, `results/hrrr_sens_summary.csv`) is post hoc.
- Buoy gust is a short-period peak (NDBC: 5 or 8 seconds) at about 5 m and is not adjusted to 10 m; the model gust is a 10 m diagnostic. The buoys are far from the shoreline band.
