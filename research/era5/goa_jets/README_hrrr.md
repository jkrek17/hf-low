# Gulf of Alaska jets seen by HRRR-Alaska (3 km model): does the 0.25 degree ERA5 proxy miss them?

HRRR values are a 3 hour MODEL forecast, not observation; ERA5 is a PROXY. Pipeline: this study's own zone maxima (zones, regimes and definitions of the first study, `README.md`), neither pipeline A nor B.
Plan: `PREREGISTRATION_hrrr.md` (committed 5b43cbc before any HRRR field was pulled; deviations logged at its end). Every 00 and 12 UTC time, October to April, 2018-19 to 2025-26: 3,394 of 3,396 times
(2 files missing, both in March 2019), 12.64 GB streamed, kept at `/mnt/project-files/goa-jets/hrrr/`. A fresh agent recomputed the point estimates and counts quoted here; it did **not** recompute the
bootstrap intervals or season ranges. "HRRR 54 kt" is the first study's 54 kt applied to HRRR's own gust, which is a different diagnostic; it is not ERA5's 54 kt and is not calibrated.

## Short answer

Away from the coast ERA5 is not badly wrong; the large difference is in a narrow band along the shoreline, where this model has no observation to check it.
- **At the seven buoys** (offshore, 0.7 to 1.7 km from the nearest HRRR cell), on times chosen by the buoy gust (35 kt or more, 123 to 244 times each), HRRR gust is 0.99 to 1.15 times the buoy gust (medians) and ERA5 gust is 0.96 to 1.14 times. Both read about 10 percent high; neither reads low. (Buoy anemometers are about 5 m up and are not adjusted.)
- **HRRR sees many more storm-force zone maxima than ERA5** (North Gulf 741 vs 226 times in 3,394; Panhandle 376 vs 205; Cook Inlet 719 vs 222), and nearly every ERA5 event is also an HRRR event (92, 88 and 98 percent). Part of that is a maximum taken over about 40 times as many cells (3 km against 0.25 degree).
- **Most of the extra events sit within 10 km of land.** For HRRR events in southeasterly flow along the North Gulf, the peak-gust cell is a median 3 km from land, 74 percent within 10 km (ERA5: 39 km, 27 percent); Cook Inlet northerly flow: 3 km, 81 percent within 10 km (ERA5 34 km, 6 percent). Panhandle peaks are far out (median 100 km).
- **Leave out cells within 25 km of land** (post hoc) and the gap shrinks: HRRR zone gust is 1.04 to 1.09 times ERA5's, and 54 kt happens at 15.6 percent of North Gulf barrier times, 16.2 Panhandle, 13.2 Cook Inlet gap, against 10.1, 11.3 and 9.2 for ERA5.
Whether the shoreline gusts are real jet cores or an artefact of the model's gust diagnostic over mixed land-sea cells cannot be told from these data. A coastal station (Middleton Island, Cape Hinchinbrook and similar) or SAR winds would settle it; none is in this study.

## H1: HRRR against ERA5 zone-maximum gust, same times

| Zone, flow | Spearman | HRRR/ERA5 median [8-season interval] | 54 kt: HRRR / ERA5 / both | 64 kt: HRRR / ERA5 |
|---|---|---|---|---|
| North Gulf, all | 0.82 | 1.19 [1.18, 1.21] | 741 / 226 / 207 | 265 / 28 |
| Panhandle, all | 0.87 | 1.12 [1.11, 1.13] | 376 / 205 / 181 | 81 / 30 |
| Cook Inlet, all | 0.90 | 1.17 [1.16, 1.18] | 719 / 222 / 217 | 241 / 20 |

By flow the ratios are 1.15 to 1.27 (North Gulf), 1.10 to 1.17 (Panhandle), 1.14 to 1.22 (Cook Inlet) (`results/hrrr_h1.csv`). HRRR sees 64 kt (strong) 3 to 12 times as often as ERA5.

## H3: buoys (results/hrrr_h3.csv)

Per buoy, buoy gust 35 kt or more: HRRR gust / buoy gust median 0.99 (46077) to 1.15 (46082), quartiles 0.90 to 1.27; HRRR wind / buoy wind 0.99 to 1.11; ERA5 gust / buoy gust 0.96 to 1.14 on the same times.
At buoy gust 45 kt or more (9 to 39 times per buoy) HRRR's zone gust was 54 kt or more at 79 to 100 percent of those times (35 kt or more: 61 to 84 percent). Each buoy's HRRR cell is water (land flag false).

## H4: chance of an HRRR zone gust of 54 kt or more (every 00/12 UTC time; exact rates; 8 season blocks)

| Zone, flow | times | HRRR 54 kt | HRRR 64 kt | ERA5 54 kt, same times |
|---|---|---|---|---|
| North Gulf, southeasterly (barrier) | 989 | 28.7 in 100 [26.3, 30.7] | 10.7 | 10.1 |
| Panhandle, southeasterly | 842 | 17.0 [14.4, 20.7] | 4.6 | 11.3 |
| Cook Inlet, northerly (gap) | 1,055 | 24.9 [19.7, 28.9] | 8.9 | 9.2 |

These include the shoreline band; the 25 km figures are in the short answer. By ingredient, top 5 percent against the lower half of the same flow (forecast-style table, same-time values):
- Panhandle coast pile-up: 86 in 100 [74, 96] (43 times) against 4.5 in 100 [2.6, 6.5]; the 80th to 95th band is 44 in 100.
- North Gulf coast pile-up: 60 in 100 [47, 78] (50 times) against 28 in 100 [25, 31]; only the top 5 percent stands out.
- Cook Inlet interior-minus-exit pressure: 64 in 100 [50, 77] (53 times) against 22 in 100 [17, 28].
Season ranges are wide (the top-5-percent cells run 0 to 100 percent across seasons because each holds 3 to 12 times), so read the bands as ordering, not precise rates. Files: `results/hrrr_h4.csv`, `hrrr_analysis.txt`.

## H2: where the peak sits (results/hrrr_h2.csv)
Distance from the peak-gust cell to the nearest land for events in barrier or gap flow, same land mask for both: North Gulf HRRR 284 events median 3.0 km (79 percent within 25 km), ERA5 100 events 38.6 km (30 percent); Panhandle HRRR 143 events 100 km (29 percent), ERA5 95 events 50 km (28 percent); Cook Inlet HRRR 263 events 3.0 km (90 percent), ERA5 97 events 34 km (28 percent).

## Not done / limits
- Strong-versus-storm odds ratios, depth and low distance (the HRRR domain covers only 66 percent of the reference-low box), and stability were not run; this is a 3 hour forecast, not a forecast skill test.
- Eight seasons; HRRR-Alaska versions changed in 2018-2026 and were not reconstructed; HRRR assimilates little over the open Gulf; the maximum over many more cells is larger for the same field, so event counts are not comparable between grids without the near-shore and cell-count caveats above.
- The near-coast sensitivity (`hrrr_sensitivity.py`, `results/hrrr_sens_summary.csv`) is post hoc.
- Buoy gust is a 2 to 3 second value at about 5 m; model gust is a diagnostic at 10 m.
