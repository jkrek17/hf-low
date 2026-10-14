# Pre-registration: calibrated outlook shading from GFS gust ("OPC agreed"), and a Replay page

Written 2026-10-14 before any GFS gust field of this study was pulled. Jason: "do an outlook type of graphic", "go ahead and build it", card choice "Replay page". **GFS gust is a model forecast, ERA5 is not used here; OPC's call is an analyst's judgement, so the target is "OPC agreed", never "verified". OPC warns on sustained wind; the GFS gust index is a proxy for that.** Existing caches checked first: the GFS fix tables (`gfs-raw-gust*`) keep no gust fields, so a daily maximum gust grid has to be pulled once (about 11 GB, under the 50 GB gate) and saved as a compact reusable product.

## Question
For a point on the map and forecast day 1 to 4, how often did OPC have a hurricane-force (or storm-force) low near that point when the GFS maximum gust within the day's radius was in a given band? This is a **neighbourhood** frequency read straight from the field, not a per-low probability (the advisor panel's warning), so storms the tracker would miss still count.

## Data
- GFS 0.25 degree `GUST:surface`, 00 UTC cycles only, 1 Oct to 30 Apr, seasons 2021-22 to 2025-26 (1,061 cycles), valid every 6 h from f006 to f096, by `.idx` byte ranges from the NOAA open-data bucket (method of `gfs_raw_gust_96h/a_fetch.py`). Decoded with our own GRIB2 template 5.3 reader (`grib2_t53.py`; eccodes cannot be installed here). Reader checked against ERA5 gust for one valid time (see README).
- Stored per cycle: **forecast-day maximum** over each day's four 6-hourly samples (day d = valid hours 24(d-1)+6 .. 24d), on a 0.5 degree grid, 75N to 20N, ocean only (ERA5 land-sea mask; land cells set to 0 before taking the 2x2 maximum), kt in 0.5 kt steps.
- Truth: OPC archive fixes (`data/hf_lows/HF_Data_-_{Atl,Pac}.csv`), 6-hourly. Event **HF-agreed** = a Category `HF` fix; event **storm-agreed** = Category `S` or `HF`. `DHF` (developing HF), `G`, `DS`, tropical categories are not events. The archive holds low centres, not wind areas.
- A day-d window includes the fixes at the four valid times above. A point's truth = at least one event fix within **324 nm (600 km, the radius of all earlier GFS-versus-OPC tables)** during that window. Predictor = the maximum of the day-d daily-maximum gust over ocean cells within the **day's radius R_d = 60, 120, 270, 270 nm for days 1 to 4** (PR 210; day 4 reuses 270, its 96 h lead measured at that radius).

## Calibration
- Evaluation points: ocean cells of the 0.5 degree grid, every other cell in both directions (1 degree spacing) to limit size; North Atlantic 100W to 10E and North Pacific 120E to 120W, 20N to 75N.
- Gust bands (kt): <40, 40-48, 48-56, 56-64, 64-72, 72-80, >=80. For each basin, day and band: the share of point-days with an event, and n (point-days) and n of cycles with any event. Pooled over basins as a second view.
- Leave-one-season-out: bands are fixed in advance; each season's table is built from the other four and scored on the held-out one with the Brier score and Brier skill against the held-out season's base rate. With **five seasons there is no meaningful bootstrap interval**: the per-season values are reported and the range between them stands in for uncertainty. Point-days within a cycle and across days are strongly dependent; the effective sample is cycles and seasons.
- Reliability: predicted band probability versus observed share, held-out.
- A band is reported as usable only if it holds at least 100 point-days and at least 5 distinct cycles with an event in every held-out fold; otherwise "too few".

## Replay cases (fixed now, before looking)
Per season: the 2 cycles with the most HF-agreed fixes in the following 96 h; plus 1 cycle with no HF fix in 96 h but the highest GFS day-2 gust band (the "false alarm" picture, language: OPC did not call one); plus 1 cycle chosen at random (seeded 20261014). 4 per season, 20 in all. Each is shown with the real OPC fixes on top, no selection of attractive cases.

## What would make the shading wrong
Day 1 and 2 shading is narrow (60/120 nm) and was sized for position error only (PR 210, upper bound); intensity error is in this calibration, so the observed shares absorb it. Gust at 6-hourly samples can miss a short peak. Version changes in GFS happened inside the five seasons. ERA5 is not involved.

## Not claimed
No claim about warnings or sustained wind, no skill claim beyond the held-out Brier numbers, nothing for June to September, nothing beyond day 4. Days 5 to 7 on the mock are layout only.

## Deviations (post hoc)
- Replay picks: the two 'most HF called' cycles of a season must be at least 4 days apart (otherwise near-duplicates of one storm); 'none called' ranks by the number of day-2 points in the two top gust bands (72 kt and above). Both settled before the cases were looked at.
- Pull: the first run stopped at a stale `.idx` (offsets off by a few hundred bytes); `fetch_dmax.py` now falls back to walking GRIB message lengths, as `a_fetch.py` does. No data differences.
- Reported intervals: per-season ranges only (five seasons), as pre-registered.

## Addendum 2026-10-14 (second batch, Jason chose "Better replay")
Added a **routine sample** before looking at any of its results: every 10th 00 UTC cycle of the 1,061 (sorted order, offset 5), 106 cycles not already among the 20 picked cases (`select_sample.py`, `results/replay_sample.csv`). Nothing is chosen for interest, so coverage over this group is the fair read. Same shading, same leave-one-season-out tables (unchanged), same 324 nm coverage definition as the page. No new data pulled (saved daily-max grids reused). `build_payload.py` now writes 126 cases (1.25 MB).
