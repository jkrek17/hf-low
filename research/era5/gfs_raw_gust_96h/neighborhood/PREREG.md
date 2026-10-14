# Pre-registration: gust search radius around the model low (HF only)

Written 2026-10-14 before any score by radius was computed. Jason's ask (thread "Notes for the other site tabs"): measure from the model low; does the model storm reach the warning category (HF) within a radius of its low.

**Data**: `../results/raw_gust_fixes_96h.csv.gz` (180,784 GFS PROXY fixes, pipeline A detector on GFS, 5 seasons 2021-22 to 2025-26, 12-hourly leads 0-96 h). Truth Y0 unchanged: an OPC archive Category HF fix within +-3 h and 600 km of the model low. OPC is an analyst call, not a measurement; the archive holds low positions, not wind areas.

**Radii** (gust maximum over ocean within R of the model low): columns `g500`, `g800`, `g1200` = 500 / 800 / 1200 km = 270 / 432 / 648 nautical miles. Jason asked for about 500 nm (926 km); that radius was not stored and the GFS fields were not kept (31 GB), so 500 nm falls between the 432 and 648 nm results and is not scored directly.

**Scoring**: per basin and lead, as `score96.py`: held-out leave-one-season-out count-matched threshold (fitted separately for each radius at lead 0, applied as G > T), plus fixed 64 and 72 kt. Report POD ("caught"), agree = hits/(hits+false alarms), HSS.

**Decision rule (fixed now)**: for each basin and lead, the gust radius has "stopped falling off" at the smallest radius whose next-larger radius adds less than 0.02 HSS (matched threshold). Summarise as the count of basin-lead cells (of 2 x 9) where 270 nm, 432 nm or 648 nm is that radius.

**Not claimed**: a rise in score with radius is partly built in (a larger search area finds more gust), so it is read as how far the model's gust placement can be off, not as added skill. No intervals; fixes in a storm are not independent; one model with version changes. Gale and storm are out of scope (the archive has no gale truth).
