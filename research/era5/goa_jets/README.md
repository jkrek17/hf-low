# Gulf of Alaska barrier jets and Cook Inlet / Shelikof gap flow: when do they get strong?

**ERA5 PROXY throughout** (reanalysis 0.25 degree gust and wind, not observation; October-April, 2004-05 to 2025-26, 00 and 12 UTC). Plan and every
definition: `PREREGISTRATION.md`, committed first (`2f74f02`) before any wind direction or pressure was read. Everything quoted here is in `results/`.
Lows are this study's own ocean pressure minimum, not pipeline A or B. "Storm force" and "strong" below are ERA5 gust classes (54 kt and 64 kt), not warning categories.

## What a forecaster can use

**1. Southeasterly flow along the coast roughly doubles the chance of a storm-force gust.** Chance that a 12-hourly time has a 54 kt gust somewhere in the zone
(`results/analysis.txt`, Q1):

| Zone, flow (zone-mean 10 m wind direction) | Reaches 54 kt | Reaches 64 kt | Storm-force episodes per season | Share of episodes that reach 64 kt |
|---|---|---|---|---|
| North Gulf coast, southeasterly (from 100-170) | about 9 in 100 | 1 in 100 | 7.5 | 16% |
| North Gulf coast, onshore (from 170-260) | 5 in 100 | 0.5 in 100 | 3.3 | 15% |
| North Gulf coast, offshore (the rest) | 6 in 100 | 0.7 in 100 | 5.8 | 17% |
| Panhandle coast, southeasterly | about 10 in 100 | 2 in 100 | 7.1 | 23% |
| Panhandle coast, onshore | 5.5 in 100 | 0.7 in 100 | 4.6 | 15% |
| Panhandle coast, offshore | 4 in 100 | 0.5 in 100 | 4.3 | 18% |
| Cook Inlet / Shelikof, northerly (from 300 through north to 060) | about 8.5 in 100 | 0.7 in 100 | 5.0 | 12% |
| Cook Inlet / Shelikof, inflow (from 090-200) | 5.5 in 100 | 0.4 in 100 | 4.2 | 12% |
| Cook Inlet / Shelikof, other | 6 in 100 | 0.5 in 100 | 6.1 | 14% |

Episodes are storm-force times linked within 36 h. Hurricane force (71.7 kt proxy) is too rare to say anything about: 8, 13 and 5 episodes in 22 seasons
for the three zones.

**2. The flow direction opens the door; how deep and how near the low is decides how strong it gets.** Chance of a 54 kt gust (and in brackets of a 64 kt gust)
by the deepest low in the Gulf region (44-62N, 175-125W), within the southeasterly flow (`results/q1_table.csv`):

| Deepest low | North Gulf coast | Panhandle coast |
|---|---|---|
| Weaker than about 978 hPa | 3.5 in 100 (0.2) | 3 in 100 (0.2) |
| 967 to 978 hPa | 13 in 100 (0.6) | 14 in 100 (2) |
| 957 to 967 hPa | 19 in 100 (3) | 22 in 100 (6) |
| Deeper than 957 hPa | 20 in 100 (5) | 19 in 100 (4) |

Band edges are percentiles of all times in that flow, so they differ a little by zone (978.8/966.9/956.6 and 978.0/966.6/955.7 hPa). The top band holds only 9 to 10 times, so its
range is wide (13 to 41 in 100 for the north Gulf). In **northerly Cook Inlet / Shelikof flow depth matters much less**: 5, 12, 12 and 11 in 100 across the same four bands.

**3. A low within about 700 nm of the zone is worth about two to four times the chance** (`q1_table.csv`, LOW_DIST). In the southeasterly flow, 13 in 100 times reach 54 kt
with the low within about 680 nm (north Gulf) and 14 in 100 within about 720 nm (Panhandle), against 4 to 9 in 100 farther out; in northerly Cook Inlet flow 14 in 100 within
about 650 nm against 1 to 4 in 100 beyond. At 64 kt: 1.7, 3.1 and 1.4 in 100 inside against at most 0.7 outside.

**4. The pressure pattern across the coast or the range.** Top 5% of times in the same flow (`q1_table.csv`, noisy: 9 to 10 times):
- Pressure higher at the coast than 300 to 400 km offshore by 9.5 hPa or more, in southeasterly flow: about 27 in 100 (north Gulf, range 15 to 53) and 24 in 100 (Panhandle, range 12 to 80) reach 54 kt, against 7 to 12 in 100 in the lower bands. The same band in offshore flow is as high (33 in 100 north Gulf), so this measures a strong cross-coast pressure gradient in general, not barrier flow in particular.
- Cook Inlet / Shelikof, northerly flow: inland sea-level pressure (62-66N, 155-145W) higher than at the strait exit (56-59N, 156-150W) by 15.7 to 21.7 hPa: 16 in 100 times reach 54 kt (range 9 to 29); by 21.7 hPa or more: 30 in 100 (range 17 to 75); below 15.7 hPa: 4 to 8 in 100. Inflow flow shows the same (15 and 43 in 100) and "other" flow in its top band (31 in 100), so again it is the pressure gradient, not the wind direction. Inland pressure is sea-level-reduced and cold air inflates it.

**5. Cold air over the warm sea matters most for the Cook Inlet gap flow (a lead, not a confirmation).** In northerly flow, with the open Gulf (50-58N, 160-145W) sea surface
warmer than the 2 m air by 3.8 to 6.1 K: 22 in 100 times reach 54 kt (range 16 to 35), by 6.1 K or more: 29 in 100 (range 18 to 45), against 3 in 100 below 1.8 K and 8 in 100 between 1.8 and 3.8 K.
Among Cook Inlet gap episodes, strong (64 kt) against storm-force ones leans 3.7 times (range 1.6 to 9.2, q 0.010) per 2 K, but only 13 episodes are strong, so this is
thin. It stands up to adding a late-winter term (3.7), to dropping any one season (3.3 to 6.2), and to a within-season shuffle (no shuffle of 1,000 matched it), and points the same way in even and odd seasons (6.2 and 2.1, `results/audit.txt`). It is post hoc in the sense that it is one of 15 tests and
no held-out seasons remain; it needs new winters.

![chance by ingredient](results/fig_chance_by_ingredient.png)

## Strong versus storm force among episodes (Q2)

The pre-registered held-out test could not run (it needed 30 strong episodes in each half of the seasons; the biggest class has 18 and 18 against 52 and 69 in the other
group, and the Cook Inlet class has 6 and 7). The **all-seasons estimates are the pre-registered fallback and are not held out** (`results/q2_tests.csv`, odds ratio per standard deviation, strong
against 54 to 64 kt episodes, ridge-penalised logistic, adjusted for depth, season-block bootstrap, Benjamini-Hochberg over the 15):

| Factor (one SD) | North Gulf southeasterly (27 of 166 strong) | Panhandle southeasterly (36 of 157) | Cook Inlet northerly (13 of 109) |
|---|---|---|---|
| Deeper low (11, 10, 13 hPa) | **2.1 [1.4, 3.7]** | 1.5 [1.0, 2.3] (q 0.057) | **1.35 [1.07, 1.8]** |
| Low farther away (620, 770, 450 km) | 0.84 [0.56, 1.3] | **0.60 [0.33, 0.90]** | 0.72 [0.52, 0.97] (q 0.054) |
| Coast pile-up (4.4, 4.5 hPa) or inland pressure (8.8 hPa) | 1.06 [0.73, 1.5] | 1.38 [1.01, 1.9] (q 0.064) | 1.16 [0.59, 2.7] |
| Sea warmer than air (1.6, 1.9, 2.2 K) | 1.19 [0.89, 1.9] | 1.26 [0.92, 1.7] | **3.7 [1.6, 9.2]** |
| Pressure gradient at the gust maximum | **2.5** | **4.0** | **2.6** |

**Bold** = Benjamini-Hochberg q below 0.05. The gradient row is the wind itself seen as pressure (the same gust field defines where it is measured), so it is not a clue;
it is listed because it was pre-registered. Depth and distance describe the storm at the same time as the gust; they are not an advance warning.
**Cannot tell** (intervals include both sides, power for effects below about 1.7 to 2.3 per SD is low): coast pile-up in the north Gulf, inland pressure in Cook Inlet, sea-air temperature in
either barrier zone.

## Does ERA5 under-read? (Q3, buoys, `results/buoy_*.csv`)

Six open-water NDBC buoys (46001, 46061, 46076, 46080, 46082, 46083; 5 m anemometers, no height adjustment, positions approximate) at the pulled times:
- ERA5 10 m wind is about equal to the buoy wind (median ratio 0.95 to 1.04 when the buoy reads 20 kt or more).
- ERA5 gust is **8 to 17% higher** than the buoy gust (median ratio 1.08 to 1.17). When ERA5's cell shows a 54 kt gust, the buoy's median gust was **44 to 48 kt**, and 0 to 9% of those times had a buoy gust of 54 kt or more. So the 54 kt proxy class reads like a high gale on these buoys, in line with the ASCAT result (ERA5 15 to 25% above).
- **Shelikof Strait (46077) sits in an ERA5 land cell**: ERA5 wind is 0.27 of the buoy's at 30 kt or more (median). ERA5 cannot see the strait or the inlet, so the Cook Inlet / Shelikof numbers above describe the open-ocean zone south and east of them (ocean cells only), not channel winds. Real gap winds there are underestimated by an unknown, large amount.
- This says nothing about the cores of barrier jets within a few tens of nm of the mountains, where no buoy sits.

## What was not done

Front position (a "front pinching" is only represented by the flow, the low's distance and the cross-coast pressure), low-level stability above the surface (no pressure-level temperature pulled), ASCAT, 1979-2004, hurricane force (too rare), channel winds in Cook Inlet or Shelikof, teleconnections.

## Looks and deviations

All of Q2 is one look at all 22 seasons (no held-out data). Q1 bands use 600 random controls; the top bands hold 9 to 13 controls each, so their ranges are wide and the bootstrap range can exceed 100%
(three bands, in the Cook Inlet inflow rows). Post hoc: `audit.py` (month term, halves, leave-one-season-out, placebo); the buoy comparison (stations from memory, approximate). A bug fix before any
result was read: the sea-air temperature used a plain mean that returned missing over cells without sea surface temperature; it now uses the mean of valid cells (`features.py`).

## Files

`PREREGISTRATION.md`, `stage0_counts.py` (gust counts), `stage1_pull.py` (22.5 GB, 1,800 times), `features.py`, `analyse.py`, `audit.py`, `buoy.py`, `figures.py`. Results in `results/`:
`q1_table.csv`, `q2_tests.csv`, `audit.csv`, `episodes.csv`, `features.csv`, `buoy_summary.csv`, `buoy_when_era5_54.csv`, `analysis.txt`, `audit.txt`. The pulled window of fields
(44-70N, 175-125W, u10, v10, MSLP, 2 m temperature, SST at the 1,800 times) is kept at `/mnt/project-files/goa-jets/` so nothing needs pulling again.
