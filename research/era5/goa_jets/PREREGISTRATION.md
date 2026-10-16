# Gulf of Alaska barrier jets and Cook Inlet / Shelikof gap flow: when do they get strong?

Written and committed BEFORE any wind direction, pressure or other ingredient was looked at for any time. Only
`stage0_counts.py` has run: it reads the cached ERA5 gust (a PROXY, not observation) and reports the largest ocean gust in
each zone at each time. Jason asked 2026-10-16 09:04Z ("barrier jets ... fronts pinching up against the coast, and channelled flow
... towards Kodiak; no real feel for when the effects will be significant"), and on 09:54Z chose Cook Inlet and Shelikof gap flow
off the Alaska and Aleutian Range for the second half; "Go" 11:46Z and 12:10Z (defaults: no 1979-2004 extra, buoy check yes).

## Zones (ERA5 0.25 degree OCEAN cells, land-sea mask < 0.5, inside the box)

| Zone | Box | Meant to catch |
|---|---|---|
| Z1 North Gulf coast | 56-61N, 150-138W | barrier flow along the Chugach and St. Elias ranges, Prince William Sound to Yakutat |
| Z2 Panhandle coast | 53-58N, 140-130W | the same barrier further south, Cape Spencer to Dixon Entrance |
| Z3 Cook Inlet / Shelikof / Kodiak | 55-60N, 158-148W | gap and channelled flow off the Alaska and Aleutian Range |

Zone gust = largest ERA5 10 m gust over the zone's ocean cells at a 00 or 12 UTC time, October to April, seasons 2004-05 to 2025-26
(decision 1: no earlier seasons; 9,338 times). **ERA5 0.25 degree smooths narrow straits and steep coast: gap winds in Shelikof Strait and
Cook Inlet and the core of a barrier jet within a few tens of nautical miles of the mountains will read LOW, and coastal
cells mix land into the gust. Every result is about the proxy, not about what a ship or buoy at a specific point felt.**

## Definitions (fixed now)

- **Event time**: zone gust at or above 54 kt (nominal proxy storm class, not calibrated on observed events). **Strong**: zone gust at or above
  64 kt; **hurricane-force proxy** (71.7 kt) is counted but NOT tested, because stage 0 shows 5 to 13 times per zone in 22 seasons.
  ERA5 71.7 kt is about 62 kt in ASCAT (ascat-check); read every class as "ERA5 proxy class".
- **Episode**: one zone's event times linked when 36 h or less apart, represented by the peak-gust time; strong if any time is strong.
  **The unit of analysis is the episode; the independent sample for inference is the season (22), not times or episodes.**
- **Regime** (assigned to EVERY time from the zone-mean vector 10 m wind direction over the zone's ocean cells, FROM direction, so it exists
  for non-event times too; this replaces any use of the direction at the gust maximum):
  - Z1, Z2 (coast faces southwest; the long axis runs 135/315): **BARRIER** = wind from 100-170 (southeasterly, along the coast, mountains
    on the right-hand side); **ONSHORE** = from 170-260; **OFFSHORE** = the rest (from 260 through north to 100: outflow from the ice fields and valleys).
  - Z3: **GAP** = from 300 through north to 060 (cold-air outflow down Cook Inlet and Shelikof); **INFLOW** = from 090-200 (ahead of a low);
    **OTHER** = the rest.
- **Reference low**: the minimum MSLP over OCEAN cells in 44-62N, 175-125W (so no sea-level-reduction low over mountains can be chosen).

## Ingredients (all standardised, computed at the peak time of an episode or at the time itself; none chosen after outcomes)

1. DEPTH: reference-low MSLP, sign flipped (deeper = higher).
2. LOW_DIST: great-circle distance from the reference low to the zone's gust maximum (km).
3. GRAD: mean MSLP gradient magnitude (hPa per 100 km) over ocean cells within 200 km of the zone's gust maximum.
4. COASTPILE (Z1, Z2 only): mean MSLP over ocean cells within 100 km of land minus the mean over ocean cells 300-400 km from land, within the
   zone's longitudes. Higher means pressure piled up against the coast (the barrier-jet mechanism).
5. INTERIOR_DP (Z3 only): mean MSLP in the interior box 62-66N, 155-145W minus the mean MSLP of ocean cells in 56-59N, 156-150W (the Shelikof
   and Cook Inlet exit). Higher means interior high pressure over the exit low (the push through the gap). Over land this is sea-level-reduced pressure;
   cold interior air can inflate it, which is a stated caveat, so it is reported beside INTERIOR_T2M.
6. INTERIOR_T2M (Z3 only, reported with 5, not as a separate test): mean 2 m temperature in the interior box (negative sign: colder = higher).
7. STAB: SST minus 2 m temperature (K), mean over ocean in the open Gulf box 50-58N, 160-145W, the same box for every zone. This is a surface-based
   stability proxy only; no pressure-level temperature profile is pulled, so true low-level stability and the height of the blocking layer are NOT tested.
8. WSPD: zone-mean 10 m wind speed (kt), described only (it defines the regime, so it is not a test factor).

Fronts are not diagnosed. "A front pinching against the coast" is represented only by the regime plus LOW_DIST and COASTPILE, and the write-up says so.

## Data (pull estimate, stated before pulling)

ARCO-ERA5 0.25 degree u10, v10, MSLP, 2 m temperature and SST at: (a) every event time in any zone (union), and (b) **600 CONTROL times drawn
uniformly at random from all 9,338 October-April times** (seed 16), kept whether or not they are events, so the controls are an unbiased sample of the
climate and give the base rates. Streaming cost is whole global fields: about 13.5 MB per time (5 variables x 2.7 MB, rate of the tip-jet pull).
The union count is written below after stage 0, before pulling. Window kept: 44-70N, 175-125W. No existing cache holds these variables at these times
(the 4.1 GB cache is gust only). WeatherBench2 upper air is not used.

## Questions and tests

**Q1, how often (descriptive, the forecaster table).** For each zone and regime: how many episodes there are per season, the share that reach strong,
the share of ALL times in the regime that are event times, and, using the controls for the ingredient distribution, the chance that a time is an event
(and strong) given the ingredient bin. Bins are fixed on the CONTROL distribution of each ingredient in the regime: the bottom half, 50th-80th, 80th-95th,
and top 5% percentile bands. P(event | bin) = events in bin / (9,338 x control share of the bin). Intervals by season-block bootstrap. No test, no selection:
every ingredient is shown.

**Q2, strong versus storm (tested).** Among episodes in one zone-regime, which ingredients are higher at strong than at 54-64 kt peaks? Classes
tested: Z1 BARRIER, Z2 BARRIER, Z3 GAP. Factors: DEPTH, LOW_DIST, GRAD, STAB, plus COASTPILE (Z1, Z2) or INTERIOR_DP (Z3): 5 per class, family of 15.
Univariate ridge-penalised (0.5) logistic regression of strong on one standardised factor, adjusting for DEPTH when the factor is not DEPTH, season-block
bootstrap (2,000 resamples), BH over the 15. A factor is **supported** if q < 0.05 in the fit half (even or odd seasons) and the same sign with
one-sided p < 0.10 in the held-out half in at least one direction; **replicated** if so in both. A class is tested only if each half holds at least 30 episodes of
the rarer outcome; otherwise it is described only. Minimum detectable odds ratio per SD at 80% power is reported either way. A well-powered null is a result.
Because the 64 kt class may be too rare for the half-split rule, the **all-seasons estimates with season-block intervals are pre-registered as the
fallback report** and labelled "not held out"; they are leads, not confirmations.

**Q3, does the proxy under-read (descriptive).** At event times, ERA5 10 m wind at the grid cell nearest each available buoy in or near the zones
(NDBC standard meteorological history, `ndbc.noaa.gov/data/historical/stdmet`, hourly wind speed and gust, taken at 00 and 12 UTC; candidates 46001, 46061, 46080, 46076, 46082, 46083, 46077, whichever have data) against the buoy's reported wind,
reported as the ratio ERA5/buoy by regime. Buoy winds are at anemometer height, ERA5 at 10 m; no height adjustment; the ratio is a bias indicator, not a correction.
If the buoys are not reachable this is reported as not done.

Not tested: NAO or other teleconnections, anything before 2004-05, upper-air or low-level stability, front position, land-fast ice effects, individual
fjord or channel winds.

## Reporting rules

Name the pipeline: gust is the shared cache proxy; lows come from this study's own ocean-MSLP minimum, not pipeline A or B. Say proxy wherever it appears.
Nautical miles in text (km in tables). Plain rates for forecasters ("about N in 10"). Every quoted number recomputed by a fresh Sonnet agent; any
surprising positive gets hf-surprise-audit. Figures in plain-words text must come from a result file (PR 244). Any later change is logged below.

## Stage 0 counts (written before the pull)

`results/stage0_out.txt`: 9,338 October-April times. Zone gust at or above 38 / 50 / 54 / 64 / 71.7 kt: Z1 3,999 / 1,132 / 625 / 68 / 9;
Z2 3,659 / 1,015 / 582 / 87 / 13; Z3 4,382 / 1,191 / 625 / 51 / 5. **Union of the three zones at 54 kt or more: 1,282 times.** With 600 random controls
(about 12% of them also events) the pull is about 1,810 times x 13.5 MB = **about 24 GB, under the 50 GB gate**. Strong (64 kt) times are 51 to 87 per zone, so
episodes will be few dozen per zone and the half-split rule is likely to fail; the all-seasons fallback above is the planned report.

## Post hoc deviations

- 2026-10-16, before any result was read: `features.py` took a plain mean of SST minus 2 m temperature that came out missing (SST is missing at some coastal cells); changed to the mean of valid cells. No outcome had been read.
- 2026-10-16 (POST HOC): the pre-registered half-split tests did not run (strong episodes per half: 11/16 North Gulf, 18/18 Panhandle, 6/7 Cook Inlet; 30 needed in the rarer group of each half). The pre-registered all-seasons fallback is reported, not held out.
- 2026-10-16 (POST HOC): `audit.py` (late-winter term, even/odd halves, leave-one-season-out, within-season placebo) was added after the Q2 results, because the Cook Inlet sea-air temperature result was the largest and rests on 13 strong episodes.
- 2026-10-16: GRAD (pressure gradient near the gust maximum) is pre-registered and reported but flagged as the wind itself; not a clue.
- 2026-10-16: Q3 buoy stations were first put in from memory; the claim reviewer found 46077 about 160 km from its NDBC position, which put it in a land cell and produced a false "ERA5 cannot see the strait" statement. Positions now match the NDBC station pages (`buoy.py check`), the buoy results were rerun, and the README was corrected. Logged in the failure library.
