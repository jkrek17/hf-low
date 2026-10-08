# Seasonal cycle of hurricane-force lows, and the Pacific midwinter question (pre-registration)

Written 2026-10-08, before any monthly count, share, latitude or eddy-activity number in this
directory was computed. Agenda item 4 (`/mnt/project-files/climatology-agenda/hf-climatology-agenda.md`),
ledger item 5 under "Future tasks". Anything changed after the first computation is logged in
`README.md` under "Departures from this plan".

**What has already been seen** (so the reader can judge how open the test is):
- The agenda's archive shares of first HF fixes by month: Pacific Nov 13.2%, Dec 20.6%, Jan 22.1%, Feb 14.6%;
  Atlantic Dec 18.5%, Jan 20.5%, Feb 19.2%. These are shares of counts, not per-day rates, and not
  uncertainty-tested. So the archive's January peak is a retrospective hypothesis check, not a prediction.
- The map thread's month-of-year view appeared to show the Pacific belt sliding south in January. A visual
  hint only.
- Nothing from pipeline A by month, nothing about eddy activity, and none of the decompositions below.

Everything ERA5 here is **pipeline A** (`research/era5/hf_history`, 800 km ocean gust index, HF-equivalent at
71.7 kt) and a **proxy**; "archive" is the OPC hurricane-force CSVs on `main`. Eddy fields are ERA5 reanalysis
(a different use of ERA5, not the gust index).

## Question

The archive's Pacific HF count peaks in January. Nakamura (1992, *J. Atmos. Sci.*) found that Pacific
baroclinic eddy activity is suppressed in midwinter, when the jet is strongest. Is the HF-low seasonal cycle
inconsistent with that, and if the two really differ, which of three things reconciles them: (R1) fewer
cyclones but a higher HF share in midwinter, (R2) storms forming and peaking farther south, (R3) storm-track
eddy activity and extreme events are different things.

## Literature prediction, per month

Nakamura (1992) and the work citing it (e.g. Afargan and Kaspi 2017): in the Pacific, 2-6 day eddy activity
(v'^2 at 250 hPa, v'T' at 850 hPa) is high in Nov and in Mar-Apr, and lower through midwinter (Dec-Feb,
lowest Jan-Feb) although baroclinicity and the jet are strongest then. In the Atlantic there is no
suppression: eddy activity peaks Dec-Feb. (From memory, not re-read for this plan; the sign of the contrast
below is the prediction, not any magnitude.) Expected signs of the two contrasts below (D for Jan, M for
midwinter against shoulders):

| quantity | Pacific eddy activity | Atlantic eddy activity |
|---|---|---|
| D_Jan (Jan vs mean of Dec, Feb) | < 0 | >= 0 |
| M (Dec-Feb vs mean of Nov, Mar) | < 0 | > 0 |

If HF lows simply followed eddy activity the same signs would hold for Pacific HF counts. The archive's
January share says they do not (seen, above), so the competing explanations below are what the test is for.

## Definitions

- Season: 1 June to 31 May, labelled by its starting year. Months compared: Nov, Dec, Jan, Feb, Mar (the
  contrasts); all twelve months are tabulated for description.
- Tracks: all 75,087 pipeline A tracks in `research/era5/hf_history/results/all_tracks.csv.gz` (a low below
  1010 hPa for at least 24 h, with at least two fixes in a basin domain), basin as given. HF = `gust800_kt >= 71.7`.
  The domain is fixed at 20-75N; storms crossing its edge are not handled (as in `freq_split`).
- Genesis month: month of the track's first fix (`start`). Primary assignment of every pipeline A count.
  A cohort is "all tracks starting in month m", so HF count = cyclone count x HF share exactly.
- HF onset month (secondary, matches the archive's convention): month of the first 6-hourly fix with
  `g800 >= 71.7` in `era5_hf_catalog_tracks.csv`. Onset latitude is that fix's latitude.
- Archive: month of an event's first fix with Category `HF` (DHF, the forecast category, does not count),
  `data/hf_lows/HF_Data_-_{Atl,Pac}.csv`, seasons 2004-05 to 2025-26. Latitude is that fix's latitude.
- Real cyclones: the deep stage of `freq_split` (`minp <= 1000 hPa and n_fix >= 8`, or HF).
- Depth-strong: `minp <=` the count-matched cut already fixed by `freq_split` (`results/meta.json`: Atlantic
  966.2 hPa, Pacific 965.0 hPa; matched to the HF count over Oct-Apr 2004-05 to 2025-26). Pressure depth did
  not drift before 2001, so this series uses all 47 seasons (1979-80 to 2025-26).
- Rates: per day, sum of counts over seasons divided by the summed days of that calendar month in the same
  seasons (leap Februaries included).
- Gust-based quantities use 2004-05 onward only (decision 1). The 1979-2000 gust-based shape is shown
  descriptively and flagged "within-era only, not tested".

## Eddy activity (a small ERA5 pull; well under 50 GB)

Source: ARCO-ERA5, 1.5 degree, 6-hourly, `1959-2022-6h-240x121_equiangular_with_poles_conservative.zarr`
(it ends 2021-12-31, so the last complete season is 2020-21). Fields:
- v-wind at 250 hPa (the Nakamura measure), seasons 2004-05 to 2020-21 (17 seasons), 20-75N.
  Chunks hold 8 times x 13 levels (about 10.5 MB each, measured on two chunks); about 3,220 chunks, **estimate
  about 34 GB**, no cache covers any of it. The full 1979 start would be about 82 GB and is not pulled.
- Mean sea level pressure, seasons 1979-80 to 2020-21 (42 seasons), 20-75N; about 0.55 MB per chunk, about 4.2 GB.
  This is a cheaper surface-based measure and is a secondary quantity, not a substitute for v250.

Band-pass: 4th-order Butterworth, zero-phase, 2-6 days (passband 1/6 to 1/2 cycle per day), applied to the
continuous 6-hourly series at each grid point. Eddy activity = mean of the squared band-passed field, per
calendar month and season, cos(lat)-weighted mean over a fixed box:
- Pacific 30-60N, 150E-150W. Atlantic 30-60N, 70W-10W.
Descriptive only: the latitude of the sector-mean maximum (storm-track axis) by month.

## Estimators

For a monthly quantity x, `D_Jan = ln x(Jan) - 0.5 [ln x(Dec) + ln x(Feb)]` and
`M = mean(ln x(Dec), ln x(Jan), ln x(Feb)) - mean(ln x(Nov), ln x(Mar))`. For a count, x is the pooled
per-day rate; for the share, HF count / cyclone count; for eddy activity, the season-mean of monthly means.
For latitudes the contrast is in degrees without logs (x = pooled mean latitude).
D and M are additive, so `D_HF = D_cyclones + D_share` and the share's fraction of the effect is
`f = D_share / D_HF` (only interpreted when D_HF is clearly non-zero).

Uncertainty: 95% percentile interval and two-sided p from a 2,000-draw season-block bootstrap (seasons
resampled with replacement, all months and all quantities in a comparison resampled jointly), seed 20261009.
Leave-one-season-out range for each primary D. Detectable effect at 80% power ~ 2.8 x bootstrap SE.
Effective n is the number of seasons (22 gust-based, 47 depth-based, 17 or 42 eddy), not the number of storms.

## Tests (the whole family, fixed now)

Quantities, each for both basins, each with D_Jan and M:

| id | quantity | seasons |
|---|---|---|
| Q1 | archive HF count by first-HF-fix month | 2004-05 to 2025-26 |
| Q2 | pipeline A HF count by genesis month | same |
| Q3 | pipeline A all-cyclone count by genesis month | same |
| Q4 | HF share = Q2 / Q3 | same |
| Q5 | pipeline A real-cyclone count (genesis month) | same |
| Q6 | pipeline A HF count by HF onset month | same |
| Q7 | pipeline A depth-strong count (genesis month) | 1979-80 to 2025-26 |
| Q8 | depth-strong share = Q7 / Q3 | same |
| Q9 | v250 eddy activity, box | 2004-05 to 2020-21 |
| Q10 | MSLP eddy activity, box | 1979-80 to 2020-21 |
| Q11 | archive mean latitude of the first HF fix | 2004-05 to 2025-26 |
| Q12 | pipeline A mean HF onset latitude | same |
| Q13 | pipeline A all-track mean `peak_lat` | same |
| Q14 | paired contrast `D(Q2) - D(Q9)` (and the same for M), on the 17 seasons both cover | 2004-05 to 2020-21 |

14 quantities x 2 basins x 2 contrasts = 56 tests; Benjamini-Hochberg q is across all 56. The primary set,
Pacific D_Jan for Q1, Q2, Q3, Q4, Q9, Q11 and Q14 (7 tests), also gets Bonferroni x7. Everything else is
secondary and labelled so. Not in the family (descriptive, no test): the monthly tables, the 1979-2000
gust-based shape, the storm-track axis latitude, season leave-one-out ranges.

## Hypotheses and decision rules

- **H0-peak (is there a January peak at all?).** Pacific D_Jan for Q1 and Q2. "Peak established" if both
  intervals lie above 0; "no peak" if both include 0; otherwise "archive only" or "proxy only" and said so.
  If the peak is not established, the reconciliation questions are not interpreted, and the result is reported
  with its power ("a dip or peak smaller than X% cannot be excluded").
- **H-eddy (does our reanalysis reproduce the lull?).** Pacific Q9 D_Jan < 0 and M < 0 (interval below 0).
  If it does not, the comparison with Nakamura is not a fair test of the HF cycle and this is the first
  thing reported.
- **R1 (cyclone count dips, HF share rises).** Supported if, given a peak in Q2, D(Q4) has an interval above 0
  and f > 0.5. Also read with Q5, Q8.
- **R2 (storms farther south in midwinter).** Supported if D(Q11) (archive) and D(Q12) (pipeline A) both have
  intervals below 0 in the Pacific. One of the two only: "partial".
- **R3 (eddy activity and extremes differ).** Supported if H-eddy holds and D(Q14) has an interval above 0
  (the HF cycle departs upward from the eddy cycle on the same seasons).
- These can hold together; each is reported separately. If none is supported but a peak is established,
  the answer is "unexplained by these three".
- The Atlantic is the comparison: eddy activity is predicted not to dip, so a Pacific-only dip in Q9 with a
  peak in HF counts is evidence that the Pacific question is a specific one.

## Pitfalls named in advance

- Months differ in length; every count is a per-day rate.
- A monthly share of the season total is not a rate and is not used for inference.
- The agenda notes eddy activity is not storm count; Q9 and Q10 are amplitudes, Q1-Q8 are counts, which is
  why R3 is a separate hypothesis.
- Month assignment: genesis month (cohort) and onset month (archive convention) differ for storms born late in
  a month; both are reported.
- ERA5 1.5 degree data damp small-scale variance; the seasonal contrast of band-passed variance is the
  quantity, not its level.
- Archive detection (scatterometer coverage, forecaster practice) could depend on the month; this cannot be
  tested and is a caveat on Q1, Q11.

## Not in scope

New cyclone tracking; jet or Eady-growth diagnostics (v250 only is pulled); trends; interannual links to
teleconnections; the full 1979 v250 record (over the gate's comfortable range and of little added value).
