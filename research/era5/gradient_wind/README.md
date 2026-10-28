# Pressure gradient to wind speed, western Atlantic (ERA5 proxy, gust)

Plan: `PREREGISTRATION.md` (committed first, `ea5e057`, before any gust was matched to a gradient). Asked by Jason 2026-10-28; options A (table for all cells) and C (storm-core cut). Code: `build_pairs.py` (match), `analyse.py` (tables and tests), `ascat_s1.py` (secondary check). Results: `results/gradient_wind_tables.txt`, `results/ascat_s1.txt`. No download; 21 cold-season winters (October to April, 2004-05 to 2024-25), 641,808 cell-times at 00 and 12 UTC, 30-45N, 78-55W, ocean at least 1 degree from land, sampled every 2 degrees (every 8th grid cell). A "chance" below is per grid point and time, not the chance that somewhere in an area reaches the class; the chance for an area maximum is higher.

**Answer in plain words.** In the western Atlantic cold season the pressure gradient gives a usable chance table for ERA5 *gust* at a grid point (below), and the same pattern shows in ASCAT mean wind up to its ceiling near 45 kt. It is a lookup and a sanity check against an independent pressure measurement. It has **not been shown to improve the GFS wind forecast**: the one gain found (15%) used the true analysis pressure and is not forecast skill. The wind-to-gradient ratio falls as the gradient tightens, and the isallobaric (pressure-tendency) wind is a real part of that, but small as a forecast term.

**Everything is the ERA5 proxy.** The wind is ERA5 instantaneous 10 m **gust**, which reads 15 to 25% above ASCAT and runs high in unstable air. Gale, storm and hurricane-force are nominal gust classes (38, 54, 71.7 kt). ERA5 assimilates ASCAT and buoys. Neither gust pipeline is used (the cells are chosen from the pressure field alone). Western Atlantic cold season only: nothing here is tested for the Pacific.

## How to use it

Measure the gap between isobars near the area you care about and turn it into hPa per 100 nm (4 hPa isobars 60 nm apart is 6.7 hPa per 100 nm). The model gradient here is taken over about 1 degree (60 nm north-south, about 45 nm east-west at 40N), so measure over a similar span, not across a single tight pair of isobars. A gradient read off a sharper analysis or model will usually be larger than the smoothed ERA5 value and put you in a higher row than this table's.

**Rule of thumb check.** The textbook value for 4 hPa per 60 nm at 45N is a geostrophic wind of about 54 kt. In the western Atlantic cold season the gust at that gradient (the 6-8 hPa per 100 nm row) is about 42 kt at the median (37.5-45N), about 0.8 of 54 kt (the median of the per-cell ratios of gust to geostrophic wind, which differs because the latitude varies, is 0.72). The ASCAT wind in the same gradient range is about 34 kt, about 0.63 of the textbook 54 kt (storm-neighbourhood cells, 30-45N; the median per-cell ratio is 0.56; see the satellite section). So the geostrophic value overshoots, and more so at tight gradients: the ratio of gust to geostrophic wind falls from about 0.8 at 4-6 to about 0.3 to 0.45 at 12 and over, partly because the gust does not go to zero with the gradient (the ratio is 1.4 to 1.7 under 2) and partly because ERA5 smooths the tightest gradients.

### Gust against gradient, all cells (option A), 37.5-45N and 30-37.5N

Median ERA5 gust, and how often the gust reached each class (winter-resampled 95% range in brackets). "thin" means fewer than 200 cell-times or fewer than 8 winters.

| Gradient (hPa per 100 nm) | 37.5-45N median gust | gale (38 kt gust) | storm (54 kt gust) | hurricane force (71.7 kt gust) | 30-37.5N median gust | gale | storm | hurricane force |
|---|---|---|---|---|---|---|---|---|
| under 2 | 17 kt | 0.5% | 0 | 0 | 17 kt | 0.1% | 0 | 0 |
| 2-4 | 27 kt | 9.4% (8.8-10.1) | 0.1% | 0 | 28.5 kt | 8.6% (7.9-9.3) | 0 | 0 |
| 4-6 | 36 kt | 40.7% (39.2-42.3) | 1.8% (1.4-2.1) | 0 | 38.5 kt | 55.1% (53.6-56.8) | 1.6% (1.1-2.3) | 0 |
| 6-8 | 42.5 kt | 75.1% (73.6-76.4) | 7.6% (6.9-8.5) | 0.0% | 46 kt | 88.1% (86.5-89.6) | 12.6% (9.6-15.9) | 0 |
| 8-10 | 47.5 kt | 88.2% (86.9-89.4) | 19.6% (17.4-22.1) | 0.2% (0.1-0.4) | 49.8 kt | 91.1% (88.0-94.4) | 28.9% (23.1-34.5) | 0.5% (0.0-1.1) |
| 10-12 | 50.5 kt | 93.2% (91.4-94.8) | 36.1% (32.2-40.2) | 0.2% (0.0-0.5) | thin | thin | thin | thin |
| 12 and over | 56 kt | 95.3% (93.8-96.7) | 60.1% (56.1-64.1) | 5.5% (3.8-7.5) | 55.5 kt | 95.8% (93.0-98.0) | 60.3% (50.9-68.9) | 9.6% (4.6-16.5) |

In plain terms at 37.5-45N: about 4 in 10 times at 4-6 hPa per 100 nm reach a 38 kt gust; about 6 in 10 times at 12 and over reach a 54 kt gust; a 71.7 kt gust happens about 1 time in 20 at 12 and over and almost never below 10 (the sample is 21 winters, so the hurricane-force rows rest on few events; the numbers of cell-times are in the results file). Medians at the same gradient differ by up to 3.5 kt between the two latitude bands, but the chances differ by up to 15 points (38 kt gust at 4-6: 41% at 37.5-45N, 55% at 30-37.5N; 54 kt gust at 8-10: 20% against 29%): use the band you are in.

**These are ERA5 gust thresholds, not warning categories.** Warnings are issued on sustained wind; a 38 kt gust is roughly a 30 kt sustained wind. The only sustained-wind evidence here is the ASCAT table, which stops near 45 kt. Do not use the table in the Pacific, north of 45N, east of 55W, or in summer.

Rows (cell-times) at 37.5-45N by gradient bin: 95,422 / 103,431 / 36,399 / 9,906 / 2,771 / 902 / 761; at 30-37.5N: 254,808 / 116,708 / 17,483 / 2,378 / 436 / 164 (thin) / 239. **Hurricane-force counts:** 37.5-45N, 12 and over: 42 cell-times on 31 days in 18 winters; 8-10: 6; 10-12: 2. 30-37.5N, 12 and over: 23 cell-times on 14 days in 10 winters; 8-10: 2. The intervals show winter-to-winter sampling only, not ERA5 bias.

### Storm core (option C)

"Core" is a cell within 216 nm of a surface low below 1000 hPa found in the pressure field itself. Gust to geostrophic ratio, **non-core minus core within the same gradient bins** (post hoc, see below): +0.045 (+0.035 to +0.054), 21 winters. At 0.6 to 0.7 that is about 7%: for the same gradient, the gust in a storm core is a little lower than outside it. The tables are close (37.5-45N, 8-10: core 46.5 kt median, 18.8% reaching storm gust; all cells 47.5 kt, 19.6%) because core cells are 22% of cell-times at 6-8, 44% at 8-10, 64% at 10-12 and 87% at 12 and over. Post hoc, one band and one bin: at 30-37.5N and 8-10, non-core cells (n=207) reach a 54 kt gust 37.7% (29.4-44.9) of the time against 21.0% (15.3-26.9) in cores (n=229). It does not repeat at 37.5-45N and 8-10 (non-core 20.3% against core 18.8%), so treat it as a lead.

### Stability

Classes of sea-minus-air difference are below 0 C, 0 to below 6 C, and 6 C or more; 12 winters (2013-14 to 2024-25) have temperatures. At 37.5-45N, median gust at 2-4 / 4-6 / 6-8 hPa per 100 nm: stable air 23 / 33 / 42 kt, weakly unstable 26.5 / 36 / 42.5 kt, strongly unstable (water 6 C or more warmer) 30.5 / 40.5 / 46 kt. The chance of a 38 kt gust at 4-6 hPa per 100 nm: 25% stable, 41% weakly unstable, 62% strongly unstable. The registered ratio test at 6-10 (P1 below) found little difference, but the chance of a 54 kt gust at 6-8 is still about 2.5 times higher in strongly unstable air (16.5% against 6.3% in stable air, 37.5-45N; descriptive). The larger effect at weak and moderate gradients is the post hoc P1c below. ERA5's gust adds a turbulent part in unstable air, so part of this is the gust, not the mean wind.

## Planned tests

| Test | Result | Verdict (plan's words) |
|---|---|---|
| P1 strongly unstable minus stable, gust/geostrophic, gradient 6-10 | +0.017 (-0.001 to +0.034), 12 winters, q 0.070 | by the plan's rule "no (well powered)", reading "plus or minus 5%" as plus or minus 0.035 (my reading, 5% of a typical ratio of 0.7): the interval lies inside it, but only just (upper end 0.034), q is 0.07, and the post hoc P1b upper end (0.044) is above it, so this "no" is fragile |
| P2 non-core minus core, gust/geostrophic, gradient 6 and over | +0.107 (+0.096 to +0.120), 21 winters, q below 0.001 | "yes" by the plan's rule (interval excludes 0 in the expected direction), but the pooled statistic is mostly the mix of gradient bins (core cells sit at tighter gradients, where the ratio is smaller); see P2b |
| P3 curvature: Spearman(gust, gradient wind) minus Spearman(gust, geostrophic), core, gradient 4 and over | -0.058 (-0.078 to -0.039), 9,472 cell-times, 21 winters, q 0.0008 | no verdict word fits: the interval lies wholly on the opposite side (worse), beyond the plus or minus 0.03 "no" band. The curvature-corrected wind tracks the gust less well than the plain geostrophic wind |

- P3 reasons (not tested): the isobar curvature is computed from ERA5 pressure smoothed over 1 degree and its error is large where the gradient is small, and 3.9% of core cells with gradient 4 and over have no real gradient-wind root and are dropped (6.3% of non-core cells). The verifier found 93% of the no-root cells have a negative discriminant, that is, strongly anticyclonic isobar curvature, not grid-edge failures. **The curvature correction is not recommended for a forecaster rule on this evidence.**
- Multiplicity: the three planned tests are corrected together (Benjamini-Hochberg); the descriptive tables are not tests.

## Post hoc (added after the pooled P2 was read; labelled)

- P2b: non-core minus core within gradient bins 6-8, 8-10, 10-12, 12 and over, averaged: +0.045 (+0.035 to +0.054), 21 winters. This replaces P2 as the number to quote for the core effect, and P2 is kept as registered.
- P1b: strongly unstable minus stable averaged over bins 6-8 and 8-10: about +0.01 to +0.02 depending on the minimum points per winter and bin (the verifier got +0.012 to +0.018; my first run gave +0.022), interval includes 0, 12 winters. P1c (bins 2-4 and 4-6, where the table shows the biggest differences): +0.245 (+0.224 to +0.266). The registered range (6-10) was chosen for the gradients where the Gulf Stream work looked; the bigger effect at weak gradients was not registered.

## Satellite check (secondary, ASCAT, labelled; `results/ascat_s1.txt`)

48,810 ASCAT cells (19 winters, 2006-07 to 2024-25) within 1.5 h of 00 and 12 UTC, no rain flag, ocean at least 1 degree from land. They lie within 800 km of a tracked low, so they describe the neighbourhood of lows, not every cell. ASCAT is equivalent-neutral wind, a mean wind, not a gust, and saturates near 45 to 50 kt, so rows from 8 hPa per 100 nm up are lower bounds. The plan restricted the check to bins under 45 kt of geostrophic wind; the table also shows the higher rows as a **deviation from the plan**, for context only. The ratio of 1.09 under 2 hPa per 100 nm shows the ratio is meaningless at weak gradients: do not apply it under 4.

| Gradient (hPa per 100 nm) | n | median ASCAT wind | ASCAT as a share of geostrophic | cells at 45 kt or more |
|---|---|---|---|---|
| under 2 | 5,950 | 15.2 kt | 1.09 | 0% |
| 2-4 | 17,583 | 22.2 kt | 0.79 | 0% |
| 4-6 | 13,442 | 28.8 kt | 0.65 | 0.9% |
| 6-8 | 6,418 | 34.2 kt | 0.56 | 8.8% |
| 8-10 | 2,720 | 38.1 kt | 0.50 | 23.9% |
| 10-12 | 1,191 | 43.9 kt | 0.47 | 47.7% |
| 12 and over | 1,506 | 52.5 kt | 0.36 | 71.3% |

For a sustained-wind rule of thumb in this region: roughly two-thirds of the geostrophic wind at 4-6 hPa per 100 nm, about a half at 8-10, and under that beyond (the upper rows are lower bounds). This is less than the textbook 70 to 80% over water, possibly because these are near-storm cells (not tested; the curvature test did not support a curvature explanation). The ERA5 gust medians in the first table sit above these ASCAT winds, as expected for a gust against a mean wind and with ERA5 reading above ASCAT; the two tables cover different samples (all cells against storm-neighbourhood cells), so the ratio between them is not computed.

## What this does not show

- Only the western Atlantic, 30-45N, October to April. Pacific, higher latitudes and the open ocean east of 55W are untested.
- Gust, not sustained wind, except the ASCAT rows. A sustained-wind table for storm and hurricane force would need a 10 m wind pull (not made).
- ERA5 at 0.25 degree smooths the tightest gradients, so the upper rows describe ERA5, not the real storm. Hurricane-force rows are few events.
- The 12-hourly samples are not independent: neighbouring cells and times share a storm. The intervals resample whole winters (21, or 12 with temperatures).
- Stability uses ERA5 sea-surface and 2 m air temperature from 2013-14 on only.
- The registered gradient-wind (curvature) estimate failed its test; do not use it.

## Does it add skill to the GFS forecast? (`PREREGISTRATION_SKILL.md`, `skill_test.py`, `results/skill_test.txt`)

ASCAT wind against GFS 24 h 10 m wind, leave-one-winter-out, 407,221 box-passes, 4 winters (2021-22 to 2024-25), wind below about 45 kt. Adding the ERA5 gradient wind to the GFS wind cuts squared error 15.0% (13.5 to 16.5; RMSE 3.10 to 2.86 kt; all four winters). At GFS wind of 25 kt or more: -2.3% (-8.3 to +2.6), cannot tell; at 35 kt or more the straight-line add was worse (-21.9%). The planned secondary curved fit (M1q, not a primary test) at 25 kt or more: +5.4% (+3.3 to +7.3), four of four winters. The planned verdicts: T1 yes (interval excludes 0, gain at least 2%, all four winters), T2 cannot tell. The 35 kt row is secondary: -21.9% (-44 to -5), three of four winters worse. Adding stability (M2, secondary) did as well as the gradient: +19.7% overall, +6.6% at 25 kt or more (three of four winters).

**Do not read the 15% as forecast skill.** The gradient used is the ERA5 *analysis* at the valid time, which a forecaster does not have 24 h ahead, so the gain is likely mostly the model's position and intensity error being corrected by the true pressure field (not tested). ERA5 also assimilates ASCAT, the outcome being predicted, so part of the gain is circular. With GFS's own forecast pressure the gain should be much smaller (not tested; GFS pressure is not cached). The gradient helps only as an independent measurement, such as an observed surface analysis, buoy or ship pressure pair, or nowcast, set against the model wind.

## Isallobaric wind (`PREREGISTRATION_ISALLOBARIC.md`, `isallobaric.py`, `results/isallobaric.txt`)

The wind toward falling pressure, from the pressure tendency, the centred difference from 12 h before to 12 h after the valid time (so it uses future pressure and is not available in real time), smoothed over 1 degree. Median size by gradient (hPa per 100 nm): 2.3 kt under 2, 3.0 at 2-4, 4.0 at 4-6, 5.1 at 6-8, 6.2 at 8-10, 7.7 at 10-12, 10.3 at 12 and over; share above 5 kt 52% at 6-8 and 87% at 12 and over. Its part along the geostrophic wind runs from about -14.5 to +11.7 kt (middle 80%) at 12 and over, so it adds to or takes from the speed.
- I1: the corrected speed |Vg + Va| ranks gust better than |Vg|: Spearman +0.0271 (+0.0236 to +0.0304), 21 winters (meaningful size 0.01): yes.
- I2: adding it to the GFS baseline against adding |Vg|: +0.92% (+0.51 to +1.31) squared error, below the 2% bar, "no (well powered)"; at GFS wind of 25 kt or more (secondary, not a primary test) -0.34% (-1.25 to +0.65), also inside plus or minus 2%. Three of four winters gain. It needs future pressure, so it could not be used in real time anyway.
- I3: at 6-10 hPa per 100 nm the gust-to-geostrophic ratio is lower where pressure falls by 1 hPa per 12 h or faster than where it rises by 1 or more: -0.038 (-0.050 to -0.027). Two-sided test; the plan set no meaningful size and predicted no direction, so no yes/no verdict. Rising pressure is the windier side at a given gradient. A possible reason, not tested: rising pressure sits behind the low in colder, unstable air.
- Quasi-geostrophic form only (no friction or curvature). All of this is ERA5 and the 12-hourly sample.

## Checks

A fresh agent recomputed the tables, bootstrap intervals, P1, P2, P2b, the spot-checked gradients (within 0.3%, a length-convention difference), the gust match (exact) and the ASCAT table from the raw files; all matched (ASCAT within a few cells). A second fresh agent rechecked the skill test, the isallobaric tests, P3 (-0.0581), the no-root shares (3.9% core, 6.3% non-core), P1c and the 8-10 core/non-core intervals: all matched, except P1b (restated above). It used the stored core flag, so core membership was not recomputed from raw pressure. The first verifier also caught five wording errors, fixed here.

The claim review (a third agent, before merge) and what it changed: the sentences above on gust versus sustained wind, the P1/P2/P3 and I2/I3 verdict words, the ASCAT rows above the ceiling, the single-band core lead, the circularity of ERA5 against ASCAT, and the note that the tendency uses future pressure. The post hoc items (P2b, P1b, P1c, the single-bin core comparison) are logged in `PREREGISTRATION.md` under "Changes after seeing results".
