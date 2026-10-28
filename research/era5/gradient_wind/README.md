# Pressure gradient to wind speed, western Atlantic (ERA5 proxy, gust)

Plan: `PREREGISTRATION.md` (committed first, `ea5e057`, before any gust was matched to a gradient). Asked by Jason 2026-10-28; options A (table for all cells) and C (storm-core cut). Code: `build_pairs.py` (match), `analyse.py` (tables and tests), `ascat_s1.py` (secondary check). Results: `results/gradient_wind_tables.txt`, `results/ascat_s1.txt`. No download; 21 cold-season winters (October to April, 2004-05 to 2024-25), 641,808 cell-times at 00 and 12 UTC, 30-45N, 78-55W, ocean at least 1 degree from land.

**Everything is the ERA5 proxy.** The wind is ERA5 instantaneous 10 m **gust**, which reads 15 to 25% above ASCAT and runs high in unstable air. Gale, storm and hurricane-force are nominal gust classes (38, 54, 71.7 kt). ERA5 assimilates ASCAT and buoys. Neither gust pipeline is used (the cells are chosen from the pressure field alone). Western Atlantic cold season only: nothing here is tested for the Pacific.

## How to use it

Measure the gap between isobars near the area you care about and turn it into hPa per 100 nm (4 hPa isobars 60 nm apart is 6.7 hPa per 100 nm). The model gradient here is taken over 1 degree (60 nm), so measure over a similar span, not across a single tight pair of isobars.

**Rule of thumb check.** The textbook value for 4 hPa per 60 nm at 45N is a geostrophic wind of about 54 kt. In the western Atlantic cold season the gust at that gradient (the 6-8 hPa per 100 nm row) is about 42 kt at the median (37.5-45N), about 0.8 of 54 kt (the median of the per-cell ratios of gust to geostrophic wind, which differs because the latitude varies, is 0.72). The ASCAT wind in the same gradient range is about 34 kt (about 0.56 of geostrophic; storm-neighbourhood cells, see the satellite section). So the geostrophic value overshoots, and more so at tight gradients: the ratio of gust to geostrophic wind falls from about 0.8 at 4-6 to about 0.3 to 0.45 at 12 and over.

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

In plain terms at 37.5-45N: about 4 in 10 times at 4-6 hPa per 100 nm reach a 38 kt gust; about 6 in 10 times at 12 and over reach a 54 kt gust; a 71.7 kt gust happens about 1 time in 20 at 12 and over and almost never below 10 (the sample is 21 winters, so the hurricane-force rows rest on few events; the numbers of cell-times are in the results file). Medians at the same gradient differ by 0 to 3.5 kt between the two latitude bands.

### Storm core (option C)

"Core" is a cell within 216 nm of a surface low below 1000 hPa found in the pressure field itself. Gust to geostrophic ratio, **non-core minus core within the same gradient bins** (post hoc, see below): +0.045 (+0.035 to +0.054), 21 winters. At 0.6 to 0.7 that is about 7%: for the same gradient, the gust in a storm core is a little lower than outside it. The tables are close (37.5-45N, 8-10: core 46.5 kt median, 18.8% reaching storm gust; all cells 47.5 kt, 19.6%) because core cells are 22% of cell-times at 6-8, 44% at 8-10, 64% at 10-12 and 87% at 12 and over. At 30-37.5N and 8-10, non-core cells (n=207) reach a 54 kt gust 37.7% (29.4-44.9) of the time against 21.0% (15.3-26.9) in cores (n=229): the tight gradients away from lows are windier than the same gradient in a core (a possible reason, not tested: fronts and cold-air outbreaks).

### Stability

Classes of sea-minus-air difference are below 0 C, 0 to below 6 C, and 6 C or more; 12 winters (2013-14 to 2024-25) have temperatures. At 37.5-45N, median gust at 2-4 / 4-6 / 6-8 hPa per 100 nm: stable air 23 / 33 / 42 kt, weakly unstable 26.5 / 36 / 42.5 kt, strongly unstable (water 6 C or more warmer) 30.5 / 40.5 / 46 kt. The chance of a 38 kt gust at 4-6 hPa per 100 nm: 25% stable, 41% weakly unstable, 62% strongly unstable. The difference is large at weak and moderate gradients and small at tight ones (registered test P1 below). ERA5's gust adds a turbulent part in unstable air, so part of this is the gust, not the mean wind.

## Planned tests

| Test | Result | Verdict (plan's words) |
|---|---|---|
| P1 strongly unstable minus stable, gust/geostrophic, gradient 6-10 | +0.017 (-0.001 to +0.034), 12 winters, q 0.070 | by the plan's rule "no (well powered)": the interval lies inside plus or minus 0.035 (5% of a typical ratio of 0.7), but only just (upper end 0.034) and q is 0.07, so read it as "small or none at tight gradients" |
| P2 non-core minus core, gust/geostrophic, gradient 6 and over | +0.107 (+0.096 to +0.120), 21 winters, q below 0.001 | direction as expected, but the pooled statistic is mostly the mix of gradient bins (core cells sit at tighter gradients, where the ratio is smaller); see P2b |
| P3 curvature: Spearman(gust, gradient wind) minus Spearman(gust, geostrophic), core, gradient 4 and over | -0.058 (-0.078 to -0.039), 9,472 cell-times, 21 winters, q 0.0008 | **opposite sign**: the curvature-corrected wind tracks the gust less well than the plain geostrophic wind. Not a "yes" |

- P3 reasons (not tested): the isobar curvature is computed from ERA5 pressure smoothed over 1 degree and its error is large where the gradient is small, and 3.9% of core cells with gradient 4 and over have no real gradient-wind root and are dropped (6.3% of non-core cells). **The curvature correction is not recommended for a forecaster rule on this evidence.**
- Multiplicity: the three planned tests are corrected together (Benjamini-Hochberg); the descriptive tables are not tests.

## Post hoc (added after the pooled P2 was read; labelled)

- P2b: non-core minus core within gradient bins 6-8, 8-10, 10-12, 12 and over, averaged: +0.045 (+0.035 to +0.054), 21 winters. This replaces P2 as the number to quote for the core effect, and P2 is kept as registered.
- P1b: strongly unstable minus stable averaged over bins 6-8 and 8-10: +0.022 (-0.003 to +0.044), 12 winters. P1c (bins 2-4 and 4-6, where the table shows the biggest differences): +0.245 (+0.224 to +0.266). The registered range (6-10) was chosen for the gradients where the Gulf Stream work looked; the bigger effect at weak gradients was not registered.

## Satellite check (secondary, ASCAT, labelled; `results/ascat_s1.txt`)

48,810 ASCAT cells (19 winters, 2006-07 to 2024-25) within 1.5 h of 00 and 12 UTC, no rain flag, ocean at least 1 degree from land. They lie within 800 km of a tracked low, so they describe the neighbourhood of lows, not every cell. ASCAT is equivalent-neutral wind, a mean wind, not a gust, and under-reads near 50 kt, so rows from 8 hPa per 100 nm up are lower bounds.

| Gradient (hPa per 100 nm) | n | median ASCAT wind | ASCAT as a share of geostrophic | cells at 45 kt or more |
|---|---|---|---|---|
| under 2 | 5,950 | 15.2 kt | 1.09 | 0% |
| 2-4 | 17,583 | 22.2 kt | 0.79 | 0% |
| 4-6 | 13,442 | 28.8 kt | 0.65 | 0.9% |
| 6-8 | 6,418 | 34.2 kt | 0.56 | 8.8% |
| 8-10 | 2,720 | 38.1 kt | 0.50 | 23.9% |
| 10-12 | 1,191 | 43.9 kt | 0.47 | 47.7% |
| 12 and over | 1,506 | 52.5 kt | 0.36 | 71.3% |

For a sustained-wind rule of thumb in this region: roughly two-thirds of the geostrophic wind at 4-6 hPa per 100 nm, about a half at 8-10, and under that beyond (the upper rows are lower bounds). This is less than the textbook 70 to 80% over water because these are near-storm cells with curved, tight isobars. The ERA5 gust medians in the first table sit above these ASCAT winds, as expected for a gust against a mean wind and with ERA5 reading above ASCAT; the two tables cover different samples (all cells against storm-neighbourhood cells), so the ratio between them is not computed.

## What this does not show

- Only the western Atlantic, 30-45N, October to April. Pacific, higher latitudes and the open ocean east of 55W are untested.
- Gust, not sustained wind, except the ASCAT rows. A sustained-wind table for storm and hurricane force would need a 10 m wind pull (not made).
- ERA5 at 0.25 degree smooths the tightest gradients, so the upper rows describe ERA5, not the real storm. Hurricane-force rows are few events.
- The 12-hourly samples are not independent: neighbouring cells and times share a storm. The intervals resample whole winters (21, or 12 with temperatures).
- Stability uses ERA5 sea-surface and 2 m air temperature from 2013-14 on only.
- The registered gradient-wind (curvature) estimate failed its test; do not use it.

## Checks

A fresh agent recomputed the tables, bootstrap intervals, P1, P2, P2b, the spot-checked gradients (within 0.3%, a length-convention difference), the gust match (exact) and the ASCAT table from the raw files; all matched (ASCAT within a few cells). Not rechecked independently: P3, P1b, P1c, the no-real-root shares, and the 8-10 core/non-core intervals. The verifier also caught five wording errors, fixed here.
