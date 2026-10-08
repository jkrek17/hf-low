# Pre-registration: Kuroshio-Oyashio SST-front strength and Pacific HF-centre hours (RA-16, ERA5 proxy)

Written and committed **before any SST index was computed and before any season's box hours were joined to a predictor.**
Post-hoc changes go under "Deviations (post hoc)" at the bottom, with both versions reported.

Question (agenda RA-16, from literature Q11): with seasons as the unit, does the strength of the SST gradient in the
Kuroshio-Oyashio zone relate to hurricane-force (HF) centre hours in the 5 x 10 degree box 40-45N, 160-170E (the atlas's busiest Pacific box)?

Plain-words answer I will give: **yes** (positive, p < 0.05), **no** (not reachable here, see decision rules), or **cannot tell**, with the correlation and its
interval, and the correlation this study could have detected stated first.

## What I had looked at before writing this

- The atlas's descriptive result that the box is the busiest Pacific box: 43.1 (archive) and 47.7 (pipeline A proxy) HF-centre hours per season, averaged over 22 seasons.
  That is a mean, not a season series. I have not opened any per-season box series, any SST gradient, or the relationship between them.
- Earlier nulls: ONI against basin-wide Pacific counts (freq_split) and EMI against Kuroshio genesis (PR 38, `enso_kuroshio`); neither is a front-strength test.
- Inventories only (file columns, chunk sizes, the ERA5 SST product note below).

## Outcome (unit: season, June-May, labelled by start year)

- **Primary outcome H_s** = HF-centre hours in the box in season s: 6 h x the number of pipeline A fixes (gust800 >= 71.7 kt, in-domain, Pacific) with
  40 <= lat < 45 and 160 <= lon < 170 (longitude as degrees east, `lon % 360`; the box does not cross the dateline). Source: `hf_history/results/era5_hf_catalog_tracks.csv`
  and the catalog events, exactly as `climatology_atlas/atlas.py::load_proxy`. **Pipeline A, a proxy**; tropical-cyclone-linked events are in.
- Phase split of the same fixes: a fix is **deepening-phase** if its time is at or before the time of that track's minimum MSLP, **mature/decay-phase** if after
  (no tunable parameter). Hours H_dep, H_mat; H = H_dep + H_mat.
- Storm count C_s = distinct tracks with at least one such fix in the box.
- Secondary source: OPC archive HF-category fixes in the same box (`load_archive` in the atlas), seasons 2004-05 onward.

## Seasons

- Pipeline A gust-based outcome: **2004-05 to 2025-26 (22 seasons)**, Decision 1 (nothing earlier is used).
- ERA5's SST product changes in September 2007 (HadISST2 before, OSTIA after). Front gradients are sharper in OSTIA, so a pre-switch period would add a level break
  to the predictor. **Primary sample = seasons 2007 to 2025 (19 seasons)**: autumn 2007 onward is entirely OSTIA. The 22-season analysis with a product dummy is a registered sensitivity test (S1).

## Predictor: SST-front strength (ERA5 0.25 degree SST)

- SST sampled at 12 UTC on days 3, 6, 9, ..., 27 (9 samples a month) for September-February of each season; monthly mean per grid cell. Pull size: 22 seasons x 6 months x 9 samples = 1,188 global chunks x 1.38 MB = **about 1.6 GB** (measured chunk 1.383 MB). No cache covers it. This is far below the 50 GB gate.
- Smooth the monthly mean SST with a 1 x 1 degree box filter (5 x 5 cells, NaN-aware). Meridional front strength G(x, y) = -dSST/dy in K per 100 km (centred difference, positive where SST falls to the north).
- For each longitude in 150-170E take the maximum of G over a latitude band, then average across the longitudes. Longitudes are degrees east in 0-360; no band crosses the dateline.
  - **FRONT (primary)**: band 32-46N.
  - **FRONT_KE**: band 32-40N (Kuroshio Extension front).
  - **FRONT_OE**: band 38-46N (Oyashio / subarctic front, overlapping the box).
- **Primary predictor = FRONT averaged over September-November (SON)** of the start year, i.e. before the December-March peak: a lagged predictor, so the HF storms
  of the season being explained cannot have imprinted on it (cold wakes, mixing). The concurrent December-February mean is a labelled contrast (S4) because storms feed that field.
- Standardised over the sample seasons. The latitude bands, longitude range, smoother and months are fixed here and are not re-chosen.

## Primary test

**P1.** Spearman rank correlation rho between FRONT(SON) and H_s over the 19 seasons 2007-2025, two-sided, by permutation (50,000 permutations of seasons; with n = 19 the null is exact in the limit).
H1 is directional (rho > 0: stronger front, more HF hours), but p is two-sided. Effect size beside it: quasi-Poisson slope of fixes (H/6) on the standardised index (ratio per SD, robust interval),
and the Pearson r for comparison with the agenda's number.

## Secondary tests (all reported; Benjamini-Hochberg over P1 and S1-S10 = 11 tests)

| id | outcome | predictor | sample |
|---|---|---|---|
| S1 | H | FRONT(SON), both variables with the product-era mean removed (dummy for seasons before 2007) | 2004-2025, n = 22 |
| S2 | H | FRONT_KE(SON) | 2007-2025 |
| S3 | H | FRONT_OE(SON) | 2007-2025 |
| S4 | H | FRONT(DJF), concurrent contrast | 2007-2025 |
| S5 | H_dep | FRONT(SON) | 2007-2025 |
| S6 | H_mat | FRONT(SON) | 2007-2025 |
| S7 | C (storm count) | FRONT(SON) | 2007-2025 |
| S8 | archive box hours | FRONT(SON), outcome with a dummy for seasons from 2013-14 removed (recording practice changed) | 2007-2025 |
| S9 | H | FRONT(SON), both variables with a linear year term removed (partial Spearman) | 2007-2025 |
| S10 | rho(S5) minus rho(S6): is the effect larger in the deepening phase? | | 2007-2025, season bootstrap (10,000), two-sided bootstrap p |

All permutation tests for S1-S9 use the same scheme as P1 (on residuals where something was removed). The prediction from the literature (Kuwano-Yoshida and Minobe 2017, not rechecked by me)
is rho_dep > 0 and rho_mat near 0, i.e. S5 positive, S6 about zero, S10 positive.

## Decision rules (fixed now)

- **Yes**: P1 rho > 0 with permutation p < 0.05. It "holds up" only if its BH q is below 0.10 over the 11 tests; I will say which.
- **No (well-powered)**: would need 80% power at |rho| = 0.40. That is not reachable with 19 seasons (below), so this study **cannot return a well-powered no**. If P1 is not significant the answer is
  "cannot tell", reported with the Fisher interval for rho so that the reader can see which correlations are excluded.
- **Power**: analytic minimum detectable |rho| at 80% power, two-sided alpha 0.05, by Fisher z: **n = 19: 0.60; n = 22: 0.57** (the agenda's 0.55 for 22 seasons is slightly low). After an
  11-test correction (alpha 0.05/11) it is about 0.73 for n = 19. A planted-effect simulation using the observed marginal distribution of H (no use of the predictor) will be added after the outcomes are joined; the analytic number is what is registered.
- A significant S5 with an insignificant S6 does **not** count as the phase prediction unless S10 is also significant.

## Checks registered with the plan

- **Autocorrelation**: report lag-1 autocorrelation of the predictor and of H. The permutation test treats seasons as exchangeable; S9 (year removed) is the guard against shared trends.
- **0/360 averaging**: the predictor and outcome both stay in degrees east 0-360 inside 150-170E; no longitude mean is taken across 0/360.
- **Outcome dominated by one slow storm**: report the largest single-storm share of hours in the box and run the rank statistic (done) and the storm count (S7).
- **Land/ice**: SST is NaN over land and ice in ERA5; none lies in 150-170E, 32-46N except the Sea of Japan margin west of 142E, outside the longitudes.

## Looks at held-out blocks

This study uses all 22 seasons including 2015-25 in a whole-sample correlation; no model is fitted on one block and scored on another. By the project's convention it counts as one scoring of 2015-25
(look 21 by the agenda count), logged in `research/era5/looks/ra16_front.log`. It uses no pre-2001 season (gust-based, Decision 1).

## What this cannot show

Mechanism (deepening, latent heat, baroclinicity) at the scale of storms; any effect smaller than the detectable correlation; ERA5 SST is a reanalysis product with a 0.25 degree grid that smooths sharp fronts.

## Deviations (post hoc)

(none yet)
