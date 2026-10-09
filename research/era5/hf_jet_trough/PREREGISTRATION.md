# Pre-registration: jet-streak quadrant and upstream-trough predictors for P(HF within 24 h)

Written 2026-10-09 before any field was read for this study and before any outcome was joined to a feature. ERA5 proxy,
pipeline A (`research/era5/hf_history`); outcome = pipeline A's 800 km ocean gust index reaching 71.7 kt. Pull approved by
Jason's card tap (October to April, about 55 GB, stop at about 60 GB). Supersedes the draft in
`/mnt/project-files/near-storm-environment/PREREG_DRAFT.md`; every change from the draft is listed under "Changes from the draft".

## Question (Jason: "A deep upper trough and jet streak are important")
Does an objective, storm-relative description of (a) where the low sits relative to a jet streak and (b) the upstream
trough add skill to the PR 12 near-storm model for P(gust index reaches 71.7 kt within 24 h)? Answer in plain words:
yes (adds skill), no (well-powered null, with the smallest gain excluded), or can't tell, and by how much.

## What had been looked at (this is not a clean slate)
- PR 12 skill and coefficients; `jet250` (1.5 deg, 1000 km maximum, storm's own wind included) has about 1/20 of the gust index's HF weight.
- PR 78 `jet_trough`: jet speed x2.7-2.9 per SD for bombs; upstream Z500 trough depth OR 1.15-1.24 per SD, not significant, under-powered.
- `intensity_extra` group J (position of the 250 hPa maximum in the heading frame): gain -0.0002. Group R (trough) was registered there and never run.
- `hf_conversion` stage B (0.25 deg): distance/bearing to the 250 hPa maximum do not separate converters (MDE ~0.2 SD at t0).
- PR 146 composites: the 250 hPa maximum lies ahead of the storm in 54% (Atlantic, 48 h before onset) falling to 33% at peak, and the mean 500 hPa low starts about 2,000 km left of the track. These came from composites, not from the features below.
No quadrant, tilt or storm-relative trough feature has been computed on these fixes, and no skill was computed with one.

## Sample
- Pipeline A fixes (00/12 UTC) in `../intensity/results/fixes_2004.csv.gz`, seasons 2004-05 to 2021-22 (18; 1 June to 31 May), fixes dated 1 October to 30 April only (user choice, "Oct-Apr is fine"). Expected: 74,009 fixes on 16,874 tracks, 4,012 with hf24. Features are also extracted for 28-30 September fixes, used only as lag sources for early-October fixes (75,000 fixes, 1,946 WeatherBench2 chunks, 55.5 GB estimate at 28.53 MB per chunk).
- Fit and test 2004-05 on only (decision 1). 2022-23 on is not used (WeatherBench2 ends 2023-01-09). Transitioning tropical cyclones are in (pipeline A tracks).
- Leave-one-season-out (LOSO) over the 18 seasons; intervals resample whole seasons. The season is the independent unit.
- Base = PR 12 `full` set (STATE + HART + ENV, `../intensity/model.py`), refitted on this Oct-Apr sample; its skill differs from PR 12's and is reported.

## Fields and frame
WeatherBench2 1.5 deg ERA5, 6-hourly: u, v at 250 and 300 hPa, geopotential at 500 hPa, read at the 00/12 UTC fix times only.
Each fix gets its own features in its own frame (north-up box 81 x 81 at 100 km, +-4000 km, bilinear, great-circle placement as in `../hf_env_composites/common.py`). A lag-L feature for a fix is the same feature computed at the same track's fix L hours earlier, in that earlier fix's own frame (L = 12, 24, 48 h; fixes are 12-hourly). A missing earlier fix gives NaN, median fill and an indicator `miss{L}`. A fix with no heading (first fix of a track) gets NaN for features that need motion (T2 bearings, `couple`) and `nohead` = 1; jet features do not need motion.
Metre/second to knot: 1.943844.

### Storm's own circulation
Before any jet or divergence feature, the azimuthal mean tangential wind about the low (100 km rings, 0-1500 km, both levels) is subtracted from (u, v), with weight 1 inside 1500 km falling linearly to 0 at 2000 km. The trough search excludes r < 750 km. Sensitivity runs (removal radius 1000 km; trough exclusion 500/1000 km; no removal) are reported for the primary tests and never used to select.

### Anomalies
No climatology pull. Z500' = Z500 minus its zonal mean on the same latitude row at the same time (m). V' = vortex-removed 250 hPa speed minus the zonal-mean 250 hPa speed (not vortex-removed) on that row (kt). 3 x 3 binomial smoothing on the box before detection.

## Predictors
**J2 jet streak** (vortex-removed 250 hPa wind, smoothed). 12 columns:
1. Streak mask = V >= 60 kt and V' >= 20 kt. Streak = the 8-connected component containing the cell of largest V within 2500 km of the low. None -> `nojet` = 1, other J2 values NaN (median fill).
2. `vmaxp` = largest V' in the streak (kt). `Lhalf` = half the along-axis extent (km/1000) of streak cells with V >= half way between 60 kt and the streak maximum V. Axis = principal axis of the component weighted by V - 60 kt, oriented downstream (positive dot with the mean wind in the component). `a` = downstream unit vector.
3. (s, n) = coordinates of the low in the streak frame centred on the streak's V maximum: s = r_low . a (+ low downstream of the maximum), n = r_low . n_hat with n_hat the left normal (+ left of the flow facing downstream), all with r_low measured from the maximum. Quadrant one-hots `RE` (s<0, n<0), `LE` (s<0, n>0), `RX` (s>0, n<0), `LX` (s>0, n>0), defined when |s| <= 2000 km and |n| <= 1500 km, otherwise all zero. Continuous `s` and `n` in 1000 km.
4. `dvds` = (V(s+300 km) - V(s-300 km)) / 0.6 along the axis at the low's projection, kt per 1000 km (+ = accelerating, entrance).
5. `div300max` = maximum vortex-removed 300 hPa divergence (1e-5 /s) within 1500 km; `div300dist` = its distance (km/1000). Divergence on the box by centred differences at 100 km, tangent-plane (metric terms ignored).

**T2 upstream trough** (Z500', smoothed). 8 columns:
1. Sector = bearing within +-75 deg of the direction opposite to the motion, 500-3000 km. Trough = the deepest local minimum (8-neighbour) of Z500' in the sector. None -> `notrough` = 1, rest NaN.
2. `tdepth` = -Z500' at the minimum (m); `tdist` (km/1000); `tbear_cos`, `tbear_sin` = cosine and sine of the minimum's bearing minus the heading.
3. `tamp` = max Z500' over cells east of the trough minimum within 3000 km of the low, minus Z500' at the minimum (m).
4. `ttilt` = atan(b cos(lat_t)) in degrees, b = slope of longitude against latitude from, for each 1.5 deg latitude row within +-8 deg of the trough latitude, the longitude of the Z500' minimum within +-30 deg of the trough longitude (rows whose minimum is on the window edge are dropped; fewer than 5 rows -> NaN, median fill). Negative = axis NW-SE, positive = SW-NE.
5. `tphase` = signed perpendicular distance (km/1000) of the low from the fitted axis through the trough minimum, + when the low is east (downstream) of the axis.

**Standardisation.** Every continuous column is standardised by basin and calendar month, with mean and SD taken over all extracted fixes (predictors only, no outcome is used; a pre-registered simplification of per-fold statistics). Indicators and one-hots are not standardised. The model's own Prep (median impute, clip 0.5/99.5, standardise) is then applied inside each fold as in the base.

**Lags.** Lag set {0}: the 20 columns at the fix. Lag set {0, 24}: 20 + 20 columns + `miss24` + the 24 h changes `dvmaxp`, `dtdepth`, `dtamp` (current minus lag 24). Lag set {0,12,24,48}: all four lags, `miss12/24/48` and the three changes.
**Coupled indicator** `couple` = (`LX` or `RE`) and `tdist` <= 2 and `tdepth` above the 75th percentile of `tdepth` over all extracted fixes (outcome-free threshold).

## Detection QC (outcome-free: run before any outcome column is read; the QC script loads only track, time, basin, lat, lon, heading)
- Q1 a streak is found for at least 80% of fixes whose vortex-removed 250 hPa maximum within 2500 km is at least 60 kt; a trough for at least 70% of all fixes (both bases).
- Q2 physical sign: the mean vortex-removed 300 hPa divergence within 1000 km of the low is higher for quadrants RE and LX pooled than for LE and RX pooled, in both basins, season-block bootstrap (2000) one-sided q < 0.05 (Benjamini-Hochberg over the 2 basins). This checks that the quadrants are real; it is not a skill test.
- Q3 twelve cases (four named storms in the probability README, eight random with a fixed seed) drawn with the streak axis, V maximum, quadrant, trough minimum and axis, reviewed by eye and logged; failures are logged, not tuned away.
If Q1 or Q2 fails, the features are not used for the skill test and the failure is the result. No threshold is changed after seeing a QC result; any change is a logged deviation and Q1-Q3 are rerun.

## Skill tests (LOSO, logistic, L2, C = 1, as the base)
Gain = pooled LOSO BSS(base + group) minus BSS(base), both against the basin x month climatology fitted on the training seasons; 90% interval and sign-flip p from whole seasons (as `../intensity_extra/evaluate.py`, same code reused). Rule: gain >= 0.005 and 90% interval above 0 and q < 0.05 -> "adds skill". Upper interval end below 0.005 -> well-powered null (a gain of 0.005 BSS or more is excluded). Otherwise can't tell.
- **Primary family (BH over 3):** J2, T2, J2+T2 added to base; target hf24; lag set {0, 24}.
- **Secondary families (BH within each):** the same three groups for (i) hf24 on onset fixes only, (ii) rapid deepening (`ndr24 >= 1`), (iii) hf48; (iv) lag set {0} and (v) lag set {0,12,24,48} for hf24; (vi) `couple` alone on hf24. All reported, none promoted afterwards.
- Season checks: gain per season (count better out of 18); by basin, descriptive.
- Yes/no scores (POD, FAR, CSI, HSS at the training-count-matched cut, HSS change with its season-bootstrap interval) for any group that passes.
- Test count and FDR: report passing tests out of the total, within each family and across all.

## Descriptive companion (no model)
HF onset (first 00/12 fix with gust index >= 71.7 kt) against storm-force-only peaks (tracks whose maximum gust index is 54 to below 71.7 kt, at their peak fix), October-April, 2004-05 to 2021-22. Storm-force-only fixes are weighted to the HF count in each basin x month cell. Quantities: quadrant occupancy, `tdepth`, `ttilt`, `tphase`, `vmaxp`, `div300max` at the anchor and at lags 12-48 h. Season-block bootstrap (2000), Benjamini-Hochberg over all comparisons; labelled descriptive.

## Pull plan and gates
1,946 chunks x 28.53 MB (u 10.36 + v 10.51 + z 7.66; measured on 30 chunks) = 55.5 GB; the lag times sit in the same chunks. Newest seasons first, resumable, features and QC drawings only are kept (no raw field is saved). The script counts bytes read and stops at 60 GB, in which case a new card is posted. Nothing raw is committed; nothing over 10 MB is committed.

## Limits
Not a forecast test (perfect prognosis: ERA5 at the fix and earlier). 1.5 deg smooths streak cores (V reads low) and cannot resolve fine structure. Tangent-plane box; storm circulation removed only as an azimuthal mean. No PV (the WeatherBench2 PV array is empty). Pipeline A only; the ERA5 record is a proxy. Oct-Apr only, so nothing here applies to summer lows. 2022-23 onward and the 0.25 deg period are untouched.

## Changes from the draft
Oct-Apr scope and the 18-season, 74,009-fix sample; vortex-removal taper 1500 to 2000 km; (s, n) measured from the streak's V maximum; `ttilt`, `tphase` and `tamp` made operational; standardisation by basin and month uses all extracted fixes (outcome-free); `couple` threshold outcome-free; lags are lookups of the track's own earlier fixes (12-hourly); the descriptive companion defines the comparison group directly.

## Deviations (post hoc)
- Clarifications found by the independent re-implementation (code unchanged, wording only): the 300 hPa divergence is computed on the 3 x 3 smoothed vortex-removed wind; the streak frame is centred on the largest V in the whole connected component (not only within 2500 km); azimuthal rings are binned by rounding r/100 km with the 1500 km ring held out to the 2000 km taper; the zonal mean is the mean over the grid row.
- Post-run deviations are listed in README.md.
