# Pre-registration: the environment and the storm-scale structure of HF lows against storm-force-only lows

Committed before any ERA5 field was read for this question. **ERA5 proxy, pipeline A** (`research/era5/hf_history`: 800 km ocean gust index `g800`, HF-equivalent at
71.7 kt). Seasons 2004-05 to 2025-26 (decision 1). Transitioning tropical cyclones are in. Atlantic and Pacific separately, never pooled.

Question (Jason): what is the typical pattern for an HF low to form, and what does its composite wind field look like against a storm-force low?

Answer forms I will give: (1) *environment*: HF onset is preceded by a pattern that differs from storm-force lows' by X in region R at lead L (q < 0.05), or the
difference is not detectable (well-powered null if the minimum detectable difference is stated and small, otherwise "can't tell"); (2) *storm scale*: the HF composite
differs from the storm-force composite in the listed structure measures by X, or does not.

## What had been looked at before this plan

Only labels: group sizes, anchor times and positions, the age of each storm at its anchor, which tracks lack a heading, and the pull size. No field, no gust, no
MSLP map and no scalar built from a field. Positions were used to fix the regional boxes below. Known from this look: the median storm age at peak is 36 h for HF and
24 h for storm-force storms (so the groups are at different points of their life, which is reported, not corrected).

## Groups (pipeline A tracks, `hf_history/results/all_tracks.csv.gz`)

- **HF**: track peak `g800` >= 71.7 kt. 2,000 storms (Atlantic 1,104, Pacific 896).
- **SF (storm-force only)**: track peak `g800` in **[54.0, 71.7)** kt. 8,136 storms (Atlantic 4,032, Pacific 4,104). The bound 54 kt is 0.75 x 71.7, the ratio of the
  storm-force (48 kt) to hurricane-force (64 kt) sustained-wind thresholds, applied to the gust index. Fixed here, not tuned. **Secondary S2:** SF-high, [60.0, 71.7).
- Matching: strata are basin x calendar month of the peak time. SF is weighted to HF's stratum counts (large scale: all SF with weights N_HF / N_SF per stratum) or
  drawn 1:1 within the stratum without replacement (storm scale, seed 20261012).

## Anchors

- **HF onset**: first in-domain 6-hourly point with `g800` >= 71.7 (catalog tracks). **Peak**: pipeline A's peak time and position (HF and SF).
- **Primary comparison C1: HF onset against SF peak** (the moment of formation against the moment of the best a storm-force low managed). **C2: HF peak against SF peak.**
- Heading = bearing between the latest 00/12 UTC fix at or before t - 6 h and the earliest at or after t + 6 h on the same track
  (`intensity/results/fixes_2004.csv.gz`); one-sided at a track end. Short-lived tracks with fewer than two 00/12 fixes have no heading: they are left out of rotated
  composites and kept in north-up ones. The counts are reported.

## Part 1: large-scale environment (5.625 degree, 12 UTC daily, Northern Hemisphere 25-87N)

- Fields: Z500 (m), MSLP (hPa), U250 (m/s), SST (K), total column water vapour (kg m-2) from WeatherBench2 64x32 to 2023-01-09 and ARCO-ERA5 regridded onto the same cells
  after (method of `hemispheric/fields.py`, copied). Anomalies against a 2004-2025 calendar-day climatology smoothed over +-15 days. Composites in the geographic frame
  (basin sector), with the mean anchor position marked. All SF and all HF storms are used (weights above).
- Lags: day k before the anchor, k = 0..10, the 12 UTC field on or before anchor - k days. **Lead window L = mean of days -10 to -4. D0 = onset day (k = 0).**
- Statistic: stratified difference HF - SF = sum over strata of (N_HF,s / N_HF) x (mean_HF,s - mean_SF,s). Uncertainty: **season-block bootstrap** (2,000 resamples of the 22
  seasons within basin, seed 20261012); two-sided p from the resampled sign. Storms are not independent within a season and synoptic situations repeat across storms.
- **Primary family P1 (20 tests): box-mean difference** for 2 basins x 5 fields x {L, D0}. Box = +-15 degrees latitude x +-30 degrees longitude around the HF-onset
  centre: Atlantic 53.5N 41.1W (318.9E), Pacific 42.25N 171.5E (fixed from positions above). BH-FDR over the 20; yes if q < 0.05.
- Maps (descriptive): per-cell p from the same bootstrap, BH q < 0.05 within each map marks cells; also shown, the share of cells passing and the 95% bootstrap interval of
  the box mean. **Secondary P1b (40 tests, own BH family): box-mean difference at k = 7, 4, 2, 1.** The lead-up time series (k = 10..0) is shown with bootstrap intervals, untested beyond that.
- Power: minimum detectable difference at 80% power = 2.8 x bootstrap standard error, quoted per test, in the field's units and in storm-to-storm SD.

## Part 2: storm scale (ARCO-ERA5 0.25 degree)

- Sample: 400 HF storms per basin (random draw, seed 20261012) and 400 SF per basin matched 1:1 in basin x month. Anchors: HF onset, HF peak, SF peak (800 each; 1,894
  distinct times). Fields per time: instantaneous 10 m gust, MSLP, 10 m u and v (wind speed), 2 m dewpoint. 2 m temperature is not pulled (it would take the total past the
  50 GB gate). Size: 27.2 GB (HEAD requests on a 400-time sample; `results/sizes_storm.txt`) plus the large-scale pull (about 12-14 GB).
- Lows are re-detected with pipeline A's own rules and the 800 km owned gust is recomputed (`hf_structure/extract.py` logic); the recomputed `g800` must match the
  catalog's within 0.5 kt for >= 99% of anchors with an in-domain value, else the discrepancy is reported before any composite is read.
- Composites: 121 x 121 boxes, 25 km spacing, +-1500 km, rotated so +x is the direction of motion (right of motion is down), and north-up. Pixelwise weighted difference with the
  season-block bootstrap; BH q < 0.05 within each map. Always shown: the HF composite, the SF composite, the difference, with n.
- **Primary scalar family P2 (7 measures x 2 basins x 2 comparisons = 28 tests, BH):** central MSLP (hPa); MSLP gradient (ring-mean MSLP at 450-550 km minus central, per 100 km);
  radius of the maximum gust (km); outer radius of the 48-kt gust area (km); share of the >= 48-kt owned gust area right of motion; mean 2 m dewpoint within 500 km
  (degC); gust factor (maximum gust over maximum 10 m wind speed, within 1200 km, owned ocean). Differences in the maximum gust, in 10 m wind speed and in the HF-equivalent
  area exist by construction (HF is defined by the gust index) and are shown but excluded from the family.
- **Secondary S1: depth-matched**, SF peak against HF onset or peak weighted into 5-hPa bins of recomputed central MSLP within basin, over the bins both groups occupy; if
  fewer than 100 storms per group remain in the overlap the result is "not enough overlap". Asks whether HF differs in structure at the same depth.
  **S2** SF-high (>= 60 kt) in place of SF. **S3** drop storms whose maximum gust lies within 100 km of Greenland or Iceland (Atlantic only; terrain is not masked in pipeline A).
  Secondary tests are reported with their own BH family and labelled secondary.

## Decision rules and limits

- Differences are "significant" only at BH q < 0.05 in the stated family. A difference that is not significant is reported as such with its minimum detectable size. "No difference"
  needs 80% power at an effect stated in advance as small (box means: 0.25 storm-to-storm SD; scalars: 0.25 SD); otherwise "can't tell".
- The composites are descriptive. Nothing here is a forecast skill claim, there is no fitting, and no held-out seasons are scored, so no entry is written to `research/era5/looks/`.
- ERA5 reads low in extreme storms (gust and 10 m wind in the core), and areas in km2 are biased low; shapes and HF-minus-SF differences are more trustworthy than absolute values. Lower-confidence: Atlantic anchors north of 60N.
- Groups differ in storm age at the anchor and the SF group is much larger and more varied. The SF population is not "almost HF": a 54 kt gust index is an ordinary storm.
- The pull is resumable and uncommitted. Committed: scripts, tables under 10 MB, figures.

## Deviations (post hoc)

None yet.
