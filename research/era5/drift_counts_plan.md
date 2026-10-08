# Plan, written before the test was run: do pipeline A's thresholded counts drift before 2001 at fixed depth?

Committed before any test statistic was computed. The only things looked at in
`research/era5/hf_history/results/all_tracks.csv.gz` before writing this were its
header and two rows. The analysis follows this plan; anything added after
seeing results is labelled "post hoc" in the result file.

Pipeline A = `research/era5/hf_history/` (800 km track-maximum ocean gust,
`gust800_kt`; event = index >= 71.7 kt; June-May seasons labelled by start
year). It is an ERA5 proxy; nothing before 2001 is validated.

## Question

Decision 1 in `STATUS.md`. `drift_check-result.txt` (`c55e74d`) found that the
gust index value rises at fixed MSLP over 1979-2000 (+0.665 kt/decade in A's
own fixes) and bounded the effect on A's event counts at about +1 to +6 events
per decade, without being able to tell it from a real rise in deep storms. The
full track population (75,087 tracks, `8b109b5`) now allows a test at fixed
depth over every track, not only events and null cases.

## Hypotheses

- **H0:** at fixed track depth (minimum MSLP, `minp`) and basin, the probability
  that a track reaches 71.7 kt has no trend over seasons 1979-2000.
- **H1 (drift):** that probability rises over 1979-2000, so that depth-adjusted
  counts trend upward. The size c9dc994's +0.65 kt/decade predicts is
  about +5.7 events/decade, both basins.

## Data and units

- Tracks in `all_tracks.csv.gz`, both basins. Pre-archive period: seasons
  1979-2000 (n = 22). Reference period: seasons 2004-2025 (n = 22). Seasons
  2001-2003 are excluded from both (archive starting up; RECORD_START).
- Depth strata: basin x 5 hPa bins of `minp`, with an open bin below 945 hPa
  and one above 1005 hPa.
- Exceedance: `gust800_kt >= 71.7`.
- The sample size for every trend is seasons (22), never tracks.

## Primary test (T1): depth-standardised excess

1. Reference exceedance rate r(basin, bin) = share of 2004-2025 tracks in that
   stratum with index >= 71.7.
2. For each season s in 1979-2000: observed O_s = events; expected E_s =
   sum over its tracks of r(basin, bin); excess D_s = O_s - E_s.
3. Statistic: OLS slope of D_s on season, per decade, with 95% CI and t
   (n - 2 = 20 df).
4. Also reported: the same with reference = 1979-2000 itself, and the same
   slope over 2004-2025 (control period; the per-fix gust is flat there).

## Secondary test (T2): within-stratum trend

Logistic regression, tracks of seasons 1979-2000: exceed ~ stratum fixed
effects + season. Standard error clustered on season (22 clusters). Reported as
the season coefficient, its se, and the implied events/decade (average marginal
effect x tracks per season). Also per bin: OLS trend of the per-season share,
listed without multiplicity correction and labelled descriptive.

## Power

Effect sizes: count drifts of +1, +3 and +6 events/decade. Induced by adding
d x (season - 1989.5)/10 kt to every 1979-2000 track's index, with d chosen so
that the induced change in the T1 slope equals the effect; power is computed
from that shift against the observed T1 standard error (two-sided 5%, t with
20 df).

## Decision rule (fixed in advance)

On the T1 95% CI for the 1979-2000 excess slope (events/decade, both basins):

- **(a)** upper bound below +2: counts are era-comparable to within that bound.
  Recommend gust-based counts from 1979, quoting the CI as the residual
  uncertainty and labelling the record a proxy.
- **(b)** lower bound above 0: count drift detected. Recommend counts from
  2001-02 for trends, or a depth-standardised correction that is stated.
- **(c)** otherwise: not settled; the CI replaces the +1 to +6 band.

T2 must agree in sign with T1 for (a) or (b) to be called.

## Robustness (pre-specified)

- Bin width 2 and 10 hPa.
- Each basin separately.
- Depth distribution itself: trend over 1979-2000 in tracks deeper than 960 hPa
  and in mean `minp`, to show whether the depth strata move.

## Known limits

- Depth is ERA5's own MSLP. If ERA5's MSLP drifts too, "fixed depth" is not fixed
  in the real atmosphere. The test asks whether the gust-based count drifts
  relative to ERA5's own depth, which is the artefact c9dc994 described.
- A real change in gust at fixed depth (for example in storm size or speed) would
  look the same as an artefact. Only the independent buoy and ship comparison in
  another thread can separate those.
