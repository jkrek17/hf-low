# Analysis plan: storm-relative structure of HF-strength winds (ERA5 proxy, pipeline A)

Written 2026-10-08 ~07:00 UTC, before the full composite run, and committed so the
timestamp precedes the results.

**What had already been looked at when this was written (stated so the plan is
not overclaimed as blind):** a test run of `composite.py` on the first 2,363 of
the 5,983 fixes (seasons 2004-05 onward in extraction order, Atlantic summary
only, plus overall distributions of RMG and HF-equivalent area). The hypotheses
below come from the literature, not from that test, but the test was seen. No
Pacific, stage, Hart-phase or wind-speed numbers had been looked at.

## Sample

Every in-domain point of a pipeline A event track with gust index g800 >= 71.7 kt,
seasons 2004-05 to 2025-26: 5,311 distinct times, 5,983 fixes. ERA5 instantaneous
10 m gust and MSLP from ARCO-ERA5 at each time; ERA5 10 m wind speed
(WeatherBench2 0.25 deg) at one random HF fix per storm for 1,500 storms before
2023-01-10. Ownership and ocean masking are pipeline A's own (nearest low within
1200 km, ocean points). The independent unit for structure is the storm; season
for anything compared across years (no across-year comparison is planned).

## Hypotheses

- **H1 (Von Ahn et al. 2006, 17 QuikSCAT cases, earth-relative).** HF winds sit
  south of the centre. Test: share of the HF-equivalent area (owned ocean gust
  >= 71.7 kt) in the southern half, and share of fixes whose maximum gust lies
  south of the centre. Supported if both exceed 50% with the storm-resampled 90%
  interval excluding 50%, in each basin.
- **H2 (Browning 2004, motion-relative).** The maximum is right of the track and
  behind the centre (the cold-air side of the bent-back front). Test: share of
  HF-equivalent area right of motion and in the rear-right quadrant; same rule.
- **H3 (life cycle).** HF area is largest near minimum pressure and the radius of
  maximum gust (RMG) increases from deepening to filling as the wind field
  broadens. Test: medians by stage (hours from track minimum MSLP: deepening
  <= -12 h, mature within 6 h, filling >= +12 h), with intervals. Supported if
  mature area exceeds deepening and filling, and filling RMG exceeds deepening
  RMG, intervals not overlapping.
- **H4 (structure type).** Thermally symmetric warm-core (seclusion-like, Hart
  B < 10 m and VTL > 0) fixes have a smaller RMG and a more compact core than
  asymmetric cold-core (frontal) fixes. Test at 00/12 UTC fixes with Hart terms.
- **H5 (basin).** No prior; report Atlantic minus Pacific for RMG, area and the
  quadrant shares.

## Planned outputs

`structure.txt` (all of the above, including those that come out null),
`fixes.csv`, `composite.npz`, motion-relative and north-up composite figures.
The ERA5 under-resolution caveat applies to every absolute area and wind value.
Terrain-affected fixes (maximum gust within 100 km of Greenland or Iceland) are
flagged and the headline repeated without them; they are not investigated here.

Null or contrary results are reported in full in `structure.txt`.

## Fixed choices (set before the full run; `composite.py` at `2076a8d` implements them)

- Ownership radius 1200 km; index radius 800 km (pipeline A's own).
- Thresholds: HF-equivalent gust 71.7 kt (pipeline A calibration); also reported,
  not tested: gust 64, 48, 34 kt; core area >= 90% and >= 80% of the storm's
  own maximum gust; sustained wind 64, 51.2 (64/1.25, Jelenak et al.), 48, 34 kt.
- Quadrants: 90 degree sectors clockwise from the direction of motion (FR, RR,
  RL, FL) and from north (NE, SE, SW, NW). Motion = bearing from the point 6 h
  before to 6 h after.
- Stages: hours from track minimum MSLP, -12 / +-6 / +12 h cut points.
  Hart: B 10 m and VTL 0 cut points, 00/12 UTC only.
- Terrain flag: maximum gust within 100 km of land whose nearest land point is
  at or north of 59N between 285E and 350E.
- Composite box +-1500 km at 25 km, bilinear; frequency = share of fixes with an
  owned-ocean box point at or above 71.7 kt.
- Intervals: 90%, 1000 storm resamples, seed 7.
- Every comparison listed above is reported, whatever it shows.

## Deviations log

Any change after the full results are seen goes here, with both versions reported.

- (none yet)
