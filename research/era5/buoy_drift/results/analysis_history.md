# Analysis history: ERA5 against moored buoys

Written 2026-10-08 **after** the comparison was run, at the coordinator's request. It is not a pre-registration. It records which choices were fixed before any ERA5-minus-buoy number was seen and which were made after, so a reader can weigh each result.

## Question and hypothesis

Question (set by the coordinator's brief before any data were pulled): do ERA5's strong surface winds drift against in-situ observations over 1979-2000, compared with later years, and is any drift the size of the gust-index drift (pipeline A +0.665 kt/decade, B +0.645, about +1.0 to +1.1 % per decade)?

Hypothesis under test: the gust-index drift is in ERA5's own surface wind or gust. If so, ERA5 against buoys should rise by about +1 % per decade over 1979-2000 and be flat afterwards.

## Fixed before any ERA5-minus-buoy number was seen

- Data: ISD moored-buoy reports (NDBC and the other hosts were unreachable); ERA5 10 m wind speed and instantaneous 10 m gust at the nearest grid point; 00/06/12/18 UTC; October to March; 1979-2004. Volume chosen to stay under 50 GB.
- Buoy gust not used, because its ISD code changes between eras (seen in the raw files before matching).
- Height adjustment 5 m to 10 m for every buoy, log profile, z0 = 2e-4 m.
- Ocean screen on the land-sea mask of the 3 x 3 neighbourhood.
- Quantities r_wind, r_gust, r_gf; bins by pair mean at 5, 10, 15 and 17.5 m/s.
- Trend model: one value per station-season, station fixed effects and station-segment fixed effects, slope on season, se clustered by season, leave-one-season-out and leave-one-station-out ranges, comparison in percent per decade with the pipeline drift.
- Break search: binary segmentation of the de-commoned monthly r_wind, |t| >= 5, |step| >= 0.03, 6 months each side.

## Decided after seeing results

Each of these came after a run of `analyse.py`. Three runs gave numbers before the final one: a partial run on 1994-2004 (code check), a full run, and a run with the sections below added.

1. **ISD unit-error screen (section 1b).** The first run showed a step of ln 0.6 (factor about 1.9) at eight stations. That is the knot factor, so it was diagnosed as an ISD unit error and those station-months were dropped by a fixed rule (|monthly median ln ratio| > 0.35). The rule was set by the size of the error, not by its effect on any trend.
2. **Long-record subset (5b), season effects (7) and correlation with pipeline B (8).** Added to look at the shape of the series.
3. **1985-2000 window (9).** Chosen after section 7 showed single-station steps of 10-20 % in 1979-1985. It is labelled in the results as chosen after looking.
4. **Per-station slope inference.** Added after noticing that the season-clustered se ignores station-to-station differences in buoy history. It is reported alongside, not instead.
5. **Near deep lows (section 10).** Added to test ERA5 against the buoys inside the gust index's own domain.

## Nulls and non-findings kept in full

- No upward drift of ERA5 against the buoys over 1979-2000 in any bin up to 15 m/s, with either fixed-effect choice. Power against +1.07 % per decade is under 10 %, so this does not show that ERA5 is free of a drift of that size.
- The gale-bin rise (+4.55 % per decade, t = +2.31) is reported, with why it is not read as an ERA5 drift.
- Every bin, quantity and window computed is in `buoy_drift-result.txt`; none was dropped from the output.
