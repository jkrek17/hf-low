# ERA5 hurricane-force-equivalent event history, 1979-2026

A proxy record of hurricane-force (HF) extratropical lows for the North
Atlantic and North Pacific, built from ERA5 and calibrated against the last
five seasons of the archive. Research code: nothing here feeds `docs/` or CI.

**This is a proxy.** An event here is "a cyclone whose ERA5 gust field looks
like the ones OPC warned for as hurricane force in 2021-26". It is well
defined, but it is not the archive and cannot be validated before 2001.

**The gust index is not era-uniform.** At fixed MSLP (955-975 hPa) the 800 km
gust rises +0.67 kt/decade over 1979-2000 (t = +2.44) and is flat after 2004,
matching the drift pipeline B found in `c9dc994`. On the full track population
(`results/all_tracks.csv.gz`), a pre-registered test finds no detectable ramp
in depth-adjusted counts within 1979-2000 (+2.8 events/decade, 95% CI -3.0 to
+8.6). But a post hoc level comparison finds 1979-2000 has 7.8% fewer events
than the same depths give after 2004 (t = -2.77), and that offset accounts for
the whole 1979-2025 trend in the table below. See `../drift_counts-result.txt`
and decision 1 in `STATUS.md` on `main`. Use counts before 2001 for variation
within that era, not for levels or trends across 2001.

## Definition

- **Lows.** Every 6 h, minima of lightly smoothed ERA5 MSLP below 1010 hPa
  over 20-75N, deduplicated within 400 km, linked into tracks (<= 900 km per
  6 h about a motion-persistence guess). Tracks shorter than 24 h are dropped.
- **Domains** (the archive's footprint): Atlantic 30-67N, 98W-10E; Pacific
  27-67N, 135E-120W. A track counts in the basin where most of its fixes fall,
  and only its in-domain fixes are scored.
- **Gust index.** At each fix, the maximum ERA5 `instantaneous_10m_wind_gust`
  over ocean within 800 km of the centre, counting only grid points nearer to
  this low than to any other. A track's index is its in-domain maximum.
- **HF-equivalent event.** Track index >= **71.7 kt**.
- **Season.** June through May, labelled by the June year, matching the
  archive's `season` field.
- **Null cases.** For each event, one cyclone from the same basin and season
  that never reaches the threshold, lasts >= 48 h, deepens to <= 1000 hPa, and
  peaks within +/-7 days of the event's peak; nearest in time, each used once.
  4,154 of 4,157 events have one.

## Calibration (seasons 2021-22 to 2025-26)

469 archive events; 453 (97%) match an ERA5 track (a fix within 400 km at the
same time; 800 km for tip-jet and centreless events). The threshold is the one
where the ERA5 event count equals the archive count (bias 1), so seasonal
totals are on the archive's scale.

    AUC 0.990 over 8,138 tracks
    threshold 71.7 kt   POD 0.77  FAR 0.23  CSI 0.63  HSS 0.76  bias 0.99

Leave-one-season-out, the refitted threshold stays within 71.4-72.3 kt and the
held-out HSS is 0.68-0.78. Atlantic and Pacific score alike (HSS 0.75 / 0.76).
Full report: `results/skill.txt`.

Applied unchanged to the earlier archive, skill holds at HSS 0.65-0.69 from
2006 on, with ERA5 counting 6-19% more events than the archive. 2001-05 shows
bias 2.05 because the archive was still ramping up, not because ERA5 drifted
(the archive has 0 Atlantic events in 2001-02 and 2002-03).

Known false-alarm source: tropical cyclones crossing into the domain before
transition (a handful per season). Masking tropical-season fixes south of
33-38N did not improve skill, so they are left in.

## Probability of hurricane force

The catalog also carries `p_hf` and `p_hf_basin`: the probability that OPC
would have warned the cyclone as hurricane force, from a logistic fit of
"matched to an archive HF event" on the gust index over every ERA5 track in
the calibration seasons (8,144 tracks, 440 positives). `p_hf_basin` adds a
Pacific term. Full report: `results/probability.txt`.

    p_hf        P = 0.5 at 73.6 kt; Atlantic over-predicted 5.9%, Pacific under 6.1%
    p_hf_basin  P = 0.5 at 74.2 kt (Atlantic), 73.0 kt (Pacific); both basins sum to observed
    basin term  likelihood ratio 4.91 on 1 df (p = 0.03); held-out log-likelihood
                improves by only 1.0 over five seasons, so treat it as marginal

The probabilities are a function of the index and basin alone, so they apply to
any catalog row. Summing them over the catalog does **not** give an expected HF
count: the catalog holds events and their null cases, not every sub-threshold
cyclone. 16 of the 469 calibration archive events have no ERA5 track and are
outside the fit.

## Results

4,157 events (Atlantic 2,254, Pacific 1,903) and 4,154 null cases across 47
seasons. Atlantic averages 48.0 per season (sd 8.1), Pacific 40.5 (sd 7.6).
Atlantic shows a weak upward trend of +1.7 per decade (p = 0.05); Pacific
none (+0.7, p = 0.42). Both trends span 2001; the depth-expected count is flat
over 1979-2025 (-0.19 per decade, both basins), so read them as the gust-drift
offset above, not as a climate trend, until independent evidence says otherwise. Over 2006-2025, ERA5 and archive seasonal counts
correlate at r = 0.62 (Atlantic) and 0.70 (Pacific).

| season | Atl ERA5 | Atl archive | Pac ERA5 | Pac archive |
|---|---|---|---|---|
| 1979-80 | 50 |  | 37 |  |
| 1980-81 | 38 |  | 41 |  |
| 1981-82 | 49 |  | 30 |  |
| 1982-83 | 42 |  | 42 |  |
| 1983-84 | 44 |  | 47 |  |
| 1984-85 | 42 |  | 32 |  |
| 1985-86 | 54 |  | 43 |  |
| 1986-87 | 46 |  | 38 |  |
| 1987-88 | 40 |  | 47 |  |
| 1988-89 | 47 |  | 22 |  |
| 1989-90 | 57 |  | 29 |  |
| 1990-91 | 33 |  | 34 |  |
| 1991-92 | 51 |  | 36 |  |
| 1992-93 | 56 |  | 39 |  |
| 1993-94 | 53 |  | 36 |  |
| 1994-95 | 45 |  | 45 |  |
| 1995-96 | 43 |  | 38 |  |
| 1996-97 | 49 |  | 40 |  |
| 1997-98 | 44 |  | 48 |  |
| 1998-99 | 47 |  | 46 |  |
| 1999-00 | 45 |  | 49 |  |
| 2000-01 | 35 |  | 60 |  |
| 2001-02 | 46 | 0 | 39 | 1 |
| 2002-03 | 62 | 0 | 51 | 22 |
| 2003-04 | 32 | 11 | 38 | 27 |
| 2004-05 | 32 | 45 | 45 | 36 |
| 2005-06 | 39 | 33 | 37 | 33 |
| 2006-07 | 58 | 67 | 40 | 48 |
| 2007-08 | 47 | 42 | 50 | 51 |
| 2008-09 | 53 | 51 | 29 | 33 |
| 2009-10 | 32 | 30 | 49 | 41 |
| 2010-11 | 48 | 32 | 36 | 30 |
| 2011-12 | 52 | 40 | 52 | 46 |
| 2012-13 | 58 | 42 | 50 | 40 |
| 2013-14 | 49 | 46 | 34 | 27 |
| 2014-15 | 62 | 57 | 32 | 31 |
| 2015-16 | 53 | 51 | 49 | 50 |
| 2016-17 | 56 | 50 | 35 | 34 |
| 2017-18 | 65 | 47 | 36 | 27 |
| 2018-19 | 53 | 46 | 45 | 35 |
| 2019-20 | 55 | 42 | 30 | 33 |
| 2020-21 | 43 | 42 | 43 | 44 |
| 2021-22 | 51 | 57 | 43 | 44 |
| 2022-23 | 43 | 42 | 37 | 41 |
| 2023-24 | 51 | 41 | 50 | 50 |
| 2024-25 | 46 | 43 | 42 | 42 |
| 2025-26 | 58 | 66 | 32 | 43 |

## Life cycle

`lifecycle.py` reads the committed catalog and track files only (no ERA5
access) and writes `results/lifecycle.txt` and `results/lifecycle_events.csv`.
Proxy, seasons 2004-05 to 2025-26, tropical-cyclone-linked events left out
(1,886 events). An HF fix is an in-domain fix at or above 71.7 kt.

    genesis -> first HF fix     median 36 h (quartiles 24/48), n = 1,152 with genesis observed
    time at HF                  median 12 h (quartiles 6/24); 35% a single 6-hourly fix
    onset before min pressure   71% (Atlantic 66%, Pacific 77%)
    peak gust vs min pressure   before 54%, same fix 26%, after 20%
    ERA5 vs archive onset       median 0 h, within 12 h for 88% of 1,244 matched events

`results/lifecycle-plan.md` records the method and says it was written after
the run. Every number above was recomputed by an independent agent.

## Files

    results/era5_hf_catalog.csv          one row per event or null case
    results/era5_hf_catalog_tracks.csv   6-hourly fixes for those tracks
    results/era5_hf_counts_by_season.csv the table above, with null counts
    results/skill.txt                    calibration and transfer report
    results/threshold.json
    results/all_tracks.csv.gz            every in-domain track 1979-2025 (75,087), one row each:
                                         gust800_kt (track index), minp, peak, season, n_fix.
                                         Re-extracted 2026-10-08; reproduces the catalog's 2,254 / 1,903
                                         events. The full sub-threshold population, so counts can be
                                         compared at fixed depth.
    results/probability.txt, probability.json  the P(HF) fits
    results/lifecycle.txt, lifecycle_events.csv, lifecycle-plan.md  life cycle

Catalog columns: `role` (event / null_case), `null_for` (the event track a
null is matched to), `gust800_kt` (track index), `minp` (hPa), `peak_*` (the
fix with the highest index), `archive_events` (matched archive ids, 2001 on).
Track ids are only unique within one run.

## Reproducing

Needs `numpy numcodecs gcsfs scipy pandas`. Streams about 370 GB from the
public ARCO-ERA5 store (MSLP and gust, 6-hourly); about 45 minutes with 12
processes. Intermediate files go to `work/` (or `$ERA5_WORK`).

    python3 extract.py 1979 2025      # lows + gust summaries, one CSV per month
    python3 track.py 197906 202605 work/track_points.csv
    python3 calibrate.py              # fits and scores the threshold
    python3 apply.py work/track_points.csv 71.7 results
    python3 probability.py work/track_points.csv results/era5_hf_catalog.csv
