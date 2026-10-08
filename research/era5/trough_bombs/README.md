# Upstream trough depth and bombs, 1979-2003 (RA-24) — ERA5 proxy, pipeline A

Pre-registered in `PREREGISTRATION.md` (commit `c0fc1ca`) before any field was read. **This is a proxy.** Pipeline A
(`research/era5/hf_history`) re-run on 0.25 degree ERA5 surface pressure only, 6-hourly, June 1979 to May 2004;
a bomb is the intensity framework's rapid-deepening class (24 h normalised deepening rate of 1.0 Bergeron or more),
a pressure-defined quantity, so Decision 1 (gust drift) does not bar these seasons. Tropical cyclones are in.
Trough fields are WeatherBench2 5.625 degree ERA5 (the grid PR 78 used), not 1.5 degree.

## Answer in plain words

**Yes: a deeper upstream trough goes with more bombs, about 1.2 times the odds per standard deviation of trough depth
(climatological-anomaly definition D1: 1.21, 95% interval 1.16 to 1.25), at fixed storm state, on seasons PR 78 never used.**
That sits inside the agenda's prediction (1.10 to 1.25) and PR 78's held-out estimate (1.23 and 1.24). Three cautions come with it.

- **The two definitions disagree on size, not sign.** The zonal-eddy definition D2 gives 1.75 (1.68 to 1.83). Both pass
  the pre-registered rule, so the verdict is yes, but D2 is inflated by geography: it keeps the stationary-wave trough, so
  it partly marks storms near the climatological bombing grounds. **Post hoc**, with position-cell fixed effects D2 falls to
  1.43 and D1 rises to 1.32; with position cells and upstream jet speed too they come together at 1.22 (D2) and 1.17 (D1).
  The honest range for "the trough effect at fixed storm state and place" is therefore about 1.2 per SD.
- **It overlaps with the upstream jet.** Adding the upstream jet speed (`jet_up`, 5.625 degree) before the trough
  takes D1 from 1.21 to 1.05 (1.01 to 1.09) and D2 from 1.75 to 1.42, so much of the trough signal is shared with the jet
  (which PR 78 found to be the stronger ingredient). The jet-adjusted result differs by definition, so what is left for the trough
  alone is unsettled.
- **It is a same-time association with the next 24 h**, not a forecast signal at a longer lead: with the outcome window shifted to
  24 to 48 h ahead, D1 is 0.92 (0.85 to 0.99) and D2 is 1.17 (1.06 to 1.27), so it does not carry to day 2.

The test had power: 80% at an odds ratio of about 1.06 per SD (simulation, an upper bound because it ignores season clustering;
the real intervals have a half-width of about 0.04 on the log scale).

## Numbers (odds ratio of a bomb per reference-season SD; season-clustered 95% interval; 25 seasons)

| Test | D1 anomaly | D2 zonal eddy | note |
|---|---|---|---|
| **Primary, pooled basins** | **1.206 (1.163-1.251)** | **1.752 (1.678-1.828)** | p 1.5e-10 and 1.6e-19; q < 0.001 |
| Atlantic | 1.317 (1.232-1.408) | 1.699 (1.587-1.819) | secondary |
| Pacific | 1.151 (1.078-1.228) | 1.823 (1.667-1.992) | secondary |
| Pooled, 1979-80 to 2000-01 only | 1.214 (1.165-1.264) | 1.786 (1.712-1.863) | 22 seasons; secondary |
| Season bootstrap, pooled | 1.167-1.249 | 1.689-1.824 | agrees with the clustered interval |

8 of 8 pre-registered tests pass FDR (q < 0.001). 112,604 eligible fixes, 3,904 bomb fixes, 25 seasons; 0 fixes dropped for missing trough.

Sensitivities (not in a family; `results/sensitivities.csv`): no year term 1.203 / 1.750; basin x covariate interactions 1.220 / 1.755;
Oct-Apr only 1.185 / 1.758; reference seasons 2004-14 alone (not fresh) 1.217 / 1.715; pooled 1979-2014 1.210 / 1.740.
Spline against linear (leave-one-season-out log-loss): D1 no evidence of curvature (p 0.32), D2 curvature p 0.01.
Bomb frequency by basin and five-season block is flat (Atlantic 2.8-3.4%, Pacific 3.4-4.1%; `results/diag_bomb_rate.csv`), so there is no sign of a record-quality step.

## What was run and what it cost

- Surface pressure 1979-2003: 79.4 GB streamed (about 36,500 six-hourly times at 2.17 MB each). WeatherBench2 Z500, U250, V250 1979-2014: 15.1 GB.
  A first fields run without V250 (9.6 GB) was discarded; a three-month tracker check took about 1.1 GB. About 105 GB in all
  against roughly 85 GB approved (Deviations 1 in the pre-registration).
- **Tracker check (pre-registered), passed:** all 71,472 committed catalog points before June 2004 were re-detected at the same
  time and position and the same pressure; of the 30,388 in-domain 00/12 UTC catalog points all are among the new fixes in the same basin;
  all 4,311 catalog tracks map to exactly one re-linked track (`results/tracker_check.txt`). Three re-extracted 2004+ months
  matched the committed fixes exactly. The trough code reproduces PR 78's committed values to 0.0 m (80,320 reference fixes).
- Fix table: 179,077 forecast fixes in 40,803 tracks, 1979-80 to 2003-04. Re-run: `extract_lows.py`, `build_fixes.py`,
  `extract_fields.py`, `features.py`, `check_tracker.py`, `power.py`, `analyse.py`, `posthoc.py` (raw fields and lows are not committed;
  the two result tables are).

## Looks spent

One scoring of 1979-80 to 2003-04 (counts as one pre-2001 look; the coordinator numbers it the 8th). Nothing scored at 2015-26.
Logged in `research/era5/looks/ra24-trough-bombs.log`.

## Verification

A fresh agent that did not read the analysis code recomputed (from the committed tables and the pre-registration only):
the table counts (179,077 / 112,604 / 3,904 new; 80,320 / 50,463 / 1,679 reference), the pooled odds ratios and intervals for D1 and D2,
the per-basin and 1979-2000 ones, the jet_up-adjusted ones, the position-cell ones (with and without jet_up), and the tracker check
(30,388 points, 4,311 tracks, none unmatched or split). Every one matched to three decimals (own statsmodels code, CR1 variance, t on G-1 df).

Not independently recomputed: the power simulation, the season-bootstrap intervals, the spline test, the lead-24, Oct-Apr, no-year-term and
interaction sensitivities, and the bomb-frequency diagnostic.

## Limits

- Trough depth is measured at the fix time, so it is partly the storm's own development (a deepening low deepens the trough behind it); the model
  holds storm state fixed (pressure, 12 h change, latitude, season), which reduces but does not remove that.
- ERA5 surface pressure before the satellite era is less constrained over open ocean; the 1979-2000 estimate (1.21) equals the full-period one,
  and bomb frequency is flat, which argues against an artefact, but cannot rule it out.
- D2 is dominated by geography (above); D1 is the cleaner definition but was not chosen for that reason after the fact: the pre-registered rule needed both.
- Within-season fixes are correlated; inference clusters by season (25 clusters), so the interval is honest about that but the simulated power is not.
