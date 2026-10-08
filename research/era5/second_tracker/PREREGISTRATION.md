# Pre-registration: does a second, independent cyclone tracker find the same hurricane-force lows? (RA-15)

Written before any second-tracker fix was computed and before any ERA5 field was pulled for this study. Everything below the "Already looked at" line is fixed by this commit; changes go in "Deviations (post hoc)" at the end, with both versions reported.

## Question

Every pipeline A result rests on one detector and linker (`research/era5/hf_history`: Gaussian-smoothed 0.25 degree MSLP minima below 1010 hPa, nearest-neighbour linking with pressure-change penalty, gust within 800 km). Do two other trackers, built on different algorithms and (for one) a different variable, find the same hurricane-force (HF) events? All results are an ERA5 **proxy**. "A" = pipeline A (`hf_history`), "M" = tracker M below, "V" = tracker V below. Pipeline B is not used.

Answer will be given in plain words: yes (agreement at or above 85%), no (well below), or can't tell, with the number and its interval.

## Already looked at (before this file)

Only inventory facts about pipeline A, none of them an outcome of this study: A's catalog has 4,157 events (role "event") 1979-2025; 1,641 in seasons 2004-05 to 2021-22 (Atlantic 906, Pacific 735), of which 1,512 have their peak fix in October-April. A's catalog track file has fixes with `g800` (gust within 800 km); an HF fix has `g800 >= 71.7`. Chunk sizes of the WeatherBench2 store were measured (4.5 GB for MSLP 1979-2022 all year; 11.05 MB per 2-day chunk of vorticity, all levels). No second-tracker output exists. Jason approved the pull (about 26 GB) on the decision card, 2026-10-08 16:09Z.

## Data

ERA5 via the public WeatherBench2 store `1959-2023_01_10-6h-240x121_equiangular_with_poles_conservative.zarr` (1.5 degree, 6-hourly, ends 2023-01-10). It is coarser than A's 0.25 degrees, so pressure depths are shallower and small intense lows are smoothed; this is a stated limitation, not tuned away. Seasons 2022-23 to 2025-26 are not covered (18 of A's 22 gust-era seasons are).

- M input: `mean_sea_level_pressure`, 1979-01-01 to 2022-12-31, all year, 6-hourly, 20-75N, all longitudes.
- V input: `vorticity` at 850 hPa, 1 October to 30 April plus a 4-day margin each side, seasons 2004-05 to 2021-22, same grid and domain.

## The two trackers (parameters fixed here, literature-style, not tuned against A)

**M, Murray-Simmonds style (MSLP, different detector and linker than A).** Detect minima of MSLP on the 1.5 degree grid, refined by a biquadratic fit in a 3x3 neighbourhood; keep a minimum only if the Laplacian of MSLP (mean over the 8 neighbours, per degree squared) is at least 0.2 hPa per degree squared, the central pressure is below 1020 hPa, and it is the lowest within 5 grid points (about 750 km). No 1010 hPa cut and no pre-smoothing (A uses both). Link with a global one-to-one assignment (Hungarian) per time step on distance from the position extrapolated by the last full displacement (first step: previous position), maximum 1000 km per 6 h and 1500 km from the previous fix, cost = distance only (A adds a pressure term and half-weight persistence). A track needs at least 4 fixes (24 h).

**V, Hodges style (850 hPa relative vorticity, different variable).** Smooth the vorticity field with a Gaussian of 1 grid point; detect maxima at least 1.0e-5 per second (northern hemisphere), lowest-neighbour rule of 5 grid points. Link by Hodges' cost minimisation: the cost of a three-fix triplet is the sum of a direction-change term and a speed-change term (equal weights, Hodges 1994/1999 form, 1000 km per 6 h limit), solved per step by the Hungarian algorithm over the extrapolated positions. A track needs at least 4 fixes.

The code is committed and frozen before any match to A is computed (commit hash recorded in the README). Allowed checks on the trackers before that point do not use A: tracks per month, lifetime distribution, a handful of plotted tracks. Any parameter change after the match step exists is a deviation, reported with both versions.

## Matching rule (primary)

An A HF event is a catalog event with its HF fixes H (fixes with `g800 >= 71.7`). A tracker track T **matches** the event if T has a fix at the same time within D km of A's fix for at least half of the fixes in H (at least one), and T lasts 24 h or more by construction. D = 500 km for M and 700 km for V (the vorticity centre sits farther from the pressure centre in a deep occluding low). If several tracks qualify, the one with the most matching fixes counts. The sensitivities D/2 and 2D, and a whole-track rule (at least half of all A fixes matched), are registered secondaries.

## Primary tests (family P: 4 tests)

For each tracker (M, V) and each basin (Atlantic, Pacific, by A's basin tag): **recall** = share of A's HF events, seasons 2004-05 to 2021-22, that the tracker matches.

- M is scored on all-year events and, for like-for-like with V, on the October-April subset. V is scored on events whose first HF fix falls between 5 October and 25 April (the margin keeps tracks from starting at a window edge). The four primary tests use: M on October-April events, V on October-April events; the all-year M recall is a labelled secondary.
- Interval: 95% (Bonferroni over the 4 primary tests, so 98.75% two-sided) season-block bootstrap, 10,000 resamples of the 18 seasons.
- **Decision rule (agenda threshold 85%).** Agree if the lower interval bound is at least 85%. Disagree if the upper bound is below 85%. Otherwise can't tell. The overall answer is "yes" only if all four agree; "no" if any is disagree; otherwise can't tell.

## Secondary tests (family S, BH within family and across all tests)

- S1. Reverse direction (precision). For each tracker, season and basin, take the N tracks with the most extreme intensity (lowest MSLP for M, largest peak vorticity for V; basin by A's rule at the intensity-peak fix), where N is A's number of HF events that season and basin. Report the share of them matched (same rule, with the tracker's peak fix standing in for the HF fixes) to any A HF event. Ceiling for interpretation: the same share when A itself is ranked by minimum pressure (`minp`) in place of the tracker. A gust-based and a depth-based ranking differ even within one tracker, so S1 is read as a ratio to its ceiling.
- S2. Count agreement. Seasonal number of tracks of at least 24 h in 20-75N from each tracker against A's `all_tracks.csv.gz`, 18 seasons: Pearson and Spearman correlation and the ratio of means, by basin.
- S3. Weekly agreement. Weekly counts of A HF events with the number of matched events and the number of unmatched, and the Pearson correlation between A's weekly HF count and the tracker's weekly count of its S1 top-N tracks (season-block bootstrap interval). This is the series the PR 41 design uses.
- S4. Recall sensitivities: D/2 and 2D, whole-track rule, events with at least 2 HF fixes and at least 3 HF fixes (the duration rules from the count reconciliation), M all-year.
- S5. **Who is missed** (descriptive, only if any primary recall is below 95%). Missed vs matched events compared on: peak latitude 60N or north (the lower-confidence Atlantic group), Greenland within 100 km (A's mask check), HF fix count, track lifetime, A's minimum pressure (`minp`), transitioning tropical cyclone link. Fisher or Mann-Whitney test per factor, BH within S5. This tells which earlier results lean on the missed subset.

## What the answer would mean for earlier results

Decided now so it is not argued after the numbers:

- All four primary recalls agree: the HF subset does not depend on the detector, linker or (V) the tracked variable; PR 14, 41, 47 and the others built on A's events keep their standing on this count and need no re-run. They are still a proxy and still carry the gust-drift rule (Decision 1).
- A basin or tracker disagrees: list the unmatched events by S5. Results most exposed are those that lean on that subset: Atlantic north of 60N and Greenland barrier flow (PR 29, PR 58, PR 71) if the misses concentrate there; the HF share results (PR 14, 47, 64) if misses concentrate in short or shallow events; the Pacific counts and the Kuroshio null (PR 38, which reused A's own detector) if Pacific misses are broad. A re-run of the exposed result on matched events only would be proposed, not done here.
- V disagrees but M agrees: the detector and linker are not the issue; the pressure-centred definition of "a low" is, and results that use vorticity-like structure should be checked.

## Power

Recall with about 800 Atlantic and 600 Pacific October-April events over 18 seasons: a binomial standard error near 1.3% and 1.5%; the season-block interval will be wider because storms cluster in seasons. The test can tell 85% from 80% or 90%; it cannot tell 85% from 86%.

## What this study cannot show

- That A or either tracker is right: agreement of two pressure-centred or vorticity-centred trackers on one reanalysis says the HF subset is not an artefact of one algorithm, not that ERA5 has the storms right (the gust field is the proxy; ERA5 reads low in extreme storms).
- Anything about 2022-23 onward, or about gust: no gust field is pulled, so a tracker's own HF counts are not formed. Matching uses A's gust label.
- Results at 0.25 degrees: M and V run at 1.5 degrees, which can drop small intense lows (polar-low-like and some barrier-flow cases). A miss for that reason is a resolution result, and is reported with the lows' pressure and size where available.
- Effect re-runs (PR 14, 41, 47) are not registered as tests here; there is no gust attached to the new tracks.

## Looks spent

This study uses no teleconnection index and no held-out outcome, so it spends no held-out look. It is the first tracker-agreement look at the 18 seasons 2004-05 to 2021-22. The project's other pre-registered looks on the same 22 gust-era seasons are not added to.

## Deviations (post hoc)

All of these were logged after the first full run of `match.py` (commit `24d0216` froze the tracker code before that run; the match script was committed before it was run). The primary results did not change between runs; later edits fixed crashes in the weekly (S3) and missed-event (S5) sections only.

1. **Chance control added (post hoc, not registered).** Because matching at 500-700 km against thousands of tracks can succeed by chance, recall was recomputed with A's HF fixes moved 15 degrees east, or 5 days later. Reported beside the primary results; it does not change any registered decision.
2. **Event count.** The inventory above (1,512 events in October-April) counted by the month of the peak fix. The registered scoring window is the first HF fix between 5 October and 25 April, which gives 1,490 events (Atlantic 817, Pacific 673). The registered rule was used.
3. **Detector wording made exact.** "Laplacian" is 4 x (mean of the 8 neighbours minus the centre) / (1.5 degrees)^2, in grid units without a latitude factor. "Biquadratic fit" is the closed-form quadratic through the 3x3 neighbourhood from central differences. Neither was tuned.
4. **S5 flags.** "Greenland within 100 km" uses the 1.5 degree WeatherBench2 land-sea mask (coarser than the 0.25 degree mask of `research/era5/highlat`). "Transitioning tropical cyclone" is the `tc` flag of A's own `lifecycle_events.csv`, for the events that have one.
5. **S3.** Only the weekly correlation was computed; the weekly matched/unmatched counts are in the per-event tables (`results/events_M.csv`, `events_V.csv`).
