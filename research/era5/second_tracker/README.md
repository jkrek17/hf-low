# Does a second, independent cyclone tracker find the same hurricane-force lows? (RA-15, ERA5 proxy)

Plan, committed before any pull or tracker output: [PREREGISTRATION.md](PREREGISTRATION.md) (commit `6c7b61d`). Tracker code frozen before any match to pipeline A: commit `24d0216`.
Code: `common.py` (grid access, detectors, linkers), `detect.py`, `track.py`, `match.py`, `smoke.py`. Results: [results/](results/) (`match.txt` is the full readable report).

All numbers are an **ERA5 proxy** at **1.5 degrees, 6-hourly** (WeatherBench2), seasons 2004-05 to 2021-22 (18 of the 22 gust-era seasons; the store ends in January 2023). "A" is **pipeline A** (`research/era5/hf_history`, 0.25 degrees, gust within 800 km at 71.7 kt). Pipeline B is not used. The hurricane-force (HF) label always comes from A; the new trackers carry no gust.

## Answer

**Yes. Two other trackers find the same HF lows.** Of A's 1,490 HF events whose first HF fix falls between 5 October and 25 April (Atlantic 817, Pacific 673), the pressure tracker M recovers 88.9% (Atlantic) and 92.3% (Pacific), and the vorticity tracker V recovers 95.8% and 97.6%. All four lower interval bounds are at or above the registered 85% rule, so the registered answer is "agree".

| tracker | what it is | basin | recall | 98.75% interval (season blocks) | decision |
|---|---|---|---|---|---|
| M | MSLP minima, Murray-Simmonds style, different detector and linker from A | Atlantic | 88.9% (726/817) | 86.2-91.2 | agree |
| M | | Pacific | 92.3% (621/673) | 89.8-94.6 | agree |
| V | 850 hPa relative vorticity, Hodges style, a different variable | Atlantic | 95.8% (783/817) | 93.7-97.8 | agree |
| V | | Pacific | 97.6% (657/673) | 96.1-99.0 | agree |

Matching rule: the tracker has a track within 500 km (M) or 700 km (V) of at least half of the event's HF fixes at the same times. The 98.75% level is Bonferroni over the four tests; 10,000 resamples of the 18 seasons.

## How much of that is real agreement (post hoc control, not registered)

Matching at several hundred kilometres against thousands of tracks can succeed by chance, so the same rule was applied to A's HF fixes moved 15 degrees east, or moved 5 days later:

| | M Atl | M Pac | V Atl | V Pac |
|---|---|---|---|---|
| actual | 88.9% | 92.3% | 95.8% | 97.6% |
| moved 15 degrees east | 6.5% | 1.2% | 30.4% | 14.3% |
| moved 5 days later | 17.3% | 12.8% | 46.5% | 37.3% |

M is far above chance. V is also above chance but its rule is looser (700 km, a dense track set), so **M is the cleaner confirmation**, and V's 96-98% should be read as "at least as high as M", not as a tighter result.

## What else was registered, and what it says

- **Stricter matching.** Requiring half of *all* of A's fixes (including the weak early and late ones) to be matched drops recall to M 82.6% (Atl) and 89.0% (Pac), V 80.5% and 93.5%. The two Atlantic values are below 85%: the trackers find the storm during its HF stage but not always during its whole life (V only tracks vorticity of 1e-5 per second or more; A starts at 1010 hPa). Tightening or loosening the distance by a factor of two moves M to 87.4-96.9% and V to 88.2-100%. Events with at least 2 HF fixes: M 93.8% / 96.4%, V 96.7% / 99.1%; at least 3: M 94.7% / 97.7%, V 96.0% / 99.7% (Atl / Pac).
- **Reverse direction is weak.** Taking the most intense tracker tracks (as many as A has HF events that season and basin) and asking whether they are A HF events: M 36.7% (Atl) and 28.8% (Pac), V 45.9% and 44.7%. The ceiling, ranking *A itself* by minimum pressure, is 61.9% and 54.4%. So a depth or vorticity ranking is not a gust ranking even within one tracker: HF in A is a gust label, and many of the deepest lows are not HF-gust lows and the reverse. This is a statement about the label, not a tracker disagreement.
- **Counts.** Tracks of 24 h or more (peak between 5 October and 25 April, by A's basin rule): M finds about half of A's (Atlantic 209 vs 422 a season, Pacific 230 vs 443; correlation across the 18 seasons 0.66 and 0.30); V finds about 1.2 times as many (503 and 542; correlation 0.36 and 0.79). Cyclone counts depend on the tracker (an expected result, resolution and threshold); the HF subset does not.
- **Weekly HF series.** Weekly counts of A HF events against weekly counts of the tracker's top-N tracks: r = 0.59 and 0.50 (M, Atl and Pac), 0.55 and 0.51 (V), with intervals about +-0.06. Same caveat as the reverse direction.
- **Who M misses (S5, Benjamini-Hochberg within the family; V was above 95% so not run).** The missed Atlantic events are shorter (median 9 fixes against 17 for those found), mostly single-HF-fix events (median 1 against 2), shallower (median minimum pressure 971.3 against 961.5 hPa), and more often north of 60N (recall 85.1% against 91.1%, q 0.017) or within 100 km of Greenland (74.4% against 89.6%, n = 39, q 0.013). The same contrasts hold in the Pacific for length, count and depth. TC-linked events are not missed more (q 0.42-1.0).

## What this means for earlier results

- Pipeline A's HF events are not an artefact of its own detector or linker, or of tracking pressure rather than vorticity. The results built on A's HF events (PR 14, 41, 47, 64 and the rest) keep their standing on this count; nothing needs re-running because of the tracker.
- The most exposed are results that lean on **short, shallower, high-latitude Atlantic** HF events: the Atlantic north of 60N and Greenland barrier flow (PR 29, 58, 71). A pressure tracker at 1.5 degrees misses about a quarter of the events within 100 km of Greenland. The vorticity tracker does not (95.8% overall; the Greenland subset was not broken out for V because V was above 95%). That is a resolution and variable statement, not evidence against these events.
- **Not tested here:** whether PR 14, 41 or 47 keep their effect sizes under a second tracker. There is no gust on the new tracks, so the share of cyclones that become HF cannot be formed from them. This is a limit of the study, not a result.

## What this does not show

- Two trackers on one reanalysis agreeing says the HF subset is not an artefact of one *algorithm*; it does not say ERA5 has the storms right (ERA5 reads low in extreme storms, and the HF label is a gust proxy).
- 1.5 degrees: both trackers can drop small intense lows; recall at 0.25 degrees would be at least as high for M. The 2022-23 to 2025-26 seasons are not covered. The fit/test rule (Decision 1) is unaffected: this is an agreement check, with no teleconnection test and no held-out look spent. It is the first tracker-agreement look at the 18 seasons.
- The reverse-direction agreement is weak for reasons given above; a reader who wants "do the trackers agree on *which lows are the strongest*" should use the ceiling comparison, not the raw share.

## Data and reproduction

    python3 research/era5/second_tracker/detect.py M 1979 2022 6   # about 4.5 GB, resumable, writes work/M/
    python3 research/era5/second_tracker/detect.py V 2004 2021 6   # about 21 GB (vorticity, all levels per chunk), writes work/V/
    python3 research/era5/second_tracker/track.py M; python3 research/era5/second_tracker/track.py V
    python3 research/era5/second_tracker/match.py research/era5/second_tracker/work research/era5/second_tracker/results

`work/` is ignored (detections 64 MB, track files). Needs numpy, pandas, scipy, numcodecs, statsmodels. Tracks M: 83,038 (1979-2022); V: 50,040 (Oct-Apr 2004-05..2021-22).

## Verification

A fresh Sonnet agent that had not seen `match.py`, the results or this README implemented the matching rule from the registered text alone and recomputed, from the committed pipeline A files and the tracker outputs: the event counts (906 Atlantic and 735 Pacific in 2004-05 to 2021-22; 817 and 673 in the 5 October-25 April set), the track counts (M 83,038, V 50,040), the four recalls (726/817, 621/673, 783/817, 657/673) and both chance controls (53, 8, 248, 96 matched for the 15-degree shift; 141, 86, 380, 251 for the 5-day shift). All agree exactly. Its own 95% season-block intervals (M Atl 86.9-90.8, Pac 90.4-94.1; V Atl 94.2-97.4, Pac 96.4-98.8) are the same width class as the registered 98.75% ones above; they are a different level and seed, so the registered intervals themselves were not independently recomputed.

**Not independently checked:** the 98.75% intervals, the stricter-match and sensitivity numbers (S4), the reverse-direction shares and the A-by-minimum-pressure ceiling (S1), the seasonal-count correlations (S2), the weekly correlations (S3), all of S5 (missed-event contrasts, Greenland distances, q values), and the counts of tracks per season quoted in the Counts bullet.

## Post-hoc deviations

See the end of [PREREGISTRATION.md](PREREGISTRATION.md): chance control (added), event count window, detector wording, S5 flags, S3 scope.
