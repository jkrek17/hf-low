# Life cycle of HF-equivalent lows: questions, planned method, and the order of work

ERA5 proxy, pipeline A (`research/era5/hf_history`). Thread "HF low climatology gaps", 2026-10-08.

## This note was written after the analysis ran

It does not pre-register anything. The method was chosen before the results were seen, but it was only written down in the docstring of `lifecycle.py`, which was committed with the results in `896f39d`. There is no earlier timestamp. The analysis was run twice:

1. The full analysis.
2. A rerun after seeing the results, to add bootstrap intervals to the medians by onset month.

No definition, threshold or sample was changed after the results were seen. Treat every number here as descriptive and exploratory, not as a confirmatory test.

## Questions

The analysis is descriptive, so it asks questions rather than testing hypotheses.

1. How long after a low is first detected does it reach HF-equivalent winds, and how long do they last?
2. Do HF winds come before, at, or after the deepest point of the storm? The expectation from the literature is "during late deepening". Von Ahn et al. (2006) tie HF winds to the sting region, which is a late-deepening feature.
3. Do the Atlantic and Pacific differ, and does the timing change through the season?
4. Do the ERA5 onset times agree with the archive's first HF fix?

## Method

The method is fixed in `lifecycle.py` and was chosen before the run.

**Sample and definitions**
- Main sample: seasons 2004-05 to 2025-26. Earlier seasons are reported separately, because the gust index drifts upward before 2001.
- An HF fix is an in-domain fix with a gust index of at least 71.7 kt. This is the same rule that makes a track an event.
- Genesis is the track's first fix. It counts as observed only if that fix is at 1000 hPa or above and north of 21N.
- Onset counts as censored if the first in-domain fix is already HF, or if an out-of-domain fix at or above 71.7 kt comes before it.
- Tropical-cyclone-linked events are reported apart from the main numbers. An event is linked if it passes within 400 km of an IBTrACS point at the same time.

**Outcomes**
- Time from genesis to onset.
- Time at HF: the number of HF fixes times 6 h.
- The number of separate HF episodes.
- The lag from onset, and from peak gust, to the time of minimum pressure.
- The lag from the middle of the fastest 24 h deepening to onset.
- Latitude at onset.

**Uncertainty**
- 95% season-block bootstrap intervals: seasons resampled within basin, 2,000 draws.
- Contrasts between basins are reported as differences with intervals.
- No significance tests were run. None were planned.

**Check against the archive**
- One-to-one matched events: ERA5 onset compared with the archive's first HF fix.

## Null results kept

Several results are null and stay in full:
- The basins do not differ in the time from genesis to onset (0 h difference, interval [0, 6] h).
- They do not differ in the time spent at HF (0 h difference, interval [0, 0] h).
- Seasons 1979-80 to 2000-01 give the same medians as 2004-05 onward. This does not show the gust drift is harmless for timing.

The full output is in `lifecycle.txt`.
