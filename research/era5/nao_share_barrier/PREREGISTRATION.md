# Is the Atlantic NAO share effect storm wind or barrier flow? (pre-registration)

Written 2026-10-08 before any outcome was related to any predictor in this analysis. Everything ERA5 is a
**proxy** (pipeline A, `research/era5/hf_history`; 800 km ocean gust index, HF-equivalent at 71.7 kt).
Anything changed after the first outcome is read goes under "Deviations (post hoc)" in `README.md`, with both
versions shown.

## What was looked at before writing this

- The published results of hf-low PR 14 (`freq_split/README.md`: Atlantic NAO, lagged, RR(HF) 1.096, RR(all)
  0.974, RR(share) 1.125, f 1.29), PR 11, PR 29 (`highlat/FINDINGS.md`), PR 58 (`greenland_jets/README.md`,
  including NAO-vs-high correlation -0.59 and "NAO adds nothing once the high is in").
- Column names and row counts of `all_tracks.csv.gz` and `highlat/gustloc_fixes.csv`, the match between their
  track ids (1,103 of 1,104 Atlantic HF tracks 2004-05 on have fix-level rows), and **how many HF tracks each
  exclusion rule below removes** (for the month-defined Oct-Apr Atlantic sample: terrain-type 110 of 1,007,
  100 km mask 57, 300 km mask 214, centre-only 142). No count or fit was related to NAO, PNA or the Greenland
  high.
- Not yet pulled: the daily Greenland-high series (below).

## Question

PR 14: the Atlantic NAO raises HF counts through the share of cyclones that reach HF. PR 11: NAO does not
change storm intensity. PR 29: about a third of Atlantic HF-strength fixes north of 60N take their gust from
terrain-type barrier flow. **Is the share effect real storm wind, or barrier flow?** Plain-words answer to give:
how much of the NAO share effect survives when barrier-type fixes are removed (a fraction, with an interval),
and whether the Greenland high accounts for it.

## Sample and model (identical to PR 14's primary, `freq_split/split.py`)

- Atlantic pipeline A tracks (lows below 1010 hPa, at least 24 h), genesis day = first fix. Window 1 Oct plus
  210 days, seasons 2004-05 to 2025-26 (Decision 1: gust-based, 2004-05 on). Transitioning tropical cyclones are
  in, as in PR 14.
- Daily genesis counts, Poisson, log link, NAO + calendar-month effects + linear season trend. NAO is the CPC
  daily index averaged over days -10..-4 before genesis, standardised over the analysis days. Share =
  stacked model, `b_share = b(HF) - b(all)`. Season-clustered SEs; season-block bootstrap (2,000 draws, seed
  20261008) for intervals; season-block permutation (2,000) for p of a single coefficient. The effective n is 22
  seasons x about 5.3 independent NAO values per season.
- The original HF (`hf`) is a track whose peak 800 km gust is at or above 71.7 kt. The PR 14 value must be
  reproduced first (RR(share) 1.125, 1.05-1.20); if it is not, stop and explain before going on.

## Barrier exclusion rules

Fix-level data: `highlat/gustloc_fixes.csv` (every in-domain HF-strength fix on a pipeline A Atlantic event
track, 2004-05 on, with the maximum's position and the index recomputed with terrain removed).

A track is **HF without barrier** if it has at least one HF-strength fix that survives the rule:

| Code | Rule | Source |
|---|---|---|
| **T (primary)** | the fix is not terrain-type: its gust maximum is not both more than 400 km from the centre and within 300 km of Greenland (`max_dkm > 400 and max_dgl < 300` removes it) | PR 29's classification |
| G100 | recomputed index with ocean within 100 km of Greenland removed is at or above 71.7 kt (`g800_gl100`) | PR 29 |
| G300 | same with 300 km removed (`g800_gl300`) | PR 29 |
| R400 | recomputed index from points within 400 km of the centre only is at or above 71.7 kt (`g800_r400`) | PR 29 |

A track reclassified as not-HF **stays in the denominator** as a cyclone (it is a cyclone whose storm wind did
not reach HF). Variant D drops the reclassified tracks from both numerator and denominator.

North-of-60N cut (agenda wording): **N60** removes every track with `peak_lat >= 60` from both sides (all
cyclones and HF), consistently. `peak_lat` is the position of the track's peak gust, the only position in
`all_tracks.csv.gz` for all tracks.

## The Greenland-high mediator

GH = mean MSLP (hPa) over the Greenland ice-sheet land points 64-78N, 55W-30W, **the same 0.25 degree ARCO-ERA5
field and land set as PR 58's stage 1**, but sampled once a day at **12 UTC** for every day from 21 Sep to 28
Apr of seasons 2004-05 to 2025-26 (about 4,850 fields, about 11 GB streamed, under the 50 GB gate). It is lagged
exactly like NAO (mean of days -10..-4 before genesis) and standardised over the analysis days. Check, before it
is used: at the 12 UTC times shared with PR 58's `stage1_times.csv`, my GH must equal PR 58's to rounding
(correlation at least 0.999); if not, stop. GH is an MSLP index, a reduction to sea level over the ice, not a
pressure that exists.

## Tests

Direction note, written before the run: PR 58 found a stronger high raises barrier gusts, and the high is
anticorrelated with NAO (-0.59, at case times). NAO+ goes with a weaker high, so a barrier-flow story in which
NAO+ raises the share needs a route other than the high. I therefore register every test **two-sided**, and name
a rise in the NAO coefficient when GH enters as *suppression* (barrier flow works against the NAO share effect),
which is not the same as "no change".

**Family F1 (headline, 3 tests)**

- **P1.** `b_share(NAO)` with HF replaced by rule T. Reported with its permutation p and bootstrap interval.
- **P2.** Paired difference `b_share(T) - b_share(original)` (same seasons resampled, same everything else).
  This is the test of whether the reclassified tracks carry the effect.
- **P3.** Mediation: with NAO and GH in the same model (original HF), `r = b_share(NAO | GH) / b_share(NAO)`;
  the test is the paired difference `b_share(NAO|GH) - b_share(NAO)`. Prediction (agenda RA-3): RR(share) falls
  from 1.125 toward 1.06, that is r at or below 0.5.

**Family F2 (secondary: variants of P1 and P2, 10 tests).** For each of G100, G300, R400, N60 and D (rule T with
tracks dropped): `b_share(NAO)` for the variant (permutation p) and the paired difference against the original
HF fitted on the same sample (bootstrap p). For N60 the original is refitted on the reduced sample.

**Family F3 (decomposition, 3 tests).** Split HF into storm-HF (rule T) and barrier-only HF (HF but not T).
NAO log RR for each count; and their difference. These show where the NAO effect on HF counts sits, not only
the share.

**Family F4 (the high, 4 tests).** `b_share(GH)` in a GH-only model for original HF and for rule T, and their
paired difference; and `b_share(NAO | GH)` with rule T (what remains of the storm-wind share effect after the
high).

Multiplicity: Benjamini-Hochberg q within each family and across all 20 tests. Anything not listed is post hoc.

## Decision rules (thresholds fixed now)

Let `q = log RR_share(T) / log RR_share(original)`, point estimate and 95% season-block bootstrap interval.

| q (point) | Read |
|---|---|
| 1.25 or more | the share effect **grows** without barrier fixes (barrier flow works against it) |
| 0.75 to 1.25 | the effect **survives**: it is storm wind south of Greenland or elsewhere |
| 0.25 to 0.75 | **partly** barrier flow |
| 0.25 or less | **mostly barrier flow** (RR_share(T) at or below about 1.03) |

A verdict is **firm** only if the whole 95% interval of q lies inside the row; otherwise it is stated as
"leans X, not resolved". Survival additionally needs the 95% interval of `RR_share(T)` to exclude 1. Mediation
(P3): supported if `r <= 0.5` and the 95% interval of r has an upper bound below 1; "leans, not resolved" if
`r <= 0.5` and the interval reaches 1; **prediction fails** if `r > 0.5`; `r > 1.2` is read as suppression.
Rules for variants are the same table on that variant's `q`.

A null is a full result: if q's interval lies in 0.75-1.25, the answer is "no, the barrier explanation does not
account for it". Power: expect P3 to be weak, since GH and NAO are collinear (NAO's SE grows by about 1/sqrt(1-
0.59^2) = 1.24 when GH enters), and the agenda already expects the 60N cut to be inconclusive unless the point
estimate collapses. The minimum detectable shift in a paired difference is reported as 2.8 x its bootstrap SE.

## What the study cannot show

It cannot separate proxy error from archive omission near Greenland (PR 29), cannot say what an observation of
hurricane-force wind there would show (no station sees one), and the fix-level classification covers only
HF-strength fixes of tracks that are events, so a track's "barrier-only" label depends on the 00/06/12/18 UTC
fixes pipeline A keeps. There is no terrain mask in pipeline A, and Atlantic fixes north of 60N are
lower-confidence throughout. 22 seasons is the sample; a held-out look is not used.

## Deviations (post hoc)

None yet.
