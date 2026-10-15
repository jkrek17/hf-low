# Greenland tip jets: typical setup, and what leans hurricane force over storm force

Written and committed BEFORE any wind direction, pressure or ingredient was looked at for any event. Only
`stage0_counts.py` has run: it reads the cached ERA5 gust (a PROXY, not observation) and counts times by box
gust. Counts from it: 6,654 Nov-Mar times (2004-05 to 2025-26, 00/12 UTC); box gust >= 54 kt at 2,214, >= 64 kt at
777, >= 71.7 kt at 218.

Related work (hf-low PR 58, `research/era5/greenland_jets/`) asked a different question: which ingredients raise
the odds of a coastal gust >= 71.7 kt anywhere within 300 km of Greenland, barrier flow included. This study
isolates the tip jet at Cape Farewell, separates three setups, and asks storm force versus hurricane
force among tip-jet events. Ice-sheet pressure (PR 58's one supported ingredient) is re-tested here, not assumed.

## Definitions (fixed now)

- **Tip box**: ERA5 0.25 degree ocean cells within 400 km of Cape Farewell (59.77N, 43.9W). **Box gust** = the
  largest 10 m gust in the box at a 00 or 12 UTC time, Nov-Mar, seasons 2004-05 on (decision 1; no earlier seasons).
- **Event time**: box gust >= 54 kt. **Class**: hurricane force (HF) if >= 71.7 kt (the project's count-matched
  proxy threshold); storm force (S) if 54 to <71.7 kt. 54 kt is the nominal proxy storm class used by the
  gust atlas; it is not calibrated on observed storm-force events. ERA5 71.7 kt corresponds to about 62 kt
  ASCAT (ascat-check), and no station near Cape Farewell has reported 64 kt: the HF class is a statement about the
  proxy only.
- **Setup class** (Jason's three common setups, 2026-10-15 04:57Z; replaces the earlier forward/reverse split, changed
  before any pull or ingredient was seen). Bearing of the box-gust maximum from the cape, and 10 m wind direction (the
  direction the wind blows FROM) averaged over the 3x3 cells around it:
  **S** = south of the tip (bearing 135-225) in westerly flow (from 225-315);
  **E** = east of the tip (bearing 45-135) in northerly flow (from 315-45);
  **W** = west of the tip (bearing 225-315) in easterly flow (from 45-135), the setup ahead of a front.
  Anything else is "other", described only. Class E overlaps with the barrier flow of PR 58 by construction: it is
  reported as such, not renamed a tip jet.
- **Episode**: event times linked when 36 h or less apart; the episode is represented by its peak-gust time and
  is HF if any of its times is HF; type from the peak time. **The unit of analysis is the episode, and the
  sample for inference is seasons (22), not episodes or times.**
- **Controls** for composites only: 700 Nov-Mar 00/12 UTC times drawn at random (month-stratified, seed 15) from
  ALL Nov-Mar times (no condition), so anomalies are against the true climate. Anomaly = field minus the control mean.
- **Mask**: lows over the ice sheet are never drawn or counted. The reference low is the minimum MSLP
  over OCEAN cells 50-70N, 60W-10W (no track file, so the ice-sheet reduction artefact cannot create a low).

## Data (pull estimate, stated before pulling)

ARCO-ERA5 0.25 degree u10, v10, MSLP, 2 m temperature, SST at the event times (2,214) plus 700 controls = 2,914
times x 5 variables x about 2.7 MB (stage-1 rate of PR 58) = about 39 GB, under the 50 GB gate. No existing
cache holds these variables at these times (the 4.1 GB cache is gust only; PR 58's work files were not kept).
Upper air (Z500, T850) from WeatherBench2 1.5 degree for composites only (small; ends 2023-01-09, so upper-air
composites use the events up to then and say so).

## Questions and tests

**Q1, typical setup (descriptive, no test).** For the three classes, HF and S separately: composite MSLP
and 10 m wind anomaly maps; the low's median position, depth, and the box wind direction and distance of the maximum
from the cape; month and time-of-day shares; the share of HF episodes with an OPC hurricane-force or storm low
fix within 600 nm in the 12 h either side (OPC has low fixes only, no wind areas, so this is "OPC agreed" through the
parent low, not a check of the tip jet itself).

**Q2, HF versus storm (tested).** Among episodes of one type, which ingredients are higher at HF than at S?
Candidate factors, each standardised, computed at the peak time, none chosen after seeing outcomes:

1. LOW_DEPTH: reference-low MSLP (hPa, sign flipped so deeper is higher).
2. LOW_DIST: distance from the reference low to Cape Farewell (km).
3. GRAD: mean MSLP gradient magnitude (hPa per 100 km) over ocean within 200 km of the gust maximum.
4. GH: mean MSLP over Greenland land points 64-78N, 55W-30W (ice-sheet high, as PR 58).
5. COLD: SST minus 2 m temperature (K), mean over ocean in the upstream box 56-63N, 60W-46W for all three classes
   alike (fixed in advance).

Primary family: 5 factors x 3 classes (S, E, W) = 15 tests, univariate logistic regression of HF
(vs S) on one factor, adjusting for LOW_DEPTH when the factor is not LOW_DEPTH. Even seasons fit, odd seasons held
out, then swapped (as PR 58). 2,000 season-block bootstrap resamples. Benjamini-Hochberg over the 15. A factor is
**supported** if q < 0.05 in the fit half and the same sign with one-sided p < 0.10 in the held-out half, in at least one
direction of the swap; **replicated** if so in both. A class is tested only if each half holds at least 30
episodes of the rarer class; otherwise it is described only, and the minimum detectable odds ratio per SD at 80%
power is reported either way. A well-powered null is a result. Secondary (not in the family, reported as such):
T850 minus SST (stability proxy, WeatherBench2, seasons to 2021-22), and a joint model with all five.

**Stability caveat fixed in advance:** the ERA5 pressure-level temperature profile is not pulled, so low-level
stability is only represented by T850 - SST and by COLD. If those cannot carry it, "stability" is reported as not
tested.

Not tested: NAO (null in PR 58 pre-registered family), low speed and direction (unresolved in PR 58), anything
before 2004-05.

## Reporting rules

Name pipeline: gust is the cache proxy, lows come from this study's own MSLP minimum, not pipeline A. Say proxy
wherever it appears. Nautical miles in text (km in tables). Plain rates for forecasters ("about N in 10"). Every
quoted number recomputed by a fresh Sonnet agent; surprising positives get hf-surprise-audit. Any later change is
logged below as post hoc.

## Post hoc deviations

- 2026-10-15 04:57Z, before the pull: event types changed from forward/reverse to Jason's three setups S, E, W (family 10 -> 15 tests). Not post hoc in the outcome sense: no wind direction had been read.
- 2026-10-15, after the pull (POST HOC): (1) the half-split rule (30 of the rarer class in each half) was met by no class (hurricane-force episodes per class: S 9, E 18, W 26 of 36, 73, 97), so the pre-registered tests did not run; `results/tests.csv` is kept as the record and its numbers must not be read. (2) The pre-registered S sector left the biggest group in `other`: 122 of 308 'other' episodes have the gust maximum east of the cape (bearing 45-135) in westerly wind, the classic forward tip jet. S was widened to bearing 45-225 with westerly wind (S2), chosen from the class cross-table and the class HF shares (about 25% in every class), not from the ingredient results; both versions are reported. (3) All-seasons tests (2,000 season-block bootstraps, ridge 0.5, BH over 15 per version, then 5 pooled) in `posthoc.py`; consistency and shuffle checks in `audit.py`.
- 2026-10-15, after the pull: the control set is 469 times, not 700 (stage1_pull.py drops comparison times that are also events), so the composite baseline excludes event times; the colour bar now says so. Also `extras.py` (monthly counts, SDs, month-adjusted refit in the widened westerly class), added after the claim review.
