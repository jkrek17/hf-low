# Do the share results hold for sustained HF lows? (RA-22, pre-registration)

Written 2026-10-08 before any 2-fix or 3-fix outcome was related to any predictor. Everything ERA5 is a **proxy**
(pipeline A, `research/era5/hf_history`, 800 km ocean gust index, HF-equivalent at 71.7 kt); the archive is the OPC HF
archive on `main`. Deviations after the first 2-fix or 3-fix result go under "Deviations (post hoc)" at the end, with
both versions shown.

## Question

About a quarter of archive HF events are one 6 h HF fix; requiring at least 2 fixes takes the 2004-05 to 2025-26 counts
from 46.0 to 33.8 (Atlantic) and 38.9 to 29.1 (Pacific) per season (`count_reconcile`). Do the headline ratios of
PR 14, PR 41, PR 63 and PR 64 keep sign and size when "HF" means at least 2, or at least 3, HF fixes? And does the
teleconnection or hemispheric pattern raise the number of sustained HF lows more than the number of brief ones?

Plain-words answers to give: for each headline, yes / no (well powered) / can't tell, and by how much; and overall
"favours sustained", "favours brief", or "no clear preference".

## What was looked at before writing this

- The published results of PR 14, 41, 63, 64 (their READMEs) and the count reconciliation.
- Label bookkeeping only: the number of HF tracks surviving each rule. Oct-Apr 2004-05 to 2025-26, pipeline A: Atlantic
  1,007 HF tracks of 9,656; at least 2 fixes 664, at least 3 fixes 439, exactly 1 fix 343. Pacific 825 of 10,040; 542, 374 and 283.
  Archive (class `low`, season 2004-05 to 2025-26): Atlantic 970 events, 720 with at least 2 HF fixes, 476 with at least 3;
  Pacific 859, 641, 452 (`hfN` in `docs/data/hf-lows.json`).
- **Controls (k = 1) were run and reproduce the originals**, because the plan needs them to be exact before k = 2, 3 are run:
  archive weekly RR per SD 1.236 / 1.193, proxy weekly RR(HF) 1.291 / 1.165 and RR(all) 1.016 / 1.033, daily RR(HF) 1.096 / 1.072 and
  RR(all) 0.974 / 1.028 (Atlantic / Pacific). Those are the numbers already published; no 2-fix or 3-fix outcome has been joined to
  any predictor.

## Definitions

- **HF fix.** An in-domain 6-hourly fix with 800 km gust index at or above 71.7 kt (proxy; `lifecycle.py` rule). Archive: a fix of category `HF`.
- **HF_k** (k = 1, 2, 3). A pipeline A track whose peak gust is at least 71.7 kt **and** that has at least k HF fixes in total
  (`n_hf`, `results/track_nhf.csv`, from `era5_hf_catalog_tracks.csv`, agreeing with `lifecycle_events.csv` for all 4,157 events);
  for the archive, an event of class `low` with `hfN >= k`. k = 1 is the original definition, except that the archive k = 1
  keeps the 4 of 1,829 class-`low` events in 2004-2025 that have no HF-category fix, as PR 41 did (so it reproduces PR 41 exactly).
- **Brief** = exactly 1 HF fix. **Sustained** = at least 3 HF fixes (12 h or more at HF-equivalent). Exactly 2 fixes is the middle class.
- A track that fails HF_k **stays a cyclone in every denominator**: only the HF flag changes. Cyclone counts (RR(all), T1) do not depend on k.
- Duration is counted in total HF fixes, not consecutive ones (archive `hfN` is also a total). A consecutive-run version is secondary S1.
- Not repeated, because the depth cut has no duration: the Pacific depth variants (PR 14 S4, PR 64 S1) and the 1979-2000 within-era checks (PR 41 S5, PR 64 S7).
  The OPC area rule is left out (its edges were recalled, not sourced).

## Re-runs (the originals' own code, unchanged except which tracks are called HF)

| Study | Setup, unchanged | Command | Headline ratios |
|---|---|---|---|
| PR 14 (`freq_split`) | daily genesis counts, Oct-Apr, lagged CPC NAO (Atlantic) / PNA (Pacific), days -10..-4, 22 seasons | `split.decompose`, primary only | RR(HF_k), RR(share_k) = RR(HF_k)/RR(all), f |
| PR 41 (`hemispheric`) | archive weekly counts by first-fix date, leave-one-season-out pattern index (`oos_index_*.csv`), as `hem_channels` S5 | `contrast.py archive` | RR per SD of the index for HF_k |
| PR 63 (`nao_share_barrier`) | Atlantic NAO share, barrier-type fixes removed, Greenland-high mediator | `barrier_share.py` on the k-tracks file | RR(share_k); P1 RR(share_k without terrain-type fixes); P3 RR(NAO share_k given the Greenland high); F4 RR(Greenland-high share_k per SD) |
| PR 64 (`hem_channels`) | weekly pipeline A counts, same pattern index | `channels.py` on the k-tracks file | T3 RR(HF_k), T2 RR(share_k), f, share part of the change |

PR 41's pattern **cannot be refitted** (it needs the 17 GB of WeatherBench2 fields, and no pull is allowed here). The index used
is the one fitted to archive HF counts at k = 1 and is out of sample for each season. A weaker response at k = 2 or 3 could therefore
partly reflect that the index was fitted to the k = 1 label; this is a stated limit, not tested.

New contrast (`contrast.py`, the same three setups): for each basin, log RR(sustained) minus log RR(brief) per SD of the index, plus
the paired differences log RR(HF_k) minus log RR(HF_1) for k = 2, 3. Season-block bootstrap interval (2,000 draws), season-block
permutation p (2,000 draws, two-sided), season-clustered SE; detectable RR at 80% power = exp(2.8 x SE).

## Predictions (written before the outcomes)

- **P1.** Every headline ratio keeps its k = 1 sign at k = 2 and k = 3 (agenda prediction).
- **P2.** The pattern effect is at least as large for sustained HF as for brief: the contrast log RR(sustained) - log RR(brief) is not negative.
- **P3.** The Greenland-high mediation of the Atlantic NAO share effect (PR 63) survives: RR(NAO share_k given the high) is closer to 1 than RR(NAO share_k).

## Decision rules (thresholds, fixed now)

For each headline cell (study x basin x k in {2, 3}):
1. **Holds** if the sign equals the k = 1 sign and BH q < 0.05 over the whole family below. **Lost, well powered** if q >= 0.05 and the detectable RR (exp(2.8 SE)) is no larger than the k = 1 estimate, so the original effect would have been seen. **Can't tell** otherwise.
2. **Size relative to k = 1:** the paired difference log RR(HF_k) - log RR(HF_1) with its bootstrap interval: **larger** if the interval is above 0, **smaller** if below 0, **same within noise** otherwise (interval half-width reported).

For the question "sustained or brief?" (one verdict per study-basin, six cells: PR 14 x 2, PR 64 x 2, PR 41 archive x 2):
3. **Favours sustained** if the contrast is above 0 with q < 0.05 (BH over the six). **Favours brief** if below 0 with q < 0.05. Otherwise **no clear preference**, with the detectable contrast stated.
4. **Overall:** favours sustained if at least 4 of 6 cells favour sustained and none favours brief; the mirror for brief; otherwise no clear preference, with "leans sustained / leans brief" only if at least 4 of 6 point estimates share a sign.
5. Brief lows have lower peak gust by construction, so a contrast above 0 may be intensity, not duration. **S2** (below) is the control. If S2 gives the same contrast with a count-matched gust cut, the answer is "stronger storms", and is reported that way.

## Multiplicity

- **Family A (headline cells, BH):** PR 14 RR(HF_k) and RR(share_k) x 2 basins x 2 k = 8; PR 64 T3 and T2 x 2 x 2 = 8; PR 41 archive RR x 2 basins x 2 k = 4; PR 63 P1, P3, F4_GH_share_orig x 2 k = 6 (Atlantic only; `orig_share` is the PR 14 Atlantic share and is not counted twice). **26 tests.** p values are the originals' own season-block permutation p (PR 63 P3 uses its bootstrap p, as in the original).
- **Family B (contrasts, BH):** the six sustained-minus-brief contrasts. The paired differences HF_k - HF_1 (24 estimates) are descriptive, reported with intervals and p, and enter a third family (BH over those 24).
- q over all tests (A + B) is reported next to the within-family q. Decision rules use the A and B within-family q as stated above.
- Secondary S1-S3 are labelled secondary and carry their own BH family; none changes a primary verdict.

## Secondary checks

- **S1 consecutive run:** HF_k by longest run of consecutive HF fixes (`all_tracks_run{k}`), PR 14 and PR 64 headline ratios.
- **S2 count-matched gust cut:** the same number of tracks as HF_k, taken as the highest peak gust in each basin over 2004-05 to 2025-26 (`all_tracks_gcut{k}`; the cut depends on counts only). PR 14 and PR 64 headline ratios. If duration matters beyond intensity, HF_k and S2 differ.
- **S3 middle class:** exactly 2 HF fixes, in the contrast tables.

## Power and what this cannot show

- Sustained events are 44% (Atlantic) and 45% (Pacific) of the original HF tracks, so intervals widen about 1.3 to 1.5 times at k = 3 (sqrt of the count ratio). Detectable RR at 80% power is reported for every cell; a "no" needs the detectable RR to be no larger than the original effect.
- The independent sample is 22 seasons (effective n about 22 x 5 index values per season), not tracks.
- Not a cause, and not a new search: this is a robustness check of four existing findings, using the same 22 gust-era seasons and archive weeks again. It adds no new data.
- ERA5 proxy: the proxy's HF duration is not the archive's (proxy: 66% of HF tracks have at least 2 fixes; archive 74%).
- Duration and peak intensity are linked; S2 separates them only partly.
- The Atlantic fixes north of 60N are a lower-confidence group in the originals and remain in the k tracks. PR 63's terrain-type rule is applied as in the original (a track needs one non-terrain HF-strength fix), not re-defined per k.

## Looks at the data (to be logged in `hemispheric/results/heldout_looks.log`)

This thread spends **one look at each of the 22 gust-era seasons for pipeline A** (PR 14, 63 and 64 designs, all with the same
k labels) **and one at the 2004-05 to 2025-26 archive weekly counts** (PR 41 design). The 2015-16 to 2025-26 seasons are
among them, so it is counted as look 7 at the 2015-25 seasons (leave-one-season-out index; no model is refitted and no held-out
split is scored). Pre-2004 seasons are not used. The k = 2 and k = 3 results are each run once; a rerun forced by an error, with the
outcome definition unchanged, will be logged. Nothing is re-picked after the results.

## Deviations (post hoc)

None yet.
