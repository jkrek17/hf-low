# RA-12: are close HF pairs parent-daughter, or a shared environment?

Plan: `PREREGISTRATION.md` (commit `8fe1f08`, before any p-value). Code: `common.py`, `engine.py`, `count_pairs.py`, `power.py`, `run.py`. Output: `results/pair_families.json`, `results/run.log`. Everything is **pipeline A** (`research/era5/hf_history`), an ERA5 **proxy**, 2004-05 to 2025-26, transitioning tropical cyclones left in. Reproduce: `PYTHONPATH=. python3 -s count_pairs.py; PYTHONPATH=. python3 -s run.py 2000` (about 6 minutes).

## Answer in plain words

**Parent-daughter pairs do not dominate, in either basin; a smaller share cannot be excluded, and for the Atlantic the data cannot tell a mixture from none.** Among HF lows born within 1,500 km of a live HF low and within 48 h of its peak, the share that start in the trailing sector (within 60 degrees of directly behind) is 0.44 in the Atlantic (70 of 160 pairs) against 0.42 expected when storms are shifted in time within their season-month (p = 0.29, BH q = 0.59). The Pacific share is 0.30 (18 of 60) against 0.50 expected: if anything fewer daughters sit behind the parent than chance (two-sided p = 0.0005, not a pre-registered direction). The Atlantic excess of rear-sector pairs is +3 (95% upper bound +17), against PR 37's 42-pair spatial excess; the Pacific excess is -12 (upper bound -4), against PR 37's 24. So families explain at most about 40% of the Atlantic persistence and about none of the Pacific, and the remainder is either shared environment or something this design cannot see. It is not a demonstration of shared environment (see the plan's decision rule).

**Power.** The test reaches 80% power only at +32 extra rear pairs in the Atlantic (76% of PR 37's excess) and +24 in the Pacific (all of it). It therefore rules out "families dominate", not a mixture of the size the agenda expected (a quarter to a third). Only 15 Atlantic and 9 Pacific candidate pairs have an observed genesis, so no claim rests on that subset.

## Results

| | Atlantic | Pacific |
|---|---|---|
| HF tracks | 1,104 | 896 |
| Candidate pairs / rear sector | 160 / 70 | 60 / 18 |
| Rear fraction, observed / null mean (null SD) | 0.438 / 0.418 (0.034) | 0.300 / 0.502 (0.050) |
| Primary p, one-sided; BH q (family of 2) | 0.295; 0.59 | 1.00; 1.00 |
| Two-sided p (not decision-bearing) | 0.58 | 0.0005 |
| Excess rear pairs (95% upper bound) | +3.1 (+16.8) | -12.1 (-4.1) |
| Minimum detectable excess, 80% power | 32 pairs | 24 pairs |
| Observed-genesis daughters: rear / n; null mean fraction | 10 / 15; 0.53 | 5 / 9; 0.56 |
| Same, one-sided p / two-sided p (secondary, uncorrected) | 0.014 / 0.033 | 0.55 / 0.92 |

Tests: 2 primary (0 pass FDR). With the 2 secondary one-sided tests added to the family (4 tests), none passes: the Atlantic observed-genesis p = 0.014 (15 pairs) has q = 0.058. Two-sided p-values are descriptive. The conditional NAO/PNA test in the plan was not run because no primary test passed.

## Caveats

- **Orientation was seen before the plan was written** (rear counts 70/160 and 18/60 were printed by an early power run). The primary result is a pre-registered computation on a definition fixed after that glimpse, so it is exploratory-leaning, not a clean test. No other sector, distance or time definition was tried.
- The pair-count null is biased: shifted tracks give 197 (Atlantic) and 99 (Pacific) candidate pairs against 160 and 60 observed, probably because the tracker, or real storms, keep new lows from forming beside live ones. Pair counts are therefore not tested; only orientation.
- First detection (first 00/12 fix below 1010 hPa) stands in for genesis; only about a quarter of candidate daughters have an observed genesis. A secondary that stays on the parent's tracker track never appears as a second storm.
- The Pacific rear deficit is a post hoc direction; read it as "no sign of rear daughters", not as a repulsion effect.
- Looks at the 2015-25 block: none (nothing fitted or scored); no looks-log entry.

## Post hoc: heading convention

The plan did not say which two fixes give the parent's heading; the code used the fix at the daughter's time and the one before it (backward). A verifier found every daughter's first fix falls on a parent fix time, so the choice matters: rear-sector pairs Atlantic / Pacific are 70 / 18 (backward, as run), 55 / 17 (forward), 61 / 20 (centred), 67 / 22 (the file's own heading column), of 160 / 60. Only the backward version was tested against the null. The other versions give an Atlantic rear fraction of 0.34 to 0.42 (at or below the 0.42 null) and a Pacific one of 0.28 to 0.37 (below 0.50), so the answer does not depend on it, but the +3 Atlantic excess is the most favourable of the four.

## Verification

See `verify/VERIFICATION.md` (fresh agent recomputing the counts and the null from the committed inputs).

Recomputed independently (`verify/`): tracks, candidate pairs, rear counts, rear fractions, observed-genesis counts and the null mean fraction (within permutation noise) in both basins. Not independently checked: the 2,000-permutation p-values, BH q, two-sided p, power and MDE, the 16 percent statement, and PR 37's excess figures of 42 and 24 (quoted from PR 37's README, 251 vs 209.0 and 135 vs 111.0).
