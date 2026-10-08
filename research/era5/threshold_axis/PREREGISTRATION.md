# Do the share results move with the gust threshold? (RA-23, pre-registration)

Written 2026-10-08, before any outcome was joined to a predictor at 68 or 75 kt. Everything is the ERA5 **proxy**, pipeline A
(`research/era5/hf_history`, 800 km ocean gust index; HF-equivalent at 71.7 kt in the originals). Departures after the first 68 or 75 kt
outcome go under "Deviations (post hoc)" at the end.

## Question

Pipeline A's cut of 71.7 kt was calibrated on five seasons. If the share results of PR 14, 63 and 64 are sharp functions of that cut, they
describe the threshold, not the storms. Also (asked by the coordinator): is RA-22's lean toward sustained lows just a higher gust cut?
Plain-words answers: for each headline, "holds at 68 and 75 kt" / "collapses or flips" / "can't tell", and whether the effect grows with the cut.

## Fixed in advance

- **Thresholds:** 68.0, 71.7 (control, reproduces the originals) and 75.0 kt on the track's peak gust index `gust800_kt`. No other values will be run or reported in place of these.
- **HF_T:** a pipeline A track with peak gust at or above T. Tracks below T stay cyclones in every denominator (as in RA-22). Duration is not used (k = 1).
- **Setups (the originals' code, unchanged except which tracks are HF):** PR 14 daily NAO (Atlantic) / PNA (Pacific), lagged days -10..-4, Oct-Apr 2004-05..2025-26;
  PR 64 weekly pattern index (leave-one-season-out, PR 41), same seasons; PR 63 Atlantic NAO share, terrain-type exclusion (P1), mediation by the Greenland high (P3), Greenland high on share (F4).
- **Not run, and why:** the archive (PR 41) has no threshold. PR 63 P1 is run only at 71.7 and 75 kt: its terrain classification comes from `highlat/gustloc_fixes.csv`, which holds only fixes at 71.7 kt or above, so at 68 an HF track would have no rows and be kept as unclassified HF, which biases P1 towards the original share. P(HF) as a soft label is not run (per-cyclone P(HF) exists only for catalog events). Pre-2004 seasons are barred (Decision 1).
- **Bookkeeping already seen (counts only):** all-months 2004-05..2025-26 tracks at or above 68 / 71.7 / 75 kt: Atlantic 1,589 / 1,104 / 772, Pacific 1,399 / 896 / 621, out of 16,992 and 18,464 tracks. Control to be run: T = 71.7 must reproduce PR 14 (1.096, 1.125, 1.072, 1.043), PR 64 (T3 1.291 / 1.165, T2 1.271 / 1.129) and PR 63 (P1 1.141, P3 0.923, F4 0.870) exactly.

## Headline tests (family A, BH over 21)

- PR 14: RR(HF_T) and RR(share_T) x 2 basins x T in {68, 75} = 8 (p: the original's season-block permutation p).
- PR 64: T3 RR(HF_T), T2 RR(share_T) x 2 basins x {68, 75} = 8.
- PR 63 (Atlantic): P3 (mediation) and F4 at {68, 75} = 4; P1 at 75 = 1.

## Predictions and decision rules

- **P1 (sign and rank).** Every headline ratio keeps its 71.7 sign at 68 and 75, and the Atlantic share RR stays above the Pacific one in PR 14 and PR 64 at each T. **Holds** per cell if same sign and q < 0.05. **Collapse or flip** if the sign changes, or q >= 0.05 while the detectable RR (exp(2.8 SE)) is no larger than the 71.7 estimate. Otherwise **can't tell**.
- **P2 (growth with the cut).** log RR(share_75) - log RR(share_68) is above 0. Judged on the paired bootstrap interval in PR 14 and PR 64 (family C: paired differences t75-t68, t75-t717, t717-t68 x 2 basins x 2 setups, HF and share are identical differences because RR(all) does not depend on T, so 12 distinct; BH over those 12). "Grows" if the interval is above 0, "flat" if it contains 0, "shrinks" if below.
- **P3 (the RA-22 lean).** If the RR of HF_T rises with T by about as much as it rose with duration in RA-22 (HF_1 to HF_3: PR 64 Atlantic 1.291 to 1.385, Pacific 1.165 to 1.264; PR 14 1.096 to 1.137 and 1.072 to 1.152), the lean toward sustained lows is consistent with an intensity gradient and nothing more. Compared descriptively with the RA-22 count-matched gust cut results (already committed); the gust cut in kt that matches the HF_3 counts is reported.
- Power: detectable RR at 80% (exp(2.8 x clustered SE)) is stated for every cell. Counts fall about 30% at 75 kt (Oct-Apr similarly), so intervals widen about 1.2 times.

## Looks

One look at the 22 gust-era seasons (pipeline A, PR 14 / 63 / 64 setups; leave-one-season-out index, nothing refitted), counted as look 14 at the 2015-16..2025-26 seasons by the agenda thread's count (13 scorings before this one). Each of 68 and 75 is run once. It is a robustness check of existing findings, not a new search, and its seasons are the same ones as before: it is not an independent confirmation. Effective n is 22 seasons.

## What it cannot show

Not a cause. ERA5 reads low in extreme storms, so the "true" HF cut may lie above 71.7 kt in ERA5 units; a stable effect across 68 to 75 kt says the mechanism acts on storm strength broadly, not that 71.7 is right. The thresholds span only about 10% of the index range. Peak-gust cuts say nothing about duration (the fix-level gust table covers events at 71.7 kt or above only).

## Deviations (post hoc)

None yet.
