# Does the NAO x PNA shift of Atlantic storm longitude replicate on 1979-2000? (RA-8; ERA5 proxy)

Plan, committed before any 1979-2000 position was related to an index: [PREREGISTRATION.md](PREREGISTRATION.md) (commit `08dbc0d`).
Code: `common.py`, `power.py`, `run.py`. Results: `results/` (`summary.txt`, `results.csv`, `plan_tests.csv`, `power.txt`).
Reproduce (about 2 minutes, no pull): `python3 power.py results/power.txt 400; python3 run.py results`.

Everything is an **ERA5 proxy**. Cyclones are the **MSLP-only tracker run of PR 38**, which reuses **pipeline A's** detector and linker (not pipeline B). The position is the fix of **minimum pressure**,
not the gust peak the lead used. The table covers **Dec-Mar only** (the lead used Oct-Apr).

## Answer, in plain words

**It does not replicate.** On the 22 winters 1979-80..2000-01 the interaction goes the other way: **-0.68 deg of longitude per SD x SD** (95% interval -1.54 to +0.16), against the lead's +1.11.
The interval excludes the lead at face value and also excludes the shrunken expectation of +0.55. It is **not a well-powered null**, though: the smallest shift the sample could detect at 80% power is 1.29 deg,
so a shift of up to 1.5 deg either way cannot be ruled out, and the sample cannot say the effect is zero. In the pre-registered words: *does not replicate at face value, and cannot tell if the true effect is small.*

Read it with the bridge: on the source era with this same position definition and season (2001-2013, Dec-Mar, not independent of the lead), the interaction is **+0.63 deg** (interval -0.14 to +1.49, one-sided p 0.055).
So the pressure-minimum position carries about half of the lead's size in the seasons where it was found, and the pre-2001 sign is opposite. Taken together: the lead is weak where it was found, and absent or reversed before 2001.

## Primary result (P1) and secondaries

Per SD x SD of lagged NAO (days -10..-4) times lagged PNA, degrees. Interval: winter-pairs bootstrap. p: sign-flip score test with winters as clusters; one-sided where the sign was predicted.
q: Benjamini-Hochberg over the 7 planned tests.

| id | what | tracks | winters | gamma (deg) | 95% interval | p | q |
|---|---|---|---|---|---|---|---|
| **P1** | Atlantic longitude, 1979-2000 | 4,364 | 22 | **-0.68** | -1.54 to +0.16 | 0.92 (one-sided, wrong way) | 0.96 |
| S1 | Atlantic latitude | 4,364 | 22 | -0.15 | -0.39 to +0.23 | 0.35 (two-sided) | 0.81 |
| S2 | joint lon+lat | 4,364 | 22 | | | 0.34 | 0.81 |
| S3 | basin taken at the first fix, longitude | 4,405 | 22 | -0.70 | -1.48 to -0.00 | 0.96 (one-sided) | 0.96 |
| S4 | same-time NAO/PNA (days -3..+3), longitude | 4,364 | 22 | -0.39 | -1.31 to +0.55 | 0.81 (one-sided) | 0.96 |
| S5 | **bridge** 2001-2013, longitude (source era; not independent) | 2,743 | 13 | +0.63 | -0.14 to +1.49 | 0.055 (one-sided) | 0.39 |
| S6 | Pacific control, longitude, 1979-2000 | 5,740 | 22 | -0.23 | -0.84 to +0.40 | 0.47 (two-sided) | 0.82 |

0 of 7 pass FDR. The one that comes nearest is the bridge (q 0.39).

## Power

Season-block permutation of the index series (400 permutations, outcome marginals only; `results/power.txt`): null SE of gamma (longitude) is **0.46 deg for P1**, so the smallest shift detectable at 80% power is
**1.29 deg** (bridge 1.57, Pacific 1.08). A well-powered null (interval inside +/-0.55) was not achievable at the planned sample; this was learned from `power.py` before P1 was fitted. The effective sample is 22 winters,
not 4,364 tracks. The lead's nominal +1.11 was therefore detectable at about 67% power, and a shrunken +0.55 at about 22% (normal approximation, 5% two-sided).

## What this does and does not show

- **A different position and a different season window.** A failure is partly ambiguous (agenda pitfall). The bridge says the pressure-minimum definition retains about half the lead's size in 2001-2013 but with an interval that includes 0, so the
  position change alone does not explain the reversal; it is also not enough to prove the original lead is real.
- **Oct, Nov and Apr are not in the committed table**, and 2014-2025 was left out of the bridge on purpose (no 2015-25 look spent).
- **Pre-2001 pressure fields carry no gust drift** (Decision 1 concerns gust), so the pressure position is allowed back to 1979. ERA5 was less constrained by observations early on, which can blur track positions.
- No tropical-cyclone transitions were removed; no hurricane-force label is used (barred before 2004).
- **Looks:** 6th look at 1979-2000. No 2015-25 season scored.

## Numbers independently recomputed

(see the verification note below, filled in after the check)
