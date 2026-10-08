# Is late HF wind in the high-latitude Atlantic barrier flow behind the low? (RA-26)

ERA5 **proxy**, **pipeline A** (`research/era5/hf_history`, 800 km ocean gust index, threshold 71.7 kt). Seasons 2004-05 to 2025-26, non-TC events, Atlantic events whose
onset-fix gust maximum is at or north of 60N (n = 283: 151 late, 132 non-late). Plan committed first: `PREREGISTRATION.md` (`eb0d5a5`). Nothing was pulled. No terrain mask: terrain-flagged
fixes (gust maximum within 100 km of Greenland, 145 of the 283) are **kept** in the primary run and removed in S5.

## Answer in plain words

**No, and the test could have found it.** In late high-latitude Atlantic storms the gust maximum is *not* more often in the "barrier" position (within 300 km of land and more than 400 km from the centre):
**38.4% of late against 39.4% of non-late storms** (difference -1.0 points, 95% interval -11.9 to +12.0). The agenda predicted at least 60% against under 35%. The study had 98% power
to see the predicted 25-point gap and could detect 17 points at 80% power, so this is a well-powered null for *this pattern*, not an underpowered miss.
The slow-deepener alternative is not favoured either: the maximum lies within 400 km of the centre in 61% of late and 58% of non-late storms (+2.6 points, [-10.1, +13.4]).

So the threefold odds of lateness north of 60N (RA-7, PR 86) is **not explained by where the gust maximum sits relative to the coast and centre** in the proxy. What does differ is mild:
nearly all maxima in this latitude band are within 300 km of land (98% of late, 89% of non-late; +8.6 points, nominal p 0.008, q 0.08, fails FDR), so the coast criterion does not discriminate
and the pattern's share is set by the centre distance, which is the same in both groups.

## All tests (one family of 10, BH-FDR; 0 of 10 pass q < 0.05)

Shares are late / non-late; difference with a 2,000-draw season-block bootstrap interval (22 seasons); `results/tests.csv`, `results/late_highlat.txt`.

| test | n late / non | late | non-late | difference [95%] | p | q |
|---|---|---|---|---|---|---|
| **P1 barrier pattern (primary)** | 151 / 132 | 0.384 | 0.394 | -0.010 [-0.119, 0.120] | 0.85 | 0.94 |
| S1 coast <= 300 km | 151 / 132 | 0.980 | 0.894 | +0.086 [0.027, 0.149] | 0.008 | 0.08 |
| S2 centre distance > 400 km | 151 / 132 | 0.391 | 0.417 | -0.026 [-0.134, 0.089] | 0.61 | 0.92 |
| S3 wind from the north (315-45 deg) at the maximum | 150 / 132 | 0.540 | 0.545 | -0.005 [-0.120, 0.113] | 0.94 | 0.94 |
| S4 barrier, LATE6 (strictly after the minimum, >= 6 h) | 66 / 217 | 0.455 | 0.369 | +0.086 [-0.048, 0.220] | 0.20 | 0.67 |
| S5 barrier, terrain-flagged fixes removed | 68 / 70 | 0.338 | 0.357 | -0.019 [-0.197, 0.170] | 0.83 | 0.94 |
| S6 barrier at the peak-gust fix | 151 / 132 | 0.404 | 0.348 | +0.055 [-0.069, 0.188] | 0.38 | 0.92 |
| S7 maximum within 400 km of centre | 151 / 132 | 0.609 | 0.583 | +0.026 [-0.101, 0.134] | 0.64 | 0.92 |
| S8 barrier, Atlantic south of 60N (contrast) | 202 / 554 | 0.129 | 0.065 | +0.064 [0.011, 0.137] | 0.017 | 0.085 |
| S9 adjusted odds ratio, barrier (marginal storm, speed, month) | 283 | | | OR 0.85 [0.49, 1.47] | 0.56 | 0.92 |

Decision rule (registered): P1 interval upper bound 0.12 < 0.15 and power at a 25-point gap 0.98, so **"no (well powered)"**. Power is the normal approximation from the bootstrap standard error (not a
planted-effect simulation). Minimum detectable difference at 80%: 0.17 (P1), 0.19 (S4, the smallest sample), 0.27 (S5).

## What this does and does not say

- It rejects "late means barrier flow behind the departed low" as the explanation of PR 86's OR 3.1, in the proxy and in the way I defined barrier flow (coast of any land, centre distance, wind from the north).
  It does not say what the explanation is. PR 86's joint model leaves OR 1.82 for the 60N group after marginal storms, slow translation and deepening are accounted for.
- S1 and S8 lean the same way (late storms have the maximum nearer land), but neither passes FDR and S1 is nearly saturated.
- Caveats: ERA5 under-resolves the gust peak; the proxy near Greenland has the highest false-alarm ratio (0.44); the coast is the nearest land of any kind; the wind direction is the 10 m wind at the gust maximum, 6-hourly,
  so a short-lived jet is sampled coarsely. 63% of LATE storms are same-fix ties with the pressure minimum (S4 is the strict version, underpowered for a small gap).
- No post hoc analyses; the "Deviations" section of the plan is empty.

## Looks and verification

- Held-out looks: one descriptive use of all 22 seasons, including 2015-25, counted as **look #16** at 2015-25 (agenda count); zero looks at pre-2001. Entry in `../looks/late_highlat.log` (own file, per the coordinator, to avoid log conflicts).
- Verification: see `verify/VERIFICATION.md`.

## Reproduce

`python3 -I late_highlat.py` (needs numpy, pandas, scipy, statsmodels; reads committed files; bootstrap seed 7).
