# Is PR 86's high-latitude "late HF" odds ratio a domain-edge artefact? (RA-26 follow-up)

ERA5 **proxy**, **pipeline A** (800 km ocean gust index, 71.7 kt), Atlantic, non-TC events, 2004-05 to 2025-26. Plan committed first: `PREREGISTRATION.md` (`d57904c`). No pull.
Terrain-flagged fixes are kept in the primary run and removed in E6.

## Answer in plain words

**No, it is not an edge artefact: the effect stays.** Among Atlantic tracks that never leave the domain through the north or east edge, storms whose onset gust maximum is at or north of 60N are
late (HF at or after the pressure minimum) **54.5%** of the time against **28.8%** south of 60N: **odds ratio 2.97 (95% interval 2.04 to 4.32)**, against 3.14 [2.25, 4.38] in the full sample here
(PR 86 reported 3.10; the 0.04 gap was not traced, probably RA-7's month terms, and the verifier also got 3.14). Tracks that do touch the edge give 3.95 [2.49, 6.26]; the two do not differ detectably (interaction OR 1.33 [0.83, 2.13], q 0.24), so
edge contact may inflate the effect somewhat (point estimate a third higher, interval from none to about double), but not detectably. The registered rule (OR at least 1.76 and the interval above 1) is met with room to spare.

It remains unexplained: RA-26 (PR 94) already ruled out barrier flow, and this rules out the domain edge.

## Why the agenda's mechanism could not act on the minimum itself

`lifecycle.py` takes the lowest pressure over the whole track, in or out of domain, so the minimum is not truncated at 67N. The edge could act only on onset (HF counts only in domain), which is what
the touch definition covers. TOUCH = any fix on the full 6-hourly track that is out of domain at 50N or higher (north or east exit, or the 75N detector limit): 343 of 1,039 Atlantic events (129 of 283 in the
high-latitude group). The buffered version (also any in-domain fix north of 64N or between 7E and 10E) gives 522 (212).

## All tests (one family of 7, BH-FDR; 6 of 7 pass, the seventh is the interaction)

| test | n (high-lat) | late hi / lo | OR [95% CR1] | boot 95% | q | MDE80 | power at OR 2 |
|---|---|---|---|---|---|---|---|
| **E1 not touching (primary)** | 696 (154) | 0.545 / 0.288 | **2.97 [2.04, 4.32]** | [2.05, 4.34] | <0.001 | 1.71 | 0.95 |
| E2 touching | 343 (129) | 0.519 / 0.215 | 3.95 [2.49, 6.26] | [2.62, 6.58] | <0.001 | 1.93 | 0.84 |
| E3 interaction (touch minus not) | 1,039 (283) | | 1.33 [0.83, 2.13] | | 0.24 | 1.96 | 0.82 |
| E4 not touching, adjusted for marginal storm, speed, month | 696 (154) | | 2.21 [1.37, 3.57] | [1.34, 3.52] | 0.002 | 1.98 | 0.81 |
| E5 not touching, LATE6 (strictly after, >= 6 h) | 696 (154) | 0.221 / 0.111 | 2.28 [1.41, 3.66] | [1.44, 3.61] | 0.001 | 1.97 | 0.81 |
| E6 not touching, terrain-flagged fixes removed | 582 (65) | 0.554 / 0.269 | 3.38 [2.23, 5.10] | [2.31, 5.21] | <0.001 | 1.80 | 0.91 |
| E7 buffered definition of touching | 517 (71) | 0.521 / 0.303 | 2.51 [1.40, 4.48] | [1.47, 4.63] | 0.002 | 2.29 | 0.65 |

Full-sample check (T2, all 1,039 Atlantic events): OR 3.14 [2.25, 4.38]. The unadjusted E1 and the adjusted E4 (2.21) bracket PR 86's joint OR 1.82, which had more covariates.
E7 is underpowered by construction (71 high-latitude events). Season-clustered Wald intervals (22 clusters) with a season-block bootstrap beside them; power and minimum detectable OR from the cluster-robust standard error, not a planted-effect simulation.

## Caveats

- This shows the effect is not made by tracks leaving through the north or east edge. It does not make the effect physical: ERA5 6-hourly ties, the false-alarm ratio near Greenland (0.44), and no terrain mask all remain.
- Removing touching tracks removes about half the high-latitude events, so the clean group may differ in storm type.
- No post hoc analyses; the "Deviations" section of the plan is empty.

## Looks and verification

Look #17 at 2015-25 by the agenda count (entry in `../looks/late_edge.log`), none pre-2001. Verification: `verify/VERIFICATION.md`.

## Reproduce

`python3 -I late_edge.py` (numpy, pandas, scipy, statsmodels; bootstrap seed 7).
