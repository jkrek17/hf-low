# Greenland tip jets: typical setups, and what goes with hurricane force rather than storm force

**ERA5 proxy throughout.** Gust and winds are reanalysis, not observation. Gust: the shared ERA5 gust cache (`gust_nh_12h`), ocean cells within about
216 nm (400 km) of Cape Farewell. Lows: this study's own ocean MSLP minimum (50-70N, 60-10W), neither pipeline A nor B. November to March,
2004-05 to 2025-26 (decision 1), 00 and 12 UTC. Plan: `PREREGISTRATION.md`, committed first (`36bbf57`; setups changed to Jason's three before the pull).
Pull: 33.6 GB, 2,683 times (2,214 events plus 469 comparison times), saved at `/mnt/project-files/tip-jets/fields_2683_times.tar`.

**Event:** the largest ERA5 gust in that area is 54 kt or more ("storm-force proxy") or 71.7 kt or more ("hurricane-force proxy", the project's count-matched
threshold; ERA5 71.7 kt is about 62 kt in ASCAT, `ascat-check`). 6,654 times, of which 2,214 reach 54 kt and 218 reach 71.7 kt, linked within 36 h into
**514 episodes (126 reach the hurricane-force proxy, about 1 in 4)** across all 22 seasons. No station check was made here; in the earlier Greenland jets
check (PR 29) no Prins Christian Sund report reached 64 kt. "Hurricane force" below always means this proxy class.

## The typical setups (composites: `results/composites_s2.png`)

| Setup (Jason's three; the first was widened post hoc, see below) | Episodes | Hurricane-force proxy share | Composite low, hurricane-force cases |
|---|---|---|---|
| Westerly flow, gust maximum east to south of the cape (the forward tip jet; 122 of the 158 are east of the cape) | 158 | 31% (49) | about 63N 36W, northeast of the cape; westerlies below it |
| Northerly flow, east of the tip | 73 | 25% (18) | about 58N 39W, southeast of the cape; north wind on its west side |
| Easterly flow, west of the tip (ahead of a front) | 97 | 27% (26) | about 57N 55W, in the Labrador Sea; easterlies ahead of it |
| None of the three | 186 | 18% (33) | mostly north of the tip in northerly flow, or south in northerly flow |

Hurricane-force composites have a low near 967-974 hPa against 983-985 for storm force. About 3 in 10 episodes in each of the three setups reach the
hurricane-force proxy against about 2 in 10 outside them (differences between setups not tested). The three setups together are spread over November to March
(59 to 72 episodes a month, `results/extras.txt`).

## What goes with hurricane force rather than storm force (post hoc, all seasons, same-time associations)

Odds ratio per standard deviation (SD), hurricane-force against storm-force proxy episodes, adjusted for how deep the low is (except depth itself), 95%
season-block intervals, `results/posthoc_s2.txt`. **Bold** = BH q < 0.05 across the 15. This is a post hoc lead, not a pre-registered result.
One SD: depth 13 to 14 hPa of central pressure, gradient about 1.3 hPa per 100 km, ice-sheet pressure about 10 hPa, cold-air term 2.3 to 2.9 K.

| Factor | Westerly | Northerly (east) | Easterly (west) |
|---|---|---|---|
| Deeper low | **5.3 [3.7, 9.6]** | **2.0 [1.2, 4.3]** | **2.9 [1.8, 5.6]** |
| Sharper pressure gradient at the wind maximum | **7.4 [4.1, 17]** | **1.8 [1.2, 3.1]** | **4.0 [2.3, 8.6]** |
| Higher pressure over the Greenland ice sheet | **0.66 [0.41, 0.95]** | **2.1 [1.5, 3.1]** | 1.8 [1.0, 3.4], q 0.066 |
| Colder air than sea upstream (SST minus air temperature) | **2.4 [1.4, 4.5]** | 0.9 [0.3, 1.8] | 1.0 [0.3, 1.8] |
| Low farther from the cape | 0.70 [0.44, 1.05] | 1.1 [0.7, 1.9] | 0.8 [0.5, 1.3] |

In plain terms (all post hoc, all about the hurricane-force proxy):
- **A deeper low goes with the hurricane-force proxy in every setup**, about 2 to 5 times the odds per 13 hPa. Depth and gradient are measured at the same time as the gust, so they describe the storm; they are not an advance warning. The gradient is largely the wind itself (rank correlation with gust 0.81 in the westerly setup).
- **Higher ice-sheet pressure goes with the proxy in the northerly setup** (2.1 per 10 hPa) and possibly the easterly one (1.8, q 0.066, cannot tell). In the westerly setup it goes the other way (0.66). The ice-sheet pressure may partly reflect the low's own field when the low sits close (not tested; the sea-level reduction over the ice is another possibility). The northerly result agrees in direction with the earlier Greenland jets study (PR 58), which found about 4 per SD for coastal cases that were mostly barrier flow.
- **Colder air over the sea upstream** goes with the proxy in the westerly setup (2.4; 3.0 [1.7, 5.7] with calendar month added; both halves of the seasons agree in direction, 3.8 and 1.7). It appears only with the setup widened (1.2, not significant, in the pre-registered class). In the other two setups there is no detectable lean (0.9 and 1.0; effects below about 2 cannot be detected). The upstream box (56-63N, 60-46W) may include sea-ice cells; not checked.
- **Distance of the low from the cape**: cannot tell (0.70 [0.44, 1.05] in the westerly setup, 0.43 [0.17, 0.97] in the original class).

## OPC link (`results/opc_link_s2.txt`)

An OPC hurricane-force low fix (Atlantic archive) lies within 600 nm of the gust maximum within 12 h in **92% of westerly, 94% of northerly and 96% of
easterly hurricane-force proxy episodes**, against 24%, 36% and 54% for storm-force proxy episodes. This is "OPC agreed" through an OPC low (usually the parent;
not checked that it is the same low), not a check of the tip jet itself; OPC archives low fixes only.

## What could not be tested, what is post hoc, what was looked at

- The pre-registered fit-and-hold-out tests **did not run**: no class had 30 of the rarer class in each half (hurricane-force episodes in the original classes: 9, 18, 26 in total; per half 3/6, 6/12, 17/9; even the widened westerly class, 23/26 per half, would fail). `results/tests.csv` is kept as the record and **must not be quoted**.
- Everything above is post hoc. The westerly class was **widened after the pull** (bearing 45-225 instead of 135-225) because the original left the classic forward tip jet (122 episodes) in "none of the three". The choice used the class cross-table and the class hurricane-force shares (outcome, about 25% each), and `analyse.py` had already written half-split estimates for the original classes. Looks: original-class half-split (both directions), original-class all seasons, widened all seasons, widened even/odd with a within-season shuffle, one month-adjusted refit. No held-out data remain; confirmation needs 2026-27 seasons or later. Both versions are shown (`results/posthoc_orig_classes.txt`).
- Detectable size (analytic, episodes treated as independent, no adjustment for depth, so a lower bound): about 1.6 (westerly), 2.1 (northerly), 1.9 (easterly) times per SD at 80% power.
- Comparison times: 469 (not the 700 planned, because comparison times that were also events were dropped); the composite anomaly is against their mean.
- Low-level stability: not tested (no pressure-level temperature pulled). Upper-air composites: not built. The whole sample is in the high-latitude Atlantic near steep terrain, where the proxy is least certain (PR 29); northerly-setup maxima often sit next to land.
- Even/odd season split (`results/audit.txt`): depth, gradient and northerly Greenland-high effects keep their direction in both halves; the westerly cold-air effect is 3.8 and 1.7; the easterly Greenland high is 1.3 and 2.4.

## Checks

Recomputed by a fresh Sonnet agent from the committed files: the 6,654 / 2,214 / 218 counts, 514 / 126 episodes, class counts, wind direction and bearing for 15 episodes from the raw fields, OPC shares, the odds ratios in the table for depth, gradient, Greenland high and cold air, and the 0.81 correlation: all matched. **Not independently checked:** bootstrap intervals and q values, detectable sizes, the even/odd and shuffle audit, composite low positions, the monthly counts and SDs in `extras.txt`, the month-adjusted refit, the 33.6 GB figure. A claim review (`hf-reviewer`) corrected eleven wording and figure errors in an earlier version of this file.
