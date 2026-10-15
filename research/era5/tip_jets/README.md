# Greenland tip jets: typical setups, and what leans hurricane force over storm force

**ERA5 proxy throughout** (gust and winds are reanalysis, not observation; ASCAT agreement is a consistency check, ERA5 assimilates ASCAT).
November to March, 2004-05 to 2025-26 (decision 1), 00 and 12 UTC. Plan: `PREREGISTRATION.md`, committed first (`36bbf57`, setups
changed to Jason's three before the pull). Pull: 33.6 GB, 2,683 times, saved at `/mnt/project-files/tip-jets/fields_2683_times.tar`.

**Event:** the largest ERA5 gust over ocean within 400 km (about 215 nm) of Cape Farewell is 54 kt or more (storm-force proxy) or 71.7 kt or more
(hurricane-force proxy, the project's count-matched threshold). 6,654 times; 2,214 reach 54 kt, 218 reach 71.7 kt. Linked within 36 h into
**514 episodes (126 hurricane-force proxy, 24.5%)**, all 22 seasons. Nothing reaches a station-confirmed 64 kt there: this is a statement about the proxy.

## The typical setups (composites: `results/composites_s2.png`)

| Setup (Jason's three) | Episodes | Hurricane-force share | Where the composite low is (hurricane-force cases) |
|---|---|---|---|
| South/south-east of the tip, westerly flow | 158 | 31% (49) | about 63N 36W, northeast of the cape; wind runs west to east below it |
| East of the tip, northerly flow | 73 | 25% (18) | about 58N 39W, south-east of the cape; north wind on its west side |
| West of the tip, easterly flow | 97 | 27% (26) | about 57N 55W, in the Labrador Sea; easterlies ahead of it |
| None of the three | 186 | 18% | mostly north of the tip in northerly flow, or south in northerly flow |

Hurricane-force composites have a mean low near 967-974 hPa against 983-985 for storm force. Rates are the same in all three setups, so the setup does not
tell you whether a tip jet reaches hurricane force; the low does. The three setups together are spread evenly over November to March (36 to 48 episodes a month) and S and W occur in all 22 seasons (E in 20).

## What leans hurricane force over storm force (post hoc, all seasons; see "What could not be tested")

Odds ratio per standard deviation, hurricane force against storm force, adjusted for how deep the low is (except for depth itself). 95% season-block
intervals. `results/posthoc_s2.txt`. "Supported" = BH q < 0.05 across the 15.

| Factor | South (westerly) | East (northerly) | West (easterly) |
|---|---|---|---|
| Deeper low | 5.3 [3.7, 9.6] | 2.0 [1.2, 4.3] | 2.9 [1.8, 5.6] |
| Sharper pressure gradient at the wind maximum | 7.4 [4.1, 17] | 1.8 [1.2, 3.1] | 4.0 [2.3, 8.6] |
| Strong high over the Greenland ice sheet | **0.66 [0.41, 0.95]** | **2.1 [1.5, 3.1]** | 1.8 [1.0, 3.4], q 0.066 |
| Cold air over the sea upstream (SST minus air temperature) | **2.4 [1.4, 4.5]** | 0.9 [0.3, 1.8] | 1.0 [0.3, 1.8] |
| Low farther from the cape | 0.70 [0.44, 1.05] | 1.1 | 0.8 |

In plain terms:
- **A deeper low is the main lever in every setup**, about 2 to 5 times the odds per standard deviation (about 11 hPa in the low's centre).
- **The pressure gradient at the wind maximum** is strongly tied to hurricane force, but it is largely the same thing as strong wind (rank correlation with gust 0.81 in the south setup): use it as a description of the jet, not as an advance warning.
- **The Greenland high matters only for the east (northerly) setup**, where it raises the odds about 2 times per standard deviation. This matches the earlier Greenland jets study (PR 58), whose coastal cases were mostly this barrier flow. In the south setup it goes the other way (0.66), but the ice-sheet pressure there partly reflects the low's own field when the low sits close, so do not read it as a cause.
- **Cold air over the sea** leans the south (westerly) setup toward hurricane force (2.4, both halves of the seasons agree in direction: 3.8 and 1.7) and does nothing in the other two.
- **Distance of the low from the cape** does not separate them.

## OPC link (parent low, "OPC agreed", `results/opc_link_s2.txt`)

An OPC hurricane-force low fix (Atlantic archive) lies within 600 nm of the gust maximum within 12 h in **92% of south, 94% of east and 96% of west
hurricane-force proxy episodes**, against 24%, 36% and 54% for storm-force episodes. OPC archive has low fixes only, so this is through the parent low, not a check of the tip jet itself.

## What could not be tested, and what is post hoc

- The pre-registered fit-and-hold-out tests **did not run**: no setup had 30 of the rarer class per half (hurricane-force episodes: 9, 18, 26 in the original classes). `results/tests.csv` is kept as the record and **must not be quoted**.
- The tables above are **post hoc**, all 22 seasons, ridge-penalised logistic fits, season-block bootstrap, BH over 15. The south setup was **widened after the pull** (bearing 45-225 instead of 135-225) because the original left the classic forward tip jet (122 episodes) in "none of the three"; the original-class results are in `results/posthoc_orig_classes.txt` (south: 36 episodes, 9 hurricane-force; same signs for depth and gradient, Greenland high 1.2, cold air 1.2, neither significant).
- Detectable size: about 1.6 (south), 2.1 (east), 1.9 (west) times per standard deviation at 80% power, treating episodes as independent. The effective sample is seasons, and episodes in a season share a regime.
- Low-level stability: not tested (no pressure-level temperature pulled). Upper-air composites: not built.
- Even/odd season split (`results/audit.txt`): depth, gradient and east-setup Greenland-high effects keep their direction in both halves; the south-setup cold-air effect is 3.8 and 1.7; the west-setup Greenland high is 1.3 and 2.4.

Checks: see the pull request; not independently recomputed until a fresh agent reports there.
