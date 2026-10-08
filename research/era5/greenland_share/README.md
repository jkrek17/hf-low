# Why does a strong Greenland high lower the Atlantic HF share? (EXPLORATORY; ERA5 proxy, pipeline A)

**Exploratory.** All 22 gust-era seasons (2004-05 to 2025-26) were already used by hf-low PR 58 and PR 63, and PR 63
reported the effect explained here. No held-out data exist, so nothing below is confirmatory. The comparisons were
written down before any outcome was related to the high: [PREREGISTRATION.md](PREREGISTRATION.md) (`7a74004`); all 21
are reported. Everything is **pipeline A** (`research/era5/hf_history`; 800 km ocean gust index, HF-equivalent at 71.7
kt), a **proxy**, Atlantic Oct-Apr, transitioning tropical cyclones in, 9,636 tracks of which 1,003 are HF. The
Greenland high (GH) and NAO are PR 63's lagged series (days -10..-4 before genesis); effects are per +1 SD.

Code: [sg.py](sg.py) (data), [run.py](run.py) (tests, 2,000 season-block draws, seed 20261008), [composite.py](composite.py),
[posthoc.py](posthoc.py). Results: [results/](results/). Reproduce (about 10 minutes, no ERA5 pull):

    python3 research/era5/greenland_share/run.py research/era5/hf_history/results/all_tracks.csv.gz <cpc_indices dir> . \
        research/era5/nao_share_barrier/gh_daily.csv research/era5/intensity/results/fixes_2004.csv.gz \
        research/era5/intensity/results/env_2004.csv.gz research/era5/greenland_share/results 2000

## Answer in plain words

**The data lean toward "the high weakens the storms", not "the high steers them away", and they cannot rule out that
the high is partly standing in for the NAO. Which of the two is right is not settled by 22 seasons.**

1. **How far the two indices can be separated: only partly.** The correlation is -0.69 over days and over storms and
   **-0.89 over the 22 season means**, which is the sample that carries the independent information. With both in
   the model the GH share ratio is 0.893 (0.818-0.981) and the NAO share ratio 1.039 (0.934-1.148); alone they are
   0.870 and 1.125 (PR 63 reproduced exactly). Variance inflation is 1.92, and the smallest GH effect detectable
   given NAO is 13% in the share (0.128 log), against 8.4% alone. The joint 95% region of the two coefficients
   contains "GH only" comfortably (Mahalanobis d^2 0.54) and "NAO only" at its edge (5.97 against a limit of 5.99).
2. **The share falls because fewer cyclones reach HF, not because there are more cyclones.** Per SD of GH, HF
   cyclones x0.895 (0.848-0.950), all cyclones x1.029 (1.010-1.051); with NAO in the model x0.910 (0.835-1.005) and
   x1.020 (0.999-1.043).
3. **Steering: real but not the explanation.** With NAO held fixed, storms start and mature farther east and leave
   the Irminger Sea and Denmark Strait box (55-67N, 50-15W): first position +1.96 deg longitude (0.94-2.95), the
   chance the minimum-pressure fix lies in the box falls 12% in odds per SD (odds ratio 0.876, 0.816-0.952), peak-gust position
   0.56 deg south. But taking position out of the share model changes nothing (the GH coefficient keeps 116%
   of its size, 97% to 149%). **STEER is not favoured.**
4. **Weakening: favoured by the pre-registered rule, with wide intervals.** With NAO fixed, a strong high goes with a
   shallower track minimum, +1.0 hPa per SD (0.25-1.64, q 0.037), and weaker Eady growth at first position
   (-0.016, q 0.030; the sample mean is 0.83); the maximum deepening rate (-0.013, q 0.14), peak gust (-0.55 kt, q 0.11) and
   250 hPa jet (-0.6 m/s, q 0.27) point the same way without passing. With minimum pressure and deepening rate in the
   model the GH coefficient retains 30% of its size; with all mediators 5%. **But** the interval for the retained
   fraction runs from -182% to +140% (C2) and -283% to +99% (C4), no C test passes q < 0.05, and minimum pressure
   sits close to the outcome (gust is set mainly by pressure gradient), so this "explains" the effect partly by
   construction. Read it as: the high's association with HF share is compatible with shallower storms, not shown to be
   caused by them.
5. **Stand-in for the NAO: not supported by the pre-registered rule, not excluded.** GH given NAO keeps 82% of the
   GH-alone coefficient (interval 20% to 134%) and its own interval excludes 0 (q 0.046); NAO given GH keeps 33% of its
   own size (-125% to +90%). Within every tercile of NAO the HF share still falls with GH (below), which is
   what a stand-in would not do. Against that, the two indices move together across seasons (r -0.89), so a
   season-level NAO effect cannot be separated from a GH effect.
6. **Reconciling PR 58 (barrier gusts up 4x) with PR 63 (share down).** Not shown here. At fixed minimum pressure and
   deepening rate, the GH coefficient on cyclone HF is -0.04 (-0.20 to +0.10) per SD: consistent with zero, with no
   positive sign. PR 58's outcome (a 0.25 degree coastal-band gust, among times with a deep low already near
   Greenland) and this one (peak 800 km gust index of a whole cyclone) differ, so the two are not contradictory, but this
   study cannot say why the same high both raises the barrier gust and lowers the share.

## What the composites show (descriptive)

[results/composites.txt](results/composites.txt), [results/composite_density.png](results/composite_density.png). Terciles
over storm-level indices (cut points: GH z -0.42, 0.41; NAO z -0.40, 0.47).

| | HF share % | n storms | first-position lon | min pressure hPa | peak gust kt |
|---|---|---|---|---|---|
| GH tercile 1 (weak) | 13.4 | 3,210 | -50.0 | 985.0 | 53.0 |
| GH tercile 2 | 10.6 | 3,205 | -48.8 | 988.0 | 50.6 |
| GH tercile 3 (strong) | 7.3 | 3,207 | -45.2 | 989.3 | 48.9 |

HF share by GH tercile within NAO tercile (cells with under 200 storms in brackets): NAO low 18.3 [153], 10.4, 7.2; NAO
middle 12.6, 10.4, 8.0; NAO high 13.3, 11.0, 6.2 [177]. At fixed GH, NAO has no consistent effect on the share
(GH tercile 3: 7.2, 8.0, 6.2). The density map shows fewer low fixes in and around the Irminger box and more over the
east Atlantic at strong GH, in the all-storm and the middle- and low-NAO panels (the low-NAO panel compares 2,253 strong-high storms with only 153 weak-high ones). The strata are not independent of the season
structure: 2,253 of the 3,207 strong-high storms are also in the low-NAO tercile.

## Numbers

Share coefficients (log, per SD; ratio in brackets), season-block bootstrap 95% interval. [results/tests.csv](results/tests.csv) has p, q and leave-one-season-out ranges for all 21 tests; headline numbers in [results/headline.json](results/headline.json).

| Test | Estimate (95%) | p | q family | q all 21 | Min. detectable |
|---|---|---|---|---|---|
| A2a share, GH alone | -0.139 (-0.195, -0.079) [0.870] | 0.0005 | 0.002 | 0.005 | 0.084 |
| A2b share, GH given NAO | -0.114 (-0.201, -0.019) [0.893] | 0.020 | 0.027 | 0.046 | 0.128 |
| A2c share, NAO alone | 0.118 (0.046, 0.178) [1.125] | 0.002 | 0.004 | 0.007 | 0.095 |
| A2d share, NAO given GH | 0.038 (-0.068, 0.138) [1.039] | 0.45 | 0.45 | 0.47 | 0.147 |
| B1 first-position lat (deg N) | -0.32 (-0.71, 0.11) | 0.17 | 0.19 | 0.20 | 0.59 |
| B2 first-position lon (deg E) | +1.96 (0.94, 2.95) | 0.0005 | 0.006 | 0.005 | 1.41 |
| B3 min-pressure-fix lat | -0.45 (-0.87, 0.01) | 0.062 | 0.089 | 0.108 | 0.65 |
| B4 min-pressure-fix lon | +1.42 (0.58, 2.29) | 0.001 | 0.006 | 0.007 | 1.19 |
| B5 peak-gust lat | -0.56 (-1.00, -0.09) | 0.022 | 0.038 | 0.046 | 0.65 |
| B6 peak-gust lon | +1.29 (0.44, 2.15) | 0.002 | 0.006 | 0.007 | 1.17 |
| B7 track min pressure (hPa) | +1.00 (0.25, 1.64) | 0.014 | 0.028 | 0.037 | 1.00 |
| B8 peak 800 km gust (kt) | -0.55 (-1.05, 0.03) | 0.067 | 0.089 | 0.108 | 0.78 |
| B9 max 24 h deepening rate | -0.013 (-0.027, 0.002) | 0.10 | 0.12 | 0.14 | 0.021 |
| B10 min-pressure fix in box (log odds) | -0.132 (-0.203, -0.050) | 0.002 | 0.006 | 0.007 | 0.111 |
| B11 250 hPa jet at first position (m/s) | -0.58 (-1.57, 0.42) | 0.24 | 0.24 | 0.27 | 1.44 |
| B12 Eady growth at first position | -0.016 (-0.027, -0.005) | 0.010 | 0.024 | 0.030 | 0.016 |
| C1 GH coef change, + first position | -0.022 (-0.043, 0.002) | 0.073 | 0.18 | 0.11 | 0.032 |
| C2 + minimum pressure, deepening | +0.095 (-0.031, 0.212) | 0.16 | 0.20 | 0.20 | 0.18 |
| C3 + first-position environment | +0.035 (-0.012, 0.077) | 0.14 | 0.20 | 0.18 | 0.064 |
| C4 + all mediators | +0.129 (0.010, 0.237) | 0.030 | 0.15 | 0.057 | 0.17 |
| C5 GH coef with C2 mediators | -0.040 (-0.201, 0.103) | 0.60 | 0.60 | 0.60 | 0.22 |

**10 of 21 pass q < 0.05 across all 21** (A2a-A2c, B2, B4, B5, B6, B7, B10, B12; none of the C tests). The storm-level
unadjusted GH coefficient in the C models is -0.135 (-0.240, -0.026) log odds per SD, on 9,622 storms with fix rows.
Fraction of that retained: C1 116% (97-149), C2 30% (-182 to 140; leave-one-season-out 2-53), C3 74% (-59 to 111;
59-82), C4 5% (-283 to 99; -26 to 25). Bootstrap one-sided probability that the C2 coefficient is at or above 0: 30%.
Minimum detectable effects are 2.8 x bootstrap se. A "null" here is inconclusive unless its minimum detectable
effect is small against the size of interest.

## Verdict against the rules fixed in advance

| Explanation | Rule | Result |
|---|---|---|
| STEER | C1 retains under 50% and a position test passes | **No.** C1 retains 116%; B2, B4, B5, B6, B10 pass but do not mediate |
| WEAKEN | C2 or C3 retains under 50%, one of B7-B9, B11, B12 passes, C5 holds | **Met on point estimates** (C2 30%; B7, B12 pass; C5 holds); intervals are very wide and no C test passes |
| STAND-IN | GH given NAO keeps under 50% of GH alone with interval including 0, NAO given GH keeps 75% of NAO alone | **Not met** (82%, interval excludes 0; NAO keeps 33%), but GH-vs-NAO "NAO only" lies at the edge of the joint region |

More than one can hold; here WEAKEN (leaning) and an unexcluded NAO contribution are not mutually exclusive.

## Post hoc (not in the plan, labelled)

[results/posthoc.txt](results/posthoc.txt). Among the 1,003 Atlantic HF events (conditional on HF), a strong high
raises track minimum pressure by +1.9 hPa per SD (0.5-3.1). Using the tele_intensity background table, roughly 0.5 hPa of
that is a higher background anomaly (-0.6 to +1.4), 0.7 is depth against the surroundings (-0.1 to +1.5) and the
rest is climatology at the new, more eastern position; none is significant alone. So even the shallower minimum is not
clearly a weaker storm rather than higher background pressure, which is why the WEAKEN reading stays a lean.

## What this cannot show

- Causation. The high is partly storm-made (the 4-10 day lag lowers this, not removes it); mediators are measured
  after genesis, so "mediation" is association.
- Separation of GH from NAO beyond what is stated: r -0.89 across season means. 22 seasons is the sample.
- Anything about proxy error near Greenland (no terrain mask; Atlantic fixes north of 60N are lower-confidence).
- Seasons outside 2004-05 to 2025-26 (Decision 1).

## Deviations and clarifications

1. **Debug run seen.** A 20-draw debug run printed point estimates (and intervals meaningless at that size) before the
   2,000-draw run. The only change made afterwards was the one in item 2, a crash/missing-data fix; no test, mediator
   set or rule was altered in response to an estimate.
2. **Missing environment values (before the full run).** `sstgrad` is missing at the first position for 27% of storms
   (2,591 of 9,636); the first debug run crashed on it. C3 and C4 fill it with the sample mean and add a missing flag.
   The 14 tracks with no fix rows are dropped from storm-level models (9,622 of 9,636). Missing deepening rate (2,233
   storms) is filled with 0 and flagged in C2 and C4. B9 uses only storms with a value.
3. **p-values** for all 21 tests are two-sided season-block bootstrap values (the plan did not name the method); PR 63
   used permutation for single coefficients.
4. **Post hoc** section above is not in the plan.

## Verification

PR 63's GH, NAO and NAO-given-GH share ratios (0.870, 1.125, 1.039) and the 9,636 / 1,003 counts are reproduced
exactly. A fresh Sonnet agent recomputed the counts, correlations, share and count ratios, tercile table, three storm-level
regressions, the logistic mediation coefficients and a 500-draw bootstrap: all matched. The box-logistic claim (B10) was
checked only by a second method of mine (the agent's script used the indicator as a covariate). Not independently
checked: 2,000-draw intervals, p and q values, minimum detectable effects, leave-one-season-out ranges, retained-fraction
intervals, the post hoc table. See [VERIFICATION.md](VERIFICATION.md).
