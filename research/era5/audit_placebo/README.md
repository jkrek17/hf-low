# Audit pack: placebo pattern indices, conversion versus depth, test register (ERA5 PROXY, pipeline A)

Plan, committed before any placebo index was fitted: [PREREGISTRATION.md](PREREGISTRATION.md) (`23f1b5d`). Origin: the panel review of 2026-10-11 (critiques
C2, C3). No ERA5 pull (0 GB): committed files only. All outcomes are **pipeline A** (800 km ocean gust index, 71.7 kt), a **proxy**. Per +1 SD, 22
seasons 2004-05 to 2025-26, Oct-Apr. Transitioning tropical cyclones are in. Look at 2015-26 logged in `../looks/audit_placebo.log`.

Code: `placelib.py` (index refit), `run_a.py` + `summarise_a.py` (Part A), `run_b.py` (Part B), `run_c.py` (Part C). Results: [results/](results/).
Reproduce (about 30 min, 4 cores; Part B about 25 min):

    python3 -I research/era5/audit_placebo/run_a.py research/era5/audit_placebo/results 2000 2000 100
    python3 -I research/era5/audit_placebo/summarise_a.py research/era5/audit_placebo/results
    python3 -I research/era5/audit_placebo/run_b.py research/era5/audit_placebo/results 2000
    python3 -I research/era5/audit_placebo/run_c.py research/era5/audit_placebo/results

## Answer

**1. The "pattern acts through conversion" localisation is not specific to the HF-trained pattern index.** Indices refit on other cyclone outcomes also raise
the share of cyclones that reach HF, and the raise sits in the conversion step (P(HF given deepening)), not in whether cyclones deepen. The decomposition
is therefore partly a property of the denominator and the label, not evidence that the HF-trained pattern picks out a conversion mechanism. What does
stand out is size in the Atlantic.

| index (trained on) | basin | share of cyclones reaching HF (A4) | P(HF given deepening) (A2) | f_conv | reading (fixed in advance) |
|---|---|---|---|---|---|
| P0 archive HF counts (refit of PR 41) | Atl | x1.288 (1.219-1.361) | x1.273 (1.202-1.344) | 0.95 (0.84-1.09) | conversion |
| P1 all cyclones | Atl | x0.954 (0.895-1.009) | x0.992 | not defined | no resolved share effect |
| P2 deepening cyclones | Atl | x1.059 (1.003-1.125) | x1.092 (1.048-1.141), q 0.028 | 1.54 (1.06-3.55) | conversion |
| P3 strong non-HF (gust 55-71.7 kt) | Atl | x1.062 (1.028-1.101) | x1.060 (1.023-1.103), q 0.155 | 0.96 (0.68-1.18) | conversion |
| P0 | Pac | x1.132 (1.063-1.202) | x1.143 (1.070-1.219) | 1.08 (0.89-1.31) | conversion |
| P1 all cyclones | Pac | x1.080 (1.014-1.158) | x1.120 (1.053-1.198), q 0.004 | 1.48 (1.13-3.06) | conversion |
| P2 deepening cyclones | Pac | degenerate (no PCs chosen in any fold) | | | no placebo possible |
| P3 strong non-HF | Pac | x1.104 (1.031-1.165) | x1.097 (1.021-1.163), q 0.022 | 0.94 (0.58-1.28) | conversion |

Four of the five definable placebos have a resolved share effect with a conversion fraction whose interval stays above 0.5. The Atlantic HF-trained index
has the largest share effect (x1.29 against x1.06 for the placebos); in the Pacific the HF-trained index (x1.13) is no larger than the placebos (x1.08, x1.10).

**Shuffle null (100 season-shuffled P0 indices per basin).** About half of the shuffled indices are constant (cross-validation chose no PCs in any fold: 47 of 100 Atlantic, 53 of 100 Pacific), so S is defined for 53 and 47 draws. The share of the absolute effect that sits in conversion, S, is 0.93 for the Atlantic P0 index against a median of 0.65 and 95th percentile of 0.91 among the valid shuffles (1 of 53 at or above; the verifier's own draws gave 3%); 0.78 for the Pacific against 0.61 and 0.87 (11 of 47 at or above; the verifier 7%). An earlier version of this file gave 1% and 11% of 100 draws, dividing by all draws including the undefined ones; that was wrong. The Atlantic value is above most chance indices and the Pacific value is not distinguishable from them, but the comparison is weak and the exceedance rates move from draw to draw.

**Interval caveat.** f_conv is a ratio whose denominator A4 is near zero for the placebos, so its bootstrap interval is heavy-tailed and its upper bound is not reproducible (Atl P2 upper 3.55 here, 5.17 for the verifier; Pac P1 3.06 against 3.47 to 3.35). The lower bound for Atl P2 depends on dropping draws with A4 <= 0 (the plan's `frac_boot_ok` rule): 1.06 here, 0.83 to 1.00 on the verifier's draws. The reading "conversion dominates" (lower bound above 0.5) holds on all of these. Atl P1 "no resolved share effect" holds only marginally (A4 upper limit 1.009 here, 1.010 for the verifier). The A4 and A2 rate-ratio intervals reproduce within bootstrap noise.

**2. Conversion is partly "reached deeper" and partly gust at a given depth.** Among deepening tracks (committed `oos_index`, track-level logistic), adding the track's
minimum pressure cuts the pattern's log-odds on HF by 0.31 in the Atlantic (0.01 to 0.60) and 0.51 in the Pacific (0.10 to 1.25). Neither meets the rule for "mostly depth"
(lower bound above 0.25 with point at least 0.5) or "mostly gust at a given depth", so both are **mixed**. Adding a 5 by 10 degree position cell does not reduce the effect further
(Atlantic 0.11, Pacific 0.29 attenuation, intervals include 0): position is not what carries it. minp is itself moved by the pattern, so this is a decomposition, not an adjustment for confounding.
Atlantic M2 q 0.080, Pacific M2 q 0.036 (BH over the four attenuation intervals).

**3. Test register (machine-readable only).** 114 result files on the research branch carry a p or q column, across 42 result directories, with 4,225 rows; 1,258 are marked q < 0.05
by their own file. Many rows are sensitivity variants and sub-family rows of the same tests, so this is rows, not independent tests, and it excludes anything reported only in prose or JSON. It
supports the review's point that multiplicity is controlled within families and never across them; it does not estimate a global false discovery rate. (Nominal passes expected
if every row were a global null: about 211.)

## Deviations from the plan, and what did not reproduce
- **Reproduction gate not met.** P0 correlates 0.938 (Atlantic) and 0.852 (Pacific) with the committed `oos_index`; the plan required 0.90 in both. By the rule fixed in advance, Part A compares
  placebos against the P0 refit only. P0's own decomposition is close to PR 85 (A2 x1.273 / x1.143 against x1.252 / x1.137), but the index is an approximation (no SST, whole-record EOFs).
- **P2 Pacific is degenerate:** cross-validation chose no PCs in every fold, so the index is constant. The initial output gave it p = 0.0005 (a NaN artefact). It is excluded, and BH is over the 15
  defined placebo tests, not the 18 in the plan (deviation, made after seeing the artefact).
- Part B's logistic M1 slope (0.287 Atlantic, 0.169 Pacific per SD) differs from PR 85's weekly Poisson x1.252 (0.225 log) and x1.137 (0.128), because the unit and model differ (tracks versus weeks, no previous-week term). It is the reference for the attenuation only.

## Limits
- The index is not refit inside the bootstrap, so intervals understate the index's own noise (stated in the plan).
- "Strong non-HF" (gust 55-71.7 kt) is this audit's label, not the project's storm-force class.
- Placebo share effects are small; with 22 seasons, the conversion fraction intervals are wide (P2 Atlantic 1.06-3.55). The result is "conversion is not specific", not "a mechanism is absent".
- Not in this audit: noisy-ingredient mediation check; archive-label conversion; a nested bootstrap. Permutation p-values and the BH correction were not independently recomputed; point estimates and bootstrap intervals for A4, A2 and f_conv were.
- Independent recompute: see VERIFICATION.md (two fresh agents).
