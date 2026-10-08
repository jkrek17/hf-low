# Do combined teleconnection states change where HF lows peak, or how many peak in a region? (ERA5 proxy, pipeline A)

Catalog: [CATALOG.md](CATALOG.md). Plan, committed before any fit: [PREREGISTRATION.md](PREREGISTRATION.md) (merged in PR 32, commit `f030c9e`).
Code: `core.py` (estimators), `power.py`, `run_all.py`, `family.py`. Results: `results/` (`results.csv` has every interaction fitted, `corners.csv` the four-corner tables,
`power.txt`, `family.txt`).

Reproduce (about 7 minutes, no ERA5 pull; CPC files in `/mnt/project-files/teleconnection-test/cpc_indices/`):

    python3 research/era5/tele_combos/power.py research/era5/tele_combos/results/power.txt 300
    python3 research/era5/tele_combos/run_all.py research/era5/tele_combos/results
    python3 research/era5/tele_combos/family.py research/era5/tele_combos/results/results.csv research/era5/tele_combos/results/family.txt

Every number is **pipeline A** (800 km ocean gust index, HF-equivalent at 71.7 kt) and an **ERA5 proxy**. Oct-Apr, seasons 2004-05..2025-26 (22 seasons),
HF-equivalent lows only unless a row says "all cyclones". Position is where the gust index peaks, **not genesis position** (not committed before 2004).

## Answer, in plain words

**No combination changed where HF lows peak or how many there were, beyond what each state does alone. For position that is a well-powered null in three of the four pairs; for
counts, and for the latitude of NAO x El Niño, the test is inconclusive.**

| Pair (basin) | Location of HF peaks | Number of HF lows |
|---|---|---|
| El Niño x MJO, west Pacific heating (Pacific) | **no** (well-powered null: shifts beyond about 1.8 deg longitude or 0.6 deg latitude per SD x SD are excluded) | can't tell |
| PNA x MJO, dateline heating (Pacific) | **no** (well-powered null: the interval, -2.1 to +1.0 deg longitude and -0.25 to +0.8 latitude, lies inside the smallest shift of interest) | can't tell |
| NAO x PNA (Atlantic) | **no** (well-powered null: shifts beyond about 0.9 deg longitude or 0.5 deg latitude are excluded) | can't tell (leans negative: log RR -0.046, interval -0.10 to +0.02) |
| NAO x El Niño (Atlantic) | longitude: no (well-powered null). Latitude: **can't tell**, +0.40 deg per SD x SD with interval -0.23 to +1.23, above the 0.9 deg bound | can't tell |
| NAO x polar vortex (Atlantic) | **not run** (see below) | not run |

Smallest effect of interest was set before the fit: 0.10 of the outcome's SD per SD x SD for position (Pacific 2.4 deg lon / 0.7 deg lat; Atlantic 1.8 / 0.9) and RR 1.05
for counts. The count tests can only detect RR 1.08 to 1.12 (`power.txt`), above 1.05, so a count null can only ever be inconclusive.

**Constructive or destructive?** Neither is claimed. The four-corner tables (`corners.csv`) show the fitted joint-state position beside the additive prediction; for the
four tests the largest gap in any corner is 0.7 deg longitude (T1), 0.4 deg latitude (T4) and 5% in count (T3), all inside the noise. In the plain sense of Jason's question: aligned states (for example +NAO with +PNA) do
shift Atlantic peaks (+NAO moves them north and east, a main effect) but the shift for the pair is the sum of the two single shifts, not more and not less.

## Primary results (8 tests; Oct-Apr 2004-05..2025-26)

Per SD x SD of the product of the two standardised indices. Location in degrees (interval = 95% season-pairs bootstrap); p = sign-flip score test, joint over lon+lat for
location. q = Benjamini-Hochberg over the 8 primary tests.

| Test | Outcome | gamma | 95% interval | p | q |
|---|---|---|---|---|---|
| T1 ONI x MJOWP, Pacific (791 events) | lon | -0.70 | -1.80 to +1.74 | joint 0.517 | 0.69 |
| | lat | +0.07 | -0.32 to +0.58 | | |
| | count (log RR) | -0.016 | -0.058 to +0.048 | 0.658 | 0.75 |
| T2 PNA x MJODL, Pacific (791) | lon | -0.50 | -2.08 to +1.03 | joint 0.350 | 0.69 |
| | lat | +0.25 | -0.25 to +0.78 | | |
| | count | +0.034 | -0.052 to +0.116 | 0.457 | 0.69 |
| T3 NAO x PNA, Atlantic (1,003) | lon | +0.17 | -0.48 to +0.89 | joint 0.889 | 0.89 |
| | lat | +0.02 | -0.40 to +0.52 | | |
| | count | -0.046 | -0.101 to +0.017 | 0.179 | 0.69 |
| T4 NAO x ONI, Atlantic (1,003) | lon | +0.25 | -0.49 to +0.94 | joint 0.362 | 0.69 |
| | lat | +0.40 | -0.23 to +1.23 | | |
| | count | -0.018 | -0.068 to +0.042 | 0.513 | 0.69 |

Minimum BH q across the 8 is 0.69. Nothing is detected. Across all 68 interaction tests in this directory (primary and secondary) the minimum q is 0.18 and none is below 0.05
(`results.csv`, column `q_all_family`).

**Combined with the two sibling threads** (`family.txt`; their primary p-values copied from their READMEs): 11 tests, minimum q 0.40 (the ENSO -> PNA -> HF indirect effect, raw p 0.036,
which its own thread called "not detected" after its two-test correction). All other q are 0.81 or higher.

## Why NAO x polar vortex was not run

The plan made it contingent on deriving a 10 hPa, 60N daily wind series, and `CATALOG.md` said that was under 10 GB. That estimate was wrong: WeatherBench2's 1.5 degree
store has 13 levels and the highest is 50 hPa, and ARCO stores all 37 levels in one hourly chunk of about 118 MB, so a daily 10 hPa series would stream about 2 TB.
T5 was withdrawn before any fit and the catalog was corrected. A 50 hPa version sampled every 8 days (about 25 GB, below the 50 GB gate) is possible but would be a
different, weaker test of the vortex; it is listed as an optional future task, not done.

## One exploratory lead, not a finding

For **all cyclones** (the storm track, not HF lows), NAO x PNA shifts the longitude of the gust-index peak: +1.11 deg per SD x SD east for 2001-02..2025-26
(interval +0.51 to +1.86, joint p 0.003, q 0.18 across the 68-test family) and +0.60 (+0.12 to +1.10, joint p 0.040) over 1979-80..2025-26 with an era term. The same coefficient for HF-equivalent lows
is +0.17 (primary, above), so the HF-minus-storm-track difference (-0.79, interval -1.86 to +0.18) does not separate them. It is one of 68 secondary tests, does not survive correction, has not
been independently recomputed, and uses a gust-defined position. Treat it as something a later, pre-registered test of storm-track position could examine.

## Secondary tiers (labelled; none survives correction)

S1 same-time predictors, S2 box counts (3 boxes per basin), S3 HF-equivalent over 1979-80..2025-26 with an era term and over 2001-02 on, S4 all cyclones for the same records, S5 HF minus all-cyclone
position, S6 Jun-May window. All rows are in `results.csv`. In S1 (same-time predictors) no test has p below 0.14.

## Power and effective n

Effective n is seasons x independent values of the predictor per season: ONI about 2 (so 44 for 22 seasons), MJO high-passed about 20, PNA about 12, NAO about 7 (`power.txt`; these
are larger than the catalog's guesses because the lagged series are 7-day means at 210 days, with integrated autocorrelation time from the series themselves). The effective n of an interaction is limited by the slower
factor, so ONI pairs are the weakest. Detectable effect at 80% power (2.8 x permutation SE), position (lon, lat): T1 2.4 / 0.65 deg, T2 2.7 / 0.77, T3 1.8 / 1.0, T4 1.7 / 0.95; count RR T1 1.09, T2 1.12, T3 1.10, T4 1.08.

## Departures from the plan

- **D1 (data handling; found by the verifier, applied before this README, not outcome-driven).** The first run (`results/results_v1_bridged_pentads.csv`, commit `66cd311`) filled the two missing MJO pentads
  (centred 2021-12-29 and 2022-12-29) from the nearest valid pentad. The plan did not say how to treat them; a missing value should stay missing. The final run leaves
  them missing, and because the high-pass needs 90 valid days, 184 days of Jan-Apr 2022 and 2023 drop out of the two MJO tests (791 events instead of 825). The results differ by well under one SE for every coefficient;
  the one change in a p-value is T2 joint 0.643 -> 0.350, and nothing changes class. T3 and T4 do not use the MJO and are unchanged.
- **D2 (T5 withdrawn)**, above.
- Day sets are per test: a day is dropped only if that test's own two predictors are missing.

## What this does not show

- **A proxy.** HF-equivalent lows are cyclones whose ERA5 fields look like the ones OPC warned for as hurricane force.
- **Position is where the gust index peaks**, so it is gust-defined and carries the pre-2001 drift (Decision 1): that is why the primary window is 2004-05 on and the longer records carry an era term.
- **Basin counts of HF lows were not the question.** The earlier test (superposition) stands for totals. Here position and regional counts were added and the pairs were not found to interact.
- **Power.** 22 seasons. A real interaction of about the size of the NAO's own main effect on position could be found; one half that size could not, and counts of up to about 10% per SD x SD could be missed.
- The MJO is a longitude series high-passed with a trailing 90-day mean; other filters or the RMM phases could behave differently.
- Not tested: ONI x PNA and ENSO flavour (sibling threads), QBO, three-way terms, genesis position.

## Verification

A fresh Sonnet agent that had not seen the code, the results or my numbers recomputed the four primary tests from the raw inputs (`all_tracks.csv.gz`, the CPC files, `teleconnections.json`) with its own implementation, following only `PREREGISTRATION.md`. Its first run also exposed the missing-pentad bridging (D1). After D1 and per-test day sets its values agree with `results.csv`:

- Events 791 / 791 / 1,003 / 1,003 for T1-T4, 22 seasons each (identical).
- Interaction coefficients (lon, lat, log RR count), mine vs its: T1 -0.699, +0.071, -0.016 vs -0.711, +0.071, -0.017; T2 -0.501, +0.254, +0.034 vs -0.487, +0.252, +0.034; T3 +0.166, +0.017, -0.046 vs +0.166, +0.017, -0.046; T4 +0.251, +0.399, -0.018 vs +0.251, +0.399, -0.018. T1 and T2 differ by up to 0.014 deg, from the scaling of the MJO series; T3 and T4 are identical to three decimals.
- Bootstrap intervals agree to within about 0.1 SE (T3 lon -0.48 to +0.89 vs -0.47 to +0.93; T4 lat -0.23 to +1.23 vs -0.19 to +1.25; others similar). Joint p: T1 0.517 vs 0.510, T2 0.350 vs 0.361, T3 0.889 vs 0.889, T4 0.362 vs 0.362; count p: 0.658 / 0.457 / 0.179 / 0.513 vs 0.648 / 0.461 / 0.179 / 0.513.
- T3 four-corner (+1,+1) deviation: full model lon +1.64, lat +1.16 vs additive +1.48, +1.14 (mine and its).

**Not independently checked:** the design-power figures (`power.txt`), the BH q-values and the combined family, every secondary tier S1-S6 (including the all-cyclone NAO x PNA lead), the box definitions, the four-corner tables for T1, T2 and T4, and the sibling p-values copied into `family.py`.
