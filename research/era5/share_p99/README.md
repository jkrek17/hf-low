# RA-28: does scoring HF on the 99th percentile of the owned gust cells change the pattern's effect on conversion?

ERA5 **proxy**, **pipeline A** (`research/era5/hf_history`). Plan: `PREREGISTRATION.md` (committed first, `1a58bd2`; cuts committed before any outcome).
Looks entry: `../looks/share_p99.log`. A sensitivity re-run of PR 85 (still open; this PR carries its own copy of the estimation code), not a redefinition of HF.

## Answer in plain words

**No, a smoother index does not strengthen the pattern's effect on conversion, in either basin.** Per SD of the hemispheric pattern, on 00/12 UTC fixes, both labels count-matched to the archive:

| | published label (6-hourly, 71.7 kt) | maximum, 00/12 (cut 69.8 kt) | 99th percentile, 00/12 (cut 64.2 kt) | Delta (p99 - max) |
|---|---|---|---|---|
| Atlantic P(HF given deepening) rate ratio | 1.252 | 1.203 | **1.180** [1.125, 1.242] | -0.024 [-0.064, +0.015] |
| Pacific | 1.137 | 1.137 | **1.140** [1.085, 1.205] | +0.003 [-0.045, +0.055] |

- **Sign kept, both basins** (q 0.001 each). The agenda's predicted rise (toward 1.3 Atlantic, 1.2 Pacific) did not happen.
- **Pacific: keeps, well-powered.** The paired interval lies inside the pre-registered +-0.06 margin. The index choice does not move the Pacific effect.
- **Atlantic: can't tell, leaning slightly weaker.** The point estimate is 0.024 lower and the interval just misses the +-0.06 margin (lower end -0.064). Nothing points up. Minimum detectable change in the rate ratio is about 0.06 (Atlantic) and 0.07 (Pacific).
- So label noise at the cut does not look like what limits the pattern signal. That is a statement about the gust 99th percentile; RA-27's smoother-index result was for sustained wind, and nothing here shows the gust p99 is less noisy. Against the archive, p99 agrees no better than the maximum (HSS 0.73 for both, difference -0.006 [-0.025, +0.012], 5 seasons).
- **The 00/12 sampling itself matters more than the index.** In the Atlantic the 00/12 maximum label gives 1.203 against the published 6-hourly 1.252 (Delta +0.048 [+0.017, +0.082], paired): dropping the 06/18 UTC fixes costs part of the signal. This is why the comparison is paired on 00/12; a 6-hourly percentile would need an 80-100 GB pull, not done.
- Shares and mediation unchanged: the deepening factor stays null (A1 1.009 / 0.982 under every label) and the conversion share f_conv stays 0.90-1.11. Mediated fraction through the eight near-environment ingredients (Part 2): Atlantic 0.31 (published), 0.36 (max), 0.38 (p99), intervals overlapping [0.19, 0.61]; Pacific 0.27, 0.30, 0.25, unresolved as in PR 85.

## What was run

Pull 38.21 GB (7,135 times, gust + MSLP, 12,772 fixes with catalog g800 >= 60 kt for the population plus the calibration seasons). Re-derived g800 matched the catalog for 12,772 of 12,772 fixes (max 0.05 kt). Calibration (`calibration.txt`): forecast count = archive count in 2021-22..2025-26 on 00/12 fixes; cuts 69.8 (max), 64.2 (p99), 62.0 (p98) kt, all above the 60 kt pull floor. 107 of 461 tracks differ in label between max and p99 in those seasons.
Reproduction gates passed: under the published label `conv.py` and `med.py` return PR 85's Part 1 primary table and Part 2 primary table with zero difference.

Files: `extract.py`, `make_pull_list.py`, `calibrate.py`, `lab.py`, `conv.py` (Part 1, paired bootstrap), `med.py` (mediation), `make_report.py`; results in `results/` (`summary.txt`, `conv_tests.csv`, `paired_delta.csv`, `part2_tests.csv`, `cuts.json`, `calibration.txt`, `fix_p99.csv.gz`).

## Tests and multiplicity

Primary family (4 tests, BH): A2(p99) > 1 Atlantic p 0.0005 q 0.001, Pacific p 0.0005 q 0.001; Delta Atlantic p 0.23 q 0.31, Pacific p 0.90 q 0.90. 0 of 2 Deltas pass. Secondary families (p98, shares, mediation) are in the CSVs; p98 behaves like p99 (Atlantic 1.180, Pacific 1.151, Delta -0.024 and +0.013, neither resolved). Global BH over all 132 pre-registered p values: 88 pass, almost all of them the unchanged a/b path tests of PR 85.

## Limits

- 00/12 UTC fixes only; tracks seen only at 06/18 UTC count as not HF under the 00/12 labels. 22 seasons shared with PR 85, so this does not replicate it.
- Looks: one more scoring of the 2015-25 seasons (sensitivity re-run; number assigned by the agenda thread). No pre-2001 season.
- Seasons 2004-05 on (Decision 1). Pipeline A, ERA5 proxy; transitioning tropical cyclones in; Atlantic north of 60N not separated.

## Deviations (post hoc)

None.

## Verification

A fresh Sonnet agent (no access to my code) recomputed from the committed tables (`verify/verify_part1.py`, `verify_part1.txt`): population sizes, HF/D/HD track counts per basin and label, the 8 A2 rate ratios (all within 0.002 under the lag convention that includes tracks from the 7 days before 1 Oct, as `chanlib.prev_week_counts` does; the Atlantic max-label ratio is 1.2056 against 1.203 under the stricter population-only lag, off by 0.0026), and the max/p99 label disagreements among deepening tracks (Atlantic 218, Pacific 235). Under the published label, my code also reproduces PR 85 exactly (zero difference, Part 1 and Part 2).
**Not independently checked:** bootstrap intervals, permutation p and q values, the paired Delta intervals and MDEs, the Part 2 mediation numbers (reproduced against PR 85's table only under the old label), the calibration cuts and archive HSS, and the extraction itself beyond the g800 reproduction gate.
