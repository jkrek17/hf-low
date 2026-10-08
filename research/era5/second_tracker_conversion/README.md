# Does the pattern still turn deepening storms into HF lows when a second tracker finds the storms? (RA-30, ERA5 proxy)

Plan, committed before any tracker track met a deepening flag, an HF label or the pattern index: [PREREGISTRATION.md](PREREGISTRATION.md) (commit `5d38151`; post hoc deviations at its end).
Code: `vpressure.py` (central pressure for V fixes), `build.py` (per-track tables), `run.py` (weekly models). Results: [results/](results/) (`tests.csv`, `tracks_{M,V}.csv.gz`, `build_report.csv`).

All numbers are an **ERA5 proxy**. The cyclones come from the PR 116 trackers (**M** pressure minima, **V** 850 hPa vorticity; 1.5 degrees, 6-hourly, seasons 2004-05 to 2021-22, 18 seasons; RA-15 coverage ends January 2023). The **HF label is pipeline A's** (`research/era5/hf_history`, 800 km ocean gust at 71.7 kt), carried onto a tracker track when that track is the PR 116 match of an A HF event. "Pipeline A" rows are pipeline A's own tracks on the same 18 seasons and the same weekly design (re-run, so the comparison is like for like). The pattern is the PR 41 leave-one-season-out index (not refitted). Per +1 SD, effective n = 18 seasons, transitioning tropical cyclones in, deepening = 12 h pressure fall of at least 3.6 hPa at a 00/12 UTC fix (PR 85's rule).

## Answer

**Yes, it holds in both basins, with both trackers.** The pattern's effect on P(HF | deepening) keeps its sign and size when the deepening storms come from a different tracker. Every point estimate lies inside PR 85's pipeline A interval, every test passes FDR, and the paired difference from pipeline A includes 1.0 in all four cases (best read as no difference; the Pacific M-minus-A difference only just includes it, 1.000 to 1.111).

| rate ratio of P(HF given deepening) per SD | Atlantic | Pacific |
|---|---|---|
| PR 85, pipeline A, 22 seasons (reference) | 1.252 (1.189-1.331) | 1.137 (1.074-1.223) |
| pipeline A, these 18 seasons | 1.252 (1.185-1.337) | 1.117 (1.045-1.216) |
| **tracker V (primary)** | **1.236 (1.154-1.327), q 0.001** | **1.149 (1.069-1.245), q 0.003** |
| tracker M (secondary) | 1.214 (1.129-1.310), q 0.001 | 1.172 (1.080-1.286), q 0.001 |
| V relative to pipeline A (paired, same seasons) | 0.987 (0.944-1.032) | 1.028 (0.990-1.064) |
| M relative to pipeline A (paired) | 0.970 (0.913-1.026) | 1.050 (1.000-1.111) |
| minimum detectable effect (V / M) | 1.10 / 1.11 | 1.10 / 1.11 |

The decision rule fixed in advance ("holds": rate ratio above 1, q below 0.05, point inside PR 85's interval) is met for all four of {V, M} x {Atlantic, Pacific}. Q is Benjamini-Hochberg over the registered family of four; all 4 pass (q 0.001 to 0.003); across all 20 conversion tests run on the trackers (primary plus secondaries) the largest q is 0.018 (BH across all 40 tracker tests, deepening and conversion, column `q_all` of `tests.csv`). Intervals: 95% season-block bootstrap (2,000 draws); p by season-block permutation of the index (2,000, floor 0.0005).

## Sensitivities (registered secondaries; each its own BH family; none changes the answer)

- **Deepening threshold -2.4 hPa per 12 h** (the 1.5 degree grid is shallower than pipeline A's 0.25 degrees): V 1.250 Atlantic, 1.141 Pacific; M 1.218 and 1.164; all q below 0.01.
- **Match distance halved or doubled** (HF-label noise): V 1.244 / 1.217 Atlantic, 1.152 / 1.144 Pacific; M 1.215 / 1.219 and 1.175 / 1.186. Doubling the distance lifts M's Pacific paired difference to 1.062 (1.020-1.123), the only case where a tracker's ratio is above pipeline A's with an interval excluding 1; read as a looser label adding some correlated non-HF tracks, not as a different mechanism.
- **Atlantic north of 60N** (peak at 60N or north; M recall there was lower): V 1.229 (1.047-1.454, q 0.019), M 1.189 (1.076-1.324, q 0.016), pipeline A 1.190. South of 60N: V 1.218, M 1.222, A 1.270. The effect holds in the far north; the V interval is wide (minimum detectable 1.26), so the V north-of-60N result is a consistent estimate, not a sharp one. The difference from pipeline A is not resolved either side of 60N.
- **Whether storms deepen at all (P(deepening))**, the PR 85 well-powered null: Atlantic V 1.013 (0.995-1.032), M 1.006 (0.974-1.041); Pacific V 1.001 (0.984-1.017), M 0.977 (0.960-0.993, q 0.067). All intervals except Pacific M and pipeline A's own Pacific (0.975) include 1; PR 85's Pacific value was 0.982. Atlantic south of 60N, V only: 1.031 (1.017-1.048, q 0.010), a secondary lead not seen with M or pipeline A.

## What this means

- PR 85's central statement (the pattern acts on conversion, not on deepening) does not depend on pipeline A's detector, linker or variable: the vorticity tracker and a differently built pressure tracker give the same ratios, 1.21-1.24 Atlantic and 1.15-1.17 Pacific.
- It does not make the effect independent evidence: the HF label is still pipeline A's gust label, the pattern index and the 18 seasons are the same ones, and the trackers' agreement with each other is partly the same ERA5 fields. It is a robustness check on the denominator, not a replication.
- The size is not sharpened: the intervals are as wide as PR 85's. The Pacific rate ratio sits 1.12-1.17 across pipelines against 1.14 in PR 85; the data cannot rank them.

## Limits

- 1.5 degree fields; V pressure is the minimum 1.5 degree MSLP within 500 km of the vorticity centre (a nearby deeper low can be picked up), the V previous-week count for week 0 covers only 26-30 September.
- The numerator loses A events that no tracker track matches (PR 116 recall 89-98%) and events that share one track; HD counts are 586-629 (Atlantic) and 558-609 (Pacific) against 736 and 642 for A. This bias pushes ratios toward 1, not away from it, and the paired differences show no sign of it.
- Seasons 2022-23 to 2025-26 are not covered. Fit and test use 2004-05 onward only (Decision 1). No pre-2004 season used. Sensitivity re-run on the same seasons and pattern index, not a new pattern search; it is a look at the 2015..2021 seasons, logged in `research/era5/looks/second_tracker_conversion.log` (agenda assigns its number).
- Not done: the planted-effect power simulation and S5 (see deviations).

## Data and reproduction

    python3 research/era5/second_tracker_conversion/vpressure.py /mnt/project-files/second-tracker/V_fixes.csv.gz work/v_mslp500.csv.gz   # 1.10 GB, resumable
    python3 research/era5/second_tracker_conversion/build.py /mnt/project-files/second-tracker work/v_mslp500.csv.gz research/era5/second_tracker_conversion/results
    python3 research/era5/second_tracker_conversion/run.py research/era5/second_tracker_conversion/results research/era5/second_tracker_conversion/results 2000 2000

The tracker fix files come from the PR 116 thread (`/mnt/project-files/second-tracker/`, 11.9 MB and 10.7 MB, not committed); they are regenerable with `second_tracker/detect.py` and `track.py`. `build.py` asserts 83,038 M and 50,040 V tracks on read and that its HF labels reproduce PR 116's matched track sets exactly (1,495 and 1,481 events). The per-track tables `results/tracks_{M,V}.csv.gz` (1.8 MB, 1.2 MB) are the committed denominator. Needs numpy, pandas, scipy, numcodecs.

## Verification

See [VERIFICATION.md](VERIFICATION.md).
