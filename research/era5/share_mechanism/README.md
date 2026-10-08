# What turns cyclones into HF lows? (ERA5 proxy, pipeline A)

Plan, committed before any weekly deepening count or any ingredient met the pattern index: [PREREGISTRATION.md](PREREGISTRATION.md) (`cb4892f`).
Code: `part1.py` (RA-17), `part2.py` (RA-18 with RA-6), `power_m.py`, `make_report.py`. Results: [results/](results/) (`summary.txt`,
`part1_tests.csv`, `part2_tests.csv`, `part2_meta.json`, `power_m.csv`).

All outcomes are **pipeline A** (`research/era5/hf_history`, 800 km ocean gust index, HF-equivalent at 71.7 kt), a **proxy**, not the
archive. The pattern index is the leave-one-season-out index of hf-low PR 41, not refitted. Per +1 SD, 22 seasons 2004-05 to 2025-26,
Oct-Apr, effective n 22 seasons. Transitioning tropical cyclones are in. Environment: `intensity/results/env_2004.csv.gz` (ERA5, 1.5 degrees).

Reproduce (no ERA5 pull; about 1 min, 20 min and 15 min):

    python3 -I research/era5/share_mechanism/part1.py . research/era5/share_mechanism/results 2000 2000
    python3 -I research/era5/share_mechanism/part2.py . research/era5/share_mechanism/results 2000 2000
    python3 -I research/era5/share_mechanism/power_m.py . research/era5/share_mechanism/results 60 400
    python3 -I research/era5/share_mechanism/make_report.py research/era5/share_mechanism/results

## Answer

**1. Which step does the pattern change? The second one: whether a deepening cyclone goes on to reach HF. Not whether cyclones deepen.**

| per +1 SD of the pattern | Atlantic | Pacific |
|---|---|---|
| share of cyclones reaching HF (reproduces PR 64) | x1.271 (1.19-1.35) | x1.129 (1.06-1.21) |
| P(deepening) | x1.009 (0.988-1.029), well-powered null | x0.982 (0.964-1.002), well-powered null (q 0.11) |
| **P(HF given deepening)** | **x1.252 (1.189-1.331), q 0.0015** | **x1.137 (1.074-1.223), q 0.0015** |
| HF without a deepening fix (residual) | x1.006 (0.98-1.03) | x1.011 (0.995-1.028) |
| fraction of the log effect in P(HF given deepening) | **0.94 (0.83-1.06)** | **1.06 (0.87-1.31)** |

The agenda prediction (at least 70% in P(HF given deepening)) is **supported** in both basins by the rule fixed in advance
(whole interval above 0.5, point above 0.70). In the Pacific the point is above 1 because P(deepening) is, if anything, 2% lower. P(deepening)
is a well-powered null in both basins (interval inside 0.95-1.05; minimum detectable 1.03 and 1.03). The verdict is the same with
fixed-depth HF (Atlantic 0.82, 0.69-0.92), with the deepening threshold at -2.4 or -6.0 hPa per 12 h (0.95, 0.85 Atlantic), in both halves
of the record, and without Atlantic peaks above 60N (0.96). Secondary indices (S3; the verdict is conversion for the Greenland high, mixed for NAO):
the Greenland high (Atlantic, negative effect on the share, x0.872 per SD) acts through conversion (P(HF given deepening) x0.887, fraction 0.88, 0.65-1.14;
P(deepening) x0.985, q 0.28), the answer to the RA-19 relay: the high does not change how many storms deepen but lowers the chance a deepening
storm reaches HF. NAO: conversion fraction 0.77 (0.39-1.05, mixed), with a small deepening term (x1.025, q 0.044). PNA (Pacific): no net share effect (x1.006) so no fraction;
P(deepening) is x0.964 per SD (q 0.001 in its family of 3), offset by a non-significant rise in conversion.

**2. How much of that is carried by the near environment? About a third in the Atlantic, and the Pacific cannot say.**

| | Atlantic | Pacific |
|---|---|---|
| cyclones (deepening, with a reference fix and all eight ingredients) / HF | 3,285 / 754 | 4,247 / 691 |
| effect of the pattern on P(HF given deepening), log-odds per SD | 0.343 (0.26-0.43), q < 0.001 | 0.200 (0.12-0.30), q < 0.001 |
| **mediated fraction M (g-formula, all eight ingredients)** | **0.31 (0.19-0.47)** | **0.27 (0.02-0.72)** |
| log-odds difference method | 0.23 (0.08-0.42) | 0.28 (-0.07 to 0.77) |
| verdict by the rule | **little mediated** (point below 1/3, upper bound below 0.5) | **partly mediated / unresolved** |

Prediction P18 (at least half mediated) is **not supported** in the Atlantic (interval 0.19-0.47 excludes 0.5). The Atlantic verdict
sits on the boundary: the point (0.312) is just under 1/3, and adding the five secondary ingredients gives 0.34 (0.22-0.49), which
the rule calls "partly mediated"; leaving the seasonal cycle in gives 0.36 (0.24-0.52). So the defensible statement is "about a third (0.2-0.5), not
half or more". The Pacific interval is too wide to say anything between 0 and 0.7. So roughly two thirds of the pattern's effect on conversion sits in something these
1.5 degree near-storm fields (read 12 h before first deepening) do not hold: mesoscale structure, timing, moisture sources, or limits of the proxy. That fits PR 56's
small gains from further predictors.

What carries the third (Atlantic, single-ingredient M, unique part in brackets): Eady growth 0.23 (0.13), jet maximum 0.18 (0.04), surface heat flux
0.08 (0.00), 300 hPa divergence 0.06 (0.02); the others are about 0. Pacific: jet 0.28 (0.07), Eady 0.21 (0.11), flux 0.14, column water vapour -0.13 (it falls with the pattern but raises HF
within a month, so it opposes). The pattern moves jet, Eady, SST minus T500, flux and divergence (Atlantic a-path q 0.004-0.005; five of eight in each basin
pass q < 0.05) but only by 0.04 to 0.13 SD per SD of the pattern, small against how strongly each relates to HF (below).

**RA-6 folded in: within a month, stronger-than-normal jet and Eady growth do raise P(HF given deepening).** Odds ratios per SD with the seasonal
cycle removed (month fixed), reference fix 12 h before first deepening: Atlantic jet 2.05 (1.88-2.26), Eady 1.96 (1.81-2.12); Pacific jet 2.25 (2.06-2.50), Eady 1.63 (1.50-1.80). The prediction (above 1.05) holds in
both basins (q 0.001); 14 of 16 b-path tests pass, the exceptions being SST minus T500 (0.93 and 0.96, n.s.). These are markers: the environment is read near a low that is already organising, even with the 12 h lag.

## Tests and multiplicity

Primary tests: Part 1 six (A1-A3 x 2 basins; BH within the six: A2 passes in both, q 0.0015; the others q 0.11-0.59), Part 2 beta_tot x 2 (both pass), a-path 16
(10 pass), b-path 16 (14 pass), 40 in all of which **28 pass q < 0.05 within their own family**. M and its fractions are estimates with intervals, decided by the rules, not
p-tested. Bootstrap p floors are 0.001 (2,000 draws). Secondary variants are their own families (`part1_tests.csv`, `part2_tests.csv`).

## Looks, power and limits

- **Looks at the 2015-25 held-out seasons: 7 (Part 1) and 8 (Part 2)**, logged in `hemispheric/results/heldout_looks.log`. Both decompose an effect already found.
- Power: Part 1 minimum detectable rate ratios are 1.03 (P(deepening)), 1.09 (conversion Atlantic) and 1.10 (conversion Pacific). Part 2, planted-effect check of the decision rule (`results/power_m.csv`, 60 simulations per cell, 400 draws each; approximate): Atlantic, true M 0 / 0.3 / 0.5 / 0.8 returns "little" in 100% / 85% / 0% / 0% of runs and "mostly" in 0% / 0% / 30% / 82%; so a true M of 0.5 would usually read as "partly" (70%) and the observed Atlantic "little" is what a true M near 0.3 gives, not 0.5. Pacific: true M 0.3 gives "little" 32%, "partly" 68%; 0.5 gives "partly" 85%; 0.8 gives "partly" 82% ("mostly" never), so the Pacific cannot separate any of these.
- The deepening flag is a whole-track property, so P(HF given deepening) has a mechanical element (HF storms usually deepen); 143 HF cyclones lack a deepening fix (residual factor, near 1).
- Mediation is statistical attribution under no unmeasured confounding, which is not credible; the environment is a 1.5 degree reanalysis field. Pipeline A only; a change of HF definition
  from the "counts vs published" work is bounded by the fixed-depth variants (same verdicts). No 1979-2000 replication was possible (no per-fix pressure history before 2004).

## Deviations (post hoc)

1. `part1.py`: the bootstrap-stability guard for fractions assumed a positive share effect, so the Greenland-high secondary (negative effect) returned missing intervals on the first run. Made direction-neutral (sign of the observed effect) and rerun; no primary number changed.
2. `power_m.py`: the first run bisected the planted direct coefficient over a range where M is not monotonic (planted M came out negative); the range was restricted to non-negative coefficients and the simulation rerun.
3. The preregistration quotes 789 Pacific HF cyclones with a deepening fix (counted over genesis months Oct-Apr); the weekly window holds 788 (one starts on day 210, outside it). Window edge, not an error.
4. A 20-draw test run of `part2.py` was read before the full run (same seed, same estimates; intervals only changed).

## Independent check

See `VERIFICATION.md` (fresh Sonnet agent, own code, 300-draw bootstrap): every quoted point estimate, sample size and rate ratio matched. **Not independently checked:** permutation p and BH q values and the 28-of-40 count, secondary variants S1-S5 in both parts, single-ingredient and unique M_j, `power_m` results and minimum detectable effects, and the wording of verdicts against the rules.
