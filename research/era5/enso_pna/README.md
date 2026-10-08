# Does El Niño change what +PNA does to Pacific HF lows? (ERA5 proxy, pipeline A)

Plan, committed before any fit: [PREREGISTRATION.md](PREREGISTRATION.md) (commit `c6a085d`).
Code: [enso_pna.py](enso_pna.py) (reuses `../freq_split/split.py`). Results: [results/](results/)
(`summary.txt` readable, `results.csv` every number, `extra.json` asymmetry and phase table).

Reproduce (about 25 minutes, no ERA5 pull):

    python3 research/era5/enso_pna/enso_pna.py research/era5/hf_history/results/all_tracks.csv.gz <cpc_indices dir> . research/era5/enso_pna/results 2000

All numbers are **pipeline A** (800 km ocean gust index, HF-equivalent at 71.7 kt) and a **proxy**. Pacific basin,
Oct-Apr 2004-05..2025-26 (22 seasons, 10,040 cyclones, 825 HF), PNA lagged to days -10..-4 before genesis, ONI of
the month of day -7, both standardised. Effects are per SD (interaction: per SD x SD).

## Answer

1. **Interaction (pre-registered P1): inconclusive, centred on none.** ONI x PNA on HF-low frequency is RR 1.014
   (95% season-block bootstrap 0.952-1.089; permutation p 0.71, BH q 0.95). El Niño does not visibly amplify or damp
   PNA's effect, but the data cannot exclude an amplification up to 1.09 or a damping down to 0.95 per SD x SD
   (detectable at 80% power: 1.096; the pre-registered smallest effect of interest was 1.05). The interaction is also
   about 1.00 in every other cut (same-time index, Jun-May, 2001+, a 47-season depth version with era term, 1979-2000
   alone; HF-count q >= 0.78 in every cell) and in an El Niño-only / La Niña-only split (both 0.995).
2. **Pathway (pre-registered P2): real in sign, negligible in size.** El Niño does shift PNA positive (path a = +0.15 SD
   per SD ONI, bootstrap +0.001 to +0.235) and PNA drives HF lows (b = RR 1.071 per SD). The indirect effect
   is RR 1.010 per SD ONI (0.999-1.022; permutation p 0.036, x2 for the two primary tests = 0.072, so not "detected").
   Under the plan's rule (whole interval inside 1/1.03..1.03) this is a **well-powered null for any effect above 3%**,
   not a proof of zero: the point estimate is positive, as predicted, and secondary cuts put it at 1.006-1.028.
3. **ENSO's own effect on HF-low frequency (gust version): none detectable.** Total T = 1.018 (0.915-1.079), direct with
   PNA held = 1.007 (0.912-1.066). With no total effect to mediate, the "share carried through PNA" is **undefined** in
   the gust version, by the plan's rule. The change T - D is +0.010 in log terms, which is the same indirect effect.
4. **Count versus share.** The interaction splits as all-cyclone count 1.013 (permutation p 0.056, q 0.16), share
   1.000. Nothing is carried by either.

## A depth-version contrast that needs care (secondary, not part of either primary test)

With "strong" meaning deep (a count-matched 965.0 hPa cut, 1979-80..2025-26, era term added), ONI does have a total effect
on Pacific deep lows: RR 1.066 (1.013-1.119) per SD, and 1.107 (1.033-1.168) for 2001-02..2025-26 alone. About a quarter
(0.25; 0.16 for 2001+) of it runs through PNA; the direct effect with PNA held is 1.049 (1.004-1.098). Gust-based HF counts
over 2001-02 onward show no such total effect (1.023, 0.937-1.089). So gust and depth disagree again on a Pacific question,
as in the frequency-split work. **Total-effect (T) tests are not in the plan's BH family** and are one of many; treat this as
a lead. Pressure-defined outcomes are also exposed to circularity with pressure-defined indices (lag is the guard).

## Phase table (S4, descriptive, untested)

HF genesis per 100 days, relative to ENSO-neutral and neutral PNA (|ONI| >= 0.5; PNA beyond +/-0.5 SD):

| | -PNA | neutral PNA | +PNA |
|---|---|---|---|
| La Niña (11 seasons) | 0.96 | 1.09 | 1.16 |
| neutral (15 seasons) | 0.91 | 1.00 | 0.98 |
| El Niño (8 seasons) | 0.84 | 1.02 | **1.25** |

El Niño with +PNA has the highest HF rate (23.9 per 100 days against 14.1 in the neutral cell), but La Niña with +PNA is
nearly as high (22.5), and the ENSO-neutral cells are lowest. The El Niño/+PNA cell is 1.25 times the additive
prediction (z +1.25), not distinguishable from additive. Only 8 El Niño seasons contribute. So "El Niño plus +PNA gives
more HF lows than neutral" is true in the raw table, but it is not evidence that the two reinforce each other.

## Power and effective n

ONI carries about one independent value per season, so n is 22 seasons (47 for the depth version), not the 825 events.
Detectable at 80% power: interaction RR 1.096 (primary), 1.057 (47-season depth); indirect effect RR 1.017. An
interaction as large as PNA's own main effect (about 1.07 per SD) is below the detectable size, so it could be missed.

## Multiplicity

BH family: every gamma and indirect permutation p (HF, all-cyclone and share) in all 9 cells, plus the asymmetry test:
55 tests. P1 and P2 also get Bonferroni x2. Lowest q in the family: 0.099 (indirect effect, several cells). Nothing
passes q < 0.05. The phase table and T/D estimates are outside the family.

## What this does not show

- Linear, one-product interaction on standardised indices. A threshold or phase-dependent mechanism (PNA only matters
  once the jet is extended by El Niño, say) could be missed; the phase table is the only look at that and is thin.
- Daily-genesis Poisson counts: the response is where/when storms start, not their whole life.
- Lag choice (days -10..-4) is the primary; the same-time index gives larger PNA effects (b 1.188) and an indirect
  effect of 1.028 (1.002-1.050; secondary), partly because storms feed their own index.
- Pacific only. The Atlantic NAO x ONI pair is in the additive test (null).
- 14 index days (CPC missing) were filled with the standardised mean so every season keeps its full calendar.

## Departures from the plan

None changes a pre-registered estimate. Logged because they happened after the plan was committed.

1. The plan's permutation statistic for the share terms is the point estimate (not a Wald); for the HF and all-cyclone
   terms it is the clustered Wald, as planned.
2. A 20-draw smoke run of the full script was done before the 2,000-draw run, and its output was not read beyond
   confirming it ran and wrote files.
3. The plan's item S7 ("not run") is a stray line; there is no S7.
4. S3 (asymmetry) and S4 (phase table) were run on the primary window only.

## Verification

A fresh Sonnet agent that had not seen this code recomputed with its own implementation, from `all_tracks.csv.gz`, the CPC
files and the ONI json, and matched: the primary-window counts (10,040 cyclones, 825 HF, 4,620 days, 14 days with a
missing PNA lag); the HF interaction RR 1.014 (se 0.0327); ONI-alone RR 1.018, direct ONI 1.007, PNA 1.071 (ses
0.0332 / 0.0307 / 0.0329); the A-path +0.150 (se 0.049, season-clustered only); the indirect effect (log +0.0103); the
all-cyclone interaction 1.013 (se 0.0053); the Pacific depth cut (965.0 hPa, 831 tracks; the count depends on the
genesis-day window, six tracks tie at 965.0) and the 47-season depth ONI-only 1.066, direct 1.049, PNA 1.131,
interaction 1.006; and the El Nino/+PNA, El Nino/neutral and neutral/neutral cell sizes (569/136, 618/114, 519/73).

**Not independently checked:** bootstrap intervals, all permutation p-values and BH q-values, leave-one-season-out
ranges, detectable-effect figures, the share (stacked) terms, S1, S2, S3, the 1979-2000, 2001-2025 and gust-era-term
cells, the 'excess over additive' columns of the phase table, and the 'share mediated' values (0.25, 0.16).
