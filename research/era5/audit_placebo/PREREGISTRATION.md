# Audit pack: placebo pattern indices, conversion versus depth, test register (ERA5 PROXY, pipeline A)

Written and committed **before any placebo index was fitted and before any depth-attenuation model was run**. It answers three
items of the panel review (artifact "HF Research Panel Review", 2026-10-11) with committed files only. No ERA5 pull (0 GB).
Outcomes are **pipeline A** (800 km ocean gust index, HF-equivalent at 71.7 kt), a **proxy**, not the archive. Transitioning
tropical cyclones are in. Seasons 2004-05 to 2025-26 (label 2004..2025), Oct-Apr, 30 weeks, 660 weeks per basin.

## What had and had not been looked at
- Looked at (structure only): columns of `hemispheric/results/weekly_table.csv.gz`, `hf_history/results/all_tracks.csv.gz`,
  `oos_index_{atl,pac}.csv`, the code of `hemispheric/attribute.py`, `hemlib.py`, `share_mechanism/part1.py`.
- Known (background): PR 85 results (P(HF | deepening) x1.252 Atl, x1.137 Pac; f_conv 0.94 / 1.06).
- **Not** looked at: any index refit on a non-archive outcome; any HF-on-index model with a pressure control.
- Held-out seasons: no new look is claimed for a pattern test. Every analysis uses leave-one-season-out over all 22 seasons, so it
  is a further look at 2015-26 (counted in `looks/`, as PR 85 did). This is an audit of an existing decomposition, not a new search.

## Part A. Placebo indices (critique C3)
Question: does the share decomposition (D-N, HD-D, H-HD) put the effect in conversion for *any* index with a share effect, or
only for the index trained on archive HF counts?

Index construction (approximation to PR 41, stated as such): from `weekly_table.csv.gz` rows with season >= 2004, the
unsupervised PC scores z1-z10 and u1-u10 (SST channels are left out; PR 41 showed SST adds nothing). For each outer
leave-one-season-out fold: penalised Poisson (`hemlib.fit_model`: month dummies, log1p previous-week count of the training
outcome, PCs penalised) with lambda chosen by inner leave-one-season-out deviance over `hemlib.LAMBDAS`; the index is the PC
linear predictor on the held-out season (`hemlib.parts`), standardised over the 660 weeks, as `chanlib.index`. The EOFs behind
the PCs are computed without an outcome, so using them whole-record does not leak an outcome.

Indices, per basin:
- **P0** trained on the archive weekly HF counts (`y_atl`, `y_pac`), the PR 41 outcome. Reproduction check, below.
- **P1** trained on weekly counts of all pipeline A cyclones (N).
- **P2** trained on weekly counts of deepening cyclones (D; dp12 <= -3.6 hPa/12 h, as PR 85).
- **P3** trained on weekly counts of strong non-HF cyclones (gust800 in [55, 71.7) kt). "Strong non-HF" is this audit's label;
  it is not the project's storm-force class.
- **P4** season-shuffled: P0's training outcome (and its previous-week count) with whole seasons permuted among seasons;
  100 draws, seed 20261011.

Decomposition: `share_mechanism/part1.py` machinery on the pipeline A weekly counts N, D, HD, H (primary HF definition),
per +1 SD, 22 seasons: A1 = D-N, A2 = HD-D, A3 = H-HD, A4 = H-N, f_conv = A2/A4. For P0-P3: 2,000 season-bootstrap draws
(no refit of the index inside the bootstrap; this understates the index's own noise and is stated in the report) and 2,000
permutations, verdict by `part1.verdict`. For P4 only point estimates.

Statistic for the shuffle null: S = |A2| / (|A1| + |A2| + |A3|), the share of the absolute effect that sits in conversion.
Compare S of P0 with the P4 distribution (basin-wise).

**Reproduction gate (decided before the run).** P0 is accepted as a faithful stand-in if its Pearson r with the committed
`oos_index` is >= 0.90 in both basins and its A2 rate ratio lies inside PR 85's published interval (Atlantic 1.189-1.331,
Pacific 1.074-1.223). If not, comparisons are made against the P0 refit only, the committed index is not used in Part A, and the
report says the approximation did not reproduce.

Readings fixed in advance (per placebo index P1-P3 and basin):
1. A4 resolved (bootstrap interval excludes 0 and `frac_boot_ok` >= 0.95) and f_conv lower bound > 0.5: **conversion dominates for a
   non-HF-trained index too.** The "pattern acts through conversion" localisation is then partly a property of the denominator and
   the label, and is reported as not specific.
2. A4 not resolved: the non-HF index carries no share effect. The share effect belongs to the HF-trained index, which supports
   specificity (but not mechanism).
3. A4 resolved and f_conv lower bound <= 0.5: mixed; reported as is.
Shuffle: if S(P0) is below the 95th percentile of the P4 distribution, the "most of the effect is in conversion" statement is
not distinguishable from what chance indices give.
Tests: 2 basins x 3 placebo indices x (A1, A2, A3) = 18 permutation p-values, BH over those 18; P4 and the gate are not p-tested.

## Part B. Conversion: depth reached versus gust at a given depth (critique C3)
Tracks with a deepening fix (dp12 <= -3.6), week of first fix, committed P0-equivalent index = the committed `oos_index` (as
PR 85), per basin. Logistic regression of HF (gust800 >= 71.7) on: M1 month + x; M2 + minp and minp^2; M3 + minp, minp^2 and a
5 degree latitude by 10 degree longitude cell of the peak position (fixed effects). Report the log-odds per SD on x in each, and
attenuation 1 - b(M2)/b(M1), 1 - b(M3)/b(M1), with a season-cluster bootstrap (2,000 draws). minp is itself moved by the pattern,
so this is a decomposition of the conversion effect, not an adjustment for confounding.
Readings fixed in advance: attenuation by M2 >= 0.5 with lower bound > 0.25: conversion is mostly "reached deeper". Attenuation
by M2 <= 0.2 or upper bound < 0.5: mostly gust at a given depth (mixing, geography, parametrisation). Otherwise mixed. The M3 minus M2
change says how much is position. Tests: 2 basins x 2 attenuation measures; reported with intervals, BH over the 4 intervals' two-sided
bootstrap p-values.

## Part C. Test register (critique C2)
Machine-readable only: every committed results CSV on the research branch with a p-value or q-value column, counted by directory.
Report the number of rows with a p-value, the number the file itself marks as passing q < 0.05, the number of distinct result
directories, and the expected number of nominal passes at q = 0.05 under a global null as a crude upper-level note. It does not
include tests reported only in prose or JSON, so it is a lower bound on the number of tests. It does not estimate a global FDR.

## Not done / limits (stated in advance)
- The index is not refit inside the bootstrap (nested bootstrap not done).
- P0 is an approximation to PR 41's index (no SST, whole-record EOFs).
- The ingredient error question (critique on mediation, planted-effect with noisy ingredients) is not in this audit.
- The archive-label conversion check is not in this audit.
