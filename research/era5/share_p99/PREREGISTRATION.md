# RA-28: does a smoother HF index change the pattern's effect on conversion? Pre-registration

Written and committed **before any pattern index, weekly count or HF label from the new index met an outcome**. ERA5 proxy throughout,
**pipeline A** (`research/era5/hf_history`: lows from 6-hourly MSLP, gust index = maximum ERA5 instantaneous 10 m gust over the ocean cells a
low owns within 800 km). A proxy for the archive, not the archive.

Question (agenda RA-28, Jason's thread "HF threshold noise" and "What turns cyclones into HF lows"): RA-27 (PR 98) found the
single-cell maximum is a noisy crossing statistic and that the 99th percentile of the owned cells agrees better between two analyses
(HSS 0.80 against 0.61, for **sustained 10 m wind**). PR 85 found the hemispheric pattern (PR 41) raises P(HF | deepening) by x1.252
(Atlantic) and x1.137 (Pacific) per SD. **If HF is scored on the 99th percentile of the owned gust cells instead of their maximum, does that
effect keep its sign and rise, as label-noise attenuation would predict?**

This is a **sensitivity re-run, not a redefinition of HF**. The working definition (gust index >= 71.7 kt) is unchanged; any change of it goes
to Jason.

Answer will be given in plain words per basin: strengthens, keeps (well-powered), weakens, or can't tell, and by how much.

## What had and had not been looked at

- **Looked at (nothing against a pattern index or an HF outcome):** the column layouts of `all_tracks.csv.gz`, `fixes_2004.csv.gz`; PR 85's
  published Part 1 table (the headline numbers above, read to write this plan; the code is reused unchanged, see "Reproduction gate"); the number
  of 00/12 UTC fixes and distinct times at several fix floors (sizing); a **pilot** of 60 random Oct-Apr 00/12 fixes with catalog g800 60-85 kt
  (`results/pilot_fixes.csv`, `results/pilot_out.csv`): the re-derived g800 returned the catalog value for 60 of 60 fixes (max difference
  0.05 kt), and the owned-cell percentiles relative to the maximum were p99/g800 = 0.928 (sd 0.040, range 0.79-0.98), p98 0.899, p95 0.841.
  The pilot sets only the pull floor below.
- **Not looked at:** any percentile index of any track outside the 60 pilot fixes; any HF label from the new index; any weekly count of it;
  any of it against the pattern index; any ingredient.
- **Held-out seasons.** The pattern index is PR 41's leave-one-season-out index over all 22 seasons, which includes 2015-16 to 2025-26 (already
  scored many times; `../looks/*.log`, `../hemispheric/results/heldout_looks.log`). This run is **one more look** at that block, flagged as a
  sensitivity re-run of an existing decomposition (no new pattern is searched). Its entry goes to `../looks/share_p99.log` when the weekly
  design is first scored. It scores no season before 2001.

## Data, size and why only 00/12 UTC

The percentiles need the per-cell gust values of every low at every fix, which are not stored anywhere (pipeline A keeps only the maximum,
`g800`). **Re-extraction is needed.** The fix table that exists for 2004-05 on (`intensity/results/fixes_2004.csv.gz`) has the 00 and 12 UTC
in-domain fixes of every track, with the catalog g800 at each. Pipeline A's own 6-hourly fixes (06 and 18 UTC) are not stored, and a full
6-hourly pull of the same population would be about 80-100 GB, over the 50 GB gate. So:

- **The comparison is made on 00/12 UTC fixes, for both indices.** The maximum is re-derived on the same fixes (from `fixes_2004.csv.gz`, no new
  pull), recalibrated by the same rule, so the only thing that differs between the two labels is maximum against 99th percentile. The published
  6-hourly label (`gust800_kt >= 71.7`) is run alongside as the reference and must reproduce PR 85 (below). A track whose in-domain fixes are all
  at 06/18 UTC has no 00/12 label and counts as not HF under the 00/12 labels; the count of such tracks is reported.
- **Pull.** Fixes with catalog g800 >= **60 kt** (a percentile of the same owned cells cannot exceed g800; the pilot's percentile cut is expected
  near 66 kt, and the run stops if a calibrated cut is below 60 kt), for tracks in the share-mechanism population (first fix in Oct 1 + 210 d,
  seasons 2004..2025, as `chanlib.assign_week`) plus all months of the calibration seasons 2021-22..2025-26. `make_pull_list.py` lists
  12,772 fixes at 7,135 distinct 00/12 times. Sizing: gust 3.19 MB + MSLP 2.17 MB per time (mean of 40 HEAD requests) = 5.36 MB x 7,135 =
  **38.2 GB**, under 50 GB; no existing cache covers it. MSLP is needed to re-detect the lows and so reproduce pipeline A's ownership. The
  re-derived g800 must return the catalog value (gate below). Extraction is resumable and writes only small per-time CSVs.

## Labels (all pipeline A, ERA5 proxy)

For each track, the index is the maximum over its in-domain 00/12 UTC fixes of: **(M)** the maximum owned-cell gust, **(P99)** the 99th percentile
of the owned ocean cells within 800 km (unweighted by area, exactly the RA-27 S7 definition), **(P98)** the 98th percentile (secondary). Owned
cells are those pipeline A's ownership assigns to the low (nearest low, within 1200 km of one, ocean, within 800 km of this one). Fixes with
g800 < 60 kt are not pulled and enter P99/P98 as 0 (they cannot reach the cut).

**Calibration, once, before any outcome** (`calibrate.py`). Same recipe as the 71.7 kt cut: match archive events (`docs/data/hf-lows.json`) to
ERA5 tracks with `hf_history/calib.match` (400 km, 800 km for tip jets and centreless events), here on the 00/12 UTC fixes; pick the cut on a
0.1 kt grid whose forecast count equals the observed count in the calibration seasons 2021-22..2025-26 (bias closest to 1, the lowest such
cut on ties). Done separately for (M) and (P99) (and (P98)). The cuts are written to `results/cuts.json` and committed **before the weekly
design is run**. The calibration also reports POD, FAR, CSI, HSS, and bias for each index against the archive in those seasons (descriptive).
No pattern index and no cyclone-level outcome enters the calibration.

Labels in the weekly design: **L_ref** (published, 6-hourly, gust800_kt >= 71.7, as PR 85), **L_M** (00/12, max, its cut), **L_P99** (00/12, p99, its
cut), **L_P98** (secondary). The deepening flag D (`dp12 <= -3.6` hPa per 12 h, PR 47/85), the cyclone population, the pattern index and the
weekly Poisson design are PR 85 Part 1's, unchanged.

## Reproduction gate

The code is `conv.py`, a copy of `part1.py` at PR 85's head `67ddfe6` with the label column made a parameter and a paired bootstrap added (PR 85
is still open, so this PR does not depend on it). Gate: under L_ref the A1-A4 estimates, intervals and fractions must reproduce
`share_mechanism/results/part1_tests.csv` (primary) to four digits, else stop. The pull's gate: re-derived g800 returns the catalog g800 within
0.5 kt for at least 99% of fixes (RA-27's gate), else stop.

## Estimands and tests

Per basin (Atlantic, Pacific; never pooled), per SD of the pattern index (lagged days -10..-4 as PR 41), weekly Poisson with month effects, the
log1p count of the same outcome in the previous 7 days and the index; season-clustered sandwich, joint SE for contrasts; 22 seasons.

- **A2(L)**: log rate ratio of P(HF | deepening) under label L = b_HD - b_D. **Headline quantity.**
- **Delta** = RR_A2(L_P99) - RR_A2(L_M) (difference of rate ratios), and the log ratio, from **one set of 2,000 paired season-block bootstrap
  draws** (the same resampled seasons for both labels; seed 20261011, as PR 85). Percentile 95% interval, and two-sided bootstrap p = 2 min(
  P(Delta* <= 0), P(Delta* >= 0)) with the +1 convention.
- Secondary quantities per label (as PR 85 Part 1): A1, A3, A4, f_conv, f_deep, f_res with the same intervals; the number of tracks that flip
  between L_M and L_P99 among deepening cyclones, and agreement (kappa) between them.

**Primary family (4 tests, BH within it):** A2(L_P99) against 1 in each basin (season-block permutation of the index, 2,000 permutations, floor
1/2001), and Delta in each basin (paired bootstrap p). **Secondary families (own BH, labelled secondary, none changes a primary):**
(S1) the same four tests with L_P98; (S2) the A1-A3 x 2 basins x {L_M, L_P99} tests, 6 per label as PR 85; (S3) mediation, below.
**Across all pre-registered tests** a global BH q is also reported, with the count passing out of the total.

**S3: mediation (PR 85 Part 2 unchanged, label swapped).** Total effect beta_tot of the pattern on logit P(HF | deepening) and the g-formula
mediated fraction M through the eight near-environment ingredients (reference fix 12 h before first deepening; month fixed; complete-case sample;
season-block bootstrap, 2,000 draws), for L_ref, L_M, L_P99, each basin. `med.py` is `part2.py` at `67ddfe6` with the label parameterised. Gate: under
L_ref it reproduces PR 85's primary beta_tot and M.

**S4 (descriptive): agreement with the archive.** From the calibration: POD, FAR, CSI, HSS, bias of L_M and L_P99 against the archive, 2021-22 to
2025-26, with a season-block bootstrap interval of the HSS difference (5 seasons: coarse). Not part of any decision rule.

## Predictions and decision rules (written before any label is built)

Agenda prediction: same sign in both basins; the Atlantic RR rises from 1.252 toward 1.3 or more and the Pacific from 1.137 toward 1.2; if it stays
within about 0.03 of the headline, label noise is not what limits the pattern signal; a fall below the headline would say the headline is partly
produced by the extremes of the index.

Per basin, with Delta_RR = RR_A2(L_P99) - RR_A2(L_M) and its 95% bootstrap interval [lo, hi]:

- **Strengthens**: lo > 0. **Weakens**: hi < 0.
- **Keeps (well-powered)**: [lo, hi] lies inside [-0.06, +0.06] (twice the agenda's 0.03 pointer, so that the interval and not the point judges it).
- **Can't tell**: anything else; the minimum detectable Delta (2.8 x the bootstrap SE of Delta, in RR) and the point estimate are reported.
- Also reported, not decided on: whether RR_A2(L_P99) >= 1.30 (Atlantic) and >= 1.20 (Pacific); the point difference to the published headline
  (1.252, 1.137) and to L_M; whether the point estimate is within 0.03 of the headline.
- **Sign**: RR_A2(L_P99) > 1 in both basins and q < 0.05 within the primary family, else the sign claim is not made.
- The two basins are judged separately; no pooled verdict.

## Power and what this cannot show

- Power is the bootstrap precision of Delta (MDE = 2.8 x SE of Delta); a "keeps" verdict needs the whole interval within +-0.06, which needs
  SE(Delta) below about 0.03. If SE is larger, the verdict is "can't tell" however close the point estimate is to zero.
- Delta is a paired difference over the same storms, so its SE is smaller than the SEs of the two RRs, but the labels differ only for storms near the
  cut, so Delta may be small and hard to resolve.
- The 99th percentile of gust is not shown to be less noisy than the maximum: RA-27 showed it for sustained wind between two analyses, and the IFS HRES
  t0 analyses carry no gust. S4 asks only whether it agrees better with the archive. A rise of the RR would be consistent with, not proof of, label
  noise being attenuated.
- 00/12 UTC fixes only; seasons 2004-05 onward (Decision 1); Oct-Apr; transitioning tropical cyclones in; the Atlantic north of 60N is a
  lower-confidence group and is not separated here. The 22 seasons are one sample, shared with PR 85, so this cannot replicate PR 85.

## Deviations (post hoc)

None yet.
