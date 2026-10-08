# What turns cyclones into HF lows? Pre-registration (agenda RA-17, RA-18 with RA-6)

Written and committed **before any cyclone-level outcome was split by deepening or joined to a pattern index**. ERA5 proxy
throughout; outcomes are **pipeline A** (`research/era5/hf_history`, 800 km ocean gust index, HF-equivalent at 71.7 kt), a
proxy for the archive, not the archive.

Questions (Jason's thread "What turns cyclones into HF lows"):

1. **RA-17.** The share of cyclones reaching HF is P(deepening) x P(HF | deepening). PR 47 found the seasonal rise sits in the
   second factor. Which factor carries the share effect of the PR 41 hemispheric pattern found in PR 64?
2. **RA-18 (RA-6 folded in).** How much of the pattern's effect on P(HF | deepening) runs through the storm's near
   environment (jet, Eady growth, SST gradient and the other PR 47 ingredients)?

Answers will be given in plain words per basin: yes, no (well-powered), or can't tell, and by how much.

## What had and had not been looked at

- Looked at (structure and marginals only, nothing against any index): the column layouts of `all_tracks.csv.gz`,
  `fixes_2004.csv.gz`, `env_2004.csv.gz`; for the 19,928 Oct-Apr cyclones of 2004-05 to 2025-26, the two-way count of
  (HF, has a deepening fix) by basin (Atlantic: HF 901 with, 106 without a deepening fix; not HF 3,661 with, 5,111 without.
  Pacific: HF 789 with, 37 without; not HF 4,198 with, 5,125 without), and the index (row number) of each cyclone's first
  deepening fix (1,252 at the first fix, 9,343 at the second, then a tail). That is why the HF-without-deepening-fix group
  is a named third factor below.
- Known from other threads: PR 64 (weekly counts, channels, f 0.94 / 0.79), PR 47 (seasonal rise in P(HF | deepening):
  Atlantic 13.3% to 27.6%, Pacific 14.3% to 20.2%; ingredient attribution 0.30 Atlantic, Pacific collinear), PR 56
  (eight further predictor groups add +0.0045 BSS), PR 63 (Greenland high).
- **Not** looked at: weekly counts of deepening cyclones, or of HF-and-deepening cyclones, against any index; any
  ingredient against the pattern index; any ingredient against HF within a month.
- **Looks at the held-out seasons.** The pattern index is the leave-one-season-out index of PR 41 over all 22 seasons
  (`hemispheric/results/oos_index_*.csv`), which includes 2015-16 to 2025-26 (already scored six times, the sixth being PR
  64; see `hemispheric/results/heldout_looks.log`). Part 1 is **look 7** and Part 2 **look 8** (same index and seasons,
  different outcomes, counted separately to be conservative). Both are decompositions of an effect already found, not a
  search. They are appended to the log when run.

## Data and units (both parts)

- **Seasons.** 2004-05 to 2025-26 (label 2004..2025), Oct-Apr (Decision 1: fit and test 2004-05 onward). Weekly design as
  PR 41/64: 30 weeks of 7 days from 1 October, 660 weeks per basin, effective n 22 seasons.
- **Cyclones.** Every pipeline A track in `all_tracks.csv.gz`, Atlantic or Pacific, placed in the week of its first fix.
  Transitioning tropical cyclones are in. Basins are never pooled.
- **HF lows.** `gust800_kt >= 71.7` (gust primary; fixed depth 966.2 / 965.0 hPa is S1).
- **Deepening cyclone** (D): the track has at least one 00/12 UTC in-domain fix in `fixes_2004.csv.gz` with `dp12 <= -3.6`
  hPa (the 12 h pressure change at the fix, the PR 47 rule, fixed there before its frequency was seen). This needs a fix and
  its predecessor, not a 24 h future, so rapid-decay truncation does not apply, but a track with a single in-domain fix
  cannot be deepening and counts as not deepening. The flag is a property of the whole track, so a cyclone that deepens
  only after it is already HF is still D; this is stated as a limit, not corrected.
- **The three counts.** D (deepening), HD (HF and deepening), H (HF), N (all). Cyclones that are HF without a deepening fix
  (H - HD) are the residual third factor.
- **Pattern index.** `idx` of `oos_index_{atl,pac}.csv`, standardised over the 660 weeks of each basin; effects are per
  +1 SD. Not refitted; no new pattern searched.

## Part 1 (RA-17): which factor carries the share effect

Weekly Poisson models exactly as PR 64 (month effects Oct..Apr, plus the log1p count of the same outcome in the previous 7
days, plus idx; season-clustered sandwich; the stacked joint SE of `chanlib.joint_se`; `chanlib.py` reused unchanged).
Outcomes per basin: N, D, HD, H. With b_X the idx coefficient of outcome X:

- **A1** log RR of P(deepening) = b_D - b_N.
- **A2** log RR of P(HF | deepening) = b_HD - b_D.
- **A3** log RR of the residual (HF without a deepening fix, as a ratio to HF with one) = b_H - b_HD.
- **A4** log RR of the share = b_H - b_N = A1 + A2 + A3 (PR 64 T2; reproduction check: N and H must give PR 64's T1 and
  T2 to four digits, else stop).
- **Fractions** f_conv = A2/A4, f_deep = A1/A4, f_res = A3/A4 (sum to 1), reported for each basin with 95% season-block
  bootstrap intervals (2,000 draws, seed 20261011). Defined only when A4's interval excludes 0 (it did in PR 64).
- **Tests.** Primary family: A1, A2, A3 x 2 basins = 6 tests, p by season-block permutation of idx (2,000 permutations,
  PR 14/64 scheme, floor 1/2001), Benjamini-Hochberg q within the 6. A4 is a reproduction, not a test.

**Prediction P17 (agenda).** At least 70% of the log effect is in P(HF | deepening): f_conv >= 0.70.

**Decision rule, per basin.** conversion carries most if the whole interval of f_conv is above 0.5; deepening carries most
if the whole interval of f_deep is above 0.5; otherwise mixed / unresolved. P17 is *supported* if the verdict is
"conversion carries most" and the point estimate of f_conv is >= 0.70, *half right* if the verdict holds but the point
estimate is below 0.70, *falsified* if the interval of f_conv lies wholly below 0.5. P(deepening) is a **well-powered null**
only if the 95% interval of exp(A1) lies within 0.95-1.05; otherwise inconclusive, with the minimum detectable effect
(2.8 x clustered SE of A1, reported as a rate ratio).

**Secondary (own BH families, labelled secondary, none changes a primary):**

- S1 fixed-depth HF cuts in place of gust (A1-A3, fractions).
- S2 deepening thresholds -2.4 and -6.0 hPa per 12 h in place of -3.6 (A1-A3, fractions).
- S3 other indices, same weekly design, window days -7..-1 before the week (PR 41 "lag1") computed from
  `docs/data/teleconnections.json` (NAO, PNA) and `nao_share_barrier/gh_daily.csv` (Greenland high, the PR 63 series):
  Atlantic NAO and Greenland high, Pacific PNA. These are separate indices, not the PR 41 pattern, and each is its own
  look at 2004-2025 seasons with the stated limits of PR 14/63. Fractions as above; no verdict is drawn for an index whose
  A4 interval includes 0.
- S4 halves: 2004-14 and 2015-25 (not independent of the leave-one-season-out index).
- S5 Atlantic: leave out tracks with peak latitude above 60N.
- **Not possible without a new pull: a 1979-2000 replication.** The pre-2004 deepening flag needs per-fix pressure histories;
  `all_tracks.csv.gz` has only track minima. Stated as a limitation, not attempted.

## Part 2 (RA-18, with RA-6): mediation by the near environment

**Population.** Deepening cyclones (D = 1) in the weekly window, per basin. Outcome: HF (gust). Pattern value x = idx of the
cyclone's first-fix week. Month = month of the first fix (7 dummies). Ingredients z (standardised within basin on the
analysis sample), from `env_2004.csv.gz` at a reference fix:

- **Primary eight (PR 47):** jet250, eady, sstgrad, sst_t500, flux, tcwv, div300, vadv500. Secondary set (S): lat, B, VTL,
  VTU, sst.
- **Reference fix.** *Primary (R_B, lagged):* the 00/12 UTC fix immediately before the first deepening fix, 12 h earlier,
  so the environment is not the storm's own imprint at the point of deepening (the guard asked for by PR 47 and the
  agenda). Cyclones whose first deepening fix is their first (no earlier fix) are dropped from R_B. *S1 (R_A):* the first
  deepening fix itself (PR 47's primary), which keeps the whole population and is expected to lean towards more mediation.
  Every mediated fraction is computed on one common complete-case sample (all eight ingredients present), total and
  direct effect on the same cyclones. The sample size and dropped counts are reported.

**Quantities (per +1 SD of the pattern, risk scale).**

- **Total effect, beta_tot:** coefficient of x in logit P(HF | D) = month + x (complete-case sample); p by season-block
  permutation of idx, 2,000 permutations; one BH family over (beta_tot x 2 basins).
- **Path a_j:** coefficient of x in OLS z_j = month + a_j x (j over the eight). 16 tests, bootstrap p, own BH family.
- **Path b_j (RA-6):** odds ratio per SD of z_j in logit P(HF | D) = month + z_j (month fixed, so the seasonal cycle is
  removed and this is a within-month anomaly effect; x is not in the model), 16 tests, bootstrap p, own BH family.
  Prediction (agenda RA-6): Atlantic jet250 and eady have OR above 1.05.
- **Mediated fraction M (primary), g-formula.** Outcome model: logit P(HF | D) = month + z(8) + x, logistic with a tiny
  ridge. Natural indirect effect NIE = mean_i [ p(z_i + Delta, x_i) - p(z_i, x_i) ], with Delta the vector of the eight
  fitted a_j (the environment shift that +1 SD of the pattern brings); natural direct effect NDE = mean_i [ p(z_i, x_i + 1)
  - p(z_i, x_i) ]. M = NIE / (NIE + NDE). Also reported: the log-odds difference method 1 - beta_dir / beta_tot
  (beta_dir from the outcome model above, beta_tot from the month+x model on the same sample), which is not collapsible and
  is a secondary view.
- **Single-ingredient and unique parts:** M_j with only z_j in the outcome model; M_unique,j = M(all 8) - M(all but j).
  Parts do not add up (ingredients correlated); reported with intervals, not tested.
- Intervals: 95% season-block bootstrap (2,000 draws, seed 20261011) of M, M_j, M_unique,j and of every a_j, b_j (a
  draw refits all models).

**Prediction P18.** M >= 0.5 (point estimate). Also fixed now: if M < 1/3 the effect lives in something these fields do not
hold (mesoscale structure, timing, moisture sources, or the proxy's own limits), consistent with the small gains of PR 56.

**Decision rule, per basin.**

0. If beta_tot has q >= 0.05 (permutation), the effect on P(HF | D) is not resolved in that basin; M is reported but not
   interpreted.
1. Otherwise: **mostly mediated** if M >= 0.5 and the lower interval bound is above 1/3; **little mediated** if M < 1/3 and
   the upper bound is below 0.5; **partly mediated / unresolved** otherwise (state M and its interval).
2. A bootstrap draw with beta_tot below 0.01 in log-odds is flagged; if more than 5% of draws are flagged the fraction is
   called "unstable" and no verdict is drawn.

**Secondary (own BH families):**

- S1 reference fix R_A (the first deepening fix).
- S2 secondary ingredient set (lat, B, VTL, VTU, sst) added to the primary eight.
- S3 fixed-depth HF in place of gust (Pacific check of PR 52; both basins).
- S4 Atlantic leaving out peaks above 60N.
- S5 no month dummies, instead the seasonal cycle left in (an upper bound for mediation; descriptive).

## Power and limits

- Part 1: report A4 interval first; minimum detectable effect for each factor as 2.8 x clustered SE (rate-ratio scale).
- Part 2: a planted-effect check on the decision rule: simulate HF outcomes from the fitted outcome model with the direct
  coefficient scaled so that the true M is 0, 0.3, 0.5 or 0.8, refit, and report how often each verdict is returned
  (`power_m.py`, 60 simulations per cell, 400 bootstrap draws per simulation; approximate).
- **What this cannot show.** Cause: the pattern is a one-week-ahead association in reanalysis and ingredients are
  reanalysis fields at 1.5 degrees read near an existing low. Mediation here is a statistical attribution under no
  unmeasured mediator-outcome confounding, which is not credible; "carried by the environment" means "accounted for by",
  not "caused through". The reference-fix lag guards against the storm's own imprint but not against the pattern's
  influence on the storm's later evolution. Pipeline A outcomes only; the proxy is not the archive. The HF definition may
  change if "Our HF counts vs published" changes it; S1 (fixed depth) bounds that sensitivity.
- The deepening flag is a whole-track property, so P(HF | deepening) carries a mechanical element (HF storms usually
  deepen). It is the PR 47 convention and is kept for comparability.

## Deviations (post hoc)

Logged after the plan was committed. None yet.
