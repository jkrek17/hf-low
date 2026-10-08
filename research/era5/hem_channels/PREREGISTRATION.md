# Which channel carries the hemispheric pattern's effect on HF lows? Pre-registration

Written and committed **before any cyclone-level outcome was joined to the pattern index**. ERA5 proxy throughout; the
outcomes are **pipeline A** (`research/era5/hf_history`, 800 km ocean gust index, HF-equivalent at 71.7 kt), a proxy for
the archive, not the archive.

Question (agenda RA-1, from Jason's "how does the hemispheric pattern work"): the lagged hemispheric Z500 and 250 hPa jet
pattern of hf-low PR 41 predicts weekly HF counts about a week ahead in both basins. Does it do so by (1) producing more
cyclones, (2) turning a larger share of them into HF lows, or (3) bringing more storms into the basin from upstream?
Method: the frequency split of hf-low PR 14 (`research/era5/freq_split`), on a weekly design.

## What had and had not been looked at

- Looked at (structure only, no outcome against the index): column layout of `all_tracks.csv.gz`,
  `intensity/results/fixes_2004.csv.gz`, `hemispheric/results/oos_index_*.csv`, `weekly_table.csv.gz`; the basin boxes in
  `hf_history/track.py`; the distribution of the age at the first in-domain fix (used to define "entrant" below; it does not
  involve the index).
- Known from other threads (background): PR 41 results (archive weekly HF counts, pipeline A HF counts S4, depth counts
  1979-2000 S5), PR 14 (NAO, PNA share/count split), PR 52 (gust vs depth), PR 29 (Atlantic fixes north of 60N).
- **Not** looked at: pipeline A counts of all cyclones, of entrants, or of HF lows split by origin against the pattern index;
  any position of cyclones against the pattern index.
- **Held-out seasons.** The index used here is the leave-one-season-out index over all 22 seasons (each season's index
  from a fit on the other 21). That includes 2015-16 to 2025-26, which PR 41 already scored at least five times
  (`hemispheric/results/heldout_looks.log`). This analysis is a further look at them (at least the sixth). Each look is
  appended to that log and counted in the README. It is a decomposition of an effect already found, not a new search, and it
  says so.

## Data and units

- **Seasons.** 2004-05 to 2025-26 (label 2004..2025), Oct-Apr, 30 weeks of 7 days from 1 October (as PR 41; 210 days).
  22 seasons x 30 weeks = 660 weeks per basin. Fit and test on 2004-05 on (decision 1).
- **Cyclones.** Every pipeline A track in `hf_history/results/all_tracks.csv.gz` with `basin` Atlantic or Pacific (a low
  below 1010 hPa for at least 24 h with at least two fixes in the basin box; Atlantic 30-67N, 262E to 10E; Pacific 27-67N,
  135E-240E). A track sits in the week of its first fix (`start`), the rule PR 41 used for its proxy counts. Transitioning
  tropical cyclones are **in** (not removed), as in PR 14 and PR 41. Atlantic tracks with a peak north of 60N are a
  lower-confidence group (PR 29); a sensitivity (S3) leaves them out.
- **HF lows.** `gust800_kt >= 71.7`. Pacific gust and depth disagree on which storms count (PR 52); gust is primary
  (archive-like counts, held-out HSS 0.61 against 0.48), fixed depth is S1 (965.0 hPa Pacific, 966.2 hPa Atlantic, the
  count-matched cuts from `freq_split`, not re-chosen).
- **Origin of a cyclone** (from `intensity/results/fixes_2004.csv.gz`, the first 00/12 UTC in-domain fix of the track,
  `age` = hours since the track's first fix): **entrant** if `age >= 12` (the low existed at least 12 h outside the basin
  box before it first appeared in it, so it came in from upstream or from the other side of a box edge); **local**
  otherwise (it was first detected in the box). The 73 tracks missing from the fixes table are counted as local and
  reported. Local tracks are split in S2 into observed genesis (first fix at or above 1000 hPa) and already mature at
  detection (below 1000 hPa); that split is descriptive.
- **Pattern index.** `hemispheric/results/oos_index_{atl,pac}.csv` column `idx` (the pattern model's linear predictor,
  produced out of sample for each season), standardised to mean 0 and SD 1 over the 660 weeks of each basin. All effects
  are **per SD of the index**. The index was fitted to the archive's weekly HF counts of the same basin; it is not refitted
  here and no new pattern is searched for.

## Models (all per basin; Atlantic and Pacific never pooled)

For a count outcome Y_i in week i: log E[Y_i] = month effects (Oct..Apr, month of the week's midpoint) + d x log(1 +
previous-week count of the same outcome) + b x idx_i. Previous week = the 7 days before the week starts (for week 0 these are
24-30 September, from the full track list). This is PR 41's baseline B1. Poisson, log link, season-clustered sandwich SEs.
The month effects and previous-week term are not penalised; there is no season trend (PR 41 had none; a season trend is S6).

- **RR(x)** = exp(b) per SD for each outcome: all cyclones N, HF lows H, entrants N_e, local N_l, HF entrants H_e, HF local H_l.
- **Share channel.** RR(share) = RR(H)/RR(N), estimated exactly as in PR 14: H and N stacked in one Poisson model with
  every term interacted with an HF indicator, so log RR(share) = b_H - b_N and its clustered SE accounts for H being a subset
  of N. Likewise the share within entrants (H_e vs N_e) and within local storms (H_l vs N_l).
- **Count-versus-share fraction.** f = log RR(share) / log RR(H) (PR 14), defined only when RR(H) is distinguishable from 1.
- **Entry contrast.** log RR(N_e) - log RR(N_l) by the same stacking: does the pattern change the entrant fraction.
- **Position.** Track-level OLS of `peak_lat` and `peak_lon` (Atlantic longitude unwrapped: values <= 10 get +360) on idx,
  month effects (month of the track's first fix) and no previous-week term, season-clustered SEs, for all tracks and for HF
  tracks. Degrees per SD.
- **Absolute decomposition (the "by how much").** Fit N_k and H_k (k = entrant, local) with the same Poisson model. Per
  week, s_k = mu(H_k)/mu(N_k). Shift idx by +1 SD in every week and sum over the 660 weeks:
  DH = sum_k [dN_k x s_k0 (count) + N_k0 x ds_k (share) + dN_k x ds_k (interaction)].
  Report each term as a fraction of total DH, and the entrant fraction of DH against the entrants' fraction of baseline HF.

## Tests and multiplicity

Primary family, **24 tests** (12 per basin): T1 RR(N); T2 RR(share); T3 RR(H); T4 RR(N_e); T5 RR(N_l); T6 RR(share within
entrants); T7 RR(share within local); T8 entry contrast; T9-T12 peak latitude and longitude of all tracks and of HF tracks.
p-values: season-block permutation, two-sided, 2,000 permutations (each season's 30-week index series moved to another
season, week-of-season kept; PR 14's scheme), floor 1/2001. Benjamini-Hochberg q within the six core tests (T1-T3 x 2
basins) and across all 24. "Supported" means q < 0.05 across all 24. Every one of the 24 appears in the README with its p
and both q values, as "x of 24 pass".

Intervals: 95% season-block bootstrap (seasons resampled with replacement, 2,000 draws, seed 20261008) for RR, f and the
decomposition fractions, plus leave-one-season-out ranges for T1-T3.

**Secondary (own BH family, labelled secondary, none changes a primary):**

- S1 fixed-depth HF cut instead of gust (T1-T3 and f).
- S2 local tracks split into observed genesis and already mature (descriptive RR(N)); previous-week term dropped (T1-T3).
- S3 leave out tracks with peak latitude above 60N (T1-T3), both basins.
- S4 halves: seasons 2004-14 and 2015-25 separately (T1-T3, f). Not independent (the index is leave-one-season-out), a
  consistency check, not a replication.
- S5 archive reference: the same model on the archive's weekly HF counts (`hemispheric/results/weekly_table.csv.gz`,
  `y_atl`, `y_pac`). The archive has no cyclone denominator, so only the total is available.
- S6 add a linear season trend to T1-T3.
- S7 within-era replication on seasons never used to fit the pattern: 1979-80 to 2000-01, frozen primary-split pattern
  (PR 41, fitted on 2004-14), index standardised over those seasons as in PR 41 S5, depth-based HF (fixed cuts) only
  because pre-2001 gust values drift (decision 1), T1-T3 and f, within-era use only. The second look at pre-2004 seasons
  (PR 41 S5 was the first) and logged as such. Check first that the S5 rate ratios per SD of PR 41 are reproduced for the
  depth-HF total (1.245 Atlantic, 1.158 Pacific); if not, S7 is reported as failed-to-reproduce and not interpreted.

## Predictions (falsifiable, fixed now)

- **P1, Atlantic: the share channel carries most.** RR(share) above 1.10 per SD and RR(N) between 0.95 and 1.05, as PR 14 found
  for NAO. Supported if the whole 95% interval of f lies above 0.5 and T3 passes. Falsified if the whole interval lies below 0.5
  or RR(N) lies outside 0.95-1.05 with the interval excluding that band.
- **P2, Pacific: both channels contribute** (f between 0.2 and 0.8), the PR 14 pattern for PNA, and the Pacific pattern is
  PNA-like. This is a weak prior; "share only" or "count only" is not excluded.
- **P3, entry, Atlantic: entrants respond more than local storms** (T8 > 0), because the Atlantic pattern is a North American
  trough with an extended jet that should export lows. No prediction for the Pacific. The entrant channel is called
  disproportionate if the bootstrap interval of (entrant fraction of DH minus entrants' baseline HF fraction) excludes 0.
- **P4, position: no prediction**; two-sided tests; power reported as degrees per SD.

## Decision rules

Per basin, with T3 (HF lows, pipeline A) as the total effect:

1. T3 does not pass (q >= 0.05): "total effect not resolved on the proxy"; no channel is assigned; the channels are reported
   anyway.
2. T3 passes. Use the 95% interval of f: wholly above 0.5 -> **share carries most**; wholly below 0.5 -> **count carries
   most**; otherwise **both contribute / not resolved**. The absolute decomposition then gives the amounts and the
   entrant versus local split. If the two methods point to a different majority channel, both are reported and the f rule
   stands.
3. A channel is a **well-powered null** only if its 95% interval lies within 0.95-1.05 (counts and shares) or within
   +-1 degree (position). Otherwise a non-significant channel is **inconclusive**, with its minimum detectable effect.

## Power

Planted-effect simulation (`power.py`): weekly counts drawn from the fitted baseline (idx effect set to the planted RR,
Poisson noise; PR 41 found the dispersion near 1.0), the same stacked model and a season-clustered test (CR1 with t on 21
df), 400 simulations per RR in {1.03, 1.05, 1.08, 1.10, 1.15, 1.20}. Reported: the RR with 80% power for each of T1-T8
and degrees for T9-T12 (OLS on simulated tracks). Cluster-robust power is an approximation to the permutation test.
**What the study cannot show:** cause (the pattern is a one-week-ahead association in reanalysis; storms feed back on the
fields); anything about the archive's own cyclone population (it does not exist); whether the pattern or the HF definition
(gust) is what moves the share (S1 and PR 52 bear on it); that the three channels are physically separate (an entrant
is an upstream cyclone, whose own genesis may be favoured by the same pattern). The index is fitted to the archive, so
effects on pipeline A are a transfer, not a refit.

## Deviations (post hoc)

None yet.
