# Does El Niño change what +PNA does to Pacific HF lows? (pre-registration)

Written 2026-10-08, before any model in this directory was fitted. ERA5 **proxy**, **pipeline A**
(`research/era5/hf_history`: 800 km ocean gust index, HF-equivalent at 71.7 kt). Not observations.

What has been looked at before this was written: column headers and first rows of
`research/era5/hf_history/results/all_tracks.csv.gz`, and the span and gaps of the index series
(ONI has no gaps 1979-2026; NAO/PNA lag series are missing on 14 days). What was *already known*
from earlier work and could bias choices: in the Pacific the lagged PNA effect on HF-low
frequency is RR 1.072 per SD and the lagged ONI effect is about 1.00-1.02 (`freq_split`, S5);
the archive's same-time effects are PNA 1.217 and ONI 0.989 (`teleconnection-test/REPORT.md`).
So a null *total* ONI effect is expected, and the plan below is built so that is not a problem
for the pathway question. ONI x PNA was never fitted by the additive test (pairs there: NAO x PNA,
NAO x MJO, NAO x ONI, PNA x MJO, ONI x MJO).

## Questions (Jason, 2026-10-08)

"Does +PNA together with El Niño produce more HF lows?"

- **(a) Interaction.** Does ONI change the size of PNA's effect on Pacific HF-low frequency? Is the
  combined effect bigger (constructive) or smaller (destructive) than the product of the separate effects?
- **(b) Pathway.** Does ENSO's influence run through PNA (El Niño shifts PNA positive; PNA drives HF
  lows)? ENSO's total effect, its effect with PNA held fixed, and the share carried through PNA.

## Hypotheses, stated before results

- **H1 (constructive interaction).** Jason's prior: El Niño amplifies the +PNA effect, so the ONI x PNA
  coefficient gamma is positive. The test is **two-sided**; a negative gamma (destructive) is reported as such.
- **H2 (pathway).** ENSO acts on HF lows through PNA: the ONI -> PNA path a is positive (El Niño favours
  +PNA; physically established) and PNA -> HF is positive (b about log 1.07 per SD), so the indirect
  effect a x b is positive. Prior expectation, stated honestly: a of 0.3-0.5 SD PNA per SD ONI would give an
  indirect RR near 1.02-1.04 per SD ONI, which is small. These hypotheses could fail, and a null is a
  full result.
- Against both: lagged ONI has shown no direct effect, and ONI carries about one independent value per
  season (effective n about 22), so power is low. An inconclusive result is a likely outcome and will be
  called that.

## Data and variables (fixed)

- Tracks: every pipeline A track in `all_tracks.csv.gz`, Pacific basin (majority in-domain basin), placed on
  its genesis day (first fix). HF = `gust800_kt >= 71.7`. "All cyclones" = every track.
- Primary window and seasons: Oct 1 + 210 days, seasons 2004-05..2025-26 (22 seasons). This is the rule that
  gust-based quantities are fitted and tested only from 2004-05 (decision 1).
- **Lagged indices (primary).** PNA = mean of daily CPC index over days -10..-4 before the genesis day (as
  `freq_split`); ONI = monthly ONI of the calendar month of day -7. Each standardised over the analysis days,
  so effects are per SD. Same-time indices are secondary (PNA days -3..+3, ONI of the genesis month).
- ONI source: `docs/data/teleconnections.json` plus `oni_1978_2002.txt`, as in `freq_split`.
- Code reuses `freq_split/split.py` helpers (index reading, windows, daily counts, Newton Poisson solver).

## Models (fixed)

Daily genesis counts per basin, Poisson GLM, log link, calendar-month fixed effects, linear season trend,
season-clustered (sandwich) SEs. Outcome y = HF count unless stated.

- **M0** y ~ ONI. Total effect T = b_ONI.
- **M1** y ~ ONI + PNA. Direct effect D = b_ONI, PNA effect b = b_PNA.
- **M2** y ~ ONI + PNA + ONI x PNA. Interaction gamma = coefficient on the product of the two standardised
  series (product not re-standardised). Per SD x SD.
- **A-path** PNA_lag ~ ONI_lag + month FE + trend, OLS on the same days, season-clustered. a = slope (SD per SD).
- Indirect effect (primary pathway estimand) = a x b (log RR per SD ONI). Difference method T - D reported
  beside it. Share mediated = a x b / T, **interpreted only if the bootstrap CI of T excludes 0**; otherwise
  reported as "undefined: no total effect to mediate".
- Count/share split: the same models for the all-cyclone count and the stacked share term, as `freq_split`
  (log RR(HF) = log RR(all) + log RR(share)). Done for gamma and for a x b. Cheap, so included; labelled secondary.

## Primary tests (two)

- **P1** gamma in M2, HF count, Pacific, lagged indices, 2004-05..2025-26.
- **P2** indirect effect a x b, HF count, same setting.

Inference: season-clustered z; 2,000 season-block permutations (P1: each season's block of both indices moved
together to another season against the counts, day-of-season kept, Wald statistic of the product term, the
additive test's method; P2: ONI blocks moved against PNA and counts, statistic a x b); 2,000 season-block
bootstrap draws (seasons resampled jointly for all quantities) for 95% intervals; leave-one-season-out.

## Decision rules (fixed now)

Smallest effect of interest (SESOI): RR 1.05 per SD x SD for gamma (about 70% of PNA's own log effect), RR 1.03 per
SD ONI for the indirect effect (the prior expectation above).

- **Detected:** permutation p x 2 (Bonferroni over P1, P2) < 0.05 AND the bootstrap interval excludes the null.
- **Well-powered null:** the whole 95% bootstrap interval of the RR lies inside (1/SESOI, SESOI).
- **Inconclusive:** anything else. Said in those words.
- Power is stated as the detectable effect at 80% power (2.8 x clustered SE, RR scale). Effective n = seasons x
  independent index values per season (ONI about 1, PNA about 5.7), not the cyclone count.

## Secondary and exploratory analyses (fixed now; labelled so)

1. **S1 same-time indices**, primary window.
2. **S2 Jun-May window**, 2004-05..2025-26.
3. **S3 El Niño versus La Niña asymmetry.** ONI split into N+ = max(ONI,0) and N- = min(ONI,0) (standardised
   after splitting); M2 with PNA x N+ and PNA x N- as separate terms. Because the question is about El Niño specifically.
4. **S4 ENSO-by-PNA phase table (descriptive, untested).** ENSO phase by |ONI| >= 0.5 (El Niño / neutral / La Niña),
   PNA phase by lagged standardised PNA beyond +/-0.5 SD (+PNA / neutral / -PNA). Poisson cell rates with month FE and
   trend relative to the neutral/neutral cell, with the additive prediction beside each. Reports seasons and days per
   cell, because several cells will be thin.
5. **S5 depth version, 1979-80..2025-26.** HF replaced by the count-matched pressure cut (Pacific, as `freq_split` S4),
   with an **era term** (season <= 2000) added to month FE and trend. Also 1979-80..1999-2000 and 2001-02..2025-26 on
   their own. Pressure depth did not drift before 2001; NAO/PNA are pressure-defined, so circularity is the exposure.
6. **S6 gust HF, 2001-02..2025-26** alone, and 1979-80..2025-26 with the era term. Within-era use only (decision 1);
   no level or trend claim. 2004-05 onward remains the primary.
7. **S7 ONI x PNA with MJO longitude held out?** Not run: out of scope.

All of P1, P2 and S1-S6 apply to the HF count; the all-cyclone count and share components are secondary to each.

## Multiplicity

Every fitted interaction and every pathway estimate is counted. Reported: raw p, Bonferroni x2 for the two primary
tests, and Benjamini-Hochberg q across all gamma and a x b tests (HF, all and share components; family size
stated in the results). S4 is descriptive and not in the family.

## Not in scope

New ERA5 pulls; other basins (the Atlantic NAO x ONI pair was in the additive test); other index pairs; awips-tools;
intensity given HF.
