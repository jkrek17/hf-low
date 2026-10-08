# What limits the P(HF within 24 h) skill? Pre-registration

Written before any error table, residual, oracle fit or archive comparison was computed. Looked at beforehand: only the headline
skill in `research/era5/intensity/results/skill.txt` (hf24 BSS 0.423, onset-only 0.291, HSS 0.59) and the PR 56 and PR 47 summaries.

**Scope.** ERA5 proxy, pipeline A (`research/era5/hf_history`, gust index g800 >= 71.7 kt). PR 12's `hf24` model (`full` predictor set, L2 logistic,
leave-one-season-out), 159,430 fixes, 22 seasons 2004-05 to 2025-26 (Decision 1: nothing earlier). The target is pipeline A's own label, so
"label noise" here means two different things and they are kept apart: (i) the ERA5 gust index is a knife-edge threshold on a continuous,
imperfect quantity; (ii) the proxy label differs from the archive. Effective n is 22 seasons; every interval resamples whole seasons (1,000 draws, seed 20261010).

## Questions, predictions and decision rules

**Q1 Where do the errors sit?** Strata fixed now: hf_now (already HF) vs onset (not HF at t); basin; calendar month group (Oct-Nov, Dec-Feb, Mar-Apr, other);
Atlantic north of 60N; distance of the current g800 below the cut (<55, 55-65, 65-71.7 kt); 24 h realised deepening class (PR 12's five classes); storm age (<24 h, 24-72 h, >72 h);
tropical-cyclone-linked (400 km of an IBTrACS point, as in the life-cycle thread, if the flag is in `lifecycle_events.csv`; otherwise dropped and said so).
For each stratum: share of all fixes, share of Brier loss, POD, FAR, HSS at the PR 12 count-matched cut.
- **H1 (knife-edge).** Among misses and false alarms at the cut, at least 60% have the 24 h peak g800 within 6 kt of 71.7. Reported with a season-bootstrap interval; supported if the lower bound is >= 0.50.
- **H2 (onset).** Onset fixes (not HF at t, 98% of fixes) hold at least 80% of the total Brier loss.

**Q2 Is it the storm's own future (evolution), not a missing static ingredient?** Oracle fits: the same `full` model plus (a) the realised 12 h msl change t to t+12, (b) the
realised 24 h normalised deepening rate (ndr24), added as predictors, leave-one-season-out. These use the future and are a diagnostic, not a forecast.
- **H3.** The oracle (b) raises BSS by at least 0.10 over 0.423. If so, most of the unexplained skill is how much and how fast the storm deepens, and a further environment field can only help through that route.
Also report: how much of ndr24 the environment explains (R2 and RPSS from PR 12) so the two steps are on one page.

**Q3 Knife-edge ceiling.** For onset fixes, regress the 24 h peak g800 on the `full` predictors (linear, LOSO). Report the residual SD in kt, then the BSS and HSS that a forecaster would score
if its forecast of the peak gust had Gaussian error of SD 2, 4, 6, 8, 10 kt (simulated on the observed peak-gust distribution, 20 draws). This says how accurate the gust amount has to be, not what is achievable.

**Q4 Lead.** Same model for HF at exactly t+12, t+24, t+36, t+48 h among not-HF-at-t fixes (LOSO BSS vs basin-month climatology and HSS). Prediction: skill falls with lead; no threshold is set.

**Q5 Archive label.** Match archive HF fixes (category HF, 6-hourly; DHF excluded) to pipeline A track points at the same time within 600 km, nearest. Archive label at an ERA5 fix = matched
track has an archive HF fix in (t, t+24 h]. Report (i) agreement of the ERA5 hf24 label with the archive label (POD, FAR, HSS: the ceiling for a forecaster that knew the ERA5 label exactly), (ii) the PR 12 model scored against the archive label,
(iii) archive HF fixes that match no track (missed by the tracker). Descriptive; no prediction.

**Q6 Resolution.** At matched archive HF fixes, ERA5 central pressure minus archive pressure, by archive pressure band. Descriptive. Sting jets, frontal low-level jets and mesoscale structure cannot be tested from ERA5 itself; stated as untestable here.

## Multiplicity and reporting
H1-H3 are the confirmatory tests (Holm over 3). Stratum HSS contrasts (each stratum vs the rest, season bootstrap) form one Benjamini-Hochberg family, all reported. Everything else is descriptive and labelled so.
Any change after the first outcome is logged below as a post-hoc deviation, and both versions reported. A well-powered null is a full answer.

## Deviations log
All logged after the first outcome was seen, in the order they arose.
1. **Peak gust window (before outcomes, a data limit).** The committed fix table holds 00 and 12 UTC only, so "24 h peak g800" is the larger of the t+12 and t+24 fixes. 28% of hf24 positives reach 71.7 kt only at a 06 or 18 UTC point and have a 12-hourly peak below the cut. H1 is therefore reported two ways (as is; those positives placed on the cut).
2. **H1 reading.** The rule text was ambiguous (60% point estimate, lower bound 0.50). Both are reported; the verdict below uses the stricter: the point estimate must reach 0.60.
3. **Oracle (b) adds ndr24 and its square, oracle (a) dmsl12 and its square**, and rows are restricted to fixes whose track survives 24 h with a class available (the base is refitted on the same rows: 0.436).
4. **Q3 post hoc.** The homoscedastic Gaussian version failed its own consistency check at k = 1 (BSS +0.087 against the logistic onset BSS +0.291). Q3b repeats it with a residual SD by predicted-gust band; k = 1 still sits at +0.175. Read the Q3 curves as an illustration of sensitivity to gust-forecast error, not as a ceiling.
5. **Q2c concurrent diagnosis** (HF now from the same fields without g800) and **diagnosis_extra.py** (boosted trees; tier-1 fields from PR 56) were not in the plan. Added after Q2 showed the oracle recovered only part of the gap. Exploratory. One run of default settings; no tuning; 11-fold season-grouped CV, not leave-one-season-out.
6. **TC stratum dropped** (the tropical-cyclone flag exists only for the 4,157 catalog events, not for all tracks).
7. **FDR.** The stratum contrasts are 24 overlapping strata (some are complements of others), BH q over all 24; descriptive.

## Amendment 1 (10:25Z, before any moisture result was computed): latent heating and moisture flux (Jason's question)
Proxies on hand (tier-1 table, PR 56, 18 seasons 2004-05..2021-22, 130,271 fixes): tcwv, surface flux, airsea, ivt500, ivtmax, precip6 (surface proxy for latent heating). Latent heating aloft (warm-conveyor ascent) is NOT available and is what Tier 2 would add.
Moisture group = {tcwv, flux, airsea, ivt500, ivtmax, precip6}; base = PR 12 full predictors. All LOSO (18 folds), BSS gain with season-bootstrap 90% CI.
- **H4 (route).** Moisture group improves prediction of the 24 h deepening rate (linear, R2 gain) by at least 0.01, but improves hf24 given the realised deepening (oracle rows) by less than 0.005 BSS: i.e. it matters through deepening, not through wind.
- **H5 (misses).** Among fixes with realised 24 h deepening above 1 Bergeron (rapid), misses differ from hits in the moisture group: jointly by a logistic of miss vs hit on the six standardised proxies, LOSO AUC; supported if AUC lower 5% bound > 0.55.
Both are reported even if null; two tests, Holm.
