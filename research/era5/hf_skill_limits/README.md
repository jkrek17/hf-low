# What limits the skill of P(HF within 24 h)?

Plan: [PREREG.md](PREREG.md) (committed before any error table; its deviations log was written afterwards and says so).
Code: [analysis.py](analysis.py) (Q1-Q6, about 6 min), [diagnosis_extra.py](diagnosis_extra.py) (post hoc, about 4 min). Both read committed tables only, no ERA5 pull:

    python3 -I research/era5/hf_skill_limits/analysis.py . research/era5/hf_skill_limits/results
    python3 -I research/era5/hf_skill_limits/diagnosis_extra.py . research/era5/hf_skill_limits/results

Results: [results/skill_limits.txt](results/skill_limits.txt), [results/strata.csv](results/strata.csv), [results/diagnosis_extra.txt](results/diagnosis_extra.txt).

**Pipeline A, ERA5 proxy.** PR 12's `full` hf24 model, leave-one-season-out, 159,430 fixes, 22 seasons 2004-05 to 2025-26. HF means pipeline A's gust index >= 71.7 kt. Intervals resample seasons.
"About 50%" is read two ways: the PR 12 forecast (POD 0.60, FAR 0.40, HSS 0.59, onset-only BSS 0.29) and PR 47's 0.30 (see the end).

## Answer

The P(HF within 24 h) score is near the ceiling of what the resolved large-scale fields can say about the gust index, and only a small part of the gap is a missing ingredient.

| What | Score | Reading |
|---|---|---|
| PR 12 forecast, logistic | BSS 0.423, HSS 0.59 | the number being asked about |
| Same predictors, boosted trees (post hoc, one default run) | BSS 0.456, HSS 0.61 | model form is worth about +0.03 |
| Add PR 56's eight tier-1 field groups (post hoc, trees) | BSS 0.455, HSS 0.62 | nothing more, as PR 56 found |
| **Diagnose HF right now from the same fields, no gust, no forecasting** | BSS 0.437, HSS 0.60 | the forecast is already as good as a same-time diagnosis |
| Diagnose HF right now, tier-1 fields included (850 hPa wind, ring gradient...) | BSS 0.56, HSS 0.69 | best a nowcast of the gust from 1.5 degree-class fields gets |
| Forecast with the storm's realised 24 h deepening given (oracle, uses the future) | BSS 0.506 vs 0.436 on the same rows | knowing the future pressure fall is worth +0.07 |
| ERA5 hf24 label scored against the archive's HF fixes (same fixes, 2004-05 on) | HSS 0.66 | even a perfect ERA5-label forecast would be 0.66 against the archive |
| PR 12 forecast scored against the archive label | HSS 0.54 | |

Read the HSS column as nested limits, not a sum. Against the ERA5 label the forecast is at 0.59, a same-time diagnosis from these fields is at 0.60-0.70, and the label itself agrees with the archive at 0.66. The forecast loses roughly 0.1 HSS to storm evolution between now and 24 h (the best nowcast is 0.69-0.70; the oracle agrees in BSS terms), and the label difference against the archive costs about 0.05 more. The bands overlap and do not add.

1. **Most of the shortfall is that the gust index is not a function of the fields we can see.** ERA5 gust is a local maximum over 800 km that depends on structure the 1.5 degree environment table (and the tier-1 groups) do not hold. Diagnosing the current state, which needs no forecast at all, tops out at HSS 0.60 (0.69-0.70 with the tier-1 wind and gradient fields). The 24 h forecast is at 0.59 (0.62 with trees). The residual SD of a linear forecast of the 24 h peak gust is 12.3 kt overall and 8-9 kt for the fixes near the cut. This is the main limit. Whether the missing information is finer resolution of our predictors (the gust is computed on the native 0.25 degree grid, the predictors are 1.5 degree) or ERA5's own representation of the strongest winds was not separated; that would need a full-resolution pull.
2. **Storm evolution is a real but smaller piece.** Giving the model the realised 24 h deepening raises BSS by +0.070 (0.064 to 0.076), below the 0.10 the plan set. Environment explains 48% of the variance in the 24 h deepening rate (linear, LOSO). Score falls with lead for onset: BSS 0.232 / 0.163 / 0.108 / 0.070 at +12 / +24 / +36 / +48 h; HSS 0.41 / 0.33 / 0.25 / 0.21.
3. **Missing static ingredients are not where the shortfall is.** PR 56's eight groups add +0.0045; the tier-1 fields add 0.09 HSS to the nowcast but 0.004 BSS to the forecast, because the forecast already carries the current gust. The case for the Tier 2 pull (stability, trough depth, warm-conveyor ascent; about 88 GB) is weaker than before: the whole evolution slice they could act on is +0.07 BSS even with perfect hindsight, and tier 1 captured a few thousandths of it.
4. **The label adds a second, separate ceiling against the archive.** At the fix level, ERA5's own label agrees with the archive at HSS 0.66 (0.61 for fixes not already HF), so the 0.59 forecast of the ERA5 label becomes 0.54 against the archive. That is consistent with the event-level calibration (HSS 0.76 in the calibration seasons, 0.65-0.69 on 2006-2021). ERA5 minus archive central pressure is small (+0.4 hPa overall, +1.8 hPa in archive storms below 950 hPa), so the proxy damps only the deepest storms a little in pressure.

## Where the errors are (Q1; count-matched cut; 2,145 misses, 2,140 false alarms)

- **Onset is the problem, not persistence.** Fixes not yet HF are 98.2% of fixes and carry 81.4% of the Brier loss (90% CI 0.804-0.824; H2 met). Their BSS is 0.291; fixes already HF have BSS 0.681.
- **Near-threshold fixes carry the loss.** Fixes whose current gust is 55-71.7 kt are 12.9% of fixes and about 54% of the Brier loss; fixes already above the cut are 1.8% of fixes and 18.6% (storms that peak and decay within 24 h).
- **The misses are of two kinds.** about 41% of misses are storms below 55 kt now (POD 0.06 there): explosive onsets. 80% of misses peak within 6 kt of the cut, i.e. they cross it narrowly. Rapid deepeners (24 h fall above 1 Bergeron) are 2.2% of fixes and 33% of misses, 25% of the loss, POD 0.49.
- **The false alarms are different.** Their median peak gust is 63 kt, only 34% are within 6 kt of the cut, and 56% are on storms that never reach HF at any time (8% reach it 24-48 h later). They are deepening storms that stall short of the cut, not a threshold artefact.
- **Pre-registered tests:** H1 (at least 60% of errors within 6 kt of the cut): 0.52 (0.51-0.53) as is, 0.57 (0.56-0.58) with the 06/18 UTC cases placed on the cut: not met. H2 (onset at least 80% of the loss): 0.814, met. H3 (oracle gain at least 0.10): +0.070, not met. Holm over three leaves only H2.
- **Where skill is weak, with q across 24 strata (all reported, `strata.csv`):** Atlantic north of 60N HSS 0.51 vs 0.60 elsewhere (q 0.002), the Greenland barrier group; March-April 0.52 (q 0.002); storms in their first 24 h 0.50 (q 0.002); rapid deepening 0.42 (q 0.002). Basins do not differ (Atlantic 0.58, Pacific 0.60, q 0.33). May-September is better (0.66) with few events.

## Limits

- The 24 h peak gust uses only 00/12 UTC fixes; 28% of positives reach the cut at a 06/18 UTC point, which is why H1 is shown two ways.
- The knife-edge simulation (Q3) is an illustration, not a ceiling: its k = 1 point under-performs the logistic model (BSS +0.175 vs +0.291 on onset fixes). Read it as "halving the gust error roughly triples the onset skill in this simplified frame", and no more.
- Archive matching is nearest-point within 600 km (96% of archive HF fixes matched); a wrong neighbour lowers the 0.66. The ceiling is for fixes in pipeline A tracks, not for archive events the tracker never found (224 fixes).
- Boosted trees: one default run, 11-fold season-grouped, not LOSO. The diagnosis results are exploratory.
- Sting jets and frontal low-level jets cannot be examined with ERA5 alone. The fix-level archive pressure comparison above is the only resolution evidence here.

## Latent heating and moisture flux (Jason's question; Amendment 1, [moisture.py](moisture.py), [results/moisture.txt](results/moisture.txt))

Short answer: with the surface proxies we have, **no**: moisture matters for how fast a storm deepens, barely, and adds nothing to the wind. Warm-conveyor ascent aloft is untested.
- Proxies: column water vapour, surface heat flux, air-sea difference, integrated vapour transport (at 500 km and the maximum), 6 h precipitation as a surface stand-in for latent heating (tier-1 table, 130,271 fixes, 18 seasons, ERA5 proxy).
- H4 (matters through deepening, not through wind): moisture raises the R2 of the 24 h deepening rate by +0.005 (+0.005 to +0.006), under the 0.01 bar; it raises hf24 skill given the realised deepening by +0.001 (-0.000 to +0.002) and the plain forecast by +0.0005. Prediction half met (wind gain below 0.005); the deepening half not met, so H4 is not supported as written. Among rapid deepeners alone (2,834 fixes, HF rate 0.40) moisture changes BSS by -0.001 (-0.003 to +0.000).
- H5 (do misses among fast deepeners differ in moisture): as registered the test is confounded and uninformative. Misses have lower flux, IVT and precipitation than hits (standardised differences -0.6 to -0.7, AUC 0.905), but the hit/miss split is made by the model's own score, which already tracks those fields; with the score as a covariate the AUC is 1.000 by construction. The usable version is the rapid-deepener BSS gain above (none).
- What this does not say: precipitation and flux are surface proxies. Latent heating aloft (warm-conveyor ascent from 3-D vertical velocity) needs the Tier 2 pull (about 88 GB). Its upper bound: even perfect knowledge of the 24 h deepening is worth +0.07 BSS, and the proxies above recover a tiny share of that.

## Relation to PR 47 (the HF share environment, 0.30)

That 0.30 answers a different question: of the *seasonal rise* in the Atlantic HF share (log-odds +0.87 from October-November to December-February), the eight ingredients' seasonal changes account for 0.30 (0.13 to 0.50). It is not a skill score. Each ingredient is strongly tied to HF within a period (0.5-0.8 log-odds per SD) but moves only 0.1-0.4 SD between periods, so it can explain only a slice of the rise. The unexplained part is not evidence of a missing ingredient: the result of this note says the same fields cannot pin down the gust index even at a fixed time, so a joint fraction well below 1 is expected. This note did not re-run PR 47.

## Verification

A fresh Sonnet agent, without reading this directory, recomputed from the committed tables: sample counts, hf24 BSS 0.4228 and onset BSS 0.2915, misses 2,145 and false alarms 2,140 with POD/FAR/HSS, the onset and already-HF shares of Brier loss (0.8136 / 0.1864), the false-alarm track share (0.438), concurrent diagnosis (BSS 0.437, HSS 0.595), lead BSS at 12 and 24 h, the oracle (0.433 to 0.503 against the quoted 0.436 to 0.506: it differs by 0.003 in both, from its own sample handling), and the archive-label figures (HSS 0.657 and 0.543, 5,118 and 3,513 fixes, 96.0% matched). Two small differences were corrected here: 54% (not 53.5%) for the 55-71.7 kt band, and about 41% of misses below 55 kt (0.4075 recomputed). **Not independently checked:** the strata HSS intervals and q values, H1, the oracle's 12 h variant, Q3/Q3b, the Q6 pressure bias, the boosted-tree and tier-1 diagnosis (`diagnosis_extra.py`), and all of `moisture.py`.
