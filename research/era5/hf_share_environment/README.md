# What in the environment tracks the midwinter rise in the hurricane-force share?

Plan, committed before any ingredient was compared with any outcome: [PREREGISTRATION.md](PREREGISTRATION.md) (commit `33e7414`; departures: none).
Code: [analysis.py](analysis.py) (all tests, 2,000 season-block draws, seed 20261010), [make_report.py](make_report.py) (tables, verdicts, figure).
Results: [results/](results/): `summary.txt` (readable, every variant), `tests.csv` (every estimate, interval, p and q), `monthly_profile.csv`,
`ingredients.png`, `boot_primary_*.csv.gz` (primary bootstrap draws), `populations.json`.

Reproduce (about 4 min, no ERA5 pull; reads committed tables only):

    python3 -I research/era5/hf_share_environment/analysis.py . research/era5/hf_share_environment/results 2000
    python3 -I research/era5/hf_share_environment/make_report.py research/era5/hf_share_environment/results

**Pipeline A and a proxy.** Cyclones, events and the share are pipeline A (`research/era5/hf_history`, HF-equivalent gust index >= 71.7 kt). The environment is
read from the intensity framework's near-storm table (`research/era5/intensity/results/env_2004.csv.gz`; ERA5 1.5 degree). Seasons 2004-05 to 2025-26 (22),
genesis month October-April, 19,928 cyclones (19,905 have a 00/12 UTC in-domain fix and enter; 23 do not). The shares reproduce PR 45 (Atlantic January 15.2%,
Pacific December 11.4%). Transitioning tropical cyclones are in.

## Answer

The HF share rises because a larger fraction of the cyclones that **do** start deepening go on to reach HF, not because more of them start deepening, and the
ingredients tested **account for only part of that rise in the Atlantic and cannot be separated in the Pacific**.

- **Stage.** Atlantic share 7.1% (October-November) to 15.3% (December-February): +0.87 in log-odds. The chance that a cyclone has a deepening fix barely moves
  (+0.07, -0.00 to +0.15); the chance that a deepening cyclone reaches HF carries it (+0.91, 0.74 to 1.11; 13.3% to 27.6%). Pacific 7.8% to 11.0%, +0.38 (0.22 to 0.55):
  deepening +0.10 (-0.01 to +0.20), HF given deepening +0.42 (0.26 to 0.60).
- **Atlantic: baroclinic supply lines up, about a third of the rise in total.** Shares of the rise (explained fraction, 95% season-block interval):
  Eady growth rate 0.32 (0.24 to 0.43), SST gradient 0.21 (0.15 to 0.29), 250 hPa jet maximum 0.18 (0.10 to 0.26), 500 hPa vorticity advection 0.17 (0.10 to 0.25), surface heat flux
  0.10 (0.03 to 0.20). All eight ingredients together: 0.30 (0.13 to 0.50), so roughly half or more of the Atlantic rise is not in these fields. The seasonal change of each
  is small, +0.14 SD (jet), +0.39 (Eady), +0.24 (SST gradient), while storm-to-storm each is strongly tied to HF (0.6 to 0.8 log-odds per SD), so each accounts for a modest slice.
- **Atlantic: things that do not.** 300 hPa divergence does not change with season (+0.02 SD, -0.03 to +0.07). Column water vapour **falls** into winter (-0.69 SD) and
  air-sea instability (SST minus T500) **rises** while it relates negatively to HF within a period, so those two oppose the rise (F = -0.46 and -0.09); SST itself also falls (F -0.28, secondary).
- **Pacific: the baroclinic ingredients each look large but the joint model is uninformative.** One at a time: jet maximum 0.66 (0.39 to 1.15), SST gradient 0.66 (0.38 to 1.46), flux
  0.48 (0.22 to 0.99), vorticity advection 0.42 (0.28 to 0.70), Eady 0.41 (0.27 to 0.73). Together: 0.02 (-0.55 to 0.53). The reason is that column water vapour falls into
  winter (-0.51 SD) yet predicts HF within a period (F = -1.65 alone), and the ingredients are correlated, so their joint split is poorly determined. The jet maximum has the largest unique part (0.34, 0.18 to 0.61).
  The seasonal change in the jet is only +0.08 SD (0.01 to 0.14); its large share comes from a steep storm-to-storm link (1.19 log-odds per SD).
- **Pacific: sensitivity changes too (S7, descriptive).** The same environment converts to HF more readily in December-February than in October-November for SST minus T500 (+0.37, q 0.02),
  column water vapour (+0.21, q 0.02), Eady (+0.22, q 0.06) and jet (+0.24, q 0.09). That is a change in response, not in the fields, and a model with fixed response cannot absorb it. Atlantic: no such change (all q >= 0.42).
- **Position (secondary).** Pacific cyclones are 0.23 SD (0.13 to 0.31) farther south at their first deepening fix in December-February; latitude alone has F 0.35 (-0.02 to 0.96), consistent with PR 45's southward shift of first-HF latitude. Atlantic latitude: no (F -0.03).
  Jet latitude and the storm's place relative to the jet exit were not examined (not in the tables; no new pull).
- **The fall (S4, December-February against March-April)** is explained about as much: all eight together 0.36 (0.22 to 0.54) Atlantic, 0.32 (0.19 to 0.45) Pacific. Jet maximum 0.28 / 0.33; Pacific flux 0.49.

## Is it the environment, or the storm? (Limits)

The environment is read near a low that is already deepening, so moisture, flux, vorticity advection and divergence are partly the storm's own doing. The within-period
associations (0.4 to 1.2 log-odds per SD) are therefore predictive markers, not proof of cause, and are probably inflated. The earlier reading at each cyclone's first fix (S1) gives the same picture
(Atlantic joint 0.24, 0.11 to 0.38; Pacific joint -0.39, -1.07 to 0.10; the Pacific joint is unstable in sign), and a reading at the fix of greatest deepening (uses the track's future, S2) agrees (0.32 and 0.14). Fix-level onset model (S3, 500 draws): joint 0.38 Atlantic, 0.26 Pacific (-0.19 to 0.61).
"Explained fraction" is a statistical attribution; the ingredients are correlated, so the unique parts are smaller than the single-ingredient ones and do not add up.

## Verdicts under the pre-registered rule (primary, first deepening fix)

"Lines up" needs q < 0.05 on both seasonal change and within-period association, a positive product and an F interval above 0.

| ingredient | Atlantic | Pacific |
|---|---|---|
| jet250 | lines up (F 0.18) | lines up (0.66) |
| eady | lines up (0.32) | lines up (0.41) |
| sstgrad | lines up (0.21) | lines up (0.66) |
| flux | lines up (0.10) | lines up (0.48) |
| vadv500 | lines up (0.17) | lines up (0.42) |
| div300 | does not (no seasonal change) | cannot tell |
| sst_t500 | does not (opposes, F -0.09) | does not (opposes, -0.35) |
| tcwv | does not (opposes, -0.46) | does not (opposes, -1.65) |

For tcwv and sst_t500 the rule's label "does not line up" means they run the wrong way for the rise, not that they are unrelated to HF.

## Tests, multiplicity, power

Primary family: 32 two-sided p values (8 ingredients x 2 basins x {seasonal change, within-period association}); one Benjamini-Hochberg q. Every primary association has q <= 0.001
(floor of the 2,000-draw bootstrap); seasonal-change q <= 0.028 except div300 (Atlantic 0.43, Pacific 0.76). Secondary variants are separate BH families (`tests.csv`, column `q`);
fall-contrast seasonal changes (16 tests) and the interaction check (16 tests) are their own families. Explained fractions are reported with intervals, not tested. Effective n is 22 seasons; the intervals resample whole seasons.
F is a ratio and unstable when the denominator (the Pacific rise, 0.38 log-odds) is small: Pacific intervals run past 1 and the joint interval spans -0.55 to 0.53. The Pacific has about half the Atlantic's rise, so it
is the less informative basin.

## Departures and caveats

- No departures from the plan. The -3.6 hPa deepening rule was fixed before its frequency was seen (about 50% of cyclones reach it).
- S3 used 500 draws (p and q floors 0.002 and 0.004).
- Missing ingredient values: SST-based fields 12.5% (SST gradient 21%) of reference fixes (ice and coast). Primary models use complete cases per model; S6a on the common sample gives the same verdicts except Atlantic jet (seasonal change q 0.07: cannot tell) and unchanged joint.
- ERA5 is a proxy; nothing here is the archive. The share, not the archive's counts, is pipeline A's.
- Not examined: jet position, storm position relative to the jet, static stability profile beyond SST minus T500, 1979-2003 (gust-based, excluded by decision 1).
