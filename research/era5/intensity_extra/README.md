# Do more predictors improve the P(HF) model? (ERA5 proxy, pipeline A)

Pre-registered in `PREREG.md` (commit `b00e8f7`) before any feature was
extracted or scored. Code: `extract.py` (WeatherBench2 pull), `evaluate.py`
(leave-one-season-out). Output: `results/`.

**Proxy.** Outcome = pipeline A's ERA5 gust index reaching 71.7 kt within 24 h.
Not the archive, not a forecast-model test.

## Answer

**Nothing in Tier 1 adds as much as half a point of Brier skill, alone. All
eight together add about 0.0045 (0.42 to 0.424).** The PR 12 model is close to
what these fields can give it.

Sample: 130,271 fixes, 29,820 tracks, 18 seasons 2004-05 to 2021-22, 4,377 reach
HF within 24 h. Base (the PR 12 full model refitted on this sample) BSS 0.4192.

| Group | Gain in BSS, hf24 [90% over seasons] | q (primary family) | Verdict (rule: gain >= 0.005, interval above 0, q < 0.05) |
|---|---|---|---|
| P pressure environment (neighbouring high, ring gradient) | +0.0020 [+0.0006, +0.0034] | 0.10 | does not help |
| L low-level wind (850 hPa maximum) | +0.0013 [+0.0008, +0.0019] | 0.01 | real but tiny; does not help |
| K land fraction | +0.0005 [0.0000, +0.0010] | 0.15 | does not help |
| D diabatic (6 h precipitation) | +0.0003 [-0.0001, +0.0007] | 0.24 | does not help |
| A air-sea temperature difference | 0.0000 [-0.0002, +0.0002] | 0.93 | does not help |
| J jet geometry (250 hPa maximum position) | -0.0002 [-0.0004, 0.0000] | 0.15 | does not help |
| H history (gust and pressure tendency) | -0.0003 [-0.0006, 0.0000] | 0.15 | does not help |
| M moisture transport (IVT) | -0.0004 [-0.0007, 0.0000] | 0.15 | does not help |
| all eight, no selection | +0.0045 [+0.0028, +0.0059] | 0.01 | pre-registered combination; below the 0.005 bar by its own estimate, interval straddles it |

"Does not help" means the upper end of the interval is below +0.005, so a gain
of half a point of BSS or more is excluded (a well-powered null at this size).
The full table with the onset-only and rapid-deepening targets, all p and q, is
`results/skill.txt` and `results/gains.csv`.

Other targets (all in `results/`):

- **Onset only** (fixes not yet HF; base BSS 0.286): all eight together +0.0039
  [+0.0014, +0.0059]; no single group has q below 0.12.
- **Rapid deepening** (base BSS 0.281): all eight +0.0064 [+0.0042, +0.0087];
  moisture transport is the largest single (+0.0035, q 0.20), jet geometry
  +0.0016. No single group clears FDR.
- **Yes/no skill at the count-matched cut** (HSS, hf24, base 0.5855): all eight
  +0.0053 [+0.0017, +0.0088]; no single group's interval excludes zero. For
  rapid deepening the land fraction and jet geometry groups lower HSS
  (-0.006, -0.004, intervals below zero) while the Brier gain is positive or
  zero, so the HSS cut is noisy there.
- **By basin** (descriptive): the pressure environment group helps the Pacific
  (+0.0033, p 0.02) more than the Atlantic (+0.0010, p 0.33); 6 h precipitation
  helps the Atlantic (+0.0010, p 0.001) and not the Pacific (-0.0005). Neither
  survives a multiplicity correction across the 18 basin tests, and neither was
  an expected result.

## What this means

- **Surface and moisture quantities do not add.** Air-sea difference, moisture
  transport and precipitation add nothing out of sample, though the base
  already holds SST, SST minus 500 hPa temperature, the heat flux and column
  water vapour. The information they would carry is already there.
- **The gains that exist are small and sit in the current wind field.** The
  low-level wind and the pressure environment (how strong the background
  gradient is around the storm) add 0.001-0.002. These are close to the
  current gust index in what they measure, so they refine "how windy is it
  now" rather than add a precursor.
- **Jet geometry, history and land fraction do not help.**
- **Combined, the eight add 1 % of the base's skill** (0.4192 to 0.4236),
  better in 15 of 18 seasons. Real, but not a change in the forecast.
- Tier 2 (stability, trough depth, warm-conveyor ascent) has not been run: it
  needs about 88 GB, above the 50 GB gate, and waits for Jason's go-ahead.

## Checked

A fresh Sonnet verifier rebuilt the baseline and the LOSO from the committed
tables with its own code and reproduced the sample counts, the base BSS
(0.4192), the gains for L, P, D, K and all eight (to 4 decimals, intervals to
within 0.0001) and the onset-only figures. **Not independently checked:** the
p and q values, the groups H, A, M, J (only through the all-eight total), the
rapid-deepening and HSS tables, and the by-basin table.

## Not run, and limits

- `sel1` and `sel1_nested` were not built, as PREREG says, because no group met
  the "helps" rule.
- 2022-23 onward is not in the sample (WeatherBench2 ends 2023-01-09). Skill on
  those seasons is unknown.
- `airsea` is missing at 16,001 fixes (discs with no ocean cell); the median is
  imputed.
- IVT is WeatherBench2's integrated vapour transport, whose units are not
  documented in the store; it is standardised before fitting.
- The potential-vorticity candidate was dropped before the run (the array is
  empty); trough depth and jet geometry stand in for it.
- Features are on a 1.5 degree grid; a fine-scale quantity (fronts, sting jets)
  would not show here whatever the true value.
- Gain intervals describe season-to-season sampling of the skill, not
  uncertainty in the definitions. The quoted FDR q's treat the tests as
  independent; the candidates are correlated.

## Deviations and post hoc

- Before the run: PV removed (above, in `PREREG.md`).
- No candidate was added or changed after seeing results. The basin split and
  HSS table are descriptive and were listed in PREREG.

## Side finding

`../intensity/hart.py` could not be imported on Python 3.11 to 3.13: its header
added a second string literal ahead of `from __future__ import annotations`.
The header is now a comment (no other change). PR 12's committed tables were
produced before the header was prepended.

## Reproduce

```
python3 extract.py ../intensity/results/fixes_2004.csv.gz FEATDIR 6   # about 42 GB, 12 min
python3 evaluate.py FEATDIR results groups                            # about 8 min on 4 cores
```
`results/features_tier1.csv.gz` is the extracted table (`dg12` and `dpacc` are
computed in `evaluate.py` from the fixes file). The per-fix probabilities
(10 MB) are not committed.
