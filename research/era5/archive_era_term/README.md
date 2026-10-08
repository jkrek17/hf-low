# Do results that score archive HF by fixes per event survive an era term? (RA-29)

Plan: [PREREGISTRATION.md](PREREGISTRATION.md) (committed first, `13f060e`; no deviations). Code: `era_refit.py`, `ta_summary.py`, `era_splits.py`. Results: `results/`.
Reproduce (no pull, about 6 minutes; needs statsmodels): `python3 era_refit.py {none|step|asc|trend} results/ta_<mode>.csv 2000 2000; python3 ta_summary.py; python3 era_splits.py`.
Data: archive (`docs/data/hf-lows.json`, `data/hf_lows/`), seasons 2004-05..2025-26 (decision 1). ERA5 is a **proxy**; pipeline A (`research/era5/hf_history`) wherever it appears. Era step = 2009-11-23, the end of QuikSCAT (RA-11, PR 105): about 5.5 seasons before, 16.5 after.

## Answer

**The archive half of PR 79 does not change.** With a before/after-November-2009 term added, all 16 archive cells (RR per SD of the hemispheric pattern index for HF lows at 1, 2 and 3 fixes, the brief and middle classes, and the paired contrasts) keep their sign, move by at most 0.032 in log (at most 0.53 standard errors), stay inside their original 95% intervals, and pass or fail FDR exactly as before (6 of 6 headline cells pass, both before and after; 9 of 16 tests pass across the whole list, both before and after). This is a well-powered result for these cells: intervals widen by at most 1.10 times and the detectable ratio (80% power) moves by at most 0.012. The reason is plain: the step lowers the number of sustained lows in a week, but it is nearly uncorrelated with the pattern index (r +0.08), so it shifts the level, not the slope.

| archive RR per SD of the pattern index | original (95%) | with era term (95%) | shift (log) |
|---|---|---|---|
| Atlantic, at least 2 fixes | 1.268 (1.181-1.378) | 1.274 (1.184-1.389) | +0.005 |
| Atlantic, at least 3 fixes | 1.351 (1.263-1.461) | 1.364 (1.268-1.488) | +0.010 |
| Pacific, at least 2 fixes | 1.215 (1.108-1.373) | 1.233 (1.124-1.388) | +0.015 |
| Pacific, at least 3 fixes | 1.275 (1.146-1.461) | 1.302 (1.167-1.491) | +0.021 |
| Atlantic, sustained / brief | 1.148 (0.979-1.328) | 1.159 (0.984-1.355) | +0.010 |
| Pacific, sustained / brief | 1.108 (0.968-1.330) | 1.144 (1.003-1.370) | +0.032 |

The sustained-versus-brief contrast is still not resolved (q 0.086 in both basins after, 0.151 and 0.203 before), so PR 79's verdict stays "no clear preference, leaning sustained". The two other sensitivity forms agree: a step at ASCAT-B (2012-09) and a season-linear trend give 16 of 16 unchanged (largest shift 0.040 and 0.015).

**The era term does see the step, in the counts.** Same pattern index, weekly rate after versus before: sustained lows (at least 3 fixes) -18% Atlantic [-36%, +4%] and -25% Pacific [-38%, -10%]; at least 2 fixes -10% [-23%, +3%] and -19% [-32%, -5%]; all HF lows -7% [-21%, +11%] and -9% [-23%, +11%]. Only the two Pacific sustained and two-fix intervals exclude zero. So any **level** of at-least-2 or at-least-3 counts in the archive needs the era term; the pattern's **effect per SD** does not.

**Other results that use archive fixes per event (descriptive era splits, 9 registered contrasts, BH: 8 pass, `results/splits_tests.csv`):**
- **Count reconciliation: the "at least 3 fixes" rule counts differently in the two eras.** Events per year with at least 3 fixes: Atlantic 24.3 before, 21.8 after (ratio 0.90 [0.69, 1.20]); Pacific 22.3 before, 19.9 after (0.89 [0.71, 1.15]). Events with at least 1 fix are flat (ratio 1.03 [0.83, 1.38] and 1.01 [0.84, 1.27]), as RA-11 said. The share of events with at least 3 fixes falls 7.0 points (Atlantic 0.541 to 0.471, q 0.034) and 6.8 points (Pacific 0.578 to 0.510, q 0.013). Von Ahn's published means (Atlantic 20.0, Pacific 19.7 a year, QuikSCAT years) are therefore matched by the at-least-3-fixes count to within +21% and +13% in the QuikSCAT-era archive and +9% and +1% in the later archive. The rule is still not uniquely identified; the point is that the rule's count is era-dependent by about 10%, and "about 20" should be read with that.
- **Atlas hours at HF per event.** The archive mean falls about 2 hours after the step (Atlantic 19.6 to 17.6 h, -1.98 [-3.36, -0.55]; Pacific 19.6 to 17.6 h, -1.97 [-3.43, -0.57]); the proxy's mean does not (+1.2 [-0.6, +3.6]; +0.3 [-1.5, +2.0]). The archive-minus-proxy mean gap goes from +2.6 and +1.9 h before to -0.6 and -0.4 h after (change -3.2 [-5.4, -1.4], -2.3 [-3.8, -0.9]). So the atlas's "archive medians 12 and 18 h against the proxy's 12" is not a lasting difference in duration: the archive and proxy means agree after 2009 and differ only before it. The Pacific archive median stays 18 h in both eras (proxy 12 h); that is a property of the 6-hour steps in a skewed distribution, not of the mean.
- **Probability-track listed share (proxy events matched to an archive HF event).** Not lower after the step, which is what a fixes-per-event artefact would predict: Atlantic 0.577 to 0.597 (+0.019 [-0.039, +0.073]), Pacific 0.575 to 0.654 (+0.079 [+0.030, +0.126], q 0.007). The overall 61% hides a Pacific rise of eight points; whether that is the archive getting more complete, or the proxy, is not tested.
- **Skill limits: the archive-label ceiling is higher after the step, opposite to my prediction.** ERA5's hf24 label scored against the archive's "any HF fix in the next 24 h": HSS 0.599 before (40,133 fixes) and 0.675 after (119,297), difference +0.076 [+0.041, +0.107], q 0.006. The PR 12 forecast against the archive label goes 0.519 to 0.551 (+0.032 [-0.002, +0.062], descriptive). The headline 0.66 and 0.54 are pooled over both eras and are correct as pooled; the per-event step does not lower them. The step cannot be the cause (the share of archive-label positives is almost the same, 3.3% and 3.2%), and the data cannot say what is: anything that changed in the archive or the proxy near 2009 is confounded with it, and the pre-step sample is 5.5 seasons.

## Scorecard against the pre-registered predictions

Held: T-A all cells unchanged, interval widening under 1.3 (max 1.10), era coefficients negative for at least 2 and 3 fixes; T-B counts at least 1 flat and at least 3 below 1 (not resolved); T-C archive mean about 2 h lower, proxy not lower, gap narrower; T-E not lower. **Missed:** T-D predicted a lower HSS after the step and found a higher one. **Not predicted:** the archive mean falls below the proxy mean after the step (I expected it to stay above); the all-HF era coefficient was predicted about 0 and is -7% and -9% with intervals containing 0; the Pacific listed-share rise.

## What was checked and excluded

The list and reasons are in PREREGISTRATION.md ("Audit"). Re-fitted: PR 79 / RA-22 archive side (T-A), count reconciliation (T-B), atlas hours at HF (T-C), skill-limits archive label (T-D), probability-track listed share (T-E). Excluded: everything on counts or total fixes (PR 41 headline, weekly and basin counts, pair counts, clustering, seasonal cycle, `front_hours` S8) and everything on proxy fix counts (the proxy half of PR 79, lifecycle durations, tele_intensity HF duration, late wind, threshold axis, gust scaling), and ERA5-labelled P(HF) work (PR 12, 56, 68, 76). The audit read `research/era5/` for scripts that touch archive fixes and is not claimed to be exhaustive.

## Looks

None. No model was fitted on one block and scored on another; T-D reuses the stored leave-one-season-out forecasts and a fixed label. No entry in any looks log (2015-25 count stays 21).

## What this does not show

Why the step happened (RA-11 left it open) or whether the era term absorbs other changes near 2009-11. The pre-step sample is 5.5 seasons, so the era coefficients are wide; the pattern-effect estimates are not, because the era column is nearly orthogonal to the index. Interaction of era and index was not fitted (underpowered). The split designs (T-B..T-E) have unit = season-by-era segment, 6 before and 17 after; bootstrap p is two-sided with a floor of about 0.001-0.002.

## Verification

A fresh Sonnet agent recomputed every point estimate quoted for T-A to T-E (36 rows) and all matched; see [verify/VERIFICATION.md](verify/VERIFICATION.md). **Not independently checked:** all intervals, p and q values, the ASCAT-B and trend rows, the paired-contrast values, MDE and SE ratios, and the era-coefficient intervals. The malformed date `20241101018` (Pacific ID 2024202506) is already on the sheet-fix list from RA-11.
