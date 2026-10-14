# Storm-centred 500 nm pattern score and analogs

ERA5 proxy, **pipeline A** (research/era5/hf_history tracks and gust index), WeatherBench2 1.5 degree fields. Perfect-prognosis: ERA5 patches in, ERA5 gust outcomes out. Plan: `PREREGISTRATION.md` (committed before any field was pulled; deviations logged there). Run 2026-10-14.

Question: does the whole weather pattern in a 500 nm box around a low (500 hPa height, sea-level pressure and 250 hPa wind anomalies) add skill to the 24 h gale / storm / hurricane-force outlook beyond the existing 21-input model, and can past storms with similar patterns be used as analogs?

Population: 32,288 fixes (00/12 UTC) on 7,379 tracks, Oct-Apr of 2014-15 to 2021-22, 1,330 hurricane-force (HF) outcomes on 513 tracks. **8 seasons are the independent sample.** Leave-one-season-out throughout. Pull: 24.8 GB + 1.2 GB climatology.

## Results (outlook RPSS vs basin-month climatology; 90% interval by resampling the 8 seasons)

| Model | RPSS |
|---|---|
| BASE (21-input model, refit on the 8 seasons) | +0.395 |
| BASE + 15 pattern components (PC15) | +0.419 |
| BASE + PC15, patch with inner 300 km removed (S1) | +0.419 |
| BASE + composite projection (S2) | +0.398 |
| BASE + analog features (S3) | +0.403 |
| BASE + permuted patches (placebo S4) | +0.394 |

- **T1 (primary): PC15 over BASE = +0.0239 [+0.0220, +0.0258]. Pre-registered rule: "pattern adds skill".** It beats the prior I wrote down (+0.00 to +0.015). Gain is positive in all 8 seasons (+0.019 to +0.028), in both basins (Atl +0.026, Pac +0.022) and is largest for hurricane force (Brier skill +0.318 to +0.351, gain +0.033).
- S1 environment only: +0.0239 [+0.0217, +0.0259]. The gain does **not** come from the core circulation the state numbers already hold; it is in the surrounding pattern. FDR q 0.001.
- S2 composite projection +0.0033 (q 0.001), S3 analog features +0.0078 (q 0.001): real but much smaller than PC15. Placebo -0.0005, so the 15 extra columns cost nothing.
- All three secondaries pass FDR, but the intervals only reflect season-to-season variation of a pooled score; 8 seasons is a small sample and these 22 seasons have been looked at before (`heldout_looks.log`).
- **Post hoc check (descriptive):** the state model has latitude but not longitude. Adding position terms alone gives +0.0125 over BASE, so about half of T1 is "where the storm is". The pattern still adds **+0.0208 [+0.0184, +0.0229]** beyond BASE + position (post hoc, not in the FDR family).
- In plain rates: this is a modest gain in a score where BASE is already +0.395. Not a different class of skill, but not noise either.

## Analogs

- Stand-alone k = 30 analog probability of HF in the next 24 h: Brier skill +0.107 vs climatology (BASE logistic +0.318). Calibration slope 1.12; slightly underconfident at the top: of 364 fixes where 25-35% of the 30 closest past fixes reached HF, 34% did (`results/analog_reliability.csv`). So "6 of the 20 closest storms reached HF" does mean roughly 30-35%, but the logistic model's probability is the better number.
- The pattern similarity is about structure more than outcome: across the 1,330 HF-outcome fixes the five closest other-season analogs reached HF 16%, storm 42%, gale 26%, below gale 16% of the time (base rate for HF is 4%).
- `results/analogs_hf_fixes.csv.gz`: five closest other-season analogs for every HF-outcome fix. `results/case_*.png`: six lowest-pressure storms (chosen by pressure before looking at analogs) beside their five analogs.

## Audit (post hoc, requested after the result; `audit.py`, `results/audit.txt`)

The gain survives every check. Existing patches and the local gust cache only; no pull.

- **Fold hygiene:** scaling, PCA, composite, analog library, imputation and the climatology baseline are all fitted on the 7 training seasons inside each fold. The one external piece is the WeatherBench2 1990-2019 climatology behind the anomalies (it includes six of the eight test seasons). Using raw z500 and sea-level pressure (no climatology) with the PCA refit per fold gives +0.0227 [+0.0209, +0.0244] vs +0.0239.
- **Neighbour lows:** 37.8% of fixes have another tracked low within 800 km at t or t+24 h (own position at that time). Excluding them (20,092 fixes, 1,100 HF outcomes), retrained: +0.0226 [+0.0210, +0.0246]; outer pattern only +0.0230. On the excluded fixes alone +0.0238. Not tested: a tight gradient against a high with no tracked low.
- **Label radius (gust cache; the 800 km label replicates at 95% class agreement):** peak gust within 300 km (HF outcomes 985): +0.0277 [+0.0262, +0.0292], and +0.0268 with neighbour lows excluded. The gain does not come from far-away winds.
- **Where the RPSS gain sits:** 59% (+0.0140) at the gale / below-gale boundary, 30% (+0.0071) at storm vs gale, 12% (+0.0028) at hurricane force vs storm. Most of the score gain is at the gale boundary, not HF. Hurricane-force Brier skill still goes +0.318 to +0.351.
- **Intervals:** per-season gains +0.019 to +0.028 (SD 0.0035); leave-one-season-out pooled +0.0233 to +0.0246. The reported interval resamples the 8 seasons 2000 times; a t interval on the 8 season gains gives [+0.0216, +0.0262]. The interval is tight because the gain is nearly the same every season; it treats seasons as independent, takes the 15-PC / k = 30 design as fixed, and says nothing about the ERA5 proxy or a different 8 seasons.
- Fresh agent recomputed the per-threshold split, the neighbour flags and gains, and the season table; not independently checked: raw-patch variant and the 300 km label.

## Limits

ERA5 proxy, not validated against observations; ERA5 patches in, ERA5 gust out (a live product needs GFS patches and a tracker; the earlier full-model test kept ~85% of ERA5 gain on GFS, not tested for patches). 1.5 degree grid: synoptic pattern only. 8 seasons, one more look at data used before. North-up patches (not rotated to motion). The 15 PCs are not interpreted physically.

## Files

`common_patch.py` geometry; `select_fixes.py`; `clim_patch.py` and `extract_patches.py` (pull); `score_patch.py` (T1, S1-S4, calibration); `extras.py` (post hoc checks, analog display, product). Results in `results/` (`skill.txt`, `extras.txt`, `tests.csv`). Patch product, not in git (65 MB): `/mnt/project-files/storm-patches/` (see its README).
