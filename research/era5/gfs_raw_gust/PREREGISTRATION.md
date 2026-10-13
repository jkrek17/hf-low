# gfs_raw_gust: pre-registration (written 2026-10-13, before any GFS gust is joined to an OPC label)

Question (Jason, "1, 2, 3", item 3): how well does GFS's own gust field near a low, taken raw with no ERA5-trained model, match OPC's
hurricane-force designations, and is it usable on its own next to the ERA5 gust proxy (71.7 kt, pipeline A) and the model probability?

Plain-words answer I will give: "usable on its own" (yes), "not usable" (no, well powered) or "cannot tell", per basin, from the criteria below.

## What had and had not been looked at
Looked at: the column layout, row counts and category names of `results/gfs_tracks_<season>.csv.gz`, `data/hf_lows/HF_Data_-_{Atl,Pac}.csv`
and `intensity/results/fixes_2004.csv.gz`; the finished `gfs_test` and `gfs_outlook` write-ups (basin-level A1 results, A2 model scores).
Not looked at: any per-fix join of GFS `g800` to an archive label, any HSS/POD/FAR/AUC of a raw gust, any threshold value.
No new data are pulled. Everything reuses committed GFS tracks (pipeline A detector and linker on GFS PRMSL and GUST; a PROXY detector, the GFS fix
is a forecast, not an observation). Pipeline A (ERA5) is a proxy record. GFS seasons 2021-22 to 2025-26 are also the seasons that calibrated the 71.7 kt
ERA5 index, so nothing here is independent of that calibration; stated with every result.

## Data and units
- Fixes: exactly `a_score.gfs_features()` (valid 00 and 12 UTC, leads 0, 12, 24, 36, 48 h, basin boxes, Oct-Apr cycles, seasons 2021-2025; June-May seasons).
  Unit = cycle x track x lead. Fixes of one low recur in many cycles, so the independent sample is **5 seasons**; all intervals are season-block bootstraps
  (2,000 resamples, seed 20261013, 95% percentile).
- Raw predictor G = `g800` (kt): maximum ocean gust within 800 km of the fix, GFS `GUST:surface`. Same quantity as the pipeline A gust index (`g800`)
  and as the state model's gust term. Sensitivities `g500`, `g1200` are reported but are not tests.
- Truth (OPC archive, Category == "HF" only; DHF/S/G excluded; malformed-date rows dropped as in `a_score.load_archive`):
  - **Y0 (primary)**: an archive HF fix within +-3 h of the valid time and within 600 km of the GFS fix.
  - **Y24 (secondary, for the model comparison)**: an archive HF fix in (t, t+24 h] within 600 km (the A2 secondary label).
- Basins separately (Atlantic, Pacific); pooled numbers are descriptive.

## Thresholds
- **Nominal**: G > 64 kt (no calibration), applied at every lead.
- **Count-matched, held-out seasons**: for each held-out season s, T_s is set on the other four seasons at lead 0 so that the number of lead-0 fixes
  with G > T equals the number of lead-0 fixes with Y0 = 1 in those four seasons (bias = 1), per basin. T_s is then applied to season s at the same
  lead 0 and, unchanged, at leads 12-48 (a forecaster sets the number on the analysis and reads it on forecasts). The matching seasons never include the scoring season.
  The five T_s values are reported. 71.7 kt (the ERA5 proxy value) applied to GFS is reported as a reference, not a test.

## Primary tests (12; Benjamini-Hochberg q across all 12)
Contingency scores: POD, FAR, bias, Heidke skill score (HSS), as `a_score.cat_scores`; season-block bootstrap on pooled counts.
- **R1** (4 tests): nominal 64 kt, Y0, lead in {0, 24}, per basin. Criterion "usable": HSS >= 0.30 and bootstrap lower bound >= 0.30 (the A1.3 rule).
- **R2** (4 tests): count-matched held-out T_s, Y0, lead in {0, 24}, per basin. Same criterion.
- **R3** (2 tests, per basin): GFS lead 0 (count-matched, R2) against the ERA5 proxy at 71.7 kt, both against Y0, same seasons, same months (Oct-Apr),
  same basin boxes, ERA5 pipeline A fixes at 00/12 UTC from `intensity/results/fixes_2004.csv.gz` (`g800`). Statistic: HSS(GFS) - HSS(ERA5), season-block
  bootstrap with the same resamples. Criterion "same": the 95% interval lies inside +-0.10. Caveat fixed now: 71.7 kt was count-matched on these same ERA5 seasons
  (in-sample for ERA5), GFS threshold is held out, so any tie favours GFS slightly less than it appears.
- **R4** (2 tests, per basin): at lead 24, label Y24, AUC of the frozen state model probability `p` (`gfs_model.json`, as scored in `a_score.py`) minus AUC of raw G,
  paired on the same fixes. Criterion "raw gust is not meaningfully worse": the upper 95% bound of (AUC_model - AUC_raw) <= 0.05. Note in advance that the
  model contains `g800` as a predictor, so the model cannot be worse in expectation; the test asks how much the model adds.

Decision rule per test: pass if the criterion holds. "Not usable" only if the upper bound of the interval is below the criterion (R1, R2) or the interval excludes it
the other way (R3, R4). Otherwise "cannot tell". With 5 seasons the intervals will be wide; a power statement is the observed interval half-width, and
a well-powered null is reported as a full result.

## Secondary, descriptive (not in the FDR family)
- S1 profile of POD / FAR / bias / HSS over leads 0-48 for both thresholds; the same with a lead-specific matched threshold (sensitivity).
- S2 AUC of raw G by lead, per basin, on Y0 and Y24; state model p AUC on Y24; HSS of model p at the same count-matched rule (matched on Y24 positives, held-out seasons).
- S3 event level (OPC storms whose first HF fix falls in Oct-Apr 2021-25): share with a GFS fix within 600 km at the 00/12 UTC valid time nearest the first HF fix (coverage,
  i.e. tracker misses count as misses); of those, share with G > T_s; and the earliest forecast lead (<= 48 h, any cycle) at which a fix valid at that time exceeded T_s. This is a loose
  analogue of the P >= 0.5 lead in the OPC-versus-model comparison (PR 161), which used ERA5 tracks, 2004-05 to 2025-26 and a different match rule, so the two are not the same experiment.
- S4 GFS minus ERA5 `g800` bias at matched fixes (existing: -0.5 kt, `gfs_test/README.md`).

## Deviations (post hoc)
1. The plan fixed the criteria for R3 and R4 but not the p-value that feeds the BH step. First scoring run (same script, same data) used a two-sided p (diff != 0) for R3 and a wrong-direction p for R4. Corrected after seeing it, before any write-up: R3 uses a two one-sided tests p (max of P(diff <= -0.10) and P(diff >= +0.10) over the bootstrap), R4 a non-inferiority p (P(diff >= +0.05)). The two-sided R3 p stays in `raw_gust.json` as `p_diff_outside` (0.001 in both basins). The decision per test follows the registered interval criteria and did not change.

## What this cannot show
Five seasons, one forecast model with undocumented version changes, a tracker that must first find the low, and an archive whose HF label is an analyst's
judgement. It says nothing about GFS gusts beyond 48 h, other models, or the open ocean away from lows. Not a test of the 64 kt value as a physical wind speed.
