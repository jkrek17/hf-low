# GFS test of the storm-centred pattern terms

Plan: `PREREGISTRATION_GFS.md` (committed before any GFS field was read; frozen ERA5 PCA and models in `results/frozen_*`). ERA5 pipeline A is a **proxy**: the training record and the truth. Question: does the pattern gain of PR 199 (ERA5 patches in, +0.024 outlook RPSS) survive when the patch comes from a GFS forecast? Code: `fit_frozen.py` (before the pull), `gfs_fixes.py`, `extract_gfs.py`, `score_gfs.py`; results `results/gfs_skill.txt`, `gfs_tests.csv`.

## What was done
PR 185's 18,167 matched lead 24 h GFS fixes (2,122 cycles, 5 seasons 2021-22 to 2025-26). GFS f024 HGT 500, PRMSL, UGRD/VGRD 250 (byte ranges, 6.39 GB pulled against 6.3 GB sized; 0 missing cycles), patches at the GFS tracker fix, anomalies against the same ERA5 1990-2019 climatology. ERA5-trained, **frozen** 15-component PCA and logistic models (FULL7 vs FULLPC7, trained on ERA5 2014-15 to 2020-21, nothing refitted on GFS). Patches saved at `/mnt/project-files/storm-patches/gfs/` (30 MB, not in git).

## Results (lead 24, RPSS vs basin-month climatology, GFS-fed; 90% season-block interval over 5 seasons)
- FULL7 +0.357; FULLPC7 +0.384; PR 185's FULL (2004-2020 fit) +0.351.
- **T1: pattern terms add +0.0273 [+0.0259, +0.0286]; rule "survives"** (>= +0.010, lower bound > 0). Prior written in advance: +0.005 to +0.020, so the gain came out larger than expected, and larger than the ERA5 leave-one-season-out gain (+0.024).
- S1 outer pattern only (inner 300 km removed): +0.0271 [+0.0258, +0.0283]. q 0.0005 for both (5 seasons: p values bottom out).
- By season +0.025 to +0.030; Atlantic +0.026, Pacific +0.029; on the 78% of fixes with no other GFS tracker fix (Atlantic/Pacific basin boxes, same cycle and lead) within 800 km +0.0285; the verifying agent counted every GFS track and got 24.5% flagged (13,723 kept), gain +0.0286. Share of the gain by threshold: 55% at gale vs below gale (+0.0149), 36% storm vs gale (+0.0099), 9% HF vs storm (+0.0025). HF Brier skill +0.318 to +0.345.
- ERA5-fed comparison, 2021-22 only (3,473 fixes, one season, illustration): ERA5 patch +0.0201, GFS patch +0.0259 on the same fixes. The GFS-fed gain is not smaller than the ERA5-fed one. A plausible, untested reason: GFS-fed state and environment inputs are noisier, leaving more for the pattern to add. GFS and ERA5 anomaly patches agree closely (correlation 0.97 to 0.98; mean bias -2.5 m in Z500, +0.4 hPa in SLP).
- Descriptive placebo (permuted GFS patches, frozen coefficients): -0.078, i.e. the frozen terms hurt badly when the pattern does not belong to the storm.

## Limits
ERA5 proxy truth. Five GFS seasons, all in the block that calibrated the 71.7 kt cut; at most 126 distinct season resamples, so the interval and p describe season-to-season variation here, not a future season. GFS versions changed over 2021-26. The patch is centred on the GFS tracker fix. The lead 24 fix is the forecast state at the valid time, outcome is the next 24 h in ERA5, so this tests forecast fields as inputs, not a 48 h forecast. Matched storms only (storms the proxy also tracked). Not a validated product; a GFE use needs a tracker.
