# Life-cycle time-lapse: when do HF lows first look different from storm-force-only lows?

**ERA5 proxy, pipeline A** (`hf_history`), seasons 2004-05 to 2025-26, Atlantic and Pacific separately, transitioning tropical cyclones in. Pre-registered in `PREREGISTRATION.md` (commit 0686ff7) before any number was computed. Hindsight composite: groups are defined by the outcome at the peak, so this says when the two groups diverge, not forecast skill.

## Answer in plain words

Yes, and earlier than the wind shows. Relative to storm-force-only lows (peak gust index 54 to 71.7 kt), HF lows (peak >= 71.7 kt) have a steeper pressure gradient first and deeper pressure later.
- **Pressure gradient** (0-500 km): separates from about 48 h before the peak in the Atlantic (+0.45 hPa/100 km, 0.4 SD, q 0.01) and 36 h in the Pacific (+0.33, 0.4 SD, q 0.04, borderline), growing to +1.8 / +1.9 hPa/100 km (1.3 / 1.5 SD) at the peak.
- **Central pressure**: not different at -48 and -36 h; separates at -24 h in the Atlantic (-7.4 hPa, 0.6 SD) and at -12 h in the Pacific (-9.4 hPa, 0.8 SD), reaching -15 hPa at the peak. Pacific -24 h: -2.5 hPa, q 0.12, can't tell (MDE 0.41 SD).
- **Wind-field shape**: the maximum-gust radius is smaller from -24 h; the 48-kt area radius is larger only from -12 h (Atlantic, Pacific) onward. So the early lead is in the pressure field, not the wind.
- Track-table version (stage A, all storms, 00/12 UTC): 24 of 28 tests pass; central pressure separates from -36 h Atlantic and -24 h Pacific; the 12 h pressure change (deepening rate) already differs at -48 h in the Atlantic (99% of resamples).

## Numbers (stage B, 0.25 degree, 250 HF + 250 SF per basin; `results/stage_b_tests.csv`)
56 pre-registered tests, **39 pass BH q<0.05**. First persistent lag (earliest lag from which q<0.05 holds through 0 h; `stage_b_first_lag.json`): central pressure Atlantic -24 h (75% of resamples), Pacific -12 h (59%, -24 h in 40%); pressure gradient Atlantic -48 h (67%), Pacific -36 h (44%; -24 h in 36%, and -24 h in the verifier's resamples); 48-kt radius -12 h both; max-gust radius -24 h both. Pixel maps (gust, MSLP at -48, -24, 0 h): `results/figs/lifecycle_maps_*.png`, descriptive.

## Limits
- **Survivorship at early lags.** A storm only has a lag if its pipeline A track already existed: at -48 h 78 of 250 Atlantic HF and 48 of 250 SF storms (50 and 31 Pacific) qualify, so -48 h and -36 h rest on older storms; power there is low (MDE 0.4 to 0.8 SD). Coverage table: `results/stage_b_coverage.csv`.
- Groups differ in age at the peak (SF peaks earlier in life) and latitude; not corrected. Test 1 (`hf_conversion`) asks the forward-looking version.
- ERA5 under-resolves extremes; areas are biased low. ERA5 is a proxy, never validated before the archive begins.
- Held-out looks: none scored, no model fitted.

## Reproduce
`stage_a.py REPO` (track tables); `select_lc.py REPO`, `extract_lc.py 8` (16.7 GB pull, uncommitted caches in `$ERA5_WORK/hf_lifecycle`), `hf_vs_storm/export_followup_tables.py` (writes `results/lc_values.csv.gz`), `stage_b.py REPO WORK`, `maps_lc.py REPO WORK`, `hf_vs_storm/figs_followup.py REPO`.

## Verification (fresh Sonnet agent, from the pre-registration and committed tables only; `verify/`, `VERIFICATION.md`)
Recomputed and matched exactly: all counts, group means, differences, coverage, g800 reproduction. Bootstrap p within tolerance on most rows; q pass/fail agrees on 83 of 84 tests (Pacific pressure gradient at -36 h flips: q 0.040 vs 0.053-0.059). Not checked: extraction, heading and rotation, maps, bootstrap resample distributions of the first lag, CI columns.
