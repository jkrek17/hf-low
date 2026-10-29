# GEFS + ECMWF ENS blend of 10 m wind chances: pre-registration

Registered 2026-10-29, before any 2025-26 field is pulled or any blend score is computed. Asked for by Jason ("go", 2026-10-29) as the follow-up to `research/ens_wind/` (PR 285). Sustained 10 m wind from the models, in knots; no ERA5.

## Question
Does a weighted blend of the calibrated GEFS chance and the calibrated ECMWF ENS chance beat the calibrated GEFS chance alone, on weeks the blend was not fitted on? If not, report the null and stop (no site layer is proposed).

## Data (frozen)
- Winter 2025-26 only (the only winter with ENS): 00 UTC runs every 5th day from 15 Oct 2025 to 15 Apr 2026 (37 dates), leads 24/48/72/96 h, GEFS 31 members and ECMWF ENS 51 members, same cell maxima (1 degree, ocean points) as `ens_fetch.py`. Re-pulled because the PR 285 maps were kept in scratch, not saved. Sized by probe before the pull: GEFS 51 MB and ENS 86 MB per lead, so 148 leads each is about 20 GB plus about 1 GB ASCAT, under the 50 GB gate. AIFS ENS is dropped (PR 285: 0.17 to 0.33 worse, never reaches 64 kt).
- Truth as in PR 285: 48 kt = ASCAT cell max (RSS v02.1, cells with at least half their ocean sub-cells seen within +-3 h; no COARE, stated with every number); 64 kt = OPC hurricane-force low fix within 324 nm (600 km) in +-6 h, neighbourhood chance. "OPC agreed" wording, nautical miles.
- A run-lead is used only if both models posted it (paired cases).
- 16 cells: 2 basins x 2 thresholds x 4 leads.

## Forecasts compared
- **GEFS alone (baseline):** the committed all-winter isotonic table (`results/gefs_forecaster_tables.txt` source, `cal_table` in `gefs_primary.json`), applied to the GEFS member count. That table was fitted on all six winters including 2025-26, so it has seen the test weeks. This favours the baseline, which makes the test harder for the blend. Stated with the result.
- **Blend:** p = w * p_GEFS + (1 - w) * p_ENS. p_GEFS is the same committed table. p_ENS is a logistic calibration of the ENS fraction fitted on training weeks only: logit(y) = a + b * logit(clip(k/51, 0.5/51, 1 - 0.5/51)) over k >= 1 (damped Newton, as `stats.slope`), with k = 0 taking its own training rate. Two parameters, not isotonic, because there is one winter. w is chosen on training weeks only from 0, 0.1, ..., 1 by lowest Brier score (this is "weight by the 2025-26 scores"; one w per cell).
- No other blend form is tried before the primary is read.

## Held-out design
Leave-one-ISO-week-out over the run dates (the weeks of PR 285's `compare.py`): for each ISO week, fit a, b, k = 0 rate and w on the other weeks, score the held-out week. Scores are summed over held-out weeks.

## Primary test (one test, no multiplicity)
D = mean over the 16 cells of [BSS(blend) - BSS(GEFS alone)], BSS against the cell's base rate over the 2025-26 cases. Interval: 2,000 resamples of ISO weeks, jointly across cells (seed 20261029); the resample re-uses the held-out predictions (does not refit). **Rule: the blend "beats GEFS alone" if the lower end of the 95 % interval of D is above 0.** Otherwise the result is "cannot tell" or "worse", reported as a null and stopped.

## Secondary (labelled secondary; BH over the 16 per-cell differences at q = 0.05)
1. Per-cell BSS difference with the same week bootstrap.
2. Fixed w = 0.5 (no fitted weight).
3. Stricter baseline: GEFS recalibrated with the same leave-one-week-out logistic fit as the ENS term (so both are fitted on 2025-26 only); blend minus this.
4. Fitted w by cell, to say which model carries which lead and threshold.
5. Reliability of the blend by chance band (0-5, 5-20, 20-50, 50-80, 80+ %), pooled by threshold.
Cells with fewer than 30 events: reported as too few, excluded from D (count stated).

## Known weaknesses stated up front
- One winter, 37 run dates, every 5th day: wide intervals, and weights and the ENS calibration are rough. Plan: refit in spring 2027 and after winter 2026-27 (the ENS archive then holds two winters).
- PR 285 found ENS about equal to GEFS (ENS minus GEFS BSS between -0.07 and +0.005, intervals crossing 0). A blend gain would come from averaging two roughly equal forecasts with different errors, not from either being better.
- ASCAT sees about half of an area and saturates near 50 kt; 64 kt truth is an analyst-decided neighbourhood event.
- Weights and the ENS calibration are fitted inside the same winter that tests them, hence leave-one-week-out; adjacent run dates (5 days) share storms, so neighbouring weeks are not fully independent. The week bootstrap does not remove that.
- Marine wind project (GFE blend) numbers are reported there, not rechecked here.

## Deviations
(none yet; any change after a score is read is logged here as post hoc)
