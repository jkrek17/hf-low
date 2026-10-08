# Do archive HF fixes per event drift, or step at a scatterometer change? (RA-11)

Plan: [PREREGISTRATION.md](PREREGISTRATION.md) (committed first, `5e16d30`; addendum A `708f1b2`; deviations listed at its end).
Code: `analyse.py`, `power.py`, `late_onset.py`; results in `results/` (`fixes_per_event-result.txt`, `power.txt`, `late_onset-result.txt`, `tests.csv`).
Reproduce: `cd research/era5/fixes_per_event && python3 analyse.py && python3 power.py && python3 late_onset.py` (no pull; about 3 minutes).
Data: archive CSVs, 22 seasons 2004-05..2025-26 (decision 1: nothing earlier). ERA5 is a **proxy**, pipeline A (`research/era5/hf_history`), used only for the S7 contrast.

## Answer

**Yes, the decline is real in the archive, and it looks like one step near the end of QuikSCAT (Nov 2009), not a gradual drift.**
Fixes per event fall from about 3.27 (2004-09) to 2.94 (2010-25), a drop of 0.33 (about 10%). Events and total HF fixes per season are flat, so event counts are unaffected.

- Linear drift: -0.180 per decade (se 0.072, t -2.51; season-bootstrap 95% [-0.32, -0.07]); q 0.035 in the 5-test primary family.
- Step at QuikSCAT's end (2009-11-23, verified): -0.348 fixes per event [-0.48, -0.21], q 0.005. AIC prefers it to linear by 4.7, so by the pre-registered rule this is a **step at a sensor date**.
- Step at ASCAT-B launch (2012-09): -0.257, q 0.005, but AIC ties with linear (0.8), and it is not separable from the QuikSCAT step (the decline is one drop; AIC prefers the QuikSCAT date by 3.9). ASCAT-A (2006-10) and Metop-C (2018-11) steps: not detected (q 0.31, 0.28). The ASCAT-A test is uninformative: only 3 seasons precede it and its bootstrap rejects 17% under the null.
- Data-chosen best split (secondary, scan): seasons 2010 on, step -0.336, permutation p 0.016, q 0.044. Chosen on this data, so it supports the date region, not the exact date.
- Post hoc: inside 2004-09 and inside 2010-25 the slopes are flat (+0.05 and +0.00 per decade). With a QuikSCAT step in the model the remaining linear term is -0.017 per decade (p 0.87). That is, one drop and no continuing drift.
- Both basins show the QuikSCAT step (Atlantic -0.348 q 0.044, Pacific -0.348 q 0.042; both -0.348 to the third decimal by coincidence of rounding). Atlantic first HF fix south of 60N -0.43 (q 0.042); north of 60N -0.19 (p 0.52, only 61 events before; difference not detected, p 0.38).
- Shape (PR 79 link): the share of one-fix events rises +3.1 points per decade (q 0.050) and the share of events with 3 or more fixes falls 4.1 points per decade (q 0.044). Archive-side sustained-HF counts (k = 2, 3) therefore carry the step; PR 79's archive contrast should include an era term. The proxy side is unaffected (below).
- Not the storms (as far as the archive's pressures show): event minimum pressure trend -1.45 hPa/decade (p 0.15); adjusting for it, basin and month leaves the slope at -0.238 per decade (unadjusted -0.213, both p below 0.01 by season-bootstrap).
- ERA5 pipeline A proxy: fixes per HF track do not fall; they rise +0.22 per decade (t +2.46; archive minus proxy -0.40 [-0.56, -0.27], q 0.014). Different definition (gust at least 71.7 kt in domain), so this says the storms as the proxy sees them did not shorten, consistent with recording practice. It does not prove the sensor is the cause.

What it does not show: why. A step at QuikSCAT's end fits forecasters having less direct wind confirmation, but a coincident change in practice or in the spreadsheet is not excluded. ASCAT operational-use dates (not launch dates) were not verified; Metop-C's launch date is from memory.

## Power (planted-effect simulation, `results/power.txt`)

- Linear: slope MDE 0.25 per decade (80%); power for the observed -0.18 is 71%. The observed linear and QuikSCAT sizes sit at or just under the 80%-power minimum (QuikSCAT step MDE 0.40; power 0.91 at 0.40 and 0.57 at 0.25), so the sizes are likely somewhat inflated (winner's curse), not confirmed values.
- Step MDEs: ASCAT-B 0.30, Metop-C 0.25, ASCAT-A above 0.80 (uninformative).
- Under a planted 0.40 step, AIC picks the step form in 74-84% of runs and linear in 4-5%. Under a planted linear decline AIC picks linear over a QuikSCAT step in 58% and a step in 8%. So an AIC lean to the step is informative but not decisive; the flat within-segment slopes are post hoc support.
- Under no effect the step test rejects 7-8% for P3-P5 (slightly liberal) and 17% for P2. The reported P3 and P4 p-values (0.001, 0.002) are far below this, but steps near the middle of the record share the same decline, so the two significant steps are not independent confirmations.

## Duplicate archive IDs and data handling

Primary uses raw IDs (so -0.18 reproduces the earlier -0.181 and -0.180). Merging the 21 listed pairs (19 apply inside 2004-2025) leaves 1842 events; linear -0.16 per decade (-0.160 to -0.170 depending on which duplicate row is kept), QuikSCAT step -0.318, ASCAT-B -0.231; the conclusion does not change. Counting identical timestamps once: -0.174, QuikSCAT -0.334. The HF row `pac:2024202506` with stamp 20241101018 is read as 2024-11-01 18Z (belongs on the sheet-fix list, with the other invalid dates). Only Category `HF` counts, not `DHF`.

## Addendum A: late HF category north of 60N in the archive (separate, secondary)

PR 86 found the strongest wind coming at least 6 h after minimum pressure about 3 times as often for Atlantic storms with the gust maximum at or north of 60N (ERA5 proxy). The archive has no wind, so this tests the HF **category** instead, using rows already loaded:
- Archive HF onset is almost never late: 96% of Atlantic events (963) have the first HF fix at or before the lowest-pressure fix (mean 15 h before), 3.9% are 6 h or more after (proxy: 10.9%).
- North of 60N versus south: 4.5% vs 3.8%, OR 1.22 [0.43, 2.51], q 0.67; at 12 h or more 2.0% vs 1.1%, OR 2.01 [0.87, 5.50], q 0.29. Power for a planted OR of 3.1 at this base rate is 0.90, so the interval above (not reaching 3.1) means the archive does **not** show the proxy's pattern in the HF category. It is not a contradiction: the quantities differ (a forecaster's HF category against an ERA5 gust maximum), and the base rate is a quarter of the proxy's.

## Looks

None: nothing was fitted on one block and scored on another. No entry in any looks log. Pre-2001 untouched.

## Verification

A fresh Sonnet agent recomputed, from the CSVs and `lifecycle_events.csv` with its own code and without reading the scripts: event counts and mean; the linear slope, se, t, Durbin-Watson; all four step differences and their event splits; the 22 season means and segment means; the one-fix and 3-fix shares and slopes; the basin slopes and QuikSCAT steps; the north/south 60N groups and steps; the proxy mean and slope; the pressure mean and trend; the adjusted and unadjusted event-level slopes; the joint linear-plus-step fit; and the late-onset counts, shares and odds ratios. All match to the rounding shown. The merged-ID slope matched only under one duplicate-row rule (-0.160 keeping HF rows, -0.164 keeping the first; reported as a range).
**Not independently checked:** every bootstrap interval and p, permutation p, BH q, the AIC differences, the power simulation, the 60N difference interval, the proxy minus archive interval, and the Pacific/Atlantic q values. The sensor dates come from web sources (QuikSCAT end: CIMSS blog; Metop-A and -B launch months: NOAA OSPO); Metop-C is from memory.

## Post hoc and corrections

Listed in PREREGISTRATION.md: the S6 bootstrap relabelling error (fixed before anything was reported), the post hoc segment block, and the reading of the malformed date.
