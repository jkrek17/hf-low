# Ensemble 10 m wind chances: pre-registration

Registered 2026-10-28, before any field is pulled or any score is computed. Asked for by Jason (2026-10-28): "10m winds from the GFS ... do we calibrate ensemble probabilities?" Earlier outlook numbers in this repository use GFS gust, not sustained wind, and no ensemble probability has been calibrated. This study does that for sustained 10 m wind.

All winds are model 10 m sustained wind (sqrt(u^2+v^2)), knots. No ERA5 is involved, so the proxy caveat does not apply; the truth sources have their own limits (below).

## Question
When k of N ensemble members reach 48 kt (or 64 kt) sustained 10 m wind, how often does it happen? Is the raw member fraction reliable, and what table should a forecaster read instead?

## Data (frozen choices)
- GEFS 0.25 degree (pgrb2sp25), 31 members (c00, p01-p30), 00 UTC, leads 24, 48, 72, 96 h (days 1-4). Archive exists from 2020-10. Seasons 2020-21 to 2025-26.
- Run dates: every 5th calendar day counted from 15 Oct each season, to 15 Apr. (Spacing is for volume; adjacent days share storms.)
- ECMWF ENS (51 members) and AIFS ENS (51), same dates and leads, season 2025-26 only (archive starts Mar 2025 / Aug 2025). Comparison only; never pooled with GEFS.
- Domains: North Atlantic 30-67N, 262-10E; North Pacific 27-67N, 135-240E (pipeline A boxes). Ocean points only (GEFS LAND mask, plus ASCAT's own land flag).

## Forecast unit
1 degree cell (4x4 points of the 0.25 degree grid), valid time = run + lead. Member "reaches T" when the cell-max 10 m wind >= T. Chance = members reaching T / N. Evaluated for T = 48 and 64 kt.

## Truth (two, reported separately)
1. ASCAT (RSS v02.1, Metop-A/B/C), 12.5 km wind converted to 10 m neutral-equivalent -> 10 m stress-equivalent with COARE as in the Gulf Stream work (PR 271), cell-max within +-3 h of valid time. A cell is scored only if ASCAT covers >= 50 % of its area in that window (coverage is decided by orbit, not by wind). ASCAT saturates near 50 kt, so T = 48 is the top usable threshold; 64 kt cannot be scored against ASCAT and the 64 kt row uses truth 2. Rain-flagged cells dropped; the count is reported.
2. OPC "agreed", never "verified": the OPC archive's low fixes (HF = hurricane-force category; storm-force category). For T = 64: forecast = any member reaching 64 kt within 324 nm of the cell centre; truth = OPC HF low fix within 324 nm in the 12 h centred on valid time. For T = 48: same with the storm-or-stronger category. This is a neighbourhood event, a different quantity from truth 1, and is labelled that way. OPC warns on sustained wind, so this is the right yardstick for 64 kt, but it is analyst judgement.

## Primary analysis (GEFS, 6 seasons)
Per basin x threshold x lead (2 x 2 x 4 = 16 cells, family for BH-FDR at q = 0.05):
- Reliability table: bins 0, (0,10 %], (10,30 %], (30,60 %], (60,100 %]; observed rate with a season-block bootstrap interval (2,000 resamples, seed 20261028; seasons are the independent unit).
- Test H1: calibration slope of logit(observed) on logit(clipped raw fraction) is below 1 (members overconfident). Block-bootstrap one-sided p.
- Brier skill score of raw fraction and of calibrated chance against the season-pooled climatology; calibration (isotonic, monotone) is fit leave-one-season-out and scored on the held-out season. Reported gain = calibrated minus raw BSS with interval. "Cannot tell" when the interval spans 0.
- Cases with fewer than 30 events in a cell are reported as "too few to calibrate", not fitted.
- Output the plain table: "when about 1 in 5 / 1 in 2 members show 48 kt, it happened about N in 10."

## Secondary (labelled as such)
- Season 2025-26 only: GEFS, ENS, AIFS ENS on the identical cases; paired BSS difference by 1-week blocks of run dates. Many cells will be too thin at 64 kt; stated as such.
- Lead dependence: does calibration shrink with lead (slope vs lead)?

## What counts as a result
Nulls are full results. No threshold, bin edge or neighbourhood radius is tuned after looking. Any change is logged below as a dated post hoc deviation.

## Known weaknesses stated up front
- ASCAT undercounts above ~50 kt and sees about half of any area; "did not happen" really means "ASCAT did not see it".
- 64 kt GEFS events at 0.25 degree are rare; the 64 kt table may be too thin in the Pacific.
- GEFS model versions: v12 throughout; any change inside 2020-26 is listed if found.
- Run dates every 5th day: fewer cases, so wider intervals; adjacent days not sampled.
- Marine wind project (GFE blend) numbers (about 15 % of ASCAT 48 kt+ peaks at about 4 % false alarm; 64 kt nothing skilful) are reported there, not rechecked here.

## Deviations
(none yet)
