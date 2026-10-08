# Status

The ledger for this repository: what is in flight, what exists, what does not, and what may not be relied on yet. Read it at the start of a session and update it at the end, in the same pull request as the work. When this file and git disagree, git is right and this file gets fixed.

Last updated: 2026-10-08, after the gust-drift findings and the reconciliation of `c9dc994` were recomputed from committed data (`c55e74d`, merged into the integration branch as `8c95de7`). Earlier the same day this repository (`jkrek17/hf-low`) was split from `jkrek17/awips-tools`. Commit hashes below are this repository's. Pull request numbers refer to `jkrek17/awips-tools` unless they say otherwise.

## Threads

| Branch | Head | Was | What it holds |
|---|---|---|---|
| `main` | | | The archive CSVs, the site, the build and publish tools, the site tests, and the session rules. No research code. |
| `claude/exciting-fermat-8vcvgq` | `cfc6475` | integration branch | **All of the research:** precursor recovery, ERA5 pipelines A and B, per-cyclone P(HF), the gust-drift finding (`c9dc994`), its recomputation (`drift_check.py`, `c55e74d`, hf-low PR 8) and the fixed-depth count test (`drift_counts.py`, plan `fe0d0cc`, result `69bc06b`, hf-low PR 18), the near-storm intensity framework (`research/era5/intensity/`, hf-low PRs 7 and 12), pipeline A's full track population, the HF life cycle (`hf_history/lifecycle.py`, hf-low PR 16), teleconnections and HF-low intensity (`research/era5/tele_intensity/`, hf-low PR 11, merged as `d23546b`), collision review. hf-low PR 13 proposes merging this branch into `main` (decision 3). Has `main` merged in, so sessions here load the rules. |
| `claude/hf-lows-qc-mode` | `0800e90` | new during the split | Everything on the integration branch, plus a QC mode for the page that writes corrections to the spreadsheet. |
| `claude/era5-hf-history` | `b964cd1` | PR 80, with PR 82 merged in | Superseded: fully contained in the integration branch. |
| `claude/era5-hf-probability` | `71013e5` | PR 82 | Superseded: fully contained in the integration branch. |
| `claude/hf-lows-cleanup` | `b6e8f52` | PR 81 | Superseded: fully contained in the integration branch. |

**Work on the integration branch** unless the task is the QC mode. The three superseded branches predate `CLAUDE.md`, hold nothing the integration branch lacks, and can be deleted once Jason is satisfied; do not start work on them.

The merges on 2026-10-08 (`2b9af76`, `0b479db`, `0800e90`) had no conflicts. Afterwards the build check, the four node suites, and the Python tests passed on both live branches, and the 25 QC tests passed on the QC-mode branch.

## Decisions waiting on Jason

1. **How far back a gust-based record holds.** The numbers are now in the repository: `research/era5/drift_check.py` and `drift_check-result.txt` on the integration branch (`c55e74d`). What they show:
   - The index *value* drifts before 2001 in both pipelines. At fixed depth (955-975 hPa) over 1979-2000, B's 500 km gust rises +0.645 kt/decade (t = +2.40) and A's own 800 km gust rises +0.665 (t = +2.44). Both are flat from 2004.
   - Whether thresholded *counts* drift is **not settled either way**. The walk-back rested on pipeline A's per-basin trends over 1979-1996 (t = +0.82 and -0.02). Those tests had 12-15% power against the bias the drift implies, so they cannot show that counts are clean. Over all 22 pre-archive seasons, A's counts rise +5.1 per decade (t = +1.98). The drift predicts +5.7. Deep-storm counts rise by a similar amount (B tracks below 960 hPa: +5.8, t = +1.96).
   - With depth held fixed at track level, the drift estimates are +0.1 to +0.6 kt/decade, each with an se of 0.2 to 0.5. That bounds the count bias at about +1 to +6 events per decade over 1979-2000, the same size as the possible real rise.
   - So neither `c9dc994`'s "cannot be carried back" nor the walk-back's "counts are fine" is supported. A trend in A that spans 2001 may be anywhere from a fifth to all artefact.

   **The fixed-depth count test on the full track population** (`research/era5/drift_counts-result.txt`, hf-low PR 18; plan committed first in `fe0d0cc`):
   - **Pre-registered result (not settled, rule (c)).** No ramp in depth-adjusted counts is detectable within 1979-2000: +2.77 events per decade, 95% CI -3.0 to +8.6. The logistic version agrees (+3.6, z = +1.57). Power is 53% at +6 per decade and 17% at +3, so a drift inside the old band can be neither found nor excluded.
   - **Post hoc result (labelled as such, not replicated).** The *level* differs. At the same ERA5 depth, 1979-2000 has 7.8% fewer gust-based events per season than the 2004-2025 rates give (Welch t = -2.77), and both basins show it. Over 1979-2025, pipeline A's observed events rise +2.39 per decade (t = +2.35), while the depth-expected count is flat (-0.19). The whole trend sits in that offset, in the direction `c9dc994`'s per-fix drift predicts. A real change in gust at fixed depth would look the same, and the buoy and ship comparison is the test that can tell them apart.

   - **Buoy test** (`research/era5/buoy_drift/`, hf-low PR 17, 21 moored buoys from NOAA ISD, October-March 1979-2004; ERA5 is a proxy). It does not clear or confirm the gust field. ERA5 gust against buoy wind at 15 m/s and above has no trend over 1979-2000: -1.80 %/decade (se 1.92, t = -0.94, leave-one-season-out -3.19 to -1.21), against an index drift of about +1.1 %/decade, and the test has under 10 % power against a drift that size. The level moves the wrong way for the post hoc offset: ERA5 gust relative to the buoys is 2.9 % lower in 2001-04 than over 1979-2000 (se 1.3, t = -2.31, but only 4 seasons after, and ISD has no moored buoys after 2004). Over 1985-2000 ERA5 falls against all 11 stations with enough seasons (wind -6.83 %/decade, t = -6.59; gust -8.10, t = -6.31), a window chosen after seeing the data. Whether that fall is ERA5 or buoy hardware needs NDBC's hull history, which the work containers cannot reach. A general upward drift in ERA5's surface wind does not explain the gust-index drift; its cause (the storms, or ERA5's storm cores away from the buoys) is still open.
   - **Buoy thread's view: it agrees with the recommended answer below.** The buoys add a second, opposite-sign era dependence in 1985-2000, so levels and trends that cross 2001 stay out and within-era use stays in.

   **Recommended answer (Claude, pending Jason):** use pipeline A's gust-based counts before 2001 for variation *within* that era: rankings, interannual correlations, and teleconnection tests with an era term. Do not use them for levels or trends that cross 2001. Report no 1979-2025 trend in hurricane-force-equivalent counts as a climate signal unless the buoy and ship comparison clears the gust field. The depth-only count is the cross-era statement, and it is flat. This is option (a) of the earlier choice, extended to allow within-era use.
2. **Two ERA5 pipelines.** A and B were built in parallel and differ in almost every definition (table below). Keep both with distinct names and purposes, or retire one?
3. **Research into `main`.** None of it is there yet. The integration branch now holds all of it and could go in as one pull request. Before decision 1 is settled, or after?
4. **Validation floor year.** Files say the proxy "cannot be validated before 2001" in some places and "before 2004" in others. The pipeline B session's position is that `RECORD_START = 2004` must gate every window, because the archive was still starting before then. Confirm 2004 as the statement of record?
5. **QC mode and the live sheet.** `claude/hf-lows-qc-mode` adds a path that writes to the spreadsheet. Should a work session ever use it against the real sheet, or only against the mocks in `tests/qc/`?
6. **Cutover** from `awips-tools` (section at the end).

## Assumptions (Coordinator, pending Jason's review)

Working defaults that sessions follow until Jason confirms or changes them. Jason has said the sessions lead the research and he steers (2026-10-08).

- **Ledger-only pull requests into `main`** (changing nothing but `STATUS.md`) are routine. The thread that opens one merges it once the numbers are verified. Research code going into `main` (decision 3) and anything that deploys still wait for Jason.
- **Working default for decision 1**, until it is settled:
  - Trend and count claims use pipeline A's gust-based counts from 2001-02 on.
  - Counts for 1979-2000 may be shown and used within that era, labelled a proxy. Not for levels or trends across 2001: at fixed depth they run 7.8% below the post-2004 rate (post hoc, t = -2.77), and a ramp of up to +8.6 events per decade cannot be excluded.
  - Pressure-depth counts are the cross-check before 2001.
  - Fitting and testing stay at 2004-05 and later.
  - The buoy test (`research/era5/buoy_drift/`, hf-low PR 17) did not settle decision 1. It found no upward drift of ERA5's surface wind against buoys before 2001, but it cannot resolve a drift as small as the gust-index drift. This default stands.
- **Life-cycle definitions** (thread "HF low climatology gaps", 2026-10-08):
  - An HF fix is an in-domain fix at or above 71.7 kt, the rule that makes a track an event.
  - Genesis counts as observed only when the first fix is at 1000 hPa or above and north of 21N.
  - Onset counts as censored when the first in-domain fix is already HF, or when an out-of-domain HF fix precedes it.
  - Events within 400 km of an IBTrACS point at the same time are reported apart from the main numbers.
  - The main numbers use 2004-05 on.

- **Intensity framework working choices** (thread "Cyclone phase-space intensity framework"; details in `research/era5/intensity/README.md`):
  - Population: every pipeline A low below 1010 hPa in domain at 00/12 UTC, not only catalog events.
  - A class needs the track to survive 24 h, so rapid decay is under-counted. |NDR| > 3 is dropped as a tracker relink.
  - Transitioning tropical cyclones are left in.
  - L2 logistic with C = 1, never tuned. Outcomes are ERA5's own (perfect prognosis).
- **Teleconnection-intensity working choices** (thread "Teleconnection and intensity questions"; details in `research/era5/tele_intensity/README.md`):
  - Pipeline A catalog events only, so the results describe intensity given HF. The all-cyclone version is question 2, queued for the full track population.
  - Index averaged over days -10..-4 before the track's first fix. The same-time index is reported only as a contrast.
  - Background pressure has two forms: the calendar-month MSLP climatology, from every 4th WeatherBench2 chunk over 1979-2022 at 1.5°, and the 900-1100 km ring mean at the deepest fix.
  - Pressure outcomes use 1979+. Gust outcomes use 2004-05 on.
  - Transitioning tropical cyclones are included.

## What exists

### Archive and site (`main`)

- Two basin CSVs under `data/hf_lows/`, exported from the spreadsheet.
- `python3 tools/build_hf_lows.py --check` on `main` reports 1,932 lows, 8,073 fixes, 25 seasons, 3 dropped rows.
- Site tests pass on `main` and on the integration branch (checked 2026-10-08 after the split): composite 39, playback 40, regress 50, teleconnect 37, plus the Python teleconnections test.
- `flat/` regenerates from `docs/` with `python3 tools/publish.py --no-fetch --deploy flat --flat --yes` and matches the committed copy apart from build stamps.
- On `claude/hf-lows-qc-mode` only: opened with `?qc`, the page lists 194 data-quality flags and lets a fix be edited. Edits post to a new `doPost` in `web/HFArchiveExport/Code.gs`, which needs its own `QC_TOKEN`, checks the cell still holds what the page showed, writes it, and logs before and after to a "QC log" tab. Tests with mocks are in `tests/qc/`.

### Precursor recovery (integration branch)

- Fetch, parse, and track stages with tests. `data/hf_lows/precursors.csv` is committed: 6,658 rows (2,416 high confidence, 4,242 medium), covering 1,423 distinct events.
- Sidecars committed: `recovered_chain_hours.csv`, `archive_position_suspects.csv`, `collision_pairs.csv`, with `tools/review_collisions.py`.
- The build attaches precursors as a separate series: 696 of 1,928 hurricane-force events have a usable window.
- Validation numbers that may be quoted, and their limits, are in the `tools/track_hsf.py` docstring. Current statement: coverage 32% (16% at high confidence only); wrong-storm rate 3.5% on usable fixes; recall 50% on the fastest deepeners against about 72% for the rest.

### ERA5 pipeline A: `research/era5/hf_history/` (integration branch)

- A threshold on a gust index, calibrated on seasons 2021-22 to 2025-26, applied 1979 to 2025. Threshold 71.7 kt.
- Skill in `results/skill.txt`: AUC 0.990, POD 0.77, FAR 0.23, CSI 0.63, HSS 0.76, bias 0.99; leave-one-season-out threshold 71.4 to 72.3 kt.
- Catalog committed: 4,157 events and 4,154 matched null cases over 47 seasons.
- Per-cyclone P(HF) (`probability.py`, `results/probability.txt`): a logistic fit on the same gust index, P = 0.5 at 73.6 kt; a Pacific term is reported alongside as marginal.
- Tested against the gust drift (`research/era5/drift_check-result.txt`, section 6): its 800 km gust at fixed MSLP 955-975 hPa rises +0.665 kt/decade over 1979-2000 (t = +2.44), flat from 2004 (t = +0.61). Density at the threshold, 8.82 events per kt per season, sets how far a drift moves the count. The README now carries the caveat.
- Reproduction streams about 370 GB and needs go-ahead.
- Life cycle (`lifecycle.py`, `results/lifecycle.txt`, hf-low PR 16). It reads committed files only.
  - Sample: seasons 2004-05 on, leaving out the 114 events linked to tropical cyclones. That leaves 1,886 events.
  - Genesis to the first fix at or above 71.7 kt: median 36 h, quartiles 24/48, over the 1,152 events whose genesis is observed.
  - Time at HF: median 12 h, quartiles 6/24.
  - Onset comes before minimum pressure in 71% of events (Atlantic 66%, Pacific 77%).
  - ERA5 onset is within 12 h of the archive's first HF fix for 88% of the 1,244 matched events.
  - Every number was recomputed independently.
  - The method note, `results/lifecycle-plan.md`, was written after the run and says so. Read the result as descriptive, not pre-registered.

### Near-storm intensity framework: `research/era5/intensity/` (integration branch)

- Hart phase space (B, -V_T lower, -V_T upper) plus nine environment predictors at every pipeline A low fix at 00/12 UTC, on a 1.5 degree grid. Logistic models give P(24 h intensity class: rapid decay, decay, steady, deepening, rapid deepening, in Bergerons) and P(pipeline A gust index reaches 71.7 kt within 24 h and 48 h). ERA5 proxy throughout.
- Fitted and tested on 2004-05 to 2025-26 only (22 seasons, 159,430 fixes), leave-one-season-out against basin-month climatology. Full model: class RPSS 0.305 (storm state alone 0.242); rapid deepening HSS 0.47 (state 0.34, state plus Hart 0.45); HF within 24 h BSS 0.423, POD 0.60, FAR 0.40, CSI 0.43, HSS 0.59, bias 1.00; HF onset within 24 h BSS 0.291 (state 0.221); rapid decay weak, HSS 0.18. Trained 2004-14 and tested 2015-25: RPSS 0.309, BSS 0.426. All quoted values recomputed independently (hf-low PR 12 comment).
- Committed: code, `results/skill.txt`, `coefficients.csv`, `model.json`, phase diagrams, and the 2004-2025 input tables that reproduce the numbers. Pre-registration and post-hoc changes are logged in its README ("Analysis history").

### Teleconnections and HF-low intensity: `research/era5/tele_intensity/` (integration branch)

- Question: do NAO and PNA change how intense pipeline A events get, once the index's own pressure signature and the storm-position shift are taken out? ERA5 proxy; intensity given the HF-equivalent threshold. The question list it comes from is `/mnt/project-files/science-questions/teleconnection-intensity-questions.md` (question 1).
- Per SD of the index over days -10..-4 before genesis, with NAO and PNA fitted together. Month fixed effects and a season trend are included. p comes from 1,000 season-block permutations, and q is Benjamini-Hochberg over 28 primary tests.
  - Atlantic, NAO:
    - Central pressure is -0.92 hPa (q 0.017), but only -0.07 against the month climatology at the deepest fix. NAO+ moves storms +1.46° lat and +1.90° lon.
    - Against the 900-1100 km ring, storms are -0.72 hPa deeper (p ≤ 0.002, q 0.009).
    - Deepening rate, gust index and HF hours (2004+) show nothing.
  - Pacific, PNA:
    - Central pressure is -0.90 hPa. All of it is background (-1.03 hPa, q 0.009); against the ring the storm is +0.12 hPa.
    - Gust is +0.52 kt (q 0.13) and HF hours +0.87 h (q 0.29), both 2004+. This is a lead, not a finding.
  - Neither basin shows a change in deepening rate.
- Nulls exclude effects above about 0.7 hPa (depth), 0.04 Bergeron and 0.65 kt per SD. Season ICC of every outcome is -0.01 to 0.05.
- All quoted numbers were recomputed by a fresh agent from the committed files. `results/HYPOTHESIS_AND_ORDER.md` records the order of the exploratory look, the plan and the run; it was written after the run and says so.

### Frequency split: more cyclones or a larger HF share? (`research/era5/freq_split/`, hf-low PR 14, merged into the integration branch as `46097c8`)

Pipeline A, ERA5 proxy. Plan committed first (`933e68f`); results `95809d1`. Oct-Apr 2004-05..2025-26, NAO (Atlantic) and PNA (Pacific) lagged to days -10..-4 before genesis, per SD.
- Atlantic NAO: RR(HF) 1.096 = RR(all cyclones) 0.974 x RR(share) 1.125; the share carries it (f = 1.29, 95% bootstrap 1.07-1.78).
- Pacific PNA: RR(HF) 1.072 = 1.028 x 1.043; not resolved (f = 0.61, -0.65 to 0.79). With a count-matched pressure cut the Pacific share is 1.111 (f = 0.79, 0.60-0.89) and holds in 1979-2000.
- The same-time index gives 1.210 and 1.186 (the published effects); the lag removes more than half of the log effect.
- ONI, Pacific NAO: null. Atlantic PNA on the share: 0.88 (exploratory, 1 of 52 distinct tests).
- Verified by a fresh agent: the headline decomposition, odds ratios, same-time and depth values, the two main permutation p-values and the Atlantic leave-one-season-out range. Not independently checked: bootstrap intervals, the Jun-May, chain, other-index and joint rows, power figures and FDR q-values (listed in the README).

### ERA5 pipeline B: `event_fields.py`, `criterion.py`, `series.py` (integration branch)

- A per-moment probability from four features, fitted on early (2004-05), late (2021-25), and full (2004 to 2025) windows.
- The ocean-masked extraction is complete for all 47 seasons (`c9dc994`: 51,372 steps, none failed, about 190 GB). `hf_probability.csv.gz` on the integration branch holds 174,700 rows for 1979 to 2025.
- Skill at a count-matched cut, leave-one-season-out over 22 seasons (`c9dc994`): POD 0.727, FAR 0.274, HSS 0.705; gust alone gives HSS 0.659.
- Masking transitioning tropical cyclones (IBTrACS within 400 km) does not change skill.
- A basin term removes the pooled Atlantic and Pacific bias at no cost to skill. It is not in the committed series.
- `criterion-result.txt` holds the transfer test between windows; `series-result.txt` carries the stationarity flags.
- The copy of `hf_probability.csv.gz` on the superseded `claude/hf-lows-cleanup` is an older 13-season, pre-mask snapshot.

### ERA5 against moored buoys: `research/era5/buoy_drift/` (hf-low PR 17, merged into the integration branch as `cfc6475`)

- ERA5 10 m wind and instantaneous 10 m gust (the pipelines' gust variable) at 21 moored buoys (NDBC and Canadian; from NOAA ISD on AWS), 00/06/12/18 UTC, October to March, 1979-2004. About 31 GB streamed by reading two blosc blocks per Zarr chunk. Results in `results/buoy_drift-result.txt`; per station-season values in `results/station_seasons.csv`.
- ISD holds a knots unit error for eight buoys in 1992-2001 (speeds too low by 1.94). The script drops those station-months (27,189 of 183,560 pairs).
- Numbers and their limits are under decision 1 and in the module README. Choices fixed before the results and ones made after are in `results/analysis_history.md` (written after the run, not a pre-registration). A fresh agent recomputed the trend numbers quoted here from the committed files and the raw downloads and matched all of them. It found the pairs-dropped figure depends on an unstated minimum of 10 valid pairs per station-month; the result file now says so. Not independently checked: the ISD station list and the 5 m height assumption. The buoys are assimilated by ERA5, so drift away from buoys could be larger than drift at them.

### How A and B differ

| | A (`hf_history`) | B (`event_fields`) |
|---|---|---|
| Output | Event or not, by threshold; catalog with null cases | Probability per 6-hourly cyclone moment |
| Gust feature | Maximum within 800 km, points nearest this low only | Maxima within 300 and 500 km, plus area and gradient terms |
| Low detection | Below 1010 hPa | Below 1005 hPa and 4 hPa deeper than surroundings; scored at 990 hPa or below |
| Atlantic domain | 30-67N, 98W-10E | 30-70N, 80W-10E |
| Pacific domain | 27-67N, 135E-120W | 30-67N, 160E-120W |
| Months | June to May | September to May |
| Positive label | Track matched to an archive event, any fix category | Moment within 400 km of an archive HF fix of class "low" |
| Calibration | 2021-22 to 2025-26 | 2004-05, 2021-25, and 2004 to 2025 |
| Seasons committed | 47 | 47 |
| Ocean mask | From the start | Since `c9dc994` |

## What does not exist

- Any research code or data on `main`.
- The basin term in the committed pipeline B series.
- Validation of either ERA5 pipeline before the archive begins. There is nothing to validate against.
- Repair of the archive position errors the tracker found. They are listed, not fixed.
- Tests for `tools/review_collisions.py` or anything under `research/`.
- Moored-buoy records after 2004, and NDBC's per-station hull and anemometer history. Neither is reachable from the work containers (NDBC, MEDS, ICOADS and Copernicus hosts are refused by the network policy; ISD on AWS has no moored buoys after 2004).

## Gates in force

Do not build on these without closing the gate or stating the dependence.

- **Backfilled statistics across seasons.** The build reports "COVERAGE TREND: Bf not safe to compare across seasons". Headline coverage must come from recovered rows alone.
- **Recovered subset.** Not a random sample of storms; a rate from it describes it, not the archive.
- **Low-confidence precursor tier.** Refused.
- **ERA5 record is a proxy.** Label it so everywhere. No validation before the archive.
- **Gust index values before 2001.** Recomputed in `drift_check-result.txt`. At fixed storm depth over 1979-2000, the gust rises in both pipelines: B +0.645 kt/decade (t = +2.40), A +0.665 (t = +2.44). It is flat over the archive period. The index value is not comparable across eras.
- **Gust-based event counts before 2001.** In pipeline A at fixed depth there is no detectable ramp within 1979-2000 (+2.8 per decade, 95% CI -3.0 to +8.6, pre-registered). The level is 7.8% below the post-2004 rate (post hoc, t = -2.77), and that offset carries the whole 1979-2025 trend. Use these counts within the era only, not for levels or trends across 2001. Do not quote the subject line of `c9dc994` or the walk-back as a finding. The depth-expected count is flat over 1979-2025 (-0.19 per decade).
- **ERA5 surface wind against buoys, 1985-2000.** ERA5 falls against 11 moored buoys at -6.83 %/decade in winds of 15 m/s and above (t = -6.59; `research/era5/buoy_drift`, section 9, a window chosen after seeing the data). The fall may be in ERA5 or in buoy hardware. Any ERA5 wind-speed quantity compared across 1985-2000 should state this dependence.
- **Pipeline B series against the archive.** `series-result.txt` flags a trend in the series minus the archive over the 22 overlap seasons (+30.4 per decade, t = 2.63 against a critical 2.09). The commit notes that Mann-Kendall does not confirm it and that part of the divergence is the archive's own labelling. Do not read the series sum as a trend estimate.
- **Pipeline B `p_full`.** In-sample for seasons from 2004.
- **Pipeline B at high latitude.** Over-predicts by about 55% at 60 to 71N and in the most land-affected quartile (`c9dc994`).
- **Windows before 2004.** `RECORD_START = 2004`. The 2001-2006 transfer row in pipeline A's `skill.txt` (bias 2.05, HSS 0.44) reflects the archive starting up, not a failure of the threshold: the 2001-02 season holds 1 archive event against 85 ERA5 events in pipeline A's own table. Do not use seasons before 2004-05 to fit or to test.

## Reported by a session, now confirmed in the repository

The `c9dc994` session (https://claude.ai/code/session_01H34U5Bp9jBSLYgoGUR4VD8) reported three findings on 2026-10-08. All three were recomputed from committed files in `research/era5/drift_check-result.txt` (`c55e74d`). A fresh agent independently recomputed seven of the quoted numbers.

- **The archive's fixes per event drift.** HF fixes per event fall -0.181 per decade (t = -2.54) from the payload, -0.180 (t = -2.51) from the raw CSVs keyed by basin and ID, and -0.217 (t = -2.96) for class "low" only. Events and total HF fixes are flat. From 2006 the trend weakens to t = -2.12. **Confirmed.** It explains why pipeline B's per-moment comparison with the archive diverges over the overlap (+30.4 per decade) while A's per-event comparison does not (-2.2, t = -0.54 from 2004).
- **`RECORD_START = 2004`.** Seasons 2001-02 and 2002-03 hold 1 and 22 archive events against 85 and 113 in pipeline A. **Confirmed.** Already a gate.
- **ERA5 gust rises at fixed storm depth before 2001.** +0.645 kt per decade (t = +2.40). **Confirmed**, and it holds for pipeline A's index too.
- **The claim that thresholded counts are nevertheless era-comparable is not confirmed.** The reconciliation's own test reproduces (Atlantic t = +0.82, Pacific t = -0.02 over 1979-1996), but it had 12-15% power. See decision 1.

The same session's handoff document (https://claude.ai/code/artifact/e6afbb4c-269b-4237-bdd2-cfe70937491d) was corrected on 2026-10-08. A new section near its top, "Correction, 8 October 2026", gives the three findings, the reconciliation with its numbers, and the drift check, each with steps to reproduce. Its 1979 start-year gate and its overlap-flag item are marked superseded. The rest of the document predates `c9dc994` and the move to this repository; where it disagrees with this ledger, the ledger wins.

## Statements in the repository that are stale or conflict

Fix these as the files are next touched, in the same commit.

- `research/era5/README.md` says there is no calibrated gust index, no extended record, and that pressure gradients run 11% stronger before 1979. Pipelines A and B exist, and `stationarity2` reversed the gradient result. Its script list omits both pipelines.
- Commit `d34dd3d` gives 1979 as a safe start year. Commit `c9dc994` says in its subject and body that a gust-based criterion cannot be carried back before 2001. Its author has since said that conclusion was wrong at event level. Commit messages cannot be edited, so this entry is the correction: read `c9dc994` for its measurements, not for its headline.
- `research/era5/hf_history/track.py` docstring names two output files; the code writes one.
- Validation floor: 2001 in `criterion-result.txt`, `research/era5/README.md`, and `hf_history/README.md`; 2004 in `series.py` and `series-result.txt`.
- Event totals: 1,928 in the tracker docstring and build message, 1,932 lows in the build count and later commits. These may count different things (hurricane-force events against all lows); nobody has written down which.
- Commit `c48a8a8` says precursors are attributed to 1,860 events; the committed file has 1,423 distinct event ids.
- Track counts: 8,138 in `hf_history/results/skill.txt`, 8,144 in `results/probability.txt`, although the later commit says the calibration reproduces exactly.
- High Seas cache size: about 70 MB in `.gitignore`, about 450 MB in commit `1c70427`.
- Left over from `awips-tools`: `tools/publish.py` and `.gitignore` mention `docs/cps/`, `docs/img/`, and other paths that are not in this repository; `.claude/skills/clasp/SKILL.md` describes a different Apps Script app; `docs/README.md` omits several files from its list; `flat/README.md` is a copy of `docs/README.md`; the publish manifest is still named `.awips-publish-manifest.json`; User-Agent strings still say `awips-tools`.

## Open work

- Decision 1: the fixed-depth count test ran (hf-low PR 18). The ramp is not detectable, and the post hoc level offset needs independent confirmation. The buoy comparison (PR 17) did not provide it: it has under 10 % power against a drift that size and found the level moving the other way over 2001-04 (4 seasons). A comparison with the power needed would need NDBC's hull history and its post-2004 records, which the work containers cannot reach; ship reports (ICOADS) were not tried because that host is unreachable too. A pre-registered replication of the level test would be a second one, for example on pipeline B tracks or with 2001-03 held out.
- Teleconnections and intensity: question 2 is answered (frequency split, entry above). Open from it: the Pacific share disagrees between the gust and the pressure-depth definitions of "strong"; the cyclone count could change with storms moving across the fixed domain edges (not checked). The Pacific PNA gust and HF-duration lead from question 1 needs an independent check.
- Intensity framework: skill by basin is not broken out; untested on operational (GFS) analyses; gale and storm thresholds deferred by Jason; the fitted model has not been applied to 1979-2003.
- Pipeline B: put the basin term in the committed series; look at the over-prediction at high latitude.
- Precursors: resolve the collision worklist (39 pairs: 31 sequential, 8 concurrent); recall on fast deepeners is not at parity and tuning was stopped deliberately; tropical and post-tropical systems are not parsed.
- The unexplained +2.6 residual in the explosive share (commit `a9be84c`).
- Hard-coded `/home/user/awips-tools/...` paths in `research/era5/` pilot scripts (`compare.py`, `envdemo.py`, `matchmonth.py`, `mismatch.py`, `sample.py`, `tip.py`). They will not run from a clone of this repository as written.

## Future tasks

Work nobody has started. A thread picking one up changes its status here in its ledger PR. Items marked "carried over" were already logged as open or deferred elsewhere and are copied here so Jason can prune them. Work already running (teleconnection question 1, high-latitude Atlantic, HF wind structure, ERA5 wind drift vs buoys, the decision-1 count-drift test) is not listed. Started 2026-10-08.

| # | Task | From | Data, and whether it needs a go-ahead | Status |
|---|---|---|---|---|
| 1 | **Test the near-storm intensity framework on GFS analyses**, to see whether its skill holds on operational data. The model (`research/era5/intensity/`, hf-low PR 12) is perfect-prognosis on ERA5, so its skill on a forecast model's analysis is unknown. | **Requested by Jason**, 2026-10-08 | GFS analyses; source, period and size not yet checked. Not ERA5, but size it first and treat a pull over 50 GB as needing Jason's go-ahead. Pre-register the comparison before scoring. | Not started |
| 2 | Break the near-storm framework's skill out by basin. | Carried over: `intensity/README.md` | Committed 2004-2025 input tables. No pull. | Not started |
| 3 | Gale and storm-force targets for the near-storm framework. | Carried over: Jason deferred these ("later") | Committed inputs; needs a gust-index threshold for each class. No pull. | Not started |
| 4 | Apply the fitted near-storm model to 1979-2003. | Carried over: Open work above | Pipeline A tracks and environment fields exist. Gust outcomes before 2001 carry the decision-1 caveat. No pull. | Not started |
| 5 | Seasonal cycle, and whether HF lows share the Pacific midwinter suppression. | Carried over: climatology agenda question 4 | Archive and pipeline A tracks. No pull. | Not started |
| 6 | Where HF lows form and intensify: ocean fronts and moisture. | Carried over: climatology agenda question 5 | Near-storm environment fields; a pull only if the non-HF population lacks them (size it then). | Not started |
| 7 | Deep-storm record since 1979 and HF trend since 2001. | Carried over: climatology agenda question 6 | Pipeline A tracks. No pull. Depends on decision 1 for anything gust-based before 2001. | Not started |
| 8 | Sea state under HF lows. | Carried over: climatology agenda question 7 | ERA5 wave fields, estimated under 10 GB for 2004-05 on (unverified; size with HEAD requests first). May use the buoy thread's data. | Not started |
| 9 | Whether HF lows come in clusters. | Carried over: climatology agenda question 8 | Archive and pipeline A catalog. No pull. | Not started |
| 10 | Share of HF lows that come from tropical cyclones. | Carried over: climatology agenda question 9 | Pipeline A catalog and IBTrACS (`research/era5/tc_candidates.csv`). No pull. Shares a TC definition with item 15. | Not started |
| 11 | Does the MJO's heating longitude change Pacific HF-low intensity 10 days later? One pre-registered test. | Carried over: science question 4 | CPC MJO pentads on the repo's frozen EOFs; pipeline A catalog. No pull. | Not started |
| 12 | Does ENSO change the seasonal intensity distribution (tail counts below 950 hPa)? Expected null. | Carried over: science question 5 | Monthly ONI; pipeline A catalog, depth outcomes only. No pull. | Not started |
| 13 | Do the archive and the ERA5 proxy agree on the intensity-teleconnection relation? | Carried over: science question 6 | Archive CSVs on `main`, 2004-05 on, with a season trend term for the fixes-per-event drift. No pull. | Not started |
| 14 | Does the environment (jet, Eady growth, SST gradient) carry the teleconnection signal? | Carried over: science question 7 | Near-storm thread's environment fields and fitted outputs (hf-low PR 12). No pull. | Not started |
| 15 | Do teleconnection states change how often transitioning tropical cyclones become HF lows? Low power. | Carried over: science question 8 | `research/era5/tc_candidates.csv`, pipeline A catalog. No pull. | Not started |
| 16 | A properly powered buoy drift test: does ERA5 wind drift against buoys before 2001, to help settle decision 1? | Carried over from the buoy drift thread | NDBC hull and anemometer metadata, plus moored-buoy records after 2004. NDBC, MEDS, ICOADS and CDS are blocked by this environment's network policy, so the network allowlist must be widened first. Size any pull before running. | Not started; blocked on network access |

Sources: `research/era5/intensity/README.md` (integration branch); the buoy drift thread (item 16); `/mnt/project-files/climatology-agenda/hf-climatology-agenda.md`; `/mnt/project-files/science-questions/teleconnection-intensity-questions.md` (question 3 was answered with question 1).

## Fixes wanted in the sheet

Corrections that belong in the spreadsheet, because a CSV-only fix is overwritten by the next fetch.

- 38 suspected position errors in `data/hf_lows/archive_position_suspects.csv` (integration branch), including 8 first hurricane-force fixes the tracker refuses to anchor on. The QC mode on `claude/hf-lows-qc-mode` is a route for making such fixes in the sheet once it is merged and deployed.

## Cutover from awips-tools

`awips-tools` still holds its copy of these files and still publishes `/hf-lows/` on the public site. Until the steps below are done, **`awips-tools` is the copy that publishes** and changes made only here do not reach the public site.

**Do all hurricane-force work here from now on.** Sessions were still committing to `awips-tools` while the split was under way. Three commits that landed there were carried over (`c9dc994`, `5f05669`, `747d87d`). Anything committed to `awips-tools` after 2026-10-08 05:01 UTC is not in this repository; check before relying on either copy.

1. Done 2026-10-08: the setup pull request (number 1 in this repository) is merged.
2. Here: turn on GitHub Pages by hand, Settings > Pages > Build and deployment > Source = "GitHub Actions", then re-run "Deploy archive site to Pages". Its first run failed at the step that tries to enable Pages automatically; the data rebuild, the page validation, and the flat-sync check all passed.
3. In `awips-tools`: remove `hf-lows` from `OWN_FOLDERS` and the archive build from `site_publish.yml`, so that workflow copies `/hf-lows/` forward as it does other projects' folders.
4. Here: set the secret `WEB_TOKEN` and the variable `PUBLIC_SITE_REPO`. `publish-site.yml` then publishes `/hf-lows/`. Do not do this before step 3, or two repositories write the same folder.
5. In `awips-tools`: remove the moved files and leave a pointer to this repository. Close pull requests 80, 81, and 82 there with a link to the branches here.
6. Production: the NOAA copy is made by hand with `tools/publish.py` from a checkout. Point that checkout at this repository.

One deliberate difference: `publish-site.yml` leaves `data/qc-report.txt` out of the public copy, as `tools/publish.py` does for production. The `awips-tools` workflow copied it.
