# Pre-registration, test 3: life-cycle time-lapse (when do HF and storm-force-only lows first become distinguishable?)

Committed before any life-cycle number was computed and before any field was pulled for it. **ERA5 proxy, pipeline A**, seasons 2004-05 to 2025-26, Atlantic and Pacific separately,
transitioning tropical cyclones in. Plan from `hf_vs_storm/PLAN_THREE_TESTS.md` (test 3); groups, matching and bootstrap as in `hf_vs_storm/PREREGISTRATION.md`
(HF = track peak gust index >= 71.7 kt, 2,000 storms; SF = peak in [54, 71.7), 8,136; strata basin x calendar month of the peak; SF weighted to HF's strata; season-block bootstrap of 22 seasons,
2,000 resamples, seed 20261012; BH-FDR; MDE = 2.8 x SE).

Question (Jason's test 3): at what lead before peak do the two groups first look different, i.e. what lead time could a forecaster use?

**What this is not.** The groups are defined by the outcome at the peak, so a hindsight composite will always separate near the peak. The lead time found here is "when the groups diverge", not forecast skill: a storm
that looks like an HF storm early is not thereby an HF storm (test 1 asks that question, conditioning on the state at a lead). Stated before any result.

What had been looked at: group sizes and the existing 00/12 UTC table's columns. No life-cycle value was computed.

## Stage A: track tables only (no pull)

- Source: `intensity/results/fixes_2004.csv.gz` (00/12 UTC fixes of every pipeline A low: central MSLP `msl`, 12 h pressure change `dp12` (negative is deepening), gust index `g800`, `lat`, `lon`).
- Anchor: the track's 00/12 UTC fix with the highest `g800` (ties: the earlier one). Lags: -48, -36, -24, -12, 0, +12, +24 h on the same track (7 lags).
- **Missing-data rules (fixed now):** a lag is missing when the track has no fix at anchor + lag (before its start, after its end) or the field is NaN (`dp12` at a track's first fix). Storms are used at every lag where they have a value (primary, n given per lag).
  Sensitivity: the complete-lag subset (storms with all 7 lags). Tracks with fewer than two 00/12 fixes are dropped, with the count reported. Anchors whose season file does not exist cannot occur (all seasons 2004-2025).
- **Primary family (28 tests):** HF minus SF stratified difference of `msl` and `dp12`, 2 basins x 7 lags, BH over 28. Descriptive, outside the family: `g800` (the group definition is a gust-index peak), storm age, latitude, speed.
- **Headline:** the **first persistent lag**: the earliest lag, scanning from -48 h forward, at which the `msl` test (and, separately, the `dp12` test) has BH q < 0.05 and stays at q < 0.05 at every later lag through 0 h. Reported per basin with the share of 2,000 season bootstraps in which that first lag is -48, -36, -24, -12 or 0 or none. The bootstrap interval is the 2.5 and 97.5 percentile of the resampled first lag; the q used inside each resample is the BH of that resample's p from a fixed-sample permutation-free t-statistic of the stratified difference (se from an inner 200-resample), a rule fixed here; the point estimate uses the main 2,000-resample p.
- **Power and MDE:** per test, MDE80 = 2.8 x SE, in hPa and in storm-to-storm SD. A lag with q >= 0.05 and MDE above 0.25 SD is "can't tell", not "the same".
- **Confounds to report with numbers:** storm age at the anchor (SF peaks earlier in life: 24 h vs 36 h median), anchor latitude, owned-ocean fraction (stage B), tracks that start inside their lag window. Pressure at the anchor differs by construction (peak of a stronger storm); the early-lag question is the one that matters.

## Stage B: 0.25 degree composites at 12 h steps (pull)

- Sample: 250 HF and 250 SF per basin drawn at random (seed 20261012) from the 400 + 400 per basin already pulled for `hf_vs_storm`, so lag 0 reuses stored fields. Lags -48 to +24 h every 12 h around the **6-hourly peak time** (pipeline A `peak_time`), 7 lags. 6 new lags x 1,000 storms = about 6,000 times.
- Fields per time: instantaneous 10 m gust and MSLP only (5.4 MB per time). Position at a lag: linear interpolation of the 00/12 UTC track positions to the lag time; snapped to the nearest low re-detected with pipeline A's rules within 150 km; otherwise missing (counted, reported per lag and group). Rotation uses the heading at the anchor. Pipeline A's owned 800 km gust must reproduce the table's `g800` within 0.5 kt for >= 99% of fixes where the table has a value at 00/12 UTC.
- **Primary family (56 tests):** stratified HF minus SF difference of central MSLP (min within 100 km), MSLP gradient (ring 450-550 km minus centre, per 100 km), outer radius of the 48-kt owned gust area, radius of the maximum gust; 4 measures x 2 basins x 7 lags, BH over 56. Descriptive: maximum gust, 48-kt area, composites with per-pixel BH maps at lags -48, -24, 0.
- **Headline:** first persistent lag as in stage A for central MSLP and MSLP gradient (the pre-specified pair); the other two are secondary.
- **Pull size (to be sized by HEAD requests before pulling):** about 6,000 times x 5.4 MB, about 32 GB. The thread's total across tests 1 to 3 was stated to Jason as about 94 GB and he said go ahead; if any run turns out much larger than estimated, stop and say so.
- Dry run on a few times with no all-NaN output before the full run.

## Deviations (post hoc)

None yet.
