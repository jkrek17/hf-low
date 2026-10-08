# HF lows against storm-force-only lows: environment and storm-scale composites

Plan: [PREREGISTRATION.md](PREREGISTRATION.md) (committed first, `1d11f8c`; two deviations logged at its end). **ERA5 proxy, pipeline A** (`research/era5/hf_history`), seasons 2004-05 to 2025-26,
transitioning tropical cyclones in, Atlantic and Pacific separately. HF = gust index >= 71.7 kt (2,000 storms); SF = peak gust index in [54, 71.7) kt (8,136), matched on basin x month.
Published page: https://claude.ai/artifact/3PCaev43n6ATE2AoAKE8eL . Follow-up plan (not started): [PLAN_THREE_TESTS.md](PLAN_THREE_TESTS.md).

## Answer in plain words
- Before formation: HF onsets follow a slightly deeper trough and stronger 250 hPa jet than storm-force peaks. Atlantic box (53.5N, 41W): Z500 -7.2 m [-11.8, -2.2] and jet +1.12 m/s over days -10 to -4, -22.0 m and +1.94 m/s on the onset day. Pacific (42N, 171E): -5.8 m and +0.82 m/s, -9.1 m and +1.31 m/s. SST: no difference. The difference builds from about day -7 (Atlantic).
- At the storm: central MSLP 971 vs 984 hPa (Atlantic), 974 vs 985 (Pacific); pressure gradient within 500 km about 60% steeper; radius of the maximum gust 302 vs 407 km (Atlantic), 233 vs 333 km (Pacific); gust factor 1.51 vs 1.47 and 1.52 vs 1.48; wider 48-kt gust area at onset (+158 and +154 km), which shrinks to +42 km (Atlantic) and none (Pacific) at the same central pressure; moister air within 500 km (partly latitude: HF anchors are 3 to 6 degrees farther south). Share of the 48-kt area right of motion: no difference detected (about 0.2 SD detectable).
- 12 of 20 pre-registered environment tests and 23 of 28 structure tests pass BH q < 0.05. Gust and 10 m wind differences are true by construction (HF is defined by the gust index) and are outside the family.

## Files
`select_times.py` (storms and anchors, labels only), `size_pull.py`, `fields_large.py` (copy of `hemispheric/fields.py` plus water vapour, seasons 2004-05 on, 20 Aug to 31 May),
`extract_storm.py` (helpers copied from `hf_structure/extract.py`), `strat.py` (stratified difference, season-block bootstrap, BH), `analyse_large.py`, `analyse_storm.py`, `make_figs.py`, `build_page.py`,
`results/` (`large/boxtests.csv`, `series.csv`, `maps_*.npz`; `storm/scalars*.csv`, `anchors.csv`, `comp_*.npz`; `figs/`), `verify/`.

## Reproduce (about 44 GB streamed: 16.7 GB large scale, 27.2 GB storm scale; raw fields are not committed)
    python3 -I research/era5/hf_vs_storm/select_times.py . 400
    python3 -I research/era5/hf_vs_storm/fields_large.py $ERA5_WORK/large
    ERA5_WORK=$ERA5_WORK python3 -I research/era5/hf_vs_storm/extract_storm.py 8
    python3 -I research/era5/hf_vs_storm/analyse_large.py . $ERA5_WORK/large 2000
    ERA5_WORK=$ERA5_WORK python3 -I research/era5/hf_vs_storm/analyse_storm.py . 2000 1000

## Limits
Different storm age at the anchor (HF onset 34-37 h after track start, SF peak 29-33 h) and latitude; the SF group is large and ordinary; 54 kt is a choice; pixel maps are tested pixelwise and spatially correlated; 5.6 degree large-scale fields are coarse;
storms in the same week share an environment (season-block bootstrap reduces this); no heading for short-lived storms (10% HF, 12% SF); 2 m temperature not pulled; ERA5 reads low in extreme storms. No held-out scoring, no looks logged.

## Verification
A fresh Sonnet agent recomputed group sizes, anchor counts, the gust-index reproduction (max 0.1 kt), all 36 structure-measure means and differences for C1, and the Z500 and jet box-mean differences (within 0.09 m, 0.01 m/s; its storm counts were 2 to 14 lower from a different missing-day rule). **Not independently checked:** bootstrap p, intervals, q, MDE, C2, right-of-motion share, MSLP/SST/water-vapour tests, secondary families, maps, figures, the extraction, heading, and the random draw. See `verify/VERIFICATION.md`.
