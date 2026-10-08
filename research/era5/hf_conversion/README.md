# What separates deepening storms that become HF from those that do not?

**ERA5 proxy, pipeline A**, seasons 2004-05 to 2025-26, Atlantic and Pacific separately, transitioning tropical cyclones in. Pre-registered in `PREREGISTRATION.md` (commit 4c312ee) before any converter and non-converter were compared. Converter = fix with gust index < 71.7 kt and normalised 24 h deepening rate ndr24 >= 1 whose gust index reaches 71.7 kt in the next 24 h; non-converter = same, but it does not. Matched 1:1 within basin x month on ndr24 (within 0.25).

## Answer in plain words

The converters sit in a stronger dynamical setting, but the structure of the near-storm environment at 0.25 degree is mostly not what distinguishes them.
- **Stage A (environment table at 1.5 degree, 788 pairs: 374 Atlantic, 414 Pacific):** 19 of 48 tests pass BH q<0.05 (17 of 48 when matched on latitude instead). The largest: **surface heat flux** about 0.4 to 0.5 SD higher for converters at t0 and 12 h before in both basins and 24 h before in the Atlantic (Atlantic +54 W/m2 at every lag, Pacific +44 to +46 at t0 and 12 h), **500 hPa vorticity advection** (+0.47 SD Atlantic at t0, +0.41 Pacific), **300 hPa divergence** (+0.2 to +0.4 SD), **Eady growth rate** (+0.29 SD Atlantic, still there 12 h earlier), column water vapour Atlantic at 12 and 24 h before, jet at t0 Atlantic. The SST gradient and SST-minus-T500 instability do not separate (MDE about 0.16 to 0.22 SD at t0: a well-powered null at t0). The 24 h change in each variable does not separate (0 of 16).
- **Stage B (0.25 degree, 250 pairs per basin):** 1 of 30 tests passes: Pacific 850 hPa baroclinicity at t0 (+0.32 SD). Not separating, with MDE about 0.2 SD at t0: distance to the 250 hPa jet maximum, where the jet sits relative to the motion, 850 hPa theta-e, near-surface stability (Atlantic stability +0.26 SD, q 0.06, borderline). At 12 and 24 h before t0 power is lower (MDE 0.3 to 0.8 SD): can't tell.
- **Big caveat for stage A:** converters are already deeper and older at t0 than their matched non-converters (central pressure -6.7 hPa Atlantic, -6.2 Pacific, about 0.8 SD; storm age +4.9 h and +6.8 h), because matching is on the deepening rate, basin and month only. Heat flux, vorticity advection and divergence are all things a deeper, more mature storm produces; stage A does not show that they cause conversion. The stage B result (that the structure at 0.25 degree adds little) is consistent with the storm itself being where the difference lives.

## Numbers
Population: 3,316 deepening fixes on 2,274 storms (Atlantic 1,475 / 1,023, Pacific 1,841 / 1,251); converting 1,287 fixes on 966 storms; one fix per storm per group; 149 tracks appear in both. Pair counts: 788 matched (178 converters unmatched), mean |ndr24 difference| 0.029. Stage B pull: 1,576 times, 33.8 GB (sized at 33.9 GB before the pull, N = 250 pairs per basin). Files: `results/stage_a_tests.csv`, `stage_a_change.csv`, `stage_a_latmatched.csv`, `stage_a_confounds.csv`, `stage_b_tests.csv`, figures `results/figs/conv_stageA.png`, `conv_stageB.png`, `conv_maps_*.png` (250 hPa wind, 850 hPa theta-e, |grad T850| rotated to motion; pixel stipple descriptive).

## Limits
- Hindsight: the label uses the next 24 h; a difference at t0 is an association that could be a forecast signal because the variable is known at t0, not a forecast skill claim (no model fitted, no held-out look).
- The 800 km owned gust is not recomputed in stage B (no gust pulled); the 25 km position check stands in.
- Sample 250 pairs per basin at 0.25 degree; lags -12 and -24 h have fewer pairs (lag fixes missing for young storms).
- ERA5 proxy, never validated before the archive. Counts near q = 0.05 move with the resample stream: "19 of 48" and "1 of 30" are good to about plus or minus 2 (deviation 7).

## Reproduce
`stage_a.py REPO` (no pull), `select_conv.py REPO 250`, `extract_conv.py 8` (`HFVS_SUB=hf_conversion_pull`, uncommitted), `hf_vs_storm/export_followup_tables.py` (`results/conv_values.csv.gz`), `stage_b.py REPO WORK`, `hf_vs_storm/figs_followup.py`.

## Verification (fresh Sonnet agent; `verify/`, `VERIFICATION.md`)
All 154 rows compared. Counts, group means, differences, SD, coverage and pair properties match exactly; p, MDE and four q flips differ within resampling noise (see `VERIFICATION.md`). Not checked: the pair draw, extraction, rotation, 25 km and 150 km position rules, maps, figures.
