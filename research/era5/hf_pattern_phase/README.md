# Does the hemispheric pattern change the structure of the HF lows that form?

**ERA5 proxy, pipeline A**, seasons 2004-05 to 2025-26, Atlantic and Pacific separately, transitioning tropical cyclones in. Pre-registered in `PREREGISTRATION.md` (commit 4c312ee) before any structure scalar was compared across index values. Index: the PR 41 leave-one-season-out pattern index (days -10 to -4 lead), standardised per basin, tercile of the week of each storm's first fix.

## Answer in plain words

Mostly no, with one Atlantic hint. Top-tercile weeks produce more HF lows (PR 41/85 already showed that); this asks whether they differ in structure once formed.
- **Pacific: no detectable difference** in any of 11 structure measures (MDE about 0.2 SD, median 0.24).
- **Atlantic:** HF onsets in top-tercile weeks have a larger storm-force wind area (+117,000 km2, +23%, 0.30 SD, q 0.011), a larger outer 48-kt radius (+66 km, q 0.011) and more of that area right of the track (+0.065, q 0.015). Central pressure is 2.1 hPa lower (0.17 SD) but not significant (q 0.21).
- **Fragile:** when strata are matched on latitude (S1) none of the 22 tests passes (Atlantic area +77,000 km2, q 0.13); at the HF peak anchor (S4, 161 vs 71 storms) none passes but the signs agree. A continuous per-SD slope (S2) passes for the same three Atlantic measures plus central pressure (-1.05 hPa per SD). So it's a lead worth carrying, not a finding.
- Gust and 10 m wind speed do not differ (the Atlantic peak gust is 0.3 kt lower in the top tercile; MDE 0.11 SD), so more index means a wider storm-force area, not stronger peak wind.

## Numbers
1,832 in-scope HF storms (Atlantic 1,007, Pacific 825; 168 out of scope with first fix May to September). Tercile counts: Atlantic top 432, middle 317, bottom 258; Pacific 353, 249, 223. Primary 22 tests: **3 pass BH q<0.05** (all Atlantic). S1 0 of 22; S2 4 of 22; S3 (Atlantic, no Greenland/Iceland maximum gusts) 3 of 11, same three; S4 0 of 22. Files: `results/pp_primary.csv`, `pp_S1..S4.csv`; maps (gust, MSLP, dewpoint, rotated to motion, pixelwise BH): `results/figs/pp_maps_*.png`, `pp_forest.png`. Low-pressure difference in the maps is northeast of the centre (about 4000 pixels pass), descriptive.

## Limits
- Held-out look: further look at the PR 41 index across all 22 seasons (no refit), logged in `research/era5/looks/hf_pattern_phase.log`, number provisional (25th at 2015-25 by the agenda count).
- The index and latitude/season move together; the latitude-matched variant removes the Atlantic result. Storm age at the onset is not matched.
- ERA5 under-resolves extremes; areas are biased low. Proxy only; no claim before 2004.

## Reproduce
`select_storms.py REPO` (tercile assignment, `times_pp.csv`), `hf_vs_storm/extract_storm.py` with `HFVS_SUB=hf_pattern_phase HFVS_TIMES=.../times_pp.csv` (15.3 GB pull), `hf_vs_storm/export_followup_tables.py` (`results/pp_values.csv.gz`), `analyse_pp.py WORK`, `maps_pp.py WORK`, `hf_vs_storm/figs_followup.py`.

## Verification (fresh Sonnet agent; `verify/`, `VERIFICATION.md`)
99 rows compared (primary, S1 to S4): counts, means, differences, SD, n and q pass/fail match on all rows; bootstrap p within tolerance on 69, off by more than 0.02 on 30 (resampling noise; differences are identical). Not checked: the extraction, tercile cut points, maps and figures, the resample stream.
