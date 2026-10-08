# Pre-registration, test 2: HF lows forming under a high or low hemispheric-pattern index

Committed before any HF storm's structure was joined to the pattern index. **ERA5 proxy, pipeline A**, seasons 2004-05 to 2025-26, Atlantic and Pacific separately, transitioning tropical cyclones in.
Plan from `hf_vs_storm/PLAN_THREE_TESTS.md` (test 2).

Question (Jason's test 2): does the hemispheric pattern (PR 41/85 index) change the *structure* of the HF lows that form, or only how many storms convert to HF? Answer forms: yes (q < 0.05, size given), no (well-powered, MDE given), or can't tell.

What had been looked at: the 2,000 HF storms and their labels, the index table's layout and the structure scalars' names. No structure scalar was compared across index values, and no HF count was split by the index here (PR 41 and PR 85 did the count questions).

## Data
- Index: `hemispheric/results/oos_index_{atl,pac}.csv` column `idx`, the leave-one-season-out pattern index of PR 41 (fitted to archive weekly HF counts, not refitted), 30 weeks of 7 days from 1 October, 22 seasons. It already carries its lead (fields in days -10 to -4 before the week). Standardised per basin to mean 0, SD 1 over the 660 weeks. Per SD of this index.
- Storms: pipeline A HF (peak `g800` >= 71.7 kt) with first fix in October to April (the weeks the index covers); each assigned the index of the week containing its first fix (the PR 85 / `hem_channels` rule). HF storms with first fix May to September are out of scope (counted).
- Anchor: HF onset (first 6-hourly in-domain point at or above 71.7 kt), as in `hf_vs_storm`. Storm-scale structure from `hf_vs_storm/extract_storm.py` (same code, re-detection, ownership, 121 x 121 boxes, 25 km). The 800 HF onsets already extracted (400 per basin, random draw of `hf_vs_storm`) are used as they are; the rest of the in-scope onsets are pulled so the primary sample is **every in-scope HF storm** (expected about 1,860 of 2,000). Pull sized by HEAD requests before pulling; 5 fields x about 14.4 MB per time, about 17 GB.
- **Missing-data rules (fixed now):** a storm has no heading if fewer than two 00/12 UTC fixes exist (rotated composites drop it; scalars that do not need a heading keep it); a scalar missing for a storm drops that storm from that test only; `g48_right` needs a heading and a 48-kt area. Re-detected low must lie within 25 km of the catalog position, else dropped (counted).

## Groups and tests
- **Top tercile vs bottom tercile** of the weekly index (cut points are the 33rd and 67th percentiles of the 660 weekly values per basin), HF storms in those weeks; middle tercile shown, not tested. Stratified by calendar month (strata basin x month), season-block bootstrap (22 seasons, 2,000 resamples, seed 20261012), BH-FDR, MDE = 2.8 x SE. Counts per group are reported (the index predicts counts, so the top tercile holds more storms).
- **Primary family P (22 tests):** top minus bottom tercile difference of 11 structure scalars x 2 basins: central MSLP, MSLP gradient 0-500 km, maximum gust, maximum 10 m wind, 48-kt gust area, radius of maximum gust, outer radius of 48-kt gust, right-of-motion share of 48-kt area, 2 m dewpoint within 500 km, gust factor, area of gust >= 71.7 kt. (Gust and area are not by construction here: every storm is already HF.) BH over 22.
- **Secondary families (own BH):** S1, latitude-matched: strata = basin x 5 degree latitude band of the anchor instead of month (the index and latitude may both move with the season). S2, continuous: slope of each scalar on the index per SD, month fixed effects, season-block bootstrap (all in-scope storms, not terciles). S3, drop storms whose maximum gust is within 100 km of Greenland or Iceland (Atlantic). S4, anchor at HF peak instead of onset (uses stored peak fields).
- **Maps (descriptive):** top and bottom tercile composites (rotated, 121 x 121) of gust, wind speed, MSLP, dewpoint, their difference with pixelwise BH q < 0.05 stipple.
- Power: per scalar MDE in SD of the scalar. Expected group sizes about 400 to 650 (Atlantic) and 300 to 500 (Pacific); a difference smaller than 0.2 SD with q >= 0.05 is "no difference" only if MDE <= 0.2 SD, else "can't tell".

## Held-out looks
The index is out of sample for every season, and PR 41 scored 2015-16 to 2025-26 at least five times (`hemispheric/results/heldout_looks.log`). This is a further look at those seasons, not a new pattern search and no model is refitted; logged in `research/era5/looks/hf_pattern_phase.log`.

## Deviations (post hoc)

1. The g800 reproduction check is vacuous here: the extractor writes the catalog g800 into the table it compares with, and the verifier found g800 = g800_cat on all 1,832 onsets and 742 peaks (max difference 0.0 kt). It checks that the re-detected low is the catalog low (match <= 25 km held for all), not an independent gust recomputation.
2. The S2 table's n column counts all storms with an index value (1,007 Atlantic, 825 Pacific), not the 928 and 768 that have a right-of-motion share (needs a heading and a 48-kt area). Slopes are unaffected.
3. The 'bottom' column in the result tables is the top mean minus the stratified difference (reweighted onto the strata where both terciles are present), not the raw bottom-tercile mean.
4. The pixel maps use 1,000 resamples, not 2,000 (descriptive).
5. Held-out look numbering is provisional (see the looks log).
