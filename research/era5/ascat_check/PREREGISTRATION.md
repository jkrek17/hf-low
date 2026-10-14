# Does the ERA5 gust proxy agree with scatterometer winds? (ASCAT check, pre-registration)

Written 2026-10-13, before any scatterometer value was read or joined to any storm. **No scatterometer data had been
reachable when this was written** (the container's network policy blocks every scatterometer host), so nothing below
can have been tuned to outcomes. ERA5 **proxy**, **pipeline A** (`research/era5/hf_history`; transitioning tropical
cyclones are in). Requested by Jason as item 2 of the 2026-10-13 list ("the review panel's top criticism").

**Looked at before this file was committed:** column names and counts of `hf_history/results/all_tracks.csv.gz`
(35,456 tracks from 2004-05 on; 17,552 with a peak time at 00 or 12 UTC); the cache layout in `era5-cache/README.md`;
the definition of the index in `hf_history/extract.py`. No scatterometer wind, and no cache gust at a storm, was looked at.

## Question and the answer we will give

Does pipeline A's gust index (strongest ERA5 gust over ocean within 800 km of a low) agree with independent satellite
winds around the same lows, and does it sort lows into the right wind classes? Answer: **yes** (agrees well enough to
use as a class proxy), **no** (does not), or **can't tell**, with the ASCAT wind that corresponds to the 71.7 kt
index level.

## Design decisions that follow from scatterometer limits (fixed now)

1. **Sustained wind vs gust.** Scatterometers retrieve an equivalent-neutral 10 m wind averaged over a 12.5-25 km cell;
   the index is an instantaneous gust. The primary comparisons are therefore **rank-based or threshold-free**, not
   anchored to an assumed gust factor. The conversion is *estimated* (the ratio of ERA5 gust to scatterometer wind at
   matched quantiles), never assumed; for display only a nominal gust factor of 1.3 is shown, labelled nominal.
2. **Saturation and low bias at high wind.** Standard ASCAT winds saturate/under-read above roughly 50 kt and
   QuikSCAT Ku-band is rain-contaminated. So 71.7 kt-equivalent cannot be tested directly. Level tests use a
   registered ASCAT level of **45 kt** ("storm-force-equivalent"), and the class test is whether ERA5 index >= 71.7 kt
   tells apart footprints with ASCAT max >= 45 kt. Winds above 60 kt are reported but never used as truth.
3. **Rain / quality.** Cells with the product's own rain, ambiguity or quality flag raised are dropped; the number
   dropped is reported by basin. If the product has no rain flag, the QuikSCAT part is not run (ASCAT's C-band is
   far less rain-sensitive; it is kept with a caveat).
4. **Swath gaps and timing.** ERA5 cache holds 00 and 12 UTC only. A scatterometer cell counts if its time is within
   **+-1.5 h** of the cache time. ERA5 gust is bilinearly interpolated to each cell. Comparison is on the **same cells**:
   ERA5 is sampled only where the scatterometer observed, so swath gaps cannot lower the scatterometer maximum
   relative to ERA5.
5. **Sample.** Pipeline A tracks, seasons 2004-05 to 2025-26, whose peak-index fix time is 00 or 12 UTC, with
   scatterometer coverage of >= 40 % of the ocean cells inside the 800 km disc at that fix (coverage reported; the
   rule is applied before any outcome is computed). Unit of independence for intervals: **season block bootstrap**
   (22 seasons), and storms within season for the pairing tests. No pre-2004 season is fitted or tested (Decision 1).
   QuikSCAT (1999-2009) is used only for 2004-05 to 2009-10.
6. **Land and edges.** Ocean cells only (cache land-sea mask < 0.5); the index definition ignores cells with gust < 17 m/s.
   Atlantic fixes north of 60 N are a separate lower-confidence group.

## Pre-run check on the ERA5 side (outcome-free, can run now)

R0: recompute the index from the cache for every sample track (disc max over ocean cells, 800 km, ignoring
ownership by other lows). Pass if the recomputed value is >= the stored `gust800_kt` minus 1 kt and within +-1 kt of it
for >= 90 % of tracks. If not, stop and report; the comparison is then to the recomputed disc maximum, named as such.

## Hypotheses and tests (primary family, 5 tests; BH-FDR q < 0.05 across all five)

Footprint = scatterometer-observed ocean cells within 800 km of the fix, at the accepted pass.

- **P1 (agreement of maxima).** Spearman rho between ERA5 gust max on the footprint and scatterometer max, per basin
  (2 tests). Decision threshold: rho >= 0.6 with the lower 95 % bound >= 0.5 counts as agreement.
- **P2 (classes).** AUC of the ERA5 footprint max in separating footprints with scatterometer max >= 45 kt from the
  rest, per basin pooled (1 test). Agreement if AUC >= 0.85, lower bound >= 0.80.
- **P3 (the HF class).** Among footprints, scatterometer max in those with ERA5 footprint max >= 71.7 kt versus a
  depth-matched set (same basin, minimum pressure within 5 hPa, ERA5 max < 71.7 kt): paired difference in
  scatterometer max (1 test). The proxy is supported if the difference is > 0 with q < 0.05.
- **P4 (stability by sensor/era).** Ratio of scatterometer max to ERA5 footprint max in ERA5 >= 54 kt footprints,
  ASCAT-A/B/C and (if kept) QuikSCAT, tested for a common value across sensors (heterogeneity, 1 test). A drift
  here would reopen Decision 1's question for 2004 on; stated before seeing it.
- Not tested, descriptive only: the matched-quantile gust:wind conversion at 71.7 kt, pointwise slope, bias by
  latitude and by distance from the fix, coverage by basin. These are labelled descriptive; none enters a
  verdict.

**Verdict rule.** "Yes": P1 and P2 meet their thresholds in both basins, P3 passes FDR. "No": P3 fails with >= 80 %
power at +5 kt, or P2 lower bound < 0.7. Otherwise "can't tell". A well-powered null is a full result.

## Power

Planted-effect simulation (`test_synthetic.py`) on the real ERA5 side (storm sample and footprints from the cache),
with a synthetic scatterometer built from the ERA5 field (wind = gust / 1.3 + 3 kt noise, saturating near 50 kt,
random swath gaps), run before any real data. It reports the P1/P2/P3 detection rates at the actual sample
size. If the achieved n is far below the simulated one, power is recomputed and reported.

## Pull plan

Sources (any one that becomes reachable): NOAA CoastWatch/ERDDAP ASCAT swath or RSS daily 0.25 degree ASCAT grids
(with time of day per cell); PO.DAAC QuikSCAT. Only passes within +-1.5 h of 00/12 UTC on the sample dates are
needed (about 17,500 storm-times). Estimate (unverified): 2-10 GB for ASCAT, < 5 GB for QuikSCAT. Over 20 GB: stop
and ask Jason; over 50 GB is his gate. Nothing is re-pulled from ERA5. A compact storm-matched table
(`scat_footprints.csv.gz` plus cell file) goes to `/mnt/project-files/ascat-check/` with a README.

## Deviations (post hoc)

Logged in order, each tagged.

1. **Before any outcome.** R0 on a random 300 of the 17,552 sample tracks (seed 0; ERA5 side only): 77 % of tracks
   recover the stored index within 1 kt and 100 % are >= stored - 1 kt. The misses are two kinds: 20 % of tracks have
   a stored index of 0 (no ocean cell of that low's own territory at >= 17 m/s at the stored peak fix), and among
   tracks with a stored index > 0, 95.1 % recover within +-1 kt, the other 4.9 % are higher (a neighbouring low's cells in
   the disc, median +3.3 kt). Resolution: the sample is restricted to tracks with stored `gust800_kt` > 0 (R0 passes at
   95 % of those); the disc maximum ignoring ownership is what the footprint ERA5 value is, and is named so.
2. **Before any outcome.** P3 classes the storm by pipeline A's **full-disc index** (`gust800_kt` >= 71.7), not by the
   footprint ERA5 maximum, because the footprint maximum is partial by construction. Scatterometer response is still
   measured on the footprint.
3. **Before any outcome (2026-10-14, source fixed once scatterometer hosts became reachable).** Source is RSS ASCAT v02.1
   daily 0.25 degree bytemaps (`data.remss.com/ascat`), one file per satellite per UTC day, each cell carrying its own
   observation time (6 min steps), so the +-1.5 h rule is applied per cell. Metop-A (2007-03 to 2021-11), Metop-B
   (2012-12 on), Metop-C (2019-07 on). Metop-A uses the `bmaps_v02.1` folder (the reprocessed `v02.1.1` folder holds no files
   before 2016), so all three are on the same product version. Seasons 2004-05 to 2006-07 therefore have no scatterometer
   pass (the registered QuikSCAT part is **not run**: it was optional, "if kept", and ASCAT-A/B/C already give three sensors
   for P4); the sample is storms from 2007-03 on, and the season bootstrap resamples the seasons that have data. Inventory
   13,059 daily files, 20.0 GB; 11,974 were opened (days with no sample fix in a window were skipped). Only cells within
   800 km and +-1.5 h of a sample fix were kept (15.26 million cell rows).
4. **Before any outcome.** Quality rule for "the product's own rain flag": cells with the **scatterometer rain flag** (bit 0 of
   the product's rain byte) raised are dropped; the product carries no ambiguity flag. The radiometer-collocated rain code
   and the sum-of-squares residual are **not** used as filters (no registered threshold, and radiometer rain would remove the
   heaviest-rain storm sectors preferentially). A descriptive sensitivity with radiometer rain > 0 also dropped is reported
   and enters no verdict. Wind is the product's 10 m speed (0.2 m/s steps) x 1.9438 kt per m/s; codes above 250 (land, ice,
   no data) are not winds.
5. **Before any outcome.** The shared cache decoder needed `numcodecs`, which cannot be installed here (no PyPI). `blosc_nolib.py`
   decodes it with the system libzstd; R0 on the same 300 tracks reproduces `results/r0_recovery_300.csv` exactly (max
   absolute difference 0.0 kt). Looked at before running: file counts and sizes, one day of one satellite to check the layout
   (cell count only), nothing about wind values against ERA5.
6. **Before any outcome statistic (2026-10-14, bug fix).** The first collocation gave 0 of 12,326 footprints at >= 40 % coverage
   (maximum 0.24). Cause: the scatterometer grid has cell centres on x.125 degrees and the ERA5 grid on x.0, and the cell
   counting rounded half-integers to even, merging neighbouring cells, while the denominator counted ERA5 nodes.
   `collocate.py` now counts scatterometer-grid cells on both sides (denominator: scatterometer-grid ocean cells with centres
   within 800 km; ocean = the four surrounding ERA5 land-sea-mask nodes average < 0.5, symmetric about the half-cell offset).
   Only coverage and cell counts had been looked at; no wind value had been compared. The synthetic test originally placed
   its fake cells on ERA5 nodes and so could not reveal this; it now uses scatterometer-grid centres (1,637 footprints,
   1,542 accepted, mean cover 0.62, same recovery as before). The 40 % rule is unchanged. Result: 12,327 footprints, 5,791
   with coverage >= 40 %.
