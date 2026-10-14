# Independent recomputation of the ASCAT check numbers (2026-10-14)

Done by a fresh Sonnet agent given the footprint table, the raw cells, the ERA5 cache and the definitions, but not the analysis code or README. The statistics (part A) were recomputed from the analysis's own footprint table; only 40 footprints (part B) were rebuilt from the raw cells.
Scripts: `verify_part_a_stats.py` (statistics from `../results/scat_footprints.csv.gz`), `verify_part_b_rebuild.py`
(rebuilds 40 random footprints from the raw cells and ERA5 cache; imports `blosc_nolib` from this folder's parent for the cache decode).
Part B needs the shared-folder cells and ERA5 cache (not committed).

**Matched:** storm counts (4,252; Atlantic 2,352, Pacific 1,900; 430 with index >= 71.7, 265 and 165; 20 seasons; sensors 1,886 / 1,632 / 734);
Spearman 0.9460 and 0.9236; AUC 0.9713; P3 +12.003 kt (n 423); class table; matched quantiles 43.54 / 54.82 / 62.20 kt; per-sensor ratios 0.8312 /
0.8383 / 0.8243 (n 797 / 758 / 357); bootstrap intervals for both correlations and AUC (to 0.001). 40 of 40 rebuilt footprints match (cell count exactly, ASCAT max and ERA5 max to 0.01 kt).

**Difference:** the P3 interval. The analysis code fixes the depth matching once and resamples seasons of the matched differences: [11.1, 13.1]. The verifier redid the matching in each draw: [10.4, 13.8]. Both are reported; the conclusion does not change.

**Recomputing the stored index** from the cache gives a value within 1 kt of the stored one for 37 of 40 (the misses are higher, a neighbouring low's cells).

**Not independently checked:** pointwise cell-pair table (slope, bias by latitude and by distance), rain-flag loss shares, season-ratio series and slope, tercile and partial correlations, the
radiometer-rain sensitivity run, the synthetic test, the coverage fix's effect on counts other than the 40 rebuilt rows.
