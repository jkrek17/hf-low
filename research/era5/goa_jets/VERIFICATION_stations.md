# Verification record, coastal-station test (2026-10-17)

A fresh Sonnet agent, not shown the analysis scripts, recomputed from `results/stations_pairs.csv`, `hrrr_features.csv` and `stage0_zone_gust.csv`:
- H1 per station (n at 35 kt, median HRRR and ERA5 ratios, both Spearman values) and the count of stations where HRRR is closer to 1: all matched (4 of 9).
- H2 class counts and hits per 100 for HRRR-only, both and neither in the three zones: all matched (494/9.1, 191/14.1, 2390/0.9; 190/11.6, 167/27.5, 2744/0.5; 467/21.0, 199/50.8, 2463/1.9). It noted the Z2 HRRR-only class holds about 22 hits and the Z1 neither class about 21, so the ratios are noisy.
- Five station values (Amatuli 2021-02-12 00 UTC 51.1 kt, Augustine 2023-02-13 12 UTC 52.9, Cape Decision 2025-03-19 12 UTC 45.1, Cape Spencer 2020-01-18 12 UTC 42.0, Flat Island 2023-11-17 12 UTC 43.0) against the raw NDBC historical files: matched to 0.1 kt.
Not checked by it: the HRRR and ERA5 columns against their sources, the zone-to-station mapping, bootstrap intervals, minimum detectable ratios, and the post hoc file `results/stations_posthoc.txt`.
A reviewer (hf-reviewer) separately re-ran the registered bootstrap and obtained identical intervals; that is a re-run of the same logic, not an independent recomputation.
