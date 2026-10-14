# Verification of tools/storm_review_facts.csv

Recomputed by `verify_numbers.py` (own implementation, numpy/pandas only, run with `python3 -I`) from the pipeline A ERA5-proxy 6-hourly fixes CSV (138,481 fixes) and the 30-storm list. The ERA5 record is a proxy.

Run: `python3 -I verify_numbers.py FIXES_CSV FACTS_CSV [STORMS_CSV]`

## Matched (30 of 30 storms each)
t_minp, minp, lat_minp, lon_minp, t_gust, gust_idx, first_hf, last_hf, n_hf, first_fix, last_fix, p_first, dp24_max, t_dp24, lat_dp24, bomb_thr. Zero mismatches.

Also: the track set, basin and rank in the facts match the storm list (30/30); within each basin, recomputed minp rises with rank (atl ranks 1-15, pac ranks 1,4-11,13,14,16-19).

hf_lat/hf_lon were not defined in the brief. Tested candidates: the first HF fix matches 29/29 storms that have an HF fix (within 0.01 deg); last, max-gust, mean, median and min-pressure HF fix each match at most 4/29. So hf_lat/hf_lon = position of the first HF fix (inferred, not an instructed definition).

## Tolerances
Times and n_hf exact. msl, p_first, dp24_max, gust_idx within 0.05; lat/lon within 0.01 deg; lat_dp24 and bomb_thr within 1e-6. Facts values came out equal to float noise in practice. HF fix = g800 >= 71.7; g800 <= 0 treated as missing for the gust maximum. Ties for lowest msl: first fix used (none tied).

## Not checked
- sw_max (needs ERA5 fields).
- Columns in the storm list but not in the facts (anom_clim, depth_ring, tc, name_inferred, peak_*, n_dom, season, date_minp) and facts text columns name and gust_flag.
- The fixes CSV itself (its extraction from ERA5) and the Pacific rank gaps (ranks 2, 3, 12, 15 are absent from the list; the reason is not given in the inputs).

## Things that look odd (values match, so these are properties of the data, not errors)
- Track 330264 (atl rank 10): no HF fix (peak gust index 70.7 < 71.7), so first_hf, last_hf, hf_lat, hf_lon are empty and n_hf = 0.
- Track 388856 (atl rank 13): the track is only 24 h long (2020-02-14 18Z to 2020-02-15 18Z, 5 fixes), so dp24 exists at one fix only (the last) and it is a rise: dp24_max = -7.1 hPa. It also has just 5 fixes for the lat_mean.
- dp24_max below the bomb threshold (24 h fall smaller than the latitude-scaled bomb value): 100785 (13.1 vs 24.1), 388856 (-7.1 vs 23.4), 252423 (8.6 vs 20.5). These storms were not explosively deepening by this measure.
