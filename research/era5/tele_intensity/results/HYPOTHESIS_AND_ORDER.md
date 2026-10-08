# Hypothesis, planned test, and the order things happened

This note was written **after** the main test had run. The project's coordinator asked for the hypothesis and planned test to be on record before results; for this question that request arrived once the fit was already done. What follows is an honest timeline, so a reader can judge how far the analysis was fixed in advance. Times are UTC, 2026-10-08.

## Hypothesis

H1: NAO (Atlantic) and PNA (Pacific) change the intensity of hurricane-force lows, beyond their known effect on how often they occur.

H0, the alternative made explicit before the test: the apparent intensity effect is the index's own pressure signature and the shift in where storms sit, so it vanishes when depth is measured against the background and the index is taken before the storm exists.

## Timeline

1. **About 06:20. Exploratory look, not a test.** This used pipeline A catalog events and the same-time 5-day index ending at peak gust, with plain OLS and a season bootstrap. It found NAO+ and PNA+ storms about 2 hPa deeper, no change in deepening rate, and a large NAO northward shift. A lagged index halved the depth effect. These numbers are in `/mnt/project-files/science-questions/teleconnection-intensity-questions.md`, last written 06:26:54, together with the confounds they raised.
2. **06:26:54. Plan written down in that file, question 1.** The plan named these choices:
   - an index from days -10 to -4 before genesis
   - depth against the climatological background and against a 1000 km ring
   - deepening rate, peak gust and HF duration as outcomes
   - latitude and longitude controls, to measure how much of the effect is the storm moving
   - season-block inference
   - pressure outcomes on 1979+ and gust outcomes on 2004+
3. **06:27. Coordinator's defaults.** The coordinator added: pipeline A, the lagged index, the climatological anomaly, the same season windows, leave-one-season-out ranges, and season-clustered or permutation inference.
4. **06:28 to 06:33. Background fields built** (`background.py`). Climatology and 1000 km ring pressure for every event. No outcome had been fitted against an index yet.
5. **About 06:35. Smoke test** of `analyse.py`, using 20 bootstraps to check the code ran. Its point estimates were visible. **Nothing in the model, outcomes, windows, or predictors changed after it.** The only differences from the plan in step 2, all fixed before the smoke test:
   - NAO and PNA are fitted together in each basin, not one at a time.
   - Month fixed effects and a linear season term were added, as in the earlier frequency test.
   - Benjamini-Hochberg runs over the 28 primary rows.
6. **06:34 to 06:42. Full run**, 1,000 bootstraps and permutations: `summary.txt`, `results.csv`.
7. **06:43. Commit** `6a4069f`.
8. **After the commit:**
   - A fresh agent recomputed every quoted number from the committed files; all agree. It found that `event_table.csv` had been written at 4 significant figures, so `analyse.py` now writes it at 6. Results are unchanged.
   - The Atlantic NAO ring-depth permutation p is 1/1001 in this run and 2/1001 in the verifier's. So it is quoted as p ≤ 0.002.

## What was tested, and what counts as primary

The primary results are 28 tests, BH-corrected:
- Indices: NAO and PNA, both in each basin.
- Basins: Atlantic and Pacific.
- Pressure outcomes on 1979+: `minp`, `anom_clim`, `depth_ring`, `bg`, `ndr_max`.
- Gust outcomes on 2004+: `g800max`, `hf_hours`.

Everything else in `summary.txt` is secondary or a contrast:
- the 2004+ pressure rows
- the position-controlled fits
- the same-time index

No outcome, window, lag or control was added or dropped after the results were seen.

## Results, kept in full including nulls

See `README.md` and `summary.txt`. In short:
- **Atlantic NAO:** depth against climatology shows nothing (-0.07 hPa per SD). The raw effect (-0.92) is the storm-position shift. Ring-relative depth is -0.72 hPa.
- **Pacific PNA:** the depth effect (-0.90) is all background (-1.03). Storm-relative depth shows nothing (+0.12).
- **Deepening rate:** nothing in either basin.
- **Gust and HF duration:** nothing in the Atlantic. In the Pacific they rise but do not pass FDR (gust q 0.13, duration q 0.29).

## Power of the null results

At 80% power and two-sided α = 0.05, the smallest effect each null could have detected is about 2.8 × its season-bootstrap SE (`summary.txt`, lagged rows). All values are per SD of the index:

| Null result | SE | Detectable effect |
|---|---|---|
| Atlantic NAO depth against climatology | 0.25 | ≈ 0.7 hPa |
| Pacific PNA storm-relative (ring) depth | 0.26 | ≈ 0.7 hPa |
| Deepening rate, Atlantic | 0.015 Bergeron | ≈ 0.04 Bergeron |
| Deepening rate, Pacific | 0.013 Bergeron | ≈ 0.04 Bergeron |
| Atlantic gust index, 2004+ | 0.24 kt | ≈ 0.65 kt |
| Atlantic HF duration, 2004+ | 0.46 h | ≈ 1.3 h |

For scale, the HF-low events in this sample have:
- a deepest-pressure SD of about 11 hPa against climatology
- a median fastest deepening of about 1.3 Bergerons

So the nulls exclude effects larger than about 6% of an SD in depth and about 3% of the median in deepening rate. They do not exclude smaller ones.
