# gfs_raw_gust: raw GFS gusts against OPC hurricane-force designations

Pre-registration `PREREGISTRATION.md` (commit fcfc9eb, before any join of GFS gust to an archive label; one scoring run plus the p-value fix in Deviation 1). `score.py` writes `results/raw_gust.json` and `results/raw_gust_fixes.csv.gz` (per-fix G, model p and labels; 2 MB). No new pull: gust is `g800` in the committed GFS tracks of `gfs_test` (pipeline A detector on GFS PRMSL and GUST, a PROXY detector; the GFS fix is a forecast, not an observation). Truth is the OPC archive, Category HF only. Seasons 2021-22 to 2025-26 (Oct-Apr cycles, 00/12 UTC, leads 0-48 h), 90,131 fixes, 2,899 with an archive HF fix within +-3 h and 600 km (label Y0). The independent sample is 5 seasons; these seasons also calibrated the 71.7 kt ERA5 index, so nothing here is independent of that calibration.

Raw predictor G = maximum ocean gust within 800 km of the low in GFS `GUST:surface` (kt). No ERA5-trained model is involved.

## Results

| | Atlantic | Pacific |
|---|---|---|
| Nominal 64 kt, lead 0 (R1): HSS [95% CI] | 0.470 [0.392, 0.538] | 0.426 [0.367, 0.495] |
| ... POD / FAR / bias | 0.92 / 0.66 / 2.67 | 0.92 / 0.70 / 3.03 |
| Count-matched, held-out seasons, lead 0 (R2): HSS | 0.524 [0.488, 0.550] | 0.567 [0.540, 0.610] |
| ... POD / FAR / bias | 0.54 / 0.45 / 0.98 | 0.58 / 0.42 / 0.99 |
| Count-matched, lead 24: HSS | 0.492 [0.482, 0.507] | 0.542 [0.516, 0.568] |
| Nominal 64 kt, lead 24: HSS (bias) | 0.419 [0.344, 0.480] (2.96) | 0.404 [0.355, 0.450] (3.22) |
| Threshold set on other seasons (kt), 5 folds | 71.5 to 72.3 | 71.5 to 72.5 |
| ERA5 proxy at 71.7 kt, same truth, seasons, months (R3): HSS (POD / FAR) | 0.618 (0.67 / 0.40) | 0.634 (0.65 / 0.36) |
| GFS lead 0 matched minus ERA5 proxy: HSS difference [CI] | -0.094 [-0.111, -0.073] | -0.067 [-0.086, -0.048] |
| Frozen state model p vs raw G, AUC on Y24 at lead 24 (R4): model minus raw [CI] | -0.001 [-0.012, +0.006] | +0.014 [+0.011, +0.017] |
| AUC of raw G (Y0), leads 0 to 48 | 0.964 to 0.969 | 0.968 to 0.977 |

Decisions (12 registered tests, BH q): R1 and R2 pass the HSS >= 0.30 criterion in both basins at both leads (8 of 8, q = 0.0006). R4 passes in both basins (the model adds at most 0.017 AUC over raw G; 2 of 2, q = 0.0006). R3: Pacific passes the "within 0.10" criterion (q = 0.001), Atlantic does not (interval reaches -0.111; q = 0.19). 11 of 12 pass FDR at 0.05 in this reading; the HSS >= 0.30 bar inherited from A1.3 is lenient, so the pass count says more about the bar than about forecast quality.

In plain words:
- **Raw GFS gust near a low is usable on its own, but not at 64 kt.** At 64 kt it catches 92% of OPC HF fixes but about two thirds of its yes calls are false (bias near 3). Setting the number so yes-calls equal OPC's count gives bias 1, POD 0.54 to 0.58, FAR 0.42 to 0.45, HSS 0.52 to 0.57.
- **The matched number is the ERA5 one.** Matched on other seasons, the GFS threshold is 71.5 to 72.5 kt in every fold, i.e. the same 71.7 kt used for the ERA5 proxy. GFS and ERA5 `g800` are alike on average (GFS minus ERA5 -0.5 kt, from `gfs_test`).
- **Lead barely matters out to 48 h.** HSS with the lead-0 number applied at 12-48 h stays 0.49 to 0.52 (Atlantic) and 0.52 to 0.57 (Pacific); AUC of G flat. This repeats the A1.2 finding of no measurable decay. The lead-specific matched number gives the same picture (post-run sensitivity in `raw_gust.json`, `lead_specific_matched`).
- **Compared with the ERA5 proxy** at 71.7 kt on the same archive and months, GFS lead 0 is lower by 0.07 (Pacific) to 0.09 (Atlantic) HSS: the ERA5 proxy is the better match to OPC (POD 0.65 to 0.67 against 0.54 to 0.58). The ERA5 71.7 kt value was count-matched on these same seasons (in-sample for ERA5); the GFS number is held out. Using 71.7 kt on GFS gives the same HSS as the held-out value (0.533 / 0.568).
- **Compared with the state model probability**, raw G is almost as good a ranker at 24 h: the model adds 0.00 (Atlantic) and 0.014 (Pacific) AUC. The model contains `g800`, so it cannot be worse in expectation; most of its skill at these labels is the gust. (PR 12/178 model BSS and RPSS were against different labels and were not recomputed here.)
- **By storm (S3, descriptive):** of 234 OPC storms whose first HF fix fell in Oct-Apr 2021-25, 228 (97%) had a GFS track fix within 600 km at the nearest 00/12 UTC. At that time 105 of 228 (46%) had a lead-0 fix above the matched threshold (the first HF fix is the weakest point of the event); 145 of 228 (64%) had some forecast, valid at that time with lead up to 48 h, above it. Of those 145, the earliest flagging lead was 48 h for 95, 36 h for 22, 24 h for 17, 12 h for 4 and 0 h for 7. This is a loose analogue of the P >= 0.5 lead in the OPC-vs-model comparison (PR 161, median 12 h, ERA5 tracks, 2004-25), not the same experiment, and 48 h is the cap of the forecast leads scored.

## What this cannot show
Five seasons, so the season-block intervals are crude (the matched-threshold HSS intervals at lead 24 are very narrow because the thresholds barely change across folds; they are not a statement about other seasons). One forecast model with version changes inside the window. The tracker must find the low first (tracker misses count as misses in S3 only). The HSS depends on the base rate of the fix sample, which differs between the GFS and the ERA5 lists, so the R3 gap is best read with POD and FAR beside it. Gust beyond 48 h, other models, and storms in May to September are not covered. Pipeline A (ERA5) is a proxy.

## Verification
(to be filled by the independent check)
