# Verification (fresh Sonnet agent, 2026-10-09)
Recomputed independently from fixes_2004.csv.gz, comp_rot.npz and the per-storm boxes; did not read analyse.py or select.py.
- MATCH at printed precision: population (1,223; 664/559; missing headings 27/18; fix counts), all scalars S1, S2, S4 (values and positions), S3 shares and median distances, raw Z500 minimum positions.
- Independent re-extraction of 3 storm-times from one WB2 chunk (~30 MB): Z500 and Z500 anomaly agree to <1e-3 m; u/v250 exact. ws250 matches only as grid speed interpolated (not speed of interpolated u,v); stated in README.
- Five sentences were unsupported and were corrected in README: Pacific stacking headline; 'at the composite minimum' wording; ridge size; 'convergence behind it' (weak, only to 12 h before onset); 500 hPa divergence 'order of magnitude' (now 2 to 4 times).
- NOT checked: bootstrap intervals and BH stippling (agree with scalars.csv but not recomputed); the 16 figures' pixel content; ws500 box re-extraction; the climatology source itself.
