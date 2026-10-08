# Verification (RA-29)

A fresh Sonnet agent recomputed with its own code (`verify_numbers.py`, raw output `VERIFICATION_raw.txt`), from the committed inputs and the definitions only, not the scripts or results. **All 36 claim rows matched to the rounding shown, no mismatches:**
- T-A: class totals (902/675/453 Atlantic, 811/609/430 Pacific); RR per SD without and with the era column for h1, h2, h3 in both basins; the six era coefficients; the index-era correlation; Pacific brief-class era coefficient and RR.
- T-B: per-season counts, exposure (5.48 / 16.52 years), events per year before and after at 1, 2, 3 fixes, shares at 2 and 3 fixes.
- T-C: mean hours at HF and event counts by era and source; single-fix shares.
- T-E: listed shares by era.
- T-D: fix counts, positives, and HSS (ceiling and PR 12 forecast) overall, before, after.

Notes from the verifier: the malformed date in `HF_Data_-_Pac.csv` (ID 2024202506, 20241101018) was dropped and the T-D figures still reproduce; h1 means all class-`low` lows including 4 per period with hfN 0 (as PR 41); 224 archive HF fixes found no track in the T-D matching.

**Not independently checked:** every bootstrap interval, p, q and the BH results; the ASCAT-B and trend sensitivity rows; the sustained-versus-brief and other paired-contrast values; the T-A MDE and SE-ratio figures; the era-coefficient intervals; the differences in T-C and T-D (they follow from the checked means and HSS, but the arithmetic was not rerun).
