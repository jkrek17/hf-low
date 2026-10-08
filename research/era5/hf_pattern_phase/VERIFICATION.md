# Verification, pattern-phase test (fresh Sonnet agent)

Inputs: `PREREGISTRATION.md`, `hf_index.csv`, `pp_values.csv.gz`, result CSVs for comparison; no code, README or notes. Script and output: `verify/verify_numbers.py`, `verify/out.txt`.

**Matched:** all counts (1,007 Atlantic and 825 Pacific in scope, tercile counts, 168 out of scope), group means, stratified differences, SD, d_std on all 99 rows (primary 22, S1 22, S2 22, S3 11, S4 22); q pass/fail agrees on every row; SE and MDE within 1 to 5%.

**Differences:** bootstrap p off by more than 0.02 on 30 of 99 rows (largest 0.064, pac gust_factor 0.753 vs 0.817); the only one crossing 0.05 is S4 pac gmax (0.075 vs 0.047; q 0.165 vs 0.115, both above 0.05). The committed 'bottom' column is the top mean minus the difference (reweighted onto strata present in both); the S2 n column counts all storms (1,007/825), not those with a right-of-motion share (928/768).

**Not independently checked:** the resample stream, tercile cut points (from the 660 weekly index values, not an input given to the verifier), the extraction of the structure scalars, maps, figures. The g800 check is circular here (g800 equals the catalog value by construction in the table); see deviation 1.
