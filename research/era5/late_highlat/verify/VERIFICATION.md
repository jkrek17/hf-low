# Verification (RA-26)

A fresh Sonnet agent, given the pre-registration and committed inputs but not my script, README or results, reimplemented the RA-7 filters and recomputed the numbers. Its script and output: `verify_numbers.py`, `verify_output.txt`.

**Matched exactly** (counts and shares): group n 283, late 151, non-late 132, LATE6 66, terrain-flagged 145; P1 barrier share 0.3841 (58/151) vs 0.3939 (52/132); coast <= 300 km 0.980 vs 0.894; centre distance > 400 km 0.391 vs 0.417 and <= 400 km 0.609 vs 0.583 (S7); northerly 0.540 (81/150) vs 0.5455 (72/132); S4 0.4545 (30/66) vs 0.3687 (80/217); S5 0.3382 (23/68) vs 0.3571 (25/70); S8 0.1287 (26/202) vs 0.0650 (36/554).
The P1 season-block bootstrap interval was reproduced independently at [-0.119, +0.120] (seed 7), and is stable to about +-0.005 across seeds 1, 2 and 123.

**Not independently checked:** S6 (barrier at the peak-gust fix), the other bootstrap intervals and p values, BH q values, the S9 adjusted logistic model, the power at a 0.25 difference and the minimum detectable differences.

Ambiguities the verifier noted, none changing a number: "onset fix" as earliest row or `groupby.first()` give identical values; bootstrap resamples the 22 seasons within the group; northerly boundaries inclusive.
