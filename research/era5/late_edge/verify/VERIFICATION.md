# Verification (RA-26 follow-up, edge artefact)

A fresh Sonnet agent, given the pre-registration and committed inputs but not my script, README or results, recomputed the numbers (`verify_numbers.py`, `verify_output.txt`).

**Matched exactly:** all counts (Atlantic 1,039; high-latitude 283; touching 343, 129 of them high-latitude; buffered-touching 522, 212 high-latitude; non-touching 696 with 154 high-latitude and 542 south); late shares 0.5455 (84/154) vs 0.2878 (156/542);
odds ratios and Wald intervals for the full sample 3.14 [2.25, 4.38], E1 2.97 [2.04, 4.32], E2 3.95 [2.49, 6.26], E4 2.21 [1.37, 3.57], E5 2.28 [1.41, 3.66], E6 3.38 [2.23, 5.10], E7 2.51 [1.40, 4.48]; E3 interaction 1.33 [0.83, 2.13], p 0.238.

**Not independently checked:** the season-block bootstrap intervals, BH q values, power and minimum-detectable-OR figures.

Notes from the verifier, none changing a number: the full-sample OR is 3.14, not RA-7's 3.10 (cause not traced; both of us use no covariates beyond T2); E6 drops events whose onset fix is terrain-flagged; the buffered group has 212 high-latitude events including the 129 touching ones (my pre-registration guessed about 70 as the *non-touching* remainder; the actual remainder is 71).
