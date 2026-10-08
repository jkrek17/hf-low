# Verification, conversion test (fresh Sonnet agent)

Inputs: `PREREGISTRATION.md`, `fixes_2004.csv.gz`, `env_2004.csv.gz`, `pairs.csv`, `conv_sample.csv`, `conv_values.csv.gz`, result CSVs; no code, README or notes. Script and outputs (2,000 and 100,000 resamples): `verify/`.

**Matched exactly:** eligible 3,316 fixes on 2,274 storms; converting 1,287 on 966; Atlantic 1,475 / 1,023 and 620 / 468; Pacific 1,841 / 1,251 and 667 / 498; 788 pairs (374 + 414), each pair same basin and month, |ndr24 difference| < 0.25 (max 0.249), no non-converter reused, 149 tracks in both groups; n, group means, stratified differences, SD and coverage on every row of stage A primary (48), change (16), latitude-matched (48), confounds (10) and stage B (30; coverage 12 of 12).

**Differences:** bootstrap p and MDE (SE) differ within 2,000-resample noise (rows: 10, 6, 14, 1, 13); four q flips: stage A atl jet250 t0 (committed 0.048, verifier 0.053 at 2,000 and 0.049 at 100,000), stage A atl eady t0-24 h (0.057 vs 0.049), latitude-matched atl jet250 t0 (0.058 vs 0.038), stage B atl stab lag 0 (0.060 vs 0.0435). Standardising SD is ddof 0 in stage A and ddof 1 in stage B (under 1% effect).

**Not independently checked:** the random draw and matching algorithm, stage B extraction and heading rotation, the 25 km and 150 km position rules, stage B scalar values (taken from `conv_values.csv.gz`), maps, figures, CI columns.
