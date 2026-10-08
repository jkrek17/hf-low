# Verification (RA-25 stage 1)

A fresh Sonnet agent recomputed the numbers below from `weekly_table.csv.gz` and the PR 41 index files, with its own code
(`verify/verify_numbers.py`), having read only PREREGISTRATION.md (not `run.py`, `posthoc.py`, the results or the README).

**Matched (all):** primary r (ARCH -0.0933, PRX -0.1061, DEP -0.0970) and their permutation p (0.0179, 0.0069, 0.0148 from 200,000 shuffles);
S1 month-dummy r; S2a season-total r; S3 NAO/PNA-removed r; index r -0.600 and ratio 0.155; halves (ARCH -0.136/-0.050, PRX -0.151/-0.055);
share of weeks with both counts zero (0.129, 0.121, 0.124); post hoc 5- and 10-week block r; ARCH variance/mean 0.953 Atlantic, 0.872 Pacific.

**Not independently checked:** bootstrap intervals; the planted-effect power and the minimum detectable r; BH q values and the 13-of-16 count;
p values of S2a, S2b, S3 and S4 (r values for S2a, S3 matched; S2b and S4 r not recomputed); the lag profile and its envelope; S7 r after dropping
both-zero weeks; post hoc blocks of 1, 2, 3 and 6 weeks and their intervals; variance/mean for the proxy samples and DEP; the DEP halves.
The verifier assumed the order of the index files (season x week) from the brief and did not test it; `run.py` asserts the same order for
the archive column only through the shared table.
