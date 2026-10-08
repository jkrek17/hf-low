# Verification, life-cycle test (fresh Sonnet agent)

Inputs given: `PREREGISTRATION.md`, `storms.csv`, `fixes_2004.csv.gz`, `lc_sample.csv`, `lc_values.csv.gz` and the result CSVs for comparison. No code, README or notes read. Script and output: `verify/verify_numbers.py`, `verify/out.txt`.

**Matched exactly:** stage A 28 tests (n, means, differences; q pass/fail on all 28; p within tolerance on 26); stage B 56 tests (n, means, differences on all 56; 49 fully inside tolerance); coverage 28 of 28 rows; first persistent lag point estimates for all four stage A cases and seven of eight stage B cases; g800 reproduction 2,140 of 2,140 rows within 0.5 kt (max 0.00).

**Differences:** bootstrap p at large p (eight rows, all above 0.05, differences 0.02 to 0.07); Pacific pressure gradient at -36 h: committed q 0.040, verifier 0.053 to 0.059 on three seeds and a pass on one, so the Pacific pressure-gradient first lag is -36 h (committed) or -24 h (verifier). Atlantic 48-kt radius at -12 h and Atlantic max-gust radius at +12 h also flip on other seeds.

**Not independently checked:** extraction, heading and rotation, pixel maps and figures, the resample distribution of the first lag, CI columns, the 84 'complete' sensitivity rows.

**Input notes:** two storms in `storms.csv` have no fixes; basin label differs between `storms.csv` and `fixes_2004` for some tracks (storms.csv used); reference SF means are reweighted onto HF strata (the pre-registered estimator).
