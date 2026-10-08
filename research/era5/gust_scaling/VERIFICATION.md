# Verification (fresh Sonnet agent, plan and committed inputs only; it did not read `scaling.py`, `diag.py`, `results/` or this README)

Script `verify/verify_numbers.py`, output `verify/verify_output.txt`.

**Matched exactly (4 decimals):** kept fixes, storms, seasons and the 49 dropped for Dp < 5 (32 Atlantic, 17 Pacific);
all ten M2 truncated-normal slopes and both sigmas (0.1282, 0.1241); all M3 slopes; the OLS contrast slopes
(depth 0.0835 / 0.0700); the five correlations per basin; the KS distances (0.0104, 0.0158). Its own sigma profile and joint
fit land on the same maximum.

**Bootstrap (500 draws, other seed, versus 2,000 here):** SEs agree within about 10% (for example Atlantic depth 0.0325 vs
0.0320, Pacific latitude 0.0619 vs 0.0619). Same 6 of 14 tests pass q < 0.05. The Pacific PNA-to-radius test sits on the line:
q 0.0502 (verifier) vs 0.0524 (here). It is reported as a lead that just misses, not as a finding either way.

**Not independently checked:** the 2,000-draw intervals and p values themselves, the minimum detectable effects, the ONI
shrink difference, all secondary models S1-S10 (including the GPD fits and the Req size sign flip), the profile table, the
percentage and kt conversions in the README (my arithmetic), q across all 193 rows.

Differences in handling: 2 missing daily PNA values (13 fixes) were averaged over the available days by the verifier; this
analysis' window mean also skips missing days (nanmean). Standardising SD: verifier ddof 0, this analysis ddof 1 (changes
ONI and PNA slopes by under 0.1%).
