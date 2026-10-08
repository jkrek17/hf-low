# Verification (RA-7)

A fresh Sonnet agent, given only `PREREGISTRATION.md` and the committed inputs (not my code or results), recomputed the numbers below.
Script and output: `verify_numbers.py`, `verify_output.txt`.

**Matched:** sample (1,886; 1,039 / 847); late shares 0.293 / 0.340 / 0.235, LATE6 0.109, LATE12 0.047; predictor prevalences T1-T4;
odds ratios T1-T7 (1.41, 3.10, 2.96, 0.45, 0.35, 0.49, 0.59); 552 late storms, single HF fix 0.612, onset on peak-gust fix 0.804,
onset on the minimum fix 63% of late, onset pressure within 2 hPa of the minimum 0.784; median pressure gap 8.2 hPa (early) and 2.7 hPa
(late by 6 h or more); joint leave-one-season-out AUC 0.764 (0.568 basin and month only), Brier skill 0.144.

**Differences, by convention:** T5 prevalence is 64.3% / 77.4% only with the 38 events lacking maxdeep in the denominator (65.9% / 78.7% over
defined events); the README does not quote it. The basin-month standardised late shares match to within 0.005 when weighted the verifier's way (largest gap T4 without: 0.407 against 0.412);
the raw shares are within 0.02. Brier skill 0.144 is against a basin-month reference; against the overall mean it would be 0.166.

**Not independently checked:** confidence intervals (clustered and bootstrap), p and q values, power and minimum detectable effects, the S1-S9
environment tests, the Atlantic-only joint model and the post hoc joint intervals, heterogeneity and sensitivity tables, the late share by month group and by number of HF fixes,
the peak-gust-minus-minimum percentiles, the 18.4% same-fix share (it is in the verifier's h_on_minp==0 count only indirectly).
