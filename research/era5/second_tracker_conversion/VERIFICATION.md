# Verification (RA-30)

A fresh Sonnet agent that had not seen `build.py`, `run.py`, `vpressure.py`, `results/tests.csv`, the per-track tables or this README implemented the design from `PREREGISTRATION.md` and the raw inputs (tracker fix files, `events_{M,V}.csv` of PR 116, pipeline A tracks and fixes, the PR 41 index) and recomputed the quoted point estimates. Its script and output: `verify/verify_numbers.py`, `verify/out.txt`.

**Matched exactly:**
- Track counts on read: 83,038 (M) and 50,040 (V).
- Weekly cyclone counts N, deepening D, and HF-and-deepening HD, for pipeline A, V and M in both basins (for example Atlantic V 9,317 / 4,828 / 629; Pacific M 4,262 / 2,486 / 558; pipeline A Atlantic 7,830 / 3,699 / 736).
- exp(A1), the deepening rate ratio, to three decimals in all six cases (for example Atlantic V 1.0125 against 1.0130).
- 30 of 30 randomly drawn V central pressures (`v_mslp500`) from the cached MSLP chunks, within 0.01 hPa.

**Matched to a stated cause:** exp(A2) differs by 0.002 to 0.006 (largest: pipeline A Atlantic 1.2470 against 1.2516; V Atlantic 1.2324 against 1.2356; M Atlantic 1.2082 against 1.2139). The verifier put each week's month dummy on the month of the week's start; `run.py` uses the month of the week's mid-point (start plus 3 days), the convention of PR 64/85. Re-running pipeline A Atlantic with start-of-week months gives 1.2470, the verifier's value, to four digits. Neither choice was fixed by the plan; the quoted numbers use the PR 85 convention so they compare with PR 85. The difference is far inside every interval.

**Not independently checked:** the bootstrap intervals, permutation p-values, Benjamini-Hochberg q, clustered standard errors and minimum detectable effects; the paired differences; all secondary variants (-2.4 hPa, distance halved or doubled, north and south of 60N); the label construction at doubled and halved distance; the 79 / 66 basin-less HF tracks; the PR 116 reproduction check in `build.py` (it was only checked by `build.py` itself against `events_{M,V}.csv`).
