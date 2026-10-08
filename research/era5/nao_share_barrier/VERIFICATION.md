# Independent recomputation

A fresh Sonnet agent that had not seen `barrier_share.py`, the README or the results recomputed the numbers below from
`all_tracks.csv.gz`, `gustloc_fixes.csv`, `gh_daily.csv`, the CPC NAO file and `PREREGISTRATION.md`, with its own
Poisson fits (`verify/verify_numbers.py`). All matched.

**Matched:** sample counts (9,636 tracks, 1,003 HF, T 894, G100 947, G300 791, R400 862, barrier-only 109, N60 6,207/642,
D 9,527/894, 4,606 days); RR(HF) 1.096, RR(all) 0.974, RR(share) 1.125; RR(share) for T 1.141, G100 1.129, G300 1.132,
R400 1.150, N60 1.091, D 1.141; storm-HF count 1.111 and barrier-only 0.979; corr(NAO, GH) -0.693; NAO|GH share 1.039;
GH share 0.870 (original), 0.862 (T); bootstrap SE of the T-minus-original difference 0.0125 (reported 0.0124), of the
NAO change when GH enters 0.0315 (0.0318); original-HF interval 1.047-1.194.

**Scope note:** the 9,636 and 1,003 are tracks whose genesis day has a non-missing NAO index; the window holds 9,656
tracks and 1,007 HF (T 898, G100 951, G300 794, R400 866). One HF track has no fix rows and counts as HF in all variants.

**Not independently checked:** permutation p-values, BH q-values, leave-one-season-out ranges, the interval for r,
bootstrap intervals other than those listed, the minimum detectable effects, and the GH check against PR 58 (done by
`gh_check.py`, same author).
