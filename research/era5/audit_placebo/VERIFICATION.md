# Verification

A fresh agent (not the one that produced the numbers; told not to read notes) recomputed from committed files with its own Poisson and logistic fits.

| Item | Claimed | Recomputed | Verdict |
|---|---|---|---|
| Atl exp(A2), exp(A4), committed oos_index | 1.252, 1.27 | 1.2517, 1.2711 | pass |
| Pac exp(A2), exp(A4) | 1.137, 1.13 | 1.1374, 1.1287 | pass |
| Part B Atl n, HF, b_M1, b_M2 | 4521, 901, 0.287, 0.199 | same | pass |
| Part B Pac n, HF, b_M1, b_M2 | 4947, 788, 0.169, 0.083 | same | pass |
| Part C files / directories / rows / q<0.05 | 114 / 42 / 4,225 / 1,258 | same | pass |

**Not independently recomputed:** the placebo indices and their decomposition (Part A: P0-P3 index values, correlations with the committed index, A1-A4 and f_conv, the shuffle null),
all bootstrap intervals, permutation p-values, q-values, Part B's M3 and the attenuation intervals. The verifier was not given the index builder (`placelib.py`), so `r_committed` (0.938, 0.852) is unchecked.
