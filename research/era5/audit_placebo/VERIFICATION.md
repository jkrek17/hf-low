# Verification

A fresh agent (not the one that produced the numbers; told not to read notes) recomputed from committed files with its own Poisson and logistic fits.

| Item | Claimed | Recomputed | Verdict |
|---|---|---|---|
| Atl exp(A2), exp(A4), committed oos_index | 1.252, 1.27 | 1.2517, 1.2711 | pass |
| Pac exp(A2), exp(A4) | 1.137, 1.13 | 1.1374, 1.1287 | pass |
| Part B Atl n, HF, b_M1, b_M2 | 4521, 901, 0.287, 0.199 | same | pass |
| Part B Pac n, HF, b_M1, b_M2 | 4947, 788, 0.169, 0.083 | same | pass |
| Part C files / directories / rows / q<0.05 | 114 / 42 / 4,225 / 1,258 | same | pass |

## Second verifier: Part A (own weekly counts, decomposition, bootstrap and shuffle code; imported only `hemlib`'s fitting functions)
- Pass: P0-P3 point estimates of A4, A2, f_conv in both basins; A4 and A2 95% intervals within bootstrap noise; r with the committed index (0.938, 0.852); Pacific P2 constant in all 22 folds.
- Did not reproduce: f_conv upper bounds for Atl P2, Pac P1, Pac P3 (heavy-tailed, A4 near 0); Atl P2 f_conv lower bound (0.83-1.00 against 1.06; the claimed value drops draws with A4 <= 0). The pre-registered reading is the same on all of them.
- Shuffle: medians and 95th percentiles within 0.05. The exceedance rates in the first version of the README (1% Atl, 11% Pac of 100 draws) were wrong: 47 and 53 of the 100 shuffled indices are constant and S is undefined for them. Corrected to 1 of 53 and 11 of 47; verifier 3% and 7%.

**Still not independently recomputed:** permutation p-values, BH q-values, Part B's M3 and attenuation intervals.
