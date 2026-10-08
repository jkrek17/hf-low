# Independent verification of pair-family numbers

Script: `verify_numbers.py` (written from PREREGISTRATION.md only; common.py, engine.py, run.py, power.py, count_pairs.py and results/ not read). Run: `python3 -I verify_numbers.py atl 300`.

## Matched
| Quantity | Atl claim / mine | Pac claim / mine |
|---|---|---|
| HF tracks | 1104 / 1104 | 896 / 896 |
| Candidate pairs | 160 / 160 | 60 / 60 |
| Rear pairs | 70 / 70 | 18 / 18 |
| Rear fraction | 0.4375 / 0.4375 | 0.300 / 0.300 |
| Obs-genesis daughters | 15 / 15 | 9 / 9 |
| Obs-genesis rear | 10 / 10 | 5 / 5 |
| Null mean f (300 perms, seed 20261008) | 0.418 / 0.4136 | 0.502 / 0.4994 |

One-sided permutation p (300 perms): Atlantic 0.229, Pacific 1.000. Null means differ from the claims by 0.004 and 0.003, within permutation noise (different RNG stream and 300 vs 2,000 draws).

## Caveat found: heading convention
Every daughter first fix falls exactly on a 00/12 UTC fix time of the parent, so "the two fixes bracketing" is ambiguous. Using the interval (t_i, t_i+1] with the fix at the daughter's time and the one before it (backward heading) reproduces 70 and 18 exactly. The forward interval [t_i, t_i+1) gives 55 (Atl) and 17 (Pac); a centered heading gives 61 and 20; the file's own `heading` column gives 67 and 22. The Atlantic rear count therefore depends on this choice (55 to 70), and PREREGISTRATION.md does not fix it. Recommend stating the backward convention in the plan.

## Not checked
- The 2,000-permutation run, BH correction, two-sided p, the conditional NAO/PNA test, power.py and the minimum detectable effect.
- The "16 percent of tracks start >12 h before first fix" statement and PR 37 excess figures (42, 24).
- The observed-genesis definition: I used lifecycle_events `gen_obs` and first fix within 12 h of `start` (all_tracks); the parenthetical "1000 hPa or above, north of 21N" was not applied separately (counts matched anyway).
- Null mechanics follow my reading (shift whole tracks by the permuted peak-time difference within basin x season x calendar month of peak); null-mean agreement supports it but is not proof of identical implementation.
