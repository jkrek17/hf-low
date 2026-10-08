# Independent verification (own implementation, verify_numbers.py)

Readings used: deepening point = midpoint of max-B 24 h window chosen over the whole track, then box (25-67N,120-240E) and DJF window (midpoint 03 Dec 00Z to before 01 Mar 00Z) applied. Winter-level WLS weights = storms/winter; z-scores over all 47 winters. Permutation = Freedman-Lane on the t-statistic (20000 perms, two seeds); bootstrap 5000 winter pairs. Power: null surrogate = reduced-model fit + permuted residuals, plus delta*EMI_z, 300 plants, 299 perms.
Alternative reading (restrict windows to DJF/box before taking max-B) gives bombs/winter 33.5 (21-49), also not exactly 34/22.
S12 = winter-mean CP minus EP (unweighted -5.25; bomb-weighted -4.91); the permutation test was not run.

# Verification output (raw, verify_numbers.py)

```
tracks total 42923 with deepening window 35773
bombs/winter mean,min,max 33.15 21 49 | rapid mean 70.34 55 79 | total bombs 1558
HF-reaching (DJF, Pacific box, deepening pt) per winter 2004-2025 mean 21.95 min 12 max 34 | HF tracks total in match 1209
check: DJF EMI from raw vs file max abs diff 0.0

## Point estimates (deg per SD)
P1 -1.088
P2 -0.07
S1 -0.116
S2 -0.027
S3 -0.056
S4 -0.066
S5 -1.63
S6 -1.081
S7 -1.063
S8 0.401
S9 -0.462
S10 -0.801
S11 -4.032
S13 -1.603
S14 -0.791
S15 -0.69
S16 -2.373
S12 El Nino winters 13 CP 4 EP 9
S12 CP-EP winter-mean (unweighted) -5.254 ; bomb-weighted -4.908 ; pooled storms -4.908
P1 FL perm p=0.1746 (rerun other seed 0.1754) t=-1.379; boot 95% CI [-2.473, 0.372]
P2 FL perm p=0.7518 (rerun other seed 0.7527) t=-0.316; boot 95% CI [-0.600, 0.352]

## Planted-effect power for P1 (null surrogate = reduced-model fit + permuted residuals)
delta 0 power 0.04
delta 1.5 power 0.383
delta 2.0 power 0.723
delta 2.5 power 0.88
delta 3.0 power 0.973

## Item 5
SD winter-mean bomb lon: unweighted 3.829 weighted 3.715 | lat unweighted 1.037 weighted 1.009
weighted ddof-corrected lon 3.755 lat 1.020
corr EMI-N34 0.732 ; EMI-PNA 0.275 ; EMI-SON EMI 0.879 ; EMI-N4mN3 0.244
```

## Note from the author (added after reading this report)

The verifier's circular longitude midpoint and its un-boxed lowest-pressure fix exposed two errors in the first run of `analysis.py`/`deepening_points.py` (README, "Post-hoc deviations"). After the fixes the author's results match the values above except for Monte-Carlo error in p and intervals. Not checked by the verifier: see README "Verification".
