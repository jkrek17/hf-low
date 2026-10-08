"""Planted-effect power for the within-time (conditional) test, using the real risk sets (labels only; no field read).

A planted feature z is N(0,1) for controls and N(delta,1) for cases, independent of everything else. The test is the
conditional-logit score test (sum over case times of [case z - n_case * risk-set mean z], standardised).
alpha 0.05 / 6 groups (Bonferroni-style stand-in for the FDR family) and 0.05 uncorrected are both reported.
Run: python3 -I research/era5/explosive_onset/power.py <repo_root>
"""
import sys, os, numpy as np, pandas as pd
from scipy import stats
root = sys.argv[1]
f = pd.read_csv(os.path.join(root, "research/era5/intensity/results/fixes_2004.csv.gz"))
E = f[(f.g800 < 55) & (~f.hf_now)]
ct = E[E.time.isin(E[E.hf24].time.unique())]
sets = [(int(d.hf24.sum()), len(d)) for _, d in ct.groupby("time")]
nc = np.array([a for a, b in sets]); n = np.array([b for a, b in sets])
rng = np.random.RandomState(1)
def run(delta, nsim=1500, alpha=(0.05, 0.05 / 6)):
    hit = np.zeros(len(alpha))
    for _ in range(nsim):
        S = 0.0; V = 0.0
        for a, b in zip(nc, n):
            z = rng.randn(b); z[:a] += delta       # first a are cases
            m = z.mean(); S += z[:a].sum() - a * m
            V += a * (b - a) / (b * (b - 1)) * (np.sum((z - m) ** 2)) if b > 1 else 0
        p = 2 * stats.norm.sf(abs(S / np.sqrt(V)))
        hit += np.array([p < al for al in alpha])
    return hit / nsim
L = [f"risk sets (case times) {len(sets)}, cases {nc.sum()}, controls {(n-nc).sum()}"]
for d in [0.05, 0.10, 0.15, 0.20, 0.25, 0.30]:
    p = run(d, 400)
    L.append(f"delta {d:.2f} SD: power {p[0]:.2f} at alpha 0.05, {p[1]:.2f} at alpha 0.05/6")
open(os.path.join(root, "research/era5/explosive_onset/results/power.txt"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
