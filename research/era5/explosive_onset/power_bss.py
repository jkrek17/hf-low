"""Planted-effect power for the primary metric (stratum BSS gain over the PR 76 tree baseline), labels and baseline
probabilities only; no ERA5 field is read.

A feature z ~ N(0,1) independent of the baseline is planted in the truth: logit p* = logit(p_M1) + beta z + c (c keeps
the mean probability). Outcomes y* are redrawn for all 135,663 stratum fixes, the pull design is applied to y* (census of
times with a case, 600 sampled other times, weight = n_other/600), z is "measured" without error, and the corrective model
logit p = logit(p_M1) + a + b z is fitted leave-one-season-out with weights. Detection = sign-flip p < 0.05 on the 22 season
Brier differences (BSS metric) or on the 22 season held-out log-likelihood differences (within-time metric), positive point estimate. Reported: mean DeltaBSS (versus the stratum climatology) and power.
Run: python3 -I research/era5/explosive_onset/power_bss.py <repo_root> [nsim]
"""
import sys, os, numpy as np, pandas as pd
root = sys.argv[1]; NSIM = int(sys.argv[2]) if len(sys.argv) > 2 else 100
f = pd.read_csv(os.path.join(root, "research/era5/intensity/results/fixes_2004.csv.gz"), usecols=["track", "time", "season", "g800", "hf_now", "hf24"])
p = pd.read_csv(os.path.join(root, "research/era5/hf_boosted/results/loso_probs.csv.gz"), usecols=["track", "time", "M1"])
d = f.merge(p, on=["track", "time"])
d = d[(d.g800 < 55) & (~d.hf_now)].reset_index(drop=True)
pb = d.M1.clip(1e-5, 1 - 1e-5).to_numpy(); lb = np.log(pb / (1 - pb)); tm = d.time.to_numpy(); se = d.season.to_numpy()
ut, inv = np.unique(tm, return_inverse=True)
seasons = np.unique(se); rng = np.random.RandomState(7)
expit = lambda x: 1 / (1 + np.exp(-x))
def fit(lbo, z, y, w, iters=12):
    X = np.c_[np.ones(len(z)), z]; b = np.zeros(2)
    for _ in range(iters):
        q = expit(lbo + X @ b); g = X.T @ (w * (y - q)) - np.r_[0, 1.0 * b[1]]
        H = (X * (w * q * (1 - q))[:, None]).T @ X + np.diag([0, 1.0]); b = b + np.linalg.solve(H, g)
    return b
def one(beta):
    z = rng.randn(len(d)); c = 0.0
    if beta > 0:
        lo, hi = -3, 3
        for _ in range(30):
            c = (lo + hi) / 2
            if expit(lb + beta * z + c).sum() > pb.sum(): hi = c
            else: lo = c
    ys = (rng.rand(len(d)) < expit(lb + beta * z + c)).astype(float)
    ncase_t = np.bincount(inv, weights=ys, minlength=len(ut)) > 0
    other = np.where(~ncase_t)[0]; samp = rng.choice(other, 600, replace=False)
    wt = np.zeros(len(ut)); wt[ncase_t] = 1; wt[samp] = len(other) / 600
    w = wt[inv]; keep = w > 0
    diffs = []; sb = 0.0; sa = 0.0
    for s in seasons:
        tr = keep & (se != s); te = keep & (se == s)
        b = fit(lb[tr], z[tr], ys[tr], w[tr])
        qa = expit(lb[te] + b[0] + b[1] * z[te])
        bsb = (w[te] * (pb[te] - ys[te]) ** 2).sum(); bsa = (w[te] * (qa - ys[te]) ** 2).sum()
        diffs.append(bsb - bsa); sb += bsb; sa += bsa
    diffs = np.array(diffs); pi = (w * ys).sum() / w.sum() if False else ys.mean()
    clim = pi * (1 - pi) * (w.sum())
    dbss = (sb - sa) / clim
    flips = rng.choice([-1, 1], size=(4000, len(diffs)))
    pval = (np.abs((flips * diffs).sum(1)) >= abs(diffs.sum())).mean()
    # within-time (Breslow-type conditional) metric on case-time risk sets: held-out log-likelihood gain per case
    ci = keep & ncase_t[inv]; g = inv
    d2 = []; ncs = 0
    for s in seasons:
        tr = ci & (se != s); te = ci & (se == s)
        b = 0.0
        for _ in range(10):
            q = np.exp(lb[tr] + b * z[tr]); S0 = np.bincount(g[tr], q, len(ut)); S1 = np.bincount(g[tr], q * z[tr], len(ut)); S2 = np.bincount(g[tr], q * z[tr] ** 2, len(ut))
            nc = np.bincount(g[tr], ys[tr], len(ut)); zc = np.bincount(g[tr], ys[tr] * z[tr], len(ut)); ok = S0 > 0
            U = (zc[ok] - nc[ok] * S1[ok] / S0[ok]).sum() - b; I = (nc[ok] * (S2[ok] / S0[ok] - (S1[ok] / S0[ok]) ** 2)).sum() + 1.0
            b += U / I
        def ll(bb):
            q = np.exp(lb[te] + bb * z[te]); S0 = np.bincount(g[te], q, len(ut)); nc = np.bincount(g[te], ys[te], len(ut)); ok = nc > 0
            return (ys[te] * (lb[te] + bb * z[te])).sum() - (nc[ok] * np.log(S0[ok])).sum()
        d2.append(ll(b) - ll(0.0)); ncs += ys[te].sum()
    d2 = np.array(d2); flips2 = rng.choice([-1, 1], size=(4000, len(d2)))
    p2 = (np.abs((flips2 * d2).sum(1)) >= abs(d2.sum())).mean()
    return dbss, (pval < 0.05) and dbss > 0, d2.sum() / max(ncs, 1), (p2 < 0.05) and d2.sum() > 0
L = [f"stratum fixes {len(d)}, planted-effect simulations {NSIM} per beta (beta = log odds ratio per SD of an independent feature)"]
for beta in [0.0, 0.1, 0.2, 0.3, 0.45]:
    r = [one(beta) for _ in range(NSIM)]
    L.append(f"beta {beta:.2f}: BSS metric: mean DeltaBSS {np.mean([x[0] for x in r]):+.4f} (sd {np.std([x[0] for x in r]):.4f}), detected {np.mean([x[1] for x in r]):.2f}; within-time metric: held-out log-lik gain per case {np.mean([x[2] for x in r]):+.4f}, detected {np.mean([x[3] for x in r]):.2f}")
    print(L[-1], flush=True)
open(os.path.join(root, "research/era5/explosive_onset/results/power_bss.txt"), "w").write("\n".join(L) + "\n")
