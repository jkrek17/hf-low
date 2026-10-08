"""RA-5 planted-effect power for P1..P4 (PREREGISTRATION.md, "Power"). Pacific series is shuffled in season blocks inside
every simulation, exactly as in the real test (1,000 shuffles). One-sided; reports p < 0.05 and p < 0.05/6 (the
Bonferroni level for the six primary tests, which BH can only loosen).
P1/P2: Poisson counts drawn from the fitted null (no Pacific term) times exp(log(RR) * standardised Pacific index).
P3/P4: null fit plus planted slope for a target partial r, null residuals resampled in season blocks.
usage: power.py [NSIM] [NPERM]"""
import json
import os
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd

import leadlib as L
import run as R

NSIM = int(sys.argv[1]) if len(sys.argv) > 1 else 400
NPERM = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
GRID = {"pois": [1.03, 1.05, 1.08, 1.10, 1.15, 1.20], "ols": [0.03, 0.05, 0.10, 0.15]}
DS = L.load_2004("archive")
SPEC = {s.name: s for s in R.PRIMARY[:4]}
PERMS = [np.random.default_rng(L.SEED + 700 + i).permutation(DS.ns) for i in range(NPERM)]


def cols(sp):
    y, X = L.design(DS, sp)
    xs = [L.design(DS, sp, perm=p)[1][:, -1] for p in PERMS]
    return y, X[:, :-1], X[:, -1], np.array(xs)


def sim_one(args):
    name, g, seed = args
    sp = SPEC[name]
    rng = np.random.default_rng(seed)
    y, Xc, xt, xp = cols(sp)
    if sp.kind == "pois":
        b = L.C.pois(Xc, y)
        mu = np.exp(Xc @ b)
        hit = np.zeros(2)
        for _ in range(NSIM):
            ys = rng.poisson(mu * np.exp(np.log(g) * xt))
            obs = L.C.pois(np.column_stack([Xc, xt]), ys)[-1]
            nul = np.array([L.C.pois(np.column_stack([Xc, v]), ys)[-1] for v in xp])
            p = (1 + np.sum(nul >= obs)) / (1 + len(nul))
            hit += [p < 0.05, p < 0.05 / 6]
        return name, g, hit / NSIM
    b, *_ = np.linalg.lstsq(Xc, y, rcond=None)
    fit0, res0 = Xc @ b, (y - Xc @ b)
    P = Xc @ np.linalg.pinv(Xc)  # projection on controls
    rx = xt - P @ xt
    delta = g / np.sqrt(1 - g * g) * res0.std() / rx.std()
    R0 = res0.reshape(DS.ns, -1)
    RX = xp - xp @ P.T  # residualised permuted columns (nperm, n)
    hit = np.zeros(2)
    for _ in range(NSIM):
        ys = fit0 + delta * xt + R0[rng.integers(0, DS.ns, DS.ns)].ravel()
        ry = ys - P @ ys
        obs = (rx @ ry) / (rx @ rx)
        nul = (RX @ ry) / np.einsum("ij,ij->i", RX, RX)
        p = (1 + np.sum(nul >= obs)) / (1 + len(nul))
        hit += [p < 0.05, p < 0.05 / 6]
    return name, g, hit / NSIM


def dispersion():
    out = {}
    for nm in ("P1", "P2"):
        y, Xc, xt, _ = cols(SPEC[nm])
        mu = np.exp(Xc @ L.C.pois(Xc, y))
        out[nm] = float(np.sum((y - mu) ** 2 / mu) / (len(y) - Xc.shape[1]))
    return out


if __name__ == "__main__":
    jobs = [(n, g, L.SEED + 1000 * i + j) for i, n in enumerate(SPEC) for j, g in enumerate(GRID[SPEC[n].kind])]
    with Pool(4) as pool:
        res = pool.map(sim_one, jobs)
    T = pd.DataFrame([dict(test=n, planted=g, power_p05=h[0], power_bonf6=h[1]) for n, g, h in res])
    T.to_csv(os.path.join(R.OUT, "power.csv"), index=False)
    json.dump(dict(dispersion_null_fit=dispersion(), nsim=NSIM, nperm=NPERM), open(os.path.join(R.OUT, "power_meta.json"), "w"))
    print(T.round(3).to_string()); print(dispersion())
