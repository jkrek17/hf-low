"""Planted-effect power for the weekly channel tests (PREREGISTRATION.md, "Power"). ERA5 PROXY, pipeline A.

Weekly counts are simulated hierarchically from the fitted baseline with the pattern slope set to zero:
    N_k ~ Poisson(mu_Nk x RR_count_k),  H_k ~ Binomial(N_k, min(0.95, s_k x RR_share_k)),  k = entrant, local,
with the previous-week term held at its observed value. The planted RR per SD multiplies one channel; the test is the same
clustered statistic as the real analysis (t on 21 df). Reports the RR with 80% power for T1-T8 and, for the position
tests, 2.94 x the clustered SE (t on 21 df at 80% power, analytic).

usage: power.py OUTDIR [NSIM]
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

import chanlib as C
import channels as CH

GRID = [1.03, 1.05, 1.08, 1.10, 1.15, 1.20]
# scenario -> (tests it powers, RR_count_e, RR_count_l, RR_share_e, RR_share_l) as functions of the planted RR
SCEN = {
    "T1": (["T1"], lambda r: (r, r, 1, 1)),
    "T2": (["T2", "T3"], lambda r: (1, 1, r, r)),
    "T4": (["T4"], lambda r: (r, 1, 1, 1)),
    "T5": (["T5"], lambda r: (1, r, 1, 1)),
    "T6": (["T6"], lambda r: (1, 1, r, 1)),
    "T7": (["T7"], lambda r: (1, 1, 1, r)),
    "T8": (["T8"], lambda r: (r, 1, 1, 1)),
}
NAMES = {"T1": ["N"], "T2": ["N", "H"], "T3": ["N", "H"], "T4": ["Ne"], "T5": ["Nl"], "T6": ["Ne", "He"], "T7": ["Nl", "Hl"],
         "T8": ["Ne", "Nl"]}


def simulate(base, rr, rng, x):
    rce, rcl, rse, rsl = [float(v) ** x for v in rr]
    Ne = rng.poisson(base["muNe"] * rce); Nl = rng.poisson(base["muNl"] * rcl)
    He = rng.binomial(Ne, np.minimum(0.95, base["se"] * rse)); Hl = rng.binomial(Nl, np.minimum(0.95, base["sl"] * rsl))
    return dict(N=(Ne + Nl).astype(float), H=(He + Hl).astype(float), Ne=Ne.astype(float), Nl=Nl.astype(float),
                He=He.astype(float), Hl=Hl.astype(float))


def main():
    outdir = sys.argv[1]
    nsim = int(sys.argv[2]) if len(sys.argv) > 2 else 300
    rng = np.random.default_rng(C.SEED + 2)
    Tfull = CH.prev_all_tracks()
    M = C.month_dummies()
    allidx = list(range(len(C.SEASONS)))
    rows = np.arange(len(C.SEASONS) * C.NWEEK)
    cl = np.repeat(np.arange(len(C.SEASONS)), C.NWEEK)
    tcrit = stats.t.ppf(0.975, len(C.SEASONS) - 1)
    out = []
    for b in C.BASINS:
        x, _ = C.index(b)
        Y, LP, win = CH.build(b, "hf", None, Tfull)
        D = CH.Design(M, x, Y, LP)
        beta = D.betas(rows, x, ["Ne", "Nl", "He", "Hl"])
        mu = {}
        for n in beta:
            bb = beta[n].copy(); bb[-1] = 0.0
            mu[n] = np.exp(np.clip(D.X(n, rows, x) @ bb, -25, 12))
        base = dict(muNe=mu["Ne"], muNl=mu["Nl"], se=mu["He"] / mu["Ne"], sl=mu["Hl"] / mu["Nl"])
        for sc, (tests, fn) in SCEN.items():
            for r in GRID:
                hits = {t: 0 for t in tests}
                for _ in range(nsim):
                    Ys = simulate(base, fn(r), rng, x)
                    # the x effect is planted through the rates, so the index column is the real x
                    Ds = CH.Design(M, x, Ys, LP)
                    s = Ds.slopes(rows, x, NAMES[tests[0]] if sc != "T2" else ["N", "H"])
                    tv = CH.tests_from_slopes({**{k: 0.0 for k in CH.OUTS_FULL}, **s})
                    for t in tests:
                        se = CH.se_of(Ds, rows, cl, t)
                        if abs(tv[t]) / se > tcrit:
                            hits[t] += 1
                for t in tests:
                    out.append(dict(basin=b, test=t, planted_rr=r, power=hits[t] / nsim))
        print(f"{b} power done", flush=True)
    P = pd.DataFrame(out)
    P.to_csv(os.path.join(outdir, "power_grid.csv"), index=False)
    # MDE: linear interpolation of power over the grid at 0.8; grid starts at 1.03, power at RR=1 is ~0.05
    rows_ = []
    for (b, t), g in P.groupby(["basin", "test"]):
        g = g.sort_values("planted_rr")
        r = np.r_[1.0, g.planted_rr.values]; p = np.r_[0.05, g.power.values]
        mde = np.interp(0.8, p, r) if p.max() >= 0.8 else np.nan
        rows_.append(dict(basin=b, test=t, mde_rr_80=mde, power_at_1p05=float(g.power[g.planted_rr == 1.05].iloc[0]),
                          power_at_1p10=float(g.power[g.planted_rr == 1.10].iloc[0])))
    pd.DataFrame(rows_).to_csv(os.path.join(outdir, "power_mde.csv"), index=False)
    T = pd.read_csv(os.path.join(outdir, "tests.csv")) if os.path.exists(os.path.join(outdir, "tests.csv")) else None
    if T is not None:
        pos = T[T.test.isin(["T9", "T10", "T11", "T12"])]
        pd.DataFrame(dict(basin=pos.basin, test=pos.test, mde_degrees_80=2.94 * pos.se_log)).to_csv(
            os.path.join(outdir, "power_position.csv"), index=False)


if __name__ == "__main__":
    main()
