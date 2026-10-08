"""RA-30 step 3: weekly Poisson decomposition of the pattern effect (PR 85 Part 1 design, A2 = P(HF | deepening)) with the
cyclone denominator from tracker V or M, against pipeline A on the same 18 seasons. ERA5 PROXY. Plan: PREREGISTRATION.md.
usage: python3 run.py OUTDIR TRACKTABLE_DIR [NPERM] [NBOOT]
"""
import json, os, sys, time
from multiprocessing import Pool
import numpy as np, pandas as pd

OUT, TT = sys.argv[1], sys.argv[2]
NPERM = int(sys.argv[3]) if len(sys.argv) > 3 else 2000
NBOOT = int(sys.argv[4]) if len(sys.argv) > 4 else 2000
ERA = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, f"{ERA}/hem_channels")
import chanlib as C
C.SEASONS = list(range(2004, 2022))                  # 18 seasons; RA-15 coverage ends January 2023
NS, SEED = len(C.SEASONS), 20261011
os.makedirs(OUT, exist_ok=True)


def index(basin):
    o = pd.read_csv(C.OOS.format(basin)).iloc[: NS * C.NWEEK]
    assert (o.season.values == np.repeat(C.SEASONS, C.NWEEK)).all()
    x = o.idx.values
    return (x - x.mean()) / x.std()


def tables():
    A = C.load_tracks(first=2003, last=2021)
    f = pd.read_csv(C.FIXES, usecols=["track", "dp12"])
    out = {"A": {}}
    A["_d36"] = A.track.isin(set(f.track[f.dp12 <= -3.6])).astype(int)
    A["_d24"] = A.track.isin(set(f.track[f.dp12 <= -2.4])).astype(int)
    out["A"] = A
    for t in ("M", "V"):
        d = pd.read_csv(f"{TT}/tracks_{t}.csv.gz", parse_dates=["date"])
        out[t] = d[d.basin.isin(C.BASINS)].reset_index(drop=True)
    return out


def weekly(T, basin, hf, deep, mask=None):
    t = T[T.basin == basin]
    if mask is not None:
        t = t[mask(t)]
    h, dp = t[hf] == 1, t[deep] == 1
    masks = dict(N=pd.Series(True, index=t.index), D=dp, HD=h & dp)
    win = C.assign_week(t)
    Y, LP = {}, {}
    for k, m in masks.items():
        Y[k] = C.weekly_counts(win, m.reindex(win.index))
        LP[k] = np.log1p(C.prev_week_counts(t, m, basin))
    return Y, LP, win


def rows_of(idx):
    return np.concatenate([np.arange(i * C.NWEEK, (i + 1) * C.NWEEK) for i in idx])


def Xmat(M, LP, x, n, rows):
    return np.column_stack([M[rows], LP[n][rows], x[rows]])


def slopes(M, D, x, rows):
    Y, LP = D
    return {n: C.pois(Xmat(M, LP, x, n, rows), Y[n][rows])[-1] for n in ("N", "D", "HD")}


def fac(s):
    return {"A1": s["D"] - s["N"], "A2": s["HD"] - s["D"]}


def analyse(job):
    label, basin, DS, x, nperm, nboot, seed = job          # DS: {name: (Y, LP)}
    M = C.month_dummies()
    rng = np.random.default_rng(seed)
    sidx = list(range(NS))
    rows = rows_of(sidx)
    cl = np.repeat(np.arange(NS), C.NWEEK)
    obs = {n: fac(slopes(M, DS[n], x, rows)) for n in DS}
    xs = x[rows].reshape(NS, C.NWEEK)
    null = {n: {k: np.empty(nperm) for k in ("A1", "A2")} for n in DS}
    for i in range(nperm):
        xp = x.copy(); xp[rows] = xs[rng.permutation(NS)].reshape(-1)
        for n in DS:
            f = fac(slopes(M, DS[n], xp, rows))
            for k in f:
                null[n][k][i] = f[k]
    bt = {n: {k: np.empty(nboot) for k in ("A1", "A2")} for n in DS}
    for i in range(nboot):
        r = rows_of(rng.integers(0, NS, NS))
        for n in DS:
            f = fac(slopes(M, DS[n], x, r))
            for k in f:
                bt[n][k][i] = f[k]
    res = []
    for n in DS:
        Y, LP = DS[n]
        for k, contr in (("A1", ("N", "D")), ("A2", ("D", "HD"))):
            se = C.joint_se([(Xmat(M, LP, x, c, rows), Y[c][rows]) for c in contr], cl, [-1.0, 1.0])
            lo, hi = np.percentile(bt[n][k], [2.5, 97.5])
            p = (1 + int((np.abs(null[n][k]) >= abs(obs[n][k]) - 1e-12).sum())) / (1 + nperm)
            res.append(dict(variant=label, basin=basin, dataset=n, test=k, rr=np.exp(obs[n][k]), rr_lo=np.exp(lo), rr_hi=np.exp(hi),
                            log=obs[n][k], se_log=se, p=p, mde_rr=float(np.exp(2.8 * se)),
                            nN=int(Y["N"].sum()), nD=int(Y["D"].sum()), nHD=int(Y["HD"].sum())))
    if "A" in DS:                                         # paired difference (same resampled seasons)
        for n in DS:
            if n == "A":
                continue
            d = bt[n]["A2"] - bt["A"]["A2"]
            res.append(dict(variant=label, basin=basin, dataset=n + "-A", test="A2diff", rr=np.exp(obs[n]["A2"] - obs["A"]["A2"]),
                            rr_lo=np.exp(np.percentile(d, 2.5)), rr_hi=np.exp(np.percentile(d, 97.5)), log=obs[n]["A2"] - obs["A"]["A2"]))
    return res


def main():
    t0 = time.time()
    TB = tables()
    idx = {b: index(b) for b in C.BASINS}
    jobs = []

    def ds(b, hf=("hf", "hf", "hf"), deep=("_d36", "deep", "deep"), mask=None, only=("A", "V", "M")):
        out = {}
        for n, h, d in zip(("A", "V", "M"), hf, deep):
            if n not in only:
                continue
            Y, LP, _ = weekly(TB[n], b, h, d, mask)
            out[n] = (Y, LP)
        return out

    for b in C.BASINS:
        jobs.append(("primary", b, ds(b), idx[b], NPERM, NBOOT, SEED))
        jobs.append(("S2_dp12_-2.4", b, ds(b, deep=("_d24", "deep24", "deep24")), idx[b], NPERM, NBOOT, SEED))
        jobs.append(("S4_dist_half", b, ds(b, hf=("hf", "hf_half", "hf_half")), idx[b], NPERM, NBOOT, SEED))
        jobs.append(("S4_dist_double", b, ds(b, hf=("hf", "hf_dbl", "hf_dbl")), idx[b], NPERM, NBOOT, SEED))
    for lab, m in (("S1_atl_ge60N", lambda t: t.peak_lat >= 60), ("S1_atl_lt60N", lambda t: t.peak_lat < 60)):
        jobs.append((lab, "atl", ds("atl", mask=m), idx["atl"], NPERM, NBOOT, SEED))
    with Pool(4) as P:
        res = P.map(analyse, jobs, chunksize=1)
    T = pd.DataFrame([r for rr in res for r in rr])
    T["q"] = np.nan
    pri = (T.variant == "primary") & (T.test == "A2") & T.dataset.isin(["V", "M"])
    T.loc[pri, "q"] = C.bh(T.p[pri].values)                           # registered family of 4
    for v in T.variant.unique():                                      # own family per secondary variant
        if v == "primary":
            continue
        m = (T.variant == v) & T.p.notna() & T.dataset.isin(["V", "M"])
        T.loc[m, "q"] = C.bh(T.p[m].values)
    allp = T.p.notna() & T.dataset.isin(["V", "M"])
    T.loc[allp, "q_all"] = C.bh(T.p[allp].values)
    T.to_csv(f"{OUT}/tests.csv", index=False)
    print(T[T.variant == "primary"].to_string(), f"\n{time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
