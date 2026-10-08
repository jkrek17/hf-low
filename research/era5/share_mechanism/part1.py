"""RA-17: which factor carries the hemispheric pattern's share effect? ERA5 PROXY, pipeline A.
Plan, committed before any weekly deepening count met an index: PREREGISTRATION.md (cb4892f).

usage: python3 -I part1.py REPO OUTDIR [NPERM] [NBOOT]      (no ERA5 access; reads committed files only)
"""
import json, os, sys, time
from multiprocessing import Pool

import numpy as np
import pandas as pd

REPO, OUT = sys.argv[1], sys.argv[2]
NPERM = int(sys.argv[3]) if len(sys.argv) > 3 else 2000
NBOOT = int(sys.argv[4]) if len(sys.argv) > 4 else 2000
sys.path.insert(0, os.path.join(REPO, "research", "era5", "hem_channels"))
sys.path.insert(0, os.path.join(REPO, "research", "era5", "hemispheric"))
import chanlib as C  # noqa: E402
from indices import Indices  # noqa: E402

SEED = 20261011
os.makedirs(OUT, exist_ok=True)
NS = len(C.SEASONS)


# ------------------------------------------------------------------ data
def deep_tracks(thr):
    f = pd.read_csv(C.FIXES, usecols=["track", "dp12"])
    return set(f.track[f.dp12 <= thr].unique())


def load(thr=-3.6):
    T = C.add_origin(C.load_tracks(first=2003, last=2025))
    T["deep"] = T.track.isin(deep_tracks(thr)).astype(int)
    return T


def weekly(Tfull, basin, hf_col, mask=None):
    t = Tfull[Tfull.basin == basin]
    if mask is not None:
        t = t[mask(t)]
    hf, dp = t[hf_col] == 1, t.deep == 1
    masks = dict(N=pd.Series(True, index=t.index), D=dp, HD=hf & dp, H=hf)
    win = C.assign_week(t)
    Y, LP = {}, {}
    for k, m in masks.items():
        Y[k] = C.weekly_counts(win, m.reindex(win.index))
        LP[k] = np.log1p(C.prev_week_counts(t, m, basin))
    return Y, LP, win


def index_series(name, basin):
    """weekly window -7..-1 before the week start, standardised over the 660 weeks (S3 only)."""
    import datetime as dt
    if name == "gh":
        g = pd.read_csv(os.path.join(REPO, "research/era5/nao_share_barrier/gh_daily.csv"), parse_dates=["date"]).set_index("date").GH.asfreq("D")
        v = [g.loc[pd.Timestamp(C.week_start(s, k) - dt.timedelta(days=7)):pd.Timestamp(C.week_start(s, k) - dt.timedelta(days=1))].mean()
             for s in C.SEASONS for k in range(C.NWEEK)]
    else:
        I = Indices(os.path.join(REPO, "docs/data/teleconnections.json"))
        v = [I.window(C.week_start(s, k) - dt.timedelta(days=7), C.week_start(s, k) - dt.timedelta(days=1))[name]
             for s in C.SEASONS for k in range(C.NWEEK)]
    v = np.array(v, float)
    assert not np.isnan(v).any()
    return (v - v.mean()) / v.std()


# ------------------------------------------------------------------ estimation
NAMES = ["N", "D", "HD", "H"]
CONTR = {"A1": [("N", -1), ("D", 1)], "A2": [("D", -1), ("HD", 1)], "A3": [("HD", -1), ("H", 1)], "A4": [("N", -1), ("H", 1)]}


def Xmat(M, LP, x, name, rows):
    return np.column_stack([M[rows], LP[name][rows], x[rows]])


def slopes(M, LP, Y, x, rows):
    return {n: C.pois(Xmat(M, LP, x, n, rows), Y[n][rows])[-1] for n in NAMES}


def factors(s):
    return {"A1": s["D"] - s["N"], "A2": s["HD"] - s["D"], "A3": s["H"] - s["HD"], "A4": s["H"] - s["N"], "bN": s["N"], "bH": s["H"]}


def rows_of(idx):
    return np.concatenate([np.arange(i * C.NWEEK, (i + 1) * C.NWEEK) for i in idx])


def analyse(job):
    label, basin, Y, LP, x, sidx, nperm, nboot, seed = job
    M = C.month_dummies()
    rng = np.random.default_rng(seed)
    rows = rows_of(sidx)
    ns = len(sidx)
    cl = np.repeat(np.arange(ns), C.NWEEK)
    obs = factors(slopes(M, LP, Y, x, rows))
    xs = x[rows].reshape(ns, C.NWEEK)
    null = {k: np.empty(nperm) for k in ("A1", "A2", "A3")}
    for i in range(nperm):
        xp = x.copy(); xp[rows] = xs[rng.permutation(ns)].reshape(-1)
        f = factors(slopes(M, LP, Y, xp, rows))
        for k in null:
            null[k][i] = f[k]
    bt = {k: np.empty(nboot) for k in ("A1", "A2", "A3", "A4", "bN", "bH")}
    for i in range(nboot):
        pick = rng.integers(0, ns, ns)
        r = rows_of([sidx[j] for j in pick])
        f = factors(slopes(M, LP, Y, x, r))
        for k in bt:
            bt[k][i] = f[k]
    out = dict(label=label, basin=basin, n_seasons=ns)
    for k in ("A1", "A2", "A3"):
        contr = [c for _, c in CONTR[k]]
        se = C.joint_se([(Xmat(M, LP, x, n, rows), Y[n][rows]) for n, _ in CONTR[k]], cl, contr)
        p = (1 + int((np.abs(null[k]) >= abs(obs[k]) - 1e-12).sum())) / (1 + nperm)
        lo, hi = np.percentile(bt[k], [2.5, 97.5])
        out[k] = dict(est=obs[k], se=se, p=p, lo=lo, hi=hi, rr=np.exp(obs[k]), rr_lo=np.exp(lo), rr_hi=np.exp(hi), mde_rr=float(np.exp(2.8 * se)))
    lo, hi = np.percentile(bt["A4"], [2.5, 97.5])
    out["A4"] = dict(est=obs["A4"], lo=lo, hi=hi, rr=np.exp(obs["A4"]), rr_lo=np.exp(lo), rr_hi=np.exp(hi))
    out["bN"] = dict(est=obs["bN"], rr=np.exp(obs["bN"]), lo=np.exp(np.percentile(bt["bN"], 2.5)), hi=np.exp(np.percentile(bt["bN"], 97.5)))
    out["bH"] = dict(est=obs["bH"], rr=np.exp(obs["bH"]))
    ok = np.sign(obs["A4"]) * bt["A4"] > 0.01          # direction-neutral guard (negative-effect indices, S3 Greenland high)
    out["frac_boot_ok"] = float(ok.mean())
    for nm, k in (("f_conv", "A2"), ("f_deep", "A1"), ("f_res", "A3")):
        fb = np.where(ok, bt[k] / np.where(ok, bt["A4"], 1), np.nan)
        out[nm] = dict(est=obs[k] / obs["A4"], lo=float(np.nanpercentile(fb, 2.5)), hi=float(np.nanpercentile(fb, 97.5)))
    return out


def verdict(r):
    if not (r["A4"]["lo"] > 0 or r["A4"]["hi"] < 0) or r["frac_boot_ok"] < 0.95:
        return "share effect not resolved or unstable; no verdict"
    if r["f_conv"]["lo"] > 0.5:
        v = "conversion P(HF|deepening) carries most"
    elif r["f_deep"]["lo"] > 0.5:
        v = "deepening P(deepening) carries most"
    else:
        v = "mixed / unresolved"
    if v.startswith("conversion"):
        v += "; P17 " + ("supported (f_conv >= 0.70)" if r["f_conv"]["est"] >= 0.70 else "half right (f_conv < 0.70)")
    elif r["f_conv"]["hi"] < 0.5:
        v += "; P17 falsified"
    return v


def bh(p):
    return C.bh(np.asarray(p))


def main():
    t0 = time.time()
    jobs = []
    Tbase = load(-3.6)
    idx_all = {b: C.index(b)[0] for b in C.BASINS}
    allidx = list(range(NS))
    h1, h2 = list(range(0, 11)), list(range(11, 22))
    for b in C.BASINS:
        Y, LP, win = weekly(Tbase, b, "hf")
        jobs.append(("primary", b, Y, LP, idx_all[b], allidx, NPERM, NBOOT, SEED))
        Yd, LPd, _ = weekly(Tbase, b, "hfd")
        jobs.append(("S1_fixed_depth", b, Yd, LPd, idx_all[b], allidx, NPERM, NBOOT, SEED))
        jobs.append(("S4a_2004-14", b, Y, LP, idx_all[b], h1, NPERM, NBOOT, SEED))
        jobs.append(("S4b_2015-25", b, Y, LP, idx_all[b], h2, NPERM, NBOOT, SEED))
    Ya, LPa, _ = weekly(Tbase, "atl", "hf", mask=lambda t: t.peak_lat <= 60)
    jobs.append(("S5_atl_le60N", "atl", Ya, LPa, idx_all["atl"], allidx, NPERM, NBOOT, SEED))
    for thr, lab in ((-2.4, "S2a_dp12_-2.4"), (-6.0, "S2b_dp12_-6.0")):
        Tt = load(thr)
        for b in C.BASINS:
            Y, LP, _ = weekly(Tt, b, "hf")
            jobs.append((lab, b, Y, LP, idx_all[b], allidx, NPERM, NBOOT, SEED))
    for nm, b in (("nao", "atl"), ("gh", "atl"), ("pna", "pac")):
        Y, LP, _ = weekly(Tbase, b, "hf")
        jobs.append((f"S3_{nm}", b, Y, LP, index_series(nm, b), allidx, NPERM, NBOOT, SEED))
    print("jobs", len(jobs), f"{time.time()-t0:.0f}s", flush=True)
    with Pool(4) as P:
        res = P.map(analyse, jobs, chunksize=1)
    rows = []
    for r in res:
        r["verdict"] = verdict(r)
        for k in ("A1", "A2", "A3"):
            rows.append(dict(variant=r["label"], basin=r["basin"], test=k, rr=r[k]["rr"], rr_lo=r[k]["rr_lo"], rr_hi=r[k]["rr_hi"],
                             log=r[k]["est"], se_log=r[k]["se"], p=r[k]["p"], mde_rr=r[k]["mde_rr"],
                             f=r[{"A1": "f_deep", "A2": "f_conv", "A3": "f_res"}[k]]["est"],
                             f_lo=r[{"A1": "f_deep", "A2": "f_conv", "A3": "f_res"}[k]]["lo"],
                             f_hi=r[{"A1": "f_deep", "A2": "f_conv", "A3": "f_res"}[k]]["hi"],
                             A4_rr=r["A4"]["rr"], A4_lo=r["A4"]["rr_lo"], A4_hi=r["A4"]["rr_hi"], RR_N=r["bN"]["rr"], RR_N_lo=r["bN"]["lo"], RR_N_hi=r["bN"]["hi"],
                             frac_boot_ok=r["frac_boot_ok"], verdict=r["verdict"]))
    T = pd.DataFrame(rows)
    T["q"] = np.nan
    for v in T.variant.unique():                      # one BH family per variant (primary: 6 tests)
        m = T.variant == v
        T.loc[m, "q"] = bh(T.p[m].values)
    T.to_csv(os.path.join(OUT, "part1_tests.csv"), index=False)
    json.dump(res, open(os.path.join(OUT, "part1_raw.json"), "w"), indent=1, default=float)
    print("done", f"{time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
