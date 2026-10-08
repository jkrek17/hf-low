"""RA-28: PR 85 Part 1 (the weekly Poisson decomposition of the hemispheric pattern's share effect) re-run with the HF label swapped.
ERA5 PROXY, pipeline A. Plan: PREREGISTRATION.md (committed before any label or count). Adapted from share_mechanism/part1.py at PR 85's head 67ddfe6
(PR 85 is open, so this PR carries its own copy): the estimation functions are unchanged; the label is a parameter; a paired bootstrap is added.

usage: python3 -I conv.py REPO OUTDIR [NPERM] [NBOOT]      (reads committed files and results/fix_p99.csv.gz; no ERA5 access)
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
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lab import labels  # noqa: E402
from indices import Indices  # noqa: E402

SEED = 20261011
HERE_ = os.path.dirname(os.path.abspath(__file__))
os.makedirs(OUT, exist_ok=True)
NS = len(C.SEASONS)


# ------------------------------------------------------------------ data
def deep_tracks(thr):
    f = pd.read_csv(C.FIXES, usecols=["track", "dp12"])
    return set(f.track[f.dp12 <= thr].unique())


def load(thr=-3.6):
    T = C.add_origin(C.load_tracks(first=2003, last=2025))
    T["deep"] = T.track.isin(deep_tracks(thr)).astype(int)
    L, _ = labels(REPO)
    for k, v in L.items():
        T["hf_" + k] = T.track.map(v).fillna(0).astype(int)
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


def paired(job):
    """Same bootstrap picks for every label: draws of RR_A2 per label (and of A1, A3 shares), so Delta = RR(L) - RR(L_M) has a paired interval."""
    basin, Ys, LPs, x, sidx, nboot, seed = job
    M = C.month_dummies()
    rng = np.random.default_rng(seed)
    ns = len(sidx)
    draws = {k: np.empty(nboot) for k in Ys}
    for i in range(nboot):
        pick = rng.integers(0, ns, ns)
        r = rows_of([sidx[j] for j in pick])
        for k in Ys:
            s = slopes(M, LPs[k], Ys[k], x, r)
            draws[k][i] = s["HD"] - s["D"]
    return basin, draws


def main():
    t0 = time.time()
    Tbase = load(-3.6)
    idx_all = {b: C.index(b)[0] for b in C.BASINS}
    allidx = list(range(NS))
    LABS = {"ref": "hf", "M": "hf_M", "P99": "hf_P99", "P98": "hf_P98"}
    jobs, store = [], {}
    for b in C.BASINS:
        Ys, LPs = {}, {}
        for lab, col in LABS.items():
            Y, LP, _ = weekly(Tbase, b, col)
            Ys[lab], LPs[lab] = Y, LP
            jobs.append((lab, b, Y, LP, idx_all[b], allidx, NPERM, NBOOT, SEED))
        store[b] = (Ys, LPs)
    pj2 = []
    for b, (Ys, LPs) in store.items():
        pj2.append((b, Ys, LPs, idx_all[b], allidx, NBOOT, SEED))
    print("jobs", len(jobs), f"{time.time()-t0:.0f}s", flush=True)
    with Pool(4) as P:
        res = P.map(analyse, jobs, chunksize=1)
        pres = P.map(paired_wrap, pj2, chunksize=1)
    rows = []
    for r in res:
        r["verdict"] = verdict(r)
        for k in ("A1", "A2", "A3"):
            fk = {"A1": "f_deep", "A2": "f_conv", "A3": "f_res"}[k]
            rows.append(dict(label=r["label"], basin=r["basin"], test=k, rr=r[k]["rr"], rr_lo=r[k]["rr_lo"], rr_hi=r[k]["rr_hi"],
                             log=r[k]["est"], se_log=r[k]["se"], p=r[k]["p"], mde_rr=r[k]["mde_rr"],
                             f=r[fk]["est"], f_lo=r[fk]["lo"], f_hi=r[fk]["hi"],
                             A4_rr=r["A4"]["rr"], A4_lo=r["A4"]["rr_lo"], A4_hi=r["A4"]["rr_hi"], RR_N=r["bN"]["rr"], RR_N_lo=r["bN"]["lo"],
                             RR_N_hi=r["bN"]["hi"], frac_boot_ok=r["frac_boot_ok"], verdict=r["verdict"]))
    T = pd.DataFrame(rows)
    T["q_family"] = np.nan
    for lab in LABS:
        m = T.label == lab
        T.loc[m, "q_family"] = bh(T.p[m].values)          # 6 tests per label, as PR 85
    T.to_csv(os.path.join(OUT, "conv_tests.csv"), index=False)
    json.dump(res, open(os.path.join(OUT, "conv_raw.json"), "w"), indent=1, default=float)
    # paired differences
    prow = []
    for basin, draws in pres:
        rr = {k: np.exp(v) for k, v in draws.items()}
        for lab in ("ref", "P99", "P98"):
            d = rr[lab] - rr["M"]
            dl = draws[lab] - draws["M"]
            lo, hi = np.percentile(d, [2.5, 97.5])
            p = min(1.0, 2 * (1 + min(int((d <= 0).sum()), int((d >= 0).sum()))) / (len(d) + 1))
            est = {r["label"]: r for r in res if r["basin"] == basin}
            prow.append(dict(basin=basin, label=lab, rr_label=est[lab]["A2"]["rr"], rr_M=est["M"]["A2"]["rr"],
                             delta_rr=est[lab]["A2"]["rr"] - est["M"]["A2"]["rr"], lo=lo, hi=hi, se_delta_rr=float(d.std(ddof=1)),
                             mde_delta=float(2.8 * d.std(ddof=1)), log_ratio=est[lab]["A2"]["est"] - est["M"]["A2"]["est"],
                             lr_lo=np.percentile(dl, 2.5), lr_hi=np.percentile(dl, 97.5), p=p))
    PD = pd.DataFrame(prow)
    PD.to_csv(os.path.join(OUT, "paired_delta.csv"), index=False)
    print("done", f"{time.time()-t0:.0f}s")


def paired_wrap(job):
    basin, Ys, LPs, x, sidx, nboot, seed = job
    return paired((basin, Ys, LPs, x, sidx, nboot, seed))


if __name__ == "__main__":
    main()
