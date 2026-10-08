"""RA-28 S3: PR 85 Part 2 (mediation of the pattern's effect on P(HF | deepening) by the near environment) with the HF label swapped.
ERA5 PROXY, pipeline A. Plan: PREREGISTRATION.md. Adapted from share_mechanism/part2.py at PR 85's head 67ddfe6; estimation code unchanged,
the label is a parameter (hf_col in ref, M, P99).

usage: python3 -I med.py REPO OUTDIR [NPERM] [NBOOT]
"""
import json, os, sys, time
from multiprocessing import Pool

import numpy as np
import pandas as pd

SEED = 20261011
PRIMARY = ["jet250", "eady", "sstgrad", "sst_t500", "flux", "tcwv", "div300", "vadv500"]
SECOND = ["lat", "B", "VTL", "VTU", "sst"]
DEEP = -3.6
MONTHS = [10, 11, 12, 1, 2, 3, 4]


# ------------------------------------------------------------------ models
def logit_fit(X, y, it=40, ridge=1e-6):
    b = np.zeros(X.shape[1])
    b[:] = 0.0
    m = y.mean()
    # start from intercept-only through the first dummy columns: plain zeros are fine for IRLS here
    for _ in range(it):
        eta = np.clip(X @ b, -30, 30)
        p = 1 / (1 + np.exp(-eta)); w = np.maximum(p * (1 - p), 1e-9)
        s = np.linalg.solve(X.T @ (X * w[:, None]) + ridge * np.eye(len(b)), X.T @ (y - p))
        b = b + s
        if np.abs(s).max() < 1e-8:
            break
    return b


def sig(e):
    return 1 / (1 + np.exp(-np.clip(e, -30, 30)))


def one_pass(M, x, Z, y, parts=True):
    """All quantities for one sample. M: month dummies (n x 7 or n x 1 intercept); x: pattern value; Z: n x k ingredients."""
    n, k = Z.shape
    Xt = np.column_stack([M, x])
    bt = logit_fit(Xt, y)
    out = {"beta_tot": bt[-1]}
    # path a
    Ma = np.column_stack([M, x])
    coef = np.linalg.lstsq(Ma, Z, rcond=None)[0]
    a = coef[-1]
    out["a"] = a
    XZ = np.column_stack([M, Z, x])

    def gform(cols):
        Xo = np.column_stack([M, Z[:, cols], x])
        bo = logit_fit(Xo, y)
        p0 = sig(Xo @ bo)
        Xi = Xo.copy(); Xi[:, M.shape[1]:M.shape[1] + len(cols)] += a[cols][None, :]
        Xd = Xo.copy(); Xd[:, -1] += 1.0
        nie = (sig(Xi @ bo) - p0).mean(); nde = (sig(Xd @ bo) - p0).mean()
        tot = nie + nde
        return nie, nde, (nie / tot if abs(tot) > 1e-12 else np.nan), bo[-1]
    nie, nde, Mfrac, bdir = gform(list(range(k)))
    out.update(NIE=nie, NDE=nde, M=Mfrac, beta_dir=bdir, M_diff=1 - bdir / bt[-1])
    if parts:
        out["M_single"] = np.array([gform([j])[2] for j in range(k)])
        out["M_unique"] = np.array([Mfrac - gform([i for i in range(k) if i != j])[2] for j in range(k)])
        out["b"] = np.array([logit_fit(np.column_stack([M, Z[:, j]]), y)[-1] for j in range(k)])
    return out


# ------------------------------------------------------------------ data
def build(ROOT, ref="B", hf_col="hf", mask=None, ingredients=PRIMARY, month_fixed=True, thr=DEEP):
    sys.path.insert(0, os.path.join(ROOT, "research", "era5", "hem_channels"))
    import chanlib as C
    T = C.load_tracks(first=2004, last=2025)
    T = T[T.basin.isin(C.BASINS)]
    if mask is not None:
        T = T[mask(T)]
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from lab import labels
    L, _ = labels(ROOT)
    for k, v in L.items():
        T["hf_" + k] = T.track.map(v).fillna(0).astype(int)
    T = C.assign_week(T)
    f = pd.read_csv(C.FIXES)
    e = pd.read_csv(os.path.join(ROOT, "research/era5/intensity/results/env_2004.csv.gz"))
    f = f.merge(e, on=["track", "time"], how="left", validate="1:1").sort_values(["track", "time"]).reset_index(drop=True)
    f["tt"] = pd.to_datetime(f.time.astype(str), format="%Y%m%d%H")
    f["i"] = f.groupby("track").cumcount()
    f["dt_prev"] = f.groupby("track").tt.diff().dt.total_seconds() / 3600
    first = f[f.dp12 <= thr].groupby("track").head(1).set_index("track")
    f2 = f.set_index(["track", "i"])
    cols = list(dict.fromkeys(PRIMARY + SECOND))
    rows = []
    ref_rows = {}
    for tr, r in first.iterrows():
        if ref == "A":
            ref_rows[tr] = r[cols].values.astype(float)
        else:
            if r.i >= 1 and r.dt_prev == 12:
                ref_rows[tr] = f2.loc[(tr, r.i - 1), cols].values.astype(float)
    R = pd.DataFrame.from_dict(ref_rows, orient="index", columns=cols)
    D = T.set_index("track").join(R, how="inner")                 # deepening cyclones with a reference fix
    D = D[D.hf.notna()]
    meta = dict(n_window_cyclones=int(len(T)), n_deepening_with_ref=int(len(D)))
    meta["n_deepening"] = int(T.track.isin(first.index).sum())
    ok = D[ingredients].notna().all(1)
    meta["n_complete_case"] = int(ok.sum())
    D = D[ok].copy()
    D["y"] = D[hf_col].astype(float)
    D["mon"] = D.date.dt.month
    D["x"] = np.nan
    out = {}
    for b in C.BASINS:
        xb, _ = C.index(b)
        d = D[D.basin == b].copy()
        d["x"] = xb[d.row.values]
        Zr = d[ingredients].astype(float)
        Z = ((Zr - Zr.mean()) / Zr.std(ddof=0)).values
        M = (np.array([[1.0 if m == mm else 0.0 for mm in MONTHS] for m in d.mon]) if month_fixed else np.ones((len(d), 1)))
        out[b] = dict(M=M, x=d.x.values, Z=Z, y=d.y.values, sea=d.sea.values, row=d.row.values, xweek=xb, n=len(d), hf=int(d.y.sum()))
    return out, meta


# ------------------------------------------------------------------ analysis
def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = np.empty(n); prev = 1.0
    for rank, i in zip(range(n, 0, -1), o[::-1]):
        prev = min(prev, p[i] * n / rank); q[i] = prev
    return q


def run_job(job):
    label, basin, d, names, nperm, nboot, seed = job
    rng = np.random.default_rng(seed)
    M, x, Z, y, sea = d["M"], d["x"], d["Z"], d["y"], d["sea"]
    pt = one_pass(M, x, Z, y)
    # permutation p of beta_tot: season blocks of the weekly index moved between seasons
    ns_all = 22
    xw = d["xweek"].reshape(ns_all, 30)
    row = d["row"]
    null = np.empty(nperm)
    Xp = np.column_stack([M, x])
    for i in range(nperm):
        xp = xw[rng.permutation(ns_all)].reshape(-1)[row]
        Xp[:, -1] = xp
        null[i] = logit_fit(Xp, y)[-1]
    p_tot = (1 + int((np.abs(null) >= abs(pt["beta_tot"]) - 1e-12).sum())) / (1 + nperm)
    seas = np.unique(sea)
    idx = {s: np.where(sea == s)[0] for s in seas}
    bt = {k: [] for k in ("beta_tot", "M", "M_diff", "NIE", "NDE")}
    ba, bM1, bMu, bb = [], [], [], []
    for _ in range(nboot):
        pick = rng.choice(seas, len(seas))
        ii = np.concatenate([idx[s] for s in pick])
        r = one_pass(M[ii], x[ii], Z[ii], y[ii])
        for k in bt:
            bt[k].append(r[k])
        ba.append(r["a"]); bM1.append(r["M_single"]); bMu.append(r["M_unique"]); bb.append(r["b"])
    bt = {k: np.array(v) for k, v in bt.items()}
    ba, bM1, bMu, bb = map(np.array, (ba, bM1, bMu, bb))
    pct = lambda v: np.nanpercentile(v, [2.5, 97.5])
    flagged = float((bt["beta_tot"] < 0.01).mean())
    res = dict(label=label, basin=basin, n=d["n"], hf=d["hf"], p_tot=p_tot, flagged=flagged, names=names)
    res["beta_tot"] = (pt["beta_tot"], *pct(bt["beta_tot"]))
    ok = bt["beta_tot"] >= 0.01
    for k in ("M", "M_diff"):
        res[k] = (pt[k], *np.nanpercentile(np.where(ok, bt[k], np.nan), [2.5, 97.5]))
    res["NIE"] = (pt["NIE"], *pct(bt["NIE"])); res["NDE"] = (pt["NDE"], *pct(bt["NDE"]))
    two = lambda v: min(1.0, 2 * (1 + min((v >= 0).sum(), (v <= 0).sum())) / (len(v) + 1))
    res["a"] = [(pt["a"][j], *pct(ba[:, j]), two(ba[:, j])) for j in range(len(names))]
    res["b"] = [(float(np.exp(pt["b"][j])), *np.exp(pct(bb[:, j])), two(bb[:, j]), float(pt["b"][j])) for j in range(len(names))]
    res["M_single"] = [(pt["M_single"][j], *np.nanpercentile(np.where(ok, bM1[:, j], np.nan), [2.5, 97.5])) for j in range(len(names))]
    res["M_unique"] = [(pt["M_unique"][j], *np.nanpercentile(np.where(ok, bMu[:, j], np.nan), [2.5, 97.5])) for j in range(len(names))]
    return res


def verdict(r):
    if r["p_tot"] >= 0.05:      # BH over 2 basins applied by caller via q_tot
        pass
    return None


def main():
    ROOT, OUT = sys.argv[1], sys.argv[2]
    nperm = int(sys.argv[3]) if len(sys.argv) > 3 else 2000
    nboot = int(sys.argv[4]) if len(sys.argv) > 4 else 2000
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    variants = {"ref": dict(hf_col="hf"), "M": dict(hf_col="hf_M"), "P99": dict(hf_col="hf_P99")}
    jobs, metas = [], {}
    for lab, kw in variants.items():
        D, meta = build(ROOT, **kw)
        metas[lab] = {**meta, **{b: dict(n=D[b]["n"], hf=D[b]["hf"]) for b in D}}
        names = kw.get("ingredients", PRIMARY)
        for b in D:
            jobs.append((lab, b, D[b], names, nperm, nboot, SEED))
    print("jobs", len(jobs), f"{time.time()-t0:.0f}s", flush=True)
    with Pool(4) as P:
        res = P.map(run_job, jobs, chunksize=1)
    rows = []
    for r in res:
        base = dict(variant=r["label"], basin=r["basin"], n=r["n"], hf=r["hf"])
        rows.append({**base, "quantity": "beta_tot", "ingredient": "", "estimate": r["beta_tot"][0], "lo": r["beta_tot"][1], "hi": r["beta_tot"][2], "p": r["p_tot"]})
        for k in ("M", "M_diff", "NIE", "NDE"):
            rows.append({**base, "quantity": k, "ingredient": "ALL8", "estimate": r[k][0], "lo": r[k][1], "hi": r[k][2], "p": np.nan, "flagged": r["flagged"]})
        for j, nm in enumerate(r["names"]):
            a, b, ms, mu = r["a"][j], r["b"][j], r["M_single"][j], r["M_unique"][j]
            rows.append({**base, "quantity": "a_path", "ingredient": nm, "estimate": a[0], "lo": a[1], "hi": a[2], "p": a[3]})
            rows.append({**base, "quantity": "b_path_OR", "ingredient": nm, "estimate": b[0], "lo": b[1], "hi": b[2], "p": b[3]})
            rows.append({**base, "quantity": "M_single", "ingredient": nm, "estimate": ms[0], "lo": ms[1], "hi": ms[2], "p": np.nan})
            rows.append({**base, "quantity": "M_unique", "ingredient": nm, "estimate": mu[0], "lo": mu[1], "hi": mu[2], "p": np.nan})
    T = pd.DataFrame(rows)
    T["q"] = np.nan
    for v in T.variant.unique():
        for qty in ("a_path", "b_path_OR"):
            m = (T.variant == v) & (T.quantity == qty)
            if m.any():
                T.loc[m, "q"] = bh(T.p[m].values)
        m = (T.variant == v) & (T.quantity == "beta_tot")
        T.loc[m, "q"] = bh(T.p[m].values)
    T.to_csv(os.path.join(OUT, "part2_tests.csv"), index=False, float_format="%.5g")
    json.dump(metas, open(os.path.join(OUT, "part2_meta.json"), "w"), indent=1)
    print("done", f"{time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
