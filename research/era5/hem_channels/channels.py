"""Which channel carries the hemispheric pattern's effect on HF lows? ERA5 PROXY, pipeline A.
Plan, committed before any outcome met the index: PREREGISTRATION.md.

usage: channels.py OUTDIR [NPERM] [NBOOT]      (about 15 minutes, no ERA5 access; reads committed files only)

Writes OUTDIR/tests.csv (all 24 primary tests), secondary.csv, decomposition.csv, summary.txt, meta.json.
"""
import json
import os
import sys
import time

import numpy as np
import pandas as pd

import chanlib as C

OUTS_FULL = ["N", "H", "Ne", "Nl", "He", "Hl"]
OUTS_CORE = ["N", "H"]


# ------------------------------------------------------------------ data
def build(basin, hf_col="hf", mask=None, Tfull=None):
    """Weekly arrays for one basin. mask: boolean over the track table (e.g. peak latitude <= 60)."""
    t_all = Tfull[Tfull.basin == basin]
    if mask is not None:
        t_all = t_all[mask(t_all)]
    ent = t_all.entrant == 1
    hf = t_all[hf_col] == 1
    masks = dict(N=pd.Series(True, index=t_all.index), H=hf, Ne=ent, Nl=~ent, He=hf & ent, Hl=hf & ~ent)
    win = C.assign_week(t_all)
    Y, LP = {}, {}
    for k, m in masks.items():
        Y[k] = C.weekly_counts(win, m.reindex(win.index))
        LP[k] = np.log1p(C.prev_week_counts(t_all.assign(_m=m), m, basin))
    return Y, LP, win


def prev_all_tracks():
    T = C.add_origin(C.load_tracks(first=2003, last=2025))
    return T


# ------------------------------------------------------------------ estimation
def X_of(M, lp, E, x, use_prev=True):
    cols = [M]
    if use_prev:
        cols.append(lp[:, None])
    if E is not None:
        cols.append(E)
    cols.append(x[:, None])
    return np.hstack(cols)


class Design:
    def __init__(self, M, x, Y, LP, E=None, use_prev=True):
        self.M, self.x, self.Y, self.LP, self.E, self.use_prev = M, x, Y, LP, E, use_prev

    def X(self, name, rows, x=None):
        x = self.x if x is None else x
        return X_of(self.M[rows], self.LP[name][rows], None if self.E is None else self.E[rows], x[rows], self.use_prev)

    def slopes(self, rows, x, names):
        return {n: C.pois(self.X(n, rows, x), self.Y[n][rows])[-1] for n in names}

    def betas(self, rows, x, names):
        return {n: C.pois(self.X(n, rows, x), self.Y[n][rows]) for n in names}


def tests_from_slopes(s):
    out = {"T1": s["N"], "T3": s["H"], "T2": s["H"] - s["N"]}
    if "Ne" in s:
        out.update(T4=s["Ne"], T5=s["Nl"], T6=s["He"] - s["Ne"], T7=s["Hl"] - s["Nl"], T8=s["Ne"] - s["Nl"])
    return out


CONTRASTS = {  # test -> list of (outcome, coefficient)
    "T1": [("N", 1)], "T3": [("H", 1)], "T2": [("N", -1), ("H", 1)],
    "T4": [("Ne", 1)], "T5": [("Nl", 1)], "T6": [("Ne", -1), ("He", 1)], "T7": [("Nl", -1), ("Hl", 1)],
    "T8": [("Ne", 1), ("Nl", -1)],
}


def se_of(D, rows, cl, test):
    outs = [(D.X(n, rows), D.Y[n][rows]) for n, _ in CONTRASTS[test]]
    return C.joint_se(outs, cl, [c for _, c in CONTRASTS[test]])


def analyse(D, seasons_idx, names, nperm, nboot, rng, tests=None, loso=True):
    """Estimates, clustered SE, season-block permutation p (two-sided) and bootstrap interval for the weekly tests."""
    rows = np.concatenate([np.arange(i * C.NWEEK, (i + 1) * C.NWEEK) for i in seasons_idx])
    cl = np.repeat(np.arange(len(seasons_idx)), C.NWEEK)
    s_obs = D.slopes(rows, D.x, names)
    t_obs = tests_from_slopes(s_obs)
    tests = tests or list(t_obs)
    ns = len(seasons_idx)
    xs = D.x[rows].reshape(ns, C.NWEEK)
    # permutation
    null = {k: np.empty(nperm) for k in t_obs}
    for i in range(nperm):
        xp = xs[rng.permutation(ns)].reshape(-1)
        xfull = D.x.copy(); xfull[rows] = xp
        t = tests_from_slopes(D.slopes(rows, xfull, names))
        for k in t_obs:
            null[k][i] = t[k]
    res = {}
    for k in tests:
        p = (1 + int((np.abs(null[k]) >= abs(t_obs[k]) - 1e-12).sum())) / (1 + nperm)
        res[k] = dict(est=t_obs[k], se=se_of(D, rows, cl, k), p=p)
    # bootstrap
    bt = {k: np.empty(nboot) for k in t_obs}
    for i in range(nboot):
        pick = rng.integers(0, ns, ns)
        r = np.concatenate([np.arange(seasons_idx[j] * C.NWEEK, (seasons_idx[j] + 1) * C.NWEEK) for j in pick])
        t = tests_from_slopes(D.slopes(r, D.x, names))
        for k in t_obs:
            bt[k][i] = t[k]
    for k in tests:
        res[k]["lo"], res[k]["hi"] = np.percentile(bt[k], [2.5, 97.5])
    with np.errstate(divide="ignore", invalid="ignore"):
        f_obs = t_obs["T2"] / t_obs["T3"]
        fb = bt["T2"] / bt["T3"]
    res["f"] = dict(est=f_obs, lo=np.nanpercentile(fb, 2.5), hi=np.nanpercentile(fb, 97.5),
                    share_gt0_5=float(np.nanmean(fb > 0.5)))
    if loso:
        lo_ = []
        for j in range(ns):
            keep = [seasons_idx[i] for i in range(ns) if i != j]
            r = np.concatenate([np.arange(i * C.NWEEK, (i + 1) * C.NWEEK) for i in keep])
            lo_.append(tests_from_slopes(D.slopes(r, D.x, names)))
        for k in ("T1", "T2", "T3"):
            v = [t[k] for t in lo_]
            res[k]["loso_lo"], res[k]["loso_hi"] = min(v), max(v)
    return res, bt


# ------------------------------------------------------------------ decomposition
def decomposition(D, rows, x):
    """Per-week absolute decomposition of the change in HF lows for +1 SD of the index (module docstring in the plan)."""
    b = D.betas(rows, x, ["Ne", "Nl", "He", "Hl"])
    out = {}
    tot = 0.0
    parts = {}
    base = {}
    for k, (nn, hh) in {"e": ("Ne", "He"), "l": ("Nl", "Hl")}.items():
        X0 = D.X(nn, rows, x); X1 = X0.copy(); X1[:, -1] += 1.0
        XH0 = D.X(hh, rows, x); XH1 = XH0.copy(); XH1[:, -1] += 1.0
        N0, N1 = np.exp(X0 @ b[nn]), np.exp(X1 @ b[nn])
        H0, H1 = np.exp(XH0 @ b[hh]), np.exp(XH1 @ b[hh])
        s0, s1 = H0 / N0, H1 / N1
        dN, ds = N1 - N0, s1 - s0
        parts["count_" + k] = float((dN * s0).sum())
        parts["share_" + k] = float((N0 * ds).sum())
        parts["inter_" + k] = float((dN * ds).sum())
        base[k] = float(H0.sum())
        tot += float((H1 - H0).sum())
    parts["total"] = tot
    parts["base_entrant_frac"] = base["e"] / (base["e"] + base["l"])
    parts["entrant_frac_of_change"] = (parts["count_e"] + parts["share_e"] + parts["inter_e"]) / tot if tot else np.nan
    for k in ("count", "share", "inter"):
        parts[k + "_frac"] = (parts[k + "_e"] + parts[k + "_l"]) / tot if tot else np.nan
    for k in ("count_e", "count_l", "share_e", "share_l", "inter_e", "inter_l"):
        parts[k + "_frac"] = parts[k] / tot if tot else np.nan
    return parts


def decomposition_ci(D, seasons_idx, nboot, rng):
    rows = np.concatenate([np.arange(i * C.NWEEK, (i + 1) * C.NWEEK) for i in seasons_idx])
    obs = decomposition(D, rows, D.x)
    ns = len(seasons_idx)
    draws = []
    for _ in range(nboot):
        pick = rng.integers(0, ns, ns)
        r = np.concatenate([np.arange(seasons_idx[j] * C.NWEEK, (seasons_idx[j] + 1) * C.NWEEK) for j in pick])
        draws.append(decomposition(D, r, D.x))
    df = pd.DataFrame(draws)
    lo, hi = df.quantile(0.025), df.quantile(0.975)
    ex = df.entrant_frac_of_change - df.base_entrant_frac
    obs["excess_entrant"] = obs["entrant_frac_of_change"] - obs["base_entrant_frac"]
    return obs, lo.to_dict(), hi.to_dict(), (float(ex.quantile(0.025)), float(ex.quantile(0.975)))


# ------------------------------------------------------------------ position
def position_tests(win, x, hf_only, nperm, rng, seasons_idx=None):
    """Track-level OLS of peak lat/lon on the index (month effects), season-clustered, permutation p."""
    t = win if not hf_only else win[win.hf == 1]
    Mt = C.month_of_date(t.date)  # month of first fix
    rows = t.row.values
    sea = (t.sea.values - C.SEASONS[0])
    out = {}
    ns = len(C.SEASONS)

    def resid(v):
        beta, *_ = np.linalg.lstsq(Mt, v, rcond=None)
        return v - Mt @ beta

    xt = x[rows]
    for name, y in (("lat", t.peak_lat.values), ("lon", t.lon_u.values)):
        yr = resid(y.astype(float))
        xr = resid(xt)
        b = (xr @ yr) / (xr @ xr)
        e = yr - b * xr
        g = np.zeros(ns)
        np.add.at(g, sea, xr * e)
        se = np.sqrt((g ** 2).sum() * ns / (ns - 1)) / (xr @ xr)
        xs = x.reshape(ns, C.NWEEK)
        null = np.empty(nperm)
        for i in range(nperm):
            xp = xs[rng.permutation(ns)].reshape(-1)[rows]
            xpr = resid(xp)
            null[i] = (xpr @ yr) / (xpr @ xpr)
        p = (1 + int((np.abs(null) >= abs(b) - 1e-12).sum())) / (1 + nperm)
        from scipy import stats
        tcrit = stats.t.ppf(0.975, ns - 1)
        out[name] = dict(est=b, se=se, p=p, lo=b - tcrit * se, hi=b + tcrit * se, n=len(t))
    return out


# ------------------------------------------------------------------ main
def main():
    outdir = sys.argv[1]
    nperm = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
    nboot = int(sys.argv[3]) if len(sys.argv) > 3 else 2000
    os.makedirs(outdir, exist_ok=True)
    t0 = time.time()
    Tfull = prev_all_tracks()
    M = C.month_dummies()
    allidx = list(range(len(C.SEASONS)))
    rng = np.random.default_rng(C.SEED)
    rows_all, dec_rows = [], []
    meta = dict(nperm=nperm, nboot=nboot, seed=C.SEED, n_tracks_window={}, missing_fix={})
    primary = {}
    for b in C.BASINS:
        x, y_arch = C.index(b)
        Y, LP, win = build(b, "hf", None, Tfull)
        meta["n_tracks_window"][b] = int(len(win))
        meta["missing_fix"][b] = int(win.missing_fix.sum())
        D = Design(M, x, Y, LP)
        res, bt = analyse(D, allidx, OUTS_FULL, nperm, nboot, rng)
        dobs, dlo, dhi, exc = decomposition_ci(D, allidx, nboot, rng)
        pos = {}
        for hf_only in (False, True):
            r = position_tests(win, x, hf_only, nperm, rng)
            for nm, v in r.items():
                pos[("hf_" if hf_only else "all_") + nm] = v
        primary[b] = dict(res=res, dec=(dobs, dlo, dhi, exc), pos=pos, D=D, win=win)
        print(f"{b} primary done {time.time()-t0:.0f}s", flush=True)
    # tidy the 24 tests
    labels = {"T1": "RR cyclones (count channel)", "T2": "RR share reaching HF (share channel)", "T3": "RR HF lows (total)",
              "T4": "RR entrants", "T5": "RR local storms", "T6": "RR share within entrants", "T7": "RR share within local",
              "T8": "entry contrast (entrant RR / local RR)"}
    posmap = {"T9": "all_lat", "T10": "all_lon", "T11": "hf_lat", "T12": "hf_lon"}
    for b in C.BASINS:
        res, pos = primary[b]["res"], primary[b]["pos"]
        for k in [f"T{i}" for i in range(1, 9)]:
            r = res[k]
            rows_all.append(dict(basin=b, test=k, label=labels[k], scale="RR per SD", est=np.exp(r["est"]), lo=np.exp(r["lo"]),
                                 hi=np.exp(r["hi"]), se_log=r["se"], p=r["p"],
                                 loso_lo=np.exp(r["loso_lo"]) if "loso_lo" in r else np.nan,
                                 loso_hi=np.exp(r["loso_hi"]) if "loso_hi" in r else np.nan))
        for k, pk in posmap.items():
            r = pos[pk]
            rows_all.append(dict(basin=b, test=k, label=("peak " + pk.replace("_", " ") + "itude" if False else f"peak {pk.split('_')[1]} of {'HF' if pk.startswith('hf') else 'all'} tracks"),
                                 scale="degrees per SD", est=r["est"], lo=r["lo"], hi=r["hi"], se_log=r["se"], p=r["p"],
                                 loso_lo=np.nan, loso_hi=np.nan))
    T = pd.DataFrame(rows_all)
    T["q_all24"] = C.bh(T.p.values)
    core = T.test.isin(["T1", "T2", "T3"])
    T["q_core6"] = np.nan
    T.loc[core, "q_core6"] = C.bh(T.p[core].values)
    T.to_csv(os.path.join(outdir, "tests.csv"), index=False)
    # decomposition table
    for b in C.BASINS:
        dobs, dlo, dhi, exc = primary[b]["dec"]
        for k, v in dobs.items():
            if k == "excess_entrant":
                continue
            dec_rows.append(dict(basin=b, term=k, est=v, lo=dlo.get(k, np.nan), hi=dhi.get(k, np.nan)))
        dec_rows.append(dict(basin=b, term="entrant_excess_over_baseline", est=dobs["excess_entrant"], lo=exc[0], hi=exc[1]))
        f = primary[b]["res"]["f"]
        dec_rows.append(dict(basin=b, term="f_share_of_log_RR_HF", est=f["est"], lo=f["lo"], hi=f["hi"]))
    pd.DataFrame(dec_rows).to_csv(os.path.join(outdir, "decomposition.csv"), index=False)
    json.dump(meta, open(os.path.join(outdir, "meta.json"), "w"), indent=1)
    print(f"primary written {time.time()-t0:.0f}s", flush=True)
    return primary, Tfull, M, allidx, rng


if __name__ == "__main__":
    main()
