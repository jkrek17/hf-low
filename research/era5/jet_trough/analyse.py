"""Jet speed and upstream trough depth against bomb and HF onset (ERA5 proxy, pipeline A).

Runs the analysis declared in PREREGISTRATION.md. Held-out seasons (2015-16..2025-26) are scored once
per declared model in this one run; `results/looks.log` records the run.

usage: analyse.py [--sens]      writes results/*.csv, results/jet_trough-result.txt, figures
"""
import sys, itertools, time, json, datetime
import numpy as np, pandas as pd
from multiprocessing import Pool
import jtlib as J
import os
OD = "work/dev" if J.DEV else "results"
os.makedirs(OD, exist_ok=True)

BASINS = ("atl", "pac")
OUTCOMES = ("BOMB", "HFON")
PRIMARY = ("jet250", "trough_up")
SECONDARY = ("jet_up", "trough_loc")
UNITS = {"jet250": "kt", "jet_up": "kt", "trough_up": "m", "trough_loc": "m", "trough_up_eddy": "m"}
NBOOT_KNOT = 5 if J.DEV else 200
OUT = []


def w(s=""):
    OUT.append(s)
    print(s, flush=True)


# ---------------------------------------------------------------- one threshold combination
def combo(args):
    basin, outcome, c, knotboot, tag = args
    D = J.load_all()
    D = filt(D, tag)
    S, y = J.subset(D, outcome, basin)
    if tag == "lead24":
        S, y = lead24(D, outcome, basin)
    cols = [] if tag == "unadjusted" else J.cov_cols(outcome)
    r = J.three_tests(S, y, c, cols, B=2000, seed=11, with_curves=True)
    o = r.pop("_obj")
    Dg, Sf, Se = o["Dg"], o["Sf"], o["Se"]
    lo, hi, mu, sd = Dg.xs[c]
    unit = lambda z: z * sd + mu
    # g-computation curves on the held-out rows' covariates (no outcomes used)
    Ce = Dg.cov(Se)
    zf = Dg.x(Sf, c)
    grid = np.percentile(zf, np.linspace(2, 98, 41))
    cur = []
    for g in grid:
        zz = np.full(len(Se), g)
        def pdp(kind):
            if kind == "lin":
                X = np.column_stack([Ce, zz]); b = o["bL"]
            elif kind == "hinge":
                X = np.column_stack([Ce, zz, np.maximum(zz - o["k"], 0)]); b = o["bh"]
            else:
                X = np.column_stack([Ce, Dg.rcs(c, zz)]); b = o["bS"]
            return float(J.prob(X, b).mean())
        cur.append(dict(basin=basin, outcome=outcome, var=c, tag=tag, x=unit(g), lin=pdp("lin"),
                        hinge=pdp("hinge"), spl=pdp("spl")))
    cur = pd.DataFrame(cur)
    # sharpness: odds ratio across the knot +/- 0.5 SD, and probabilities at fit percentiles
    def odds_ratio(kind, k):
        res = []
        for g in (k - 0.5, k + 0.5):
            zz = np.full(len(Se), g)
            if kind == "hinge":
                X = np.column_stack([Ce, zz, np.maximum(zz - o["k"], 0)]); b = o["bh"]
            else:
                X = np.column_stack([Ce, Dg.rcs(c, zz)]); b = o["bS"]
            p = J.prob(X, b).mean(); res.append(p / (1 - p))
        return res[1] / res[0]
    pct = {q: unit(np.percentile(zf, q)) for q in (10, 25, 50, 75, 90)}
    pp = {}
    for q in (10, 25, 50, 75, 90):
        i = int(np.argmin(np.abs(grid - np.percentile(zf, q))))
        pp[q] = (cur.hinge[i], cur.spl[i], cur.lin[i])
    row = dict(basin=basin, outcome=outcome, var=c, tag=tag, unit=UNITS[c],
               n_fit=r["n_fit"], n_test=r["n_test"], ev_fit=r["ev_fit"], ev_test=r["ev_test"],
               knot=unit(o["k"]), knot_pct=r["knot_pct"],
               slope_pre=r["slope_pre"], slope_post=r["slope_post"], b_lin=r["b_lin"],
               OR_knot_hinge=odds_ratio("hinge", o["k"]), OR_knot_spl=odds_ratio("spl", o["k"]),
               or_per_sd_lin=float(np.exp(r["b_lin"])))
    for t in ("T1", "T2", "T3", "T3rev"):
        for k2 in ("gain", "lo", "hi", "p"):
            row[f"{t}_{k2}"] = r[t][k2]
    for q in (10, 25, 50, 75, 90):
        row[f"x{q}"] = pct[q]
        row[f"P{q}_hinge"], row[f"P{q}_spl"], row[f"P{q}_lin"] = pp[q]
    # knot interval from fit-season resamples
    if knotboot:
        rng = np.random.default_rng(5)
        fs = np.unique(Sf.season.values)
        ks = []
        for _ in range(NBOOT_KNOT):
            pick = rng.choice(fs, len(fs))
            rows = np.concatenate([np.where(Sf.season.values == s)[0] for s in pick])
            k, *_ = J.best_hinge(Dg, Sf, o["yf"], c, rows=rows)
            ks.append(unit(k))
        row["knot_lo"], row["knot_hi"] = np.percentile(ks, [5, 95])
    # raw deciles, fit and held-out, with season-bootstrap intervals
    edges = np.percentile(Sf[c].values, np.linspace(0, 100, 11))
    edges[0], edges[-1] = -np.inf, np.inf
    dec_f = np.clip(np.searchsorted(edges, Sf[c].values, side="right") - 1, 0, 9)
    dec_e = np.clip(np.searchsorted(edges, Se[c].values, side="right") - 1, 0, 9)
    sea = Se.season.values; us = np.unique(sea)
    rng = np.random.default_rng(3)
    bs = np.zeros((500, 10))
    pre = {(s, d): (np.sum(o["ye"][(sea == s) & (dec_e == d)]), np.sum((sea == s) & (dec_e == d))) for s in us for d in range(10)}
    for b in range(500):
        pick = rng.choice(us, len(us))
        for d in range(10):
            a = sum(pre[(s, d)][0] for s in pick); n = sum(pre[(s, d)][1] for s in pick)
            bs[b, d] = a / max(n, 1)
    raw = pd.DataFrame(dict(basin=basin, outcome=outcome, var=c, tag=tag, decile=range(1, 11),
                            x_lo=np.where(np.isinf(edges[:-1]), np.nan, edges[:-1]),
                            n_fit=np.bincount(dec_f, minlength=10), freq_fit=[o["yf"][dec_f == d].mean() for d in range(10)],
                            n_test=np.bincount(dec_e, minlength=10), freq_test=[o["ye"][dec_e == d].mean() for d in range(10)],
                            lo_test=np.percentile(bs, 5, 0), hi_test=np.percentile(bs, 95, 0)))
    return row, cur, raw


# ---------------------------------------------------------------- variants
def filt(D, tag):
    if tag == "octapr":
        m = D.month.isin([10, 11, 12, 1, 2, 3, 4])
        return D[m].reset_index(drop=True)
    return D


def lead24(D, outcome, basin):
    """Predictors at t, outcome taken from the same track's fix 24 h later (window t+24..t+48 h)."""
    t = pd.to_datetime(D.time, format="%Y%m%d%H")
    K = D[["track", "time", "bomb", "elig", "hf_now", "hf24"]].copy()
    K["t"] = t
    L = K.copy(); L["t"] = L["t"] - pd.Timedelta(hours=24)
    M = D.assign(t=t).merge(L[["track", "t", "bomb", "elig", "hf_now", "hf24"]], on=["track", "t"], how="inner", suffixes=("", "_L"))
    M = M[M.basin == basin]
    if outcome == "BOMB":
        M = M[M.elig_L]; y = M.bomb_L.values
    else:
        M = M[(~M.hf_now) & (~M.hf_now_L)]; y = M.hf24_L.astype(int).values
    return M, y


# ---------------------------------------------------------------- interaction (H3)
def interaction(args):
    basin, outcome, pair = args
    D = J.load_all()
    S, y = J.subset(D, outcome, basin)
    cols = J.cov_cols(outcome)
    fm = (S.season <= J.FIT_LAST).values
    a, b2 = pair
    Dg = J.Design(S, fm, cols, [a, b2])
    C = Dg.cov(S)
    za, zb = Dg.x(S, a), Dg.x(S, b2)
    Xadd = np.column_stack([C, Dg.rcs(a, za), Dg.rcs(b2, zb)])
    Xint = np.column_stack([Xadd, za * zb])
    ba, _ = J.fit_logit(Xadd[fm], y[fm]); bi, _ = J.fit_logit(Xint[fm], y[fm])
    la = J.logloss(J.prob(Xadd[~fm], ba), y[~fm]); li = J.logloss(J.prob(Xint[~fm], bi), y[~fm])
    g = J.season_gain(la, li, S.season.values[~fm], 2000, 21)
    return dict(basin=basin, outcome=outcome, vars=f"{a}x{b2}", coef_int=float(bi[-1]),
                OR_int_per_sd2=float(np.exp(bi[-1])), gain=g["gain"], lo=g["lo"], hi=g["hi"], p=g["p"],
                ev_test=int(y[~fm].sum()))


# ---------------------------------------------------------------- grid (3 x 3 terciles)
def grid_tables():
    D = J.load_all()
    fm_all = D.fit.values
    # framework 'full' (L2, C=1) fitted on fit seasons, binary hf24 and binary rapid deepening
    out = []
    pred = {}
    for tgt, rows in (("hf24", D.index), ("bomb", D.index[D.elig])):
        S = D.loc[rows]
        fm = S.fit.values
        cols = J.FW.SETS["full"]
        prep = J.FW.Prep(S.loc[fm, cols].values.astype(float))
        m = J.FW.fit(prep(S.loc[fm, cols].values.astype(float)), S.loc[fm, tgt].astype(int).values)
        pred[tgt] = pd.Series(m.predict_proba(prep(S.loc[~fm, cols].values.astype(float)))[:, 1], index=S.index[~fm])
    for basin in BASINS:
        B = D[D.basin == basin]
        fit = B[B.fit]
        qj = np.percentile(fit.jet250, [100 / 3, 200 / 3]); qt = np.percentile(fit.trough_up, [100 / 3, 200 / 3])
        B = B.assign(jt=np.searchsorted(qj, B.jet250), tt=np.searchsorted(qt, B.trough_up))
        for per, sub in (("fit", B[B.fit]), ("test", B[~B.fit])):
            for (a, b), g in sub.groupby(["jt", "tt"]):
                e = g[g.elig]; h = g[~g.hf_now]
                row = dict(basin=basin, period=per, jet_tercile=a + 1, trough_tercile=b + 1,
                           jet_lo=(np.nan if a == 0 else qj[a - 1]), trough_lo=(np.nan if b == 0 else qt[b - 1]),
                           n=len(g), n_elig=len(e), bomb_freq=e.bomb.mean(), n_notHF=len(h), hfon_freq=h.hf24.mean(),
                           hf24_all_freq=g.hf24.mean(), both=((e.bomb == 1) & e.hf24).mean() if len(e) else np.nan,
                           bomb_only=((e.bomb == 1) & ~e.hf24).mean() if len(e) else np.nan,
                           hf_only=((e.bomb == 0) & e.hf24).mean() if len(e) else np.nan,
                           neither=((e.bomb == 0) & ~e.hf24).mean() if len(e) else np.nan)
                if per == "test":
                    row["fw_P_hf24"] = pred["hf24"].reindex(g.index).mean()
                    row["fw_P_bomb"] = pred["bomb"].reindex(e.index).mean()
                    row["hf24_obs_test"] = g.hf24.mean()
                    row["bomb_obs_test"] = e.bomb.mean() if len(e) else np.nan
                out.append(row)
    return pd.DataFrame(out)


# ---------------------------------------------------------------- added skill (H4)
def added_skill():
    D = J.load_all()
    rows_out, tests = [], []
    base = J.FW.SETS["full"]
    variants = {"full": base, "full+trough": base + ["trough_up"],
                "full+trough+loc+jetup": base + ["trough_up", "trough_loc", "jet_up"]}
    for tgt, kind in (("hf24", "HF24"), ("bomb", "BOMB")):
        S = D if tgt == "hf24" else D[D.elig]
        S = S.reset_index(drop=True)
        fm = S.fit.values
        y = (S[tgt].astype(int)).values
        pr = {}
        for name, cols in variants.items():
            prep = J.FW.Prep(S.loc[fm, cols].values.astype(float))
            m = J.FW.fit(prep(S.loc[fm, cols].values.astype(float)), y[fm])
            pr[name] = (m.predict_proba(prep(S.loc[~fm, cols].values.astype(float)))[:, 1],
                        m.predict_proba(prep(S.loc[fm, cols].values.astype(float)))[:, 1])
        # basin-month climatology from fit seasons
        fitS = S[fm].assign(y=y[fm])
        cl = fitS.groupby(["basin", "month"]).y.agg(["sum", "count"])
        te = S[~fm]
        pc = np.array([(cl.loc[(b, m_), "sum"] + 0.5) / (cl.loc[(b, m_), "count"] + 1) for b, m_ in zip(te.basin, te.month)])
        ye = y[~fm]; sea = te.season.values
        for scope in ("pooled", "atl", "pac"):
            mk = np.ones(len(te), bool) if scope == "pooled" else (te.basin.values == scope)
            bsc = {n: (p[0][mk] - ye[mk]) ** 2 for n, p in pr.items()}
            bsc["clim"] = (pc[mk] - ye[mk]) ** 2
            for name in variants:
                us = np.unique(sea[mk]); rng = np.random.default_rng(8)
                tot = lambda v, idx: np.array([v[sea[mk] == s].sum() for s in idx]).sum()
                per = {s: (bsc[name][sea[mk] == s].sum(), bsc["clim"][sea[mk] == s].sum(), bsc["full"][sea[mk] == s].sum()) for s in us}
                bss, d_full = [], []
                for _ in range(2000):
                    pick = rng.choice(us, len(us))
                    a = sum(per[s][0] for s in pick); c0 = sum(per[s][1] for s in pick); f = sum(per[s][2] for s in pick)
                    bss.append(1 - a / c0); d_full.append(f - a)
                # yes/no at the count-matched cut
                ptr = pr[name][1][(S.basin.values[fm] == scope) if scope != "pooled" else slice(None)]
                ytr = y[fm][(S.basin.values[fm] == scope) if scope != "pooled" else slice(None)]
                cut = np.quantile(ptr, 1 - ytr.mean())
                t = J.FW.table(pr[name][0][mk] >= cut, ye[mk].astype(bool))
                rows_out.append(dict(target=kind, scope=scope, model=name, n_test=int(mk.sum()), ev_test=int(ye[mk].sum()),
                                     BSS=1 - bsc[name].sum() / bsc["clim"].sum(),
                                     BSS_lo=np.percentile(bss, 5), BSS_hi=np.percentile(bss, 95),
                                     gain_brier_vs_full=bsc["full"].sum() - bsc[name].sum(),
                                     gain_lo=np.percentile(d_full, 5), gain_hi=np.percentile(d_full, 95),
                                     p_gain=(np.sum(np.array(d_full) <= 0) + 1) / 2001,
                                     POD=t["pod"], FAR=t["far"], HSS=t["hss"], bias=t["bias"]))
    return pd.DataFrame(rows_out)


# ---------------------------------------------------------------- main
def main(sens):
    t0 = time.time()
    combos = [(b, o, c, True, "main") for b in BASINS for o in OUTCOMES for c in PRIMARY]
    combos_sec = [(b, o, c, False, "secondary") for b in BASINS for o in OUTCOMES for c in SECONDARY]
    with Pool(4) as pool:
        res = pool.map(combo, combos + combos_sec)
        inter = pool.map(interaction, [(b, o, ("jet250", "trough_up")) for b in BASINS for o in OUTCOMES])
        sensres = []
        if sens:
            jobs = ([(b, o, c, False, "lead24") for b in BASINS for o in OUTCOMES for c in PRIMARY]
                    + [(b, o, c, False, "octapr") for b in BASINS for o in OUTCOMES for c in PRIMARY]
                    + [(b, o, c, False, "unadjusted") for b in BASINS for o in OUTCOMES for c in PRIMARY]
                    + [(b, o, "trough_up_eddy", False, "eddy") for b in BASINS for o in OUTCOMES])
            sensres = pool.map(combo, jobs)
    allr = res + sensres
    T = pd.DataFrame([r[0] for r in allr]); CUR = pd.concat([r[1] for r in allr]); RAW = pd.concat([r[2] for r in allr])
    T.to_csv(f"{OD}/thresholds.csv", index=False, float_format="%.5g")
    CUR.to_csv(f"{OD}/curves.csv", index=False, float_format="%.5g")
    RAW.to_csv(f"{OD}/raw_deciles.csv", index=False, float_format="%.5g")
    I = pd.DataFrame(inter)
    I.to_csv(f"{OD}/interaction.csv", index=False, float_format="%.5g")
    G = grid_tables(); G.to_csv(f"{OD}/grid_terciles.csv", index=False, float_format="%.5g")
    A = added_skill(); A.to_csv(f"{OD}/added_skill.csv", index=False, float_format="%.5g")

    # multiplicity
    main_ = T[T.tag == "main"].copy()
    prim = []
    for r in main_.itertuples():
        for t in ("T1", "T2", "T3"):
            prim.append(dict(family="primary", basin=r.basin, outcome=r.outcome, test=f"{r.var} {t}", p=getattr(r, f"{t}_p")))
    for r in I.itertuples():
        prim.append(dict(family="primary", basin=r.basin, outcome=r.outcome, test="jet x trough interaction", p=r.p))
    for r in A[(A.model == "full+trough") & (A.scope != "pooled")].itertuples():
        prim.append(dict(family="primary", basin=r.scope, outcome=r.target, test="trough added to framework", p=r.p_gain))
    P = pd.DataFrame(prim); P["q"] = J.bh(P.p.values)
    sec = T[T.tag == "secondary"]
    S2 = []
    for r in sec.itertuples():
        for t in ("T1", "T2", "T3"):
            S2.append(dict(family="secondary", basin=r.basin, outcome=r.outcome, test=f"{r.var} {t}", p=getattr(r, f"{t}_p")))
    S2 = pd.DataFrame(S2); S2["q"] = J.bh(S2.p.values)
    pd.concat([P, S2]).to_csv(f"{OD}/fdr_table.csv", index=False, float_format="%.5g")
    w(f"Jet speed / upstream trough depth: threshold or smooth ramp (ERA5 proxy, pipeline A, framework fixes)")
    w(f"Fit seasons 2004-05..2014-15, held out 2015-16..2025-26. Run {datetime.date.today()}; {time.time() - t0:.0f} s.")
    w()
    w("== Primary tests (held-out log-loss gains, season bootstrap, BH over all %d primary tests)" % len(P))
    w(P.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    w()
    w("== Secondary variables (separate BH family of %d)" % len(S2))
    w(S2.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    with open(f"{OD}/jet_trough-result.txt", "w") as f:
        f.write("\n".join(OUT) + "\n")
    open(f"{OD}/looks.log", "a").write(f"{datetime.datetime.utcnow().isoformat()}Z full analysis run; held-out seasons scored once per declared model (sens={sens})\n")


if __name__ == "__main__":
    main("--sens" in sys.argv)
