"""RA-7: why HF onset comes at or after the minimum-pressure fix in ~29% of storms.
ERA5 PROXY, pipeline A (research/era5/hf_history), seasons 2004-05..2025-26, non-TC events.
Plan: PREREGISTRATION.md (committed before this file, 33039f0). Reads committed files only.

usage: python3 late_wind.py        -> results/late_wind.txt, results/tests.csv, results/joint.csv, results/late_wind.png
"""
import os, sys
import numpy as np, pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ERA = os.path.abspath(os.path.join(HERE, ".."))
RES = os.path.join(HERE, "results")
SEED = 7
NBOOT = 2000
NPOW = 2000
INFO = {}


# ---------- data ----------
def load():
    L = pd.read_csv(os.path.join(ERA, "hf_history/results/lifecycle_events.csv"))
    M = L[(L.era == "2004+") & (~L.tc.astype(bool))].copy()
    f = pd.read_csv(os.path.join(ERA, "hf_structure/results/fixes.csv")).sort_values(["track", "time"])
    fo = f.groupby("track").first().reset_index()
    fo["t_on_f"] = pd.to_datetime(fo.time.astype(str), format="%Y%m%d%H")
    C = pd.read_csv(os.path.join(ERA, "hf_history/results/era5_hf_catalog.csv"))[["track", "gust800_kt"]]
    M = M.merge(fo[["track", "t_on_f", "gmax_r", "gmax_rel", "gmax_lat", "terrain", "speed_kt", "time", "msl"]], on="track", how="left")
    M = M.merge(C, on="track", how="left")
    T = pd.read_csv(os.path.join(ERA, "hf_history/results/era5_hf_catalog_tracks.csv"))
    T = T[T.track.isin(M.track) & T.basin.notna() & (T.g800 >= 71.7)]
    T["t"] = pd.to_datetime(T.time.astype(str), format="%Y%m%d%H")
    M = M.merge(T.groupby("track").t.min().rename("t_on").reset_index(), on="track", how="left")
    n0 = len(M)
    M = M.reset_index(drop=True)
    bad = M[(M.t_on_f != M.t_on) | M.gmax_r.isna()]
    M = M.drop(bad.index)
    INFO["dropped"] = n0 - len(M)
    M["late"] = (M.h_on_minp >= 0).astype(int)
    M["late6"] = (M.h_on_minp >= 6).astype(int)
    M["late12"] = (M.h_on_minp >= 12).astype(int)
    M["pac"] = (M.basin == "pac").astype(float)
    mo = M.mon_on.astype(int)
    M["mg"] = np.where(mo.isin([12, 1, 2]), "DJF", np.where(mo.isin([3, 4, 5]), "MAM", "SON"))   # Jun-Aug folded into SON
    M["T1"] = (M.gmax_r > 400).astype(float)
    M["T2"] = (M.gmax_lat >= 60).astype(float)
    M["T3"] = (M.gust800_kt < 76.7).astype(float)
    M["T4"] = ((M.gmax_rel >= 90) & (M.gmax_rel < 270)).astype(float)
    M["T5"] = np.where(M.maxdeep.isna(), np.nan, (M.maxdeep >= 1).astype(float))
    M["T6"] = (M.speed_kt - M.speed_kt.mean()) / M.speed_kt.std()
    M["T7"] = M.pac
    # environment at the last 00/12 fix at or before onset, at most 6 h earlier
    E = pd.read_csv(os.path.join(ERA, "intensity/results/env_2004.csv.gz"))
    E["te"] = pd.to_datetime(E.time.astype(str), format="%Y%m%d%H")
    E = E.sort_values("te")
    M = M.sort_values("t_on")
    M = pd.merge_asof(M, E.drop(columns=["time"]), left_on="t_on", right_on="te", by="track",
                      direction="backward", tolerance=pd.Timedelta(hours=6))
    ENV = ["B", "VTL", "VTU", "jet250", "div300", "vadv500", "eady", "flux", "tcwv"]
    for i, c in enumerate(ENV, 1):
        M[f"S{i}"] = M.groupby("basin")[c].transform(lambda s: (s - s.mean()) / s.std())
    return M.sort_values("track").reset_index(drop=True)


# ---------- logistic + cluster-robust inference ----------
def design(d, xcol, basin_term=True):
    cols = [np.ones(len(d))]
    if basin_term:
        cols.append(d.pac.values)
    cols += [(d.mg == "DJF").values.astype(float), (d.mg == "MAM").values.astype(float)]
    cols.append(d[xcol].values.astype(float))
    return np.column_stack(cols)


def irls(X, y, it=30):
    b = np.zeros(X.shape[1])
    for _ in range(it):
        eta = X @ b
        p = 1 / (1 + np.exp(-eta))
        w = np.clip(p * (1 - p), 1e-9, None)
        z = eta + (y - p) / w
        XtW = X.T * w
        try:
            bn = np.linalg.solve(XtW @ X, XtW @ z)
        except np.linalg.LinAlgError:
            return b, np.nan
        if np.max(np.abs(bn - b)) < 1e-9:
            b = bn
            break
        b = bn
    return b, 0


def fit(d, xcol, y="late", basin_term=True, wald=True):
    d = d.dropna(subset=[xcol])
    X = design(d, xcol, basin_term)
    yy = d[y].values.astype(float)
    b, bad = irls(X, yy)
    out = dict(b=b[-1], n=len(d))
    if not wald:
        return out
    p = 1 / (1 + np.exp(-(X @ b)))
    w = p * (1 - p)
    A = np.linalg.inv((X.T * w) @ X)
    g = d.season.values
    G = np.unique(g)
    meat = np.zeros((X.shape[1], X.shape[1]))
    s = X * (yy - p)[:, None]
    for k in G:
        u = s[g == k].sum(axis=0)
        meat += np.outer(u, u)
    n, kk = X.shape
    c = len(G) / (len(G) - 1) * (n - 1) / (n - kk)
    V = c * A @ meat @ A
    se = np.sqrt(V[-1, -1])
    out.update(se=se, t=b[-1] / se, p=2 * stats.t.sf(abs(b[-1] / se), len(G) - 1))
    tc = stats.t.ppf(0.975, len(G) - 1)
    out.update(lo=b[-1] - tc * se, hi=b[-1] + tc * se)
    out["fit"] = (X, b)
    return out


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    m = len(p)
    q = np.empty(m)
    r = p[o] * m / (np.arange(m) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    q[o] = np.minimum(r, 1)
    return q


def boot_ci(d, xcol, basin_term, rng):
    seas = d.season.unique()
    groups = {s: d[d.season == s] for s in seas}
    bs = []
    for _ in range(NBOOT):
        pick = rng.choice(seas, len(seas))
        dd = pd.concat([groups[s] for s in pick], ignore_index=True)
        if dd[xcol].nunique() < 2:
            continue
        bs.append(fit(dd, xcol, basin_term=basin_term, wald=False)["b"])
    return np.percentile(bs, [2.5, 97.5])


def std_shares(d, xcol, basin_term, res):
    """basin-month standardised late share with the factor set to 1 and to 0 for every row."""
    X, b = res["fit"]
    X1, X0 = X.copy(), X.copy()
    X1[:, -1], X0[:, -1] = 1, 0
    f = lambda Z: float(np.mean(1 / (1 + np.exp(-(Z @ b)))))
    return f(X1), f(X0)


TESTS = [  # id, label, expected sign, basin_term, subset
    ("T1", "gust maximum > 400 km from centre", "+", True, None),
    ("T2", "Atlantic: gust maximum at or north of 60N", "+", False, "atl"),
    ("T3", "marginal storm (peak gust < 76.7 kt)", "+", True, None),
    ("T4", "gust maximum in rear half", "+", True, None),
    ("T5", "explosive deepening (>= 1 Bergeron)", "-", True, None),
    ("T6", "translation speed, per SD", "2s", True, None),
    ("T7", "Pacific (vs Atlantic)", "-", False, None),
]
ENVN = ["B", "VTL", "VTU", "jet250", "div300", "vadv500", "eady", "flux", "tcwv"]


def run_family(M, tests, rng, boot=True, y="late"):
    rows = []
    for tid, lab, sign, bt, sub in tests:
        d = M[M.basin == sub] if sub else M
        r = fit(d, tid, y=y, basin_term=bt)
        row = dict(id=tid, label=lab, expected=sign, n=r["n"], OR=np.exp(r["b"]), lo=np.exp(r["lo"]), hi=np.exp(r["hi"]),
                   p=r["p"], b=r["b"])
        if y == "late":
            if bt or tid == "T7":
                pass
            sh1, sh0 = std_shares(d, tid, bt, r) if tid not in ("T6",) and not tid.startswith("S") else (np.nan, np.nan)
            row.update(share_in=sh1, share_out=sh0)
            if tid in ("T1", "T2", "T3", "T4", "T5", "T7"):
                x = d.dropna(subset=[tid])
                row["frac_of_late"] = float(x[(x.late == 1)][tid].mean())
                row["prevalence"] = float(x[tid].mean())
        if boot:
            lo, hi = boot_ci(d, tid, bt, rng)
            row.update(blo=np.exp(lo), bhi=np.exp(hi))
        rows.append(row)
    return pd.DataFrame(rows)


def power(M, tests, rng):
    out = []
    for tid, lab, sign, bt, sub in tests:
        d = (M[M.basin == sub] if sub else M).dropna(subset=[tid]).reset_index(drop=True)
        X = design(d, tid, bt)
        b0, _ = irls(X[:, :-1], d.late.values.astype(float))
        eta0 = X[:, :-1] @ b0
        seas = d.season.values
        # centre the planted predictor so the base rate stays put
        xs = d[tid].values.astype(float)
        res = {}
        for OR in (1.5, 2.0):
            hit = 0
            for _ in range(NPOW):
                p = 1 / (1 + np.exp(-(eta0 + np.log(OR) * (xs - xs.mean()))))
                yy = (rng.random(len(d)) < p).astype(int)
                dd = d.assign(late=yy)
                r = fit(dd, tid, basin_term=bt)
                hit += (r["p"] < 0.05)
            res[OR] = hit / NPOW
        # minimum detectable OR at 80%, 400 draws per grid point
        mde = np.nan
        for OR in np.arange(1.1, 3.01, 0.1):
            hit = 0
            for _ in range(400):
                p = 1 / (1 + np.exp(-(eta0 + np.log(OR) * (xs - xs.mean()))))
                yy = (rng.random(len(d)) < p).astype(int)
                hit += fit(d.assign(late=yy), tid, basin_term=bt)["p"] < 0.05
            if hit / 400 >= 0.8:
                mde = round(float(OR), 1)
                break
        out.append(dict(id=tid, power_OR1_5=res[1.5], power_OR2_0=res[2.0], mde_OR_80=mde))
    return pd.DataFrame(out)


def auc(y, s):
    r = stats.rankdata(s)
    n1 = y.sum()
    n0 = len(y) - n1
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def joint(M, preds, basin_only=None):
    d = M if basin_only is None else M[M.basin == basin_only]
    d = d.dropna(subset=preds).reset_index(drop=True)
    base = ["pac"] if basin_only is None else []
    def X(dd, cols):
        c = [np.ones(len(dd))] + [dd[k].values.astype(float) for k in cols] + [(dd.mg == "DJF").values.astype(float), (dd.mg == "MAM").values.astype(float)]
        return np.column_stack(c)
    pj, pb = np.zeros(len(d)), np.zeros(len(d))
    for s in d.season.unique():
        tr, te = d.season != s, d.season == s
        bj, _ = irls(X(d[tr], base + preds), d.late[tr].values.astype(float))
        bb, _ = irls(X(d[tr], base), d.late[tr].values.astype(float))
        pj[te.values] = 1 / (1 + np.exp(-(X(d[te], base + preds) @ bj)))
        pb[te.values] = 1 / (1 + np.exp(-(X(d[te], base) @ bb)))
    y = d.late.values
    bsj, bsb = np.mean((pj - y) ** 2), np.mean((pb - y) ** 2)
    bfull, _ = irls(X(d, base + preds), y.astype(float))
    return dict(n=len(d), auc=auc(y, pj), auc_base=auc(y, pb), bss=1 - bsj / bsb,
                coef={k: float(np.exp(v)) for k, v in zip(["const"] + base + preds + ["DJF", "MAM"], bfull)})


def main():
    rng = np.random.default_rng(SEED)
    M = load()
    ENV = ["B", "VTL", "VTU", "jet250", "div300", "vadv500", "eady", "flux", "tcwv"]
    envtests = [(f"S{i}", f"{n} at onset, per SD", "2s", True, None) for i, n in enumerate(ENV, 1)]
    out = []
    w = out.append
    w("RA-7  Why HF onset comes at or after the minimum-pressure fix (ERA5 PROXY, pipeline A, non-TC events)")
    w("Plan: PREREGISTRATION.md (33039f0). Seasons 2004-05..2025-26. Tests: logistic, basin + month-group terms, season-clustered CR1 Wald, t on 21 df.")
    w(f"Events used {len(M)} (Atl {int((M.basin=='atl').sum())}, Pac {int((M.basin=='pac').sum())}); dropped for onset-time mismatch or missing structure: {INFO['dropped']}")
    w(f"Environment matched (00/12 UTC fix within 6 h before onset): {int(M.B.notna().sum())} of {len(M)}")
    w(f"LATE (onset at/after min pressure) share: all {M.late.mean():.3f}, Atl {M[M.basin=='atl'].late.mean():.3f}, Pac {M[M.basin=='pac'].late.mean():.3f}; "
      f"LATE6 {M.late6.mean():.3f}; LATE12 {M.late12.mean():.3f}")
    w("")
    # ---- primary ----
    T = run_family(M, TESTS, rng)
    T["q_family"] = bh(T.p)
    S = run_family(M, envtests, rng)
    S["q_family"] = bh(S.p)
    allp = np.concatenate([T.p.values, S.p.values])
    qa = bh(allp)
    T["q_all"], S["q_all"] = qa[:len(T)], qa[len(T):]
    P = power(M, TESTS, rng)
    T = T.merge(P, on="id")
    pd.concat([T, S]).drop(columns=["b"]).to_csv(os.path.join(RES, "tests.csv"), index=False)

    w("== Primary family T1-T7 (OR for LATE; share_in / share_out = basin-month standardised late share with / without the factor) ==")
    for _, r in T.iterrows():
        sh = f"  late share {r.share_in:.3f} vs {r.share_out:.3f}" if pd.notna(r.get("share_in")) else ""
        fr = f"  prevalence {r.prevalence:.3f}, share of late storms with factor {r.frac_of_late:.3f}" if pd.notna(r.get("prevalence")) else ""
        w(f"{r.id} {r.label} (expected {r.expected}) n={int(r.n)}: OR {r.OR:.2f} [{r.lo:.2f}, {r.hi:.2f}] boot [{r.blo:.2f}, {r.bhi:.2f}] "
          f"p={r.p:.4f} q={r.q_family:.4f} q_all16={r.q_all:.4f}{sh}{fr}")
        w(f"      power at OR 1.5 / 2.0 = {r.power_OR1_5:.2f} / {r.power_OR2_0:.2f}; minimum detectable OR at 80% = {r.mde_OR_80}")
    w(f"Passing q<0.05 in family: {int((T.q_family<0.05).sum())} of {len(T)}")
    w("")
    w("== Exploratory family S1-S9 (environment at onset, per SD, two-sided) ==")
    for _, r in S.iterrows():
        w(f"{r.id} {r.label}: n={int(r.n)} OR {r.OR:.2f} [{r.lo:.2f}, {r.hi:.2f}] p={r.p:.4f} q={r.q_family:.4f} q_all16={r.q_all:.4f}")
    w(f"Passing q<0.05 in family: {int((S.q_family<0.05).sum())} of {len(S)}; across all 16: {int((qa<0.05).sum())} of 16")
    w("")
    # ---- joint ----
    jp = ["T1", "T3", "T4", "T5", "T6"]
    J = joint(M, jp)
    w("== How much is explained: joint logistic of T1, T3, T4, T5, T6 + basin + month group, leave-one-season-out ==")
    w(f"n={J['n']}: AUC {J['auc']:.3f} (basin+month only: {J['auc_base']:.3f}); Brier skill vs basin-month base rate {J['bss']:.4f}")
    w("full-sample odds ratios: " + ", ".join(f"{k} {v:.2f}" for k, v in J["coef"].items()))
    JA = joint(M, jp + ["T2"], "atl")
    w(f"Atlantic-only with T2 added (T1,T3,T4,T5,T6,T2): n={JA['n']}: AUC {JA['auc']:.3f} (month only: {JA['auc_base']:.3f}); BSS {JA['bss']:.4f}")
    w("full-sample odds ratios: " + ", ".join(f"{k} {v:.2f}" for k, v in JA["coef"].items()))
    pd.DataFrame([dict(model="joint_all", **{k: v for k, v in J.items() if k != "coef"}),
                  dict(model="joint_atl_T2", **{k: v for k, v in JA.items() if k != "coef"})]).to_csv(os.path.join(RES, "joint.csv"), index=False)
    w("")
    # ---- descriptive ----
    w("== Descriptive: where the 29% comes from ==")
    M["off_pg"] = M.h_on_minp - M.h_pg_minp       # onset minus peak-gust time (<= 0)
    bins = [-1e9, -13, -7, -1, 0, 6, 12, 1e9]
    labs = ["<=-18", "-12", "-6", "0", "+6", "+12", ">=+18"]
    for b in ("atl", "pac", "all"):
        s = M if b == "all" else M[M.basin == b]
        c = pd.cut(s.h_on_minp, [-1e9, -13, -7, -1, 0, 6, 12, 1e9] if False else [-1e9, -15, -9, -3, 3, 9, 15, 1e9], labels=["<=-18", "-12", "-6", "0", "+6", "+12", ">=+18"])
        w(f"{b}: onset minus min-pressure time, h, shares " + ", ".join(f"{k} {v:.3f}" for k, v in c.value_counts(normalize=True).reindex(["<=-18", "-12", "-6", "0", "+6", "+12", ">=+18"]).items()))
    late = M[M.late == 1]
    w(f"late storms n={len(late)}: single HF fix {np.mean(late.n_hf==1):.3f} (early storms {np.mean(M[M.late==0].n_hf==1):.3f}); onset == peak-gust fix {np.mean(late.off_pg==0):.3f}; "
      f"peak gust after min pressure by >=6 h {np.mean(late.h_pg_minp>=6):.3f}; late but LATE6 false (onset at the minimum fix) {np.mean(late.late6==0):.3f}")
    w("peak-gust time minus min-pressure time, h (all events): " + ", ".join(f"{q}% {np.percentile(M.h_pg_minp, q):.0f}" for q in (5, 25, 50, 75, 95)))
    w("late share by number of HF fixes: " + ", ".join(f"{k}: {v:.3f} (n={n})" for (k, v), n in zip(M.groupby(pd.cut(M.n_hf, [0, 1, 2, 4, 100], labels=["1", "2", "3-4", "5+"]), observed=True).late.mean().items(), M.groupby(pd.cut(M.n_hf, [0, 1, 2, 4, 100], labels=["1", "2", "3-4", "5+"]), observed=True).size())))
    w("late share by month group: " + ", ".join(f"{k} {v:.3f} (n={n})" for (k, v), n in zip(M.groupby("mg").late.mean().items(), M.groupby("mg").size())))
    w("")
    # ---- heterogeneity ----
    w("== Heterogeneity (OR within subset; sign agreement only, not tests) ==")
    for tid, lab, sign, bt, sub in TESTS[:5]:
        cells = []
        for name, d in (("Atl", M[M.basin == "atl"]), ("Pac", M[M.basin == "pac"]), ("2004-14", M[M.season <= 2014]), ("2015-25", M[M.season >= 2015])):
            if tid == "T2" and name == "Pac":
                continue
            r = fit(d, tid, basin_term=(name not in ("Atl", "Pac")))
            cells.append(f"{name} {np.exp(r['b']):.2f} [{np.exp(r['lo']):.2f}, {np.exp(r['hi']):.2f}]")
        w(f"{tid}: " + "; ".join(cells))
    w("")
    # ---- sensitivity ----
    w("== Sensitivity (not counted in the families) ==")
    for y in ("late6", "late12"):
        R = run_family(M, TESTS, rng, boot=False, y=y)
        w(f"outcome {y} (share {M[y].mean():.3f}): " + "; ".join(f"{r.id} OR {r.OR:.2f} p={r.p:.3f}" for _, r in R.iterrows()))
    subs = {"events with >=2 HF fixes": M[M.n_hf >= 2],
            "Atlantic gust max north of 60N removed": M[~((M.basin == "atl") & (M.T2 == 1))],
            "Atlantic terrain-flagged removed": M[~((M.basin == "atl") & (M.terrain.astype(bool)))]}
    for k, d in subs.items():
        tt = [t for t in TESTS if t[0] != "T2"] if "north of 60N" in k else TESTS
        R = run_family(d, tt, rng, boot=False)
        w(f"{k} (n={len(d)}, late share {d.late.mean():.3f}): " + "; ".join(f"{r.id} OR {r.OR:.2f} p={r.p:.3f}" for _, r in R.iterrows()))
    w("")
    w("Not pre-registered changes: Jun-Aug onsets (rare) are folded into the Sep-Nov month group; see PREREGISTRATION.md deviations.")
    open(os.path.join(RES, "late_wind.txt"), "w").write("\n".join(out) + "\n")
    print("\n".join(out))


if __name__ == "__main__":
    main()
