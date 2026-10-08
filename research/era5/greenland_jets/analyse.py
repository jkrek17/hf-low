"""Which ingredients give more Greenland-jet gusts? The pre-registered analysis.

ERA5 proxy (ARCO-ERA5 gust and surface fields; WeatherBench2 Z500 for the GBI
check); lows are pipeline A's. Everything follows PREREGISTRATION.md. Run:

    ERA5_WORK=... python3 analyse.py [NBOOT] [--swap]

Reads results/stage1_times.csv, results/sample_times.csv, $ERA5_WORK/greenland_jets/stage2/*.csv,
the CPC daily NAO file ($CPC_DIR, default /mnt/project-files/teleconnection-test/cpc_indices) and
the intensity fix table. Writes results/results.txt (or results_swap.txt) and results/tests.csv.

Notes on points the plan left open, fixed before any model was fitted (README, "Clarifications"):
  * the fix table's speed is kt over the previous 6 h, and is missing at a track's first fix;
    those rows keep zero standardised motion and get an indicator (mmiss) in the baseline
  * z-scores use the discovery half's mean and SD weighted by the inverse sampling weight, so
    "per SD" is per SD of all low-in-region times
  * logistic fits use a ridge of 1e-6 on all terms for numerical stability only
"""
import os, sys, glob, time
import numpy as np, pandas as pd
from scipy.special import expit
from scipy.stats import chi2
from sklearn.metrics import roc_auc_score

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(os.environ.get("ERA5_WORK", os.path.join(HERE, "work")), "greenland_jets")
RES = os.path.join(HERE, "results")
CPC = os.environ.get("CPC_DIR", "/mnt/project-files/teleconnection-test/cpc_indices")
FIX = os.path.join(HERE, "..", "intensity", "results", "fixes_2004.csv.gz")
NBOOT = int(next((a for a in sys.argv[1:] if a.isdigit()), 2000))
SWAP = "--swap" in sys.argv
SEED = 20261008
KM = 111.195
CF = (60.0, -44.0)                      # Cape Farewell
out, rows = [], []


def say(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)


# ------------------------------------------------------------------ data
def load():
    s = pd.read_csv(f"{RES}/sample_times.csv")
    p2 = f"{RES}/stage2_ingredients.csv"                 # committed copy of stage2.py's per-time output
    s2 = pd.read_csv(p2) if os.path.exists(p2) else pd.concat([pd.read_csv(p) for p in glob.glob(f"{WORK}/stage2/*.csv")])
    d = s.merge(s2[["time", "STAB", "GBI"]], on="time")
    assert len(d) == len(s)
    nao = pd.read_csv(f"{CPC}/norm.daily.nao.index.b500101.current.ascii", sep=r"\s+|(?<=\d)(?=-99)", header=None,
                      names=["y", "m", "d", "NAO"], engine="python")
    nao.loc[nao.NAO <= -90, "NAO"] = np.nan             # -99.000 is missing
    nao["date"] = pd.to_datetime(dict(year=nao.y, month=nao.m, day=nao.d))
    t = pd.to_datetime(d.time.astype(str), format="%Y%m%d%H")
    d["month"] = t.dt.month
    d["date"] = t.dt.normalize()
    d = d.merge(nao[["date", "NAO"]], on="date", how="left")
    d["x"] = (d.ref_lon - CF[1]) * np.cos(np.radians(d.ref_lat)) * KM
    d["y"] = (d.ref_lat - CF[0]) * KM
    d["Dcf"] = np.hypot(d.x, d.y)
    h = np.radians(d.ref_heading)
    d["u"] = d.ref_speed * np.sin(h)
    d["v"] = d.ref_speed * np.cos(h)
    d["mmiss"] = d.ref_speed.isna().astype(float)
    return d.reset_index(drop=True)


def wstats(d, cols):
    st = {}
    for c in cols:
        ok = d[c].notna()
        w = d.w[ok]
        m = np.average(d[c][ok], weights=w)
        sd = np.sqrt(np.average((d[c][ok] - m) ** 2, weights=w))
        st[c] = (m, sd)
    return st


ZC = ["GH", "GRAD", "STAB", "NAO", "GBI", "u", "v", "ref_msl", "x", "y", "Dcf"]


def prep(d, st):
    z = d.copy()
    for c in ZC:
        m, sd = st[c]
        z[c] = (d[c] - m) / sd
    z["u"] = z.u.fillna(0.0)
    z["v"] = z.v.fillna(0.0)
    z["depth"] = z.ref_msl
    z["trend"] = (d.season - 2015) / 6.5
    z["xx"], z["yy"], z["xy"] = z.x ** 2, z.y ** 2, z.x * z.y
    for mth, nm in ((12, "m12"), (1, "m1"), (2, "m2"), (3, "m3")):
        z[nm] = (d.month == mth).astype(float)
    z["X1"] = z.GH * z.depth
    z["X2"] = z.GH * z.STAB
    z["X3"] = z.GH * z.Dcf
    z["X4"] = z.depth * z.STAB
    return z


BASE = ["m12", "m1", "m2", "m3", "depth", "x", "y", "xx", "yy", "xy", "trend", "mmiss"]
ING = ["GH", "GRAD", "STAB", "NAO"]
MOT = ["u", "v"]
XINT = ["X1", "X2", "X3", "X4"]
XDESC = {"X1": "GH x depth", "X2": "GH x STAB", "X3": "GH x distance from Cape Farewell", "X4": "depth x STAB"}


def mat(z, cols):
    return np.column_stack([np.ones(len(z))] + [z[c].values for c in cols])


# ------------------------------------------------------------------ fitting
def logit(X, y, l2=1e-6, it=40):
    b = np.zeros(X.shape[1])
    for _ in range(it):
        eta = X @ b
        p = expit(eta)
        W = np.clip(p * (1 - p), 1e-8, None)
        H = (X * W[:, None]).T @ X + l2 * np.eye(X.shape[1])
        g = X.T @ (y - p) - l2 * b
        step = np.linalg.solve(H, g)
        b = b + step
        if np.abs(step).max() < 1e-8:
            break
    return b


def wls(X, y, w):
    A = (X * w[:, None]).T @ X
    return np.linalg.solve(A + 1e-9 * np.eye(X.shape[1]), (X * w[:, None]).T @ y)


def seasons_idx(z):
    return [np.where(z.season.values == s)[0] for s in sorted(z.season.unique())]


def boot_idx(sidx, rng):
    k = len(sidx)
    return np.concatenate([sidx[i] for i in rng.integers(0, k, k)])


# ------------------------------------------------------------------ F1 style tests
def estimates(z, cols_main, y, kind="logit"):
    X = mat(z, cols_main)
    return (logit(X, y) if kind == "logit" else wls(X, z.G_T.values, z.w.values))


def run_tests(z, kind, label, nboot, seed, sign=None):
    """The 5 main-effect tests (GH, GRAD, STAB, NAO, MOT joint) and, for logit, X1..X4."""
    y = z.case.values.astype(float)
    sidx = seasons_idx(z)
    main = BASE + ING + MOT
    models = {"main": main}
    if kind == "logit":
        for k in XINT:
            models[k] = main + [k]
    names = {k: [c for c in v] for k, v in models.items()}
    est = {}
    for k, cols in models.items():
        est[k] = estimates(z, cols, y, kind)
    rng = np.random.default_rng(seed)
    bs = {k: np.zeros((nboot, len(models[k]) + 1)) for k in models}
    zz = z.reset_index(drop=True)
    for b in range(nboot):
        ix = boot_idx(sidx, rng)
        zb = zz.iloc[ix]
        yb = y[ix]
        for k, cols in models.items():
            bs[k][b] = estimates(zb, cols, yb, kind)
    res = []
    for c in ING:
        j = models["main"].index(c) + 1
        b0, bb = est["main"][j], bs["main"][:, j]
        res.append(dict(test=c, est=b0, lo=np.percentile(bb, 2.5), hi=np.percentile(bb, 97.5),
                        p=2 * min((bb >= 0).mean(), (bb <= 0).mean()) + 1 / nboot,
                        p1=(bb <= 0).mean() + 1 / nboot if b0 > 0 else (bb >= 0).mean() + 1 / nboot,
                        sd_boot=bb.std(ddof=1)))
    ju, jv = models["main"].index("u") + 1, models["main"].index("v") + 1
    bm = est["main"][[ju, jv]]
    Bm = bs["main"][:, [ju, jv]]
    S = np.cov(Bm.T)
    Si = np.linalg.inv(S)
    obs = float(bm @ Si @ bm)
    cen = Bm - Bm.mean(0)
    dist = np.einsum("ij,jk,ik->i", cen, Si, cen)
    res.append(dict(test="MOT", est=np.nan, lo=np.nan, hi=np.nan, p=(dist >= obs).mean() + 1 / nboot,
                    p1=np.nan, sd_boot=np.nan, wald=obs, bu=bm[0], bv=bm[1]))
    if kind == "logit":
        for k in XINT:
            j = models[k].index(k) + 1
            b0, bb = est[k][j], bs[k][:, j]
            res.append(dict(test=k, est=b0, lo=np.percentile(bb, 2.5), hi=np.percentile(bb, 97.5),
                            p=2 * min((bb >= 0).mean(), (bb <= 0).mean()) + 1 / nboot,
                            p1=(bb <= 0).mean() + 1 / nboot if b0 > 0 else (bb >= 0).mean() + 1 / nboot,
                            sd_boot=bb.std(ddof=1)))
    r = pd.DataFrame(res)
    r["q"] = bh(r.p.values)
    r["label"] = label
    return r, est["main"], models["main"], bs["main"]


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    n = len(p)
    q = np.empty(n)
    cm = 1.0
    for rank in range(n, 0, -1):
        i = o[rank - 1]
        cm = min(cm, p[i] * n / rank)
        q[i] = cm
    return q


# ------------------------------------------------------------------ skill
def brier(y, p, w):
    return float(np.average((y - p) ** 2, weights=w))


def skill(zd, zh, nboot, seed):
    """Fit on zd, predict zh; weights w on zh represent all low-in-region times."""
    yd = zd.case.values.astype(float)
    yh = zh.case.values.astype(float)
    f = 1.0 / np.average(zd.w[zd.case == 0])           # sampling fraction of controls in the fit half
    off = np.log(f)
    cols0, cols1 = BASE, BASE + ING + MOT
    b0 = logit(mat(zd, cols0), yd)
    b1 = logit(mat(zd, cols1), yd)
    p0 = expit(mat(zh, cols0) @ b0 + off)
    p1 = expit(mat(zh, cols1) @ b1 + off)
    cm = np.average(zd.case, weights=zd.w)
    mon = {m: np.average(zd.case[zd.month == m], weights=zd.w[zd.month == m]) for m in sorted(zd.month.unique())}
    pc = zh.month.map(mon).values
    w = zh.w.values
    sidx = seasons_idx(zh.reset_index(drop=True))
    rng = np.random.default_rng(seed)
    stat = lambda ix: dict(
        bss01=1 - brier(yh[ix], p1[ix], w[ix]) / brier(yh[ix], p0[ix], w[ix]),
        bss1c=1 - brier(yh[ix], p1[ix], w[ix]) / brier(yh[ix], pc[ix], w[ix]),
        bss0c=1 - brier(yh[ix], p0[ix], w[ix]) / brier(yh[ix], pc[ix], w[ix]),
        auc1=roc_auc_score(yh[ix], p1[ix], sample_weight=w[ix]),
        auc0=roc_auc_score(yh[ix], p0[ix], sample_weight=w[ix]),
        brier1=brier(yh[ix], p1[ix], w[ix]), brier0=brier(yh[ix], p0[ix], w[ix]), brierc=brier(yh[ix], pc[ix], w[ix]))
    allix = np.arange(len(zh))
    obs = stat(allix)
    bs = pd.DataFrame([stat(boot_idx(sidx, rng)) for _ in range(nboot)])
    ci = {k: (np.percentile(bs[k], 2.5), np.percentile(bs[k], 97.5)) for k in bs}
    return obs, ci, b1


def g_comp(z, b, cols, name, lo=0.10, hi=0.90, f=None):
    """Average predicted probability with one ingredient set to its 10th / 90th weighted percentile."""
    w = z.w.values
    v = z[name].values
    order = np.argsort(v)
    cw = np.cumsum(w[order]) / w.sum()
    q = lambda p: v[order][np.searchsorted(cw, p)]
    out = []
    for val in (q(lo), q(hi)):
        zz = z.copy()
        zz[name] = val
        for k, comp in (("X1", ("GH", "depth")), ("X2", ("GH", "STAB")), ("X3", ("GH", "Dcf")), ("X4", ("depth", "STAB"))):
            zz[k] = zz[comp[0]] * zz[comp[1]]
        eta = mat(zz, cols) @ b + np.log(f)
        out.append(float(np.average(expit(eta), weights=w)))
    return q(lo), q(hi), out[0], out[1]


# ------------------------------------------------------------------ main
def main():
    d = load()
    fit_half, test_half = ("heldout", "discovery") if SWAP else ("discovery", "heldout")
    tag = "SWAP (replication: fit on odd start years, test on even)" if SWAP else "PRIMARY (fit on even start years, test on odd)"
    say(f"Greenland jets analysis, {tag}; ERA5 proxy; outcome G_T (not pipeline A's index); {NBOOT} season-block bootstraps")
    dd = d[d.half == fit_half].reset_index(drop=True)
    dh = d[d.half == test_half].reset_index(drop=True)
    st = wstats(dd, ZC)
    zd, zh = prep(dd, st), prep(dh, st)
    zall = prep(d, st)
    say(f"sample: fit half {len(dd)} rows ({int(dd.case.sum())} cases, {dd.season.nunique()} seasons); "
        f"test half {len(dh)} rows ({int(dh.case.sum())} cases, {dh.season.nunique()} seasons)")
    say(f"weighted base rate of G_T >= 71.7 kt among low-in-region times: fit half {np.average(dd.case, weights=dd.w):.4f}, "
        f"test half {np.average(dh.case, weights=dh.w):.4f}")
    say("SD used for 'per SD' (weighted, fit half): " + ", ".join(f"{c} {st[c][1]:.2f}" for c in ["GH", "GRAD", "STAB", "NAO", "GBI"]))
    say(f"missing: NAO {int(d.NAO.isna().sum())}, motion {int(d.mmiss.sum())} rows, GBI {int(d.GBI.isna().sum())}")

    # F1
    t0 = time.time()
    f1d, bd, namesd, _ = run_tests(zd, "logit", "F1 fit half", NBOOT, SEED)
    f1h, bh_, namesh, _ = run_tests(zh, "logit", "F1 test half", NBOOT, SEED + 1)
    say(f"\nF1 primary family (9 tests, BH q within the fit half), outcome Y_T, logit scale per SD  [{time.time() - t0:.0f}s]")
    say("test      fit-half est [95% CI]   p      q   |  test-half est [95% CI]   p(2-sided)  same sign + one-sided p<0.05?  -> supported")
    sup = {}
    for a, b in zip(f1d.itertuples(), f1h.itertuples()):
        if a.test == "MOT":
            same = None
            line = (f"{a.test:5s} Wald {a.wald:.2f} (u {a.bu:+.3f}, v {a.bv:+.3f}/kt) p {a.p:.3f} q {a.q:.3f} | "
                    f"Wald {b.wald:.2f} (u {b.bu:+.3f}, v {b.bv:+.3f}) p {b.p:.3f}")
            sup[a.test] = (a.q < 0.10) and (b.p < 0.05)
        else:
            same = np.sign(a.est) == np.sign(b.est)
            p1 = (b.p1 if same else 1.0)
            sup[a.test] = bool((a.q < 0.10) and same and (p1 < 0.05))
            line = (f"{a.test:5s} {a.est:+.3f} [{a.lo:+.3f}, {a.hi:+.3f}] p {a.p:.3f} q {a.q:.3f} | "
                    f"{b.est:+.3f} [{b.lo:+.3f}, {b.hi:+.3f}] p {b.p:.3f}  same sign {same}, one-sided p {p1:.3f}")
        line += f"  -> {'SUPPORTED' if sup[a.test] else 'not supported'}"
        say(line)
        rows.append(dict(family="F1", outcome="Y_T", test=a.test, fit_est=a.est, fit_lo=a.lo, fit_hi=a.hi, fit_p=a.p,
                         fit_q=a.q, test_est=b.est, test_lo=b.lo, test_hi=b.hi, test_p=b.p, supported=sup[a.test]))
    say("odds ratio per SD, fit half: " + ", ".join(
        f"{r.test} {np.exp(r.est):.2f} [{np.exp(r.lo):.2f}, {np.exp(r.hi):.2f}]" for r in f1d.itertuples()
        if r.test in ING))
    say("odds ratio per SD, test half: " + ", ".join(
        f"{r.test} {np.exp(r.est):.2f} [{np.exp(r.lo):.2f}, {np.exp(r.hi):.2f}]" for r in f1h.itertuples()
        if r.test in ING))

    # pooled, not out of sample
    f1p, bp, namesp, _ = run_tests(zall, "logit", "pooled", max(500, NBOOT // 2), SEED + 2)
    say("\npooled all seasons (fit-half z-scores, NOT out of sample), odds ratio per SD [season-block 95% CI]; p, q over the 9:")
    for r in f1p.itertuples():
        if r.test == "MOT":
            say(f"  MOT  Wald {r.wald:.2f}, u {r.bu:+.3f}, v {r.bv:+.3f} per kt, p {r.p:.3f}, q {r.q:.3f}")
        else:
            say(f"  {r.test:5s} {np.exp(r.est):.2f} [{np.exp(r.lo):.2f}, {np.exp(r.hi):.2f}] p {r.p:.3f} q {r.q:.3f}")
        rows.append(dict(family="F1-pooled", outcome="Y_T", test=r.test, fit_est=r.est, fit_lo=r.lo, fit_hi=r.hi,
                         fit_p=r.p, fit_q=r.q))

    # baseline coefficients (context)
    say("\nbaseline terms, pooled fit (logit per z, M1): " + ", ".join(
        f"{n} {bp[i + 1]:+.2f}" for i, n in enumerate(namesp) if n in ("depth", "x", "y", "trend")))

    # skill
    obs, ci, b1 = skill(zd, zh, NBOOT, SEED + 3)
    say("\nOut-of-sample skill: M1 and M0 fitted on the fit half, scored on the test half (weights restore all low-in-region times)")
    say(f"  AUC  M1 {obs['auc1']:.3f} [{ci['auc1'][0]:.3f}, {ci['auc1'][1]:.3f}]  M0 {obs['auc0']:.3f} [{ci['auc0'][0]:.3f}, {ci['auc0'][1]:.3f}]")
    say(f"  Brier  M1 {obs['brier1']:.4f}  M0 {obs['brier0']:.4f}  month climatology {obs['brierc']:.4f}")
    say(f"  Brier skill M1 vs M0 {obs['bss01']:+.3f} [{ci['bss01'][0]:+.3f}, {ci['bss01'][1]:+.3f}];"
        f" M1 vs month climatology {obs['bss1c']:+.3f} [{ci['bss1c'][0]:+.3f}, {ci['bss1c'][1]:+.3f}];"
        f" M0 vs month climatology {obs['bss0c']:+.3f} [{ci['bss0c'][0]:+.3f}, {ci['bss0c'][1]:+.3f}]")
    rows.append(dict(family="skill", outcome="Y_T", test="BSS M1 vs M0", fit_est=obs["bss01"], fit_lo=ci["bss01"][0], fit_hi=ci["bss01"][1]))
    rows.append(dict(family="skill", outcome="Y_T", test="BSS M1 vs climatology", fit_est=obs["bss1c"], fit_lo=ci["bss1c"][0], fit_hi=ci["bss1c"][1]))

    # plain-words risk change (fit-half model, applied to test-half population)
    f = 1.0 / np.average(zd.w[zd.case == 0])
    cols1 = BASE + ING + MOT
    bfit = logit(mat(zd, cols1), zd.case.values.astype(float))
    base_p = np.average(expit(mat(zh, cols1) @ bfit + np.log(f)), weights=zh.w)
    say(f"\nPlain-words effects on the test-half population (fit-half M1; mean predicted P(Y_T) {base_p:.3f}):")
    say("  ingredient at its 10th -> 90th percentile (z), P(G_T >= 71.7 kt):")
    for c in ING:
        l, h, pl, ph = g_comp(zh, bfit, cols1, c, f=f)
        say(f"  {c:5s} z {l:+.2f} -> {h:+.2f}:  {pl:.3f} -> {ph:.3f}  (x{ph / pl:.2f})")
        rows.append(dict(family="effect", outcome="Y_T", test=c, fit_est=pl, fit_lo=ph))

    # matched sensitivity and F2..F4 run only in the primary pass
    if SWAP:
        pd.DataFrame(rows).to_csv(f"{RES}/tests_swap.csv", index=False)
        open(f"{RES}/results_swap.txt", "w").write("\n".join(out) + "\n")
        return

    # F2 classes
    say("\nF2 classes: barrier cases (tip cases removed) against controls; 5 main-effect tests, own BH")
    for cls in ("barrier",):
        zdc = zd[(zd.case == 0) | (zd[cls])].reset_index(drop=True)
        zhc = zh[(zh.case == 0) | (zh[cls])].reset_index(drop=True)
        fa, _, _, _ = run_tests_main(zdc, cls + " fit", NBOOT, SEED + 10)
        fb, _, _, _ = run_tests_main(zhc, cls + " test", NBOOT, SEED + 11)
        say(f"  {cls}: fit half {int(zdc.case.sum())} cases, test half {int(zhc.case.sum())} cases")
        for a, b in zip(fa.itertuples(), fb.itertuples()):
            if a.test == "MOT":
                say(f"    {a.test:5s} Wald p {a.p:.3f} q {a.q:.3f} | test p {b.p:.3f}")
            else:
                say(f"    {a.test:5s} OR {np.exp(a.est):.2f} [{np.exp(a.lo):.2f}, {np.exp(a.hi):.2f}] p {a.p:.3f} q {a.q:.3f} | "
                    f"test OR {np.exp(b.est):.2f} [{np.exp(b.lo):.2f}, {np.exp(b.hi):.2f}] p {b.p:.3f}")
            rows.append(dict(family="F2", outcome=cls, test=a.test, fit_est=a.est, fit_lo=a.lo, fit_hi=a.hi, fit_p=a.p,
                             fit_q=a.q, test_est=b.est, test_lo=b.lo, test_hi=b.hi, test_p=b.p))
    nt_d, nt_h = int(zd.tip.sum()), int(zh.tip.sum())
    say(f"  tip: {nt_d} cases in the fit half and {nt_h} in the test half; the fit half is under the 40-case rule when below 40 -> described, not modelled")
    for c in ["GH", "GRAD", "STAB", "NAO", "GBI"]:
        a = zall[zall.tip]
        b = zall[zall.case == 0]
        mdiff = a[c].mean() - np.average(b[c].dropna(), weights=b.w[b[c].notna()])
        say(f"    tip cases minus weighted controls, {c}: {mdiff:+.2f} SD (n tip {int(a[c].notna().sum())}; descriptive)")

    # F3 continuous
    f3d, _, _, _ = run_tests_main(zd, "F3 fit", NBOOT, SEED + 20, kind="wls")
    f3h, _, _, _ = run_tests_main(zh, "F3 test", NBOOT, SEED + 21, kind="wls")
    say("\nF3 continuous G_T (kt per SD), weighted least squares, 5 tests, own BH")
    for a, b in zip(f3d.itertuples(), f3h.itertuples()):
        if a.test == "MOT":
            say(f"  MOT  Wald p {a.p:.3f} q {a.q:.3f} | test p {b.p:.3f}")
        else:
            say(f"  {a.test:5s} {a.est:+.2f} kt [{a.lo:+.2f}, {a.hi:+.2f}] p {a.p:.3f} q {a.q:.3f} | test {b.est:+.2f} [{b.lo:+.2f}, {b.hi:+.2f}] p {b.p:.3f}")
        rows.append(dict(family="F3", outcome="G_T", test=a.test, fit_est=a.est, fit_lo=a.lo, fit_hi=a.hi, fit_p=a.p, fit_q=a.q,
                         test_est=b.est, test_lo=b.lo, test_hi=b.hi, test_p=b.p))

    # F4 GBI
    say("\nF4 GBI (Z500, WeatherBench2 1.5 degree, times before 2023-01-10): 2 tests")
    for lbl, zset in (("fit half", zd), ("test half", zh)):
        zz = zset[zset.GBI.notna()].reset_index(drop=True)
        zz2 = zz.assign(GH_orig=zz.GH)
        ya = zz.case.values.astype(float)
        sidx = seasons_idx(zz)
        res = {}
        for name, cols in (("GBI replaces GH", BASE + ["GBI", "GRAD", "STAB", "NAO"] + MOT),
                           ("GBI added to M1", BASE + ING + ["GBI"] + MOT)):
            b0 = logit(mat(zz, cols), ya)
            j = cols.index("GBI") + 1
            rng = np.random.default_rng(SEED + 30)
            bb = np.array([logit(mat(zz.iloc[ix], cols), ya[ix]) [j] for ix in (boot_idx(sidx, rng) for _ in range(NBOOT // 2))])
            res[name] = (b0[j], np.percentile(bb, 2.5), np.percentile(bb, 97.5), 2 * min((bb >= 0).mean(), (bb <= 0).mean()) + 2 / NBOOT)
        qs = bh([res[k][3] for k in res])
        for (k, v), q in zip(res.items(), qs):
            say(f"  {lbl}: {k}: OR per SD {np.exp(v[0]):.2f} [{np.exp(v[1]):.2f}, {np.exp(v[2]):.2f}] p {v[3]:.3f} q {q:.3f} "
                f"({int(zz.case.sum())} cases, {zz.season.nunique()} seasons)")
            rows.append(dict(family="F4", outcome="Y_T", test=f"{lbl}: {k}", fit_est=v[0], fit_lo=v[1], fit_hi=v[2], fit_p=v[3], fit_q=q))

    # correlations among ingredients (context for collinearity)
    cc = zall[zall.w > 0][["GH", "GRAD", "STAB", "NAO"]].corr()
    say("\ningredient correlations (sample, unweighted): " + "; ".join(
        f"{a}-{b} {cc.loc[a, b]:+.2f}" for i, a in enumerate(cc.index) for b in cc.index[i + 1:]))

    pd.DataFrame(rows).to_csv(f"{RES}/tests.csv", index=False)
    open(f"{RES}/results.txt", "w").write("\n".join(out) + "\n")
    zall.to_pickle(f"{WORK}/zall.pkl")
    pd.Series(st).to_pickle(f"{WORK}/zstats.pkl")
    np.save(f"{WORK}/bfit.npy", bfit)


def run_tests_main(z, label, nboot, seed, kind="logit"):
    """Main-effect tests only (the F2 and F3 families)."""
    y = z.case.values.astype(float)
    sidx = seasons_idx(z.reset_index(drop=True))
    cols = BASE + ING + MOT
    zz = z.reset_index(drop=True)
    est = estimates(zz, cols, y, kind)
    rng = np.random.default_rng(seed)
    bs = np.zeros((nboot, len(cols) + 1))
    for b in range(nboot):
        ix = boot_idx(sidx, rng)
        bs[b] = estimates(zz.iloc[ix], cols, y[ix], kind)
    res = []
    for c in ING:
        j = cols.index(c) + 1
        bb = bs[:, j]
        res.append(dict(test=c, est=est[j], lo=np.percentile(bb, 2.5), hi=np.percentile(bb, 97.5),
                        p=2 * min((bb >= 0).mean(), (bb <= 0).mean()) + 1 / nboot))
    ju, jv = cols.index("u") + 1, cols.index("v") + 1
    bm = est[[ju, jv]]
    Bm = bs[:, [ju, jv]]
    Si = np.linalg.inv(np.cov(Bm.T))
    obs = float(bm @ Si @ bm)
    cen = Bm - Bm.mean(0)
    dist = np.einsum("ij,jk,ik->i", cen, Si, cen)
    res.append(dict(test="MOT", est=np.nan, lo=np.nan, hi=np.nan, p=(dist >= obs).mean() + 1 / nboot))
    r = pd.DataFrame(res)
    r["q"] = bh(r.p.values)
    return r, est, cols, bs


if __name__ == "__main__":
    main()
