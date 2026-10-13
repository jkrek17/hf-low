"""Score the 24 h gale / storm / HF outlook and rapid deepening, leave-one-season-out.
Plan: PREREGISTRATION.md. ERA5 proxy, pipeline A, seasons 2004-05 to 2025-26.
usage: score.py [OUTDIR]
"""
import os, sys, json
import numpy as np, pandas as pd
from joblib import Parallel, delayed
from sklearn.linear_model import LogisticRegression

HERE = os.path.dirname(os.path.abspath(__file__))
INT = os.path.join(HERE, "..", "intensity")
sys.path.insert(0, INT)
import model as M

OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "results")
CUTS = (38.0, 54.0, 71.7)
NAMES = ["below gale", "gale", "storm", "hurricane force"]
SETS = {"state": M.STATE,
        "A": M.STATE + ["jet250", "vadv500"],
        "B": M.STATE + ["jet250", "eady"],
        "full": M.SETS["full"],
        "full_nogust": M.SETS["full_nogust"]}
NBOOT = 2000


def load():
    D = M.load(os.path.join(INT, "results", "fixes_2004.csv.gz"), os.path.join(INT, "results", "env_2004.csv.gz"))
    key = dict(zip(zip(D.track, D.time), D.g800))
    t = pd.to_datetime(D.time, format="%Y%m%d%H")
    nxt = []
    for h in (12, 24):
        tt = (t + pd.Timedelta(hours=h)).dt.strftime("%Y%m%d%H")
        nxt.append(np.array([key.get((a, b), np.nan) for a, b in zip(D.track, tt)]))
    gmax = np.fmax(nxt[0], nxt[1])
    D["gmax24"] = np.nan_to_num(gmax, nan=0.0)
    D["gcls"] = np.digitize(D.gmax24.values, CUTS)
    return D


def fit_predict(S, te, cols, target, ncls):
    tr = ~te
    Xtr = S.loc[tr, cols].values.astype(float)
    prep = M.Prep(Xtr)
    m = LogisticRegression(C=1.0, max_iter=3000).fit(prep(Xtr), S.loc[tr, target].astype(int).values)
    p = np.zeros((te.sum(), ncls))
    p[:, m.classes_] = m.predict_proba(prep(S.loc[te, cols].values.astype(float)))
    return p


def clim(S, te, target, ncls):
    Dtr = S[~te]; Dte = S[te]
    p = np.zeros((len(Dte), ncls))
    for (b, mo), g in Dte.groupby(["basin", "month"]):
        t = Dtr[(Dtr.basin == b) & (Dtr.month == mo)][target].astype(int)
        cnt = np.bincount(t, minlength=ncls) + 0.5
        p[np.flatnonzero(Dte.index.isin(g.index))] = cnt / cnt.sum()
    return p


def loso(S, target, ncls):
    seasons = sorted(S.season.unique())
    def one(s):
        te = (S.season == s).values
        r = {"clim": clim(S, te, target, ncls)}
        for k, cols in SETS.items():
            r[k] = fit_predict(S, te, cols, target, ncls)
        return s, te, r
    res = Parallel(n_jobs=-1)(delayed(one)(s) for s in seasons)
    P = {k: np.zeros((len(S), ncls)) for k in ["clim"] + list(SETS)}
    for s, te, r in res:
        for k in P:
            P[k][te] = r[k]
    return P


def rps(P, y, ncls):
    Y = np.eye(ncls)[y]
    return ((np.cumsum(P, 1) - np.cumsum(Y, 1)) ** 2).sum(1)


def brier(p, y):
    return (p - y) ** 2


def boot_gain(season, sc_a, sc_b, sc_c, rng):
    """score arrays are per-fix losses for model a, model b, climatology; gain = skill_a - skill_b"""
    u = np.unique(season)
    sa = np.array([sc_a[season == s].sum() for s in u]); sb = np.array([sc_b[season == s].sum() for s in u])
    sc = np.array([sc_c[season == s].sum() for s in u])
    pt = (sb.sum() - sa.sum()) / sc.sum()
    g = []
    for _ in range(NBOOT):
        i = rng.integers(0, len(u), len(u))
        g.append((sb[i].sum() - sa[i].sum()) / sc[i].sum())
    g = np.array(g)
    p = min(1.0, 2 * min((g <= 0).mean(), (g >= 0).mean()) + 1 / NBOOT)
    return pt, np.percentile(g, 5), np.percentile(g, 95), p


def skill(sc_m, sc_c):
    return 1 - sc_m.sum() / sc_c.sum()


def bh(p):
    p = np.asarray(p); o = np.argsort(p); n = len(p)
    q = np.empty(n); prev = 1.0
    for rank in range(n, 0, -1):
        i = o[rank - 1]; prev = min(prev, p[i] * n / rank); q[i] = prev
    return q


def main():
    rng = np.random.default_rng(20261013)
    D = load()
    lines = []
    def out(s=""):
        print(s); lines.append(s)
    out(f"Fixes {len(D)}, tracks {D.track.nunique()}, seasons {D.season.min()}-{D.season.max()} (n = {D.season.nunique()})")
    fr = np.bincount(D.gcls, minlength=4) / len(D)
    out("Outlook class frequencies (next 24 h, 00/12 UTC samples): " + ", ".join(f"{n} {f:.3f}" for n, f in zip(NAMES, fr)))
    # ---- outlook
    P4 = loso(D, "gcls", 4)
    y4 = D.gcls.values
    L4 = {k: rps(P4[k], y4, 4) for k in P4}
    out("\n== Outlook (4 ordered classes): LOSO RPSS vs basin-month climatology, season-bootstrap 90 % CI")
    for k in P4:
        pt, lo, hi, _ = boot_gain(D.season.values, L4[k], L4["clim"], L4["clim"], rng)
        out(f"  {k:12s} {skill(L4[k], L4['clim']):+.3f} [{lo:+.3f}, {hi:+.3f}]")
    # ---- rapid deepening
    R = D[D.cls >= 0].copy().reset_index(drop=True)
    R["rd"] = (R.cls == 4).astype(int)
    P2 = loso(R, "rd", 2)
    L2 = {k: brier(P2[k][:, 1], R.rd.values) for k in P2}
    out("\n== Rapid deepening (NDR >= 1 Bergeron in 24 h), n = %d fixes, base rate %.3f: LOSO BSS vs climatology" % (len(R), R.rd.mean()))
    for k in P2:
        pt, lo, hi, _ = boot_gain(R.season.values, L2[k], L2["clim"], L2["clim"], rng)
        out(f"  {k:12s} {skill(L2[k], L2['clim']):+.3f} [{lo:+.3f}, {hi:+.3f}]")
    # ---- family of 6
    fam = []
    for tgt, S_, L in (("outlook RPSS", D, L4), ("rapid-deepening BSS", R, L2)):
        for k in ("A", "B", "full"):
            pt, lo, hi, p = boot_gain(S_.season.values, L[k], L["state"], L["clim"], rng)
            fam.append((tgt, k, pt, lo, hi, p))
    q = bh([f[5] for f in fam])
    out("\n== Pre-registered family (gain over state; BH q within the 6)")
    rows = []
    for f, qq in zip(fam, q):
        verdict = "yes" if (qq < 0.05 and f[3] > 0) else ("no" if f[4] < 0.005 else "can't tell")
        out(f"  {f[0]:20s} {f[1]:5s} - state: {f[2]:+.4f} [{f[3]:+.4f}, {f[4]:+.4f}]  p {f[5]:.4f}  q {qq:.4f}  {verdict}")
        rows.append(dict(target=f[0], model=f[1], gain=f[2], lo90=f[3], hi90=f[4], p=f[5], q=qq, verdict=verdict))
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "family.csv"), index=False)
    out("\n== Secondary: full over A / B")
    for tgt, S_, L in (("outlook RPSS", D, L4), ("rapid-deepening BSS", R, L2)):
        for k in ("A", "B"):
            pt, lo, hi, p = boot_gain(S_.season.values, L["full"], L[k], L["clim"], rng)
            out(f"  {tgt:20s} full - {k}: {pt:+.4f} [{lo:+.4f}, {hi:+.4f}]")
    out("  full vs full_nogust (outlook): %+.4f" % ((L4["full_nogust"].sum() - L4["full"].sum()) / L4["clim"].sum()))
    # per-class Brier skill (full vs state), per basin, left24 sensitivity
    out("\n== Per-class Brier skill, full vs clim (and gain over state)")
    for c in range(4):
        yc = (y4 == c).astype(float)
        bm = brier(P4["full"][:, c], yc); bs = brier(P4["state"][:, c], yc); bc = brier(P4["clim"][:, c], yc)
        pt, lo, hi, _ = boot_gain(D.season.values, bm, bs, bc, rng)
        out(f"  {NAMES[c]:16s} base {yc.mean():.4f}  BSS full {skill(bm, bc):+.3f}  gain over state {pt:+.4f} [{lo:+.4f}, {hi:+.4f}]")
    out("\n== Outlook RPSS gain full - state by basin and without left24 fixes")
    for lab, msk in (("Atlantic", D.basin == "atl"), ("Pacific", D.basin == "pac"), ("track continues 24 h", ~D.left24)):
        m = msk.values
        pt, lo, hi, _ = boot_gain(D.season.values[m], L4["full"][m], L4["state"][m], L4["clim"][m], rng)
        out(f"  {lab:22s} n {m.sum():6d}  {pt:+.4f} [{lo:+.4f}, {hi:+.4f}]")
    # reliability of the full outlook
    out("\n== Reliability (full outlook, pooled held-out seasons): forecast decile bin -> observed frequency")
    rel = []
    for c in range(4):
        p = P4["full"][:, c]; yc = (y4 == c)
        edges = np.unique(np.quantile(p, np.linspace(0, 1, 11)))
        b = np.clip(np.digitize(p, edges[1:-1]), 0, len(edges) - 2)
        for j in range(len(edges) - 1):
            m = b == j
            if m.sum():
                rel.append(dict(cls=NAMES[c], bin=j, n=int(m.sum()), forecast=p[m].mean(), observed=yc[m].mean()))
        r = pd.DataFrame([x for x in rel if x["cls"] == NAMES[c]])
        slope = np.polyfit(r.forecast, r.observed, 1, w=np.sqrt(r.n))[0]
        out(f"  {NAMES[c]:16s} reliability slope {slope:.2f}  (1 = perfect); forecast {r.forecast.min():.3f}-{r.forecast.max():.3f}")
    pd.DataFrame(rel).to_csv(os.path.join(OUT, "reliability.csv"), index=False)
    # most likely class hit rate
    ml = P4["full"].argmax(1)
    out("\n== Most-likely-class table (full), rows observed, columns forecast, counts")
    cm = pd.crosstab(pd.Categorical([NAMES[i] for i in y4], NAMES), pd.Categorical([NAMES[i] for i in ml], NAMES))
    out(cm.to_string())
    o = D[["track", "time"]].copy()
    for c in range(4):
        o["p%d" % c] = np.round(P4["full"][:, c], 4)
    o["gcls"] = y4
    o.to_csv(os.path.join(OUT, "oof_outlook_full.csv.gz"), index=False)
    open(os.path.join(OUT, "skill.txt"), "w").write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
