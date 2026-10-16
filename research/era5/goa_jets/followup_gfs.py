"""Lead 2 (PREREGISTRATION_followup.md): do GFS FORECAST pressure patterns at 24 and 48 h separate ERA5-PROXY event times from non-events?

Reads results/gfs_features.csv and stage0_zone_gust.csv. Writes results/followup_gfs_t1.csv (band table), followup_gfs_t2.csv (leave-one-season-out
skill), followup_gfs.txt. Truth = ERA5 proxy zone gust >= 54 kt (event) / >= 64 kt (strong) at the valid time. 5 seasons 2021-22..2025-26.
"""
import os
import numpy as np, pandas as pd
from scipy.optimize import minimize
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, "results")
rng = np.random.default_rng(16)
out = []; P = lambda *a: out.append(" ".join(str(x) for x in a))
g = pd.read_csv(os.path.join(RES, "gfs_features.csv"), parse_dates=["time"])
s0 = pd.read_csv(os.path.join(RES, "stage0_zone_gust.csv"), parse_dates=["time"])


def regime(zn, wd):
    if zn in ("Z1", "Z2"):
        if 100 <= wd < 170: return "BARRIER"
        if 170 <= wd < 260: return "ONSHORE"
        return "OFFSHORE"
    if wd >= 300 or wd < 60: return "GAP"
    if 90 <= wd < 200: return "INFLOW"
    return "OTHER"


g["regime"] = [regime(z, w) for z, w in zip(g.zone, g.wdir)]
g["season"] = [t.year if t.month >= 10 else t.year - 1 for t in g.time]
tr = pd.concat([s0[["time", f"{z}_gust_kt"]].rename(columns={f"{z}_gust_kt": "truth_gust"}).assign(zone=z) for z in ("Z1", "Z2", "Z3")])
g = g.merge(tr, on=["time", "zone"], how="left")
assert g.truth_gust.notna().all()
g["event"] = (g.truth_gust >= 54).astype(int); g["strong"] = (g.truth_gust >= 64).astype(int)
NV = g.groupby("lead").time.nunique().to_dict()
P("valid times per lead:", NV, "; seasons", sorted(g.season.unique()))
P("event base rate (all regimes) by zone and lead, %:"); P((100 * g.groupby(["zone", "lead"]).event.mean()).round(1).to_string())
FACT = {"Z1": ["GFSGUST", "DEPTH", "LOW_DIST", "COASTPILE"], "Z2": ["GFSGUST", "DEPTH", "LOW_DIST", "COASTPILE"], "Z3": ["GFSGUST", "DEPTH", "LOW_DIST", "INTERIOR_DP"]}
BANDS = [("bottom half", 0, 50), ("50th-80th", 50, 80), ("80th-95th", 80, 95), ("top 5%", 95, 100.001)]
KEY = {"Z1": "BARRIER", "Z2": "BARRIER", "Z3": "GAP"}

# ---- T1
rows = []
for zn in ("Z1", "Z2", "Z3"):
    for lead in (24, 48):
        d = g[(g.zone == zn) & (g.lead == lead)]
        for rg in sorted(d.regime.unique()):
            dr = d[d.regime == rg]
            rows.append(dict(zone=zn, lead=lead, regime=rg, ingredient="(all times in regime)", band="all", n=len(dr), events=int(dr.event.sum()), strong=int(dr.strong.sum()),
                             pct_event=100 * dr.event.mean(), pct_strong=100 * dr.strong.mean(),
                             season_min=100 * dr.groupby("season").event.mean().min(), season_max=100 * dr.groupby("season").event.mean().max()))
            if len(dr) < 80: continue
            for fac in FACT[zn]:
                x = dr[fac]; cuts = np.nanpercentile(x, [50, 80, 95])
                for bn, lo, hi in BANDS:
                    qlo = -np.inf if lo == 0 else cuts[[50, 80, 95].index(lo)]
                    qhi = np.inf if hi > 100 else cuts[[50, 80, 95].index(hi)]
                    sel = dr[(dr[fac] >= qlo) & (dr[fac] < qhi)]
                    if len(sel) == 0: continue
                    bs = {s: x_ for s, x_ in dr.groupby("season")}; ss = list(bs); r = []
                    for _ in range(1000):
                        dd = pd.concat([bs[s] for s in rng.choice(ss, len(ss))]); sl = dd[(dd[fac] >= qlo) & (dd[fac] < qhi)]
                        if len(sl): r.append(100 * sl.event.mean())
                    ps = sel.groupby("season").event.mean() * 100
                    rows.append(dict(zone=zn, lead=lead, regime=rg, ingredient=fac, band=bn, lo=qlo, hi=qhi, n=len(sel), events=int(sel.event.sum()), strong=int(sel.strong.sum()),
                                     pct_event=100 * sel.event.mean(), pct_event_lo=np.percentile(r, 2.5), pct_event_hi=np.percentile(r, 97.5), pct_strong=100 * sel.strong.mean(),
                                     season_min=ps.min(), season_max=ps.max()))
T1 = pd.DataFrame(rows); T1.to_csv(os.path.join(RES, "followup_gfs_t1.csv"), index=False)
P("\nT1 all times in regime (forecast regime; ERA5-proxy event >= 54 kt at the valid time)")
P(T1[T1.ingredient == "(all times in regime)"][["zone", "lead", "regime", "n", "events", "strong", "pct_event", "pct_strong", "season_min", "season_max"]].round(1).to_string(index=False))
for zn in ("Z1", "Z2", "Z3"):
    for fac in FACT[zn]:
        P(f"\nT1 {zn} {KEY[zn]} forecast regime, band of {fac}: chance (%) of an ERA5-proxy event [season-block 95%], events/times")
        for lead in (24, 48):
            t = T1[(T1.zone == zn) & (T1.lead == lead) & (T1.regime == KEY[zn]) & (T1.ingredient == fac)]
            P(f"  {lead} h: " + " | ".join(f"{r.band}: {r.pct_event:.1f} [{r.pct_event_lo:.1f},{r.pct_event_hi:.1f}] {int(r.events)}/{int(r.n)}" for r in t.itertuples()))


# ---- T2
def fit(X, y, lam=0.5):
    X1 = np.column_stack([np.ones(len(y)), X])
    fn = lambda b: np.sum(np.logaddexp(0, X1 @ b) - y * (X1 @ b)) + lam * np.sum(b[1:] ** 2)
    return minimize(fn, np.zeros(X1.shape[1]), method="BFGS").x


def predict(b, X):
    return 1 / (1 + np.exp(-(np.column_stack([np.ones(len(X)), X]) @ b)))


def auc(y, p):
    from scipy.stats import rankdata
    r = rankdata(p); n1 = y.sum(); n0 = len(y) - n1
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0) if n1 and n0 else np.nan


rows = []
for zn in ("Z1", "Z2", "Z3"):
    for lead in (24, 48):
        d = g[(g.zone == zn) & (g.lead == lead) & (g.regime == KEY[zn])].copy()
        ing = [f for f in FACT[zn] if f != "GFSGUST"]
        ya = d.event.values.astype(float); pa = np.zeros(len(d)); pb = np.zeros(len(d)); pc = np.zeros(len(d))
        for s in sorted(d.season.unique()):
            te = (d.season == s).values; trn = ~te
            for cols, pp in ((["GFSGUST"], pa), (["GFSGUST"] + ing, pb)):
                mu = d.loc[trn, cols].mean(); sd = d.loc[trn, cols].std()
                Xt = ((d.loc[trn, cols] - mu) / sd).values; Xe = ((d.loc[te, cols] - mu) / sd).values
                pp[te] = predict(fit(Xt, ya[trn]), Xe)
            pc[te] = ya[trn].mean()
        ll = lambda p: -(ya * np.log(np.clip(p, 1e-6, 1)) + (1 - ya) * np.log(np.clip(1 - p, 1e-6, 1)))
        dif = ll(pb) - ll(pa)       # negative = ingredients help
        sea = {s: (dif[(d.season == s).values].sum(), int((d.season == s).sum())) for s in sorted(d.season.unique())}
        ss = list(sea); r = []
        for _ in range(2000):
            pick = rng.choice(ss, len(ss)); r.append(sum(sea[s][0] for s in pick) / sum(sea[s][1] for s in pick))
        r = np.array(r); est = dif.mean()
        rows.append(dict(zone=zn, lead=lead, regime=KEY[zn], n=len(d), events=int(ya.sum()), ll_clim=ll(pc).mean(), ll_gust=ll(pa).mean(), ll_gust_plus=ll(pb).mean(),
                         diff=est, lo=np.percentile(r, 2.5), hi=np.percentile(r, 97.5), p=2 * min((r >= 0).mean(), (r <= 0).mean()),
                         auc_gust=auc(ya.astype(int), pa), auc_gust_plus=auc(ya.astype(int), pb), seasons_improved=int(sum(v[0] < 0 for v in sea.values()))))
T2 = pd.DataFrame(rows)
p = T2.p.clip(lower=1 / 2000).values; o = np.argsort(p); q = np.empty(len(p)); m = len(p)
q[o] = np.minimum.accumulate((p[o] * m / (np.arange(m) + 1))[::-1])[::-1]; T2["q"] = np.minimum(q, 1)
T2["verdict"] = np.where((T2.hi < 0) & (T2.q < 0.05), "ADDS", "cannot tell / no gain")
T2.to_csv(os.path.join(RES, "followup_gfs_t2.csv"), index=False)
P("\nT2 leave-one-season-out log-loss (per time), forecast regime only. diff = (gust + pressure ingredients) minus (gust alone); negative = ingredients help. BH over 6")
P(T2.round(4).to_string(index=False))
open(os.path.join(RES, "followup_gfs.txt"), "w").write("\n".join(out)); print("\n".join(out))
