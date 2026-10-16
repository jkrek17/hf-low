"""Stage 3: episodes, regimes, the forecaster table (Q1) and strong-versus-storm tests (Q2). ERA5 PROXY.

Reads results/stage0_zone_gust.csv, pull_times.csv, features.csv. Writes results/episodes.csv, q1_table.csv, q2_tests.csv,
analysis.txt. Definitions and the tests are fixed in PREREGISTRATION.md (committed 2f74f02 before any field was read).
"""
import os
import numpy as np, pandas as pd
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, "results")
rng = np.random.default_rng(16)
out = []; P = lambda *a: out.append(" ".join(str(x) for x in a))
ZONES = ["Z1", "Z2", "Z3"]
NT = 9338; NCTRL = 600

s0 = pd.read_csv(os.path.join(RES, "stage0_zone_gust.csv"))
s0["time"] = pd.to_datetime(s0.time)
pt = pd.read_csv(os.path.join(RES, "pull_times.csv"))
f = pd.read_csv(os.path.join(RES, "features.csv"), parse_dates=["time"])
ctrl_times = set(pd.to_datetime(pt[pt.is_control].time))


def season(t):
    return t.year if t.month >= 10 else t.year - 1


def regime(zn, wd):
    if zn in ("Z1", "Z2"):
        if 100 <= wd < 170: return "BARRIER"
        if 170 <= wd < 260: return "ONSHORE"
        return "OFFSHORE"
    if wd >= 300 or wd < 60: return "GAP"
    if 90 <= wd < 200: return "INFLOW"
    return "OTHER"


f["regime"] = [regime(z, w) for z, w in zip(f.zone, f.wdir)]
f["season"] = [season(t) for t in f.time]
f["is_ctrl"] = f.time.isin(ctrl_times)
f["event"] = f.gust_kt >= 54
f["strong_t"] = f.gust_kt >= 64

# ---- episodes: per zone, event times linked when 36 h or less apart, peak time represents the episode
eps = []
for zn in ZONES:
    d = f[(f.zone == zn) & f.event].sort_values("time")
    # event times that were not pulled cannot exist: every event time is pulled
    if len(d) != int((s0[f"{zn}_gust_kt"] >= 54).sum()): raise SystemExit(f"{zn}: event times missing from the pull")
    grp = (d.time.diff() > pd.Timedelta(hours=36)).cumsum()
    for _, g in d.groupby(grp):
        pk = g.loc[g.gust_kt.idxmax()].copy()
        pk["n_times"] = len(g); pk["ep_max"] = g.gust_kt.max(); pk["strong_ep"] = int(g.gust_kt.max() >= 64)
        pk["hf_ep"] = int(g.gust_kt.max() >= 71.7)
        eps.append(pk)
E = pd.DataFrame(eps).reset_index(drop=True)
E.to_csv(os.path.join(RES, "episodes.csv"), index=False)
P("EPISODES (per zone, 36 h linking, peak-time regime)")
t = E.groupby(["zone", "regime"]).agg(n=("time", "size"), strong=("strong_ep", "sum"), hf=("hf_ep", "sum")).reset_index()
t["strong_share"] = (t.strong / t.n).round(3); t["per_season"] = (t.n / 22).round(2)
P(t.to_string(index=False))
P("\nepisodes per zone:", E.groupby("zone").size().to_dict(), " strong episodes:", E.groupby("zone").strong_ep.sum().to_dict())

# ---- Q1: chance that a time is an event given regime (controls give the regime shares)
FACT = {"Z1": ["DEPTH", "LOW_DIST", "GRAD", "COASTPILE", "STAB"], "Z2": ["DEPTH", "LOW_DIST", "GRAD", "COASTPILE", "STAB"],
        "Z3": ["DEPTH", "LOW_DIST", "GRAD", "INTERIOR_DP", "STAB"]}
BANDS = [("bottom half", 0, 50), ("50th-80th", 50, 80), ("80th-95th", 80, 95), ("top 5%", 95, 100.001)]
seasons = sorted(f.season.unique())
bys = {s: g for s, g in f.groupby("season")}


def rate(d, mask_fn):
    """chance (percent) that a time with the property is an event / strong: events x NCTRL / (NT x controls)"""
    c = d[d.is_ctrl & mask_fn(d)]
    e = d[d.event & mask_fn(d)]
    s = d[d.strong_t & mask_fn(d)]
    if len(c) == 0: return np.nan, np.nan
    return 100 * len(e) * NCTRL / (NT * len(c)), 100 * len(s) * NCTRL / (NT * len(c))


def boot_rate(d, mask_fn, n=500):
    g = {s: x for s, x in d.groupby("season")}; ss = list(g); r = []
    for _ in range(n):
        dd = pd.concat([g[s] for s in rng.choice(ss, len(ss))])
        r.append(rate(dd, mask_fn))
    r = np.array(r, float)
    return np.nanpercentile(r[:, 0], [2.5, 97.5]), np.nanpercentile(r[:, 1], [2.5, 97.5])


rows = []
for zn in ZONES:
    d = f[f.zone == zn]
    for rg in sorted(d.regime.unique()):
        dr = d[d.regime == rg]
        nc = int(dr.is_ctrl.sum())
        e, s = rate(dr, lambda x: np.ones(len(x), bool))
        # NB: the regime share among times is nc/NCTRL; bootstrap below uses the same estimator
        ci_e, ci_s = boot_rate(d, lambda x, rg=rg: x.regime == rg)
        rows.append(dict(zone=zn, regime=rg, ingredient="(all times in regime)", band="all", n_ctrl=nc, n_event_t=int(dr.event.sum()),
                         n_strong_t=int(dr.strong_t.sum()), pct_event=e, pct_event_lo=ci_e[0], pct_event_hi=ci_e[1],
                         pct_strong=s, pct_strong_lo=ci_s[0], pct_strong_hi=ci_s[1]))
        if nc < 40: continue
        for fac in FACT[zn]:
            x = dr[fac]
            cx = dr[dr.is_ctrl][fac].dropna()
            if len(cx) < 40: continue
            cuts = np.nanpercentile(cx, [50, 80, 95])
            # DEPTH, GRAD, COASTPILE, INTERIOR_DP: larger = stronger setup; LOW_DIST and STAB: bands are on the raw value, direction is read from the table
            for bn, lo, hi in BANDS:
                qlo = -np.inf if lo == 0 else cuts[[50, 80, 95].index(lo)]
                qhi = np.inf if hi > 100 else cuts[[50, 80, 95].index(hi)]
                fn = lambda xx, fac=fac, qlo=qlo, qhi=qhi, rg=rg: (xx.regime == rg) & (xx[fac] >= qlo) & (xx[fac] < qhi)
                e, s = rate(d, fn)
                ci_e, ci_s = boot_rate(d, fn)
                rows.append(dict(zone=zn, regime=rg, ingredient=fac, band=bn, lo=qlo, hi=qhi, n_ctrl=int((d.is_ctrl & fn(d)).sum()),
                                 n_event_t=int((d.event & fn(d)).sum()), n_strong_t=int((d.strong_t & fn(d)).sum()),
                                 pct_event=e, pct_event_lo=ci_e[0], pct_event_hi=ci_e[1], pct_strong=s, pct_strong_lo=ci_s[0], pct_strong_hi=ci_s[1]))
Q1 = pd.DataFrame(rows)
Q1.to_csv(os.path.join(RES, "q1_table.csv"), index=False)
P("\nQ1 (chance a 00/12 UTC time in the regime reaches 54 / 64 kt zone gust; regime assigned to every time from the zone-mean wind)")
P(Q1[Q1.ingredient == "(all times in regime)"][["zone", "regime", "n_ctrl", "n_event_t", "n_strong_t", "pct_event", "pct_event_lo", "pct_event_hi", "pct_strong", "pct_strong_lo", "pct_strong_hi"]].round(1).to_string(index=False))


# ---- Q2: strong versus storm among episodes in Z1 BARRIER, Z2 BARRIER, Z3 GAP
def fit(X, y, lam=0.5):
    X1 = np.column_stack([np.ones(len(y)), X])
    fn = lambda b: np.sum(np.logaddexp(0, X1 @ b) - y * (X1 @ b)) + lam * np.sum(b[1:] ** 2)
    return minimize(fn, np.zeros(X1.shape[1]), method="BFGS").x


def est(d, fac):
    cols = [fac + "_z"] + ([] if fac == "DEPTH" else ["DEPTH_z"])
    return fit(d[cols].values.astype(float), d.strong_ep.values.astype(float))[1]


def boot(d, fac, n=2000):
    g = {s: x for s, x in d.groupby("season")}; ss = list(g); r = []
    for _ in range(n):
        dd = pd.concat([g[s] for s in rng.choice(ss, len(ss))])
        if dd.strong_ep.nunique() > 1: r.append(est(dd, fac))
    return np.array(r)


def bh(p):
    p = np.asarray(p); o = np.argsort(p); q = np.empty(len(p)); m = len(p)
    q[o] = np.minimum.accumulate((p[o] * m / (np.arange(m) + 1))[::-1])[::-1]; return np.minimum(q, 1)


CLASSES = [("Z1", "BARRIER"), ("Z2", "BARRIER"), ("Z3", "GAP")]
rows = []
for zn, rg in CLASSES:
    d = E[(E.zone == zn) & (E.regime == rg)].copy()
    for fac in FACT[zn]:
        d = d.dropna(subset=[fac]) if fac in d else d
    d["half"] = np.where(d.season % 2 == 0, "even", "odd")
    nstrong = {h: int(d[d.half == h].strong_ep.sum()) for h in ("even", "odd")}
    nstorm = {h: int((1 - d[d.half == h].strong_ep).sum()) for h in ("even", "odd")}
    P(f"\n{zn} {rg}: {len(d)} episodes, {int(d.strong_ep.sum())} strong; strong per half {nstrong}, storm per half {nstorm}")
    for fac in FACT[zn]:
        d[fac + "_z"] = (d[fac] - d[fac].mean()) / d[fac].std()
    for fac in FACT[zn]:
        # DEPTH must be standardised before use as adjustment: done above
        b = boot(d, fac); e0 = est(d, fac)
        prev = d.strong_ep.mean()
        rows.append(dict(zone=zn, regime=rg, factor=fac, n=len(d), strong=int(d.strong_ep.sum()), OR=np.exp(e0), lo=np.exp(np.percentile(b, 2.5)),
                         hi=np.exp(np.percentile(b, 97.5)), p=2 * min((b <= 0).mean(), (b >= 0).mean()) if len(b) else np.nan,
                         mdor=float(np.exp(2.8 / np.sqrt(len(d) * prev * (1 - prev)))) if 0 < prev < 1 else np.nan,
                         halves_ok=min(min(nstrong.values()), min(nstorm.values())) >= 30))
Q2 = pd.DataFrame(rows)
Q2["q"] = bh(Q2.p.fillna(1).values)
Q2.to_csv(os.path.join(RES, "q2_tests.csv"), index=False)
P("\nQ2 all-seasons estimates (NOT HELD OUT; fallback report pre-registered). OR per SD, strong (>=64 kt) vs 54-64 kt episodes; BH over", len(Q2))
P(Q2.round(3).to_string(index=False))
P("\nhalf-split rule (30 of the rarer outcome in each half) met by any class:", bool(Q2.halves_ok.any()))
open(os.path.join(RES, "analysis.txt"), "w").write("\n".join(out)); print("\n".join(out))
