"""Fit and cross-validate the intensity-class and hurricane-force models.

usage: model.py FIXES_CSV ENV_CSV OUTDIR

Population: every in-domain ERA5 low fix at 00 or 12 UTC, seasons 1979-2025.
Targets
  cls   24 h normalised deepening rate (Bergerons) in five ordered classes:
        rapid decay <= -1 < decay <= -0.3 < steady < 0.3 <= deepening < 1 <= rapid deepening
  hf24  gust index reaches 71.7 kt in (t, t+24 h]
  hf48  the same in (t, t+48 h]
Models (all logistic, L2, standardised predictors)
  clim   basin x calendar-month frequencies
  state  storm state only (persistence and climatology)
  hart   state + Hart B, VTL, VTU
  full   state + Hart + environment
Validation: leave-one-season-out. Uncertainty: resampling whole seasons, so the
effective sample is the number of seasons, not the fix count.

Gate: ERA5 gusts drift upward at fixed storm depth before 2001 (hf-low commit
c9dc994), and the project rule is that no season before 2004-05 is used to fit
or to test. Every model here is fitted and scored on seasons FIT_FROM (2004)
to 2025 only. Earlier fixes may have the fitted model applied to them, with
gust-based numbers before 2001 caveated. "full_nogust" drops the current gust
predictor to show how much the class model leans on the gust index.
"""
import os, sys, json
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

CLASSES = ["rapid decay", "decay", "steady", "deepening", "rapid deepening"]
EDGES = [-np.inf, -1.0, -0.3, 0.3, 1.0, np.inf]
STATE = ["msl", "dp12", "young", "lat", "speed", "g800", "logage", "pac", "doy_c", "doy_s"]
HART = ["B", "VTL", "VTU"]
ENV = ["jet250", "div300", "vadv500", "eady", "sst", "sstgrad", "sst_t500", "flux", "tcwv", "nosst"]
SETS = {"state": STATE, "hart": STATE + HART, "full": STATE + HART + ENV,
        "full_nogust": [c for c in STATE + HART + ENV if c != "g800"]}
FIT_FROM = 2004
NDR_MAX = 3.0          # |NDR| above this is almost always a tracker relink, not a storm
NBOOT = int(os.environ.get("NBOOT", 1000))


def load(fx, ev):
    F = pd.read_csv(fx, dtype={"time": str})
    E = pd.read_csv(ev, dtype={"time": str})
    D = F.merge(E, on=["track", "time"], how="inner")
    D = D[(D.season >= FIT_FROM) & (D.season <= 2025)].copy()
    D["young"] = D.dp12.isna().astype(float)
    D["dp12"] = D.dp12.fillna(0.0)
    D["logage"] = np.log1p(D.age)
    D["pac"] = (D.basin == "pac").astype(float)
    doy = pd.to_datetime(D.time, format="%Y%m%d%H").dt.dayofyear.values
    D["doy_c"], D["doy_s"] = np.cos(2 * np.pi * doy / 365.25), np.sin(2 * np.pi * doy / 365.25)
    D["month"] = D.time.str[4:6].astype(int)
    D["nosst"] = D.sst.isna().astype(float)
    ok = D.ndr24.notna() & (D.ndr24.abs() <= NDR_MAX)
    D["cls"] = np.where(ok, pd.cut(D.ndr24, EDGES, labels=False, right=False), -1)
    # rapid decay is <= -1 and decay <= -0.3: right-closed on the decay side
    D.loc[ok & (D.ndr24 == -1.0), "cls"] = 0
    D.loc[ok & (D.ndr24 == -0.3), "cls"] = 1
    return D.reset_index(drop=True)


class Prep:
    """Median-impute, clip to the 0.5-99.5 % range, standardise; fit on training rows only."""
    def __init__(self, X):
        self.med = np.nanmedian(X, 0)
        Xi = np.where(np.isnan(X), self.med, X)
        self.lo, self.hi = np.percentile(Xi, 0.5, 0), np.percentile(Xi, 99.5, 0)
        Xc = np.clip(Xi, self.lo, self.hi)
        self.mu, self.sd = Xc.mean(0), Xc.std(0) + 1e-9

    def __call__(self, X):
        X = np.where(np.isnan(X), self.med, X)
        return (np.clip(X, self.lo, self.hi) - self.mu) / self.sd

    def to_json(self):
        return {k: getattr(self, k).tolist() for k in ("med", "lo", "hi", "mu", "sd")}


def fit(X, y):
    m = LogisticRegression(C=1.0, max_iter=3000)
    return m.fit(X, y)


def clim_probs(Dtr, Dte, target, ncls):
    p = np.zeros((len(Dte), ncls))
    for (b, mo), g in Dte.groupby(["basin", "month"]):
        t = Dtr[(Dtr.basin == b) & (Dtr.month == mo)][target]
        cnt = np.bincount(t.astype(int), minlength=ncls) + 0.5
        p[Dte.index.get_indexer(g.index)] = cnt / cnt.sum()
    return p


def events(target):
    """Outcome columns scored as yes/no events: rapid decay and rapid deepening, or HF."""
    return (0, 4) if target == "cls" else (1,)


def loso(D, target, ncls):
    """Leave-one-season-out probabilities for every model, and for each yes/no
    event the probability cut that makes the training seasons' forecast count
    equal their observed count (bias 1), applied unchanged to the held-out season."""
    rows = D.index[D[target] >= 0] if target == "cls" else D.index
    S = D.loc[rows]
    P = {k: np.zeros((len(S), ncls)) for k in ["clim"] + list(SETS)}
    CUT = {k: np.zeros((len(S), ncls)) for k in SETS}
    for s in sorted(S.season.unique()):
        te = (S.season == s).values
        tr = ~te
        Dtr, Dte = S[tr], S[te]
        P["clim"][te] = clim_probs(Dtr, Dte.reset_index(drop=True).set_index(np.arange(te.sum())), target, ncls)
        y = Dtr[target].astype(int).values
        for name, cols in SETS.items():
            prep = Prep(Dtr[cols].values.astype(float))
            m = fit(prep(Dtr[cols].values.astype(float)), y)
            P[name][te] = m.predict_proba(prep(Dte[cols].values.astype(float)))
            ptr = m.predict_proba(prep(Dtr[cols].values.astype(float)))
            for k in events(target):
                CUT[name][te, k] = np.quantile(ptr[:, k], 1 - np.mean(y == k))
    return S, P, CUT


def table(fc, ob):
    """POD, FAR, CSI, HSS, bias from yes/no forecasts and observations."""
    a = np.sum(fc & ob); b = np.sum(fc & ~ob); c = np.sum(~fc & ob); d = np.sum(~fc & ~ob)
    n = a + b + c + d
    ex = ((a + b) * (a + c) + (c + d) * (b + d)) / n
    return dict(pod=a / max(a + c, 1), far=b / max(a + b, 1), csi=a / max(a + b + c, 1),
                hss=(a + d - ex) / max(n - ex, 1), bias=(a + b) / max(a + c, 1), n_obs=int(a + c))


def categorical(S, P, CUT, target, w):
    """Forecaster-style scores at the training-count-matched cut, pooled and by season."""
    y = S[target].astype(int).values
    seas = S.season.values
    names = {0: "rapid decay", 4: "rapid deepening", 1: target}
    for k in events(target):
        w(f"  yes/no {names[k]} at the count-matched cut (bias 1 in training): "
          "POD FAR CSI HSS bias, then the held-out-season HSS range")
        for name in SETS:
            fc, ob = P[name][:, k] >= CUT[name][:, k], y == k
            t = table(fc, ob)
            hs = [table(fc[seas == s], ob[seas == s])["hss"] for s in np.unique(seas) if ob[seas == s].sum() >= 5]
            w(f"    {name:11s} POD {t['pod']:.2f} FAR {t['far']:.2f} CSI {t['csi']:.2f} HSS {t['hss']:.2f} "
              f"bias {t['bias']:.2f}  (n_obs {t['n_obs']}; season HSS {min(hs):.2f}-{max(hs):.2f})")


def rps(p, y):
    """Ranked probability score per row for ordered classes."""
    o = np.zeros_like(p); o[np.arange(len(y)), y] = 1
    return np.sum((np.cumsum(p, 1) - np.cumsum(o, 1)) ** 2, 1) / (p.shape[1] - 1)


def brier(p, y):
    return (p - y) ** 2


def season_boot(score_by_model, seasons, ref="clim", rng=np.random.default_rng(0)):
    """Skill vs ref with a 90 % CI from resampling whole seasons."""
    us = np.unique(seasons)
    idx = {s: np.where(seasons == s)[0] for s in us}
    sums = {k: np.array([v[idx[s]].sum() for s in us]) for k, v in score_by_model.items()}
    cnt = np.array([len(idx[s]) for s in us])
    out = {}
    B = rng.integers(0, len(us), (NBOOT, len(us)))
    for k in score_by_model:
        sk = 1 - sums[k].sum() / sums[ref].sum()
        bs = 1 - sums[k][B].sum(1) / sums[ref][B].sum(1)
        out[k] = (sk, *np.percentile(bs, [5, 95]))
    return out, sums, B


def gain_ci(sums, B, a, b, ref="clim"):
    """Difference in skill a - b with a 90 % season-bootstrap CI."""
    d = (sums[b].sum() - sums[a].sum()) / sums[ref].sum()
    db = (sums[b][B].sum(1) - sums[a][B].sum(1)) / sums[ref][B].sum(1)
    return d, *np.percentile(db, [5, 95])


def reliability(p, y, bins=(0, .05, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1.0001)):
    b = np.digitize(p, bins) - 1
    return [(bins[i], bins[i + 1], int((b == i).sum()), float(p[b == i].mean()), float(y[b == i].mean()))
            for i in range(len(bins) - 1) if (b == i).sum() > 0]


def coef_table(D, target, cols, labels):
    """Final fit on all seasons, standardised coefficients with season-bootstrap 90 % CIs."""
    S = D[D[target] >= 0] if target == "cls" else D
    X, y = S[cols].values.astype(float), S[target].astype(int).values
    prep = Prep(X)
    m = fit(prep(X), y)
    us = S.season.unique()
    rng = np.random.default_rng(1)
    boots = []
    for _ in range(100):
        pick = rng.choice(us, len(us))
        ii = np.concatenate([np.where(S.season.values == s)[0] for s in pick])
        boots.append(fit(prep(X[ii]), y[ii]).coef_)
    boots = np.array(boots)
    lo, hi = np.percentile(boots, [5, 95], 0)
    rows = []
    for k, lab in enumerate(labels):
        for j, c in enumerate(cols):
            rows.append(dict(target=target, outcome=lab, predictor=c, coef=m.coef_[k, j],
                             ci5=lo[k, j], ci95=hi[k, j]))
    model = dict(target=target, outcomes=labels, predictors=cols, prep=prep.to_json(),
                 coef=m.coef_.tolist(), intercept=m.intercept_.tolist())
    return pd.DataFrame(rows), model


def predict(model, rows):
    """Probabilities from a saved model.json entry for a DataFrame of predictors
    (the columns named in model["predictors"], built as in load())."""
    p = model["prep"]
    X = rows[model["predictors"]].values.astype(float)
    X = np.where(np.isnan(X), p["med"], X)
    X = (np.clip(X, p["lo"], p["hi"]) - np.array(p["mu"])) / np.array(p["sd"])
    z = X @ np.array(model["coef"]).T + np.array(model["intercept"])
    if z.shape[1] == 1:
        return 1 / (1 + np.exp(-z[:, 0]))
    z = np.exp(z - z.max(1, keepdims=True))
    return z / z.sum(1, keepdims=True)


def main(fx, ev, outdir):
    os.makedirs(outdir, exist_ok=True)
    D = load(fx, ev)
    rep = []
    w = rep.append
    w(f"Fixes: {len(D)} in-domain 00/12 UTC fixes on {D.track.nunique()} tracks, seasons "
      f"{D.season.min()}-{D.season.max()} ({D.season.nunique()} seasons)")
    w(f"Class target available for {(D.cls >= 0).sum()} fixes "
      f"({(D.ndr24.abs() > NDR_MAX).sum()} dropped as |NDR| > {NDR_MAX}, "
      f"{D.ndr24.isna().sum()} where the track ends within 24 h)")
    w("Class frequencies: " + ", ".join(f"{c} {np.mean(D.cls[D.cls >= 0] == k):.3f}" for k, c in enumerate(CLASSES)))
    w(f"HF frequencies: now {D.hf_now.mean():.3f}, within 24 h {D.hf24.mean():.3f}, within 48 h {D.hf48.mean():.3f}")
    w("")
    models = {}
    coefs = []
    rel = {}
    for target, ncls, labels in (("cls", 5, CLASSES), ("hf24", 2, ["no", "yes"]), ("hf48", 2, ["no", "yes"])):
        D[target] = D[target].astype(int) if target != "cls" else D[target]
        Dt = D
        S, P, CUT = loso(Dt, target, ncls)
        y = S[target].astype(int).values
        seas = S.season.values
        if target == "cls":
            sc = {k: rps(p, y) for k, p in P.items()}
            name = "RPSS"
        else:
            sc = {k: brier(p[:, 1], y) for k, p in P.items()}
            name = "BSS"
        sk, sums, B = season_boot(sc, seas)
        w(f"== {target}: leave-one-season-out {name} vs basin-month climatology (90 % CI over seasons), "
          f"seasons {seas.min()}-{seas.max()} (n = {len(np.unique(seas))})")
        for k, (v, a, b) in sk.items():
            extra = ""
            if target != "cls":
                extra = f"  AUC {roc_auc_score(y, P[k][:, 1]):.3f}"
            else:
                extra = (f"  AUC(rapid deepening) {roc_auc_score(y == 4, P[k][:, 4]):.3f}"
                         f"  AUC(rapid decay) {roc_auc_score(y == 0, P[k][:, 0]):.3f}")
            w(f"  {k:6s} {v:+.3f} [{a:+.3f}, {b:+.3f}]{extra}")
        for a, b in (("hart", "state"), ("full", "hart"), ("full", "state"), ("full", "full_nogust")):
            d, lo, hi = gain_ci(sums, B, a, b)
            w(f"  gain {a} over {b}: {d:+.3f} [{lo:+.3f}, {hi:+.3f}]")
        if target != "cls":
            on = ~S.hf_now.values.astype(bool)
            sc_on = {k: v[on] for k, v in sc.items()}
            sk_on, sums_on, B_on = season_boot(sc_on, seas[on])
            w(f"  onset only (not HF at t, n={on.sum()}, base rate {y[on].mean():.3f}):")
            for k, (v, a, b) in sk_on.items():
                w(f"    {k:6s} {v:+.3f} [{a:+.3f}, {b:+.3f}]  AUC {roc_auc_score(y[on], P[k][on, 1]):.3f}")
            d, lo, hi = gain_ci(sums_on, B_on, "full", "state")
            w(f"    gain full over state: {d:+.3f} [{lo:+.3f}, {hi:+.3f}]")
            rel[target] = reliability(P["full"][:, 1], y)
        else:
            rel["rapid deepening"] = reliability(P["full"][:, 4], (y == 4).astype(int))
            rel["rapid decay"] = reliability(P["full"][:, 0], (y == 0).astype(int))
        categorical(S, P, CUT, target, w)
        # drift: skill of the full model by decade of the held-out season
        dec = (seas // 10) * 10
        parts = []
        for d0 in np.unique(dec):
            m = dec == d0
            parts.append(f"{d0}s {1 - sc['full'][m].sum() / sc['clim'][m].sum():+.3f}")
        w("  full model skill by decade: " + ", ".join(parts))
        w("")
        cols = SETS["full"]
        ct, md = coef_table(Dt, target, cols, labels if target == "cls" else ["yes"])
        coefs.append(ct)
        models[target] = md
    # train-early / test-late drift check
    w("== Drift: trained 2004-2014, tested 2015-2025")
    for target in ("cls", "hf24"):
        S = D[D.cls >= 0] if target == "cls" else D
        cut = 2014
        tr, te = S[S.season <= cut], S[S.season > cut]
        if not len(tr) or not len(te):
            continue
        cols = SETS["full"]
        prep = Prep(tr[cols].values.astype(float))
        m = fit(prep(tr[cols].values.astype(float)), tr[target].astype(int).values)
        p = m.predict_proba(prep(te[cols].values.astype(float)))
        pc = clim_probs(tr, te.reset_index(drop=True), target, 5 if target == "cls" else 2)
        y = te[target].astype(int).values
        if target == "cls":
            s1 = 1 - rps(p, y).sum() / rps(pc, y).sum()
        else:
            s1 = 1 - brier(p[:, 1], y).sum() / brier(pc[:, 1], y).sum()
        w(f"  {target}: early-trained skill on {cut + 1}-2025 {s1:+.3f} (compare the leave-one-season-out decades above)")
    w("")
    w("== Reliability of the full model (leave-one-season-out): bin, n, mean forecast, observed frequency")
    for k, r in rel.items():
        w(f"  {k}")
        for lo, hi, n, pf, ob in r:
            w(f"    {lo:.2f}-{min(hi, 1):.2f}  n={n:6d}  fc {pf:.3f}  obs {ob:.3f}")
    open(os.path.join(outdir, "skill.txt"), "w").write("\n".join(rep) + "\n")
    pd.concat(coefs).round(4).to_csv(os.path.join(outdir, "coefficients.csv"), index=False)
    json.dump(models, open(os.path.join(outdir, "model.json"), "w"), indent=1)
    print("\n".join(rep))


if __name__ == "__main__":
    main(*sys.argv[1:4])
