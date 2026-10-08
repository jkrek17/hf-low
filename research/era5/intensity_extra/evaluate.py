"""Added skill of candidate predictors over the PR 12 full model (see PREREG.md).

usage: evaluate.py FEATURE_DIR OUTDIR [stage]
  stage groups   (default) leave-one-season-out for base, every Tier 1 group, all1
  stage nested   sel1 and sel1_nested, only if a group met the "helps" rule

Everything is ERA5 proxy, pipeline A, seasons 2004-05..2021-22 (18 seasons).
Gain = BSS(candidate) - BSS(base), both against basin x month climatology,
from pooled leave-one-season-out Brier scores. 90 % intervals resample whole
seasons. p = exact two-sided sign-flip test on the 18 per-season Brier
differences. q = Benjamini-Hochberg.
"""
import os, sys, glob, json, itertools
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np, pandas as pd
from joblib import Parallel, delayed

HERE = os.path.dirname(os.path.abspath(__file__))
INT = os.path.join(HERE, "..", "intensity")
sys.path.insert(0, INT)
import model as M

LAST_SEASON = 2021
NBOOT = 2000
GROUPS = {
    "H history": ["dg12", "dpacc"],
    "A air-sea": ["airsea"],
    "M moisture transport": ["ivt500", "ivtmax"],
    "D diabatic": ["precip6"],
    "J jet geometry": ["jet_along", "jet_cross"],
    "L low-level wind": ["ws850"],
    "P pressure environment": ["dphigh", "ringgrad"],
    "K land": ["landfrac"],
}
GROUPS2 = {"N stability": ["stab"], "R trough": ["trough"], "W warm conveyor": ["omega700"]}
BASE = M.SETS["full"]
RD = 4                                      # rapid deepening class index in model.py
THRESH = 0.005


def history(D):
    """dg12 and dpacc from the track's own earlier fixes (0 where not available)."""
    D = D.sort_values(["track", "time"])
    t = pd.to_datetime(D.time, format="%Y%m%d%H")
    out = {}
    for lag in (1, 2):
        same = (D.track.values == np.roll(D.track.values, lag)) & \
               ((t - np.roll(t.values, lag)).dt.total_seconds().values == 12 * 3600 * lag)
        same[:lag] = False
        out[lag] = (same, np.roll(D.msl.values, lag), np.roll(D.g800.values, lag))
    s1, m1, g1 = out[1]
    s2, m2, _ = out[2]
    D["dg12"] = np.where(s1, D.g800.values - g1, 0.0)
    prev = np.where(s1, m1 - m2, np.nan)
    now = D.msl.values - m1
    D["dpacc"] = np.where(s1 & s2, now - prev, 0.0)
    return D.sort_index()


def build(featdir, fixes, env):
    D = M.load(fixes, env)
    D = D[D.season <= LAST_SEASON].reset_index(drop=True)
    D = history(D)
    X = pd.concat([pd.read_csv(f, dtype={"time": str}) for f in sorted(glob.glob(f"{featdir}/wb2_*.csv"))])
    D = D.merge(X, on=["track", "time"], how="left", validate="one_to_one")
    D["hf24"] = D.hf24.astype(int)
    D["rd"] = np.where(D.cls >= 0, (D.cls == RD).astype(int), -1)
    return D.reset_index(drop=True)


def fold(D, te_season, models, target, rows_mask):
    """Fit every model on the other seasons, return held-out probabilities and training-count-matched cuts."""
    S = D[rows_mask]
    te = (S.season == te_season).values
    Dtr, Dte = S[~te], S[te]
    y = Dtr[target].astype(int).values
    pc = M.clim_probs(Dtr, Dte.reset_index(drop=True).set_index(np.arange(te.sum())), target, 2)[:, 1]
    out = {"clim": pc}
    cut = {}
    for name, cols in models.items():
        prep = M.Prep(Dtr[cols].values.astype(float))
        m = M.fit(prep(Dtr[cols].values.astype(float)), y)
        out[name] = m.predict_proba(prep(Dte[cols].values.astype(float)))[:, 1]
        ptr = m.predict_proba(prep(Dtr[cols].values.astype(float)))[:, 1]
        cut[name] = np.quantile(ptr, 1 - y.mean())
    return te_season, np.where(te)[0], out, cut


def loso(D, models, target, rows_mask, njobs=4):
    S = D[rows_mask]
    seasons = sorted(S.season.unique())
    res = Parallel(n_jobs=njobs)(delayed(fold)(D, s, models, target, rows_mask) for s in seasons)
    P = {k: np.zeros(len(S)) for k in ["clim"] + list(models)}
    C = {k: np.zeros(len(S)) for k in models}
    for s, ix, out, cut in res:
        for k, v in out.items():
            P[k][ix] = v
        for k, v in cut.items():
            C[k][ix] = v
    return S.reset_index(drop=True), P, C


def season_sums(brier, seasons, us):
    return np.array([brier[seasons == s].sum() for s in us])


_SIGNS = {}


def signflip_p(d):
    n = len(d)
    if n not in _SIGNS:
        _SIGNS[n] = np.array(list(itertools.product((-1, 1), repeat=n)), dtype=np.float32)
    stat = np.abs(_SIGNS[n] @ d.astype(np.float32))
    return float(np.mean(stat >= abs(d.sum()) - 1e-12))


def bh(p):
    p = np.asarray(p)
    o = np.argsort(p)
    q = np.empty(len(p))
    r = p[o] * len(p) / (np.arange(len(p)) + 1)
    q[o] = np.minimum.accumulate(r[::-1])[::-1].clip(max=1)
    return q


def gains(S, P, C, y, target_name, names, mask=None):
    """Gain over base for every name in `names`; mask restricts rows (onset only)."""
    if mask is None:
        mask = np.ones(len(y), bool)
    seas = S.season.values[mask]
    us = np.unique(seas)
    yy = y[mask]
    B = {k: (P[k][mask] - yy) ** 2 for k in P}
    sums = {k: season_sums(v, seas, us) for k, v in B.items()}
    rng = np.random.default_rng(0)
    idx = rng.integers(0, len(us), (NBOOT, len(us)))
    den = sums["clim"].sum()
    bss = {k: 1 - sums[k].sum() / den for k in P}
    rows = []
    for k in names:
        d = sums["base"] - sums[k]
        g = d.sum() / den
        gb = d[idx].sum(1) / sums["clim"][idx].sum(1)
        lo, hi = np.percentile(gb, [5, 95])
        rows.append(dict(target=target_name, model=k, bss_base=bss["base"], bss=bss[k], gain=g, lo90=lo, hi90=hi,
                         p=signflip_p(d), seasons_better=int((d > 0).sum()), n_seasons=len(us)))
    return pd.DataFrame(rows), bss


def hss_gain(S, P, C, y, name, rows=None, nb=1000):
    """Yes/no HSS at the training-count-matched cut: change from base, season-bootstrap 90 % interval."""
    ob = y == 1
    seas = S.season.values
    us = np.unique(seas)

    def h(fc, o):
        return M.table(fc, o)["hss"]
    fb, fk = P["base"] >= C["base"], P[name] >= C[name]
    d = h(fk, ob) - h(fb, ob)
    rng = np.random.default_rng(1)
    ix = {s: np.where(seas == s)[0] for s in us}
    out = []
    for _ in range(nb):
        pick = rng.choice(us, len(us))
        ii = np.concatenate([ix[s] for s in pick])
        out.append(h(fk[ii], ob[ii]) - h(fb[ii], ob[ii]))
    return d, *np.percentile(out, [5, 95]), h(fb, ob), h(fk, ob)


def run_groups(D, outdir):
    models = {"base": BASE}
    for g, cols in GROUPS.items():
        models[g] = BASE + cols
    models["all1"] = BASE + [c for cols in GROUPS.values() for c in cols]
    names = [k for k in models if k != "base"]
    tables, hss, bybasin = [], [], []
    store = {}
    for target, tname, mask in (("hf24", "hf24", np.ones(len(D), bool)), ("rd", "rapid deepening", (D.rd >= 0).values)):
        S, P, C = loso(D, models, target, mask)
        y = S[target].astype(int).values
        store[tname] = (S, P, C, y)
        t, bss = gains(S, P, C, y, tname, names)
        tables.append(t)
        if target == "hf24":
            on = ~S.hf_now.values.astype(bool)
            t2, _ = gains(S, P, C, y, "hf24 onset", names, mask=on)
            tables.append(t2)
        for k in names:
            d, lo, hi, hb, hk = hss_gain(S, P, C, y, k)
            hss.append(dict(target=tname, model=k, hss_base=hb, hss=hk, d_hss=d, lo90=lo, hi90=hi))
        if target == "hf24":
            for b in ("atl", "pac"):
                mb = (S.basin == b).values
                tb, _ = gains(S, P, C, y, f"hf24 {b}", names, mask=mb)
                bybasin.append(tb)
    T = pd.concat(tables, ignore_index=True)
    T["q_target"] = T.groupby("target").p.transform(lambda s: bh(s.values))
    fam = T[T.target.isin(["hf24", "hf24 onset", "rapid deepening"]) & (T.model != "all1")]
    T["q_family"] = np.nan
    T.loc[fam.index, "q_family"] = bh(fam.p.values)
    T.to_csv(f"{outdir}/gains.csv", index=False, float_format="%.5g")
    pd.DataFrame(hss).to_csv(f"{outdir}/hss.csv", index=False, float_format="%.4g")
    pd.concat(bybasin).to_csv(f"{outdir}/gains_by_basin.csv", index=False, float_format="%.5g")
    S, P, C, y = store["hf24"]
    np.savez_compressed(f"{outdir}/hf24_probs.npz", track=S.track.values, time=S.time.values.astype(str),
                        season=S.season.values, y=y, **{k.replace(" ", "_"): v for k, v in P.items()})
    return T, store


def helps(T):
    p = T[(T.target == "hf24") & (T.model != "all1")]
    return p[(p.gain >= THRESH) & (p.lo90 > 0) & (p.q_target < 0.05)].model.tolist()


def report(T, outdir, extra=""):
    with open(f"{outdir}/skill.txt", "w") as f:
        for tg in ("hf24", "hf24 onset", "rapid deepening"):
            t = T[T.target == tg]
            f.write(f"== {tg}: base BSS {t.bss_base.iloc[0]:+.4f}; gain = BSS(cand) - BSS(base), 90% CI over 18 seasons\n")
            for r in t.itertuples():
                f.write(f"  {r.model:24s} gain {r.gain:+.4f} [{r.lo90:+.4f}, {r.hi90:+.4f}]  BSS {r.bss:+.4f}  "
                        f"p {r.p:.4f}  q(target) {r.q_target:.3f}  q(family) {r.q_family:.3f}  seasons better {r.seasons_better}/{r.n_seasons}\n")
        f.write(extra)
    print(open(f"{outdir}/skill.txt").read())


if __name__ == "__main__":
    featdir, outdir = sys.argv[1], sys.argv[2]
    stage = sys.argv[3] if len(sys.argv) > 3 else "groups"
    os.makedirs(outdir, exist_ok=True)
    D = build(featdir, os.path.join(INT, "results/fixes_2004.csv.gz"), os.path.join(INT, "results/env_2004.csv.gz"))
    D.to_pickle(os.path.join(os.environ.get("SCRATCH", "/tmp"), "D.pkl"))
    print(len(D), "fixes", D.track.nunique(), "tracks", D.season.nunique(), "seasons; hf24", int(D.hf24.sum()),
          "; missing:", {c: int(D[c].isna().sum()) for c in sum(GROUPS.values(), [])}, flush=True)
    if stage == "groups":
        T, _ = run_groups(D, outdir)
        report(T, outdir, f"\nGroups meeting the pre-registered 'helps' rule: {helps(T)}\n")
