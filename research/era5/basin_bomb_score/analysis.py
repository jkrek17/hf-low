"""T1-T4 of PREREGISTRATION.md. usage: analysis.py   (committed files only)"""
import os, sys, json, datetime as dt
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bomblib as B
rng = np.random.default_rng(19)
NB, NP = 2000, 1000
def per_season_dev(df, y, Xs, pen_cols=None, ridge=0.0):
    """leave-one-season-out held-out deviance per season for each design in Xs"""
    out = {}
    for s in sorted(df.season.unique()):
        tr = (df.season != s).values; te = ~tr
        row = []
        for X, pen in Xs:
            b = B.poisson_fit(X[tr], y[tr], ridge=pen[0], pen_from=pen[1]); row.append(B.dev(y[te], np.exp(np.clip(X[te] @ b, -20, 8))))
        out[s] = row
    return pd.DataFrame(out, index=range(len(Xs))).T
def hf_weekly(basin, df):
    tr = pd.read_csv(os.path.join(B.ERA, "hf_history/results/all_tracks.csv.gz")); tr = tr[(tr.season >= 2004) & (tr.basin == basin) & (tr.gust800_kt >= 71.7)]
    cnt = {}
    for d in pd.to_datetime(tr.peak_time.astype(str), format="%Y%m%d%H").dt.date:
        s = d.year if d.month >= 10 else d.year - 1; k = (d - dt.date(s, 10, 1)).days // 7; cnt[(s, k)] = cnt.get((s, k), 0) + 1
    return np.array([cnt.get((s, k), 0) for s, k in zip(df.season, df.week)], float)
def pcs(df):
    w = pd.read_csv(os.path.join(B.ERA, "hemispheric/results/weekly_table.csv.gz"))
    cols = [f"{v}{i}" for v in "zus" for i in range(1, 11)]
    m = df[["season", "week"]].merge(w[["season", "week"] + cols], on=["season", "week"], how="left")
    X = m[cols].values; return X
def bootstrap_skill(ps, ncols, extra_pairs=()):
    seasons = ps.index.values; out = []
    for _ in range(NB):
        s = rng.choice(seasons, len(seasons)); sub = ps.loc[s].sum().values; out.append(1 - sub / sub[0])
    return np.array(out)
rows, extra = [], {}
for b in ("atl", "pac"):
    df = B.weekly(b); y = df.bomb.values.astype(float); idx = df.idx.values
    # T1
    df["q"] = pd.qcut(idx, 5, labels=False); mm = df.groupby("month").bomb.mean()
    def f(d):
        e = d.month.map(mm); t, bo = d[d.q == 4], d[d.q == 0]
        return (t.bomb.sum() / e[t.index].sum()) / (bo.bomb.sum() / e[bo.index].sum())
    r1 = f(df); bs = []
    for _ in range(NB):
        d = pd.concat([df[df.season == s] for s in rng.choice(sorted(df.season.unique()), 22)], ignore_index=True); bs.append(f(d))
    q5 = df.groupby("q").bomb.mean().round(2).tolist(); wk1 = df.groupby("q").bomb.apply(lambda x: (x > 0).mean()).round(2).tolist()
    rows.append(dict(basin=b, test="T1 month-adj top/bottom fifth", est=r1, lo=np.percentile(bs, 2.5), hi=np.percentile(bs, 97.5), p=(np.sum(np.array(bs) <= 1) + 1) / (NB + 1)))
    extra[b] = dict(bombs_per_week_by_fifth=q5, share_weeks_with_bomb=wk1, total_bombs=int(y.sum()))
    # T2/T3 : designs
    Xb = B.design(df); X0 = B.design(df, [idx]); P = pcs(df)
    # S1 ridge selected by inner season-block CV within each outer fold
    seasons = sorted(df.season.unique()); d_base, d_s0, d_s1 = {}, {}, {}
    best = {}
    for s in seasons:
        tr = (df.season != s).values; te = ~tr
        Ptr = P[tr]; mu_, sd_ = Ptr.mean(0), Ptr.std(0) + 1e-9; Z = (P - mu_) / sd_
        Xs1 = np.hstack([Xb, Z]); pf = Xb.shape[1]
        sc = {}
        for lam in (3.0, 30.0, 300.0):
            tot = 0.0
            for s2 in [x for x in seasons if x != s][::3]:   # every third inner season, for speed
                itr = tr & (df.season != s2).values; ite = (df.season == s2).values
                bb = B.poisson_fit(Xs1[itr], y[itr], ridge=lam, pen_from=pf); tot += B.dev(y[ite], np.exp(np.clip(Xs1[ite] @ bb, -20, 8)))
            sc[lam] = tot
        lam = min(sc, key=sc.get); best[s] = lam
        b1 = B.poisson_fit(Xs1[tr], y[tr], ridge=lam, pen_from=pf)
        bb0 = B.poisson_fit(Xb[tr], y[tr]); bb2 = B.poisson_fit(X0[tr], y[tr])
        d_base[s] = B.dev(y[te], np.exp(np.clip(Xb[te] @ bb0, -20, 8))); d_s0[s] = B.dev(y[te], np.exp(np.clip(X0[te] @ bb2, -20, 8)))
        d_s1[s] = B.dev(y[te], np.exp(np.clip(Xs1[te] @ b1, -20, 8)))
    ps = pd.DataFrame({"base": d_base, "s0": d_s0, "s1": d_s1})
    bs = bootstrap_skill(ps, 3); sk = 1 - ps.sum() / ps.sum().base
    # permutation p for S0 (idx season blocks shuffled)
    blocks = {s: idx[(df.season == s).values] for s in seasons}; ge = 0; obs = sk.s0; db = ps.base.sum()
    for _ in range(NP):
        order = rng.permutation(seasons); new = np.concatenate([blocks[o] for o in order]); Xn = B.design(df, [new]); d1 = 0.0
        for s in seasons:
            tr = (df.season != s).values; te = ~tr; bb = B.poisson_fit(Xn[tr], y[tr]); d1 += B.dev(y[te], np.exp(np.clip(Xn[te] @ bb, -20, 8)))
        ge += (1 - d1 / db) >= obs
    rows.append(dict(basin=b, test="T2 skill S0 (HF-trained index) over month+prev week", est=sk.s0, lo=np.percentile(bs[:, 1], 2.5), hi=np.percentile(bs[:, 1], 97.5), p=(ge + 1) / (NP + 1)))
    rows.append(dict(basin=b, test="T3 skill S1 (bomb-refit)", est=sk.s1, lo=np.percentile(bs[:, 2], 2.5), hi=np.percentile(bs[:, 2], 97.5), p=np.nan))
    dd = bs[:, 2] - bs[:, 1]
    rows.append(dict(basin=b, test="T3 S1 minus S0", est=sk.s1 - sk.s0, lo=np.percentile(dd, 2.5), hi=np.percentile(dd, 97.5), p=np.nan))
    extra[b]["ridge_chosen"] = pd.Series(best).value_counts().to_dict()
    # T4 HF on same weeks
    yh = hf_weekly(b, df); dfh = df.copy(); dfh["bomb"] = yh; dfh["prev"] = dfh.groupby("season").bomb.shift(1).fillna(yh.mean())
    Xbh, X0h = B.design(dfh), B.design(dfh, [idx]); dh = per_season_dev(dfh, yh, [(Xbh, (0, None)), (X0h, (0, None))])
    bh = bootstrap_skill(dh, 2); skh = 1 - dh.sum() / dh.sum().iloc[0]
    rows.append(dict(basin=b, test="T4 skill S0 for HF counts (same weeks)", est=skh.iloc[1], lo=np.percentile(bh[:, 1], 2.5), hi=np.percentile(bh[:, 1], 97.5), p=np.nan))
R = pd.DataFrame(rows)
# BH over T1 and T2 in both basins
m = R.test.str.startswith(("T1", "T2")); p = R.loc[m, "p"].values; o = np.argsort(p); n = len(p); q = np.empty(n); prev = 1.0
for i in range(n - 1, -1, -1): prev = min(prev, p[o[i]] * n / (i + 1)); q[o[i]] = prev
R.loc[m, "q"] = q
R.round(4).to_csv(os.path.join(B.HERE, "results/tests.csv"), index=False)
json.dump(extra, open(os.path.join(B.HERE, "results/extra.json"), "w"), indent=1)
print(R.round(4).to_string()); print(json.dumps(extra))
