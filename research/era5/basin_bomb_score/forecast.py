"""F1-F4 of PREREGISTRATION_FORECAST.md. GEFS-fed (and ERA5 z+u reference) index against weekly bomb and HF counts. Committed files only."""
import os, sys, json, datetime as dt
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bomblib as B
rng = np.random.default_rng(23); NB, NP = 2000, 1000
def hf_weekly(basin, df):
    tr = pd.read_csv(os.path.join(B.ERA, "hf_history/results/all_tracks.csv.gz")); tr = tr[(tr.season >= 2004) & (tr.basin == basin) & (tr.gust800_kt >= 71.7)]
    cnt = {}
    for d in pd.to_datetime(tr.peak_time.astype(str), format="%Y%m%d%H").dt.date:
        s = d.year if d.month >= 10 else d.year - 1; k = (d - dt.date(s, 10, 1)).days // 7; cnt[(s, k)] = cnt.get((s, k), 0) + 1
    return np.array([cnt.get((s, k), 0) for s, k in zip(df.season, df.week)], float)
def season_dev(df, y, idx):
    Xb = B.design(df); X1 = B.design(df, [idx]); out = {}
    for s in sorted(df.season.unique()):
        tr = (df.season != s).values; te = ~tr
        r = []
        for X in (Xb, X1):
            b = B.poisson_fit(X[tr], y[tr]); r.append(B.dev(y[te], np.exp(np.clip(X[te] @ b, -20, 8))))
        out[s] = r
    return pd.DataFrame(out, index=["base", "idx"]).T
def boot_skill(ps):
    ss = ps.index.values; o = []
    for _ in range(NB):
        sub = ps.loc[rng.choice(ss, len(ss))].sum(); o.append(1 - sub.idx / sub.base)
    return np.array(o)
def ratio(df, y, idx, mm):
    q = pd.qcut(idx, 5, labels=False); e = df.month.map(mm).values
    t, b = q == 4, q == 0
    return (y[t].sum() / e[t].sum()) / (y[b].sum() / e[b].sum())
rows = []
for b in ("atl", "pac"):
    base = B.weekly(b); gi = pd.read_csv(os.path.join(B.ERA, "gfs_test/results/index_zu.csv"))
    df = base.drop(columns=["idx"]).merge(gi[["season", "week", f"{b}_gefs", f"{b}_era5"]], on=["season", "week"], how="inner").reset_index(drop=True)
    assert len(df) == 450
    for outcome in ("bomb", "HF"):
        y = df.bomb.values.astype(float) if outcome == "bomb" else hf_weekly(b, df)
        d2 = df.copy(); d2["prev"] = pd.Series(y).groupby(d2.season.values).shift(1).fillna(y.mean()).values
        mm = pd.Series(y).groupby(d2.month.values).mean(); seasons = sorted(d2.season.unique())
        res = {}
        for name, col in (("GEFS", f"{b}_gefs"), ("ERA5zu", f"{b}_era5")):
            idx = d2[col].values; ps = season_dev(d2, y, idx); sk = 1 - ps.sum().idx / ps.sum().base; bs = boot_skill(ps)
            r1 = ratio(d2, y, idx, mm); rb = []
            for _ in range(NB):
                pick = rng.choice(seasons, len(seasons)); sub = pd.concat([d2[d2.season == s] for s in pick], ignore_index=True)
                ysub = np.concatenate([y[(d2.season == s).values] for s in pick]); rb.append(ratio(sub, ysub, sub[col].values, mm))
            p = None
            if name == "GEFS":
                blocks = {s: idx[(d2.season == s).values] for s in seasons}; ge = 0; db = ps.sum().base
                for _ in range(NP):
                    new = np.concatenate([blocks[o] for o in rng.permutation(seasons)]); Xn = B.design(d2, [new]); d1 = 0.0
                    for s in seasons:
                        tr = (d2.season != s).values; te = ~tr; bb = B.poisson_fit(Xn[tr], y[tr]); d1 += B.dev(y[te], np.exp(np.clip(Xn[te] @ bb, -20, 8)))
                    ge += (1 - d1 / db) >= sk
                p = (ge + 1) / (NP + 1)
            res[name] = (sk, bs, ps)
            rows.append(dict(basin=b, outcome=outcome, index=name, ratio=r1, ratio_lo=np.percentile(rb, 2.5), ratio_hi=np.percentile(rb, 97.5),
                             skill=sk, skill_lo=np.percentile(bs, 2.5), skill_hi=np.percentile(bs, 97.5), p_perm=p))
            # F4: 2015-16..2018-19 only, trained on the rest (leave-one-season-out within those four held out)
            late = ps.loc[[s for s in seasons if s >= 2015]].sum(); sk4 = 1 - late.idx / late.base
            bs4 = boot_skill(ps.loc[[s for s in seasons if s >= 2015]])
            rows.append(dict(basin=b, outcome=outcome, index=name + " F4 2015-18", skill=sk4, skill_lo=np.percentile(bs4, 2.5), skill_hi=np.percentile(bs4, 97.5)))
        # F3 difference
        dg, de = res["GEFS"][2], res["ERA5zu"][2]; diff = []
        for _ in range(NB):
            s = rng.choice(seasons, len(seasons)); a, c = dg.loc[s].sum(), de.loc[s].sum(); diff.append((1 - a.idx / a.base) - (1 - c.idx / c.base))
        rows.append(dict(basin=b, outcome=outcome, index="GEFS minus ERA5zu", skill=res["GEFS"][0] - res["ERA5zu"][0], skill_lo=np.percentile(diff, 2.5), skill_hi=np.percentile(diff, 97.5)))
R = pd.DataFrame(rows)
m = (R["index"] == "GEFS"); pv = R.loc[m, "p_perm"].values; o = np.argsort(pv); n = len(pv); q = np.empty(n); prev = 1.0
for i in range(n - 1, -1, -1): prev = min(prev, pv[o[i]] * n / (i + 1)); q[o[i]] = prev
R.loc[m, "q"] = q
R.round(4).to_csv(os.path.join(B.HERE, "results/forecast_tests.csv"), index=False); print(R.round(3).to_string())
