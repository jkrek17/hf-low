"""Gale/storm/HF track counts by PR 41 pattern-index quintile. ERA5 proxy, pipeline A. See PREREGISTRATION.md.
usage: quintiles.py   (reads committed files only)"""
import os, json, datetime as dt
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
ERA = os.path.dirname(HERE)
tr = pd.read_csv(os.path.join(ERA, "hf_history/results/all_tracks.csv.gz"))
tr = tr[tr.season >= 2004].copy()
pt = pd.to_datetime(tr.peak_time.astype(str), format="%Y%m%d%H")
tr["date"] = pt.dt.date
CL = {"gale": (40, 54), "storm": (54, 71.7), "HF": (71.7, 1e9)}
rng = np.random.default_rng(41)
rows, tab = [], []
for b in ("atl", "pac"):
    ix = pd.read_csv(os.path.join(ERA, f"hemispheric/results/oos_index_{b}.csv"))
    ix = ix[ix.season >= 2004].copy()
    ix["start"] = [dt.date(s, 10, 1) + dt.timedelta(days=7 * k) for s, k in zip(ix.season, ix.week)]
    ix["month"] = [(d + dt.timedelta(days=3)).month for d in ix.start]
    t = tr[tr.basin == b]
    for c, (lo, hi) in CL.items():
        tc = t[(t.gust800_kt >= lo) & (t.gust800_kt < hi)]
        cnt = {}
        for d in tc.date:
            k = (d - dt.date(d.year if d.month >= 10 else d.year - 1, 10, 1)).days // 7
            s = d.year if d.month >= 10 else d.year - 1
            cnt[(s, k)] = cnt.get((s, k), 0) + 1
        ix[c] = [cnt.get((s, k), 0) for s, k in zip(ix.season, ix.week)]
    ix["all"] = ix[list(CL)].sum(axis=1)
    # total tracks of any class in window for shares
    ix["q"] = pd.qcut(ix.idx, 5, labels=False)
    seasons = sorted(ix.season.unique())
    def stats(d):
        out = {}
        for c in list(CL):
            m = d.groupby("month")[c].transform("mean")  # placeholder, replaced below
        return out
    def ratios(d, c):
        mm = ix.groupby("month")[c].mean()  # month mean over all weeks (fixed)
        exp = d.month.map(mm)
        top, bot = d[d.q == 4], d[d.q == 0]
        r = top[c].mean() / bot[c].mean()
        ra = (top[c].sum() / exp[top.index].sum()) / (bot[c].sum() / exp[bot.index].sum())
        return r, ra
    for c in CL:
        pt_ = ratios(ix, c)
        bs = []
        for _ in range(2000):
            samp = rng.choice(seasons, len(seasons))
            d = pd.concat([ix[ix.season == s] for s in samp], ignore_index=True)
            try: bs.append(ratios(d, c))
            except Exception: pass
        bs = np.array(bs)
        p = (np.sum(bs[:, 1] <= 1) + 1) / (len(bs) + 1)
        rows.append(dict(basin=b, cls=c, ratio=pt_[0], ratio_lo=np.percentile(bs[:, 0], 2.5), ratio_hi=np.percentile(bs[:, 0], 97.5),
                         ratio_month_adj=pt_[1], adj_lo=np.percentile(bs[:, 1], 2.5), adj_hi=np.percentile(bs[:, 1], 97.5), p_one_sided=p))
    for qq in range(5):
        d = ix[ix.q == qq]
        r = dict(basin=b, quintile=qq + 1, weeks=len(d), idx_mean=d.idx.mean())
        for c in CL: r[c + "_per_week"] = d[c].mean()
        r["tracks_per_week"] = d["all"].mean()
        # share of weeks with at least one
        for c in CL: r[c + "_weeks_with_one"] = (d[c] > 0).mean()
        tab.append(r)
R = pd.DataFrame(rows)
# BH over the 6 month-adjusted ratios
p = R.p_one_sided.values; o = np.argsort(p); m = len(p)
q = np.empty(m); prev = 1
for i in range(m - 1, -1, -1):
    prev = min(prev, p[o[i]] * m / (i + 1)); q[o[i]] = prev
R["q"] = q
def verdict(r):
    if r.adj_lo > 1 and r.ratio_month_adj >= 1.5 and r.q < 0.05: return "separates"
    if r.adj_lo > 1: return "weak"
    return "can't tell"
R["verdict"] = R.apply(verdict, axis=1)
R.round(3).to_csv(os.path.join(HERE, "results/ratios.csv"), index=False)
pd.DataFrame(tab).round(3).to_csv(os.path.join(HERE, "results/by_quintile.csv"), index=False)
print(R.round(2).to_string()); print(pd.DataFrame(tab).round(2).to_string())
