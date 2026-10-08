"""RA-11 addendum A: late HF category north of 60N in the archive (Atlantic, 2004-05..2025-26). Plan: PREREGISTRATION.md.
usage: python3 late_onset.py   (writes results/late_onset-result.txt)"""
import math, os
import numpy as np, pandas as pd
import analyse as A

HERE = os.path.dirname(os.path.abspath(__file__))
d = pd.read_csv(os.path.join(A.ARCH, "HF_Data_-_Atl.csv"), dtype=str)
d["ID"] = d.ID.str.strip(); d["Category"] = d.Category.str.strip()
d["date"] = d.date.where(d.date.str.len() == 10, d.date.str[:8] + d.date.str[9:])
d["t"] = pd.to_datetime(d.date, format="%Y%m%d%H", errors="coerce")
d["lat"] = pd.to_numeric(d.Latitude, errors="coerce"); d["p"] = pd.to_numeric(d.Pressure, errors="coerce")
d["season"] = d.ID.str[:4].astype(int)
d = d[d.season.between(A.FIRST, A.LAST) & d.t.notna()]
rows = []
for k, g in d.groupby("ID"):
    hf = g[g.Category == "HF"].sort_values("t")
    gp = g.dropna(subset=["p"])
    if hf.empty or gp.empty: continue
    tmin = gp.sort_values(["p", "t"]).iloc[0].t
    lag = (hf.iloc[0].t - tmin).total_seconds() / 3600
    rows.append(dict(ID=k, season=int(k[:4]), lat0=hf.iloc[0].lat, lag=lag))
e = pd.DataFrame(rows)
e["north"] = (e.lat0 >= 60).astype(int)
buf = []
def P(s): print(s); buf.append(s)
P("RA-11 addendum A: archive Atlantic events, first HF fix >=6 h after the lowest-pressure fix (late HF category)")
P("events %d, north of 60N %d, south %d; lag mean %.1f h; same-fix or earlier %.1f%%, >=6 h %.1f%%, >=12 h %.1f%%" % (
    len(e), e.north.sum(), (1 - e.north).sum(), e.lag.mean(), 100 * (e.lag < 6).mean(), 100 * (e.lag >= 6).mean(), 100 * (e.lag >= 12).mean()))

def OR(df, thr):
    y = (df.lag >= thr).astype(int).values; x = df.north.values
    a = ((x == 1) & (y == 1)).sum() + .5; b = ((x == 1) & (y == 0)).sum() + .5
    c = ((x == 0) & (y == 1)).sum() + .5; dd = ((x == 0) & (y == 0)).sum() + .5
    return a * dd / (b * c)

groups = {s: e[e.season == s] for s in A.SEASONS}
rng = np.random.default_rng(A.SEED)
res = {}
for thr in (6, 12):
    o = OR(e, thr)
    bs = np.array([OR(pd.concat([groups[A.SEASONS[i]] for i in rng.integers(0, len(A.SEASONS), len(A.SEASONS))]), thr) for _ in range(2000)])
    lo, hi = np.percentile(bs, [2.5, 97.5])
    lb = np.log(bs); se = lb.std()
    z = math.log(o) / se
    p2 = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))
    sh_n = (e[e.north == 1].lag >= thr).mean(); sh_s = (e[e.north == 0].lag >= thr).mean()
    P("  late >=%2d h: share north %.1f%% vs south %.1f%%; OR %.2f [%.2f, %.2f]; log-OR z=%.2f two-sided p=%.4f" % (thr, 100 * sh_n, 100 * sh_s, o, lo, hi, z, p2))
    res[thr] = (o, lo, hi, p2, se)
q = A.bh([res[6][3], res[12][3]])
P("  BH q (2 tests): >=6 h q=%.4f, >=12 h q=%.4f" % tuple(q))
# power: planted OR 3.1 on the observed base rate (south share), north n as observed, season-block bootstrap SE from data
n1, n0 = int(e.north.sum()), int((1 - e.north).sum())
p0 = (e[e.north == 0].lag >= 6).mean()
odds = p0 / (1 - p0) * 3.1; p1 = odds / (1 + odds)
sims = 2000; hit = 0
rng2 = np.random.default_rng(3)
se_data = res[6][4]
for _ in range(sims):
    y1 = rng2.random(n1) < p1; y0 = rng2.random(n0) < p0
    a, b_, c, dd = y1.sum() + .5, (~y1).sum() + .5, y0.sum() + .5, (~y0).sum() + .5
    z = math.log(a * dd / (b_ * c)) / se_data
    hit += z > 1.645
P("  power for planted OR 3.1 (base late share south %.1f%%, north n=%d, south n=%d; SE from the season-block bootstrap): %.2f" % (100 * p0, n1, n0, hit / sims))
P("  for reference PR 86 (ERA5 proxy pipeline A, gust maximum >= 60N, different definition): OR 3.10, late 53% vs 27%")
P("  both thresholds count the same-fix tie as NOT late; archive minimum pressure is the analysed value over rows that may not span the whole life cycle")
os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
open(os.path.join(HERE, "results", "late_onset-result.txt"), "w").write("\n".join(buf) + "\n")
