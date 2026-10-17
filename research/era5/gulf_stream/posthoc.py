"""POST HOC diagnostics (not pre-registered). Usage: python3 posthoc.py work/stationtimes.csv.gz results/"""
import sys, numpy as np, pandas as pd
import samples as S, models as M
from stats_lib import Boot, pval
rng = np.random.default_rng(7); d0 = pd.read_csv(sys.argv[1], dtype={"station": str}); w = S.wind_sample(d0); L = []
def P(*a): s = " ".join(str(x) for x in a); print(s); L.append(s)
def contrast(wd, sts, nb=2000):
    d, X, y, k = M.wind_xy(wd, sts); B = Boot(X, y, d.winter.values); f = M.wind_contrast(k, "high")
    th = f(B.fit()); dr = B.boot(f, nb, rng); lo, hi = np.nanpercentile(dr, [2.5, 97.5])
    return len(d), 100 * (np.exp(th) - 1), 100 * (np.exp(lo) - 1), 100 * (np.exp(hi) - 1), pval(dr)
P("# POST HOC diagnostics of the wind result (not pre-registered; leads only)")
P("\nD1 one buoy at a time (% wind, high dT vs 0-3 C):")
for s in S.STREAM + S.SHELF:
    try: n, e, lo, hi, p = contrast(w, [s]); P(f"{s}: n={n} {e:+.1f}% ({lo:+.1f}..{hi:+.1f})")
    except Exception as ex: P(s, "failed", ex)
P("\nD2 NW flow only (wind from 270-360), by group:")
nw = w[(w.wdir >= 270) & (w.wdir <= 360)]
for lab, sts in [("Stream", S.STREAM), ("shelf", S.SHELF)]:
    n, e, lo, hi, p = contrast(nw, sts); P(f"{lab}: n={n} {e:+.1f}% ({lo:+.1f}..{hi:+.1f})")
P("\nD3 leave one winter out, Stream group and shelf group (min..max of the % effect):")
for lab, sts in [("Stream", S.STREAM), ("shelf", S.SHELF)]:
    d, X, y, k = M.wind_xy(w, sts); B = Boot(X, y, d.winter.values); f = M.wind_contrast(k, "high"); r = []
    for i in range(len(B.w)):
        c = np.ones(len(B.w)); c[i] = 0; r.append(100 * (np.exp(f(B.fit(c))) - 1))
    P(f"{lab}: {min(r):+.1f}% .. {max(r):+.1f}%  ({len(B.w)} winters)")
P("\nD4 placebo: shuffle dT among station-times of the same station, winter and Vg tercile (expect about 0):")
for lab, sts in [("Stream", S.STREAM), ("shelf", S.SHELF)]:
    e = []
    for _ in range(20):
        x = w[w.station.isin(sts)].copy(); x["vt"] = x.groupby("station").vg.transform(lambda v: pd.qcut(v, 3, labels=False, duplicates="drop"))
        x["dT"] = x.groupby(["station", "winter", "vt"]).dT.transform(lambda v: rng.permutation(v.values))
        e.append(contrast(x, sts, nb=10)[1])
    P(f"{lab}: mean {np.mean(e):+.2f}%  range {min(e):+.1f}..{max(e):+.1f}")
P("\nD5 how well the Vg control works: slope of log wind on log Vg, with dT bins; raw median ratio V10/Vg by Vg class, 0-3 C bin:")
for lo, hi in [(15, 25), (25, 35), (35, 50), (50, 200)]:
    a = w[(w.vg >= lo) & (w.vg < hi) & (w.dT >= 0) & (w.dT < 3)]; P(f"Vg {lo}-{hi}: n={len(a)} median V10/Vg = {(a.v10/a.vg).median():.2f}")
open(sys.argv[2] + "/posthoc-result.txt", "w").write("\n".join(L) + "\n")
