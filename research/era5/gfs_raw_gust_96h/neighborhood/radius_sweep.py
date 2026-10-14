"""Gust search-radius sweep per PREREG.md. usage: radius_sweep.py FIXES_CSV_GZ OUT_JSON"""
import sys, json
import numpy as np, pandas as pd
F = pd.read_csv(sys.argv[1], dtype={"cycle": str}); out = {}
RAD = {"270nm": "g500", "432nm": "g800", "648nm": "g1200"}
SE = sorted(F.season.unique())
def match_T(g, y):
    k = int(y.sum()); gs = np.sort(g)[::-1]
    return float(gs[k]) if 0 <= k < len(gs) else float("-inf")
def score(sel, hit):
    H = int((hit & (sel.Y0 == 1)).sum()); FA = int((hit & (sel.Y0 == 0)).sum()); M = int((~hit & (sel.Y0 == 1)).sum()); CN = len(sel) - H - FA - M
    N = len(sel); e = ((H + FA) * (H + M) + (M + CN) * (FA + CN)) / N
    return dict(n=N, hits=H, fa=FA, miss=M, cn=CN, pod=H / max(H + M, 1), agree=H / max(H + FA, 1), hss=((H + CN) - e) / (N - e))
for rn, col in RAD.items():
    Tm = {(b, s): match_T(*(lambda d: (d[col].values, d.Y0.values))(F[(F.basin == b) & (F.lead == 0) & (F.season != s)])) for b in ("atl", "pac") for s in SE}
    for b in ("atl", "pac"):
        for l in sorted(F.lead.unique()):
            s = F[(F.basin == b) & (F.lead == l)]; T = np.array([Tm[(b, x)] for x in s.season])
            out.setdefault(f"{b}_{int(l)}", {})[rn] = dict(matched=score(s, (s[col].values > T)), kt64=score(s, (s[col] > 64).values), kt72=score(s, (s[col] > 72).values))
cells = {}
for k, v in out.items():
    h = [v[r]["matched"]["hss"] for r in RAD]
    cells[k] = list(RAD)[0] if h[1] - h[0] < 0.02 else (list(RAD)[1] if h[2] - h[1] < 0.02 else list(RAD)[2])
json.dump(dict(scores=out, stop_radius=cells), open(sys.argv[2], "w"), indent=1)
for k, v in out.items(): print(k, {r: (round(x["matched"]["pod"], 2), round(x["matched"]["agree"], 2), round(x["matched"]["hss"], 3)) for r, x in v.items()}, cells[k])
