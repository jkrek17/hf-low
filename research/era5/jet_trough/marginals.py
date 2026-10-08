"""Plain-words tables from the held-out seasons (2015-16..2025-26): observed frequency of bomb and HF
by jet-speed and trough-depth tercile and by jet threshold, with the framework's own forecast alongside.
Terciles are from the fit seasons (2004-05..2014-15), per basin. ERA5 proxy, pipeline A. Outcome-free
choices only: no cut here was picked from the held-out outcomes.

  bomb      ndr24 >= +1 Bergeron, class-eligible fixes
  hf24      gust index reaches 71.7 kt within 24 h, all fixes (the framework's target; includes fixes already HF)
  hfon      the same among fixes not already HF
  fw_*      framework `full` model (L2, C = 1) fitted on the fit seasons, mean forecast probability
"""
import numpy as np, pandas as pd
import jtlib as J

D = J.load_all()
S = {}
for tgt, rows in (("hf24", D.index), ("bomb", D.index[D.elig])):
    X = D.loc[rows]; fm = X.fit.values
    cols = J.FW.SETS["full"]
    prep = J.FW.Prep(X.loc[fm, cols].values.astype(float))
    m = J.FW.fit(prep(X.loc[fm, cols].values.astype(float)), X.loc[fm, tgt].astype(int).values)
    S[tgt] = pd.Series(m.predict_proba(prep(X.loc[~fm, cols].values.astype(float)))[:, 1], index=X.index[~fm])
D["fw_hf24"] = S["hf24"].reindex(D.index); D["fw_bomb"] = S["bomb"].reindex(D.index)
rows = []
def summ(g, **k):
    e = g[g.elig]; h = g[~g.hf_now]
    rows.append(dict(**k, n=len(g), n_elig=len(e), n_bomb=int(e.bomb.sum()), bomb=e.bomb.mean(), fw_bomb=e.fw_bomb.mean(),
                     n_hf24=int(g.hf24.sum()), hf24=g.hf24.mean(), fw_hf24=g.fw_hf24.mean(),
                     n_hfon=int(h.hf24.sum()), hfon=h.hf24.mean()))
for basin in ("atl", "pac"):
    B = D[D.basin == basin]; fit = B[B.fit]; te = B[~B.fit]
    qj = np.percentile(fit.jet250, [100 / 3, 200 / 3]); qt = np.percentile(fit.trough_up, [100 / 3, 200 / 3])
    for a in range(3):
        g = te[np.searchsorted(qj, te.jet250) == a]
        summ(g, basin=basin, cut="jet tercile", level=a + 1, lo=(np.nan if a == 0 else qj[a - 1]), hi=(np.nan if a == 2 else qj[a]))
    for a in range(3):
        g = te[np.searchsorted(qt, te.trough_up) == a]
        summ(g, basin=basin, cut="trough tercile", level=a + 1, lo=(np.nan if a == 0 else qt[a - 1]), hi=(np.nan if a == 2 else qt[a]))
    for thr in (80, 90, 100, 110, 120, 130, 140, 150, 160):
        summ(te[te.jet250 >= thr], basin=basin, cut="jet >= kt", level=thr, lo=thr, hi=np.nan)
        summ(te[te.jet250 < thr], basin=basin, cut="jet < kt", level=thr, lo=np.nan, hi=thr)
    sd = np.clip(fit.jet250, *np.percentile(fit.jet250, [1, 99])).std()
    rows.append(dict(basin=basin, cut="jet SD (fit, clipped 1-99 pct)", level=np.nan, lo=sd))
    sd = np.clip(fit.trough_up, *np.percentile(fit.trough_up, [1, 99])).std()
    rows.append(dict(basin=basin, cut="trough SD (fit, clipped 1-99 pct)", level=np.nan, lo=sd))
M = pd.DataFrame(rows); M.to_csv("results/marginals.csv", index=False, float_format="%.5g")
with open("results/marginals.txt", "w") as f:
    f.write("Held-out seasons 2015-16..2025-26; terciles from fit seasons 2004-05..2014-15; ERA5 proxy, pipeline A.\n")
    f.write(M.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
print(open("results/marginals.txt").read())
