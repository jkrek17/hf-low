"""Pre-registered tests of PREREGISTRATION.md. Usage: python3 analyse.py work/stationtimes.csv.gz results/"""
import sys, os, numpy as np, pandas as pd
import samples as S, models as M
from stats_lib import Boot, pval, bh
rng = np.random.default_rng(20261017); NB = 10000
d0 = pd.read_csv(sys.argv[1], dtype={"station": str}); OUT = sys.argv[2]; os.makedirs(OUT, exist_ok=True)
L = []
def P(*a):
    s = " ".join(str(x) for x in a); print(s); L.append(s)
def ci(x): return np.nanpercentile(x, [2.5, 97.5])
w = S.wind_sample(d0)
vref = {}
for vg in (30, 45):
    r = w[(w.dT >= 0) & (w.dT < 3) & (abs(w.vg - vg) <= 3)]
    vref[vg] = float(r.v10.median()); 
P("# Gulf Stream effects: pre-registered results (buoys + ERA5 pressure-gradient control; ERA5 is a proxy control, not an outcome)")
P(f"Wind sample {len(w)} station-times, {w.winter.nunique()} winters. Reference 0-3 C wind at Vg 30 / 45 kt: {vref[30]:.1f} / {vref[45]:.1f} kt.")
# descriptive table: mean of V10 by dT bin at fixed Vg (observed, no model)
def desc(wd, label):
    P(f"\n## Observed 10 m wind by water-minus-air temperature, {label}")
    P("dT bin (C) | n | median V10 (kt) at Vg 25-35 | n | median V10 (kt) at Vg 40-50")
    for lo, hi in [(-9, 0), (0, 3), (3, 6), (6, 9), (9, 99)]:
        a = wd[(wd.dT >= lo) & (wd.dT < hi)]; x = a[abs(a.vg - 30) <= 5]; y = a[abs(a.vg - 45) <= 5]
        P(f"{lo}..{hi} | {len(x)} | {x.v10.median():.1f} | {len(y)} | {y.v10.median():.1f}")
desc(w, "all 7 buoys"); desc(w[w.station.isin(S.STREAM)], "Stream buoys"); desc(w[w.station.isin(S.SHELF)], "shelf buoys")
# --- family W
def fitgroup(wd, sts, yname="v10"):
    d, X, y, k = M.wind_xy(wd.assign(v10=wd[yname]) if yname != "v10" else wd, sts)
    return k, Boot(X, y, d.winter.values), d
res = {}
def draws_for(wd, yname="v10", sts=None):
    k, B, d = fitgroup(wd, sts, yname); return k, B
kp, Bp = draws_for(w, sts=S.STREAM + S.SHELF); ks, Bs = draws_for(w, sts=S.STREAM); kh, Bh = draws_for(w, sts=S.SHELF)
def boot_pair(B1, f1, B2=None, f2=None):
    out = np.empty(NB); 
    for i in range(NB):
        c = np.bincount(rng.integers(0, len(B1.w), len(B1.w)), minlength=len(B1.w)).astype(float)
        a = f1(B1.fit(c))
        if B2 is not None:
            c2 = np.bincount(rng.integers(0, len(B2.w), len(B2.w)), minlength=len(B2.w)).astype(float)
            a = a - f2(B2.fit(c2))
        out[i] = a
    return out
tests = {}
tests["W1 pooled, high-dT minus 0-3 C"] = (boot_pair(Bp, M.wind_contrast(kp, "high")), Bp.fit(), M.wind_contrast(kp, "high"))
tests["W2 Stream buoys, same contrast"] = (boot_pair(Bs, M.wind_contrast(ks, "high")), Bs.fit(), M.wind_contrast(ks, "high"))
tests["W3 shelf buoys, same contrast"] = (boot_pair(Bh, M.wind_contrast(kh, "high")), Bh.fit(), M.wind_contrast(kh, "high"))
tests["W4 Stream minus shelf"] = (boot_pair(Bs, M.wind_contrast(ks, "high"), Bh, M.wind_contrast(kh, "high")), None, None)
tests["W5 dT>=9 minus 6-9, pooled"] = (boot_pair(Bp, M.wind_contrast(kp, "top_minus_mid")), Bp.fit(), M.wind_contrast(kp, "top_minus_mid"))
est = {"W1 pooled, high-dT minus 0-3 C": M.wind_contrast(kp, "high")(Bp.fit()), "W2 Stream buoys, same contrast": M.wind_contrast(ks, "high")(Bs.fit()),
       "W3 shelf buoys, same contrast": M.wind_contrast(kh, "high")(Bh.fit()),
       "W5 dT>=9 minus 6-9, pooled": M.wind_contrast(kp, "top_minus_mid")(Bp.fit())}
est["W4 Stream minus shelf"] = est["W2 Stream buoys, same contrast"] - est["W3 shelf buoys, same contrast"]
names = list(tests); pv = [pval(tests[n][0]) for n in names]; qv = bh(pv)
P("\n## Family W (5 tests): wind change for water warmer than air, at the same pressure-gradient wind")
P("test | estimate (% wind) | 95% interval | kt at Vg 30 | kt at Vg 45 | p | q (BH, family)")
for n, p, q in zip(names, pv, qv):
    th = est[n]; lo, hi = ci(tests[n][0])
    f = lambda t, v: vref[v] * (np.exp(t) - 1)
    P(f"{n} | {100*(np.exp(th)-1):+.1f} | {100*(np.exp(lo)-1):+.1f} to {100*(np.exp(hi)-1):+.1f} | {f(th,30):+.1f} ({f(lo,30):+.1f} to {f(hi,30):+.1f}) | {f(th,45):+.1f} ({f(lo,45):+.1f} to {f(hi,45):+.1f}) | {p:.4f} | {q:.4f}")
# bins one by one, pooled and groups (descriptive of the model)
P("\n## Model bin effects (% wind relative to 0-3 C, same Vg), pooled / Stream / shelf")
for lab, (k, B) in {"pooled": (kp, Bp), "Stream": (ks, Bs), "shelf": (kh, Bh)}.items():
    b = B.fit(); bd = np.array([B.boot(lambda x, j=j: x[k + 1 + j], 2000, rng) for j in range(4)])
    P(lab, " | ".join(f"{n}: {100*(np.exp(b[k+1+j])-1):+.1f}% ({100*(np.exp(np.nanpercentile(bd[j],2.5))-1):+.1f}..{100*(np.exp(np.nanpercentile(bd[j],97.5))-1):+.1f})" for j, n in enumerate(["<0", "3-6", "6-9", ">=9"])), f"| log Vg slope {b[k]:.2f}")
# --- secondary
P("\n## Secondary (not in family; raw p, labelled secondary)")
def sec(label, wd, yname="v10"):
    k, B, d = fitgroup(wd, S.STREAM + S.SHELF, yname)
    th = M.wind_contrast(k, "high")(B.fit()); dr = B.boot(M.wind_contrast(k, "high"), 3000, rng)
    P(f"{label}: n={len(d)}  {100*(np.exp(th)-1):+.1f}% ({100*(np.exp(ci(dr)[0])-1):+.1f}..{100*(np.exp(ci(dr)[1])-1):+.1f}), p={pval(dr):.4f}")
sec("S1a flow from 270-360 deg", w[(w.wdir >= 270) & (w.wdir <= 360)]); sec("S1b other flow directions", w[~((w.wdir >= 270) & (w.wdir <= 360))])
sec("S2 Dec-Mar only", w[w.month.isin([12, 1, 2, 3])])
sec("S3 buoy gust instead of mean", w.dropna(subset=["gst"]).query("gst>=1").assign(v10=lambda x: x.gst))
sec("S4 raw (not 10 m-adjusted) wind", w.assign(v10=w.wspd_raw))
# --- family V
v = S.wave_sample(d0)
P(f"\n## Waves: sample {len(v)} steady-wind hours (>=20 kt); opposed/following hours in Stream water at Stream buoys:")
q = v[v.station.isin(S.STREAM) & v.sector.isin(["opposed", "following"]) & v.instream]
P(f"n={len(q)} hours, {q.winter.nunique()} winters; opposed {int((q.sector=='opposed').sum())}, following {int((q.sector=='following').sum())}")
d, X, y, k = M.wave_xy(v, S.STREAM); nb = k + 2 + 1 + 1
Bw = Boot(X, y, d.winter.values); bw = Bw.fit()
f1 = lambda b: b[nb] + b[nb + 1]; f2 = lambda b: b[nb + 1]
dsh, Xh, yh, kh2 = M.wave_xy(v, S.SHELF, instream_interaction=False); Bsh = Boot(Xh, yh, dsh.winter.values); f3 = lambda b: b[kh2 + 3]
vt = {"V1 opposed minus following, Stream water, Stream buoys": (f1(bw), Bw.boot(f1, NB, rng)),
      "V2 diff-in-diff vs non-Stream water": (f2(bw), Bw.boot(f2, NB, rng)),
      "V3 placebo, shelf buoys (41002 labels)": (f3(Bsh.fit()), Bsh.boot(f3, NB, rng))}
pv2 = [pval(x[1]) for x in vt.values()]; qv2 = bh(pv2)
P("test | estimate (ft) | 95% interval | p | q (BH, family)")
for (n, (e, dr)), p, qq in zip(vt.items(), pv2, qv2): P(f"{n} | {e:+.2f} | {ci(dr)[0]:+.2f} to {ci(dr)[1]:+.2f} | {p:.4f} | {qq:.4f}")
allq = bh(pv + pv2); P("q across all 8 tests:", " ".join(f"{x:.4f}" for x in allq))
# wave descriptive: mean WVHT in ft by wind bin and sector, Stream water, Stream buoys
P("\n## Observed wave height (ft), Stream buoys in Stream water: mean by wind bin and sector (hours)")
for vb in ["20-30", "30-40", "40+"]:
    r = []
    for sct in ["opposed", "cross", "following"]:
        a = v[v.station.isin(S.STREAM) & v.instream & (v.vbin == vb) & (v.sector == sct)]; r.append(f"{sct} {a.wvht_ft.mean():.1f} ft ({len(a)})")
    P(vb, "kt:", "; ".join(r))
open(f"{OUT}/gulf_stream-result.txt", "w").write("\n".join(L) + "\n")
