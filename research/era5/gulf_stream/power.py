"""Planted-effect simulation (pre-registration 'Power'); prints power and SE only, never the data's own coefficients.
Usage: python3 power.py work/stationtimes.csv.gz OUT.txt"""
import sys, numpy as np, pandas as pd
import samples as S, models as M
from stats_lib import Boot, pval
rng = np.random.default_rng(20261017)
d0 = pd.read_csv(sys.argv[1], dtype={"station": str}); OUT = sys.argv[2]
def sim(X, y_base, resid, plant, winter, fn, nout=200, nin=400):
    """y_sim = y_base(fitted without effect terms) + resid(after full fit) + plant; outer: resample winters; inner: bootstrap."""
    ys = y_base + resid + plant; W = np.unique(winter); hits = 0; ests = []
    sub = {w: np.where(winter == w)[0] for w in W}
    for _ in range(nout):
        pick = rng.choice(W, len(W), replace=True)
        idx = np.concatenate([sub[w] for w in pick]); wnew = np.concatenate([np.full(len(sub[w]), i) for i, w in enumerate(pick)])
        B = Boot(X[idx], ys[idx], wnew); est = fn(B.fit()); ests.append(est)
        dr = B.boot(fn, nin, rng); hits += (pval(dr) < 0.05) and est > 0
    return hits / nout, float(np.std(ests))
lines = ["Planted-effect power (winter resampling, nested bootstrap; 200 outer x 400 inner). Unadjusted p < 0.05, correct sign.",
         "Residuals taken after fitting WITH the effect terms (deviation from the plan, logged before any outcome)."]
# --- wind
w = M.S.wind_sample(d0)
for name, sts in [("W1 pooled", S.STREAM + S.SHELF), ("W2 Stream", S.STREAM), ("W3 shelf", S.SHELF)]:
    d, X, y, k = M.wind_xy(w, sts); wn = d.winter.values
    Bf = Boot(X, y, wn); b = Bf.fit(); resid = y - X @ b
    Xn = X[:, : k + 1]; bn = np.linalg.lstsq(Xn, y, rcond=None)[0]; base = Xn @ bn
    hi = (d.dT.values >= 6).astype(float)
    for pct in (0.03, 0.06, 0.10):
        p, se = sim(X, base, resid, np.log(1 + pct) * hi, wn, M.wind_contrast(k, "high"))
        lines.append(f"{name}: planted +{pct*100:.0f}% wind when dT>=6 -> power {p:.2f}  (SE of contrast {se:.4f} log units)")
# W4: extra effect in Stream only
d, X, y, k = M.wind_xy(w, S.STREAM + S.SHELF); wn = d.winter.values
Bf = Boot(X, y, wn); b = Bf.fit(); resid = y - X @ b
Xn = X[:, : k + 1]; base = Xn @ np.linalg.lstsq(Xn, y, rcond=None)[0]
isS = d.station.isin(S.STREAM).values.astype(float)
# interaction design: add stream x bins
Xi = np.column_stack([X] + [X[:, k + 1 + j] * isS for j in range(4)])
bi = Boot(Xi, y, wn).fit(); ri = y - Xi @ bi
f4 = lambda b: 0.5 * (b[k + 5 + 2] + b[k + 5 + 3])   # extra effect of the two high dT bins at Stream stations
for pct in (0.03, 0.06, 0.10):
    p, se = sim(Xi, base, ri, np.log(1 + pct) * (d.dT.values >= 6) * isS, wn, f4)
    lines.append(f"W4 extra Stream effect: planted +{pct*100:.0f}% only at Stream stations -> power {p:.2f} (SE {se:.4f})")
# --- waves
v = S.wave_sample(d0)
d, X, y, k = M.wave_xy(v, S.STREAM); wn = d.winter.values
b = Boot(X, y, wn).fit(); resid = y - X @ b
nb = k + 2 + 1 + 1                       # stations + 2 vbins + dpd + instream
base = X[:, :nb] @ np.linalg.lstsq(X[:, :nb], y, rcond=None)[0]
opp = X[:, nb]; ins = X[:, nb - 1]
f1 = lambda b: b[nb] + b[nb + 1]; f2 = lambda b: b[nb + 1]
for ft in (0.5, 1.0, 2.0):
    p1, s1 = sim(X, base, resid, ft * opp * ins, wn, f1)
    p2, s2 = sim(X, base, resid, ft * opp * ins, wn, f2)
    lines.append(f"V1 planted +{ft} ft opposed in Stream water -> power {p1:.2f} (SE {s1:.2f} ft); V2 (diff-in-diff) -> power {p2:.2f} (SE {s2:.2f} ft)")
open(OUT, "w").write("\n".join(lines) + "\n"); print("\n".join(lines))
