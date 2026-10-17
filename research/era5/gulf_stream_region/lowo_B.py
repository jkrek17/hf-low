"""Pre-registered B rule: with 5 winters, B is 'yes' only if the 95% interval excludes 0 in 5 of 5 leave-one-winter-out refits."""
import sys, numpy as np, pandas as pd
from region_lib import *
d_all = prep_B(pd.read_csv(sys.argv[1])); out = []; rng = np.random.default_rng(7)
for g, sel in {"B1 all boxes": lambda d: d.v_ascat > 0, "B2 in band": lambda d: d.band == "in", "B3 outside band": lambda d: d.band == "out"}.items():
    s = d_all[sel(d_all)]; res = []
    for w in sorted(s.winter.unique()):
        d = s[s.winter != w]; B, k = fit_boot(d, False); f = contrast_high(k); th = f(B.fit()); dr = B.boot(f, 4000, rng)
        lo, hi = np.nanpercentile(dr, [2.5, 97.5]); res.append((w, pct(th), pct(lo), pct(hi)))
    n_ex = sum(1 for r in res if r[2] > 0)
    out.append(f"{g}: interval excludes 0 in {n_ex} of {len(res)} refits")
    out += [f"   without winter {w}: {t:+.2f}% ({lo:+.2f} to {hi:+.2f})" for w, t, lo, hi in res]
open(sys.argv[2], "w").write("\n".join(out) + "\n"); print("\n".join(out))
