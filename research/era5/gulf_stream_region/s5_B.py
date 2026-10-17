"""Registered S5 applied to B: convert ASCAT equivalent-neutral wind to approximate actual 10 m wind (same formula as analyse_region.py S5, GFS t2m and skin temperature) and refit B1 and B2."""
import sys, numpy as np, pandas as pd
from region_lib import *
B = prep_B(pd.read_csv(sys.argv[1])); rng = np.random.default_rng(11); U = B.v_ascat / 1.944
Rb = 9.81 / B.t2m * (B.t2m - B.skt) * 10 / np.maximum(U, 1) ** 2
z = np.where(Rb < 0, Rb, Rb / (1 - 5 * np.clip(Rb, None, 0.15))); x = (1 - 16 * np.clip(z, -2, 0)) ** 0.25
psi = np.where(z < 0, 2 * np.log((1 + x) / 2) + np.log((1 + x * x) / 2) - 2 * np.arctan(x) + np.pi / 2, -5 * z)
B["y"] = np.log(B.v_ascat * (1 - psi / 10.8) / B.v_gfs); out = []
for g, sel in {"B1 all boxes": B.v_ascat > 0, "B2 in band": B.band == "in", "B3 outside band": B.band == "out"}.items():
    d = B[sel]; Bt, k = fit_boot(d, False); f = contrast_high(k); th = f(Bt.fit()); dr = Bt.boot(f, 4000, rng); lo, hi = np.nanpercentile(dr, [2.5, 97.5])
    out.append(f"S5 {g} after equivalent-neutral conversion: {pct(th):+.2f}% ({pct(lo):+.2f} to {pct(hi):+.2f}), p={pval(dr):.4f}; knots at 25 kt model wind {25*(np.exp(th)-1):+.2f}")
open(sys.argv[2], "w").write("\n".join(out) + "\n"); print("\n".join(out))
