"""Secondary F1 check: GFS native surface stress (UFLX, VFLX, 18-24 h average) -> stress-equivalent neutral wind, compared with ASCAT in the same dT bins.
Usage: python3 stress_check.py boxes2_B.csv.gz OUT.txt"""
import sys, numpy as np, pandas as pd
from round2_lib import *
from bulk import neutral_from_stress, neutral_from_actual
rng = np.random.default_rng(31); B = pd.read_csv(sys.argv[1]).dropna(subset=["v_ascat", "v_gfs", "dT", "dpt", "pres", "sst_mur", "uflx", "vflx"]); B = B[(B.v_gfs >= 5) & (B.v_ascat >= 5)].copy()
tau = np.hypot(B.uflx.values, B.vflx.values)
B["u_tau"] = neutral_from_stress(tau, B.t2m.values, B.dpt.values, B.skt.values, B.pres.values, B.lat.values, 0)
B["u_en"] = neutral_from_actual(B.v_gfs.values, B.t2m.values, B.dpt.values, B.skt.values, B.pres.values, B.lat.values, 0)
ok = np.isfinite(B.u_tau) & (B.u_tau > 3); B = B[ok]; L = [f"stress route sample {len(B)}; median GFS neutral wind from stress / from 10 m wind: {np.median(B.u_tau / B.u_en):.3f} (dT>=6: {np.median((B.u_tau / B.u_en)[B.dT >= 6]):.3f}; 0-3: {np.median((B.u_tau / B.u_en)[(B.dT >= 0) & (B.dT < 3)]):.3f})"]
for nm, m in (("all", B.v_ascat > 0), ("near front (ERA5 band)", B.band == "in"), ("away", B.band == "out")):
    for lab, col in (("from stress (18-24 h mean)", "u_tau"), ("from 10 m wind (instantaneous)", "u_en")):
        d = B[m].assign(y=np.log(B[m].v_ascat / B[m][col])); Bt, k = fit2(d, False); f = contrast_high(k); th = f(Bt.fit()); dr = Bt.boot(f, 10000, rng); lo, hi = np.nanpercentile(dr, [2.5, 97.5])
        L.append(f"B {nm}, GFS neutral wind {lab}: {pct(th):+.2f}% ({pct(lo):+.2f} to {pct(hi):+.2f}), {25*(np.exp(th)-1):+.2f} kt at 25 kt, p={pval(dr):.4f}")
open(sys.argv[2], "w").write("\n".join(L) + "\n"); print("\n".join(L))
