"""Registered S5: convert ASCAT equivalent-neutral wind to approximate actual 10 m wind (region_lib.convert_en) and refit A1, B1, B2, B3.
Supersedes the first S5 (results/s5-first-version-result.txt), which used the bulk Richardson number itself as the stability parameter."""
import sys, numpy as np, pandas as pd
from region_lib import *
rng = np.random.default_rng(11); out = []
B = prep_B(pd.read_csv(sys.argv[1])); B["y"] = np.log(convert_en(B.v_ascat, B.t2m, B.skt) / B.v_gfs)
A = prep_A(pd.read_csv(sys.argv[3])); A["y"] = np.log(convert_en(A.v_ascat, A.t2m, A.sst))
for g, d, wv, ref in [("A1 all boxes", A, True, 21.6), ("B1 all boxes", B, False, 25.0), ("B2 in band", B[B.band == "in"], False, 25.0), ("B3 outside band", B[B.band == "out"], False, 25.0)]:
    Bt, k = fit_boot(d, wv); f = contrast_high(k); th = f(Bt.fit()); dr = Bt.boot(f, 4000, rng); lo, hi = np.nanpercentile(dr, [2.5, 97.5])
    out.append(f"S5 {g} after equivalent-neutral conversion: {pct(th):+.2f}% ({pct(lo):+.2f} to {pct(hi):+.2f}), p={pval(dr):.4f}; knots at {ref} kt: {ref * (np.exp(th) - 1):+.2f}")
open(sys.argv[2], "w").write("\n".join(out) + "\n"); print("\n".join(out))
