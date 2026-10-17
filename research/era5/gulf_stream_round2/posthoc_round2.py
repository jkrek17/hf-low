"""POST HOC (not in the plan). Is the near-minus-away gap in A a latitude artefact? Compare near-front and away boxes inside the same 4-degree latitude band, so the Coriolis parameter (and the yardstick) is similar.
Usage: python3 posthoc_round2.py boxes2_A.csv.gz OUT.txt"""
import sys, numpy as np, pandas as pd
from round2_lib import *
rng = np.random.default_rng(51)
A = pd.read_csv(sys.argv[1]).dropna(subset=["v_ascat", "vg", "dT", "tadv", "curv", "d2m", "sst_mur", "msl_hpa"]); A = A[(A.vg >= 15) & (A.v_ascat >= 5)].copy()
A["lvg"] = np.log(A.vg); A["dTm"] = (A.sst_mur - A.t2m).clip(-3, 12); A["y"] = np.log(A.v_ascat)
A["tadv_bin"] = tadv_bin(A.tadv.values); q = A.curv.quantile([1 / 3, 2 / 3]).values; A["curv_bin"] = np.select([A.curv <= q[0], A.curv <= q[1]], ["lo", "mid"], "hi")
L = []
for spec, dT, bc, ex in (("round 1 spec (ERA5 front, no controls)", "dT", "band", ()), ("MUR front and dT, with air-mass controls", "dTm", "band_mur", ("tadv_bin", "curv_bin"))):
    for lo, hi in ((29, 33), (33, 37), (37, 41)):
        d = A[(A.lat >= lo) & (A.lat < hi)]; ni, no = (d[bc] == "in").sum(), (d[bc] == "out").sum()
        if min(ni, no) < 3000: L.append(f"{spec}, {lo}-{hi}N: near {ni}, away {no} box-passes, too few for the difference"); continue
        B1, k1 = fit2(d[d[bc] == "in"], True, dT, ex); B2, k2 = fit2(d[d[bc] == "out"], True, dT, ex); f1, f2 = contrast_high(k1), contrast_high(k2)
        th = f1(B1.fit()) - f2(B2.fit()); dr = boot_diff(B1, f1, B2, f2, 4000, rng); l, h = np.nanpercentile(dr, [2.5, 97.5])
        L.append(f"{spec}, {lo}-{hi}N: near {ni}, away {no}; near {pct(f1(B1.fit())):+.1f}%, away {pct(f2(B2.fit())):+.1f}%, near minus away {pct(th):+.2f}% ({pct(l):+.2f} to {pct(h):+.2f})")
open(sys.argv[2], "w").write("\n".join(L) + "\n"); print("\n".join(L))
