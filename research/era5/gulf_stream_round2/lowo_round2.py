"""Registered B rule: with 5 winters a B result is 'yes' only if its 95% interval excludes 0 in 5 of 5 refits that leave one winter out. Usage: python3 lowo_round2.py boxes2_B.csv.gz CACHE_DIR OUT.txt"""
import sys, numpy as np, pandas as pd
from round2_lib import *
rng = np.random.default_rng(41)
B = pd.read_csv(sys.argv[1]).dropna(subset=["v_ascat", "v_gfs", "dT", "dpt", "pres", "sst_mur"]); B = B[(B.v_gfs >= 5) & (B.v_ascat >= 5)].copy()
z = np.load(f"{sys.argv[2]}/conv_round2.npz"); assert len(z["b_en"]) == len(B); B["v_gfs_en"] = z["b_en"]
defs = {"B1e all boxes (like with like)": (lambda d: d.v_ascat > 0, None, "band", "en"), "B2e near front (like with like)": (lambda d: d.band == "in", None, "band", "en"),
        "B2m near front (MUR band, actual winds)": (lambda d: d.band_mur == "in", None, "band_mur", "act"), "B4m near minus away (MUR band, actual winds)": (lambda d: d.band_mur == "in", lambda d: d.band_mur == "out", "band_mur", "act"),
        "B4e near minus away (like with like)": (lambda d: d.band == "in", lambda d: d.band == "out", "band", "en")}
L = []
for nm, (s1, s2, bc, kind) in defs.items():
    res = []
    for w in sorted(B.winter.unique()):
        d = B[B.winter != w]; d = d.assign(y=np.log(d.v_ascat / (d.v_gfs_en if kind == "en" else d.v_gfs)))
        B1, k1 = fit2(d[s1(d)], False); f1 = contrast_high(k1); th = f1(B1.fit())
        if s2 is None: dr = B1.boot(f1, 4000, rng)
        else: B2, k2 = fit2(d[s2(d)], False); f2 = contrast_high(k2); th -= f2(B2.fit()); dr = boot_diff(B1, f1, B2, f2, 4000, rng)
        lo, hi = np.nanpercentile(dr, [2.5, 97.5]); res.append((w, pct(th), pct(lo), pct(hi)))
    L.append(f"{nm}: interval excludes 0 (above) in {sum(r[2] > 0 for r in res)} of {len(res)} refits, (below) in {sum(r[3] < 0 for r in res)}")
    L += [f"   without winter {w}: {t:+.2f}% ({lo:+.2f} to {hi:+.2f})" for w, t, lo, hi in res]
open(sys.argv[3], "w").write("\n".join(L) + "\n"); print("\n".join(L))
