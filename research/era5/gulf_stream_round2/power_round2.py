"""Planted-effect simulation for round 2 (PREREGISTRATION.md 'Power'). Prints detection rates only; no coefficient is read.
Remove each band group's own fitted dT effect, plant log(1+p) where dT >= 6, keep the residuals; resample winters (outer) then winters of the resample (inner), using the
per-winter sufficient statistics (the winter resample is a count vector). Usage: python3 power_round2.py boxes2_A.csv.gz boxes2_B.csv.gz CACHE_DIR OUT.txt"""
import sys, numpy as np, pandas as pd
from round2_lib import *
from bulk import neutral_from_actual
FA, FB, CACHE, OUTF = sys.argv[1:5]; rng = np.random.default_rng(20261020); NO, NI = 200, 400
A = pd.read_csv(FA).dropna(subset=["v_ascat", "vg", "dT", "tadv", "curv", "d2m", "sst_mur", "msl_hpa"]); A = A[(A.vg >= 15) & (A.v_ascat >= 5)].copy()
B = pd.read_csv(FB).dropna(subset=["v_ascat", "v_gfs", "dT", "dpt", "pres", "sst_mur"]); B = B[(B.v_gfs >= 5) & (B.v_ascat >= 5)].copy()
A["lvg"] = np.log(A.vg); A["dTm"] = (A.sst_mur - A.t2m).clip(-3, 12)
A["tadv_bin"] = tadv_bin(A.tadv.values); q = A.curv.quantile([1 / 3, 2 / 3]).values; A["curv_bin"] = np.select([A.curv <= q[0], A.curv <= q[1]], ["lo", "mid"], "hi")
z = np.load(f"{CACHE}/conv_round2.npz")                                  # written by analyse_round2.py (same sample filters)
assert len(z["b_en"]) == len(B), "run analyse_round2.py first so the conversion cache matches this sample"
B["v_gfs_en"] = z["b_en"]; ex = ("tadv_bin", "curv_bin"); L = []
def planted(d, wv, dTcol, extras, ycol, bcol, plant_sel, p):
    y = ycol(d).values.copy(); X, k0 = design2(d, wv, dTcol, extras); rem = np.zeros(len(d))
    for g in d[bcol].unique():
        m = (d[bcol] == g).values
        if m.sum() > 1000:
            Bt, k = fit2(d[m].assign(y=y[m]), wv, dTcol, extras); b = Bt.fit(); rem[m] = X[m, k0:k0 + 4] @ b[k0:k0 + 4]
    return y - rem + np.log(1 + p) * ((d[dTcol].values >= 6) & plant_sel.values)
def boots(d, wv, dTcol, extras, y, masks):
    out = []
    for m in masks:
        Bt, k = fit2(d[m].assign(y=y[m.values]), wv, dTcol, extras); out.append((Bt, contrast_high(k)))
    return out
def sim(bs):
    """bs: [(Boot, f)] one (single test) or two (difference). Returns detection rate and SE of the estimate."""
    W = len(bs[0][0].w); hits = 0; est = []
    for _ in range(NO):
        pick = rng.integers(0, W, W)
        cnt = [np.bincount(pick, minlength=W).astype(float)] * len(bs)
        th = bs[0][1](bs[0][0].fit(cnt[0])); 
        if len(bs) == 2: th -= bs[1][1](bs[1][0].fit(cnt[1]))
        dr = np.empty(NI)
        for i in range(NI):
            v = bs[0][1](bs[0][0].fit(np.bincount(rng.choice(pick, W), minlength=W).astype(float)))
            if len(bs) == 2: v -= bs[1][1](bs[1][0].fit(np.bincount(rng.choice(pick, W), minlength=W).astype(float)))
            dr[i] = v
        est.append(th); hits += (pval(dr) < 0.05) and th > 0
    return hits / NO, float(np.std(est))
def run(label, d, wv, dTcol, extras, ycol, bcol, pcs):
    inm, outm, allm = (d[bcol] == "in"), (d[bcol] == "out"), d.v_ascat > 0
    for pc in pcs:
        for where, plant in (("everywhere where dT>=6", allm), ("only near the front", inm)):
            y = planted(d, wv, dTcol, extras, ycol, bcol, plant, pc)
            tests = {"1 all boxes": [allm], "4 near minus away": [inm, outm]} if where.startswith("every") else {"4 near minus away": [inm, outm]}
            for n, masks in tests.items():
                p, se = sim(boots(d, wv, dTcol, extras, y, masks)); L.append(f"{label}: +{pc*100:.0f}% planted {where} -> {n}: power {p:.2f} (SE {se:.4f})"); print(L[-1], flush=True)
L.append(f"Planted-effect power; {NO} outer x {NI} inner winter resamples (count-vector resamples of per-winter sufficient statistics); unadjusted p < 0.05 with the expected sign.")
run("A1c/A4c (controls)", A, True, "dT", ex, lambda d: np.log(d.v_ascat), "band", (0.03, 0.06))
run("A1m/A4m (MUR front and temperature)", A, True, "dTm", (), lambda d: np.log(d.v_ascat), "band_mur", (0.03,))
run("B1e/B4e (like-with-like winds)", B, False, "dT", (), lambda d: np.log(d.v_ascat / d.v_gfs_en), "band", (0.03, 0.06))
open(OUTF, "w").write("\n".join(L) + "\n")
