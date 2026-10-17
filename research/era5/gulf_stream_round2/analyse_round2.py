"""Round 2 tests (PREREGISTRATION.md). Usage: python3 analyse_round2.py boxes2_A.csv.gz boxes2_B.csv.gz CACHE_DIR OUT_DIR"""
import sys, os, numpy as np, pandas as pd
from round2_lib import *
from bulk import neutral_from_actual, actual_from_neutral
FA, FB, CACHE, OUT = sys.argv[1:5]; os.makedirs(OUT, exist_ok=True); os.makedirs(CACHE, exist_ok=True)
rng = np.random.default_rng(20261019); NB = 10000; L = []
def P(*a): s = " ".join(str(x) for x in a); print(s, flush=True); L.append(s)
def ci(x): return np.nanpercentile(x, [2.5, 97.5])
# ---------- samples
A0 = pd.read_csv(FA); B0 = pd.read_csv(FB)
A = A0.dropna(subset=["v_ascat", "vg", "dT", "tadv", "curv", "d2m", "sst_mur", "msl_hpa"]); A = A[(A.vg >= 15) & (A.v_ascat >= 5)].copy()
B = B0.dropna(subset=["v_ascat", "v_gfs", "dT", "dpt", "pres", "sst_mur"]); B = B[(B.v_gfs >= 5) & (B.v_ascat >= 5)].copy()
P(f"# Gulf Stream round 2. A sample {len(A)} box-passes ({A.winter.nunique()} winters; round 1 had 781652), B sample {len(B)} ({B.winter.nunique()} winters; round 1 had 474396)")
A["lvg"] = np.log(A.vg); A["dTm"] = (A.sst_mur - A.t2m).clip(-3, 12)
A["tadv_bin"] = tadv_bin(A.tadv.values); q = A.curv.quantile([1 / 3, 2 / 3]).values; A["curv_bin"] = np.select([A.curv <= q[0], A.curv <= q[1]], ["lo", "mid"], "hi")
# ---------- bulk conversions (cached)
cf = f"{CACHE}/conv_round2.npz"
if os.path.exists(cf) and len(np.load(cf)["a_act"]) == len(A) and len(np.load(cf)["b_en"]) == len(B):
    z = np.load(cf); A["v_act"], A["v_act_m"], B["v_gfs_en"], B["v_asc_act"] = z["a_act"], z["a_act_m"], z["b_en"], z["b_act"]
else:
    A["v_act"], r1 = actual_from_neutral(A.v_ascat.values, A.t2m.values, A.d2m.values, A.sst.values, A.msl_hpa.values, A.lat.values, 1)
    A["v_act_m"], r2 = actual_from_neutral(A.v_ascat.values, A.t2m.values, A.d2m.values, A.sst_mur.values, A.msl_hpa.values, A.lat.values, 1)
    B["v_gfs_en"] = neutral_from_actual(B.v_gfs.values, B.t2m.values, B.dpt.values, B.skt.values, B.pres.values, B.lat.values, 0)
    B["v_asc_act"], r3 = actual_from_neutral(B.v_ascat.values, B.t2m.values, B.dpt.values, B.skt.values, B.pres.values, B.lat.values, 0)
    np.savez_compressed(cf, a_act=A.v_act.values, a_act_m=A.v_act_m.values, b_en=B.v_gfs_en.values, b_act=B.v_asc_act.values); P(f"inversion residuals (kt): {r1:.3f} {r2:.3f} {r3:.3f}")
P(f"Mean conversion (neutral/actual) A, dT>=6: {np.nanmean((A.v_ascat / A.v_act)[A.dT >= 6]):.4f}; dT 0-3: {np.nanmean((A.v_ascat / A.v_act)[(A.dT >= 0) & (A.dT < 3)]):.4f}; B GFS neutral/actual dT>=6: {np.nanmean((B.v_gfs_en / B.v_gfs)[B.dT >= 6]):.4f}")
if os.environ.get("CONV_ONLY"): sys.exit(0)          # builds the conversion cache for power_round2.py before any coefficient is read
# ---------- generic test runner. A test is (name, spec, spec2); the estimate is spec minus spec2 (spec2 None: a single contrast)
def fitspec(s):
    d = s["d"].assign(y=s["y"](s["d"])); Bt, k = fit2(d, s["wv"], s["dT"], s["extras"]); return Bt, contrast_high(k), len(d)
def sp(d, wv, dT="dT", extras=(), y=None): return dict(d=d, wv=wv, dT=dT, extras=extras, y=y)
def run(family, tests, vref):
    res = []
    for name, s1, s2 in tests:
        B1, f1, n1 = fitspec(s1); th = f1(B1.fit())
        if s2 is None: dr = B1.boot(f1, NB, rng); n = n1
        else: B2, f2, _ = fitspec(s2); th -= f2(B2.fit()); dr = boot_diff(B1, f1, B2, f2, NB, rng); n = ""
        res.append((name, th, dr, n, len(B1.w)))
    pv = [pval(r[2]) for r in res]; qv = bh(pv)
    P(f"\n## {family}\ntest | n | winters | estimate % | 95% interval | kt at reference wind {vref} kt | p | q")
    for (name, th, dr, n, nw), p, q in zip(res, pv, qv):
        lo, hi = ci(dr); kt = lambda x: vref * (np.exp(x) - 1)
        P(f"{name} | {n} | {nw} | {pct(th):+.2f} | {pct(lo):+.2f} to {pct(hi):+.2f} | {kt(th):+.2f} ({kt(lo):+.2f} to {kt(hi):+.2f}) | {p:.4f} | {q:.4f}")
    return res, pv
ya = lambda d: np.log(d.v_ascat); yb = lambda d: np.log(d.v_ascat / d.v_gfs)
def four(prefix, d, wv, y, dT="dT", extras=(), bcol="band", which="1234"):
    sub = lambda m: sp(d[m], wv, dT, extras, y); allm = d.v_ascat > 0; inm = d[bcol] == "in"; outm = d[bcol] == "out"
    t = {"1": (f"{prefix}1 all boxes", sub(allm), None), "2": (f"{prefix}2 near front", sub(inm), None), "3": (f"{prefix}3 away from front", sub(outm), None),
         "4": (f"{prefix}4 near minus away", sub(inm), sub(outm))}
    return [t[c] for c in which]
REF_A, REF_B = 21.6, 25.0
# ---------- baselines on the restricted samples (round 1 specification)
run("Baseline A on this sample (round 1 specification: ERA5 SST and front, no controls)", four("A", A, True, ya), REF_A)
run("Baseline B on this sample (round 1 specification: actual 10 m winds)", four("B", B, False, yb), REF_B)
# ---------- F1
t1 = four("B", B, False, lambda d: np.log(d.v_ascat / d.v_gfs_en)); t1 = [(n.replace("B", "B", 1) + " (like with like)", a, b) for n, a, b in t1]
t1 = [(f"B{n[1]}e" + n[2:], a, b) for n, a, b in t1]
t1 += [("A1a all boxes (ASCAT converted to actual 10 m wind)", sp(A, True, "dT", (), lambda d: np.log(d.v_act)), None)]
F1, p1 = run("Family F1: like-with-like wind (COARE 3.6)", t1, REF_B)
run("F1 secondary: B on actual winds (ASCAT converted to actual 10 m, against GFS 10 m)", four("Bact", B, False, lambda d: np.log(d.v_asc_act / d.v_gfs)), REF_B)
# ---------- F2
ex = ("tadv_bin", "curv_bin"); caa = A.tadv <= -1
t2 = four("A", A, True, ya, "dT", ex); t2 = [(f"A{n[1]}c" + n[2:], a, b) for n, a, b in t2]
t2 += [("A1-CAA cold advection (tadv <= -1 K/day)", sp(A[caa], True, "dT", ex, ya), None), ("A1-nonCAA (tadv > -1 K/day)", sp(A[~caa], True, "dT", ex, ya), None)]
F2, p2 = run("Family F2: air mass (thermal advection and curvature controls)", t2, REF_A)
P(f"share of dT>=6 box-passes with cold advection (tadv <= -1): {(caa & (A.dT >= 6)).sum() / (A.dT >= 6).sum():.3f}; near front {(caa & (A.dT >= 6) & (A.band == 'in')).sum() / ((A.dT >= 6) & (A.band == 'in')).sum():.3f}; away {(caa & (A.dT >= 6) & (A.band == 'out')).sum() / ((A.dT >= 6) & (A.band == 'out')).sum():.3f}")
# ---------- F3
t3 = four("A", A, True, ya, "dTm", (), "band_mur"); t3 = [(f"A{n[1]}m" + n[2:], a, b) for n, a, b in t3]
tb = four("B", B, False, yb, "dT", (), "band_mur", "24"); t3 += [(f"B{n[1]}m" + n[2:], a, b) for n, a, b in tb]
F3, p3 = run("Family F3: satellite SST (MUR) for the front and the temperature difference", t3, REF_A)
P(f"MUR band counts A: {A.band_mur.value_counts().to_dict()}; B: {B.band_mur.value_counts().to_dict()}; share of 'away' A box-passes south of MUR front: {(A[A.band_mur == 'out'].lat < A[A.band_mur == 'out'].front_mur).mean():.3f}")
# ---------- all fixes (secondary)
ta = four("A", A, True, lambda d: np.log(d.v_act_m), "dTm", ex, "band_mur"); ta = [(f"A{n[1]}f" + n[2:], a, b) for n, a, b in ta]
run("Secondary: all three fixes together, A (MUR dT and front, controls, true 10 m wind)", ta, REF_A)
tbf = four("B", B, False, lambda d: np.log(d.v_ascat / d.v_gfs_en), "dT", (), "band_mur", "124"); tbf = [(f"B{n[1]}f" + n[2:], a, b) for n, a, b in tbf]
run("Secondary: all fixes together, B (like-with-like winds, MUR front)", tbf, REF_B)
# ---------- currents (winters 2017-18 on)
C = A.dropna(subset=["cur_u", "cur_v"]).copy(); C = C[C.winter >= 2017]
th = np.radians(C.dir_ascat.values); ur, vr = -C.v_ascat.values * np.sin(th), -C.v_ascat.values * np.cos(th)            # wind vector toward, relative to the water
C["v_abs"] = np.hypot(ur + C.cur_u.values * 1.94384449, vr + C.cur_v.values * 1.94384449)
P(f"\ncurrent sample: {len(C)} box-passes, {C.winter.nunique()} winters; mean |current| {np.hypot(C.cur_u, C.cur_v).mean():.2f} m/s, near front {np.hypot(C[C.band=='in'].cur_u, C[C.band=='in'].cur_v).mean():.2f}, away {np.hypot(C[C.band=='out'].cur_u, C[C.band=='out'].cur_v).mean():.2f}")
run("Secondary: current sample, wind relative to the water (as ASCAT reports)", four("Crel", C, True, lambda d: np.log(d.v_ascat)), REF_A)
run("Secondary: current sample, absolute wind (relative wind plus geostrophic current)", four("Cabs", C, True, lambda d: np.log(d.v_abs)), REF_A)
# ---------- multiplicity across the 17
allp = p1 + p2 + p3; P("\nq across all 17 primary tests:", " ".join(f"{x:.4f}" for x in bh(allp)))
# ---------- placebo
PLC = {}
P("\n## Placebo: shuffle dT among boxes of the same winter, month and Vg tercile (A) / winter, month (B), 20 draws")
def shuffled(d, col, wv):
    x = d.copy()
    if wv: x["vt"] = x.groupby("winter")["vg"].transform(lambda v: pd.qcut(v, 3, labels=False, duplicates="drop")); keys = ["winter", "month", "vt"]
    else: keys = ["winter", "month"]
    out = x[col].values.copy()
    for _, ix in x.groupby(keys).indices.items(): out[ix] = rng.permutation(out[ix])
    x[col] = out; return x
for nm, d, wv, col, extras, y, bcol, tag in [("A1c", A, True, "dT", ex, ya, "band", "A"), ("A1m", A, True, "dTm", (), ya, "band_mur", "Am"), ("B1e", B, False, "dT", (), lambda d: np.log(d.v_ascat / d.v_gfs_en), "band", "B")]:
    e1, e4 = [], []
    for _ in range(20):
        x = shuffled(d, col, wv); Bt, f, _n = fitspec(sp(x, wv, col, extras, y)); e1.append(pct(f(Bt.fit())))
        if tag != "B" or True:
            Bi, fi, _n = fitspec(sp(x[x[bcol] == "in"], wv, col, extras, y)); Bo, fo, _n = fitspec(sp(x[x[bcol] == "out"], wv, col, extras, y)); e4.append(pct(fi(Bi.fit()) - fo(Bo.fit())))
    PLC[nm] = np.mean(e1)
    P(f"{nm}: placebo mean {np.mean(e1):+.2f}% (range {min(e1):+.2f}..{max(e1):+.2f}); near-minus-away placebo mean {np.mean(e4):+.2f}% (range {min(e4):+.2f}..{max(e4):+.2f})")
for nm, th, ref in (("A1c", F2[0][1], REF_A), ("A1m", F3[0][1], REF_A), ("B1e", F1[0][1], REF_B)):
    o = np.log(1 + PLC[nm] / 100); P(f"{nm}: estimate {pct(th):+.2f}%, shuffle offset {PLC[nm]:+.2f}%, {'corrected (real minus offset in log wind) ' + format(pct(th - o), '+.2f') + '%, ' + format(ref * (np.exp(th - o) - 1), '+.2f') + ' kt' if abs(PLC[nm]) > 1 else 'offset under the 1-point rule, no correction'}")
open(f"{OUT}/round2-result.txt", "w").write("\n".join(L) + "\n")
