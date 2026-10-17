"""Pre-registered tests A (ASCAT + ERA5 control, 12 winters) and B (GFS 24 h vs ASCAT, 5 winters). Usage: python3 analyse_region.py boxes_A.csv.gz boxes_B.csv.gz results/"""
import sys, os, numpy as np, pandas as pd
from region_lib import *
rng = np.random.default_rng(20261017); NB = 10000
A0 = pd.read_csv(sys.argv[1]); B0 = pd.read_csv(sys.argv[2]); OUT = sys.argv[3]; os.makedirs(OUT, exist_ok=True); L = []
def P(*a): s = " ".join(str(x) for x in a); print(s, flush=True); L.append(s)
def ci(x): return np.nanpercentile(x, [2.5, 97.5])
def nwflow(d, col): return (d[col] >= 270) & (d[col] <= 360)
def run_family(name, d_all, with_vg, vref, groups, nw_col):
    P(f"\\n## Family {name}: wind change when water is 6 C or more warmer than air (mean of 6-9 and >=9 bins vs 0-3 C)")
    fits = {}
    for g, sel in groups.items():
        d = d_all[sel(d_all)]; B, k = fit_boot(d, with_vg); fits[g] = (B, k, len(d), d.winter.nunique(), d)
    res = {}
    for g in groups:
        B, k, n, nw, d = fits[g]; f = contrast_high(k); th = f(B.fit()); dr = B.boot(f, NB, rng); res[g] = (th, dr, n, nw)
    gi = list(groups)
    th4 = res[gi[1]][0] - res[gi[2]][0]; dr4 = boot_diff(fits[gi[1]][0], contrast_high(fits[gi[1]][1]), fits[gi[2]][0], contrast_high(fits[gi[2]][1]), NB, rng)
    tests = [(f"{name}1 {gi[0]}", res[gi[0]][0], res[gi[0]][1], res[gi[0]][2]), (f"{name}2 {gi[1]}", res[gi[1]][0], res[gi[1]][1], res[gi[1]][2]),
             (f"{name}3 {gi[2]}", res[gi[2]][0], res[gi[2]][1], res[gi[2]][2]), (f"{name}4 {gi[1]} minus {gi[2]}", th4, dr4, None)]
    if len(gi) > 3: tests.append((f"{name}5 {gi[3]}", res[gi[3]][0], res[gi[3]][1], res[gi[3]][2]))
    pv = [pval(t[2]) for t in tests]; qv = bh(pv)
    P("test | n box-passes | winters | estimate % | 95% interval | kt at ref wind {} / {} kt | p | q".format(*vref))
    for t, p, q in zip(tests, pv, qv):
        lo, hi = ci(t[2]); kt = lambda x, v: v * (np.exp(x) - 1)
        nn = t[3] if t[3] is not None else ""
        P(f"{t[0]} | {nn} | {len(np.unique(d_all.winter))} | {pct(t[1]):+.2f} | {pct(lo):+.2f} to {pct(hi):+.2f} | {kt(t[1], vref[0]):+.2f} ({kt(lo, vref[0]):+.2f} to {kt(hi, vref[0]):+.2f}) / {kt(t[1], vref[1]):+.2f} | {p:.4f} | {q:.4f}")
    return tests, pv, qv, fits
# ---- A
A = prep_A(A0); P("# Gulf Stream, region-wide: ASCAT winds against ERA5 pressure-gradient wind (A) and against GFS 24 h forecasts (B)")
P(f"A sample {len(A)} box-passes, {A.winter.nunique()} winters, {A.sat.value_counts().to_dict()}; bands {A.band.value_counts().to_dict()}")
vA = [float(A[(A.dT >= 0) & (A.dT < 3) & (abs(A.vg - v) <= 5)].v_ascat.median()) for v in (30, 45)]
P(f"Reference ASCAT wind (0-3 C bin) at pressure-gradient wind 30 / 45 kt: {vA[0]:.1f} / {vA[1]:.1f} kt")
P("\\n## Observed median ASCAT wind (kt) by water-minus-air temperature at pressure-gradient wind 25-35 kt, by band")
for gname, sel in [("all", lambda d: d.v_ascat > 0), ("in band", lambda d: d.band == "in"), ("outside", lambda d: d.band == "out")]:
    r = []
    for lo, hi in [(-9, 0), (0, 3), (3, 6), (6, 9), (9, 99)]:
        a = A[sel(A) & (A.dT >= lo) & (A.dT < hi) & (abs(A.vg - 30) <= 5)]; r.append(f"{lo}..{hi}: {a.v_ascat.median():.1f} ({len(a)})")
    P(gname, " | ".join(r))
groupsA = {"all boxes": lambda d: d.v_ascat > 0, "in band": lambda d: d.band == "in", "outside band": lambda d: d.band == "out", "in band, flow 270-360": lambda d: (d.band == "in") & nwflow(d, "dir_ascat")}
tA, pA, qA, fA = run_family("A", A, True, vA, groupsA, "dir_ascat")
# ---- B
B = prep_B(B0); P(f"\\nB sample {len(B)} box-passes, {B.winter.nunique()} winters, bands {B.band.value_counts().to_dict()}")
ratio_ref = float(np.exp(B[(B.dT >= 0) & (B.dT < 3)].y.median())); P(f"ASCAT / GFS wind ratio in the 0-3 C bin (median): {ratio_ref:.3f} (a model wind bias, same in every bin)")
vB = [25.0, 35.0]
P("\\n## Observed median ASCAT/GFS ratio by model water-minus-air temperature, by band")
for gname, sel in [("all", lambda d: d.v_ascat > 0), ("in band", lambda d: d.band == "in"), ("outside", lambda d: d.band == "out")]:
    P(gname, " | ".join(f"{lo}..{hi}: {np.exp(B[sel(B) & (B.dT >= lo) & (B.dT < hi)].y.median()):.3f} ({int((sel(B) & (B.dT >= lo) & (B.dT < hi)).sum())})" for lo, hi in [(-9, 0), (0, 3), (3, 6), (6, 9), (9, 99)]))
groupsB = {"all boxes": lambda d: d.v_ascat > 0, "in band": lambda d: d.band == "in", "outside band": lambda d: d.band == "out"}
tB, pB, qB, fB = run_family("B", B, False, vB, groupsB, "dir_gfs")
allq = bh(pA + pB); P("\\nq across all 9 tests:", " ".join(f"{x:.4f}" for x in allq))
# ---- secondary and checks
P("\\n## Secondary (raw p, labelled secondary)")
def sec(label, d, with_vg, nb=3000):
    B_, k = fit_boot(d, with_vg); f = contrast_high(k); th = f(B_.fit()); dr = B_.boot(f, nb, rng)
    P(f"{label}: n={len(d)}  {pct(th):+.2f}% ({pct(ci(dr)[0]):+.2f}..{pct(ci(dr)[1]):+.2f}), p={pval(dr):.4f}")
sec("S2 A1 Dec-Mar only", A[A.month.isin([12, 1, 2, 3])], True)
for s in "ABC": sec(f"S3 A1 Metop-{s} only", A[A.sat == s], True)
sec("S4 A1 with pressure-gradient wind >= 30 kt", A[A.vg >= 30], True)
Ac = A.copy(); U = Ac.v_ascat / 1.944; Rb = 9.81 / Ac.t2m * (Ac.t2m - Ac.sst) * 10 / np.maximum(U, 1) ** 2
z = np.where(Rb < 0, Rb, Rb / (1 - 5 * np.clip(Rb, None, 0.15))); x = (1 - 16 * np.clip(z, -2, 0)) ** 0.25
psi = np.where(z < 0, 2 * np.log((1 + x) / 2) + np.log((1 + x * x) / 2) - 2 * np.arctan(x) + np.pi / 2, -5 * z)
Ac["v_ascat"] = Ac.v_ascat * (1 - psi / 10.8); Ac["y"] = np.log(Ac.v_ascat)
sec("S5 A1 after equivalent-neutral to actual 10 m wind conversion (approximate)", Ac, True)
sec("S6 B-style: A1 with GFS-free design but 5 winters 2021-25 only", A[A.winter >= 2021], True)
sec("B5 B2 in band, flow 270-360 by GFS direction", B[(B.band == "in") & nwflow(B, "dir_gfs")], False)
P("\\n## Leave one winter out (min..max of % effect)")
for nm, (B_, k, n, nw, d) in {**{"A " + g: fA[g] for g in groupsA}, **{"B " + g: fB[g] for g in groupsB}}.items():
    f = contrast_high(k); r = []
    for i in range(len(B_.w)):
        c = np.ones(len(B_.w)); c[i] = 0; r.append(pct(f(B_.fit(c))))
    P(f"{nm}: {min(r):+.2f}% .. {max(r):+.2f}%  ({len(B_.w)} winters)")
P("\\n## Placebo: shuffle dT among boxes of the same winter, month and Vg tercile (20 draws; A) / winter, month (B)")
for nm, d, wv, nmz in [("A1", A, True, "vg"), ("B1", B, False, None)]:
    e = []
    for _ in range(20):
        x = d.copy()
        if nmz: x["vt"] = x.groupby("winter")[nmz].transform(lambda v: pd.qcut(v, 3, labels=False, duplicates="drop")); keys = ["winter", "month", "vt"]
        else: keys = ["winter", "month"]
        x["dT"] = x.groupby(keys).dT.transform(lambda v: rng.permutation(v.values))
        B_, k = fit_boot(x, wv); e.append(pct(contrast_high(k)(B_.fit())))
    P(f"{nm}: placebo mean {np.mean(e):+.2f}%  range {min(e):+.2f}..{max(e):+.2f}")
open(f"{OUT}/region-result.txt", "w").write("\n".join(L) + "\n")
