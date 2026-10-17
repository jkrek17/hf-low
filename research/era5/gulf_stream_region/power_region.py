"""Planted-effect simulation for the region-wide tests (PREREGISTRATION.md 'Power'). Prints detection rates and SEs only.
Residuals after the FULL fit (as in the buoy study); planted log(1+p) added where dT >= 6. 100 outer winter resamples x 200 inner.
Usage: python3 power_region.py boxes_A.csv.gz boxes_B.csv.gz OUT.txt"""
import sys, numpy as np, pandas as pd
from region_lib import *
rng = np.random.default_rng(20261018); NO, NI = 100, 200
A = prep_A(pd.read_csv(sys.argv[1])); B = prep_B(pd.read_csv(sys.argv[2])); L = []
def resample(d, W):
    pick = rng.choice(W, len(W), replace=True); idx = np.concatenate([np.nonzero(d.winter.values == w)[0] for w in pick])
    new = np.concatenate([np.full((d.winter.values == w).sum(), i) for i, w in enumerate(pick)]); s = d.iloc[idx].copy(); s["winter"] = new; return s
def power_for(d, with_vg, targets, plant_sel, pct_):
    """targets: dict name -> (selector for fitting, selector for second group or None). Returns {name: (power, se)}."""
    Xr, k0 = design_raw(d, with_vg); W = np.unique(d.winter.values); rem = np.zeros(len(d))
    for g in ("in", "out", "mid", "undef"):                      # remove each band group's own fitted dT effect (group-specific full fit)
        m = (d.band.values == g)
        if m.sum() > 1000:
            B0, _ = fit_boot(d[m], with_vg); rem[m] = Xr[m, k0:] @ B0.fit()[k0:]
    d = d.copy(); d["y"] = d.y.values - rem + np.log(1 + pct_) * ((d.dT.values >= 6) & plant_sel(d)).astype(float)   # keep residuals, add planted
    hits = {n: 0 for n in targets}; ests = {n: [] for n in targets}
    for _ in range(NO):
        s = resample(d, W)
        for n, (s1, s2) in targets.items():
            B1, k1 = fit_boot(s[s1(s)], with_vg); f1 = contrast_high(k1); e = f1(B1.fit())
            if s2 is None: dr = B1.boot(f1, NI, rng)
            else:
                B2, k2 = fit_boot(s[s2(s)], with_vg); f2 = contrast_high(k2); e -= f2(B2.fit()); dr = boot_diff(B1, f1, B2, f2, NI, rng)
            ests[n].append(e); hits[n] += (pval(dr) < 0.05) and e > 0
    return {n: (hits[n] / NO, float(np.std(ests[n]))) for n in targets}
allb = lambda d: d.v_ascat > 0; inb = lambda d: d.band == "in"; outb = lambda d: d.band == "out"
L.append(f"Planted-effect power; {NO} outer x {NI} inner winter resamples; unadjusted p < 0.05 with the expected sign. Residuals after the full fit.")
for pc in (0.03, 0.06):
    r = power_for(A, True, {"A1 all boxes": (allb, None), "A2 in band": (inb, None), "A4 in minus outside": (inb, outb)}, allb, pc)
    for n, (p, se) in r.items(): L.append(f"A: +{pc*100:.0f}% planted everywhere where dT>=6 -> {n}: power {p:.2f} (SE {se:.4f})")
    r = power_for(A, True, {"A4 in minus outside": (inb, outb)}, inb, pc)
    for n, (p, se) in r.items(): L.append(f"A: +{pc*100:.0f}% planted only in band -> {n}: power {p:.2f} (SE {se:.4f})")
    r = power_for(B, False, {"B1 all boxes": (allb, None), "B2 in band": (inb, None), "B4 in minus outside": (inb, outb)}, allb, pc)
    for n, (p, se) in r.items(): L.append(f"B: +{pc*100:.0f}% planted everywhere where dT>=6 -> {n}: power {p:.2f} (SE {se:.4f})")
    r = power_for(B, False, {"B4 in minus outside": (inb, outb)}, inb, pc)
    for n, (p, se) in r.items(): L.append(f"B: +{pc*100:.0f}% planted only in band -> {n}: power {p:.2f} (SE {se:.4f})")
open(sys.argv[3], "w").write("\n".join(L) + "\n"); print("\n".join(L))
