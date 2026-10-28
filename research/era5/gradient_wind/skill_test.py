"""Does the ERA5 pressure-gradient wind add to the GFS 24 h wind in predicting ASCAT? Plan: PREREGISTRATION_SKILL.md.
Usage: python3 skill_test.py OUTFILE   (reads /mnt/project-files/gulf-stream/tables; no download)"""
import sys
import numpy as np, pandas as pd

T = "/mnt/project-files/gulf-stream/tables"
B = 2000
rng = np.random.default_rng(20261028)
KEY = ["day", "seg", "lat", "lon", "sat"]


def X(d, name):
    v, g = d.v_gfs.values, d.vg.values
    cols = {"M0": [v], "M1": [v, g], "M1q": [v, g, v * v, g * g, v * g], "M2": [v, g, d.dT.values], "M3": [g]}[name]
    return np.column_stack([np.ones(len(d))] + cols)


def lowo(d, name):
    pred = np.empty(len(d))
    for w in d.winter.unique():
        tr, te = (d.winter != w).values, (d.winter == w).values
        beta = np.linalg.lstsq(X(d[tr], name), d.v_ascat.values[tr], rcond=None)[0]
        pred[te] = X(d[te], name) @ beta
    return pred


def skill(d, p0, p1, m):
    e0 = (d.v_ascat.values - p0) ** 2; e1 = (d.v_ascat.values - p1) ** 2
    days = d.day.values[m]; u, inv = np.unique(days, return_inverse=True)
    s0 = np.bincount(inv, e0[m], len(u)); s1 = np.bincount(inv, e1[m], len(u)); n = np.bincount(inv, minlength=len(u)).astype(float)
    ix = rng.integers(0, len(u), (B, len(u)))
    bs = 1 - s1[ix].sum(1) / s0[ix].sum(1)
    per = {int(w): 1 - e1[m & (d.winter.values == w)].mean() / e0[m & (d.winter.values == w)].mean() for w in np.unique(d.winter.values[m])}
    return 1 - e1[m].mean() / e0[m].mean(), np.percentile(bs, [2.5, 97.5]), np.sqrt(e0[m].mean()), np.sqrt(e1[m].mean()), int(m.sum()), per, bs


def main(out):
    a = pd.read_csv(f"{T}/boxes_A.csv.gz", usecols=KEY + ["vg", "winter"]); b = pd.read_csv(f"{T}/boxes_B.csv.gz")
    d = b.merge(a[KEY + ["vg"]], on=KEY, how="inner").dropna(subset=["v_ascat", "v_gfs", "vg", "dT"])
    d = d[d.winter <= 2024].reset_index(drop=True)
    L = [f"box-passes {len(d)}, winters {sorted(d.winter.unique())}, days {d.day.nunique()}; GFS wind >=25 kt: {int((d.v_gfs>=25).sum())}, >=35 kt: {int((d.v_gfs>=35).sum())}"]
    P = {n: lowo(d, n) for n in ["M0", "M1", "M1q", "M2", "M3"]}
    allm = np.ones(len(d), bool); strong = (d.v_gfs >= 25).values; vstrong = (d.v_gfs >= 35).values
    res = {}
    for lab, m in [("T1 all", allm), ("T2 GFS wind >= 25 kt", strong), ("GFS wind >= 35 kt (secondary)", vstrong)]:
        for mod in ["M1", "M1q", "M2", "M3"]:
            s, ci, r0, r1, n, per, bs = skill(d, P["M0"], P[mod], m)
            res[(lab, mod)] = (s, ci, bs)
            L.append(f"{lab:32s} {mod:3s} vs M0: skill {100*s:+6.2f}% ({100*ci[0]:+.2f} to {100*ci[1]:+.2f}) RMSE {r0:.2f} -> {r1:.2f} kt, n {n}, by winter " +
                     ", ".join(f"{w}:{100*v:+.1f}%" for w, v in per.items()))
    p = []
    for lab in ["T1 all", "T2 GFS wind >= 25 kt"]:
        bs = res[(lab, "M1")][2]; p.append(float(2 * min((bs <= 0).mean(), (bs >= 0).mean())))
    p = np.maximum(p, 1 / B); q = [min(1, p[0] * 2), min(1, p[1])] if p[0] <= p[1] else [min(1, p[0]), min(1, p[1] * 2)]
    L.append(f"BH over T1, T2: p {p[0]:.4f}, {p[1]:.4f}; q {q[0]:.4f}, {q[1]:.4f}")
    L.append("Baseline M0 coefficients (fit on all four winters): " + str(np.round(np.linalg.lstsq(X(d, 'M0'), d.v_ascat.values, rcond=None)[0], 3)))
    L.append("M1 coefficients (all four winters): " + str(np.round(np.linalg.lstsq(X(d, 'M1'), d.v_ascat.values, rcond=None)[0], 3)))
    open(out, "w").write("\n".join(L) + "\n"); print("\n".join(L))


if __name__ == "__main__":
    main(sys.argv[1])
