"""Collect every pre-registered test (PREREGISTRATION.md), apply Benjamini-Hochberg over the whole family and within H1/H2/H3.
usage: assemble.py RESULTS_DIR"""
import sys, os
import numpy as np, pandas as pd

R = sys.argv[1]


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = np.empty(n)
    c = np.minimum.accumulate((p[o] * n / (np.arange(n) + 1))[::-1])[::-1]
    q[o] = np.minimum(c, 1.0)
    return q


d = pd.read_csv(os.path.join(R, "diagnose_tests.csv"))
d = d.rename(columns={"what": "label"})
d["kind"] = "H1/H2"
est = pd.read_csv(os.path.join(R, "tele_estimates.csv"))
con = pd.read_csv(os.path.join(R, "tele_contrasts.csv"))
pri = lambda x: x[(x.window == "octapr") & (x.seasons == "2004-2025") & (~x.era_term)]
est, con = pri(est), pri(con)
rows = []
for r in est[est.outcome.isin(["G", "D", "Dp", "arch", "share_G", "share_D", "share_Dp"])].itertuples():
    rows.append(dict(family="H3", basin=r.basin, test=f"T3 {r.idx} {r.outcome}", label=f"{r.basin} {r.idx}: RR per SD, {r.outcome}",
                     est=r.RR, lo=r.lo, hi=r.hi, p=r.p_perm, kind="RR"))
for r in con.itertuples():
    rows.append(dict(family="H3", basin=r.basin, test=f"T3 {r.idx} {r.contrast}", label=f"{r.basin} {r.idx}: ratio of RR, {r.contrast}",
                     est=r.ratio_of_RR, lo=r.lo, hi=r.hi, p=r.p_boot, kind="ratio"))
A = pd.concat([d, pd.DataFrame(rows)], ignore_index=True)
A["q_all"] = bh(A.p.values)
for f, g in A.groupby("family"):
    A.loc[g.index, "q_family"] = bh(g.p.values)
A.to_csv(os.path.join(R, "fdr_all_tests.csv"), index=False, float_format="%.5g")
print(len(A), "tests;", "q_all<0.05:", int((A.q_all < 0.05).sum()), " q_all<0.10:", int((A.q_all < 0.10).sum()))
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 70)
print(A[A.family == "H3"][["basin", "label", "est", "lo", "hi", "p", "q_all", "q_family"]].round(4).to_string())
