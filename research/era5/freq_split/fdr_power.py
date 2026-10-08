"""Benjamini-Hochberg FDR across every permutation test in results.csv, and the
effect each test could detect. ERA5 PROXY, pipeline A.

Added after the first debug run (see README "Departures from the plan"): the
pre-registration named Bonferroni x2 over the two primary decompositions; this
adds FDR across all tests run, as the project asked. Neither replaces the other.

Detectable effect: the RR per SD at 80% power, two-sided alpha 0.05, taken as
exp(2.8 x clustered SE of the log coefficient).

usage: fdr_power.py RESULTS_DIR
"""
import os
import sys

import numpy as np
import pandas as pd


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1)


def main():
    d = sys.argv[1]
    R = pd.read_csv(os.path.join(d, "results.csv"))
    rows = []
    for _, r in R.iterrows():
        label = f"{r.analysis} | {r.basin} {r.idx}{' (joint)' if isinstance(r.joint, str) else ''} {r.window} {r.seasons} {r.timing}"
        for factor, b, se, p in (("all cyclones", r.b_all, r.se_all, r.pperm_all),
                                 ("share", r.b_share, r.se_share, r.pperm_share),
                                 (f"{r.outcome} lows", r[f"b_{r.outcome}"], r[f"se_{r.outcome}"], r[f"pperm_{r.outcome}"])):
            rows.append(dict(test=label, factor=factor, rr=np.exp(b), z=b / se, p_perm=p,
                             mde_rr=np.exp(2.8 * se), primary=r.analysis == "primary"))
    F = pd.DataFrame(rows)
    # one test per distinct estimate: S3 repeats the primary all/HF/share and S4 2004-2025
    # repeats the primary all-cyclone factor, so drop exact duplicates of the estimate
    F["key"] = F.rr.round(6).astype(str) + F.factor.str.split().str[0]
    U = F.drop_duplicates("key").copy()
    U["q_bh"] = bh(U.p_perm)
    F = F.merge(U[["key", "q_bh"]], on="key")
    P = F[F.primary].copy()
    P["p_bonf2"] = np.minimum(P.p_perm * 2, 1)
    L = ["ERA5 PROXY, pipeline A. BH FDR across all distinct permutation tests "
         f"({len(U)} tests; duplicated estimates counted once). mde_rr: RR per SD detectable",
         "at 80% power (two-sided 0.05), exp(2.8 x clustered SE). Primary rows also get Bonferroni x2.", ""]
    L.append(F.drop(columns=["key", "primary"]).to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    L += ["", "Primary, Bonferroni x2:", P[["test", "factor", "rr", "p_perm", "p_bonf2", "q_bh"]].to_string(
        index=False, float_format=lambda x: f"{x:.4f}")]
    open(os.path.join(d, "fdr_power.txt"), "w").write("\n".join(L) + "\n")
    F.drop(columns=["key"]).to_csv(os.path.join(d, "fdr_power.csv"), index=False, float_format="%.5g")


if __name__ == "__main__":
    main()
