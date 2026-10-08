"""RA-5: Pacific lead on the Atlantic. Runs the pre-registered tests (PREREGISTRATION.md) and writes results/.
usage: run.py [NPERM] [NBOOT]"""
import json
import zlib
import os
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd

import leadlib as L

NPERM = int(sys.argv[1]) if len(sys.argv) > 1 else 10000
NBOOT = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
OUT = os.path.join(L.HERE, "results")
os.makedirs(OUT, exist_ok=True)
DS = {}


def S(name, kind, k, pred, extras=False, O="atl", Q="pac"):
    return L.Spec(name, kind, k, pred, extras, O, Q)


PRIMARY = [S("P1", "pois", 1, "idx"), S("P2", "pois", 2, "idx"), S("P3", "ols", 1, "idx"), S("P4", "ols", 2, "idx"),
           S("P5", "pois", 1, "cnt"), S("P6", "pois", 2, "cnt")]
S1 = [S("S1-P1", "pois", 1, "idx", True), S("S1-P2", "pois", 2, "idx", True), S("S1-P3", "ols", 1, "idx", True),
      S("S1-P4", "ols", 2, "idx", True)]
S2 = [S("S2-rev-P1", "pois", 1, "idx", False, "pac", "atl"), S("S2-rev-P3", "ols", 1, "idx", False, "pac", "atl")]
S3 = [S("S3-P1", "pois", 1, "idx"), S("S3-P2", "pois", 2, "idx"), S("S3-P3", "ols", 1, "idx"), S("S3-P4", "ols", 2, "idx")]
S4 = [S("S4-P1", "pois", 1, "idx"), S("S4-P2", "pois", 2, "idx")]
JOBS = [(s, "a") for s in PRIMARY + S1 + S2] + [(s, "r") for s in S3] + [(s, "p") for s in S4]


def work(job):
    sp, which = job
    ds = DS[which]
    rng = np.random.default_rng(L.SEED + zlib.crc32(sp.name.encode()) % 1000)
    perms = [np.random.default_rng(L.SEED + 17 + i).permutation(ds.ns) for i in range(NPERM)]
    obs, r, p_pos, p_neg, nul = L.perm_test(ds, sp, perms)
    bt = L.boot(ds, sp, rng, NBOOT)
    lo, hi = np.percentile(bt[:, 0], [2.5, 97.5])
    rlo, rhi = (np.percentile(bt[:, 1], [2.5, 97.5]) if sp.kind == "ols" else (np.nan, np.nan))
    ex = np.exp if sp.kind == "pois" else (lambda v: v)
    return dict(test=sp.name, kind=sp.kind, lag=sp.k, pred=sp.pred, data=ds.label, est=float(ex(obs)), lo=float(ex(lo)),
                hi=float(ex(hi)), partial_r=r, r_lo=rlo, r_hi=rhi, p_pos=p_pos, p_neg=p_neg,
                n_obs=int((L.NW - 2) * ds.ns), n_seasons=ds.ns, null_sd=float(np.std(nul)))


def main():
    DS["a"] = L.load_2004("archive"); DS["p"] = L.load_2004("proxy"); DS["r"] = L.load_1979_2000()
    with Pool(4) as pool:
        rows = pool.map(work, JOBS)
    T = pd.DataFrame(rows)
    fam1 = T.test.isin([s.name for s in PRIMARY])
    T.loc[fam1, "q_family1"] = L.BH(T.p_pos[fam1].values)
    T["q_family2"] = L.BH(T.p_pos.values)
    T.to_csv(os.path.join(OUT, "tests.csv"), index=False)
    print(T.drop(columns=["data", "pred", "null_sd"]).round(4).to_string())


if __name__ == "__main__":
    main()
