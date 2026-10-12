"""Part A: placebo pattern indices. usage: python3 -I run_a.py OUTDIR [NPERM NBOOT NSHUF]. Plan: PREREGISTRATION.md (23f1b5d)."""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from multiprocessing import Pool
import numpy as np, pandas as pd
import placelib as P

OUT = sys.argv[1]
NPERM = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
NBOOT = int(sys.argv[3]) if len(sys.argv) > 3 else 2000
NSHUF = int(sys.argv[4]) if len(sys.argv) > 4 else 100
os.makedirs(OUT, exist_ok=True)
sys.argv = ["part1", P.REPO, OUT]                      # part1 reads argv at import
sys.path.insert(0, os.path.join(P.REPO, "research", "era5", "share_mechanism"))
import part1 as S1  # noqa: E402
C = P.C
SEED = 20261011
T = P.load_table()
Tbase = S1.load(-3.6)
M = C.month_dummies()
ALLIDX = list(range(len(C.SEASONS)))


def outcome_counts(b):
    t = Tbase[Tbase.basin == b]
    win = C.assign_week(t)
    g = t.gust800_kt
    masks = dict(N=pd.Series(True, index=t.index), D=t.deep == 1, S=(g >= 55) & (g < 71.7))
    out = {}
    for k, m in masks.items():
        out[k] = (C.weekly_counts(win, m.reindex(win.index)), C.prev_week_counts(t, m, b))
    return out


def shuffled(args):
    b, seed, y, prev, Y, LP = args
    rng = np.random.default_rng(seed)
    perm = rng.permutation(len(C.SEASONS))
    yy = np.concatenate([y[p * P.NW:(p + 1) * P.NW] for p in perm])
    pp = np.concatenate([prev[p * P.NW:(p + 1) * P.NW] for p in perm])
    x, _ = P.index_from_outcome(T, b, yy, pp)
    rows = np.arange(len(x))
    f = S1.factors(S1.slopes(M, LP, Y, x, rows))
    return dict(basin=b, seed=seed, A1=f["A1"], A2=f["A2"], A3=f["A3"], A4=f["A4"])


def main():
    t0 = time.time()
    jobs, idxinfo, shuf = [], [], []
    for b in C.BASINS:
        Y, LP, _ = S1.weekly(Tbase, b, "hf")
        oc = outcome_counts(b)
        committed = C.index(b)[0]
        idxs = {"P0": (T[f"y_{b}"].values, T[f"prev_{b}"].values),
                "P1": oc["N"], "P2": oc["D"], "P3": oc["S"]}
        for name, (y, prev) in idxs.items():
            x, lams = P.index_from_outcome(T, b, y, prev)
            idxinfo.append(dict(basin=b, index=name, r_committed=float(np.corrcoef(x, committed)[0, 1]), lam_modal=float(pd.Series(lams).mode()[0]), train_events=float(np.sum(y))))
            jobs.append((f"{name}", b, Y, LP, x, ALLIDX, NPERM, NBOOT, SEED))
            print(b, name, f"{time.time()-t0:.0f}s", flush=True)
        shuf += [(b, SEED + 1000 * (b == "pac") + i, T[f"y_{b}"].values, T[f"prev_{b}"].values, Y, LP) for i in range(NSHUF)]
    with Pool(4) as pool:
        res = pool.map(S1.analyse, jobs, chunksize=1)
        sh = pool.map(shuffled, shuf, chunksize=2)
    rows = []
    for r in res:
        r["verdict"] = S1.verdict(r)
        for k, fk in (("A1", "f_deep"), ("A2", "f_conv"), ("A3", "f_res")):
            rows.append(dict(index=r["label"], basin=r["basin"], test=k, rr=r[k]["rr"], rr_lo=r[k]["rr_lo"], rr_hi=r[k]["rr_hi"], p=r[k]["p"],
                             f=r[fk]["est"], f_lo=r[fk]["lo"], f_hi=r[fk]["hi"], A4_rr=r["A4"]["rr"], A4_lo=r["A4"]["rr_lo"], A4_hi=r["A4"]["rr_hi"],
                             frac_boot_ok=r["frac_boot_ok"], verdict=r["verdict"]))
    D = pd.DataFrame(rows)
    D["q_placebo"] = np.nan
    m = D["index"].isin(["P1", "P2", "P3"])
    D.loc[m, "q_placebo"] = C.bh(D.p[m].values)               # BH over the 18 placebo tests
    D.to_csv(os.path.join(OUT, "partA_tests.csv"), index=False)
    pd.DataFrame(idxinfo).to_csv(os.path.join(OUT, "partA_indices.csv"), index=False)
    Sh = pd.DataFrame(sh)
    Sh["S"] = Sh[["A1", "A2", "A3"]].abs().pipe(lambda d: d.A2 / d.sum(axis=1))
    Sh.to_csv(os.path.join(OUT, "partA_shuffle.csv"), index=False)
    json.dump(res, open(os.path.join(OUT, "partA_raw.json"), "w"), indent=1, default=float)
    print("done", f"{time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
