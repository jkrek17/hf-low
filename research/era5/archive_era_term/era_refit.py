"""RA-29 T-A: PR 79 / RA-22 archive design (PR 41 weekly archive counts) with an era column. Plan: PREREGISTRATION.md.

Reuses sustained_hf/contrast.py unchanged except that the Design gets one extra column E (era step). ERA = none reproduces the
original (control); step = QuikSCAT end 2009-11-23 (primary); asc = ASCAT-B 2012-09-01 (S1); trend = season-linear (S2).
usage: era_refit.py MODE OUT_CSV NPERM NBOOT      MODE in none | step | asc | trend
"""
import datetime as dt, os, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ERA = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(ERA, "sustained_hf"))
import contrast as K      # noqa: E402
C, CH = K.C, K.CH

STEP = dt.date(2009, 11, 23)
ASC = dt.date(2012, 9, 1)


def era_column(mode):
    if mode == "none":
        return None
    starts = [C.week_start(s, w) for s in C.SEASONS for w in range(C.NWEEK)]
    if mode == "step":
        e = [1.0 if d >= STEP else 0.0 for d in starts]
    elif mode == "asc":
        e = [1.0 if d >= ASC else 0.0 for d in starts]
    elif mode == "trend":
        e = [(d - dt.date(2004, 10, 1)).days / 3652.5 for d in starts]   # per decade
    return np.array(e)[:, None]


def archive_with_era(basin, mode):
    M, n = K.archive_weekly(basin)
    E = era_column(mode)
    if E is not None:
        M.D = CH.Design(M.M, M.x, M.D.Y, M.D.LP, E=E)
        M.E = E
    return M, n


def main():
    mode, out, nperm, nboot = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    rng = np.random.default_rng(K.SEED)       # same seed as the original archive run
    rows = []
    for basin in C.BASINS:
        M, n = archive_with_era(basin, mode)
        names = ["h1", "h2", "h3", "b", "m"]
        df = K.analyse(M, names, nperm, nboot, rng)
        df.insert(0, "basin", basin); df.insert(1, "mode", mode)
        # era coefficient for each class (log rate, Poisson fit on all 22 seasons)
        if mode != "none":
            rws = M.rows(list(range(M.ns)))
            beta = {k: C.pois(M.D.X(k, rws), M.D.Y[k][rws]) for k in names}
            ie = M.D.X("h1", rws).shape[1] - 2           # era column sits just before the index column
            df["era_coef"] = [beta[k][ie] if k in beta else np.nan for k in df.stat]
            df["corr_x_era"] = float(np.corrcoef(M.x, M.E[:, 0])[0, 1])
        rows.append(df)
        print(basin, mode, n, flush=True)
    pd.concat(rows).to_csv(out, index=False, float_format="%.6g")


if __name__ == "__main__":
    main()
