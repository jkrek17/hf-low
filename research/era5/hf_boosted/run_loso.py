"""RA-21: leave-one-season-out probabilities for M0-M3, the S9 single-feature spline logistics and the S10 product-term logistics,
plus the L2 temporal split (train 2004-14, test 2015-25). Pipeline A, ERA5 proxy, PR 12's sample and folds. Settings are fixed in PREREG.md.

usage: run_loso.py REPO_ROOT WORKDIR      (WORKDIR is not committed; it holds the held-out probabilities)
Writes WORKDIR/probs.npz and WORKDIR/meta.csv. Prints timings only, no scores.
"""
import os, sys, time, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
os.environ.setdefault("OMP_NUM_THREADS", "1")
root, work = sys.argv[1], sys.argv[2]
os.makedirs(work, exist_ok=True)
sys.path.insert(0, os.path.join(root, "research/era5/intensity"))
import model as M
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import SplineTransformer
from joblib import Parallel, delayed

res = os.path.join(root, "research/era5/intensity/results")
D = M.load(os.path.join(res, "fixes_2004.csv.gz"), os.path.join(res, "env_2004.csv.gz"))
D["hf24"] = D.hf24.astype(int); D["hf_now"] = D.hf_now.astype(int)
COLS = M.SETS["full"]
CONT = ["msl", "dp12", "lat", "speed", "g800", "logage", "B", "VTL", "VTU", "jet250", "div300", "vadv500", "eady", "sst", "sstgrad", "sst_t500", "flux", "tcwv"]
PAIRS = [("g800", "dp12"), ("msl", "dp12"), ("jet250", "eady"), ("sstgrad", "flux"), ("lat", "g800")]
ci = {c: COLS.index(c) for c in COLS}
HGB = dict(max_iter=300, learning_rate=0.06, max_leaf_nodes=15, l2_regularization=1.0, early_stopping=False, random_state=0)
X = D[COLS].values.astype(float)
y = D["hf24"].values
seas = D.season.values


def clim(tr_idx, te_idx):
    """Basin-month climatology of hf24 from training rows, with PR 12's +0.5 smoothing."""
    Dt = D.iloc[tr_idx]
    g = Dt.groupby(["basin", "month"]).hf24.agg(["sum", "count"])
    key = list(zip(D.basin.values[te_idx], D.month.values[te_idx]))
    return np.array([(g.loc[k, "sum"] + .5) / (g.loc[k, "count"] + 1) for k in key])


def spline_block(Ztr, Zte, idx):
    st = SplineTransformer(n_knots=5, degree=3, knots="quantile", extrapolation="constant", include_bias=False)
    return st.fit(Ztr[:, idx]).transform(Ztr[:, idx]), st.transform(Zte[:, idx]), st


def lr(A, yy):
    return LogisticRegression(C=1.0, max_iter=3000).fit(A, yy)


def fold(tr, te, light=False):
    """All models for one train/test split. Returns dict name -> (held-out prob, train prob or None)."""
    out = {}
    Xtr, Xte, ytr = X[tr], X[te], y[tr]
    prep = M.Prep(Xtr); Ztr, Zte = prep(Xtr), prep(Xte)
    m0 = lr(Ztr, ytr)
    out["M0"] = (m0.predict_proba(Zte)[:, 1], m0.predict_proba(Ztr)[:, 1])
    m1 = HistGradientBoostingClassifier(**HGB).fit(Xtr, ytr)
    out["M1"] = (m1.predict_proba(Xte)[:, 1], m1.predict_proba(Xtr)[:, 1])
    out["clim"] = (clim(tr, te), None)
    if light:
        return out
    m2 = HistGradientBoostingClassifier(interaction_cst=[[i] for i in range(len(COLS))], **HGB).fit(Xtr, ytr)
    out["M2"] = (m2.predict_proba(Xte)[:, 1], None)
    cidx = [ci[c] for c in CONT]
    Str, Ste, _ = spline_block(Ztr, Zte, cidx)
    m3 = lr(np.hstack([Ztr, Str]), ytr)
    out["M3"] = (m3.predict_proba(np.hstack([Zte, Ste]))[:, 1], None)
    for c in CONT:
        Str, Ste, _ = spline_block(Ztr, Zte, [ci[c]])
        m = lr(np.hstack([Ztr, Str]), ytr)
        out["S9_" + c] = (m.predict_proba(np.hstack([Zte, Ste]))[:, 1], None)
    for a, b in PAIRS:
        ptr = (Ztr[:, ci[a]] * Ztr[:, ci[b]])[:, None]; pte = (Zte[:, ci[a]] * Zte[:, ci[b]])[:, None]
        m = lr(np.hstack([Ztr, ptr]), ytr)
        out["S10_%s_x_%s" % (a, b)] = (m.predict_proba(np.hstack([Zte, pte]))[:, 1], None)
    return out


def run_season(s):
    t0 = time.time()
    tr = np.where(seas != s)[0]; te = np.where(seas == s)[0]
    o = fold(tr, te)
    print("season", s, "done in %.0f s" % (time.time() - t0), flush=True)
    return s, te, o


if __name__ == "__main__":
    t0 = time.time()
    seasons = sorted(np.unique(seas))
    parts = Parallel(n_jobs=4)(delayed(run_season)(s) for s in seasons)
    names = list(parts[0][2].keys())
    P = {n: np.full(len(D), np.nan, dtype=np.float32) for n in names}
    CUT = {n: np.full(len(D), np.nan, dtype=np.float32) for n in ("M0", "M1")}
    for s, te, o in parts:
        for n in names:
            P[n][te] = o[n][0]
        for n in ("M0", "M1"):
            CUT[n][te] = np.quantile(o[n][1], 1 - y[seas != s].mean())
    # L2: temporal split
    tr = np.where(seas <= 2014)[0]; te = np.where(seas > 2014)[0]
    o = fold(tr, te, light=True)
    T = {"L2_" + n: np.full(len(D), np.nan, dtype=np.float32) for n in o}
    for n in o:
        T["L2_" + n][te] = o[n][0]
    np.savez_compressed(os.path.join(work, "probs.npz"), **P, **{"cut_" + k: v for k, v in CUT.items()}, **T)
    D[["track", "time", "basin", "season", "month", "hf24", "hf_now", "g800"]].to_csv(os.path.join(work, "meta.csv"), index=False)
    print("all done in %.0f s; models: %s" % (time.time() - t0, ", ".join(names)))
