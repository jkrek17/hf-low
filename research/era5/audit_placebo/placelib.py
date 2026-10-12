"""Audit pack helpers (ERA5 PROXY, pipeline A). Plan: PREREGISTRATION.md (23f1b5d)."""
import os, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO, "research", "era5", "hemispheric"))
sys.path.insert(0, os.path.join(REPO, "research", "era5", "hem_channels"))
import hemlib as H  # noqa: E402
import chanlib as C  # noqa: E402

PCS = [f"z{i}" for i in range(1, 11)] + [f"u{i}" for i in range(1, 11)]
SEASONS = C.SEASONS
NW = C.NWEEK


def load_table():
    t = pd.read_csv(os.path.join(REPO, "research/era5/hemispheric/results/weekly_table.csv.gz"))
    t = t[(t.season >= 2004)].reset_index(drop=True)
    assert len(t) == 660, len(t)
    return t


def index_from_outcome(t, basin, y, prev, lams=H.LAMBDAS, verbose=False):
    """Leave-one-season-out index: for each season the penalised-PC linear predictor from a fit on the other 21."""
    df = pd.DataFrame(dict(season=t.season.values, month=t.month.values, y=y.astype(float), lprev=np.log1p(prev)))
    pcs = t[PCS].values
    idx = np.zeros(len(df))
    lam_used = []
    for s in SEASONS:
        tr, te = (df.season != s).values, (df.season == s).values
        lam, _ = H.loso_lambda(df[tr], [], pcs[tr], lambdas=lams)
        m = H.fit_model(df[tr], [], pcs[tr], lam)
        _, pc = H.predict(m, df[te], pcs[te], parts=True)
        idx[te] = pc
        lam_used.append(lam)
        if verbose:
            print(" ", s, lam, flush=True)
    return (idx - idx.mean()) / idx.std(), lam_used
