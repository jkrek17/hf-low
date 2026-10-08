"""Fit and freeze the Stage A model (PREREGISTRATION.md, Stage A2). Writes results/gfs_model.json.

PR 12 'state' model (intensity/model.py: Prep = median impute, clip 0.5-99.5 %, standardise; LogisticRegression L2 C=1, never tuned)
WITHOUT the history terms a forecast track cannot have (young, logage). Training rows: ERA5 pipeline A fixes_2004.csv.gz joined to
env_2004.csv.gz exactly as model.load does, seasons 2004..2020 only, target hf24. Nothing is scored here, not even the training rows.
  dp12 model   (lead >= 12 h): msl, dp12, lat, speed, g800, pac, doy_c, doy_s; rows with dp12 missing dropped
  nodp12 model (lead 0 only) : msl, lat, speed, g800, pac, doy_c, doy_s;       all training rows
doy_c, doy_s = cos, sin of 2 pi dayofyear / 365.25 of the fix valid time; pac = 1 for the Pacific box.
Climatology: hf24 frequency by basin x calendar month over all training rows (raw rate, and the +0.5 smoothed value model.py uses).
usage: a_fit_model.py [FIXES ENV]
"""
import os, sys, json
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "intensity"))
import model as M

IR = os.path.join(HERE, "..", "intensity", "results")
FIXES = sys.argv[1] if len(sys.argv) > 1 else os.path.join(IR, "fixes_2004.csv.gz")
ENV = sys.argv[2] if len(sys.argv) > 2 else os.path.join(IR, "env_2004.csv.gz")
TRAIN_TO = 2020
P_DP = ["msl", "dp12", "lat", "speed", "g800", "pac", "doy_c", "doy_s"]
P_NO = [c for c in P_DP if c != "dp12"]


def build(cols, D):
    X = D[cols].values.astype(float)
    y = D["hf24"].astype(int).values
    prep = M.Prep(X)
    m = M.fit(prep(X), y)
    return dict(predictors=cols, prep=prep.to_json(), coef=m.coef_.tolist(), intercept=m.intercept_.tolist(),
                C=1.0, n_train=int(len(y)), n_hf24=int(y.sum()),
                n_train_by_basin={b: int((D.basin == b).sum()) for b in ("atl", "pac")},
                n_hf24_by_basin={b: int(y[(D.basin == b).values].sum()) for b in ("atl", "pac")},
                seasons=[int(D.season.min()), int(D.season.max())])


def main():
    D = M.load(FIXES, ENV)                       # seasons 2004..2025, as in model.py
    n_load = len(D)
    D = D[(D.season >= M.FIT_FROM) & (D.season <= TRAIN_TO)].reset_index(drop=True)
    D["hf24"] = D.hf24.astype(int)
    has_dp = D.young == 0
    out = dict(
        description="Frozen Stage A model. Pipeline A (ERA5, a proxy) fixes, seasons 2004-2020, target hf24. Not scored on anything.",
        use="lead 0 h: model_nodp12 ; lead >= 12 h: model_dp12. p = sigmoid(coef . standardise(clip(impute(x))) + intercept); "
            "see intensity/model.py predict(). dp12 = msl(t) - msl(t-12 h) hPa along the track as in pipeline A.",
        rows_loaded_2004_2025=int(n_load), rows_train_2004_2020=int(len(D)),
        model_dp12=build(P_DP, D[has_dp]),
        model_nodp12=build(P_NO, D),
    )
    out["model_dp12"]["rows_dropped_dp12_missing"] = int((~has_dp).sum())
    clim = []
    for (b, mo), g in D.groupby(["basin", "month"]):
        k, n = int(g.hf24.sum()), len(g)
        clim.append(dict(basin=b, month=int(mo), n=n, n_hf24=k, rate=k / n, rate_smoothed=(k + 0.5) / (n + 1.0)))
    out["climatology_hf24"] = clim
    json.dump(out, open(os.path.join(HERE, "results", "gfs_model.json"), "w"), indent=1)
    # plumbing check only (no scores): the saved JSON reproduces sklearn's own probabilities on 5 rows
    for key, S in (("model_dp12", D[has_dp]), ("model_nodp12", D)):
        mm = out[key]
        mod = dict(mm); mod["coef"] = mm["coef"]
        p = M.predict(mod, S.head(5))
        assert p.shape == (5,) and np.all((p > 0) & (p < 1))
    print(json.dumps({k: v for k, v in out.items() if not k.startswith("climatology")}, indent=1)[:6000])


if __name__ == "__main__":
    main()
