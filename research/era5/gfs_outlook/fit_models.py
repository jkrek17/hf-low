"""Fit and freeze the outlook models S, A, B (and their no-dp12 fallbacks). Writes results/outlook_models.json.
Plan: PREREGISTRATION.md. ERA5 pipeline A (a proxy) fixes, seasons 2004-2020 only, 00/12 UTC, target gcls
(four-class peak g800 over the next two 00/12 UTC fixes, cuts 38/54/71.7 kt; bomb_env/score.load()).
Multinomial logistic, L2 C = 1, never tuned; Prep = intensity/model.py (median impute, clip 0.5-99.5 %, standardise).
Predictors without young/logage: S = msl dp12 lat speed g800 pac doy_c doy_s; A = S + jet250 vadv500; B = S + jet250 eady.
dp12 models: rows with dp12 present (young == 0). Fallback models (no dp12): all rows, dp12 dropped.
Nothing is scored here. usage: fit_models.py
"""
import os, sys, json
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "bomb_env"))
import score as BS                     # bomb_env/score.py: load(), CUTS; imports intensity/model as M
M = BS.M

TRAIN_TO = 2020
S_COLS = ["msl", "dp12", "lat", "speed", "g800", "pac", "doy_c", "doy_s"]
SETS = {"S": S_COLS, "A": S_COLS + ["jet250", "vadv500"], "B": S_COLS + ["jet250", "eady"]}


def build(cols, D):
    X = D[cols].values.astype(float)
    y = D["gcls"].astype(int).values
    prep = M.Prep(X)
    m = M.fit(prep(X), y)
    assert list(m.classes_) == [0, 1, 2, 3]
    return dict(predictors=cols, prep=prep.to_json(), coef=m.coef_.tolist(), intercept=m.intercept_.tolist(),
                classes=[int(c) for c in m.classes_], C=1.0, n_train=int(len(y)),
                n_by_class=np.bincount(y, minlength=4).tolist(),
                n_train_by_basin={b: int((D.basin == b).sum()) for b in ("atl", "pac")},
                n_env_missing={c: int(D[c].isna().sum()) for c in cols if c in ("jet250", "vadv500", "eady")})


def main():
    D = BS.load()
    n_load = len(D)
    D = D[(D.season >= M.FIT_FROM) & (D.season <= TRAIN_TO)].reset_index(drop=True)
    has_dp = D.young == 0
    out = dict(
        description="Frozen outlook models. ERA5 pipeline A (a proxy) fixes, seasons 2004-2020, target gcls (4 classes). Not scored on anything.",
        use="p = softmax(coef . standardise(clip(impute(x))) + intercept). Lead 0 or dp12 missing: models with suffix _nodp12; else dp12 models.",
        rows_loaded_2004_2025=int(n_load), rows_train_2004_2020=int(len(D)), rows_dropped_dp12_missing=int((~has_dp).sum()),
        class_names=BS.NAMES, cuts_kt=list(BS.CUTS), models={})
    for k, cols in SETS.items():
        out["models"][k] = build(cols, D[has_dp])
        out["models"][k + "_nodp12"] = build([c for c in cols if c != "dp12"], D)
    clim = []
    for (b, mo), g in D.groupby(["basin", "month"]):
        cnt = np.bincount(g.gcls.astype(int), minlength=4)
        sm = (cnt + 0.5) / (cnt + 0.5).sum()
        clim.append(dict(basin=b, month=int(mo), n=int(len(g)), counts=cnt.tolist(), probs_smoothed=sm.tolist()))
    out["climatology_gcls"] = clim
    json.dump(out, open(os.path.join(HERE, "results", "outlook_models.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in out.items() if k not in ("climatology_gcls", "models")}, indent=1))
    for k, m in out["models"].items():
        print(k, m["n_train"], m["n_by_class"], m["n_env_missing"])


if __name__ == "__main__":
    main()
