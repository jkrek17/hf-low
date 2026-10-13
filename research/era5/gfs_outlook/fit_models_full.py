"""Fit and freeze the FULL outlook models (FULL and FULL_NOSF, each with a no-dp12 fallback). Writes results/outlook_models_full.json.
Plan: PREREGISTRATION_FULL.md (part 2). ERA5 pipeline A (a proxy) fixes, seasons 2004-2020 only, 00/12 UTC, target gcls (4 classes).
Same row selection, gcls, C = 1, intensity/model.py Prep as fit_models.py (imported, not copied). No young / logage.
FULL      = S + Hart (B VTL VTU) + jet250 div300 vadv500 eady sst sstgrad sst_t500 flux tcwv nosst
FULL_NOSF = FULL without sst sstgrad sst_t500 flux nosst
dp12 models: rows with dp12 present (young == 0); fallback (_nodp12): all rows, dp12 dropped. Nothing is scored here. outlook_models.json is not touched.
usage: fit_models_full.py
"""
import os, json
import numpy as np
import fit_models as FM

HERE = FM.HERE
BS, M = FM.BS, FM.M
HART = ["B", "VTL", "VTU"]
ENV = ["jet250", "div300", "vadv500", "eady", "sst", "sstgrad", "sst_t500", "flux", "tcwv", "nosst"]
SF = ["sst", "sstgrad", "sst_t500", "flux", "nosst"]
FULL = FM.S_COLS + HART + ENV
SETS = {"FULL": FULL, "FULL_NOSF": [c for c in FULL if c not in SF]}


def main():
    D = BS.load()
    n_load = len(D)
    D = D[(D.season >= M.FIT_FROM) & (D.season <= FM.TRAIN_TO)].reset_index(drop=True)
    has_dp = D.young == 0
    out = dict(
        description="Frozen FULL outlook models (part 2). ERA5 pipeline A (a proxy) fixes, seasons 2004-2020, target gcls (4 classes). Not scored on anything.",
        use="p = softmax(coef . standardise(clip(impute(x))) + intercept). Lead 0 or dp12 missing: models with suffix _nodp12; else dp12 models.",
        rows_loaded_2004_2025=int(n_load), rows_train_2004_2020=int(len(D)), rows_dropped_dp12_missing=int((~has_dp).sum()),
        class_names=BS.NAMES, cuts_kt=list(BS.CUTS), models={})
    miss_cols = HART + ENV
    for k, cols in SETS.items():
        for name, cc, DD in ((k, cols, D[has_dp]), (k + "_nodp12", [c for c in cols if c != "dp12"], D)):
            m = FM.build(cc, DD)
            m["n_env_missing"] = {c: int(DD[c].isna().sum()) for c in cc if c in miss_cols}
            out["models"][name] = m
    out["climatology_gcls"] = json.load(open(os.path.join(HERE, "results", "outlook_models.json")))["climatology_gcls"]
    out["climatology_note"] = "identical to outlook_models.json (same rows, 2004-2020, basin x month, smoothed +0.5)"
    json.dump(out, open(os.path.join(HERE, "results", "outlook_models_full.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in out.items() if k not in ("climatology_gcls", "models")}, indent=1))
    for k, m in out["models"].items():
        print(k, m["n_train"], m["n_by_class"], len(m["predictors"]), m["n_env_missing"])


if __name__ == "__main__":
    main()
