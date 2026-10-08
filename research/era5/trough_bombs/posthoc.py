"""POST HOC (not in PREREGISTRATION.md): does the size of the D2 effect, and the D1 effect, survive a control for
where the storm is? D2 (zonal eddy) keeps the stationary-wave trough, so it may partly mark geography.
Controls: fixed effects for the 5.625 degree position cell of the fix (cells with >= 100 fixes; others pooled).
Everything here is labelled post hoc in the README; the pre-registered verdict is not changed by it.

usage: posthoc.py WORK_DIR -> results/posthoc.csv
"""
import sys
import numpy as np, pandas as pd
import tblib as T

N, R = T.load_tables(sys.argv[1])
std = T.Std(R)
rows = []
for d, col in T.DEFS.items():
    S = N[N.elig & N[col].notna()].reset_index(drop=True)
    cell = (np.floor(S.lat / 5.625).astype(int).astype(str) + "_" + np.floor(S.lon / 5.625).astype(int).astype(str))
    vc = cell.value_counts()
    cell = cell.where(cell.map(vc) >= 100, "other")
    Dm = pd.get_dummies(cell, drop_first=True).values.astype(float)
    base = T.design(S, std, col)
    for lab, X in (("pre-registered model", base),
                   ("+ position-cell fixed effects", np.column_stack([base[:, :-1], Dm, base[:, -1]])),
                   ("+ position cells + jet_up", np.column_stack([base[:, :-1], Dm, std.z(S, "jet_up"), base[:, -1]]))):
        y = S.bomb.values
        r = T.summarize(S, y, X)
        r.update(definition=d, label=lab, cells=Dm.shape[1] + 1 if "cell" in lab else 0)
        rows.append(r)
        print(d, lab, round(r["OR"], 3), round(r["lo"], 3), round(r["hi"], 3), flush=True)
    # raw (unadjusted) trough-variable correlation with latitude and longitude position, descriptive
    z = std.z(S, col)
    print(d, "corr of z(trough) with lat", round(np.corrcoef(z, S.lat)[0, 1], 3), "with doy_c", round(np.corrcoef(z, S.doy_c)[0, 1], 3))
pd.DataFrame(rows)[["definition", "label", "OR", "lo", "hi", "p", "n", "events", "seasons", "cells"]].to_csv("results/posthoc.csv", index=False)
