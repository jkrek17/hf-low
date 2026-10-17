"""Sample definitions fixed in PREREGISTRATION.md."""
import numpy as np, pandas as pd
STREAM = ["41002", "41001", "41025"]; SHELF = ["41013", "44014", "41004", "44009"]
CUR_TO = {"41002": 45.0, "41001": 65.0, "41025": 30.0}   # direction the current flows toward, deg true
def angdiff(a, b): return np.abs((a - b + 180) % 360 - 180)
def wind_sample(d):
    d = d.dropna(subset=["v10", "dT", "vg"]); return d[(d.vg >= 15) & (d.v10 >= 1)].copy()
def wave_sample(d):
    d = d.dropna(subset=["v10", "wdir", "wvht_ft", "dpd", "wtmp"]).copy()
    d = d[d.v10 >= 20]
    ok = np.ones(len(d), bool)
    for k in (1, 2, 3):
        ws, wd = d[f"ws_m{k}"].values, d[f"wd_m{k}"].values
        ok &= ~np.isnan(ws) & ~np.isnan(wd) & (np.abs(ws - d.v10.values) <= 0.3 * d.v10.values) & (angdiff(wd, d.wdir.values) <= 30)
    d = d[ok]
    cur = d.station.map(lambda s: CUR_TO.get(s, CUR_TO["41002"]))   # shelf stations use 41002's labels (placebo)
    to = angdiff(d.wdir.values, cur.values)          # wind FROM within 45 deg of where the current flows TO = opposed
    d["sector"] = np.where(to <= 45, "opposed", np.where(angdiff(d.wdir.values, (cur.values + 180) % 360) <= 45, "following", "cross"))
    d["instream"] = d.wtmp >= 22
    d["vbin"] = pd.cut(d.v10, [20, 30, 40, 200], right=False, labels=["20-30", "30-40", "40+"])
    return d
