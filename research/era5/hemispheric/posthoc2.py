"""POST HOC PH5: can any nonlinear function of the named indices reproduce the out-of-sample pattern index?
Random forest on (NAO, PNA, ONI, MJO1, MJO2, AO, month), season-grouped 11-fold CV R-squared. Not pre-registered.
usage: posthoc2.py WORKDIR [OUT]"""
import os, sys, json, datetime as dt
import numpy as np, pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GroupKFold
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hemlib as H
import run as R

work = sys.argv[1]
out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(H.HERE, "results")
ctx = R.Ctx(work)
ALL = list(range(2004, 2026))
ix = ctx.index_table(ALL, "lag1")
res = {}
for b in H.BASINS:
    o = pd.read_csv(os.path.join(out, f"oos_index_{b}.csv"))
    df = pd.DataFrame(dict(season=o.season, week=o.week))
    df, _ = R.add_indices(df, ix)
    df["month"] = [(H.week_start(s, k) + dt.timedelta(days=3)).month for s, k in zip(df.season, df.week)]
    X = df[R.NAMED + ["ao", "month"]].values; y = o.idx.values
    pred = np.zeros(len(y))
    for tr, te in GroupKFold(11).split(X, y, df.season):
        rf = RandomForestRegressor(300, min_samples_leaf=10, random_state=H.SEED, n_jobs=4).fit(X[tr], y[tr])
        pred[te] = rf.predict(X[te])
    res[b] = dict(cv_r2_random_forest=float(1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()))
    print(b, res[b])
json.dump(res, open(os.path.join(out, "posthoc_rf.json"), "w"), indent=1)
