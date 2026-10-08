"""POST HOC (not in PREREG.md; logged there): is the ~0.6 HSS ceiling of concurrent diagnosis a limit of the linear form
or of the fields? Same fixes as analysis.py (pipeline A, ERA5 proxy, 2004-05..2025-26). Grouped 11-fold CV by season
(two seasons per fold), one run, default HistGradientBoosting settings, nothing tuned.

usage: diagnosis_extra.py REPO_ROOT OUTDIR
"""
import os, sys, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
root, out = sys.argv[1], sys.argv[2]
sys.path.insert(0, os.path.join(root, "research/era5/intensity"))
import model as M
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.model_selection import GroupKFold

res = os.path.join(root, "research/era5/intensity/results")
D = M.load(os.path.join(res, "fixes_2004.csv.gz"), os.path.join(res, "env_2004.csv.gz"))
D["hf24"] = D.hf24.astype(int); D["hf_now"] = D.hf_now.astype(int)
T1 = pd.read_csv(os.path.join(root, "research/era5/intensity_extra/results/features_tier1.csv.gz"), dtype={"time": str})
D = D.merge(T1, on=["track", "time"], how="left")
full = M.SETS["full"]
nog = [c for c in full if c != "g800"]
t1 = [c for c in T1.columns if c not in ("track", "time")]
lines = []
def w(s=""):
    lines.append(s); print(s, flush=True)

def cv(df, cols, target, kind):
    p = np.zeros(len(df)); pcl = np.zeros(len(df))
    gk = GroupKFold(n_splits=11)
    for tr, te in gk.split(df, groups=df.season):
        Xtr, Xte = df.iloc[tr][cols].values.astype(float), df.iloc[te][cols].values.astype(float)
        ytr = df.iloc[tr][target].values
        if kind == "logit":
            prep = M.Prep(Xtr); m = M.fit(prep(Xtr), ytr); p[te] = m.predict_proba(prep(Xte))[:, 1]
        else:
            m = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.06, max_leaf_nodes=15, l2_regularization=1.0,
                                               early_stopping=False, random_state=0).fit(Xtr, ytr)
            p[te] = m.predict_proba(Xte)[:, 1]
        trd = df.iloc[tr]
        for (b, mo), g in df.iloc[te].groupby(["basin", "month"]):
            t_ = trd[(trd.basin == b) & (trd.month == mo)][target]
            pcl[df.index.get_indexer(g.index)] = (t_.sum() + .5) / (len(t_) + 1)
    y = df[target].values
    bss = 1 - ((p - y) ** 2).sum() / ((pcl - y) ** 2).sum()
    cut = np.quantile(p, 1 - y.mean())
    t = M.table(p >= cut, y.astype(bool))
    return bss, t, M.roc_auc_score(y, p)

for label, df in (("all 159,430 fixes", D), ("fixes with the tier-1 fields (seasons 2004-05..2021-22)", D[D.ws850.notna()].reset_index(drop=True))):
    w(f"== {label}; n {len(df)}")
    for target, sets in (("hf_now", [("full minus g800 (concurrent diagnosis)", nog), ("... plus tier-1 fields (850 hPa wind, ring gradient, high, jet, IVT, land fraction, ...)", nog + t1)]),
                         ("hf24", [("full (PR 12 predictors)", full), ("... plus tier-1 fields", full + t1)])):
        for name, cols in sets:
            if cols is not full and cols is not nog and len(df) == len(D) and any(c in t1 for c in cols):
                continue
            for kind in ("logit", "boosted trees"):
                b, t, a = cv(df, cols, target, kind)
                w(f"  {target:7s} {name:78s} {kind:13s} BSS {b:+.3f} AUC {a:.3f} POD {t['pod']:.2f} FAR {t['far']:.2f} HSS {t['hss']:.2f}")
open(os.path.join(out, "diagnosis_extra.txt"), "w").write("\n".join(lines) + "\n")
