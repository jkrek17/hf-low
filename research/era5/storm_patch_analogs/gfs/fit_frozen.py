"""Fit and freeze, BEFORE any GFS field is read: the ERA5 pattern features (PCA, 15 components) and the logistic models.
Training rows: pipeline A ERA5 fixes of seasons 2014-15..2020-21, Oct-Apr, 00/12 UTC, with patches from /mnt/project-files/storm-patches (PR 199).
Season 2021-22 (the eighth patch season) is excluded because it is also a GFS test season.
Models (same form as gfs_outlook/fit_models_full.py: fit_models.build, C = 1, Prep; dp12 model and no-dp12 fallback):
  FULL7        FULL predictors (S + Hart + environment, no young/logage)
  FULLPC7      FULL7 + 15 pattern components (full 13 x 13 patch)
  FULLPC7_OUT  FULL7 + 15 components of the patch with the inner 300 km removed
usage: fit_frozen.py   -> results/frozen_pca_full.npz, frozen_pca_out.npz, frozen_models.json
"""
import os, sys, glob, json, hashlib
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "..", "gfs_outlook"))
import numpy as np, pandas as pd
import fit_models as FM
import fit_models_full as FF
from sklearn.decomposition import PCA
from frozen import Frozen, PR, NP
BS, M = FM.BS, FM.M
SD = "/mnt/project-files/storm-patches"
RES = os.path.join(HERE, "results")


def era5_rows():
    D = BS.load()
    D["key"] = D.track.astype(str) + "|" + D.time
    P = {}
    for f in sorted(glob.glob(os.path.join(SD, "patches_20*.npz"))):
        z = np.load(f, allow_pickle=True)
        for k, x in zip(z["keys"], z["x"]):
            P[str(k)] = x
    D = D[D.key.isin(P)].reset_index(drop=True)
    X = np.stack([P[k] for k in D.key]).astype(np.float32)
    return D, X


def fit_pca(Xtr, mask, path, npc=15):
    sd = Xtr.std((0, 2, 3), keepdims=True) + 1e-9
    Z = np.nan_to_num(Xtr[:, :4] / sd[:, :4], nan=0.0)[:, :, mask].reshape(len(Xtr), -1)
    mu = Z.mean(0)
    p = PCA(npc, svd_solver="full", random_state=0).fit(Z - mu)
    np.savez(path, sd=sd, mask=mask, mu=mu, comp=p.components_.astype(np.float32), pm=p.mean_.astype(np.float32))
    return Frozen(path)


def main():
    D, X = era5_rows()
    tr = ((D.season >= 2014) & (D.season <= 2020)).values
    D, X = D[tr].reset_index(drop=True), X[tr]
    print("training fixes", len(D), "tracks", D.track.nunique(), "seasons", sorted(D.season.unique()))
    out = dict(description="Frozen ERA5 pattern models for the GFS test; fitted before any GFS field was read.",
               training="ERA5 pipeline A (proxy) fixes 2014-15..2020-21 Oct-Apr", rows=int(len(D)), models={})
    has_dp = (D.young == 0).values
    cols0 = FF.FULL
    for tag, mask in (("full", np.ones((NP, NP), bool)), ("out", PR >= 300.0)):
        fr = fit_pca(X, mask, os.path.join(RES, f"frozen_pca_{tag}.npz"))
        pc = fr.pcs(X)
        for i in range(pc.shape[1]):
            D[f"pc{i+1}_{tag}"] = pc[:, i]
    pcf = [f"pc{i+1}_full" for i in range(15)]; pco = [f"pc{i+1}_out" for i in range(15)]
    for name, cols in (("FULL7", cols0), ("FULLPC7", cols0 + pcf), ("FULLPC7_OUT", cols0 + pco)):
        out["models"][name] = FM.build(cols, D[has_dp])
        out["models"][name + "_nodp12"] = FM.build([c for c in cols if c != "dp12"], D)
        print(name, out["models"][name]["n_train"], out["models"][name + "_nodp12"]["n_train"])
    # basin-month climatology of the same rows is NOT used: scoring uses outlook_models_full.json climatology (2004-2020) as in PR 185
    json.dump(out, open(os.path.join(RES, "frozen_models.json"), "w"), indent=1)
    for f in ("frozen_pca_full.npz", "frozen_pca_out.npz", "frozen_models.json"):
        print(f, hashlib.sha256(open(os.path.join(RES, f), "rb").read()).hexdigest()[:16])


if __name__ == "__main__":
    main()
