"""Descriptive companion: HF onset against storm-force-only peak, Oct-Apr 2004-05..2021-22 (no model; see PREREGISTRATION.md).

usage: desc.py OUTDIR      HF = first 00/12 fix with gust index >= 71.7 kt; storm-force-only = tracks with max gust index in [54, 71.7) kt,
at the peak fix, weighted to the HF count in each basin x month cell. Season-block bootstrap (2000, seed 20261009), BH over all comparisons.
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np, pandas as pd
from qc import bh

FIX = os.path.join(HERE, "..", "intensity", "results", "fixes_2004.csv.gz")
QS = ["RE", "LE", "RX", "LX", "nojet", "vmaxp", "tdepth", "ttilt", "tphase", "div300max"]


def main(out):
    F = pd.read_csv(FIX, usecols=["track", "time", "basin", "season", "g800"], dtype={"time": str})
    F = F[F.season <= 2021]
    X = pd.read_csv(os.path.join(HERE, "results", "features.csv.gz"), dtype={"time": str})
    X = X.drop(columns=["basin", "season", "lat", "lon", "heading"])
    t = pd.to_datetime(F.time, format="%Y%m%d%H")
    F["month"] = t.dt.month
    oa = F.month.isin([10, 11, 12, 1, 2, 3, 4])
    G = F[oa]
    mx = G.groupby("track").g800.transform("max")
    first_hf = G[G.g800 >= 71.7].sort_values("time").groupby("track").head(1)
    sf_tr = G[(mx >= 54) & (mx < 71.7)]
    sf_pk = sf_tr.loc[sf_tr.groupby("track").g800.idxmax()]
    hf = first_hf.assign(grp="HF")
    sf = sf_pk.assign(grp="SF")
    A = pd.concat([hf, sf], ignore_index=True)
    cnt = A[A.grp == "HF"].groupby(["basin", "month"]).size()
    cs = A[A.grp == "SF"].groupby(["basin", "month"]).size()
    w = (cnt / cs).rename("w")
    A = A.merge(w, left_on=["basin", "month"], right_index=True, how="left")
    A["w"] = np.where(A.grp == "HF", 1.0, A.w.fillna(0.0))
    print("HF onsets", int((A.grp == "HF").sum()), "storm-force-only peaks", int((A.grp == "SF").sum()), "weighted", round(A[A.grp == "SF"].w.sum(), 1))
    rows = []
    rng = np.random.default_rng(20261009)
    for lag in (0, 12, 24, 48):
        key = (pd.to_datetime(A.time, format="%Y%m%d%H") - pd.Timedelta(hours=lag)).dt.strftime("%Y%m%d%H)")[:0]
        k2 = (pd.to_datetime(A.time, format="%Y%m%d%H") - pd.Timedelta(hours=lag)).dt.strftime("%Y%m%d%H")
        B = A[["track", "grp", "basin", "season", "w"]].assign(time=k2.values).merge(X, on=["track", "time"], how="left")
        B = B[B.vmax2500.notna()]
        for basin in ("atl", "pac"):
            S = B[B.basin == basin]
            us = np.sort(S.season.unique())
            for q in QS:
                v = S[q].values.astype(float)
                if q in ("vmaxp", "tdepth", "ttilt", "tphase", "div300max"):
                    ok = np.isfinite(v)
                else:
                    ok = np.ones(len(v), bool)
                if q in ("RE", "LE", "RX", "LX"):
                    ok = ok & (S.nojet.values == 0)
                H = (S.grp == "HF").values & ok
                Sf = (S.grp == "SF").values & ok
                arr = []
                for u in us:
                    m = (S.season == u).values
                    arr.append([np.nansum(v[m & H]), (m & H).sum(), np.nansum((v * S.w.values)[m & Sf]), S.w.values[m & Sf].sum()])
                arr = np.array(arr, float)

                def st(a):
                    return a[:, 0].sum() / max(a[:, 1].sum(), 1) - a[:, 2].sum() / max(a[:, 3].sum(), 1e-9)
                d0 = st(arr)
                bs = np.array([st(arr[rng.integers(0, len(us), len(us))]) for _ in range(2000)])
                p = 2 * min((bs <= 0).mean(), (bs >= 0).mean())
                rows.append(dict(basin=basin, lag_h=lag, quantity=q, hf_mean=arr[:, 0].sum() / max(arr[:, 1].sum(), 1),
                                 sf_mean=arr[:, 2].sum() / max(arr[:, 3].sum(), 1e-9), diff=d0, lo90=np.percentile(bs, 5), hi90=np.percentile(bs, 95),
                                 p=max(p, 1 / 2001), n_hf=int(arr[:, 1].sum()), n_sf_w=float(arr[:, 3].sum())))
    T = pd.DataFrame(rows)
    T["q"] = bh(T.p.values)
    T.to_csv(os.path.join(out, "descriptive.csv"), index=False, float_format="%.4g")
    sig = T[T.q < 0.05]
    with open(os.path.join(out, "descriptive.txt"), "w") as f:
        f.write(f"HF onsets {int((A.grp == 'HF').sum())}; storm-force-only peaks {int((A.grp == 'SF').sum())} (weighted to HF counts by basin x month)\n"
                f"{len(T)} comparisons; q < 0.05 in {len(sig)}\n")
        f.write(sig.to_string(index=False) if len(sig) else "none pass BH q < 0.05\n")
        f.write("\n\nall at lag 0 (HF onset vs storm-force-only peak):\n" + T[T.lag_h == 0].to_string(index=False))
    print(open(os.path.join(out, "descriptive.txt")).read())


if __name__ == "__main__":
    main(sys.argv[1])
