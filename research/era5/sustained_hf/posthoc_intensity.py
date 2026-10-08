"""POST HOC (not in PREREGISTRATION.md): is the rise of the pattern effect with duration a duration effect or an intensity effect?
Paired difference log RR(HF_k by duration) - log RR(count-matched peak-gust cut, same number of tracks per basin over 2004-2025), k = 2, 3.
Same weekly (PR 64 / PR 41-index) and daily (PR 14) designs and the same seed machinery as contrast.py. One run; labelled post hoc wherever quoted.
usage: posthoc_intensity.py DESIGN OUT_CSV NPERM NBOOT      DESIGN = weekly | daily
"""
import os, sys
import numpy as np, pandas as pd
import contrast as K

C, S = K.C, K.S
design, out, nperm, nboot = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
W = os.path.join(K.HERE, "work")
G = {k: pd.read_csv(os.path.join(W, f"all_tracks_gcut{k}.csv.gz"), usecols=["track", "gust800_kt"]).set_index("track").gust800_kt > 0 for k in (2, 3)}
rng = np.random.default_rng(K.SEED + 1)
extra = [("h2", "g2"), ("h3", "g3")]
rows = []
if design == "weekly":
    Tfull = K.tracks_with_nhf(2003, 2025)
    for k in (2, 3):
        Tfull[f"g{k}"] = Tfull.track.map(G[k]).fillna(False).astype(bool)
    for basin in C.BASINS:
        x, _ = C.index(basin)
        t_all = Tfull[Tfull.basin == basin]
        win = C.assign_week(t_all)
        Y, LP = {}, {}
        for name in ("h2", "h3", "g2", "g3"):
            m = (t_all[name] if name[0] == "g" else pd.Series(K.CLASSES[name](t_all.hf.values == 1, t_all.n_hf.values), index=t_all.index))
            m = m.astype(bool)
            Y[name] = C.weekly_counts(win, m.reindex(win.index))
            LP[name] = np.log1p(C.prev_week_counts(t_all.assign(_m=m), m, basin))
        M = K.Weekly(x, Y, LP)
        df = K.analyse(M, ["h2", "h3", "g2", "g3"], nperm, nboot, rng, extra_pairs=extra)
        df.insert(0, "basin", basin); rows.append(df)
else:
    T = S.load_tracks(os.path.join(K.ERA, "hf_history/results/all_tracks.csv.gz")).merge(K.NHF, on="track", how="left")
    T["n_hf"] = T.n_hf.fillna(0)
    I = S.indices(K.CPC, K.REPO)
    for k in (2, 3):
        T[f"g{k}"] = T.track.map(G[k]).fillna(False).astype(bool)
    for basin, idx in (("atl", "NAO"), ("pac", "PNA")):
        outcomes = ["h2", "h3"]
        M = K.Daily.__new__(K.Daily)
        D = S.window_days("octapr", 2004, 2025)
        Tw = S.assign(T, "octapr", 2004, 2025); Tw = Tw[Tw.basin == basin]
        col = idx + "_lag"
        raw = I.reindex(D.date)[[col]].reset_index(drop=True)
        ok = raw.notna().all(axis=1).values
        D, raw = D[ok].reset_index(drop=True), raw[ok].reset_index(drop=True)
        z = (raw - raw.mean()) / raw.std(ddof=0)
        for o in ("h2", "h3"):
            Tw[o] = K.CLASSES[o](Tw.hf.values == 1, Tw.n_hf.values).astype(int)
        for k in (2, 3):
            Tw[f"g{k}"] = Tw[f"g{k}"].astype(int)
        M.C = S.daily_counts(Tw, D, basin, ["h2", "h3", "g2", "g3"]); M.X = S.design(M.C, z, [col])
        M.season = M.C.season.values; M.seasons = np.arange(2004, 2026); M.Xv = np.asarray(M.X, float); M.ns = 22
        df = K.analyse(M, ["h2", "h3", "g2", "g3"], nperm, nboot, rng, extra_pairs=extra)
        df.insert(0, "basin", basin); rows.append(df)
pd.concat(rows).to_csv(out, index=False, float_format="%.6g")
