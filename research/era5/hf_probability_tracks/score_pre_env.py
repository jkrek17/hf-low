"""Model N (full without the gust index) on the 1979-2003 catalog fixes, after the 111 GB WeatherBench2 environment pull.

usage: python3 -I score_pre_env.py ENV_DIR      (ENV_DIR = per-chunk CSVs written by ../intensity/env.py for results/fixes_pre2004 inputs)
Writes results/env_pre2004.csv.gz (the derived environment table, so this reproduces without the pull) and adds P_N to
results/fix_probs_pre2004.csv.gz. N is fitted on all 2004-2025 fixes exactly as in score.py. Not a test, a map/ranking input.
"""
import os, sys, glob
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import score as S
M = S.M
RES = S.RES

def main(envdir):
    out = os.path.join(RES, "env_pre2004.csv.gz")
    if envdir != "-":
        E = pd.concat([pd.read_csv(f, dtype={"time": str}) for f in sorted(glob.glob(os.path.join(envdir, "*.csv")))])
        E.round(4).to_csv(out, index=False)
    E = pd.read_csv(out, dtype={"time": str})
    P = pd.read_csv(os.path.join(RES, "fix_probs_pre2004.csv.gz"), dtype={"time": str})
    D = M.load(os.path.join(S.INT, "results", "fixes_2004.csv.gz"), os.path.join(S.INT, "results", "env_2004.csv.gz"))
    D["hf24"] = D.hf24.astype(int)
    ct = pd.read_csv(os.path.join(S.HFH, "era5_hf_catalog_tracks.csv"), dtype={"time": str})
    y_, m_ = ct.time.str[:4].astype(int), ct.time.str[4:6].astype(int)
    ct["season"] = np.where(m_ >= 6, y_, y_ - 1)
    Fx = S.FX.build(ct); Fx = Fx[Fx.season < 2004]
    Fx = S.prepare(Fx).merge(E, on=["track", "time"], how="left")
    Fx["nosst"] = Fx.sst.isna().astype(float)
    cols = S.SETS["N"]
    prep = M.Prep(D[cols].values.astype(float))
    m = M.fit(prep(D[cols].values.astype(float)), D.hf24.values)
    Fx["P_N"] = m.predict_proba(prep(Fx[cols].values.astype(float)))[:, 1]
    P = P.drop(columns=[c for c in ("P_N",) if c in P]).merge(Fx[["track", "time", "P_N"]], on=["track", "time"], how="left")
    P["P_N"] = P.P_N.round(5)
    P.to_csv(os.path.join(RES, "fix_probs_pre2004.csv.gz"), index=False)
    print(f"{len(Fx)} fixes, env missing {int(Fx.B.isna().sum())}, mean P_N events {P[P.role=='event'].P_N.mean():.3f}, nulls {P[P.role=='null_case'].P_N.mean():.3f}; "
          f"corr(P_N, P_S) {np.corrcoef(P.P_N.fillna(0), P.P_S)[0,1]:.3f}")

if __name__ == "__main__":
    main(sys.argv[1])
