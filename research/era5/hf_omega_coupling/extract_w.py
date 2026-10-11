"""Per-fix W/C features for the same Oct-Apr fixes as hf_jet_trough (see PREREGISTRATION.md). Resumable; counts bytes read.

usage: ERA5_WORK=DIR extract_w.py [NPROC] [MAX_CHUNKS]     output DIR/featw/<chunk>.csv; DIR/bytes_<chunk>.txt
Reads only WeatherBench2 vertical_velocity (11.3 MB per chunk). Stops if the pull passes MAXGB (24 GB). No raw field is saved;
outcome columns are not read.
"""
import os, sys, glob
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "hf_env_composites"))
import numpy as np, pandas as pd
from multiprocessing import Pool
from omega_features import *
from common import decode, WB2, WB2_T0, LEVELS, WORK

MAXGB = 24.0
FIX = os.path.join(HERE, "..", "intensity", "results", "fixes_2004.csv.gz")
L5, L7 = LEVELS.index(500), LEVELS.index(700)


def select():
    D = pd.read_csv(FIX, usecols=["track", "time", "basin", "season", "lat", "lon", "heading"], dtype={"time": str})
    D = D[D.season <= 2021]
    t = pd.to_datetime(D.time, format="%Y%m%d%H")
    keep = t.dt.month.isin([10, 11, 12, 1, 2, 3, 4]) | ((t.dt.month == 9) & (t.dt.day >= 28))
    D = D[keep.values].copy()
    k = ((t[keep.values].values.astype("datetime64[h]") - WB2_T0) / np.timedelta64(6, "h")).astype(int)
    D["chunk"], D["step"] = k // 8, k % 8
    return D


def job(args):
    chunk, rows = args
    path = os.path.join(WORK, "featw", f"{chunk:05d}.csv")
    if os.path.exists(path):
        return chunk, 0
    spent = sum(float(open(f).read()) for f in glob.glob(os.path.join(WORK, "bytes_*.txt")))
    if spent / 1e9 > MAXGB:
        raise SystemExit(f"stop: {spent / 1e9:.1f} GB already pulled (limit {MAXGB})")
    a, nb = decode(WB2, "vertical_velocity", f"{chunk}.0.0.0")
    recs, cache = [], {}
    for r in rows.itertuples():
        if r.step not in cache:
            cache[r.step] = time_fields_w(a[r.step, L5], a[r.step, L7])
        f = fix_features_w(cache[r.step], r.lat, r.lon, r.heading)
        f.update(track=r.track, time=r.time)
        recs.append(f)
    pd.DataFrame(recs).to_csv(path + ".tmp", index=False, float_format="%.5g")
    os.replace(path + ".tmp", path)
    with open(os.path.join(WORK, f"bytes_{chunk:05d}.txt"), "w") as fh:
        fh.write(str(nb))
    return chunk, nb


if __name__ == "__main__":
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    maxc = int(sys.argv[2]) if len(sys.argv) > 2 else None
    os.makedirs(os.path.join(WORK, "featw"), exist_ok=True)
    D = select()
    J = sorted(((c, g) for c, g in D.groupby("chunk")), key=lambda x: -x[0])
    print(len(J), "chunks,", len(D), "fixes; expected", round(len(J) * 11.28 / 1000, 1), "GB", flush=True)
    if maxc:
        J = J[:maxc]
    tot, t0 = 0, __import__("time").time()
    with Pool(nproc) as p:
        for i, (c, nb) in enumerate(p.imap_unordered(job, J)):
            tot += nb
            if i % 25 == 0:
                print(f"{i + 1}/{len(J)} chunk {c} pulled {tot / 1e9:.2f} GB {__import__('time').time() - t0:.0f}s", flush=True)
    print("done, pulled this run", round(tot / 1e9, 2), "GB")
