"""Per-fix J2/T2 features for Oct-Apr fixes (plus 28-30 Sep as lag sources). Resumable; counts bytes read.

usage: ERA5_WORK=DIR extract.py [NPROC]       output DIR/feat/<chunk>.csv (track, time, features); DIR/bytes_<chunk>.txt
Sensitivity columns (suffix _r1000, _r0: vortex removal 1000 km / none; _x500, _x1000: trough exclusion) exist only for chunks
processed after the first 304 (newest seasons; re-pulling them would pass the 60 GB limit).
Stops if the pull passes MAXGB (60 GB). Only features are kept; no raw field is saved. Outcome columns are not read.
"""
import os, sys, glob
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np, pandas as pd
from multiprocessing import Pool
from features import *
from common import decode, WB2, WB2_T0, LEVELS, WORK

MAXGB = 60.0
FIX = os.path.join(HERE, "..", "intensity", "results", "fixes_2004.csv.gz")
L = {l: LEVELS.index(l) for l in (250, 300, 500)}


def select():
    D = pd.read_csv(FIX, usecols=["track", "time", "basin", "season", "lat", "lon", "heading"], dtype={"time": str})
    D = D[D.season <= 2021]
    t = pd.to_datetime(D.time, format="%Y%m%d%H")
    keep = t.dt.month.isin([10, 11, 12, 1, 2, 3, 4]) | ((t.dt.month == 9) & (t.dt.day >= 28))
    D = D[keep.values].copy()
    D["chunk"] = ((t[keep.values].values.astype("datetime64[h]") - WB2_T0) / np.timedelta64(6, "h")).astype(int) // 8
    D["step"] = ((t[keep.values].values.astype("datetime64[h]") - WB2_T0) / np.timedelta64(6, "h")).astype(int) % 8
    return D


def job(args):
    chunk, rows = args
    path = os.path.join(WORK, "feat", f"{chunk:05d}.csv")
    if os.path.exists(path):
        return chunk, 0
    spent = sum(float(open(f).read()) for f in glob.glob(os.path.join(WORK, "bytes_*.txt")))
    if spent / 1e9 > MAXGB:
        raise SystemExit(f"stop: {spent / 1e9:.1f} GB already pulled (limit {MAXGB})")
    A, nb = {}, 0
    for var, tag in (("u_component_of_wind", "u"), ("v_component_of_wind", "v"), ("geopotential", "z")):
        a, n = decode(WB2, var, f"{chunk}.0.0.0")
        A[tag], nb = a, nb + n
    recs, cache = [], {}
    for r in rows.itertuples():
        if r.step not in cache:
            ti = r.step
            cache[ti] = time_fields(A["u"][ti, L[250]], A["v"][ti, L[250]], A["u"][ti, L[300]], A["v"][ti, L[300]], A["z"][ti, L[500]])
        f = fix_features(cache[r.step], r.lat, r.lon, r.heading)
        f.update(track=r.track, time=r.time)
        for tag, kw, cols in (("r1000", dict(rm=1000), J2 + EXTRA), ("r0", dict(rm=0), J2 + EXTRA),
                              ("x500", dict(excl=500), T2), ("x1000", dict(excl=1000), T2)):
            g = fix_features(cache[r.step], r.lat, r.lon, r.heading, **kw)
            f.update({f"{c}_{tag}": g[c] for c in cols})
        recs.append(f)
    pd.DataFrame(recs).to_csv(path + ".tmp", index=False, float_format="%.5g")
    os.replace(path + ".tmp", path)
    with open(os.path.join(WORK, f"bytes_{chunk:05d}.txt"), "w") as fh:
        fh.write(str(nb))
    return chunk, nb


if __name__ == "__main__":
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    os.makedirs(os.path.join(WORK, "feat"), exist_ok=True)
    D = select()
    J = sorted(((c, g) for c, g in D.groupby("chunk")), key=lambda x: -x[0])
    print(len(J), "chunks,", len(D), "fixes; expected", round(len(J) * 28.53 / 1000, 1), "GB", flush=True)
    tot, t0 = 0, __import__("time").time()
    with Pool(nproc) as p:
        for i, (c, nb) in enumerate(p.imap_unordered(job, J)):
            tot += nb
            if i % 25 == 0:
                print(f"{i + 1}/{len(J)} chunk {c} pulled {tot / 1e9:.1f} GB {__import__('time').time() - t0:.0f}s", flush=True)
    print("done, pulled this run", round(tot / 1e9, 2), "GB")
