"""Storm-centred 500 nm patches for every fix in results/fixes_patch.csv.

usage: ERA5_WORK=DIR extract_patches.py [NPROC]     (needs DIR/clim_patch.npz from clim_patch.py)
Per WeatherBench2 time chunk one npz in DIR/patches (resumable); stops if the running pull passes MAXGB.
Patch: 13 x 13 points, 154.33 km spacing (edge 926 km = 500 nm), north-up, bilinear from the 1.5 degree grid.
Channels CH: a_z500 (m), a_msl (hPa), a_u250, a_v250 (m/s) as anomalies from the 1990-2019 climatology, then raw z500, msl (display).
"""
import os, sys, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from multiprocessing import Pool
from common_patch import *

MAXGB = 30.0
_C = None


def clim():
    global _C
    if _C is None:
        z = np.load(os.path.join(WORK, "clim_patch.npz")); _C = {k: z[k] for k in z.files}
    return _C


def job(args):
    chunk, rows = args
    path = os.path.join(WORK, "patches", f"{chunk:05d}.npz")
    if os.path.exists(path):
        return chunk, 0
    spent = sum(float(open(f).read()) for f in glob.glob(os.path.join(WORK, "bytes_*.txt")))
    if spent / 1e9 > MAXGB:
        raise SystemExit(f"stop: {spent / 1e9:.1f} GB pulled (limit {MAXGB})")
    C, nb, A = clim(), 0, {}
    for var, tag in (("u_component_of_wind", "u"), ("v_component_of_wind", "v"), ("geopotential", "z")):
        a, n = decode(WB2, var, f"{chunk}.0.0.0"); A[tag], nb = a, nb + n
    a, n = decode(WB2, "mean_sea_level_pressure", f"{chunk}.0.0"); A["m"], nb = a, nb + n
    keys, out = [], []
    for r in rows.itertuples():
        when = np.datetime64(f"{r.time[:4]}-{r.time[4:6]}-{r.time[6:8]}T{r.time[8:10]}")
        step = int((when - WB2_T0) / np.timedelta64(6, "h")); assert step // 8 == chunk
        ti = step % 8
        h = int((when.astype("datetime64[h]") - when.astype("datetime64[D]")).astype(int)) // 6
        doy = int((when.astype("datetime64[D]") - when.astype("datetime64[Y]").astype("datetime64[D]")).astype(int)) + 1
        assert C["doy"][doy - 1] == doy and int(C["hour"][h]) == h * 6
        z = A["z"][ti, LEVELS.index(500)].T.astype(np.float64) / G
        m = A["m"][ti].T.astype(np.float64) / 100.0
        u = A["u"][ti, LEVELS.index(250)].T.astype(np.float64); v = A["v"][ti, LEVELS.index(250)].T.astype(np.float64)
        F = {"a_z500": z - C["z500"][h, doy - 1], "a_msl": m - C["msl"][h, doy - 1],
             "a_u250": u - C["u250"][h, doy - 1], "a_v250": v - C["v250"][h, doy - 1], "z500": z, "msl": m}
        la, lo = destination(r.lat, r.lon, 90.0 + PBANG, PR)
        out.append(np.stack([sample(F[c], la, lo) for c in CH]).astype(np.float16))
        keys.append(r.key)
    tmp = path + ".tmp.npz"
    np.savez_compressed(tmp, keys=np.array(keys), x=np.stack(out))
    os.replace(tmp, path)
    with open(os.path.join(WORK, f"bytes_{chunk:05d}.txt"), "w") as f:
        f.write(str(nb))
    return chunk, nb


if __name__ == "__main__":
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    os.makedirs(os.path.join(WORK, "patches"), exist_ok=True)
    D = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", "fixes_patch.csv"), dtype={"time": str})
    D["key"] = D.track.astype(str) + "|" + D.time
    t = pd.to_datetime(D.time, format="%Y%m%d%H").values.astype("datetime64[h]")
    D["chunk"] = ((t - WB2_T0) / np.timedelta64(6, "h")).astype(int) // 8
    J = sorted(((c, g) for c, g in D.groupby("chunk")), key=lambda x: -x[0])
    print(len(J), "chunks,", len(D), "rows; expected", round(len(J) * 29.08 / 1000, 1), "GB", flush=True)
    clim()
    tot, t0 = 0, __import__("time").time()
    with Pool(nproc) as p:
        for i, (c, nb) in enumerate(p.imap_unordered(job, J)):
            tot += nb
            if i % 25 == 0:
                print(f"{i + 1}/{len(J)} chunk {c} pulled {tot / 1e9:.1f} GB {__import__('time').time() - t0:.0f}s", flush=True)
    print("done, pulled this run", round(tot / 1e9, 2), "GB")
