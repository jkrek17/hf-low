"""Storm-centred boxes of Z500, 250/500 hPa wind and divergence, raw and anomaly, for every row of storm_times.csv.

usage: ERA5_WORK=DIR extract.py [NPROC]       (needs DIR/clim.npz from clim.py)
One npz per WeatherBench2 time chunk in DIR/boxes, resumable. Per row two frames, each [18, 81, 81] float32:
  rot    +x along the previous-6-h motion, +y left; the wind components are along / cross; NaN if heading is missing
  north  +x east, +y north; the wind components are u / v
Field order FIELDS: raw z500 u250 v250 u500 v500 ws250 ws500 d250 d500, then the same nine as anomalies
(ERA5 minus the WeatherBench2 1990-2019 climatology; speed anomaly = |V| - |V_clim|; divergence anomaly = divergence
of the anomaly wind). Stops if the running pull passes MAXGB (48 GB, including the climatology).
"""
import os, sys, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from multiprocessing import Pool
from common import *

HERE = os.path.dirname(os.path.abspath(__file__))
MAXGB = 48.0
BASE = ["z500", "u250", "v250", "u500", "v500", "ws250", "ws500", "d250", "d500"]
FIELDS = BASE + ["a_" + f for f in BASE]
_clim = None


def clim():
    global _clim
    if _clim is None:
        z = np.load(os.path.join(WORK, "clim.npz"))
        _clim = {k: z[k] for k in z.files}
    return _clim


def fields_at(U, V, Z, ti, when):
    """Raw and anomaly fields [121, 240] for one time of the chunk arrays ([8, 13, lon, lat])."""
    C = clim()
    h = int((when.astype("datetime64[h]") - when.astype("datetime64[D]")).astype(int)) // 6
    doy = int((when.astype("datetime64[D]") - when.astype("datetime64[Y]").astype("datetime64[D]")).astype(int)) + 1
    assert C["doy"][doy - 1] == doy and int(C["hour"][h]) == h * 6, (doy, h)
    raw, ano = {}, {}
    for lev, key in ((250, "250"), (500, "500")):
        li = LEVELS.index(lev)
        u, v = U[ti, li].T.astype(np.float64), V[ti, li].T.astype(np.float64)
        uc, vc = C["u" + key][h, doy - 1].astype(np.float64), C["v" + key][h, doy - 1].astype(np.float64)
        raw["u" + key], raw["v" + key] = u, v
        raw["ws" + key] = np.hypot(u, v)
        raw["d" + key] = divergence(u, v)
        ano["u" + key], ano["v" + key] = u - uc, v - vc
        ano["ws" + key] = np.hypot(u, v) - np.hypot(uc, vc)
        ano["d" + key] = divergence(u - uc, v - vc)
    z = Z[ti, LEVELS.index(500)].T.astype(np.float64) / G
    raw["z500"], ano["z500"] = z, z - C["z500"][h, doy - 1].astype(np.float64)
    return raw, ano


def one_frame(F, lat, lon, rot, rotate_vectors):
    la, lo = destination(lat, lon, rot + BANG, BR)
    out = np.empty((len(FIELDS), BOX_N, BOX_N), np.float32)
    for k, name in enumerate(FIELDS):
        out[k] = sample(F[name], la, lo)
    if rotate_vectors:
        th = np.radians(rot)
        for pre in ("", "a_"):
            for key in ("250", "500"):
                iu, iv = FIELDS.index(f"{pre}u{key}"), FIELDS.index(f"{pre}v{key}")
                u, v = out[iu].copy(), out[iv].copy()
                out[iu] = u * np.sin(th) + v * np.cos(th)      # along motion
                out[iv] = -u * np.cos(th) + v * np.sin(th)     # to the left of motion
    return out


def job(args):
    chunk, rows = args
    path = os.path.join(WORK, "boxes", f"{chunk:05d}.npz")
    if os.path.exists(path):
        return chunk, 0
    spent = sum(float(open(f).read()) for f in glob.glob(os.path.join(WORK, "bytes_*.txt")))
    if spent / 1e9 > MAXGB:
        raise SystemExit(f"stop: {spent / 1e9:.1f} GB already pulled (limit {MAXGB})")
    nb = 0
    A = {}
    for var, tag in (("u_component_of_wind", "u"), ("v_component_of_wind", "v"), ("geopotential", "z")):
        a, n = decode(WB2, var, f"{chunk}.0.0.0")
        A[tag], nb = a, nb + n
    keys, rot, north = [], [], []
    for r in rows.itertuples():
        when = np.datetime64(f"{r.time[:4]}-{r.time[4:6]}-{r.time[6:8]}T{r.time[8:10]}")
        step = int((when - WB2_T0) / np.timedelta64(6, "h"))
        assert step // 8 == chunk
        raw, ano = fields_at(A["u"], A["v"], A["z"], step % 8, when)
        F = dict(raw, **{"a_" + k: v for k, v in ano.items()})
        north.append(one_frame(F, r.lat, r.lon, 90.0, False))
        if np.isfinite(r.heading):
            rot.append(one_frame(F, r.lat, r.lon, r.heading, True))
        else:
            rot.append(np.full_like(north[-1], np.nan))
        keys.append(r.key)
    tmp = path + ".tmp.npz"
    np.savez_compressed(tmp, keys=np.array(keys), rot=np.stack(rot), north=np.stack(north))
    os.replace(tmp, path)
    with open(os.path.join(WORK, f"bytes_{chunk:05d}.txt"), "w") as f:
        f.write(str(nb))
    return chunk, nb


if __name__ == "__main__":
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    os.makedirs(os.path.join(WORK, "boxes"), exist_ok=True)
    D = pd.read_csv(os.path.join(HERE, "results", "storm_times.csv"), dtype={"time": str})
    t = pd.to_datetime(D.time, format="%Y%m%d%H").values.astype("datetime64[h]")
    D["chunk"] = ((t - WB2_T0) / np.timedelta64(6, "h")).astype(int) // 8
    J = sorted(((c, g) for c, g in D.groupby("chunk")), key=lambda x: -x[0])      # newest first
    print(len(J), "chunks,", len(D), "rows; expected", round(len(J) * 28.53 / 1000, 1), "GB", flush=True)
    clim()
    tot, t0 = 0, __import__("time").time()
    with Pool(nproc) as p:
        for i, (c, nb) in enumerate(p.imap_unordered(job, J)):
            tot += nb
            if i % 50 == 0:
                print(f"{i + 1}/{len(J)} chunk {c} pulled {tot / 1e9:.1f} GB {__import__('time').time() - t0:.0f}s", flush=True)
    print("done, pulled this run", round(tot / 1e9, 2), "GB")
