"""Earth-relative basin composites of HF lows: MSLP, Z500, 250 hPa wind speed (WeatherBench2 ERA5, 1.5 degrees; a pipeline A proxy).

Descriptive, no test. Population: the 1,223 HF lows of ../hf_env_composites/results/storm_times.csv (onset = first 00/12 UTC fix with
gust index >= 71.7 kt; Oct-Apr onsets, 2004-05 to 2022-23; 664 Atlantic, 559 Pacific). Two times per low: onset and onset minus 48 h
(the calendar time, whether or not the track existed then). Anomaly = ERA5 minus the WeatherBench2 1990-2019 climatology at the same
day of year and hour; wind-speed anomaly = |V| - |V_clim|.

usage: ERA5_WORK=DIR extract.py [NPROC]   (needs DIR/clim.npz from ../hf_env_composites/clim.py and DIR/clim_msl.npz from clim_msl.py)
Output DIR/basin_sums.npz: per group (basin, anchor) the count and the sums of [msl, z500, ws250, a_msl, a_z500, a_ws250]. Resumable
(checkpoint every 100 chunks); stops if the pull passes MAXGB (45 GB). No raw field is kept.
"""
import os, sys, glob
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "hf_env_composites"))
import numpy as np, pandas as pd
from multiprocessing import Pool
from common import decode, WB2, WB2_T0, LEVELS, WORK, G

MAXGB = 45.0
NAMES = ["msl", "z500", "ws250", "a_msl", "a_z500", "a_ws250"]
L250, L500 = LEVELS.index(250), LEVELS.index(500)
_C = {}


def events():
    S = pd.read_csv(os.path.join(HERE, "..", "hf_env_composites", "results", "storm_times.csv"), dtype={"time": str})
    on = S[S.anchor == "onset"]
    rows = []
    for r in on.itertuples():
        t = pd.to_datetime(r.time, format="%Y%m%d%H")
        for anchor, tt in (("onset", t), ("m48", t - pd.Timedelta(hours=48))):
            k = int((np.datetime64(tt) - WB2_T0) / np.timedelta64(6, "h"))
            rows.append(dict(track=r.track, basin=r.basin, anchor=anchor, when=np.datetime64(tt, "h"), chunk=k // 8, step=k % 8))
    return pd.DataFrame(rows)


def clim():
    if not _C:
        z = np.load(os.path.join(WORK, "clim.npz"))
        _C.update({k: z[k] for k in z.files})
        _C["msl"] = np.load(os.path.join(WORK, "clim_msl.npz"))["msl"]
    return _C


def job(args):
    chunk, rows = args
    C = clim()
    U, nb1 = decode(WB2, "u_component_of_wind", f"{chunk}.0.0.0")
    V, nb2 = decode(WB2, "v_component_of_wind", f"{chunk}.0.0.0")
    Z, nb3 = decode(WB2, "geopotential", f"{chunk}.0.0.0")
    M, nb4 = decode(WB2, "mean_sea_level_pressure", f"{chunk}.0.0")
    out = {}
    for r in rows.itertuples():
        w = np.datetime64(pd.Timestamp(r.when), "h")
        h = int((w - w.astype("datetime64[D]")).astype(int)) // 6
        doy = int((w.astype("datetime64[D]") - w.astype("datetime64[Y]").astype("datetime64[D]")).astype(int)) + 1
        assert C["doy"][doy - 1] == doy and int(C["hour"][h]) == h * 6
        u, v = U[r.step, L250].T.astype(np.float64), V[r.step, L250].T.astype(np.float64)
        uc, vc = C["u250"][h, doy - 1].astype(np.float64), C["v250"][h, doy - 1].astype(np.float64)
        z = Z[r.step, L500].T.astype(np.float64) / G
        m = M[r.step].T.astype(np.float64) / 100.0
        f = np.stack([m, z, np.hypot(u, v), m - C["msl"][h, doy - 1], z - C["z500"][h, doy - 1].astype(np.float64), np.hypot(u, v) - np.hypot(uc, vc)])
        key = (r.basin, r.anchor)
        n, s = out.get(key, (0, 0.0))
        out[key] = (n + 1, s + f)
    return chunk, out, nb1 + nb2 + nb3 + nb4


if __name__ == "__main__":
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    ck = os.path.join(WORK, "basin_ckpt.npz")
    E = events()
    chunks = sorted(((c, g) for c, g in E.groupby("chunk")), key=lambda x: -x[0])
    acc, done, spent = {}, set(), 0
    if os.path.exists(ck):
        z = np.load(ck, allow_pickle=True)
        done, spent = set(z["done"].tolist()), float(z["spent"])
        acc = {tuple(k.split("|")): (int(z["n_" + k]), z["s_" + k]) for k in z["keys"]}
    todo = [x for x in chunks if x[0] not in done]
    print(len(chunks), "chunks,", len(E), "event-times; todo", len(todo), "expected", round(len(chunks) * 29.08 / 1000, 1), "GB", flush=True)

    def save():
        keys = ["|".join(k) for k in acc]
        np.savez(ck, done=np.array(sorted(done)), spent=spent, keys=np.array(keys), **{"n_" + k: acc[tuple(k.split("|"))][0] for k in keys},
                 **{"s_" + k: acc[tuple(k.split("|"))][1] for k in keys})
    with Pool(nproc) as p:
        for i, (c, out, nb) in enumerate(p.imap_unordered(job, todo)):
            spent += nb
            done.add(c)
            for k, (n, s) in out.items():
                n0, s0 = acc.get(k, (0, 0.0))
                acc[k] = (n0 + n, s0 + s)
            if spent / 1e9 > MAXGB:
                save(); raise SystemExit(f"stop: {spent / 1e9:.1f} GB (limit {MAXGB})")
            if i % 100 == 99:
                save(); print(f"{i + 1}/{len(todo)} chunks, {spent / 1e9:.1f} GB", flush=True)
    save()
    np.savez(os.path.join(WORK, "basin_sums.npz"), **{"n_" + "|".join(k): v[0] for k, v in acc.items()}, **{"s_" + "|".join(k): v[1] for k, v in acc.items()})
    print("done", round(spent / 1e9, 2), "GB;", {k: v[0] for k, v in acc.items()})
