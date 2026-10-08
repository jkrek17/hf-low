"""Test 3 stage B extraction: gust + MSLP at the six new life-cycle lags. ERA5 proxy, pipeline A low detection and ownership (helpers imported from hf_vs_storm/extract_storm.py).
Position = interpolated track position snapped to the nearest re-detected low within 150 km, else missing. usage: HFVS_SUB=hf_lifecycle ERA5_WORK=... extract_lc.py NPROC [--dry|N]"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "hf_vs_storm"))
os.environ.setdefault("HFVS_SUB", "hf_lifecycle")
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
from multiprocessing import Pool
import extract_storm as E
WORK = E.WORK; SNAP = 150.0


def one_time(job):
    stamp, fixes = job
    spath = os.path.join(WORK, "stats", f"{stamp}.csv")
    if os.path.exists(spath):
        return stamp, 0
    t = np.datetime64(f"{stamp[:4]}-{stamp[4:6]}-{stamp[6:8]}T{stamp[8:10]}")
    ti = int((t - E.ARCO_T0) / np.timedelta64(1, "h")); nbytes = 0; F = {}
    for name, var in (("msl", "mean_sea_level_pressure"), ("gust", "instantaneous_10m_wind_gust")):
        F[name], nb = E.fetch(f"{E.ARCO}/{var}/{ti}.0.0"); nbytes += nb
    msl = F["msl"] / 100.0; gust = F["gust"] * E.KT
    lows = E.detect(msl); cen = np.array([E.XYZ[i, j] for _, i, j in lows])
    d, owner = cKDTree(cen).query(E.XYZ[E.R0:E.R1].reshape(-1, 3), k=1)
    dnear = (2 * np.arcsin(np.clip(d / 2, 0, 1)) * E.RE).reshape(E.R1 - E.R0, 1440); owner = owner.reshape(E.R1 - E.R0, 1440)
    rows, boxes = [], {}
    for fx in fixes:
        v = E.XYZ[int(round((90 - fx["lat"]) / 0.25)), int(round(fx["lon"] / 0.25)) % 1440]
        dk = np.arccos(np.clip(cen @ v, -1, 1)) * E.RE; k = int(np.argmin(dk))
        base = dict(time=stamp, track=fx["track"], anchor=fx["anchor"], basin=fx["basin"], lat=fx["lat"], lon=fx["lon"], heading=fx["heading"])
        if dk[k] > SNAP:
            rows.append(dict(base, snapped=False)); continue
        own = np.zeros((721, 1440), bool); own[E.R0:E.R1] = (owner == k) & (dnear <= 1200) & E.OCEAN[E.R0:E.R1]
        _, ci, cj = lows[k]; clat, clon = E.LAT[ci], E.LON[cj]
        G, O = gust[E.R0:E.R1], own[E.R0:E.R1]
        dist = np.arccos(np.clip(E.BAND @ E.XYZ[ci, cj], -1, 1)) * E.RE
        r1200 = O & (dist <= 1200); r800 = O & (dist <= 800); g = np.where(r1200, G, -1.0)
        im = np.unravel_index(np.argmax(g), g.shape); m48 = r1200 & (G >= 48.0)
        hd = fx["heading"]
        b = {}
        for name, fld in (("gust", gust), ("msl", msl)):
            b[name] = E.box(fld, clat, clon, hd).astype(np.float16) if np.isfinite(hd) else np.full((E.BOX_N, E.BOX_N), np.nan, np.float16)
        b["own"] = E.box(own.astype(np.float32), clat, clon, hd, 0).astype(np.uint8) if np.isfinite(hd) else np.zeros((E.BOX_N, E.BOX_N), np.uint8)
        nb_ = E.box(msl, clat, clon, 90.0)
        pc = float(np.nanmin(np.where(E.BR <= 100, nb_, np.nan))); ring = float(np.nanmean(nb_[(E.BR >= 450) & (E.BR <= 550)]))
        rows.append(dict(base, snapped=True, match_km=round(float(dk[k]), 1), g800=round(float(np.where(r800, G, 0).max()), 1), gmax=round(float(G[im]), 1), gmax_r=round(float(dist[im]), 0),
                         a_g48=round(float(E.BCELL[m48].sum()), 0), g48_rmax=round(float(dist[m48].max()), 0) if m48.any() else np.nan,
                         own_ocean_frac=round(float(r1200.sum() / max((dist <= 1200).sum(), 1)), 3), pc=round(pc, 1), msl_ring500=round(ring, 2), msl_grad=round((ring - pc) / 5.0, 3)))
        for kk, vv in b.items():
            boxes[f"{fx['track']}_{fx['anchor']}_{kk}"] = vv
    np.savez_compressed(os.path.join(WORK, "boxes", f"{stamp}.npz"), **boxes)
    pd.DataFrame(rows).to_csv(spath + ".tmp", index=False); os.replace(spath + ".tmp", spath)
    with open(os.path.join(WORK, "bytes.log"), "a") as f:
        f.write(f"{stamp},{nbytes}\n")
    return stamp, nbytes


def jobs():
    T = pd.read_csv(os.environ.get("HFVS_TIMES", os.path.join(HERE, "results", "times_lc.csv")))
    return [(str(t), g[["track", "anchor", "basin", "lat", "lon", "heading"]].to_dict("records")) for t, g in T.groupby("time")]


if __name__ == "__main__":
    J = jobs()
    if "--dry" in sys.argv:
        print(len(J), "times", sum(len(j[1]) for j in J), "anchors"); sys.exit()
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    lim = int(sys.argv[2]) if len(sys.argv) > 2 else None
    if lim: J = J[:lim]
    t0, tot = time.time(), 0
    with Pool(n) as p:
        for i, (s, nb) in enumerate(p.imap_unordered(one_time, J, chunksize=2)):
            tot += nb
            if i % 100 == 0: print(f"{i + 1}/{len(J)} {tot / 1e9:.2f} GB {time.time() - t0:.0f}s", flush=True)
    print(f"done {len(J)} times, {tot / 1e9:.2f} GB this run")
