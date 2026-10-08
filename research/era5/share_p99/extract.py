"""RA-28 extraction: owned-cell percentiles of the ERA5 instantaneous 10 m gust at pipeline A's 00/12 UTC fixes.

ERA5 PROXY, pipeline A (research/era5/hf_history). For every fix listed in a fixes CSV (columns track,time,lat,lon,g800 as in
intensity/results/fixes_2004.csv.gz) it fetches the global ERA5 instantaneous 10 m gust and MSLP at that hour from ARCO-ERA5,
re-detects lows with pipeline A's rule (hf_structure/extract.py detect, which reproduced the catalog g800 for 99.98% of fixes in RA-27)
and takes the cells the matched low owns: ocean, nearer to this low than to any other, within 1200 km of some low and within 800 km of
this one. Per fix it writes the maximum over those cells (g800_re, must return the catalog g800), and the 95th, 98th and 99th
percentiles of the same cells (unweighted by area, exactly as RA-27's S7 did for sustained wind).

One CSV per time under $ERA5_WORK/share_p99/<tag>/, resumable. Raw fields are never written.

usage: python3 -I extract.py REPO FIXES_CSV TAG [NPROC]
"""
import os, sys, time
import numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.spatial import cKDTree

root, fixes_csv, tag = sys.argv[1], sys.argv[2], sys.argv[3]
nproc = int(sys.argv[4]) if len(sys.argv) > 4 else 4
HERE = os.path.join(root, "research/era5/share_p99")
os.environ["ERA5_WORK"] = os.environ.get("ERA5_WORK", os.path.join(HERE, "work"))
WORK = os.path.join(os.environ["ERA5_WORK"], "share_p99", tag)
os.makedirs(WORK, exist_ok=True)
sys.path.insert(0, os.path.join(root, "research/era5/hf_structure"))
import extract as hs   # noqa: E402  (hf_structure/extract.py; uses ERA5_WORK for its land-sea mask cache)

KT, RE, R0, R1 = hs.KT, hs.RE, hs.R0, hs.R1
ARCO, ARCO_T0, XYZ, OCEAN = hs.ARCO, hs.ARCO_T0, hs.XYZ, hs.OCEAN
BAND = XYZ[R0:R1]
S = pd.read_csv(fixes_csv, dtype={"time": str})


def one_time(stamp):
    path = os.path.join(WORK, f"{stamp}.csv")
    if os.path.exists(path):
        return stamp, 0
    t = np.datetime64(f"{stamp[:4]}-{stamp[4:6]}-{stamp[6:8]}T{stamp[8:10]}")
    ti = int((t - ARCO_T0) / np.timedelta64(1, "h"))
    msl, b1 = hs.fetch(f"{ARCO}/mean_sea_level_pressure/{ti}.0.0")
    gust, b2 = hs.fetch(f"{ARCO}/instantaneous_10m_wind_gust/{ti}.0.0")
    msl = msl / 100.0
    gust = gust * KT
    lows = hs.detect(msl)
    cen = np.array([XYZ[i, j] for _, i, j in lows])
    d, owner = cKDTree(cen).query(XYZ[R0:R1].reshape(-1, 3), k=1)
    dnear = (2 * np.arcsin(np.clip(d / 2, 0, 1)) * RE).reshape(R1 - R0, 1440)
    owner = owner.reshape(R1 - R0, 1440)
    rows = []
    for _, fx in S[S.time == stamp].iterrows():
        vx = XYZ[int(round((90 - fx.lat) / 0.25)), int(round(fx.lon / 0.25)) % 1440]
        dk = np.arccos(np.clip(cen @ vx, -1, 1)) * RE
        k = int(np.argmin(dk))
        _, ci, cj = lows[k]
        dist = np.arccos(np.clip(BAND @ XYZ[ci, cj], -1, 1)) * RE
        r8 = (owner == k) & (dnear <= 1200) & OCEAN[R0:R1] & (dist <= 800)
        row = dict(track=int(fx.track), time=stamp, match_km=round(float(dk[k]), 1), n_cells=int(r8.sum()), g800_cat=float(fx.g800))
        if r8.any():
            g = gust[R0:R1][r8]
            row.update(g800_re=float(g.max()), p95=float(np.percentile(g, 95)), p98=float(np.percentile(g, 98)), p99=float(np.percentile(g, 99)))
        rows.append(row)
    tmp = path + ".tmp"
    pd.DataFrame(rows).to_csv(tmp, index=False)
    os.replace(tmp, path)
    with open(os.path.join(WORK, "bytes.log"), "a") as f:
        f.write(f"{stamp},{b1 + b2}\n")
    return stamp, b1 + b2


if __name__ == "__main__":
    stamps = sorted(S.time.unique())
    todo = [s for s in stamps if not os.path.exists(os.path.join(WORK, f"{s}.csv"))]
    print(f"{len(stamps)} times, {len(todo)} to do", flush=True)
    t0, done, nb = time.time(), 0, 0
    with Pool(nproc) as p:
        for st, b in p.imap_unordered(one_time, todo):
            done += 1
            nb += b
            if done % 50 == 0:
                print(f"{done}/{len(todo)} {time.time() - t0:.0f}s {nb / 1e9:.2f} GB", flush=True)
    print(f"done {time.time() - t0:.0f}s {nb / 1e9:.2f} GB", flush=True)
