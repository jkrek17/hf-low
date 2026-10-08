"""RA-27 extraction: 10 m sustained wind index in two analysis systems at the same centres.

ERA5 (ARCO 0.25 degree hourly) and the IFS HRES t0 analysis (WeatherBench2 hres_t0, 0.25 degree, 6-hourly) are read at each time in
results/sample_fixes_<floor>.csv. Centres and ownership come from ERA5 lows re-detected with hf_structure's code (pipeline A's rule,
which reproduced the catalog g800 for 99.98% of fixes), so only the wind field differs between the two systems. Per fix: the maximum
10 m wind speed (kt) over ocean cells within 800 km that the fix owns, and the 99th percentile of the same cells, in each system.
The 150 validation times (fixed by seed, inside the same run) also read the ERA5 gust and require pipeline A's g800 back. One CSV per time under $ERA5_WORK/second_analysis/,
resumable; raw fields are never written.

Run: python3 -I research/era5/second_analysis/extract.py <repo_root> <floor_kt> [n_proc]
"""
import sys, os, time, urllib.request, numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.spatial import cKDTree
import numcodecs
root, floor = sys.argv[1], int(sys.argv[2])
nproc = int(sys.argv[3]) if len(sys.argv) > 3 and not sys.argv[3].startswith("--") else 8
VALIDATE = False
HERE = os.path.join(root, "research/era5/second_analysis")
WORK = os.path.join(os.environ.get("ERA5_WORK", os.path.join(HERE, "work")), f"floor{floor}")
os.makedirs(WORK, exist_ok=True)
os.environ["ERA5_WORK"] = os.environ.get("ERA5_WORK", os.path.join(HERE, "work"))
sys.path.insert(0, os.path.join(root, "research/era5/hf_structure"))
import extract as hs
KT, RE, R0, R1 = hs.KT, hs.RE, hs.R0, hs.R1
ARCO, ARCO_T0, XYZ, OCEAN = hs.ARCO, hs.ARCO_T0, hs.XYZ, hs.OCEAN
HRES = "https://storage.googleapis.com/weatherbench2/datasets/hres_t0/2016-2022-6h-1440x721.zarr"
HRES_T0 = np.datetime64("2016-01-01T00")
_codec = numcodecs.Blosc()
BAND = XYZ[R0:R1]
S = pd.read_csv(os.path.join(HERE, f"results/sample_fixes_{floor}.csv"))
rng = np.random.default_rng(20261008)
VAL_TIMES = set(int(x) for x in rng.choice(np.array(sorted(S.time.unique())), 150, replace=False))   # 150 times also read for the gust, to check ownership


def fetch_h(url):
    for k in range(6):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                raw = r.read()
            return np.frombuffer(_codec.decode(raw), "<f4").reshape(721, 1440)[::-1], len(raw)   # HRES latitude runs south to north
        except Exception:
            if k == 5:
                raise
            time.sleep(2 ** k)


def one_time(stamp):
    path = os.path.join(WORK, f"{stamp}.csv")
    if os.path.exists(path):
        return stamp, 0
    t = np.datetime64(f"{stamp[:4]}-{stamp[4:6]}-{stamp[6:8]}T{stamp[8:10]}")
    ti = int((t - ARCO_T0) / np.timedelta64(1, "h")); th = int((t - HRES_T0) / np.timedelta64(1, "h")) // 6
    nb = 0
    eu, b = hs.fetch(f"{ARCO}/10m_u_component_of_wind/{ti}.0.0"); nb += b
    ev, b = hs.fetch(f"{ARCO}/10m_v_component_of_wind/{ti}.0.0"); nb += b
    msl, b = hs.fetch(f"{ARCO}/mean_sea_level_pressure/{ti}.0.0"); nb += b
    hu, b = fetch_h(f"{HRES}/10m_u_component_of_wind/{th}.0.0"); nb += b
    hv, b = fetch_h(f"{HRES}/10m_v_component_of_wind/{th}.0.0"); nb += b
    VAL = int(stamp) in VAL_TIMES
    if VAL:
        gust, b = hs.fetch(f"{ARCO}/instantaneous_10m_wind_gust/{ti}.0.0"); nb += b; gust = gust * KT
    wE = np.hypot(eu, ev) * KT; wH = np.hypot(hu, hv) * KT
    lows = hs.detect(msl / 100.0)
    cen = np.array([XYZ[i, j] for _, i, j in lows])
    d, owner = cKDTree(cen).query(XYZ[R0:R1].reshape(-1, 3), k=1)
    dnear = (2 * np.arcsin(np.clip(d / 2, 0, 1)) * RE).reshape(R1 - R0, 1440); owner = owner.reshape(R1 - R0, 1440)
    rows = []
    for _, fx in S[S.time == int(stamp)].iterrows():
        vx = XYZ[int(round((90 - fx.lat) / 0.25)), int(round(fx.lon / 0.25)) % 1440]
        dk = np.arccos(np.clip(cen @ vx, -1, 1)) * RE; k = int(np.argmin(dk))
        _, ci, cj = lows[k]
        dist = np.arccos(np.clip(BAND @ XYZ[ci, cj], -1, 1)) * RE
        r8 = (owner == k) & (dnear <= 1200) & OCEAN[R0:R1] & (dist <= 800)
        row = dict(track=int(fx.track), time=int(fx.time), match_km=round(float(dk[k]), 1), n_cells=int(r8.sum()))
        if r8.any():
            e, h = wE[R0:R1][r8], wH[R0:R1][r8]
            row.update(wE_max=float(e.max()), wH_max=float(h.max()), wE_p99=float(np.percentile(e, 99)), wH_p99=float(np.percentile(h, 99)))
            if VAL:
                row["g800_re"] = float(gust[R0:R1][r8].max()); row["g800"] = float(fx.g800)
        rows.append(row)
    tmp = path + ".tmp"; pd.DataFrame(rows).to_csv(tmp, index=False); os.replace(tmp, path)
    with open(os.path.join(WORK, "bytes.log"), "a") as f: f.write(f"{stamp},{nb}\n")
    return stamp, nb


if __name__ == "__main__":
    stamps = [str(t) for t in sorted(S.time.unique())]
    todo = [s for s in stamps if not os.path.exists(os.path.join(WORK, f"{s}.csv"))]
    print(f"{len(stamps)} times, {len(todo)} to do, floor {floor}", flush=True)
    t0 = time.time(); done = 0
    with Pool(nproc) as p:
        for st, nb in p.imap_unordered(one_time, todo):
            done += 1
            if done % 100 == 0: print(f"{done}/{len(todo)} {time.time()-t0:.0f}s", flush=True)
    print("done", flush=True)
