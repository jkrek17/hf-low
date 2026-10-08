"""RA-20 Arm B extraction: mesoscale surface features at native 0.25 degrees, per stratum fix, at the times in
results/times_armB.csv. ERA5 proxy, pipeline A lows re-detected with hf_structure's code (which reproduced the catalog g800 for
99.98% of fixes). Features are the 11 defined in PREREGISTRATION.md. One CSV per time under $ERA5_WORK/explosive_onset/armB/,
resumable; raw fields are never written. Needs: gust, mslp, u10, v10 per time (about 12 MB per time compressed).

Run: python3 research/era5/explosive_onset/extract_armB.py <repo_root> [n_proc] [--dry]
"""
import sys, os, time, numpy as np, pandas as pd
from multiprocessing import Pool
from scipy import ndimage
from scipy.spatial import cKDTree
root = sys.argv[1]
HERE = os.path.join(root, "research/era5/explosive_onset")
WORK = os.path.join(os.environ.get("ERA5_WORK", os.path.join(HERE, "work")), "armB")
os.makedirs(WORK, exist_ok=True)
os.environ["ERA5_WORK"] = os.environ.get("ERA5_WORK", os.path.join(HERE, "work"))
sys.path.insert(0, os.path.join(root, "research/era5/hf_structure"))
import extract as hs                      # pipeline A detection, ownership, land-sea mask (module-level setup fetches the mask only)
KT, RE, R0, R1 = hs.KT, hs.RE, hs.R0, hs.R1
ARCO, ARCO_T0 = hs.ARCO, hs.ARCO_T0
LAT, LON, XYZ, OCEAN, CELL = hs.LAT, hs.LON, hs.XYZ, hs.OCEAN, hs.CELL
DY = np.radians(0.25) * RE                            # km per 0.25 deg of latitude
VARS = {"gust": "instantaneous_10m_wind_gust", "msl": "mean_sea_level_pressure", "u": "10m_u_component_of_wind", "v": "10m_v_component_of_wind"}
F = pd.read_csv(os.path.join(root, "research/era5/intensity/results/fixes_2004.csv.gz"))
F = F[(F.g800 < 55) & (~F.hf_now)]
T = pd.read_csv(os.path.join(HERE, "results/times_armB.csv"))
BAND = XYZ[R0:R1]; BCELL = np.broadcast_to(CELL[R0:R1, None], (R1 - R0, 1440))


def quad_curv(msl, ci, cj, lat0):
    """Laplacian of MSLP (hPa per (100 km)^2) from a quadratic fit to the points within 150 km of grid point (ci, cj)."""
    n = 7
    ii = np.arange(ci - n, ci + n + 1); jj = np.arange(cj - n, cj + n + 1) % 1440
    sub = msl[np.clip(ii, 0, 720)][:, jj]
    y = (ci - ii)[:, None] * DY * np.ones((1, len(jj))); x = (np.arange(-n, n + 1) * DY * np.cos(np.radians(lat0)))[None, :] * np.ones((len(ii), 1))
    m = np.hypot(x, y) <= 150
    A = np.stack([np.ones(m.sum()), x[m], y[m], x[m] ** 2, x[m] * y[m], y[m] ** 2], 1)
    c = np.linalg.lstsq(A, sub[m], rcond=None)[0]
    return float((2 * c[3] + 2 * c[5]) * 1e4)


def one_time(stamp):
    path = os.path.join(WORK, f"{stamp}.csv")
    if os.path.exists(path):
        return stamp, 0
    ti = int((np.datetime64(f"{stamp[:4]}-{stamp[4:6]}-{stamp[6:8]}T{stamp[8:10]}") - ARCO_T0) / np.timedelta64(1, "h"))
    fld, nb = {}, 0
    for k, v in VARS.items():
        fld[k], b = hs.fetch(f"{ARCO}/{v}/{ti}.0.0"); nb += b
    msl = fld["msl"] / 100.0; gust = fld["gust"] * KT
    u, v = fld["u"], fld["v"]; ws = np.hypot(u, v) * KT
    smooth = ndimage.gaussian_filter(msl, 6, mode=["nearest", "wrap"])          # sigma 1.5 deg
    gy, gx = np.gradient(msl, axis=(0, 1))                                      # hPa per grid step
    dxkm = (DY * np.cos(np.radians(LAT)))[:, None]
    grad = np.hypot(-gy / DY, gx / dxkm) * 100.0                                # hPa per 100 km
    wy, wx = np.gradient(ws, axis=(0, 1))
    shear = np.hypot(-wy / DY, wx / dxkm) * 100.0                               # kt per 100 km
    dudx = np.gradient(u, axis=1) / (dxkm * 1e3); dvdy = -np.gradient(v, axis=0) / (DY * 1e3)
    div = (dudx + dvdy) * 1e5                                                   # 1e-5 per second
    lows = hs.detect(msl)
    cen = np.array([XYZ[i, j] for _, i, j in lows])
    d, owner = cKDTree(cen).query(XYZ[R0:R1].reshape(-1, 3), k=1)
    dnear = (2 * np.arcsin(np.clip(d / 2, 0, 1)) * RE).reshape(R1 - R0, 1440); owner = owner.reshape(R1 - R0, 1440)
    rows = []
    for _, fx in F[F.time == int(stamp)].iterrows():
        vx = XYZ[int(round((90 - fx.lat) / 0.25)), int(round(fx.lon / 0.25)) % 1440]
        dk = np.arccos(np.clip(cen @ vx, -1, 1)) * RE; k = int(np.argmin(dk))
        _, ci, cj = lows[k]
        dist = np.arccos(np.clip(BAND @ XYZ[ci, cj], -1, 1)) * RE
        own = (owner == k) & (dnear <= 1200) & OCEAN[R0:R1]
        r8 = own & (dist <= 800); G = gust[R0:R1]; W = ws[R0:R1]
        row = dict(track=int(fx.track), time=int(fx.time), match_km=round(float(dk[k]), 1), g800_re=round(float(np.where(r8, G, 0).max()), 2), g800=fx.g800)
        if not r8.any():
            rows.append(row); continue
        g = np.where(r8, G, -1.0); im = np.unravel_index(np.argmax(g), g.shape); gmax = float(G[im])
        row["A45"] = float(np.log1p(BCELL[r8 & (G >= 45)].sum()))
        row["PK"] = gmax - float(np.percentile(G[r8], 99))
        row["RG"] = float(dist[im])
        row["SHARP"] = float(msl[ci + R0, cj] - smooth[ci + R0, cj])
        row["CURV"] = quad_curv(msl, ci + R0, cj, LAT[ci + R0])
        c4 = dist <= 400
        row["GRAD"] = float(np.percentile(grad[R0:R1][c4], 90))
        row["WS99"] = float(np.percentile(W[r8], 99))
        row["GF"] = gmax / max(float(W[im]), 1.0)
        row["A34"] = float(np.log1p(BCELL[r8 & (W >= 34)].sum()))
        row["SHR"] = float(shear[R0:R1][r8].max())
        row["CONV"] = float(-div[R0:R1][c4].min())
        rows.append(row)
    tmp = path + ".tmp"; pd.DataFrame(rows).to_csv(tmp, index=False); os.replace(tmp, path)
    with open(os.path.join(WORK, "bytes.log"), "a") as f: f.write(f"{stamp},{nb}\n")
    return stamp, nb


if __name__ == "__main__":
    stamps = [str(t) for t in T.time]
    n = int(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2].isdigit() else 8
    if "--dry" in sys.argv:
        print(len(stamps), "times,", int(F.time.isin(T.time).sum()), "stratum fixes"); sys.exit()
    t0, tot = time.time(), 0
    with Pool(n) as p:
        for i, (s, nb) in enumerate(p.imap_unordered(one_time, stamps, chunksize=2)):
            tot += nb
            if i % 100 == 0: print(f"{i + 1}/{len(stamps)} {tot / 1e9:.2f} GB {time.time() - t0:.0f}s", flush=True)
    print(f"done {len(stamps)} times, {tot / 1e9:.2f} GB this run")
