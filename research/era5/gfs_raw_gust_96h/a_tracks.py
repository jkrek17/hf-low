"""Stage A tracker on GFS fields: pipeline A detector (hf_history/extract.py lows_at) and
linker (hf_history/track.py link), adapted to in-memory GFS fields.

GFS 0.25 degree global grid is 721 x 1440, lat 90 to -90, lon 0 to 359.75: the same grid as the ERA5
ARCO store used by pipeline A, so every grid-point constant (Gaussian sigma 2 points, 21x21 minimum
filter, 5x5 refinement, rows 60..280 = 75N..20N) carries over unchanged. Longitude wraps as in the original.

lows_in(msl_hpa, gust_ms, ocean)  -> list of (lat, lon, msl, g500, g800, g1200), gusts in kt.
python a_tracks.py                -> link each cycle's lows over f000..f048 and write per-season track files.
"""
import os, sys, glob
import numpy as np, pandas as pd
from scipy import ndimage
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "hf_history"))
import track as T            # pipeline A linker; importing it has no side effects

WORK = os.environ.get("GFS_WORK", "/home/claude/gfs_work")
R0, R1 = 60, 281
LAT = 90 - 0.25 * np.arange(R0, R1)
LON = 0.25 * np.arange(1440)
RE = 6371.0
KT = 1 / 0.514444
LA2, LO2 = np.meshgrid(np.radians(LAT), np.radians(LON), indexing="ij")
XYZ = np.stack([np.cos(LA2) * np.cos(LO2), np.cos(LA2) * np.sin(LO2), np.sin(LA2)], -1)


def lows_in(msl_full, gust_full, ocean_full):
    """msl_full hPa, gust_full m/s, ocean_full bool; all 721x1440 (lat 90..-90, lon 0..359.75)."""
    assert msl_full.shape == (721, 1440)
    msl = np.asarray(msl_full[R0:R1], float)
    gust = gust_full[R0:R1]
    ocean = ocean_full[R0:R1]
    sm = ndimage.gaussian_filter(msl, 2, mode=["nearest", "wrap"])
    mn = ndimage.minimum_filter(sm, size=(21, 21), mode=["nearest", "wrap"])
    r, c = np.where((sm == mn) & (sm < 1010.0))
    keep = (r > 2) & (r < len(LAT) - 3)
    r, c = r[keep], c[keep]
    if len(r) == 0:
        return []
    pts = []
    for i, j in zip(r, c):
        rr = slice(max(i - 2, 0), i + 3)
        cols = np.arange(j - 2, j + 3) % 1440
        box = msl[rr][:, cols]
        a, b = np.unravel_index(np.argmin(box), box.shape)
        pts.append((box[a, b], max(i - 2, 0) + a, cols[b]))
    pts.sort()
    sel = []
    for p, i, j in pts:
        v = XYZ[i, j]
        if all(np.arccos(np.clip(v @ XYZ[a, b], -1, 1)) * RE > 400 for _, a, b in sel):
            sel.append((p, i, j))
    tree = cKDTree(np.array([XYZ[i, j] for _, i, j in sel]))
    m = ocean & (gust >= 17.0)
    gi = np.where(m)
    g = gust[gi]
    n = len(sel)
    acc = {k: np.zeros(n) for k in (500, 800, 1200)}
    if len(g):
        d, k = tree.query(XYZ[gi], k=1)
        dkm = 2 * np.arcsin(np.clip(d / 2, 0, 1)) * RE
        for rad in acc:
            w = dkm <= rad
            np.maximum.at(acc[rad], k[w], g[w])
    return [(LAT[i], LON[j], round(float(p), 1), round(acc[500][q] * KT, 1), round(acc[800][q] * KT, 1),
             round(acc[1200][q] * KT, 1)) for q, (p, i, j) in enumerate(sel)]


def stamp(cycle, lead):
    return (pd.to_datetime(cycle, format="%Y%m%d%H") + pd.Timedelta(hours=int(lead))).strftime("%Y%m%d%H")


def link_cycle(g):
    """g: lows of one cycle (columns lead, lat, lon, msl, g500, g800, g1200). Returns tracks of >= MIN_LEN fixes
    (all links, MIN_LEN not applied, with column 'track')."""
    g = g[(g.lat >= 20) & (g.lat <= 75)].copy()
    g["time"] = g["lead"]
    g = T.link(g)
    return g.drop(columns="time")


def season_of(cycle):
    y, m = int(str(cycle)[:4]), int(str(cycle)[4:6])
    return y if m >= 6 else y - 1


def main():
    files = sorted(glob.glob(os.path.join(WORK, "lows", "*.csv")))
    by = {}
    for f in files:
        by.setdefault(season_of(os.path.basename(f)[:10]), []).append(f)
    out = os.path.join(HERE, "results")
    for s, fs in sorted(by.items()):
        parts_f, parts_all, ntr = [], [], 0
        for f in fs:
            d = pd.read_csv(f, dtype={"cycle": str})
            if not len(d):
                continue
            cyc = os.path.basename(f)[:10]
            t = link_cycle(d)
            t["cycle"] = cyc
            t["valid_time"] = [stamp(cyc, l) for l in t["lead"]]
            t["basin"] = [T.basin(a, b) for a, b in zip(t.lat, t.lon)]
            n = t.groupby("track").size()
            parts_all.append(t)
            parts_f.append(t[t.track.isin(n[n >= T.MIN_LEN].index)])
        cols = ["cycle", "lead", "valid_time", "track", "lat", "lon", "msl", "g500", "g800", "g1200", "basin"]
        for name, parts, odir in ((f"gfs_tracks_{s}.csv.gz", parts_f, out),
                                  (f"gfs_tracks_all_{s}.csv.gz", parts_all, os.path.join(WORK, "tracks_all"))):
            os.makedirs(odir, exist_ok=True)
            df = pd.concat(parts, ignore_index=True)[cols]
            df.to_csv(os.path.join(odir, name), index=False)
            if odir == out:
                print(s, "cycles", len(fs), "fixes", len(df), "tracks", df.groupby(["cycle", "track"]).ngroups,
                      "MB", round(os.path.getsize(os.path.join(odir, name)) / 1e6, 2))


def gust_table():
    """Concatenate the per-cycle basin gust files into results/gfs_basin_gust.csv.gz (columns add valid_time, season)."""
    d = pd.concat([pd.read_csv(f, dtype={"cycle": str}) for f in sorted(glob.glob(os.path.join(WORK, "gust", "*.csv")))],
                  ignore_index=True)
    d["valid_time"] = [stamp(c, l) for c, l in zip(d.cycle, d.lead)]
    d["season"] = [season_of(c) for c in d.cycle]
    d.to_csv(os.path.join(HERE, "results", "gfs_basin_gust.csv.gz"), index=False)
    print("gust rows", len(d), "MB", round(os.path.getsize(os.path.join(HERE, "results", "gfs_basin_gust.csv.gz")) / 1e6, 2))


if __name__ == "__main__":
    main()
    gust_table()
