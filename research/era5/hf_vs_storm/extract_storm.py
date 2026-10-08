"""Storm-relative ERA5 fields at HF onset / HF peak / storm-force-only peak times (hf_vs_storm).

Helpers copied from hf_structure/extract.py (low detection, ownership, boxes). Original docstring:

Storm-relative ERA5 wind fields at pipeline A's hurricane-force-strength times.

ERA5 proxy, pipeline A (research/era5/hf_history). For every event-track point
in era5_hf_catalog_tracks.csv with gust index g800 >= 71.7 kt in seasons
2004-05 on, fetch the global ERA5 instantaneous 10 m gust and MSLP at that hour
from ARCO-ERA5 and summarise the gust field the low owns:

  - Lows are re-detected exactly as hf_history/extract.py does (smoothed MSLP
    minima below 1010 hPa over 20-75N, 400 km de-duplication), and a grid point
    belongs to the low it is nearest to, within 1200 km, ocean only. This is
    pipeline A's own ownership rule, so g800 recomputed here must reproduce the
    catalog value; the stats file carries both so the check is visible.
  - Storm motion is the bearing from the track point 6 h before to the one 6 h
    after (one-sided at a track end).

For a storm-stratified subsample of times before 2023-01-10 it also fetches the
6-hourly ERA5 10 m wind speed from WeatherBench2's 0.25 degree copy (ARCO-ERA5
has only u and v, which cost twice as much), to relate sustained wind to gust.

Outputs, under $ERA5_WORK/hf_structure (ignored, resumable, never committed):
  stats/<time>.csv   one row per HF fix: native-grid statistics (see STATS)
  boxes/<time>.npz   per fix, 121 x 121 boxes at 25 km: motion-relative gust,
                     MSLP, owned-ocean mask, optional wind speed; north-up gust
                     and mask. Used by composite.py.
  bytes.log          bytes streamed per time, for the pull size

usage: extract.py [NPROC]
"""
import os, sys, time, urllib.request
import numpy as np, pandas as pd, numcodecs
from scipy import ndimage
from scipy.spatial import cKDTree
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "hf_history", "results")
OWN = os.path.join(HERE, "results")
WORK = os.path.join(os.environ.get("ERA5_WORK", os.path.join(HERE, "work")), os.environ.get("HFVS_SUB", "hf_vs_storm"))

ARCO = "https://storage.googleapis.com/gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3"
WB2 = ("https://storage.googleapis.com/weatherbench2/datasets/era5/"
       "1959-2023_01_10-wb13-6h-1440x721_with_derived_variables.zarr")
ARCO_T0 = np.datetime64("1900-01-01T00")
WB2_T0 = np.datetime64("1959-01-01T00")
WB2_END = np.datetime64("2023-01-10T00")
HF = 71.7                      # pipeline A threshold (kt), threshold.json
SEASON0 = 2004                 # RECORD_START
N_WIND = 1500                  # times in the wind-speed subsample
KT = 1 / 0.514444
RE = 6371.0
R0, R1 = 60, 281               # rows 75N .. 20N, pipeline A's detection band
LAT = 90 - 0.25 * np.arange(721)
LON = 0.25 * np.arange(1440)
CELL = (np.radians(0.25) * RE) ** 2 * np.cos(np.radians(LAT))     # km2 per cell, by row
BOX_N, BOX_D = 121, 25.0       # box points per side, km between points (+-1500 km)
_codec = numcodecs.Blosc()

LA2, LO2 = np.meshgrid(np.radians(LAT), np.radians(LON), indexing="ij")
XYZ = np.stack([np.cos(LA2) * np.cos(LO2), np.cos(LA2) * np.sin(LO2), np.sin(LA2)], -1)


def fetch(url):
    for k in range(6):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                raw = r.read()
            return np.frombuffer(_codec.decode(raw), "<f4").reshape(721, 1440), len(raw)
        except Exception:
            if k == 5:
                raise
            time.sleep(2 ** k)


def land_sea():
    path = os.path.join(WORK, "lsm.npy")
    if not os.path.exists(path):
        a, _ = fetch(f"{ARCO}/land_sea_mask/876576.0.0")         # 2000-01-01 00Z, as pipeline A
        np.save(path, a)
    return np.load(path)


def detect(msl):
    """Pipeline A's low detection (hf_history/extract.py lows_at), on full-grid rows."""
    m = msl[R0:R1]
    sm = ndimage.gaussian_filter(m, 2, mode=["nearest", "wrap"])
    mn = ndimage.minimum_filter(sm, size=(21, 21), mode=["nearest", "wrap"])
    r, c = np.where((sm == mn) & (sm < 1010.0))
    keep = (r > 2) & (r < R1 - R0 - 3)
    pts = []
    for i, j in zip(r[keep], c[keep]):
        rr = slice(max(i - 2, 0), i + 3)
        cols = np.arange(j - 2, j + 3) % 1440
        box = m[rr][:, cols]
        a, b = np.unravel_index(np.argmin(box), box.shape)
        pts.append((box[a, b], R0 + max(i - 2, 0) + a, cols[b]))
    pts.sort()
    sel = []
    for p, i, j in pts:
        v = XYZ[i, j]
        if all(np.arccos(np.clip(v @ XYZ[a, b], -1, 1)) * RE > 400 for _, a, b in sel):
            sel.append((p, i, j))
    return sel


def bearing(lat1, lon1, lat2, lon2):
    p1, p2, dl = np.radians(lat1), np.radians(lat2), np.radians(lon2 - lon1)
    x = np.sin(dl) * np.cos(p2)
    y = np.cos(p1) * np.sin(p2) - np.sin(p1) * np.cos(p2) * np.cos(dl)
    return np.degrees(np.arctan2(x, y)) % 360


def destination(lat, lon, brg, dkm):
    p, l, b, d = np.radians(lat), np.radians(lon), np.radians(brg), dkm / RE
    p2 = np.arcsin(np.sin(p) * np.cos(d) + np.cos(p) * np.sin(d) * np.cos(b))
    l2 = l + np.arctan2(np.sin(b) * np.sin(d) * np.cos(p), np.cos(d) - np.sin(p) * np.sin(p2))
    return np.degrees(p2), np.degrees(l2) % 360


_ax = (np.arange(BOX_N) - BOX_N // 2) * BOX_D
BX, BY = np.meshgrid(_ax, _ax, indexing="xy")     # BX forward (or east), BY left (or north)
BR = np.hypot(BX, BY)
BANG = np.degrees(np.arctan2(-BY, BX))             # clockwise from +x


def sample(field, lat, lon, order):
    row = (90 - lat) / 0.25
    col = lon / 0.25
    f = np.concatenate([field, field[:, :1]], 1)     # wrap one column for 359.75 -> 360
    return ndimage.map_coordinates(f, [row, col], order=order, mode="nearest")


def box(field, clat, clon, rot, order=1):
    """Field on the 121x121 box; rot = bearing of +x (heading, or 90 for north-up east-x)."""
    la, lo = destination(clat, clon, rot + BANG, BR)
    return sample(field, la, lo, order)



KT_ = KT


def one_time(job):
    stamp, fixes = job
    spath = os.path.join(WORK, "stats", f"{stamp}.csv")
    if os.path.exists(spath):
        return stamp, 0
    t = np.datetime64(f"{stamp[:4]}-{stamp[4:6]}-{stamp[6:8]}T{stamp[8:10]}")
    ti = int((t - ARCO_T0) / np.timedelta64(1, "h"))
    nbytes = 0
    F = {}
    for name, var in (("msl", "mean_sea_level_pressure"), ("gust", "instantaneous_10m_wind_gust"),
                      ("u", "10m_u_component_of_wind"), ("v", "10m_v_component_of_wind"), ("d2m", "2m_dewpoint_temperature")):
        F[name], nb = fetch(f"{ARCO}/{var}/{ti}.0.0")
        nbytes += nb
    msl = F["msl"] / 100.0
    gust = F["gust"] * KT
    ws = np.hypot(F["u"], F["v"]) * KT
    d2m = F["d2m"] - 273.15

    lows = detect(msl)
    cen = np.array([XYZ[i, j] for _, i, j in lows])
    d, owner = cKDTree(cen).query(XYZ[R0:R1].reshape(-1, 3), k=1)
    dnear = (2 * np.arcsin(np.clip(d / 2, 0, 1)) * RE).reshape(R1 - R0, 1440)
    owner = owner.reshape(R1 - R0, 1440)

    rows, boxes = [], {}
    for fx in fixes:
        v = XYZ[int(round((90 - fx["lat"]) / 0.25)), int(round(fx["lon"] / 0.25)) % 1440]
        dk = np.arccos(np.clip(cen @ v, -1, 1)) * RE
        k = int(np.argmin(dk))
        own = np.zeros((721, 1440), bool)
        own[R0:R1] = (owner == k) & (dnear <= 1200) & OCEAN[R0:R1]
        _, ci, cj = lows[k]
        clat, clon = LAT[ci], LON[cj]
        G, W, O = gust[R0:R1], ws[R0:R1], own[R0:R1]
        dist = np.arccos(np.clip(BAND @ XYZ[ci, cj], -1, 1)) * RE
        brg = bearing(clat, clon, BLAT, BLON)
        hd = fx["heading"]
        rel = (brg - hd) % 360 if np.isfinite(hd) else np.full_like(brg, np.nan)
        r1200 = O & (dist <= 1200)
        r800 = O & (dist <= 800)
        g = np.where(r1200, G, -1.0)
        im = np.unravel_index(np.argmax(g), g.shape)
        gmax = float(G[im])
        w = np.where(r1200, W, -1.0)
        iw = np.unravel_index(np.argmax(w), w.shape)
        cd, ck = COAST.query(BAND[im])
        lat_c, lon_c = LANDPTS[ck]
        cell = BCELL
        m48 = r1200 & (G >= 48.0)
        row = dict(time=stamp, track=fx["track"], anchor=fx["anchor"], basin=fx["basin"], lat=fx["lat"], lon=fx["lon"],
                   heading=fx["heading"], g800_cat=fx["g800_cat"], match_km=round(float(dk[k]), 1),
                   g800=round(float(np.where(r800, G, 0).max()), 1),
                   gmax=round(gmax, 1), gmax_r=round(float(dist[im]), 0),
                   wsmax=round(float(W[iw]), 1), wsmax_r=round(float(dist[iw]), 0),
                   gust_factor=round(gmax / float(W[iw]), 3) if W[iw] > 0 else np.nan,
                   gmax_coast_km=round(float(cd * RE), 0),
                   gmax_near_grn_ice=bool(cd * RE <= 100 and lat_c >= 59.5 and 303 <= lon_c <= 352),
                   a_hf=round(float(cell[r1200 & (G >= HF)].sum()), 0),
                   a_g48=round(float(cell[m48].sum()), 0),
                   g48_rmax=round(float(dist[m48].max()), 0) if m48.any() else np.nan,
                   g48_right=round(float(cell[m48 & (rel < 180)].sum() / cell[m48].sum()), 3) if (m48.any() and np.isfinite(hd)) else np.nan,
                   own_ocean_frac=round(float(r1200.sum() / max((dist <= 1200).sum(), 1)), 3))
        b = {}
        for rot, tag in ((hd, "m"), (90.0, "n")):
            for name, fld in (("gust", gust), ("msl", msl), ("ws", ws), ("d2m", d2m)):
                if np.isfinite(rot):
                    b[f"{name}_{tag}"] = box(fld, clat, clon, rot).astype(np.float16)
                else:
                    b[f"{name}_{tag}"] = np.full((BOX_N, BOX_N), np.nan, np.float16)
            if np.isfinite(rot):
                b[f"own_{tag}"] = box(own.astype(np.float32), clat, clon, rot, 0).astype(np.uint8)
            else:
                b[f"own_{tag}"] = np.zeros((BOX_N, BOX_N), np.uint8)
        nb_ = b["msl_n"].astype(np.float32)
        row["pc"] = round(float(np.nanmin(np.where(BR <= 100, nb_, np.nan))), 1)
        ring = (BR >= 450) & (BR <= 550)
        row["msl_ring500"] = round(float(np.nanmean(nb_[ring])), 2)
        row["msl_grad"] = round((row["msl_ring500"] - row["pc"]) / 5.0, 3)
        row["d2m_500"] = round(float(np.nanmean(b["d2m_n"].astype(np.float32)[BR <= 500])), 2)
        rows.append(row)
        for kk, vv in b.items():
            boxes[f"{fx['track']}_{fx['anchor']}_{kk}"] = vv
    np.savez_compressed(os.path.join(WORK, "boxes", f"{stamp}.npz"), **boxes)
    tmp = spath + ".tmp"
    pd.DataFrame(rows).to_csv(tmp, index=False)
    os.replace(tmp, spath)
    with open(os.path.join(WORK, "bytes.log"), "a") as f:
        f.write(f"{stamp},{nbytes}\n")
    return stamp, nbytes


def jobs():
    T = pd.read_csv(os.environ.get("HFVS_TIMES", os.path.join(OWN, "times_storm.csv")))
    S = pd.read_csv(os.path.join(OWN, "storms.csv"))
    C = pd.read_csv(os.path.join(RES, "era5_hf_catalog_tracks.csv"))
    pk = S.set_index("track").gust800_kt
    on = C.rename(columns=dict(time="time_c", g800="g800_on"))[["track", "time_c", "g800_on"]]
    T = T.merge(on, left_on=["track", "time"], right_on=["track", "time_c"], how="left")
    T["g800_cat"] = np.where(T.anchor == "HF_onset", T.g800_on, T.track.map(pk))
    T = T.drop_duplicates(["track", "anchor"])
    out = []
    for t, g in T.groupby("time"):
        out.append((str(t), g[["track", "anchor", "basin", "lat", "lon", "heading", "g800_cat"]].to_dict("records")))
    return out


os.makedirs(os.path.join(WORK, "stats"), exist_ok=True)
os.makedirs(os.path.join(WORK, "boxes"), exist_ok=True)
LSM = land_sea()
OCEAN = LSM < 0.5
_land = np.argwhere(~OCEAN)
LANDPTS = np.stack([LAT[_land[:, 0]], LON[_land[:, 1]]], 1)
COAST = cKDTree(XYZ[~OCEAN])
BAND = XYZ[R0:R1]
BLAT = np.broadcast_to(LAT[R0:R1, None], (R1 - R0, 1440))
BLON = np.broadcast_to(LON[None, :], (R1 - R0, 1440))
BCELL = np.broadcast_to(CELL[R0:R1, None], (R1 - R0, 1440))

if __name__ == "__main__":
    J = jobs()
    if "--dry" in sys.argv:
        print(f"{len(J)} times, {sum(len(j[1]) for j in J)} anchors")
        sys.exit()
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    lim = int(sys.argv[2]) if len(sys.argv) > 2 else None
    if lim:
        J = J[:lim]
    t0, tot = time.time(), 0
    with Pool(n) as p:
        for i, (s, nb) in enumerate(p.imap_unordered(one_time, J, chunksize=2)):
            tot += nb
            if i % 50 == 0:
                print(f"{i + 1}/{len(J)} {tot / 1e9:.2f} GB {time.time() - t0:.0f}s", flush=True)
    print(f"done {len(J)} times, {tot / 1e9:.2f} GB this run")
