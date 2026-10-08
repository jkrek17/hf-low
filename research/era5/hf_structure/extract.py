"""Storm-relative ERA5 wind fields at pipeline A's hurricane-force-strength times.

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
WORK = os.path.join(os.environ.get("ERA5_WORK", os.path.join(HERE, "work")), "hf_structure")

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


def one_time(job):
    stamp, fixes, want_wind = job
    spath = os.path.join(WORK, "stats", f"{stamp}.csv")
    if os.path.exists(spath):
        return stamp, 0
    t = np.datetime64(f"{stamp[:4]}-{stamp[4:6]}-{stamp[6:8]}T{stamp[8:10]}")
    ti = int((t - ARCO_T0) / np.timedelta64(1, "h"))
    msl, b1 = fetch(f"{ARCO}/mean_sea_level_pressure/{ti}.0.0")
    gust, b2 = fetch(f"{ARCO}/instantaneous_10m_wind_gust/{ti}.0.0")
    msl = msl / 100.0
    gust = gust * KT
    nbytes = b1 + b2
    ws = None
    if want_wind:
        wi = int((t - WB2_T0) / np.timedelta64(6, "h"))
        ws, b3 = fetch(f"{WB2}/10m_wind_speed/{wi}.0.0")
        ws = ws * KT
        nbytes += b3

    lows = detect(msl)
    cen = np.array([XYZ[i, j] for _, i, j in lows])
    d, owner = cKDTree(cen).query(XYZ[R0:R1].reshape(-1, 3), k=1)
    dnear = (2 * np.arcsin(np.clip(d / 2, 0, 1)) * RE).reshape(R1 - R0, 1440)
    owner = owner.reshape(R1 - R0, 1440)
    owned_any = np.zeros((721, 1440), bool)

    rows, boxes = [], {}
    for fx in fixes:
        # match the catalog point to a re-detected low
        v = XYZ[int(round((90 - fx["lat"]) / 0.25)), int(round(fx["lon"] / 0.25)) % 1440]
        dk = np.arccos(np.clip(cen @ v, -1, 1)) * RE
        k = int(np.argmin(dk))
        own = np.zeros((721, 1440), bool)
        own[R0:R1] = (owner == k) & (dnear <= 1200) & OCEAN[R0:R1]
        _, ci, cj = lows[k]
        clat, clon = LAT[ci], LON[cj]
        # statistics on pipeline A's band only: no owned point lies outside it
        G, O = gust[R0:R1], own[R0:R1]
        dist = np.arccos(np.clip(BAND @ XYZ[ci, cj], -1, 1)) * RE
        brg = bearing(clat, clon, BLAT, BLON)
        rel = (brg - fx["heading"]) % 360                      # clockwise from motion
        r1200 = O & (dist <= 1200)
        r800 = O & (dist <= 800)
        g = np.where(r1200, G, -1.0)
        im = np.unravel_index(np.argmax(g), g.shape)
        gmax = float(G[im])
        cd, ck = COAST.query(BAND[im])
        row = dict(time=stamp, track=fx["track"], basin=fx["basin"], season=fx["season"],
                   lat=fx["lat"], lon=fx["lon"], msl=fx["msl"], g800_cat=fx["g800"],
                   heading=round(fx["heading"], 1), speed_kt=round(fx["speed"], 1),
                   match_km=round(float(dk[k]), 1), low_msl=round(float(lows[k][0]), 1),
                   g800=round(float(np.where(r800, G, 0).max()), 1),
                   gmax=round(gmax, 1), gmax_r=round(float(dist[im]), 0),
                   gmax_rel=round(float(rel[im]), 0), gmax_brg=round(float(brg[im]), 0),
                   gmax_lat=LAT[R0 + im[0]], gmax_lon=LON[im[1]],
                   gmax_coast_km=round(float(cd * RE), 0),
                   gmax_coast_lat=LANDPTS[ck][0], gmax_coast_lon=LANDPTS[ck][1])
        cell = BCELL
        for thr, key in ((HF, "a_hf"), (64.0, "a_g64"), (48.0, "a_g48"), (34.0, "a_g34")):
            m = r1200 & (G >= thr)
            row[key] = round(float(cell[m].sum()), 0)
            if key == "a_hf":
                hm = m
        for frac, key in ((0.9, "a_p90"), (0.8, "a_p80")):
            row[key] = round(float(cell[r1200 & (G >= frac * gmax)].sum()), 0)
        # where the HF-equivalent gust area sits, by motion-relative and earth quadrant
        rq, bq, ca, dd = rel[hm], brg[hm], cell[hm], dist[hm]
        for name, lo_, hi_ in (("fr", 0, 90), ("rr", 90, 180), ("rl", 180, 270), ("fl", 270, 360)):
            row["q_" + name] = round(float(ca[(rq >= lo_) & (rq < hi_)].sum()), 0)
        for name, lo_, hi_ in (("ne", 0, 90), ("se", 90, 180), ("sw", 180, 270), ("nw", 270, 360)):
            row["e_" + name] = round(float(ca[(bq >= lo_) & (bq < hi_)].sum()), 0)
        row["hf_r50"] = round(float(np.median(dd)), 0) if len(dd) else np.nan
        row["hf_rmax"] = round(float(dd.max()), 0) if len(dd) else np.nan
        sf = r1200 & (G >= 48.0)
        row["g48_rmax"] = round(float(dist[sf].max()), 0) if sf.any() else np.nan
        row["hf_coast50"] = round(float((COASTD[hm] <= 50).mean()), 3) if hm.any() else np.nan
        row["own_ocean_frac"] = round(float(r1200.sum() / max((dist <= 1200).sum(), 1)), 3)
        if ws is not None:
            W = ws[R0:R1]
            w = np.where(r1200, W, -1.0)
            iw = np.unravel_index(np.argmax(w), w.shape)
            row.update(wmax=round(float(W[iw]), 1), wmax_r=round(float(dist[iw]), 0),
                       wmax_rel=round(float(rel[iw]), 0),
                       ws_at_gmax=round(float(W[im]), 1))
            for thr, key in ((64.0, "a_w64"), (51.2, "a_w51"), (48.0, "a_w48"), (34.0, "a_w34")):
                row[key] = round(float(cell[r1200 & (W >= thr)].sum()), 0)
        rows.append(row)

        b = dict(gust_m=box(gust, clat, clon, fx["heading"]).astype(np.float16),
                 msl_m=box(msl, clat, clon, fx["heading"]).astype(np.float16),
                 own_m=box(own.astype(np.float32), clat, clon, fx["heading"], 0).astype(np.uint8),
                 gust_n=box(gust, clat, clon, 90.0).astype(np.float16),
                 own_n=box(own.astype(np.float32), clat, clon, 90.0, 0).astype(np.uint8))
        if ws is not None:
            b["ws_m"] = box(ws, clat, clon, fx["heading"]).astype(np.float16)
        for kk, vv in b.items():
            boxes[f"{fx['track']}_{kk}"] = vv
    np.savez_compressed(os.path.join(WORK, "boxes", f"{stamp}.npz"), **boxes)
    tmp = spath + ".tmp"
    pd.DataFrame(rows).to_csv(tmp, index=False)
    os.replace(tmp, spath)
    with open(os.path.join(WORK, "bytes.log"), "a") as f:
        f.write(f"{stamp},{nbytes}\n")
    return stamp, nbytes


def jobs():
    P = pd.read_csv(os.path.join(RES, "era5_hf_catalog_tracks.csv"))
    C = pd.read_csv(os.path.join(RES, "era5_hf_catalog.csv"))
    P = P.merge(C[["track", "role", "season"]], on="track")
    P = P[P.role == "event"].sort_values(["track", "time"]).reset_index(drop=True)
    P["dt"] = pd.to_datetime(P.time.astype(str), format="%Y%m%d%H")
    hd, sp = np.full(len(P), np.nan), np.full(len(P), np.nan)
    for _, g in P.groupby("track", sort=False):
        ix = g.index.values
        la, lo, tt = g.lat.values, g.lon.values, g.dt.values
        for n, i in enumerate(ix):
            a, b = max(n - 1, 0), min(n + 1, len(ix) - 1)
            if a == b:
                continue
            hd[i] = bearing(la[a], lo[a], la[b], lo[b])
            dkm = np.arccos(np.clip(XYZ_pt(la[a], lo[a]) @ XYZ_pt(la[b], lo[b]), -1, 1)) * RE
            sp[i] = dkm / ((tt[b] - tt[a]) / np.timedelta64(1, "h")) / 1.852
    P["heading"], P["speed"] = hd, sp
    # in-domain points only: pipeline A's index ignores points outside its basin boxes
    H = P[(P.season >= SEASON0) & (P.g800 >= HF) & P.basin.notna()].copy()
    times = np.sort(H.time.unique())
    early = H[H.dt < WB2_END].groupby("track").time.apply(list)
    rng = np.random.default_rng(20261008)
    pick = set()
    for tr in rng.permutation(early.index.values):          # one time per storm, then repeat
        pick.add(rng.choice(early[tr]))
        if len(pick) >= N_WIND:
            break
    out = []
    for t in times:
        fx = H[H.time == t][["track", "basin", "season", "lat", "lon", "msl", "g800", "heading", "speed"]]
        out.append((str(t), fx.to_dict("records"), t in pick))
    return out


def XYZ_pt(lat, lon):
    a, b = np.radians(lat), np.radians(lon)
    return np.array([np.cos(a) * np.cos(b), np.cos(a) * np.sin(b), np.sin(a)])


os.makedirs(os.path.join(WORK, "stats"), exist_ok=True)
os.makedirs(os.path.join(WORK, "boxes"), exist_ok=True)
LSM = land_sea()
OCEAN = LSM < 0.5
_land = np.argwhere(~OCEAN)
LANDPTS = np.stack([LAT[_land[:, 0]], LON[_land[:, 1]]], 1)
COAST = cKDTree(XYZ[~OCEAN])                                # chord distance on the unit sphere
BAND = XYZ[R0:R1]
BLAT = np.broadcast_to(LAT[R0:R1, None], (R1 - R0, 1440))
BLON = np.broadcast_to(LON[None, :], (R1 - R0, 1440))
BCELL = np.broadcast_to(CELL[R0:R1, None], (R1 - R0, 1440))
COASTD = COAST.query(BAND.reshape(-1, 3))[0].reshape(R1 - R0, 1440) * RE

if __name__ == "__main__":
    J = jobs()
    if "--dry" in sys.argv:
        print(f"{len(J)} times, {sum(len(j[1]) for j in J)} fixes, {sum(j[2] for j in J)} with wind speed")
        sys.exit()
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    t0, tot = time.time(), 0
    with Pool(n) as p:
        for i, (s, nb) in enumerate(p.imap_unordered(one_time, J, chunksize=4)):
            tot += nb
            if i % 100 == 0:
                print(f"{i + 1}/{len(J)} {tot / 1e9:.2f} GB {time.time() - t0:.0f}s", flush=True)
    print(f"done {len(J)} times, {tot / 1e9:.2f} GB this run")
