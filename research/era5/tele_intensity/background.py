"""Background sea-level pressure at each pipeline A event's deepest fix (ERA5 proxy).

For every event in the pipeline A catalog, take the in-domain fix with the
lowest central pressure (its time, position and `minp`, which equals the
catalog's minp) and add two background pressures from ERA5 MSLP on the
1.5 degree grid:

  clim    calendar-month climatology at that position, bilinear. Built from
          every 4th WeatherBench2 time chunk (8 consecutive 6-hourly fields
          every 8 days) over 1979-01 .. 2023-01, so about 1,000 chunks.
  ring    mean MSLP over the 900-1100 km annulus around the fix, at the same
          time, cos(lat)-weighted. The storm's surroundings at that moment.

So minp - clim is depth against climatology at that place and month, minp -
ring is depth against the concurrent environment, and ring - clim is the
large-scale background anomaly; the first is the sum of the other two.

Sources are the ones research/era5/intensity/env.py uses: WeatherBench2's
conservative 1.5 degree regrid to 2023-01-09, then hourly 0.25 degree ARCO
coarsened with the same cell edges (coarsen() is copied from env.py).

Volume: about 2,000 climatology chunks + ~3,800 event chunks at 0.93 MB, plus
~400 ARCO MSLP fields at ~2-4 MB: about 6-7 GB streamed. Nothing is kept but
the monthly climatology (work/, ignored) and the output CSV.

usage: background.py CATALOG_CSV TRACKS_CSV OUT_CSV [NTHREADS]
"""
import os, sys, json, time
import numpy as np, pandas as pd, numcodecs, gcsfs
from concurrent.futures import ThreadPoolExecutor

WB2 = "weatherbench2/datasets/era5/1959-2023_01_10-6h-240x121_equiangular_with_poles_conservative.zarr"
ARCO = "gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3"
WB2_T0 = np.datetime64("1959-01-01T00")
ARCO_T0 = np.datetime64("1900-01-01T00")
WB2_END = np.datetime64("2023-01-10T00")
RE = 6371.0
LAT = -90 + 1.5 * np.arange(121)
LON = 1.5 * np.arange(240)
LAT2, LON2 = np.meshgrid(LAT, LON, indexing="ij")
HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")

_fs = None


def cat(key):
    global _fs
    for k in range(7):
        try:
            if _fs is None:
                _fs = gcsfs.GCSFileSystem(token="anon")
            return _fs.cat(key)
        except Exception:
            time.sleep(2 ** k)
    raise RuntimeError(f"fetch failed {key}")


_meta = None


def wb2_chunk(c):
    """MSLP for WB2 time chunk c: array [8, 121 lat (-90..90), 240 lon] in hPa (not cached)."""
    global _meta
    if _meta is None:
        _meta = json.loads(cat(f"{WB2}/.zmetadata"))["metadata"]
    za = _meta["mean_sea_level_pressure/.zarray"]
    codec = numcodecs.get_codec(za["compressor"]) if za["compressor"] else None
    raw = cat(f"{WB2}/mean_sea_level_pressure/{c}.0.0")
    if codec is not None:
        raw = codec.decode(raw)
    a = np.frombuffer(raw, za["dtype"]).reshape(za["chunks"])        # time, lon, lat
    nt = min(za["chunks"][0], za["shape"][0] - c * za["chunks"][0])
    return (a[:nt].transpose(0, 2, 1) / 100.0).astype(np.float32)


_W7 = np.array([0.5, 1, 1, 1, 1, 1, 0.5]) / 6.0


def coarsen(a):
    """0.25 degree [721, 1440] (lat 90..-90) -> 1.5 degree [121, 240] (lat -90..90). From env.py."""
    a = a.astype(np.float64)
    nan = np.isnan(a)
    w = (~nan).astype(np.float64)
    a = np.where(nan, 0.0, a)
    def smooth(x):
        xl = np.concatenate([x[:, -3:], x, x[:, :3]], axis=1)
        y = sum(_W7[k] * xl[:, k:k + 1440] for k in range(7))[:, ::6]
        yl = np.concatenate([y[3:0:-1], y, y[-2:-5:-1]], axis=0)
        z = sum(_W7[k] * yl[k:k + 721] for k in range(7))[::6]
        return z
    num, den = smooth(a), smooth(w)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.where(den > 0.5 * smooth(np.ones_like(w)), num / den, np.nan)
    return out[::-1]


def arco_field(t):
    """MSLP at hourly time t from ARCO, coarsened to 1.5 degrees, hPa, cached."""
    ti = int((t - ARCO_T0) / np.timedelta64(1, "h"))
    a = np.frombuffer(numcodecs.Blosc().decode(cat(f"{ARCO}/mean_sea_level_pressure/{ti}.0.0")), "<f4").reshape(721, 1440)
    return (coarsen(a) / 100.0).astype(np.float32)


def bilinear(F, lat, lon):
    i = (lat + 90) / 1.5
    j = (lon % 360) / 1.5
    i0, j0 = int(np.floor(i)), int(np.floor(j))
    di, dj = i - i0, j - j0
    j1 = (j0 + 1) % 240
    return ((1 - di) * (1 - dj) * F[i0, j0] + (1 - di) * dj * F[i0, j1]
            + di * (1 - dj) * F[i0 + 1, j0] + di * dj * F[i0 + 1, j1])


def ring_mean(F, lat, lon, r0=900.0, r1=1100.0):
    p1, p2 = np.radians(lat), np.radians(LAT2)
    dl = np.radians(LON2 - lon)
    d = RE * np.arccos(np.clip(np.sin(p1) * np.sin(p2) + np.cos(p1) * np.cos(p2) * np.cos(dl), -1, 1))
    m = (d >= r0) & (d <= r1)
    w = np.cos(p2[m])
    return float(np.sum(F[m] * w) / np.sum(w))


def main():
    cat_csv, tracks_csv, out_csv = sys.argv[1:4]
    nth = int(sys.argv[4]) if len(sys.argv) > 4 else 24
    os.makedirs(WORK, exist_ok=True)
    C = pd.read_csv(cat_csv)
    E = C[C.role == "event"][["track"]]
    T = pd.read_csv(tracks_csv, dtype={"time": str})
    T = T[T.basin.notna() & T.track.isin(E.track)]
    D = T.loc[T.groupby("track").msl.idxmin(), ["track", "time", "lat", "lon", "msl"]].reset_index(drop=True)
    D["t"] = pd.to_datetime(D.time, format="%Y%m%d%H").values.astype("datetime64[h]")
    print(f"{len(D)} events", flush=True)

    # climatology: every 4th WB2 chunk from 1979-01-01 to the end of WB2
    c0 = int((np.datetime64("1979-01-01T00") - WB2_T0) / np.timedelta64(48, "h"))
    c1 = int((WB2_END - WB2_T0) / np.timedelta64(48, "h"))
    cl_chunks = list(range(c0, c1, 4))
    print(f"climatology chunks {len(cl_chunks)}, ARCO fields {(D.t >= WB2_END).sum()}", flush=True)

    def clim_part(c):
        a = wb2_chunk(c)
        S = np.zeros((12, 121, 240)); N = np.zeros(12)
        for k in range(a.shape[0]):
            m = (WB2_T0 + np.timedelta64(48 * c + 6 * k, "h")).astype("datetime64[M]").astype(int) % 12
            S[m] += a[k]; N[m] += 1
        return S, N
    f = os.path.join(WORK, "clim_month.npz")
    if os.path.exists(f):
        z = np.load(f); clim, N = z["clim"], z["n"]
    else:
        S = np.zeros((12, 121, 240)); N = np.zeros(12)
        with ThreadPoolExecutor(nth) as ex:
            for s_, n_ in ex.map(clim_part, cl_chunks):
                S += s_; N += n_
        clim = S / N[:, None, None]
        np.savez(f, clim=clim, n=N)
    print("climatology done", flush=True)

    def one(r):
        t = np.datetime64(r.t, "h")
        if t < WB2_END:
            h = int((t - WB2_T0) / np.timedelta64(6, "h"))
            F = wb2_chunk(h // 8)[h % 8]
        else:
            F = arco_field(t)
        m = t.astype("datetime64[M]").astype(int) % 12
        return dict(track=r.track, tmin=r.time, lat_min=r.lat, lon_min=r.lon, minp=r.msl,
                    clim=round(bilinear(clim[m], r.lat, r.lon), 2),
                    ring=round(ring_mean(F, r.lat, r.lon), 2),
                    src="wb2" if t < WB2_END else "arco")
    with ThreadPoolExecutor(nth) as ex:
        rows = list(ex.map(one, [r for r in D.itertuples()]))
    pd.DataFrame(rows).to_csv(out_csv, index=False)
    print(f"wrote {out_csv}; climatology fields per month {N.astype(int).tolist()}")


if __name__ == "__main__":
    main()
