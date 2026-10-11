#!/usr/bin/env python3
"""Pull the ERA5 fields behind the storm-review pages (docs/storm-reviews.html).

For each storm listed in storm_review_storms.csv (pipeline A track ids from the strongest-storms table),
fetch, at 0.25 degrees:
  WeatherBench2 6-hourly MSLP, Z500, u250, v250 at peak -48, -36, -24, -12, 0, +12, +24 h  (cropped around the storm)
  ARCO-ERA5 hourly instantaneous 10 m gust, every 3 h from peak -48 to +24 h               (running maximum, kt)
Writes work/<track>.npz (ignored, never committed). Resumable. ERA5 is a proxy; nothing here is observation.

usage: storm_review_extract.py TRACKS_CSV FIXES_CSV OUTDIR
"""
import os, sys, time, urllib.request, json
import numpy as np, pandas as pd, numcodecs

WB2 = "https://storage.googleapis.com/weatherbench2/datasets/era5/1959-2023_01_10-wb13-6h-1440x721_with_derived_variables.zarr"
ARCO = "https://storage.googleapis.com/gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3"
WB2_T0 = np.datetime64("1959-01-01T00"); ARCO_T0 = np.datetime64("1900-01-01T00")
LEVELS = (50, 100, 150, 200, 250, 300, 400, 500, 600, 700, 850, 925, 1000)
KT = 1 / 0.514444
codec = numcodecs.Blosc()
BYTES = [0]
MAP_H = (-48, -36, -24, -12, 0, 12, 24)
MAP_DLAT, MAP_DLON = 22.0, 40.0           # half-widths of the map crops, degrees
SW_H = tuple(range(-48, 25, 3))


def get(url, shape):
    for k in range(6):
        try:
            with urllib.request.urlopen(url, timeout=180) as r:
                raw = r.read()
            BYTES[0] += len(raw)
            return np.frombuffer(codec.decode(raw), "<f4").reshape(shape)
        except Exception:
            if k == 5:
                raise
            time.sleep(2 ** k)


def crop(a, lat0, lon0, dlat, dlon):
    """a[..., 721, 1440] -> rows lat0+dlat .. lat0-dlat, columns lon0-dlon .. lon0+dlon (wrapping)."""
    r0 = int(round((90 - (lat0 + dlat)) / 0.25)); r1 = int(round((90 - (lat0 - dlat)) / 0.25)) + 1
    r0, r1 = max(r0, 0), min(r1, 721)
    c = (np.arange(int(round((lon0 - dlon) / 0.25)), int(round((lon0 + dlon) / 0.25)) + 1)) % 1440
    return a[..., r0:r1, :][..., c], 90 - 0.25 * np.arange(r0, r1), (c * 0.25)


def centre_at(fx, t):
    tt = fx.t.values.astype("datetime64[h]").astype(np.int64)
    x = t.astype("datetime64[h]").astype(np.int64)
    lon = np.unwrap(np.radians(fx.lon.values)); lon = np.degrees(lon)
    return float(np.interp(x, tt, fx.lat.values)), float(np.interp(x, tt, lon)) % 360


def main(tracks_csv, fixes_csv, outdir):
    os.makedirs(outdir, exist_ok=True)
    st = pd.read_csv(tracks_csv)
    fx_all = pd.read_csv(fixes_csv)
    fx_all["t"] = pd.to_datetime(fx_all.time.astype(str), format="%Y%m%d%H")
    for s in st.itertuples():
        path = os.path.join(outdir, f"{s.track}.npz")
        if os.path.exists(path):
            continue
        fx = fx_all[fx_all.track == s.track].sort_values("t").reset_index(drop=True)
        t0 = np.datetime64(fx.t[fx.msl.idxmin()].isoformat(), "h")
        out = {}
        for h in MAP_H:
            t = t0 + np.timedelta64(h, "h")
            la, lo = centre_at(fx, t)
            ti = int((t - WB2_T0) / np.timedelta64(6, "h"))
            msl = get(f"{WB2}/mean_sea_level_pressure/{ti}.0.0", (1, 721, 1440))[0] / 100
            z = get(f"{WB2}/geopotential/{ti}.0.0.0", (1, 13, 721, 1440))[0, LEVELS.index(500)] / 9.80665
            u = get(f"{WB2}/u_component_of_wind/{ti}.0.0.0", (1, 13, 721, 1440))[0, LEVELS.index(250)]
            v = get(f"{WB2}/v_component_of_wind/{ti}.0.0.0", (1, 13, 721, 1440))[0, LEVELS.index(250)]
            m, lats, lons = crop(msl, la, lo, MAP_DLAT, MAP_DLON)
            out[f"msl_{h}"] = m.astype(np.float32)
            out[f"z500_{h}"] = crop(z, la, lo, MAP_DLAT, MAP_DLON)[0].astype(np.float32)
            out[f"ws250_{h}"] = crop(np.hypot(u, v), la, lo, MAP_DLAT, MAP_DLON)[0].astype(np.float32)
            out[f"centre_{h}"] = np.array([la, lo])
            out["lats"], out["lons"] = lats, lons
        # swath: running max gust over the track's box
        tr = fx[(fx.t >= (t0 + np.timedelta64(-48, "h")).astype("datetime64[s]")) & (fx.t <= (t0 + np.timedelta64(24, "h")).astype("datetime64[s]"))]
        lat_c = float(tr.lat.mean()); lon_unw = np.degrees(np.unwrap(np.radians(tr.lon.values)))
        lon_c = float(lon_unw.mean()) % 360
        dlat = float(max(tr.lat.max() - lat_c, lat_c - tr.lat.min())) + 14
        dlon = float(max(lon_unw.max() - lon_unw.mean(), lon_unw.mean() - lon_unw.min())) + 28
        mx = None
        for h in SW_H:
            t = t0 + np.timedelta64(h, "h")
            ti = int((t - ARCO_T0) / np.timedelta64(1, "h"))
            g = get(f"{ARCO}/instantaneous_10m_wind_gust/{ti}.0.0", (721, 1440)) * KT
            c, slat, slon = crop(g, lat_c, lon_c, dlat, dlon)
            mx = c if mx is None else np.maximum(mx, c)
        out.update(sw_gust=mx.astype(np.float32), sw_lats=slat, sw_lons=slon, t0=np.array([str(t0)]))
        np.savez_compressed(path + ".tmp.npz", **out); os.replace(path + ".tmp.npz", path)
        print(s.track, str(t0), f"cumulative {BYTES[0]/1e9:.2f} GB", flush=True)
        if BYTES[0] > 30e9:
            sys.exit("stop at 30 GB")


if __name__ == "__main__":
    main(*sys.argv[1:4])
