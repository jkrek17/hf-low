"""Shared ERA5 chunk access for the ENSO / Kuroshio study (urllib only, no gcsfs)."""
import os, time, urllib.request
import numpy as np, numcodecs

ROOT = "https://storage.googleapis.com/gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3"
T0 = np.datetime64("1900-01-01T00")
CODEC = numcodecs.Blosc()
LAT = 90 - 0.25 * np.arange(721)
LON = 0.25 * np.arange(1440)
BYTES = [0]


def tindex(ts):
    return int((np.datetime64(ts) - T0) / np.timedelta64(1, "h"))


def get(var, ti):
    url = f"{ROOT}/{var}/{ti}.0.0"
    for k in range(6):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                raw = r.read()
            BYTES[0] += len(raw)
            return np.frombuffer(CODEC.decode(raw), "<f4").reshape(721, 1440)
        except Exception:
            if k == 5:
                raise
            time.sleep(2 ** k)
