"""Shared helpers: WeatherBench2 reads, grids, storm-centred boxes."""
import os, time, json
import numpy as np
import gcsfs, numcodecs
from scipy import ndimage

WORK = os.environ.get("ERA5_WORK", os.path.join(os.path.dirname(os.path.abspath(__file__)), "work"))
WB2 = "weatherbench2/datasets/era5/1959-2023_01_10-6h-240x121_equiangular_with_poles_conservative.zarr"
CLIM = ("weatherbench2/datasets/era5-hourly-climatology/"
        "1990-2019_6h_240x121_equiangular_with_poles_conservative.zarr")
WB2_T0 = np.datetime64("1959-01-01T00")
LEVELS = (50, 100, 150, 200, 250, 300, 400, 500, 600, 700, 850, 925, 1000)
G = 9.80665
RE = 6371.0
LAT = -90 + 1.5 * np.arange(121)
LON = 1.5 * np.arange(240)
BOX_N, BOX_D = 81, 100.0                      # +-4000 km
MASK_LAT = 85.0

_fs = None


def fs():
    global _fs
    if _fs is None:
        _fs = gcsfs.GCSFileSystem(token="anon")
    return _fs


def cat(key):
    for k in range(6):
        try:
            return fs().cat(key)
        except Exception:
            time.sleep(2 ** k)
    raise RuntimeError(f"fetch failed {key}")


_meta = {}


def zarray(root, var):
    k = (root, var)
    if k not in _meta:
        _meta[k] = json.loads(cat(f"{root}/{var}/.zarray"))
    return _meta[k]


def decode(root, var, key):
    za = zarray(root, var)
    raw = cat(f"{root}/{var}/{key}")
    a = np.frombuffer(numcodecs.get_codec(za["compressor"]).decode(raw), za["dtype"]).reshape(za["chunks"])
    return a, len(raw)


# ---- box geometry (same convention as hf_structure/extract.py) ----
_ax = (np.arange(BOX_N) - BOX_N // 2) * BOX_D
BX, BY = np.meshgrid(_ax, _ax, indexing="xy")      # BX forward (or east), BY left (or north)
BR = np.hypot(BX, BY)
BANG = np.degrees(np.arctan2(-BY, BX))              # clockwise from +x


def destination(lat, lon, brg, dkm):
    p, l, b, d = np.radians(lat), np.radians(lon), np.radians(brg), dkm / RE
    p2 = np.arcsin(np.sin(p) * np.cos(d) + np.cos(p) * np.sin(d) * np.cos(b))
    l2 = l + np.arctan2(np.sin(b) * np.sin(d) * np.cos(p), np.cos(d) - np.sin(p) * np.sin(p2))
    return np.degrees(p2), np.degrees(l2) % 360


def sample(field, lat, lon):
    """Bilinear sample of field[121, 240] (lat -90..90, lon 0..358.5) at points; NaN where |lat| > MASK_LAT."""
    row = (lat + 90) / 1.5
    col = (lon % 360) / 1.5
    f = np.concatenate([field, field[:, :1]], 1)
    out = ndimage.map_coordinates(f, [row, col], order=1, mode="nearest")
    return np.where(np.abs(lat) > MASK_LAT, np.nan, out)


def divergence(u, v):
    """Horizontal divergence (1e-5 /s) on the [121, 240] grid; NaN on rows within 3 degrees of the poles."""
    dl = np.radians(1.5)
    cosl = np.cos(np.radians(LAT))[:, None]
    dudx = (np.roll(u, -1, 1) - np.roll(u, 1, 1)) / (2 * dl * RE * 1e3 * cosl)
    vc = v * cosl
    dvcdy = np.full_like(v, np.nan)
    dvcdy[1:-1] = (vc[2:] - vc[:-2]) / (2 * dl * RE * 1e3)
    d = (dudx + dvcdy / cosl) * 1e5
    d[np.abs(LAT) > 88.5] = np.nan
    return d


def read_small(root, var):
    """A 1-D coordinate array stored as one chunk."""
    za = zarray(root, var)
    raw = cat(f"{root}/{var}/0")
    if za["compressor"]:
        raw = numcodecs.get_codec(za["compressor"]).decode(raw)
    return np.frombuffer(raw, za["dtype"])
