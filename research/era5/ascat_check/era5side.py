"""ERA5 side of the ASCAT check: cache index, sample, disc geometry, bilinear gust at points.
ERA5 PROXY, pipeline A. Reads the shared gust cache; never pulls."""
import glob, os, numpy as np, pandas as pd
try:
    import numcodecs
except ImportError:                       # no PyPI in the 2026-10-14 environment: decode with blosc_nolib (system libzstd)
    numcodecs = None
    import blosc_nolib

CACHE = os.environ.get("ERA5_GUST_CACHE", "/mnt/project-files/era5-cache")
TRACKS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "hf_history", "results", "all_tracks.csv.gz")
RE = 6371.0
LAT = 90 - 0.25 * np.arange(361)          # 90N .. 0N
LON = 0.25 * np.arange(1440)
_blosc = numcodecs.Blosc() if numcodecs else None


def cache_index():
    """hours_since_1900 -> (file, row)."""
    idx = {}
    for f in sorted(glob.glob(f"{CACHE}/gust_nh_12h/gust_*.npz")):
        h = np.load(f)["hours_since_1900"]
        for r, v in enumerate(h):
            idx[int(v)] = (f, r)
    return idx


def hours1900(t):
    return int((np.datetime64(t, "h") - np.datetime64("1900-01-01T00", "h")) / np.timedelta64(1, "h"))


import functools


@functools.lru_cache(maxsize=2)
def _file(f):
    z = np.load(f)
    raw = z["data"].tobytes()
    return np.frombuffer(_blosc.decode(raw) if _blosc else blosc_nolib.decode(raw), np.uint8).reshape(z["shape"])


def field(f, r):
    """one 00/12 UTC gust field in kt (361x1440); decoded files are cached, so iterate in file order."""
    return _file(f)[r].astype(np.float32) / 2.0


def lsm_ocean():
    return np.load(f"{CACHE}/lsm_nh_0p25.npy") < 0.5


def sample(min_season=2004):
    d = pd.read_csv(TRACKS, dtype={"peak_time": str})
    d = d[(d.season >= min_season) & d.peak_time.str[-2:].isin(["00", "12"])].copy()
    d["time"] = pd.to_datetime(d.peak_time, format="%Y%m%d%H")
    d["h1900"] = [hours1900(t) for t in d.time.values]
    return d.reset_index(drop=True)


def gc_km(la1, lo1, la2, lo2):
    la1, lo1, la2, lo2 = map(np.radians, (la1, lo1, la2, lo2))
    a = np.sin((la2 - la1) / 2) ** 2 + np.cos(la1) * np.cos(la2) * np.sin((lo2 - lo1) / 2) ** 2
    return 2 * RE * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def disc_cells(lat, lon, radius=800.0):
    """(row, col) index arrays of all grid cells within radius km of (lat, lon); lon in 0-360."""
    dlat = np.degrees(radius / RE)
    i0 = max(0, int(np.floor((90 - (lat + dlat)) / 0.25)))
    i1 = min(360, int(np.ceil((90 - (lat - dlat)) / 0.25)))
    dlon = min(180.0, dlat / max(np.cos(np.radians(min(abs(lat) + dlat, 89))), 0.05))
    n = int(np.ceil(dlon / 0.25))
    j0 = int(round((lon % 360) / 0.25))
    cols = (j0 + np.arange(-n, n + 1)) % 1440
    rows = np.arange(i0, i1 + 1)
    R, C = np.meshgrid(rows, cols, indexing="ij")
    d = gc_km(LAT[R], LON[C], lat, lon)
    m = d <= radius
    return R[m], C[m]


def bilinear(field, lat, lon):
    """field (361x1440) at points; lat,lon arrays (lon any, wrapped)."""
    x = (np.asarray(lon) % 360) / 0.25
    y = (90 - np.asarray(lat)) / 0.25
    y = np.clip(y, 0, 359.999)
    i, j = y.astype(int), x.astype(int)
    fy, fx = y - i, x - j
    j1 = (j + 1) % 1440
    return (field[i, j] * (1 - fy) * (1 - fx) + field[i + 1, j] * fy * (1 - fx)
            + field[i, j1] * (1 - fy) * fx + field[i + 1, j1] * fy * fx)
