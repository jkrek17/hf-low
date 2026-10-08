"""Monthly SST-front strength (Kuroshio-Oyashio zone) from ERA5 0.25 degree SST (resumable).

Plan: PREREGISTRATION.md. 9 samples a month (12 UTC, days 3,6,...,27), Sep-Feb of each season 2004-05..2025-26
(about 1.6 GB). Smooth monthly mean SST with a 1x1 degree NaN-aware box filter, G = -dSST/dy in K per 100 km, max over a latitude
band at each longitude in 150-170E, mean over longitudes. Bands: FRONT 32-46N, FRONT_KE 32-40N, FRONT_OE 38-46N.
usage: sst_front.py OUT_CSV
"""
import sys, os
import numpy as np, pandas as pd
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "enso_kuroshio"))
from common import get, tindex, LAT, LON, BYTES
from scipy.ndimage import uniform_filter

DAYS = (3, 6, 9, 12, 15, 18, 21, 24, 27)
BANDS = {"front": (32, 46), "front_ke": (32, 40), "front_oe": (38, 46)}
R = (LAT >= 25) & (LAT <= 55); C = (LON >= 145) & (LON <= 175)
LATS = LAT[R]; LONS = LON[C]


def nan_smooth(a, n=5):
    ok = np.isfinite(a)
    num = uniform_filter(np.where(ok, a, 0.0), n, mode="nearest")
    den = uniform_filter(ok.astype(float), n, mode="nearest")
    out = np.where(den > 0.5, num / np.maximum(den, 1e-9), np.nan)
    return out


def indices(sst):
    s = nan_smooth(sst)
    dy = np.gradient(s, axis=0)           # K per row; rows go north to south (LAT descending)
    # latitude descends with row index, so d/dlat = -dy / 0.25 deg; G = -dSST/dlat(north) -> positive where SST falls northward
    G = dy / 0.25 / 111.19 * 100          # K per 100 km, positive when SST rises with row index (i.e. falls northward)
    sel = (LONS >= 150) & (LONS <= 170)
    out = {}
    for k, (a, b) in BANDS.items():
        rr = (LATS >= a) & (LATS <= b)
        g = G[np.ix_(rr, sel)]
        out[k] = float(np.nanmean(np.nanmax(g, axis=0)))
    return out


def one(args):
    y, m = args
    acc = {k: [] for k in BANDS}
    for d in DAYS:
        a = get("sea_surface_temperature", tindex(f"{y}-{m:02d}-{d:02d}T12"))[np.ix_(R, C)]
        for k, v in indices(a).items():
            acc[k].append(v)
    return dict(year=y, month=m, **{k: float(np.mean(v)) for k, v in acc.items()})


if __name__ == "__main__":
    out = sys.argv[1]
    months = []
    for s in range(2004, 2026):
        months += [(s, 9), (s, 10), (s, 11), (s, 12), (s + 1, 1), (s + 1, 2)]
    with ThreadPoolExecutor(8) as ex:
        rows = list(ex.map(one, months))
    pd.DataFrame(rows).to_csv(out, index=False, float_format="%.5f")
    print(len(rows), "months", round(BYTES[0] / 1e9, 2), "GB")
