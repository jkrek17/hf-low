"""Pull 6-hourly ERA5 (1.5 degree) v-wind at 250 hPa or mean sea level pressure, 19.5-75N, resumable.

Source: ARCO-ERA5 `1959-2022-6h-240x121_equiangular_with_poles_conservative.zarr` (ends 2021-12-31).
Chunks hold 8 times (x 13 levels for v), so whole chunks are read and the level is selected after.
Output: WORK/<var>/blk_<n>.npy, float32 (time, 240 lon, 38 lat), plus WORK/<var>/times.npy.

    python3 extract_eddy.py v250  2004-06-01 2022-01-01
    python3 extract_eddy.py mslp  1979-06-01 2022-01-01
"""
import os, sys, time
import numpy as np, xarray as xr
import multiprocessing as mp

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.environ.get("ERA5_WORK", os.path.join(HERE, "work"))
P = "gs://gcp-public-data-arco-era5/ar/1959-2022-6h-240x121_equiangular_with_poles_conservative.zarr"
LAT0, LAT1 = 73, 111          # latitudes 19.5 .. 75.0 (index into -90..90 step 1.5)
BLK = 400                      # chunks (of 8 times) per saved block
_ds = None


def ds():
    global _ds
    if _ds is None:
        _ds = xr.open_zarr(P, storage_options={"token": "anon"}, chunks=None)
    return _ds


def fetch(args):
    var, a, b, out = args
    if os.path.exists(out):
        return out
    for k in range(5):
        try:
            d = ds()
            if var == "v250":
                x = d.v_component_of_wind.isel(time=slice(a, b), latitude=slice(LAT0, LAT1)).sel(level=250).values
            else:
                x = d.mean_sea_level_pressure.isel(time=slice(a, b), latitude=slice(LAT0, LAT1)).values
            np.save(out + ".tmp.npy", x.astype("float32"))
            os.replace(out + ".tmp.npy", out)
            return out
        except Exception as e:
            print("retry", a, e, flush=True)
            time.sleep(2 ** k)
    raise RuntimeError(f"fetch failed {var} {a}")


if __name__ == "__main__":
    var, t0, t1 = sys.argv[1], np.datetime64(sys.argv[2]), np.datetime64(sys.argv[3])
    t = ds().time.values
    i0 = int(np.searchsorted(t, t0)) // 8 * 8
    i1 = int(np.searchsorted(t, t1))
    d = os.path.join(WORK, var)
    os.makedirs(d, exist_ok=True)
    np.save(os.path.join(d, "times.npy"), t[i0:i1].astype("datetime64[h]"))
    jobs = [(var, a, min(a + 8 * BLK, i1), os.path.join(d, f"blk_{n:03d}.npy"))
            for n, a in enumerate(range(i0, i1, 8 * BLK))]
    # smaller units so 4 workers share the load
    jobs2 = []
    for var_, a, b, out in jobs:
        for k, s in enumerate(range(a, b, 8 * 25)):
            jobs2.append((var_, s, min(s + 8 * 25, b), out.replace(".npy", f"_{k:02d}.npy")))
    t00 = time.time()
    with mp.get_context("spawn").Pool(4) as p:
        for i, o in enumerate(p.imap_unordered(fetch, jobs2)):
            if i % 20 == 0:
                print(i, len(jobs2), f"{time.time() - t00:.0f}s", flush=True)
    print("done", time.time() - t00)
