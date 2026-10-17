"""JPL MUR daily SST (analysed_sst, degC) from NOAA CoastWatch ERDDAP, 25-45N, 85-50W, stride 5 (0.05 degree), 1 Oct-30 Apr, winters 2013-14..2025-26.
Saved per winter as int16 (K-200)*100 (-32768 = missing), shape (days, 401, 701). Resumable per winter. Usage: python3 pull_mur.py OUT [THREADS]"""
import os, sys, io, time, datetime as dt, concurrent.futures as cf, requests, numpy as np
from scipy.io import netcdf_file
OUT = sys.argv[1]; NT = int(sys.argv[2]) if len(sys.argv) > 2 else 4
URL = "https://coastwatch.pfeg.noaa.gov/erddap/griddap/jplMURSST41.nc?analysed_sst%5B({d}T09:00:00Z):1:({d}T09:00:00Z)%5D%5B(25):5:(45)%5D%5B(-85):5:(-50)%5D"
S = requests.Session()
def day(d):
    for k in range(10):
        try:
            r = S.get(URL.format(d=d), timeout=120)
            if r.status_code == 404: return np.full((401, 701), -32768, np.int16)
            r.raise_for_status()
            with netcdf_file(io.BytesIO(r.content), mmap=False) as nc:
                v = nc.variables["analysed_sst"]; a = np.array(v[:]).reshape(401, 701).astype(float)
                sc = getattr(v, "scale_factor", 1.0); of = getattr(v, "add_offset", 0.0); fv = getattr(v, "_FillValue", None)
                raw = np.array(v[:]).reshape(401, 701)
            a = raw.astype(float) * sc + of
            if fv is not None: a[raw == fv] = np.nan
            a[a < -3] = np.nan
            return np.where(np.isfinite(a), np.clip(np.rint((a + 273.15 - 200) * 100), -32767, 32767), -32768).astype(np.int16)
        except Exception:
            time.sleep(min(2 ** k, 30))
    raise RuntimeError(d)
if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    for w in range(2013, 2026):
        fn = f"{OUT}/mur_{w}.npz"
        if os.path.exists(fn): continue
        ds = []; t = dt.date(w, 10, 1)
        while t < dt.date(w + 1, 5, 1): ds.append(t.isoformat()); t += dt.timedelta(days=1)
        with cf.ThreadPoolExecutor(NT) as ex: arr = list(ex.map(day, ds))
        np.savez_compressed(fn + ".tmp.npz", day=np.array([int(d.replace("-", "")) for d in ds]), sst=np.stack(arr))
        os.rename(fn + ".tmp.npz", fn); print(w, len(ds), flush=True)
