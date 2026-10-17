"""NOAA absolute geostrophic surface currents (ugos, vgos; m/s) daily, from CoastWatch ERDDAP nesdisSSH1day, 25-45N, 85-50W (0.25 degree), 1 Oct-30 Apr, winters 2017-18..2024-25.
Saved per winter as int16 (m/s*1000, -32768 missing), shape (2, 80, 140) per day (u, v). Usage: python3 pull_currents.py OUT [THREADS]"""
import os, sys, io, time, datetime as dt, concurrent.futures as cf, requests, numpy as np
from scipy.io import netcdf_file
OUT = sys.argv[1]; NT = int(sys.argv[2]) if len(sys.argv) > 2 else 4
URL = "https://coastwatch.pfeg.noaa.gov/erddap/griddap/nesdisSSH1day.nc?ugos%5B({d}T00:00:00Z):1:({d}T00:00:00Z)%5D%5B(25.125):1:(44.875)%5D%5B(-84.875):1:(-50.125)%5D,vgos%5B({d}T00:00:00Z):1:({d}T00:00:00Z)%5D%5B(25.125):1:(44.875)%5D%5B(-84.875):1:(-50.125)%5D"
S = requests.Session()
def rd(nc, n):
    v = nc.variables[n]; raw = np.array(v[:]); a = raw.astype(float) * getattr(v, "scale_factor", 1.0) + getattr(v, "add_offset", 0.0)
    fv = getattr(v, "_FillValue", None)
    if fv is not None: a[raw == fv] = np.nan
    return a.reshape(a.shape[-2], a.shape[-1])
def day(d):
    for k in range(10):
        try:
            r = S.get(URL.format(d=d), timeout=120)
            if r.status_code == 404: return None
            r.raise_for_status()
            with netcdf_file(io.BytesIO(r.content), mmap=False) as nc: u, v = rd(nc, "ugos"), rd(nc, "vgos")
            e = lambda a: np.where(np.isfinite(a), np.clip(np.rint(a * 1000), -32767, 32767), -32768).astype(np.int16)
            return np.stack([e(u), e(v)])
        except Exception:
            time.sleep(min(2 ** k, 30))
    return None
if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    for w in range(2017, 2025):
        fn = f"{OUT}/cur_{w}.npz"
        if os.path.exists(fn): continue
        ds = []; t = dt.date(w, 10, 1)
        while t < dt.date(w + 1, 5, 1): ds.append(t.isoformat()); t += dt.timedelta(days=1)
        with cf.ThreadPoolExecutor(NT) as ex: res = list(ex.map(day, ds))
        ok = [(d, a) for d, a in zip(ds, res) if a is not None]
        np.savez_compressed(fn + ".tmp.npz", day=np.array([int(d.replace("-", "")) for d, _ in ok]), uv=np.stack([a for _, a in ok]))
        os.rename(fn + ".tmp.npz", fn); print(w, len(ok), "of", len(ds), flush=True)
