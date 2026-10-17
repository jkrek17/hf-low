"""Pull ERA5 mean sea level pressure (ARCO-ERA5, 0.25 deg) at 00 and 12 UTC, 1 Oct-30 Apr, winters 2004-05..2024-25,
cropped to 20-55N, 90-50W, stored as uint16 (Pa - 85000). Resumable per winter. Streams ~19.6 GB (global field per time).
Usage: python3 pull_msl.py OUTDIR [NPROC]"""
import os, sys, datetime as dt, numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import era5lib as L
OUT = sys.argv[1]; NP = int(sys.argv[2]) if len(sys.argv) > 2 else 12
R0, R1, C0, C1 = 140, 281, 1080, 1241   # 55N..20N, 270E..310E (90W..50W)
def times(w):
    t = dt.datetime(w, 10, 1); end = dt.datetime(w + 1, 5, 1); out = []
    while t < end:
        out += [t, t + dt.timedelta(hours=12)]; t += dt.timedelta(days=1)
    return out
def one(t):
    a = L.field("mean_sea_level_pressure", L.tidx(t.strftime("%Y%m%d%H")))[R0:R1, C0:C1]
    return np.clip(np.rint(a - 85000), 0, 65535).astype(np.uint16)
def winter(w):
    f = f"{OUT}/msl_{w}.npz"
    if os.path.exists(f): return w, 0
    ts = times(w); arr = np.stack([one(t) for t in ts])
    np.savez(f + ".tmp.npz", t=np.array([int(t.strftime("%Y%m%d%H")) for t in ts]), msl=arr)
    os.rename(f + ".tmp.npz", f); return w, L.BYTES[0]
if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    with Pool(NP) as p:
        for w, b in p.imap_unordered(winter, range(2004, 2025)): print(w, b / 1e9, "GB streamed (this proc)", flush=True)
