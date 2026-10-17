"""ERA5 2 m dewpoint temperature (00 and 12 UTC), 1 Oct-30 Apr, winters 2013-14..2024-25, cropped as pull_era5_temp.py (25-45N, 85-50W), int16 (K-200)*100.
Resumable per winter. Streams about 12 GB (global field per time). Usage: python3 pull_era5_d2m.py OUT [NPROC]"""
import os, sys, datetime as dt, numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import era5lib as L
OUT = sys.argv[1]; NP = int(sys.argv[2]) if len(sys.argv) > 2 else 4
R0, R1, C0, C1 = 180, 261, 1100, 1241
def enc(a): return np.clip(np.rint((a - 200.0) * 100), -32768, 32767).astype(np.int16)
def winter(w):
    f = f"{OUT}/era5_d2m_{w}.npz"
    if os.path.exists(f): return w, 0
    t = dt.datetime(w, 10, 1); end = dt.datetime(w + 1, 5, 1); d2, tt = [], []
    while t < end:
        for h in (0, 12):
            tm = t + dt.timedelta(hours=h); i = L.tidx(tm.strftime("%Y%m%d%H"))
            d2.append(enc(np.nan_to_num(L.field("2m_dewpoint_temperature", i)[R0:R1, C0:C1], nan=0.0))); tt.append(int(tm.strftime("%Y%m%d%H")))
        t += dt.timedelta(days=1)
    np.savez_compressed(f + ".tmp.npz", d2m_time=np.array(tt), d2m=np.stack(d2))
    os.rename(f + ".tmp.npz", f); return w, L.BYTES[0]
if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    with Pool(NP) as p:
        for w, b in p.imap_unordered(winter, range(2013, 2025)): print(w, round(b / 1e9, 2), "GB streamed (this proc)", flush=True)
