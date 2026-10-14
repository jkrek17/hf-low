"""Step 2: per cycle, forecast day and basin, count 1-degree ocean point-days by GFS gust band and by OPC outcome.
Predictor: maximum of the day's daily-max gust (kt) over ocean cells within R_d nm of the point (R_d = 60, 120, 270, 270 nm).
Truth: an OPC archive fix of Category HF (HF-agreed) / S or HF (storm-agreed) within 324 nm of the point at one of the day's four valid times.
usage: tabulate.py DMAX_DIR OPC_DIR OUT.csv [NPROC]   -> one row per cycle, day, basin, band; counts only (small)."""
import os, sys, glob
import numpy as np, pandas as pd
from multiprocessing import Pool
RD = [60, 120, 270, 270]; BANDS = [0, 40, 48, 56, 64, 72, 80, 999]; TRUTH_NM = 324.0
dm, opc_dir, out = sys.argv[1:4]; NPROC = int(sys.argv[4]) if len(sys.argv) > 4 else 4
LAT = 74.75 - 0.5 * np.arange(110); LON = 0.125 + 0.5 * np.arange(720)
LSM = np.load("/mnt/project-files/era5-cache/lsm_nh_0p25.npy")[60:280] < 0.5
OC = LSM.reshape(110, 2, 720, 2).sum(axis=(1, 3)) >= 3                      # ocean 0.5 cell: at least 3 of 4 ocean
ER, EC = np.meshgrid(np.arange(0, 110, 2), np.arange(0, 720, 2), indexing="ij")   # 1 degree evaluation points
BAS = np.where((EC * 0.5 + 0.125 >= 260) | (EC * 0.5 + 0.125 <= 10), 1, np.where((EC * 0.5 + 0.125 >= 120) & (EC * 0.5 + 0.125 <= 240), 2, 0))
OK = OC[ER, EC] & (BAS > 0)                                                  # 1 = Atlantic (100W-10E), 2 = Pacific (120E-120W)
PL, PN = LAT[ER], LON[EC]
parts = []
for tag in ("Atl", "Pac"):
    d = pd.read_csv(os.path.join(opc_dir, f"HF_Data_-_{tag}.csv"), dtype={"date": str})
    d = d[d.date.str.fullmatch(r"\d{10}", na=False)].copy(); d["t"] = pd.to_datetime(d.date, format="%Y%m%d%H", errors="coerce"); parts.append(d.dropna(subset=["t"]))
OPC = pd.concat(parts, ignore_index=True); OPC = OPC[OPC.Category.isin(["HF", "S"])].copy()
OPC["lon"] = OPC.Longitude % 360; OPC["lat"] = OPC.Latitude

_OFF = {}
def offsets(R, row):
    key = (R, row // 10)
    if key not in _OFF:
        lat = 74.75 - 0.5 * (row // 10 * 10 + 5); dx = 30.0 * np.cos(np.radians(lat))
        ky = int(R // 30); kx = int(R // dx)
        yy, xx = np.meshgrid(np.arange(-ky, ky + 1), np.arange(-kx, kx + 1), indexing="ij")
        m = ((yy * 30.0) ** 2 + (xx * dx) ** 2) <= R ** 2
        _OFF[key] = (yy[m], xx[m], ky, kx)
    return _OFF[key]

def maxwithin(G, R):
    """G (110,720) float; returns (55,360) maximum within R nm (circular in longitude, edges clipped) at the evaluation points."""
    out = np.zeros(ER.shape, np.float32); PADY, PADX = 12, 120
    P = np.pad(G, ((PADY, PADY), (0, 0))); P = np.pad(P, ((0, 0), (PADX, PADX)), mode="wrap")
    for b in range(0, 110, 10):
        sel = (ER >= b) & (ER < b + 10)
        if not sel.any(): continue
        dy, dx, ky, kx = offsets(R, b)
        r = ER[sel] + PADY; c = EC[sel] + PADX
        out[sel] = P[r[:, None] + dy[None, :], c[:, None] + dx[None, :]].max(axis=1)
    return out

def hav(la1, lo1, la2, lo2):
    la1, lo1, la2, lo2 = map(np.radians, (la1, lo1, la2, lo2))
    a = np.sin((la2 - la1) / 2) ** 2 + np.cos(la1) * np.cos(la2) * np.sin((lo2 - lo1) / 2) ** 2
    return 2 * 3440.065 * np.arcsin(np.sqrt(np.clip(a, 0, 1)))        # nautical miles

def do_cycle(path):
    day = os.path.basename(path)[:8]; t0 = pd.Timestamp(day)
    d = np.load(path)["d"]; rows = []
    for k in range(4):
        if (d[k] == 255).any(): continue
        G = d[k].astype(np.float32) / 2
        pred = maxwithin(G, RD[k])
        lo, hi = t0 + pd.Timedelta(hours=24 * k + 6), t0 + pd.Timedelta(hours=24 * k + 24)
        fx = OPC[(OPC.t >= lo) & (OPC.t <= hi)]
        hf = np.zeros(ER.shape, bool); st = np.zeros(ER.shape, bool)
        for cat, arr in (("HF", hf), ("S", st)):
            for la, lo_ in zip(fx[fx.Category == cat].lat, fx[fx.Category == cat].lon):
                arr |= hav(PL, PN, la, lo_) <= TRUTH_NM
        st |= hf
        band = np.digitize(pred, BANDS) - 1
        for bs in (1, 2):
            m = OK & (BAS == bs)
            for b in range(len(BANDS) - 1):
                s = m & (band == b); n = int(s.sum())
                if n: rows.append((day, k + 1, bs, b, n, int((s & hf).sum()), int((s & st).sum())))
    return rows

if __name__ == "__main__":
    files = sorted(glob.glob(os.path.join(dm, "cycles", "*.npz")))
    with Pool(NPROC) as p: res = p.map(do_cycle, files, chunksize=4)
    df = pd.DataFrame([r for rr in res for r in rr], columns=["cycle", "day", "basin", "band", "n", "hf", "storm"])
    df["basin"] = df.basin.map({1: "atl", 2: "pac"}); df.to_csv(out, index=False); print(len(files), "cycles", len(df), "rows")
