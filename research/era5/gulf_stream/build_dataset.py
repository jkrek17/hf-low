"""Build the station-time table for the Gulf Stream study (00 and 12 UTC, Oct-Apr, winters 2004-05..2024-25).
NDBC buoy report nearest to the synoptic hour (within 35 min) + ERA5 MSLP-derived pressure-gradient wind Vg.
Usage: python3 build_dataset.py NDBC_DIR MSL_DIR OUT.csv.gz   (Vg: ERA5 MSLP proxy, a weather control, not an outcome)"""
import sys, glob, gzip, io, os, numpy as np, pandas as pd
NDBC, MSL, OUT = sys.argv[1:4]
ST = {"41002": (33.989, -75.212), "41001": (34.738, -72.503), "41025": (35.026, -75.380), "41013": (33.436, -77.764),
      "44014": (36.603, -74.837), "41004": (32.502, -79.099), "44009": (38.460, -74.692)}
H = dict(l.strip().split(",") for l in open(f"{NDBC}/anemometer_heights.csv").read().strip().split("\n")[1:])
KT = 1.943844
def load(s):
    fr = []
    for f in sorted(glob.glob(f"{NDBC}/{s}_*.txt.gz")):
        txt = gzip.open(f, "rt").read().splitlines()
        txt = [l for l in txt if l.strip() and l.split()[0].isdigit()]
        t = pd.read_csv(io.StringIO("\n".join(txt)), sep=r"\s+", header=None, na_values=[99, 99.0, 999, 9999, 99.00, 999.0, 9999.0])
        if t.shape[1] < 15: continue
        t = t.iloc[:, :15]; t.columns = "yy mm dd hh mn wdir wspd gst wvht dpd apd mwd pres atmp wtmp".split()
        t.loc[t.yy < 100, "yy"] += 2000
        fr.append(t)
    d = pd.concat(fr)
    d["time"] = pd.to_datetime(dict(year=d.yy, month=d.mm, day=d.dd, hour=d.hh, minute=d.mn))
    # sentinels: some years use 99/999 per column; also drop impossible values
    for c, lim in [("wspd", 60), ("gst", 80), ("wvht", 30), ("dpd", 40), ("atmp", 50), ("wtmp", 50)]:
        d.loc[d[c] >= lim, c] = np.nan
    d.loc[d.wdir > 360, "wdir"] = np.nan
    return d.sort_values("time").reset_index(drop=True)
def vg(msl, la, lo, lat0, lon0):
    # msl: (141,161) Pa for 55N..20N, 270E..310E at 0.25 deg; central differences +-0.5 deg (2 cells)
    i = int(round((55 - lat0) / 0.25)); j = int(round((lon0 % 360 - 270) / 0.25))
    dy = 2 * 0.5 * 111_000; dx = dy * np.cos(np.radians(lat0))
    dpdx = (msl[i, j + 2] - msl[i, j - 2]) / dx; dpdy = (msl[i - 2, j] - msl[i + 2, j]) / dy
    f = 2 * 7.292e-5 * np.sin(np.radians(lat0))
    return np.hypot(dpdx, dpdy) / (1.25 * f) * KT
rows = []
msl = {}
for f in sorted(glob.glob(f"{MSL}/msl_*.npz")):
    z = np.load(f); w = int(os.path.basename(f)[4:8])
    for t, a in zip(z["t"], z["msl"]): msl[int(t)] = a.astype(np.float32) + 85000
for s, (la, lo) in ST.items():
    d = load(s); z = float(H[s]) if H.get(s) else 5.0
    d = d.dropna(subset=["time"]).drop_duplicates("time").set_index("time").sort_index()
    ts = sorted(msl)
    for t in ts:
        T = pd.to_datetime(str(t), format="%Y%m%d%H"); w = T.year if T.month >= 10 else T.year - 1
        sub = d.loc[T - pd.Timedelta(minutes=35): T + pd.Timedelta(minutes=35)]
        if sub.empty or np.isnan(sub.wspd.iloc[0]): continue
        r = sub.iloc[(sub.index - T).to_series().abs().argmin()]
        prev = []
        for k in (1, 2, 3):
            q = d.loc[T - pd.Timedelta(hours=k, minutes=35): T - pd.Timedelta(hours=k) + pd.Timedelta(minutes=35)]
            prev.append((q.wspd.iloc[0], q.wdir.iloc[0]) if len(q) else (np.nan, np.nan))
        rows.append(dict(station=s, time=T, winter=w, month=T.month, wdir=r.wdir, wspd_raw=r.wspd * KT, gst=r.gst * KT,
                         v10=r.wspd * KT * (10 / z) ** 0.11, wvht_ft=r.wvht * 3.28084, dpd=r.dpd, atmp=r.atmp, wtmp=r.wtmp,
                         vg=vg(msl[t], None, None, la, lo),
                         p_ws=[p[0] * KT for p in prev], p_wd=[p[1] for p in prev]))
df = pd.DataFrame(rows)
for k in range(3):
    df[f"ws_m{k+1}"] = df.p_ws.apply(lambda x: x[k]); df[f"wd_m{k+1}"] = df.p_wd.apply(lambda x: x[k])
df = df.drop(columns=["p_ws", "p_wd"])
df["dT"] = (df.wtmp - df.atmp).clip(-3, 12)
df.to_csv(OUT, index=False)
print(len(df)); print(df.groupby("station").agg(n=("v10", "size"), winters=("winter", "nunique"), wvht=("wvht_ft", "count"), dT=("dT", "count")))
