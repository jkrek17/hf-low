"""Stage 1: pull u10, v10, MSLP, 2 m temperature and SST (ERA5 PROXY, ARCO 0.25 degree) at the event and control times.

Event times = Nov-Mar 00/12 UTC with tip-box gust >= 54 kt (results/stage0_box_gust.csv); controls = 700 random Nov-Mar
times (month-stratified, seed 15) from the same list. Window 45-75N, 75W-5W is kept (float32, one npz per time) in
$ERA5_WORK/tip_jets/fields/. Resumable.
usage: stage1_pull.py probe | stage1_pull.py run [NPROC]
"""
import os, sys, time, urllib.request
import numpy as np, pandas as pd, numcodecs
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(os.environ.get("ERA5_WORK", os.path.join(HERE, "work")), "tip_jets")
OUT = os.path.join(WORK, "fields")
ARCO = "https://storage.googleapis.com/gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3"
VARS = {"u10": "10m_u_component_of_wind", "v10": "10m_v_component_of_wind", "msl": "mean_sea_level_pressure",
        "t2m": "2m_temperature", "sst": "sea_surface_temperature"}
T0 = np.datetime64("1900-01-01T00")
ROW = slice(60, 180)          # 75N .. 45.25N
COL = slice(1140, 1420)       # 285E (75W) .. 354.75E
_codec = numcodecs.Blosc()


def fetch(var, ti):
    url = f"{ARCO}/{var}/{ti}.0.0"
    for k in range(7):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                raw = r.read()
            return np.frombuffer(_codec.decode(raw), "<f4").reshape(721, 1440), len(raw)
        except Exception:
            if k == 6:
                raise
            time.sleep(2 ** k)


def times():
    d = pd.read_csv(os.path.join(HERE, "results", "stage0_box_gust.csv"))
    ev = d[d.box_gust_kt >= 54].copy(); ev["role"] = "event"
    d["month"] = pd.to_datetime(d.time).dt.month
    rng = np.random.default_rng(15)
    parts = []
    for m, g in d.groupby("month"):
        n = int(round(700 * len(g) / len(d)))
        parts.append(g.sample(n=n, random_state=int(rng.integers(1 << 30))))
    ct = pd.concat(parts); ct["role"] = "control"
    out = pd.concat([ev, ct]).drop_duplicates("time", keep="first")  # a control that is also an event stays an event
    return out[["time", "role"]].reset_index(drop=True)


def one(stamp):
    f = os.path.join(OUT, stamp + ".npz")
    if os.path.exists(f):
        return 0
    t = np.datetime64(stamp[:4] + "-" + stamp[4:6] + "-" + stamp[6:8] + "T" + stamp[8:10])
    ti = int((t - T0) / np.timedelta64(1, "h"))
    arr, n = {}, 0
    for k, v in VARS.items():
        a, b = fetch(v, ti); arr[k] = a[ROW, COL].copy(); n += b
    np.savez(f + ".tmp.npz", **arr); os.replace(f + ".tmp.npz", f)
    return n


def stamp(s):
    s = str(s); return s[:4] + s[5:7] + s[8:10] + s[11:13]


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    tt = times(); tt.to_csv(os.path.join(HERE, "results", "pull_times.csv"), index=False)
    st = [stamp(x) for x in tt.time]
    print(len(st), "times;", (tt.role == "event").sum(), "events")
    if sys.argv[1] == "probe":
        n = sum(one(s) for s in st[:3]); print("bytes per time (5 vars):", n / 3, "-> estimated GB:", n / 3 * len(st) / 1e9)
    else:
        t0 = time.time()
        with Pool(int(sys.argv[2]) if len(sys.argv) > 2 else 10) as p:
            tot = sum(p.imap_unordered(one, st, chunksize=4))
        print("streamed GB:", tot / 1e9, "seconds:", time.time() - t0)
