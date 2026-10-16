"""Stage 1: pull u10, v10, MSLP, 2 m temperature and SST (ERA5 PROXY, ARCO 0.25 degree) at event and control times.

Event times = October-April 00/12 UTC times with gust >= 54 kt in ANY zone (results/stage0_zone_gust.csv). Controls = 600 times drawn
uniformly at random (seed 16) from ALL 9,338 October-April times, kept whether or not they are events. Window 44-70N, 175-125W is kept
(float32, one npz per time) in $ERA5_WORK/goa_jets/fields/. Resumable.
usage: stage1_pull.py probe | stage1_pull.py run [NPROC]
"""
import os, sys, time, urllib.request
import numpy as np, pandas as pd, numcodecs
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(os.environ.get("ERA5_WORK", os.path.join(HERE, "work")), "goa_jets")
OUT = os.path.join(WORK, "fields")
ARCO = "https://storage.googleapis.com/gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3"
VARS = {"u10": "10m_u_component_of_wind", "v10": "10m_v_component_of_wind", "msl": "mean_sea_level_pressure",
        "t2m": "2m_temperature", "sst": "sea_surface_temperature"}
T0 = np.datetime64("1900-01-01T00")
ROW = slice(80, 185)     # 70N .. 43.75N
COL = slice(740, 941)    # 185E (175W) .. 235E (125W)
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
    d = pd.read_csv(os.path.join(HERE, "results", "stage0_zone_gust.csv"))
    g = d[[c for c in d.columns if c.endswith("gust_kt")]].max(axis=1)
    d["event"] = g >= 54
    ctrl = d.sample(n=600, random_state=16)
    ev = d[d.event].copy(); ev["role"] = "event"
    ct = ctrl.copy(); ct["role"] = "control"
    allt = pd.concat([ev, ct])
    out = pd.DataFrame({"time": d.time})
    out["event"] = d.event.values
    out["is_control"] = d.time.isin(ctrl.time).values
    return out[out.event | out.is_control].reset_index(drop=True)


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
    print(len(st), "times;", int(tt.event.sum()), "events;", int(tt.is_control.sum()), "controls;", int((tt.event & tt.is_control).sum()), "both")
    if sys.argv[1] == "probe":
        n = sum(one(s) for s in st[:3]); print("bytes per time (5 vars):", n / 3, "-> estimated GB:", n / 3 * len(st) / 1e9)
    else:
        t0 = time.time()
        with Pool(int(sys.argv[2]) if len(sys.argv) > 2 else 10) as p:
            tot = sum(p.imap_unordered(one, st, chunksize=4))
        print("streamed GB:", tot / 1e9, "seconds:", time.time() - t0)
