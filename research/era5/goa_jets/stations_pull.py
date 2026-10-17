"""Coastal-station test (PREREGISTRATION_stations.md): NDBC station gusts, HRRR f03 gust and ERA5 gust at the station points, at every 00/12 UTC time Oct-Apr 2018-19..2025-26.
modes: obs   - download NDBC historical stdmet files, keep GST (kt) in $ERA5_WORK/goa_jets/stations/obs.csv
       hrrr  - GUST record only of each HRRR f03 file (about 0.78 MB per time), nearest cell and 3x3 max at each station -> work/hrrr_st.csv (resumable)
       era5  - cached ERA5 gust at the nearest 0.25 degree cell -> work/era5_st.csv
       build - join all three into results/stations_pairs.csv (committed, small)
usage: stations_pull.py obs|hrrr [NPROC]|era5|build"""
import os, sys, gzip, urllib.request, time
import numpy as np, pandas as pd
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import hrrr_pull as HP, hrrr_grid as HG
W = os.path.join(os.environ.get("ERA5_WORK", os.path.join(HERE, "work")), "goa_jets", "stations"); os.makedirs(W, exist_ok=True)
ST = {"pila2": (59.742, -149.47, "Z1 Z3"), "mdxa2": (59.438, -146.327, "Z1"), "cdea2": (56.001, -134.136, "Z2"), "lcna2": (56.058, -132.69, "Z2"), "gixa2": (55.446, -131.881, "Z2"),
      "amaa2": (58.916, -151.952, "Z3"), "fila2": (59.331, -151.995, "Z3"), "auga2": (59.378, -153.356, "Z3"), "cspa2": (58.199, -136.64, "")}
KT = 1.9438445
LA, LO, _ = HG.grid()


def obs():
    rows = []
    for s in ST:
        for y in range(2018, 2027):
            try: t = gzip.decompress(urllib.request.urlopen(f"https://www.ndbc.noaa.gov/data/historical/stdmet/{s}h{y}.txt.gz", timeout=60).read()).decode().split("\n")
            except Exception: continue
            h = t[0].lstrip("#").split()
            if "GST" not in h: continue
            gi, mi = h.index("GST"), h.index("mm") if "mm" in h else None
            for l in t[2:]:
                f = l.split()
                if len(f) <= gi: continue
                try:
                    g = float(f[gi]); yy = int(f[0]); mn = int(f[mi]) if mi else 0
                    ts = pd.Timestamp(year=yy if yy > 100 else 2000 + yy, month=int(f[1]), day=int(f[2]), hour=int(f[3]), minute=mn)
                except ValueError: continue
                if g < 90: rows.append((s, ts, round(g * KT, 1)))
    pd.DataFrame(rows, columns=["station", "time", "gust_kt"]).to_csv(os.path.join(W, "obs.csv"), index=False); print(len(rows), "obs")


_ij = {}
def cells():
    if not _ij:
        for s, (la, lo, _) in ST.items():
            j, i = HG.ij(la, lo); _ij[s] = (int(round(float(j))), int(round(float(i))))
    return _ij


def gust_one(valid):
    u = HP.url(valid)
    try: idx = HP.get(u + ".idx")
    except HP.Missing: return valid, None
    rows = []
    for line in idx.decode().splitlines():
        p = line.split(":")
        if len(p) >= 6 and p[1].isdigit(): rows.append((int(p[1]), p[3], p[4]))
    for k, (off, v, l) in enumerate(rows):
        if (v, l) == ("GUST", "surface"):
            end = rows[k + 1][0] - 1 if k + 1 < len(rows) else ""
            try: blob = HP.get(u, f"bytes={off}-{end}"); a = HP.G.decode_message(blob)
            except HP.Missing: return valid, None
            except Exception: a = HP.G.decode_message(HP.walk(u, k))
            out = {}
            for s, (j, i) in cells().items():
                out[s + "_nn"] = float(a[j, i]) * KT; out[s + "_max3"] = float(a[j - 1:j + 2, i - 1:i + 2].max()) * KT
            return valid, out
    return valid, None


def hrrr(n):
    T = HP.vtimes(); fn = os.path.join(W, "hrrr_st.csv"); done = set()
    if os.path.exists(fn): done = set(pd.read_csv(fn, parse_dates=["time"]).time)
    T = [t for t in T if t not in done]; print(len(T), "to do", flush=True); t0 = time.time()
    with Pool(n) as p:
        for k, (v, o) in enumerate(p.imap_unordered(gust_one, T)):
            if o is None: continue
            d = pd.DataFrame([{"time": v, **o}]); d.to_csv(fn, mode="a", header=not os.path.exists(fn), index=False)
            if k % 200 == 0: print(k, f"{time.time() - t0:.0f}s", flush=True)


def era5():
    import numcodecs
    out = []; A = "/mnt/project-files/era5-cache/gust_nh_12h/"
    for f in sorted(os.listdir(A)):
        z = np.load(A + f); a = np.frombuffer(numcodecs.Blosc().decode(z["data"].tobytes()), np.uint8).reshape(z["shape"])
        t = pd.Timestamp("1900-01-01") + pd.to_timedelta(z["hours_since_1900"], unit="h")
        for k, tt in enumerate(t):
            r = {"time": tt}
            for s, (la, lo, _) in ST.items(): r[s] = a[k, int(round((90 - la) / 0.25)), int(round(((lo + 360) % 360) / 0.25))] / 2.0
            out.append(r)
    pd.DataFrame(out).to_csv(os.path.join(W, "era5_st.csv"), index=False)



def build():
    o = pd.read_csv(os.path.join(W, "obs.csv"), parse_dates=["time"]); h = pd.read_csv(os.path.join(W, "hrrr_st.csv")); h["time"] = pd.to_datetime(h.time, format="mixed")
    e = pd.read_csv(os.path.join(W, "era5_st.csv"), parse_dates=["time"])
    h = h.melt(id_vars="time"); h["station"] = h.variable.str.split("_").str[0]; h["kind"] = h.variable.str.split("_").str[1]
    h = h.pivot_table(index=["time", "station"], columns="kind", values="value").reset_index().rename(columns={"nn": "hrrr_nn_kt", "max3": "hrrr_max3_kt"})
    e = e.melt(id_vars="time", var_name="station", value_name="era5_kt")
    p = h.merge(e, on=["time", "station"])
    # largest GST within 30 minutes of the hour: bin each obs to its nearest hour, keep the max
    o["otime"] = o.time; o["time"] = o.time.dt.round("h"); o = o[(o.otime - o.time).abs() <= pd.Timedelta(minutes=30)]
    o = o.groupby(["time", "station"], as_index=False).gust_kt.max().rename(columns={"gust_kt": "obs_kt"})
    p = p.merge(o, on=["time", "station"]).sort_values(["station", "time"])
    p[["time", "station", "obs_kt", "hrrr_nn_kt", "hrrr_max3_kt", "era5_kt"]].round(1).to_csv(os.path.join(HERE, "results", "stations_pairs.csv"), index=False); print(len(p), "pairs")


if __name__ == "__main__":
    m = sys.argv[1]
    if m == "obs": obs()
    elif m == "hrrr": hrrr(int(sys.argv[2]) if len(sys.argv) > 2 else 8)
    elif m == "era5": era5()
    elif m == "build": build()
