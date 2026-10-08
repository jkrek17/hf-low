"""Exact byte sizes of the objects each arm would stream (HEAD requests only, no field is downloaded).

Arm B: ARCO-ERA5 0.25 degree hourly, one object per time per surface variable.
Arm A: WeatherBench2 1.5 degree 6-hourly, one object per 8-time chunk per variable (all 13 levels in the object).
Run: python3 -I research/era5/explosive_onset/size_pull.py <repo_root>
"""
import sys, os, datetime, numpy as np, pandas as pd, urllib.request
from concurrent.futures import ThreadPoolExecutor
root = sys.argv[1]
T = pd.read_csv(os.path.join(root, "research/era5/explosive_onset/results/times_armB.csv"))
ARCO = "https://storage.googleapis.com/gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3"
WB2 = "https://storage.googleapis.com/weatherbench2/datasets/era5/1959-2023_01_10-6h-240x121_equiangular_with_poles_conservative.zarr"
def dt(t): s = str(t); return datetime.datetime(int(s[:4]), int(s[4:6]), int(s[6:8]), int(s[8:10]))
def head(url):
    for k in range(4):
        try:
            r = urllib.request.urlopen(urllib.request.Request(url, method="HEAD"), timeout=60)
            return int(r.headers["Content-Length"])
        except Exception as e:
            err = e
    raise err
def sizes(urls):
    with ThreadPoolExecutor(16) as ex: return list(ex.map(head, urls))
B_VARS = ["instantaneous_10m_wind_gust", "mean_sea_level_pressure", "10m_u_component_of_wind", "10m_v_component_of_wind"]
A_VARS = ["temperature", "geopotential", "vertical_velocity"]
L = []
tidx = [int((dt(t) - datetime.datetime(1900, 1, 1)).total_seconds() // 3600) for t in T.time]
T["tidx"] = tidx
res = {}
for v in B_VARS:
    res[v] = sizes([f"{ARCO}/{v}/{i}.0.0" for i in tidx])
for sel, name in [(T.kind == "case_time", "case times"), (T.kind == "noncase_sample", "sampled non-case times"), (T.kind != "x", "all")]:
    idx = np.where(sel)[0]
    tot = {v: sum(res[v][i] for i in idx) for v in B_VARS}
    L.append(f"ARM B ({name}, {len(idx)} times): " + ", ".join(f"{v} {tot[v]/1e9:.2f} GB" for v in B_VARS) + f"; total {sum(tot.values())/1e9:.2f} GB; per time {sum(tot.values())/len(idx)/1e6:.2f} MB")
gustmsl = sum(res[v][i] for v in B_VARS[:2] for i in range(len(T)))
L.append(f"ARM B reduced (gust + MSLP only, no u10/v10): {gustmsl/1e9:.2f} GB")
# Arm A: WB2 covers through 2023-01-09; later times would need ARCO pressure-level objects
WB2_END = datetime.datetime(2023, 1, 10)
w = T[[dt(t) < WB2_END for t in T.time]].copy()
late = T[[dt(t) >= WB2_END for t in T.time]]
w["chunk"] = [int((dt(t) - datetime.datetime(1959, 1, 1)).total_seconds() // 3600 // 6 // 8) for t in w.time]
ch = sorted(w.chunk.unique())
L.append(f"ARM A: times in WB2 range {len(w)} (case {int((w.kind=='case_time').sum())}, non-case {int((w.kind!='case_time').sum())}); times after 2023-01-09 {len(late)} (case {int((late.kind=='case_time').sum())}); distinct WB2 8-time chunks {len(ch)}")
tot = {}
for v in A_VARS:
    s = sizes([f"{WB2}/{v}/{c}.0.0.0" for c in ch]); tot[v] = sum(s)
    L.append(f"  {v}: {sum(s)/1e9:.2f} GB over {len(ch)} chunks ({np.mean(s)/1e6:.1f} MB per chunk)")
L.append(f"ARM A total (WB2, temperature + geopotential + vertical_velocity): {sum(tot.values())/1e9:.2f} GB")
wc = w[w.kind == "case_time"]; chc = sorted(wc.chunk.unique())
L.append(f"ARM A if only case-time chunks (no non-case sample): {len(chc)} chunks, {sum(tot.values())*len(chc)/len(ch)/1e9:.2f} GB (proportional estimate)")
a = sizes([f"{ARCO}/{v}/{tidx[0]}.0.0.0" for v in ["temperature", "geopotential", "vertical_velocity"]])
L.append(f"ARCO pressure-level object per time per variable (all 37 levels), one sample time: " + ", ".join(f"{x/1e6:.0f} MB" for x in a))
L.append(f"  covering the {len(late)} post-2023-01-09 times with those 3 variables: {len(late)*sum(a)/1e9:.1f} GB (not planned)")
L.append(f"TOTAL ARM B {sum(sum(res[v]) for v in B_VARS)/1e9:.2f} GB + ARM A {sum(tot.values())/1e9:.2f} GB = {(sum(sum(res[v]) for v in B_VARS)+sum(tot.values()))/1e9:.2f} GB")
open(os.path.join(root, "research/era5/explosive_onset/results/sizes.txt"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
