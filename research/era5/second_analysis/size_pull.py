"""Exact byte sizes (HEAD requests only) of the objects RA-27 would stream, for each sample floor.
ERA5: ARCO 0.25 deg hourly, 10m u and v, one object per time per variable.
HRES t0 analysis: WeatherBench2 hres_t0 0.25 deg 6-hourly, 10m u and v (and 10m_wind_speed for reference).
Run: python3 -I research/era5/second_analysis/size_pull.py <repo_root>
"""
import sys, os, datetime, pandas as pd
from concurrent.futures import ThreadPoolExecutor
import urllib.request
root = sys.argv[1]
ARCO = "https://storage.googleapis.com/gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3"
HRES = "https://storage.googleapis.com/weatherbench2/datasets/hres_t0/2016-2022-6h-1440x721.zarr"
def dt(t): s = str(t); return datetime.datetime(int(s[:4]), int(s[4:6]), int(s[6:8]), int(s[8:10]))
def head(url):
    for k in range(4):
        try: return int(urllib.request.urlopen(urllib.request.Request(url, method="HEAD"), timeout=60).headers["Content-Length"])
        except Exception as e: err = e
    raise err
def sizes(urls):
    with ThreadPoolExecutor(16) as ex: return list(ex.map(head, urls))
L = []
for floor in (55, 60, 65):
    S = pd.read_csv(os.path.join(root, f"research/era5/second_analysis/results/sample_fixes_{floor}.csv"))
    times = sorted(S.time.unique())
    ia = [int((dt(t) - datetime.datetime(1900, 1, 1)).total_seconds() // 3600) for t in times]
    ih = [int((dt(t) - datetime.datetime(2016, 1, 1)).total_seconds() // 3600 // 6) for t in times]
    tot = {}
    for v in ("10m_u_component_of_wind", "10m_v_component_of_wind"):
        tot["ERA5 " + v] = sum(sizes([f"{ARCO}/{v}/{i}.0.0" for i in ia]))
        tot["HRES " + v] = sum(sizes([f"{HRES}/{v}/{i}.0.0" for i in ih]))
    L.append(f"floor {floor} kt: {len(times)} times, {S.track.nunique()} tracks. " + "; ".join(f"{k} {x/1e9:.2f} GB" for k, x in tot.items()) + f"; TOTAL {sum(tot.values())/1e9:.2f} GB")
open(os.path.join(root, "research/era5/second_analysis/results/sizes.txt"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
