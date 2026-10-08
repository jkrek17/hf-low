"""Pull size by HEAD requests (no field downloaded). usage: size_pull.py REPO_ROOT [N_SAMPLE]
Storm scale: ARCO-ERA5 0.25 degree, one object per time per surface variable; large scale: sized separately in the README."""
import sys, os, datetime, numpy as np, pandas as pd, urllib.request
from concurrent.futures import ThreadPoolExecutor
root = sys.argv[1]; n = int(sys.argv[2]) if len(sys.argv) > 2 else 400
T = pd.read_csv(sys.argv[3] if len(sys.argv) > 3 else os.path.join(root, "research/era5/hf_vs_storm/results/times_storm.csv"))
B = "https://storage.googleapis.com/gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3"
VARS = ["instantaneous_10m_wind_gust", "mean_sea_level_pressure", "10m_u_component_of_wind", "10m_v_component_of_wind", "2m_dewpoint_temperature", "2m_temperature"]
if len(sys.argv) > 4:
    VARS = sys.argv[4].split(",")      # variable subset, e.g. a custom list of ARCO variable names
def ti(t):
    s = str(t); return int((datetime.datetime(int(s[:4]), int(s[4:6]), int(s[6:8]), int(s[8:10])) - datetime.datetime(1900, 1, 1)).total_seconds() // 3600)
def head(url):
    for k in range(4):
        try:
            return int(urllib.request.urlopen(urllib.request.Request(url, method="HEAD"), timeout=60).headers["Content-Length"])
        except Exception as e: err = e
    raise err
times = np.sort(T.time.unique()); rng = np.random.default_rng(1)
samp = rng.choice(times, min(n, len(times)), replace=False)
out = {}
with ThreadPoolExecutor(16) as ex:
    for v in VARS:
        out[v] = np.mean(list(ex.map(head, [f"{B}/{v}/{ti(t)}.0.0" for t in samp])))
tot = sum(out.values())
print(f"{len(times)} distinct times; sample {len(samp)}")
for v, m in out.items(): print(f"  {v}: {m/1e6:.2f} MB per time")
print(f"per time, 6 variables: {tot/1e6:.2f} MB; whole pull {tot*len(times)/1e9:.1f} GB; without t2m {(tot-out[VARS[-1]])*len(times)/1e9:.1f} GB")
