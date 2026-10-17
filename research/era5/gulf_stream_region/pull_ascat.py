"""Download RSS ASCAT v02.1 daily bytemaps (Metop-A/B/C), October-April 2013-10..2026-04, crop to 25-45N, 85-50W,
keep (segment, time, wind, direction, rain) bytes. Output OUT/ascat_<sensor>_<winter>.npz (resumable). Total about 9.2 GB streamed.
Usage: python3 pull_ascat.py INVENTORY.csv OUT [WORKERS<=5]"""
import sys, os, gzip, time, urllib.request, concurrent.futures as cf
import numpy as np, pandas as pd
inv, OUT = sys.argv[1], sys.argv[2]; NW = min(5, int(sys.argv[3]) if len(sys.argv) > 3 else 5)
LAT0, LAT1, LON0, LON1 = 25.0, 45.0, 275.0, 310.0
LATS = -89.875 + 0.25 * np.arange(720); LONS = 0.125 + 0.25 * np.arange(1440)
I = np.nonzero((LATS >= LAT0) & (LATS <= LAT1))[0]; J = np.nonzero((LONS >= LON0) & (LONS <= LON1))[0]
def get(url, tries=20):
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r: return r.read()
        except urllib.error.HTTPError as e:
            if e.code == 404: return None
            err = e
        except Exception as e: err = e
        time.sleep(min(2 ** k, 30))
    raise RuntimeError(f"{url}: {err}")
d = pd.read_csv(inv, dtype={"date": str}); d = d[(d.date >= "20131001") & (d.date <= "20260430") & d.date.str[4:6].isin(["10", "11", "12", "01", "02", "03", "04"])]
d["winter"] = np.where(d.date.str[4:6].astype(int) >= 10, d.date.str[:4].astype(int), d.date.str[:4].astype(int) - 1)
os.makedirs(OUT, exist_ok=True)
def one(r):
    raw = get(r.url)
    if raw is None: return None
    a = np.frombuffer(gzip.decompress(raw), np.uint8)
    if a.size != 2 * 5 * 720 * 1440: raise ValueError(r.url)
    a = a.reshape(2, 5, 720, 1440)[:, :, I][:, :, :, J]
    return r.date, a[:, [0, 1, 2, 3]]          # time, wind, direction, rain
for (sen, w), g in d.groupby(["sensor", "winter"]):
    fn = f"{OUT}/ascat_{sen}_{w}.npz"
    if os.path.exists(fn): continue
    t0 = time.time()
    with cf.ThreadPoolExecutor(NW) as ex: res = [x for x in ex.map(one, list(g.itertuples())) if x is not None]
    np.savez_compressed(fn + ".tmp.npz", dates=np.array([x[0] for x in res]), data=np.stack([x[1] for x in res]), lat=LATS[I], lon=LONS[J])
    os.rename(fn + ".tmp.npz", fn); print(sen, w, len(res), "days", round(time.time() - t0), "s", flush=True)
