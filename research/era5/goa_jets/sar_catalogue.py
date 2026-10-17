"""SAR wind availability (asked by Jason 2026-10-17): how many Sentinel-1 wide-swath GRD scenes (IW_GRDH_1S or EW_GRDM_1S) overlap each Gulf of Alaska zone box
within +-6 h and on the same UTC day as each ERA5-proxy 54 kt episode time, 2014 on. Catalogue COUNTS only (Copernicus Data Space OData); nothing is downloaded.
Writes results/sar_catalogue.csv and sar_catalogue.txt. Zone boxes from hrrr_analyse.ZONES. Control rows of episodes.csv (is_ctrl) are left out of the summary. usage: sar_catalogue.py [--reuse] (reuse: summarise the saved csv without querying)"""
import os, sys, time, json, urllib.request, urllib.parse
from concurrent.futures import ThreadPoolExecutor
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from hrrr_analyse import ZONES
URL = "https://catalogue.dataspace.copernicus.eu/odata/v1/Products"
def count(zn, t0, t1):
    a, b = (count1(zn, t0, t1, p) for p in ("IW_GRDH_1S", "EW_GRDM_1S"))
    return None if a is None or b is None else a + b
def count1(zn, t0, t1, pt):
    s, n, w, e = ZONES[zn]; poly = f"POLYGON(({w} {s},{e} {s},{e} {n},{w} {n},{w} {s}))"
    f = (f"Collection/Name eq 'SENTINEL-1' and ContentDate/Start ge {t0:%Y-%m-%dT%H:%M:%S}.000Z and ContentDate/Start lt {t1:%Y-%m-%dT%H:%M:%S}.000Z "
         f"and OData.CSC.Intersects(area=geography'SRID=4326;{poly}') and Attributes/OData.CSC.StringAttribute/any(a:a/Name eq 'productType' and a/OData.CSC.StringAttribute/Value eq '{pt}')")
    q = URL + "?" + urllib.parse.urlencode({"$filter": f, "$count": "true", "$top": "1"}, quote_via=urllib.parse.quote)
    for k in range(4):
        try: return json.load(urllib.request.urlopen(q, timeout=60))["@odata.count"]
        except Exception as ex: time.sleep(2 * (k + 1)); err = ex
    return None
E = pd.read_csv(os.path.join(HERE, "results", "episodes.csv"), parse_dates=["time"])
E = E[(E.time >= "2014-10-01")].copy()
def job(r):
    t = r.time; d0 = t.normalize()
    return (count(r.zone, t - pd.Timedelta(hours=6), t + pd.Timedelta(hours=6)), count(r.zone, d0, d0 + pd.Timedelta(days=1)))
CSV = os.path.join(HERE, "results", "sar_catalogue.csv")
if "--reuse" in sys.argv:
    E = E[["time", "zone", "is_ctrl"]].merge(pd.read_csv(CSV, parse_dates=["time"]), on=["time", "zone"])
else:
    with ThreadPoolExecutor(4) as ex: res = list(ex.map(job, E.itertuples()))
    E["scenes_pm6h"] = [a for a, b in res]; E["scenes_day"] = [b for a, b in res]
E[["time", "zone", "is_ctrl", "gust_kt", "regime", "strong_t", "scenes_pm6h", "scenes_day"]].to_csv(CSV, index=False)
E = E[~E.is_ctrl.astype(bool)]
out = [f"episodes from 2014-10 on: {len(E)}; queries failed: {int(E.scenes_pm6h.isna().sum() + E.scenes_day.isna().sum())}"]
for zn in ZONES:
    x = E[E.zone == zn].dropna(subset=["scenes_pm6h", "scenes_day"])
    out.append(f"{zn}: episodes {len(x)}; with >=1 scene within +-6 h {int((x.scenes_pm6h >= 1).sum())} ({(x.scenes_pm6h >= 1).mean():.2f}); same UTC day {int((x.scenes_day >= 1).sum())} ({(x.scenes_day >= 1).mean():.2f}); "
               f"strong (64 kt) episodes {int(x.strong_t.sum())}, with a scene within +-6 h {int(((x.scenes_pm6h >= 1) & x.strong_t).sum())}")
open(os.path.join(HERE, "results", "sar_catalogue.txt"), "w").write("\n".join(out) + "\n"); print("\n".join(out))
