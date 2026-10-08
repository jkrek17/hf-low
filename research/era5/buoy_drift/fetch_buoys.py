"""Fetch moored-buoy wind reports from NOAA's Integrated Surface Database (ISD).

Why ISD: the NDBC, MEDS, ICOADS and CDS hosts are not reachable from the work
containers this project runs in, but NOAA's ISD "global-hourly" copy on AWS
(s3://noaa-global-hourly-pds) is. ISD carries the GTS reports (FM-18 BUOY and
FM-13 SHIP) of NDBC and Canadian moored buoys under USAF ids 99xxxx, from the
late 1970s to the end of 2004. ISD holds no moored-buoy reports after 2004, so
this record cannot reach the 2004-onward archive period except for its first
months; see README.md.

What is kept, per report: time (UTC), wind direction, wind speed (m/s, as
reported: the buoy's averaging-period mean, 8 min for NDBC), the ISD quality
flags, the report position, and the supplementary-wind (OA1-OA3) groups as
raw strings so the gust-code change can be documented. Nothing is adjusted
here; height adjustment and screening are in analyse.py.

    python3 research/era5/buoy_drift/fetch_buoys.py      # -> work/isd/<usaf>.csv.gz
Resumable: a station whose output exists is skipped. About 1.5 GB is streamed
and about 30 MB is kept (in work/, which is ignored).
"""
import concurrent.futures as cf
import csv
import gzip
import io
import os
import urllib.error
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "work", "isd")
URL = "https://noaa-global-hourly-pds.s3.amazonaws.com/{y}/{sid}99999.csv"
YEARS = range(1978, 2006)

# ISD USAF id -> (label, basin). Offshore moored buoys in the two archive
# basins with ISD reports in both the 1980s/early 1990s and 2000-2004. Labels
# are the WMO ids of the NDBC buoys at those positions; the Canadian (MEDS)
# buoys carry ISD's station name, since their WMO ids were not checkable here.
# The ocean screen (land-sea mask) is applied in analyse.py, not here.
STATIONS = {
    "992030": ("46001", "Pac"), "992160": ("46002", "Pac"), "992170": ("46003", "Pac"),
    "992190": ("MidNomad", "Pac"), "992210": ("46005", "Pac"), "992200": ("46006", "Pac"),
    "992820": ("46035", "Pac"), "992930": ("SNomad", "Pac"), "996210": ("NNomad", "Pac"),
    "992070": ("44004", "Atl"), "992150": ("41002", "Atl"), "992360": ("44011", "Atl"),
    "992790": ("44008", "Atl"), "992760": ("44005", "Atl"), "992400": ("44009", "Atl"),
    "992420": ("44013", "Atl"), "992780": ("44007", "Atl"), "996300": ("44014", "Atl"),
    "996420": ("44025", "Atl"), "995480": ("EScotia", "Atl"), "995490": ("SWGrBank", "Atl"),
    "995500": ("Banquero", "Atl"), "996250": ("LaHave", "Atl"), "996260": ("LaurFan", "Atl"),
    "992980": ("TailBank", "Atl"), "995340": ("Nickersn", "Atl"),
}
KEEP = ["DATE", "REPORT_TYPE", "LATITUDE", "LONGITUDE", "WND", "OA1", "OA2", "OA3"]


def fetch(sid):
    path = os.path.join(OUT, sid + ".csv.gz")
    if os.path.exists(path):
        return sid, "cached"
    rows = 0
    tmp = path + ".part"
    with gzip.open(tmp, "wt", newline="") as fo:
        w = csv.writer(fo)
        w.writerow(KEEP)
        for y in YEARS:
            try:
                raw = urllib.request.urlopen(URL.format(y=y, sid=sid), timeout=120).read()
            except urllib.error.HTTPError as e:
                if e.code in (403, 404):
                    continue
                raise
            for r in csv.DictReader(io.StringIO(raw.decode("utf-8", "replace"))):
                w.writerow([r.get(k, "") for k in KEEP])
                rows += 1
    os.replace(tmp, path)
    return sid, rows


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    with cf.ThreadPoolExecutor(8) as ex:
        for sid, n in ex.map(fetch, STATIONS):
            print(sid, STATIONS[sid][0], n, flush=True)
