"""Download NDBC historical standard-met files and anemometer heights for the Gulf Stream study.
Output: /mnt/project-files/gulf-stream/ndbc/<station>_<year>.txt.gz and anemometer_heights.csv"""
import os, re, sys, urllib.request, time
OUT = sys.argv[1] if len(sys.argv) > 1 else "/mnt/project-files/gulf-stream/ndbc"
STATIONS = ["41002", "41001", "41025", "41013", "44014", "41004", "44009"]
YEARS = range(2004, 2026)
os.makedirs(OUT, exist_ok=True)
def get(url):
    for k in range(3):
        try:
            with urllib.request.urlopen(url, timeout=60) as r: return r.read()
        except Exception as e:
            if k == 2: return None
            time.sleep(2)
rows = []
for s in STATIONS:
    for y in YEARS:
        f = f"{OUT}/{s}_{y}.txt.gz"
        if os.path.exists(f): continue
        b = get(f"https://www.ndbc.noaa.gov/data/historical/stdmet/{s}h{y}.txt.gz")
        if b: open(f, "wb").write(b)
        else: print("missing", s, y)
    p = get(f"https://www.ndbc.noaa.gov/station_page.php?station={s}")
    t = p.decode("latin1") if p else ""
    m = re.search(r"Anemometer height:</b>\s*([\d.]+)\s*m", t)
    rows.append((s, m.group(1) if m else ""))
open(f"{OUT}/anemometer_heights.csv", "w").write("station,anemometer_m\n" + "\n".join(f"{a},{b}" for a, b in rows) + "\n")
print(rows)
