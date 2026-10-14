#!/usr/bin/env python3
"""Daily job for the experimental live outlook page: pull the newest GFS 00 UTC run, shade it like the Outlook replay page, write the small data file.

    python3 tools/outlook_live.py [--docs docs] [--flat flat] [--cycle YYYYMMDD] [--force] [--check-against replay_payload.json]

What it does (same method as research/era5/outlook_replay/, calibrated on five past seasons, October-April):
  * GFS 0.25 degree GUST:surface at f006..f096 every 6 h (NOAA open-data bucket, .idx byte ranges, our own GRIB2 reader tools/grib2_t53.py).
  * Day d (1..4) maximum = the four valid times 24(d-1)+6 .. 24d h, on a 0.5 degree ocean grid (land zeroed, 2x2 maximum), kept in 0.5 kt steps.
  * Maximum within R_d = 60/120/270/270 nm of each 1-degree ocean point, then the gust band (<40, 40-48, 48-56, 56-64, 64-72, 72-80, >=80 kt).
  * Writes docs/outlook-live-data.js and flat/outlook-live-data.js (identical): run time, when it was made, four run-length band strings, and the
    context overlays of tools/outlook_overlays.py (lows, trough axes, jet cores, Coupled flags) for forecast hours 24/48/72/96 (null where a lead failed).
    The percentages come from the page (tables from data/outlook_replay/replay_payload.json), not from this file.
Does nothing (exit 0) outside October-April, or when the newest complete run is already in the data file; exits 1 if a pull fails.
GFS gust is a model forecast and a proxy for the sustained wind OPC warns on; the page's percentages are shares of past point-days where OPC had an
HF (or storm-force) low fix within 324 nm, an analyst's judgement, not verification. Standard library plus numpy."""
import argparse, json, os, re, struct, sys, time, urllib.error, urllib.request
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta, timezone
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); sys.path.insert(0, HERE)
import grib2_t53 as G

BUCKET = "https://noaa-gfs-bdp-pds.s3.amazonaws.com"; UA = "hf-low-outlook-live/0.1 (github.com/jkrek17/hf-low)"
LEADS = list(range(6, 97, 6)); KT = 1.9438445
RD = [60, 120, 270, 270]; BANDS = [0, 40, 48, 56, 64, 72, 80, 999]
SEASON_MONTHS = (10, 11, 12, 1, 2, 3, 4)
DATA_NAME = "outlook-live-data.js"

# --- geometry, identical to research/era5/outlook_replay/tabulate.py (checked against replay_payload.json by --check-against) ---
_z = np.load(os.path.join(ROOT, "data", "outlook_live", "ocean_0p25_75N_20N.npz"))
LSM = np.unpackbits(_z["ocean"], axis=1)[:, :1440].astype(bool)                     # 220 rows 75N..20.25N, 1440 cols, ocean True
LAT = 74.75 - 0.5 * np.arange(110); LON = 0.125 + 0.5 * np.arange(720)
OC = LSM.reshape(110, 2, 720, 2).sum(axis=(1, 3)) >= 3
ER, EC = np.meshgrid(np.arange(0, 110, 2), np.arange(0, 720, 2), indexing="ij")
BAS = np.where((EC * 0.5 + 0.125 >= 260) | (EC * 0.5 + 0.125 <= 10), 1, np.where((EC * 0.5 + 0.125 >= 120) & (EC * 0.5 + 0.125 <= 240), 2, 0))
OK = OC[ER, EC] & (BAS > 0)
FLAT = np.nonzero(OK.ravel())[0]

_OFF = {}
def offsets(R, row):
    key = (R, row // 10)
    if key not in _OFF:
        lat = 74.75 - 0.5 * (row // 10 * 10 + 5); dx = 30.0 * np.cos(np.radians(lat))
        ky = int(R // 30); kx = int(R // dx)
        yy, xx = np.meshgrid(np.arange(-ky, ky + 1), np.arange(-kx, kx + 1), indexing="ij")
        m = ((yy * 30.0) ** 2 + (xx * dx) ** 2) <= R ** 2
        _OFF[key] = (yy[m], xx[m], ky, kx)
    return _OFF[key]

def maxwithin(g, R):
    out = np.zeros(ER.shape, np.float32); PADY, PADX = 12, 120
    P = np.pad(g, ((PADY, PADY), (0, 0))); P = np.pad(P, ((0, 0), (PADX, PADX)), mode="wrap")
    for b in range(0, 110, 10):
        sel = (ER >= b) & (ER < b + 10)
        if not sel.any(): continue
        dy, dx, ky, kx = offsets(R, b)
        r = ER[sel] + PADY; c = EC[sel] + PADX
        out[sel] = P[r[:, None] + dy[None, :], c[:, None] + dx[None, :]].max(axis=1)
    return out

def rle(a):
    s = []; prev = a[0]; n = 1
    for v in a[1:]:
        if v == prev: n += 1
        else: s.append(f"{prev}{n},"); prev = v; n = 1
    s.append(f"{prev}{n},"); return "".join(s)

# --- fetching ---
class Missing(Exception):
    pass

def get(u, rng=None, tries=6):
    h = {"User-Agent": UA}; want = None
    if rng:
        h["Range"] = rng; a, b = rng.split("=")[1].split("-"); want = int(b) - int(a) + 1
    err = None; n404 = 0
    for k in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(u, headers=h), timeout=120) as r:
                b = r.read(); cl = r.headers.get("Content-Length")
            if (cl is not None and len(b) != int(cl)) or (want is not None and len(b) != want):
                raise IOError("short read")
            return b
        except urllib.error.HTTPError as e:
            if e.code in (403, 404):
                n404 += 1
                if n404 >= 2: raise Missing(u)
            err = e
        except Exception as e:
            err = e
        time.sleep(min(2 ** k, 20))
    raise RuntimeError(f"fetch failed {u}: {err}")

def walk(u, i):
    off = 0
    for k in range(i + 1):
        h = get(u, f"bytes={off}-{off + 15}")
        if h[:4] != b"GRIB": raise RuntimeError("walk lost sync")
        L = struct.unpack(">Q", h[8:16])[0]
        if k == i: return get(u, f"bytes={off}-{off + L - 1}")
        off += L

def gust_field(day, lead):
    u = f"{BUCKET}/gfs.{day}/00/atmos/gfs.t00z.pgrb2.0p25.f{lead:03d}"
    idx = get(u + ".idx"); rows = []
    for line in idx.decode().splitlines():
        p = line.split(":")
        if len(p) >= 6 and p[1].isdigit(): rows.append((int(p[1]), p[3], p[4]))
    for i, (off, v, l) in enumerate(rows):
        if v == "GUST" and l == "surface":
            end = rows[i + 1][0] - 1
            blob = get(u, f"bytes={off}-{end}")
            try:
                return G.decode_message(blob)[60:280] * KT
            except (AssertionError, struct.error, IndexError, KeyError):
                return G.decode_message(walk(u, i))[60:280] * KT        # stale .idx: walk GRIB message lengths
    raise Missing(u + " no GUST")

def available(day):
    """A run is usable once its last lead (f096) has been written."""
    try:
        get(f"{BUCKET}/gfs.{day}/00/atmos/gfs.t00z.pgrb2.0p25.f096.idx", tries=3); return True
    except Missing:
        return False                    # not written yet; any other error propagates and fails the job

def day_maxima(day, workers=4):
    """(4,110,720) uint8, 0.5 kt steps, 255 where a sample of that day is missing; plus the list of missing leads."""
    def one(lead):
        try: return lead, gust_field(day, lead)
        except Missing: return lead, None
    with ThreadPoolExecutor(workers) as ex: res = dict(ex.map(one, LEADS))
    samples = {}
    for lead, g in res.items():
        if g is None: continue
        g = np.where(LSM, g, 0.0); samples[lead] = g.reshape(110, 2, 720, 2).max(axis=(1, 3))
    d = np.full((4, 110, 720), 255, np.uint8)
    for k in range(4):
        ls = [24 * k + 6 * j for j in (1, 2, 3, 4)]
        if all(l in samples for l in ls):
            d[k] = np.clip(np.round(np.max([samples[l] for l in ls], axis=0) * 2), 0, 254).astype(np.uint8)
    return d, [l for l in LEADS if l not in samples]

def bands_from(d):
    out = []
    for k in range(4):
        if (d[k] == 255).any(): out.append(""); continue
        pred = maxwithin(d[k].astype(np.float32) / 2, RD[k])
        out.append(rle((np.digitize(pred, BANDS) - 1).ravel()[FLAT].astype(int).tolist()))
    return out

def current_run(path):
    try:
        t = open(path, encoding="utf-8").read(); return re.search(r'"run":"(\d{10})"', t).group(1)
    except (OSError, AttributeError):
        return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--docs", default=os.path.join(ROOT, "docs")); ap.add_argument("--flat", default=os.path.join(ROOT, "flat"))
    ap.add_argument("--cycle", help="YYYYMMDD: build this 00 UTC run instead of the newest (testing)")
    ap.add_argument("--force", action="store_true", help="write even if this run is already in the data file")
    ap.add_argument("--check-against", help="replay_payload.json: compare this run's band strings with the same cycle there, write nothing")
    a = ap.parse_args()
    now = datetime.now(timezone.utc)
    if a.cycle: cands = [a.cycle]
    else: cands = [(now - timedelta(days=k)).strftime("%Y%m%d") for k in (0, 1)]
    day = next((c for c in cands if available(c)), None)
    if day is None:
        print("::notice::no complete GFS 00 UTC run found for " + ", ".join(cands) + "; nothing written"); return 0
    if int(day[4:6]) not in SEASON_MONTHS and not a.cycle:
        print(f"::notice::{day} is outside October-April (the calibrated season); nothing written"); return 0
    run = day + "00"; path = os.path.join(a.docs, DATA_NAME)
    if not a.force and not a.check_against and current_run(path) == run:
        print(f"run {run} is already in {path}; nothing new"); return 0
    t0 = time.time(); d, miss = day_maxima(day); days = bands_from(d)
    print(f"run {run}: pulled and shaded in {time.time() - t0:.0f}s; missing leads {miss}; days available {[bool(x) for x in days]}")
    if miss or not all(days):
        print(f"::error::run {run} is incomplete (missing leads {miss}); nothing written, the page keeps the previous run"); return 1
    if a.check_against:
        P = json.load(open(a.check_against)); c = next((x for x in P["cases"] if x["cycle"] == day), None)
        if c is None: print("cycle not in payload"); return 2
        ok = [x == y for x, y in zip(days, c["days"])]; print("band strings equal to the replay payload, days 1-4:", ok)
        return 0 if all(ok) else 3
    made = now.strftime("%Y-%m-%dT%H:%MZ")
    import outlook_overlays as V                                          # context overlays; a failure here leaves a day without them, never without shading
    ov = V.for_cycle(day)
    print(f"overlays for days 1-4: {[None if x is None else (len(x['lows']), len(x['jet']), len(x['trough'])) for x in ov]} (lows, jet runs, trough lines)")
    text = "window.HF_LIVE=" + json.dumps({"run": run, "made": made, "days": days, "ov": ov}, separators=(",", ":")) + ";\n"
    for folder in (a.docs, a.flat):
        os.makedirs(folder, exist_ok=True)
        with open(os.path.join(folder, DATA_NAME), "w", encoding="utf-8") as f: f.write(text)
    print(f"wrote {DATA_NAME} ({len(text)} bytes) for run {run}"); return 0

if __name__ == "__main__":
    sys.exit(main())
