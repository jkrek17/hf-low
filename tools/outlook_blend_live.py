#!/usr/bin/env python3
"""Daily job, second layer of the live outlook page: the GFE-style blend of GFS, ECMWF and AIFS 10 m wind for the same 00 UTC run as the gust shading.

    python3 tools/outlook_blend_live.py [--docs docs] [--flat flat] [--force] [--cache DIR]

Reads the run in docs/outlook-live-data.js (written by tools/outlook_live.py), pulls that run's GFS, IFS, AIFS, GEFS (31), ENS (50) and AIFS ENS (51)
10 m wind at +12 to +96 h every 12 h (about 1.8 GB streamed, about 12 minutes), blends them (tools/outlook_blend.py) and writes docs/outlook-live-blend.js
and flat/outlook-live-blend.js, about 190 KB each. It never mixes cycles: the file carries its own run, and the page uses it only when it equals the gust run.
Does nothing (exit 0) when there is no gust run, the run is outside October-April, the blend file already holds this run, or a source has not posted its +96 h
yet (the ECMWF ensemble is usually out about 8 h after the cycle, so the 06:10 UTC job may wait for the 08:40 one). Exits 1 when a pull fails.

Wire format of window.HF_BLEND (append-only): run 'YYYYMMDDHH', made (UTC), res (degrees per cell, 1), lat0 (north edge, 67), fields (names, in blob order), scale (stored = round(value x scale), uint8),
basins.atl / basins.pac = {lon0 (west edge), ny, nx, days: four base64 zlib blobs (or null) holding the fields back to back, field-major, row-major, north row first}. Cell value = highest 0.25 degree
value in the 1 degree cell and, for each day, the higher of its two 12-hourly forecast times. Sustained 10 m wind in knots; p34/p48/p64 are percent chances that the wind occurs, not OPC-agreed rates."""
import argparse, json, os, re, sys, time
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); sys.path.insert(0, HERE)
import outlook_blend as O

DATA_NAME, BLEND_NAME = "outlook-live-data.js", "outlook-live-blend.js"
SEASON_MONTHS = (10, 11, 12, 1, 2, 3, 4)


def run_of(path):
    try:
        return re.search(r'"run":"(\d{10})"', open(path, encoding="utf-8").read()).group(1)
    except (OSError, AttributeError):
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--docs", default=os.path.join(ROOT, "docs")); ap.add_argument("--flat", default=os.path.join(ROOT, "flat"))
    ap.add_argument("--force", action="store_true", help="write even if this run is already in the blend file")
    ap.add_argument("--cache", help="folder for the 0.25 degree step files (testing); reused when present")
    a = ap.parse_args()
    run = run_of(os.path.join(a.docs, DATA_NAME))
    if run is None:
        print("::notice::no gust run in " + DATA_NAME + "; blend not built"); return 0
    if int(run[4:6]) not in SEASON_MONTHS:
        print(f"::notice::run {run} is outside October-April; blend not built"); return 0
    if not a.force and run_of(os.path.join(a.docs, BLEND_NAME)) == run:
        print(f"run {run} is already in {BLEND_NAME}; nothing new"); return 0
    ymd, hh = run[:8], run[8:]
    if not O.ready(ymd, hh):
        print(f"::notice::run {run}: a source has not posted +{max(O.STEPS)} h yet; the blend waits for the next job"); return 0
    t0 = time.time()
    try:
        res = O.run_cycle(ymd, hh, cache=a.cache)
        P = O.pack_live(res, run, datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%MZ"))
    except Exception as e:
        print(f"::error::blend for run {run} failed: {e!r}; nothing written"); return 1
    if len(P["basins"]) != 2 or not all(all(d for d in b["days"]) for b in P["basins"].values()):
        print(f"::error::blend for run {run} is incomplete; nothing written"); return 1
    text = "window.HF_BLEND=" + json.dumps(P, separators=(",", ":")) + ";\n"
    for folder in (a.docs, a.flat):
        os.makedirs(folder, exist_ok=True)
        with open(os.path.join(folder, BLEND_NAME), "w", encoding="utf-8") as f: f.write(text)
    print(f"wrote {BLEND_NAME} ({len(text)} bytes) for run {run} in {time.time() - t0:.0f}s"); return 0


if __name__ == "__main__":
    sys.exit(main())
