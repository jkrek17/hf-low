#!/usr/bin/env python3
"""Write docs/outlook-live.html: the experimental page that shades the newest GFS run by how often OPC agreed a hurricane-force or storm-force low was near.

    python3 tools/build_outlook_live.py [--payload data/outlook_replay/replay_payload.json]

The page is static. What changes every day is only docs/outlook-live-data.js (written by tools/outlook_live.py, about 8 KB): the run time and four
run-length band strings. The page inlines the calibration (tables_all of the Outlook replay payload: five seasons, October-April 2021-22 to 2025-26),
the 1-degree ocean grid and a trimmed coastline, so it must be rebuilt when that payload changes.

Wire format of window.HF_LIVE_STATIC (append-only; same names as window.HF_REPLAY in build_outlook_replay.py):
  grid {lat, lon, basin}, radius_nm, truth_nm, bands_kt, tables_all[basin+day] = {hf, st, n, c}, coast.
Wire format of window.HF_LIVE (outlook-live-data.js, append-only):
  run 'YYYYMMDDHH' (GFS cycle), made (UTC time the file was written), days = four run-length strings of band indices (digit, count, comma) over grid order;
  an empty string means that forecast day could not be built.
"""
import argparse
import json
import os

import build_outlook_replay as R

ROOT = R.ROOT


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--payload", default=os.path.join(ROOT, "data", "outlook_replay", "replay_payload.json"))
    ap.add_argument("--coast", default=os.path.join(ROOT, "docs", "data", "coastlines.js"))
    ap.add_argument("--template", default=os.path.join(ROOT, "tools", "outlook_live.template.html"))
    ap.add_argument("--out", default=os.path.join(ROOT, "docs", "outlook-live.html"))
    a = ap.parse_args()
    p = json.load(open(a.payload))
    d = {k: p[k] for k in ("grid", "radius_nm", "truth_nm", "bands_kt", "tables_all")}
    d["coast"] = R.coast(a.coast)
    t = open(a.template, encoding="utf-8").read()
    assert t.count(R.MARK) == 1
    html = t.replace(R.MARK, json.dumps(d, separators=(",", ":")))
    open(a.out, "w", encoding="utf-8").write(html)
    print(a.out, round(len(html) / 1e6, 2), "MB")


if __name__ == "__main__":
    main()
