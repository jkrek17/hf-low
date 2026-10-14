#!/usr/bin/env python3
"""Write docs/outlook-replay.html: replays of past GFS runs shaded by how often OPC agreed a hurricane-force or storm-force low was near.

    python3 tools/build_outlook_replay.py [--payload data/outlook_replay/replay_payload.json]

Input is research/era5/outlook_replay/results/replay_payload.json (copied to data/outlook_replay/ so this build runs on main). GFS forecast
gust (a proxy for the sustained wind OPC warns on) against OPC archive fixes, five seasons (October-April 2021-22 to 2025-26), leave-one-season-out
tables. The page is self-contained (the flat/ copy has no data/ folder for pages other than index.html), so the payload and a trimmed
coastline are inlined from tools/outlook_replay.template.html at the marker below.

Wire format of window.HF_REPLAY (append-only):
  grid = {lat, lon (east, 0..360), basin (1 Atlantic 100W-10E, 2 Pacific 120E-120W)}: the 1-degree ocean evaluation points, in order.
  radius_nm = shading radius per forecast day 1..4; truth_nm = radius of the OPC truth; bands_kt = the seven gust bands.
  tables[season][basin+day] = {hf, st (percent OPC hurricane-force / storm-force-or-stronger agreed per band), n (point-days), c (runs with an HF event)},
    built without that season; tables_all uses all five.
  cases[] = {cycle YYYYMMDD (00 UTC), season, why, days: 4 run-length strings of band indices (digit, count and a comma) over grid order,
    fixes: [lat, lon east, 'HF'|'S', hours after the run]}.
  coast = polygons of [lon, lat] (Natural Earth 110m via docs/data/coastlines.js), trimmed to 15N-80N.
"""
import argparse
import json
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MARK = "/*@@DATA@@*/"


def coast(path):
    txt = open(path, encoding="utf-8").read()
    obj = json.loads(txt[txt.index("{"): txt.rindex("}") + 1].rstrip(";"))
    out = []
    for poly in obj["polygons"]:
        if not any(15 <= p[1] <= 80 for p in poly):
            continue
        out.append([[round(p[0], 1), round(max(10, min(80, p[1])), 1)] for p in poly])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--payload", default=os.path.join(ROOT, "data", "outlook_replay", "replay_payload.json"))
    ap.add_argument("--coast", default=os.path.join(ROOT, "docs", "data", "coastlines.js"))
    ap.add_argument("--template", default=os.path.join(ROOT, "tools", "outlook_replay.template.html"))
    ap.add_argument("--out", default=os.path.join(ROOT, "docs", "outlook-replay.html"))
    a = ap.parse_args()
    d = json.load(open(a.payload))
    d["coast"] = coast(a.coast)
    t = open(a.template, encoding="utf-8").read()
    assert t.count(MARK) == 1
    html = t.replace(MARK, json.dumps(d, separators=(",", ":")))
    open(a.out, "w", encoding="utf-8").write(html)
    print(a.out, round(len(html) / 1e6, 2), "MB,", len(d["cases"]), "cases,", len(d["coast"]), "coast polygons")


if __name__ == "__main__":
    main()
