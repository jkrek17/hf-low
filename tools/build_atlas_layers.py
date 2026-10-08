#!/usr/bin/env python3
"""Build docs/data/atlas.js: the climatology atlas layers for the map.

Derived data only. The method is research/era5/climatology_atlas (atlas.py,
historic.py), which stays on the research branch; this tool redoes just the
three table builds the map needs from the same committed inputs and checks the
result against the numbers the atlas reports.

    python3 tools/build_atlas_layers.py \
        --catalog <hf_history>/results/era5_hf_catalog.csv \
        --tracks  <hf_history>/results/era5_hf_catalog_tracks.csv \
        --motion  <climatology_atlas>/results/motion_steps.csv.gz \
        --historic <climatology_atlas>/results/historic_storms.csv \
        --source  <research-branch commit>

Two sources, always named: ARCHIVE (docs/data/hf-lows.json, HF-category fixes,
seasons 2004-05 to 2025-26) and PROXY (ERA5 pipeline A, 800 km gust index at
or above 71.7 kt, same seasons). Season = 1 June to 31 May. Tropical-cyclone-
linked events are included, as in the atlas.

Wire format (append-only):
  HF_ATLAS.box[src][basin] = {lat0, lon0, dlat, dlon, nlat, nlon,
        cells: [[i, j, c1..c12], ...]}      HF-centre FIXES per calendar month in
        the 5 x 10 degree box (row i from lat0, column j from lon0); one fix =
        6 hours, so hours per season = fixes * 6 / HF_ATLAS.meta.seasons
  HF_ATLAS.motion[src][basin] = [[box lat0, box lon0, mean u kt*10 (east),
        mean v kt*10 (north), n steps], ...]   boxes with >= 25 steps only
  HF_ATLAS.historic = [[basin, ranked_by, track, season, date YYYYMMDD,
        min pressure hPa*10, gust800 kt*10, lat*4, lon*4, tc(0/1)], ...]
        within-era ERA5 proxy candidates 1979-80 to 2003-04: NOT confirmed
"""
import argparse
import csv
import gzip
import json
import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
THR = 71.7
FIRST, LAST = 2004, 2025
NSEAS = LAST - FIRST + 1
BASINS = ("atl", "pac")
LAT_EDGES = np.arange(20, 76, 5)
LON_EDGES = {"atl": np.arange(-100, 31, 10), "pac": np.arange(120, 251, 10)}
MIN_STEPS = 25


def lonfix(lon, basin):
    lon = np.asarray(lon, float)
    return ((lon + 180) % 360) - 180 if basin == "atl" else lon % 360


def archive_fixes():
    d = json.load(open(os.path.join(ROOT, "docs", "data", "hf-lows.json")))
    rows = []
    for raw in d["lows"]:
        l = dict(zip(d["lowFields"], raw))
        if l["season"] < FIRST:
            continue
        f = pd.DataFrame(l["fixes"], columns=d["fixFields"] if len(l["fixes"][0]) == len(d["fixFields"]) else None)
        f.columns = ["date", "lat", "lon", "cat", "pres"][:len(f.columns)]
        f["t"] = pd.to_datetime(f.date.astype("int64").astype(str), format="%Y%m%d%H")
        f = f.drop_duplicates("t")
        h = f[f.cat == "HF"]
        for _, r in h.iterrows():
            rows.append((l["basin"], r.t.month, r.lat, lonfix(r.lon, l["basin"])))
    return pd.DataFrame(rows, columns=["basin", "month", "lat", "lon"])


def proxy_fixes(catalog, tracks):
    C = pd.read_csv(catalog)
    E = C[(C.role == "event") & (C.season >= FIRST) & (C.season <= LAST)].set_index("track")
    T = pd.read_csv(tracks)
    T = T[T.track.isin(E.index)].copy()
    T["t"] = pd.to_datetime(T.time.astype(str), format="%Y%m%d%H")
    T = T[T.basin.fillna("") != ""]
    T = T[T.g800 >= THR]
    T["b"] = T.track.map(E.basin)
    out = []
    for b in BASINS:
        x = T[T.b == b]
        out.append(pd.DataFrame({"basin": b, "month": x.t.dt.month.values, "lat": x.lat.values,
                                 "lon": lonfix(x.lon.values, b)}))
    return pd.concat(out, ignore_index=True)


def box_cells(fx, basin):
    f = fx[fx.basin == basin]
    le = LON_EDGES[basin]
    cells = {}
    i = np.floor((f.lat.values - LAT_EDGES[0]) / 5).astype(int)
    j = np.floor((f.lon.values - le[0]) / 10).astype(int)
    ok = (i >= 0) & (i < len(LAT_EDGES) - 1) & (j >= 0) & (j < len(le) - 1)
    for a, b, m in zip(i[ok], j[ok], f.month.values[ok]):
        c = cells.setdefault((int(a), int(b)), [0] * 12)
        c[int(m) - 1] += 1
    total = sum(sum(c) for c in cells.values())
    return [[a, b] + c for (a, b), c in sorted(cells.items())], total


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--tracks", required=True)
    ap.add_argument("--motion", required=True)
    ap.add_argument("--historic", required=True)
    ap.add_argument("--source", required=True)
    ap.add_argument("--out", default=os.path.join(ROOT, "docs", "data", "atlas.js"))
    a = ap.parse_args()

    box, busiest = {}, {}
    for src, fx in (("archive", archive_fixes()), ("proxy", proxy_fixes(a.catalog, a.tracks))):
        box[src] = {}
        for b in BASINS:
            cells, total = box_cells(fx, b)
            box[src][b] = {"lat0": 20, "lon0": int(LON_EDGES[b][0]), "dlat": 5, "dlon": 10, "nlat": len(LAT_EDGES) - 1,
                           "nlon": len(LON_EDGES[b]) - 1, "cells": cells}
            top = max(cells, key=lambda c: sum(c[2:]))
            busiest[(src, b)] = (sum(top[2:]) * 6.0 / NSEAS, 20 + 5 * top[0], int(LON_EDGES[b][0]) + 10 * top[1], total * 6.0 / NSEAS)

    mo = pd.read_csv(a.motion)
    motion = {}
    for src in ("archive", "proxy"):
        motion[src] = {}
        for b in BASINS:
            m = mo[(mo.source == src) & (mo.basin == b)].copy()
            m["bi"] = (m.latc // 5) * 5
            m["bj"] = (m.lonc // 10) * 10
            g = m.groupby(["bi", "bj"]).agg(u=("u", "mean"), v=("v", "mean"), n=("u", "size")).reset_index()
            g = g[g.n >= MIN_STEPS]
            motion[src][b] = [[int(r.bi), int(r.bj), round(r.u * 10), round(r.v * 10), int(r.n)] for r in g.itertuples()]

    hist = []
    with open(a.historic, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            hist.append([r["basin"], r["ranked_by"], int(r["track"]), int(r["season"]), int(r["date"]),
                         round(float(r["minp"]) * 10), round(float(r["gust800_kt"]) * 10),
                         round(float(r["lat"]) * 4), round(float(r["lon"]) * 4), 1 if r["tc_linked"] == "True" else 0])
    assert all(h[3] < FIRST for h in hist), "historic candidates must be before 2004-05"

    meta = {"what": "ERA5 proxy + OPC archive climatology atlas layers (research/era5/climatology_atlas); boxes 5 x 10 degrees",
            "source": a.source, "seasons": NSEAS, "first_season": FIRST, "last_season": LAST, "min_steps": MIN_STEPS,
            "motion_steps": int(len(mo)), "historic": len(hist),
            "fields": ["box[src][basin].cells: [i, j, fixes per month 1..12]",
                       "motion[src][basin]: [box lat0, box lon0, u kt*10, v kt*10, n]",
                       "historic: [basin, ranked_by, track, season, date, minp*10, gust*10, lat*4, lon*4, tc]"]}
    with open(a.out, "w", encoding="utf-8") as fh:
        fh.write("/* generated by tools/build_atlas_layers.py - do not edit. ERA5 PROXY and OPC archive, descriptive only. */\n")
        fh.write("window.HF_ATLAS = " + json.dumps({"meta": meta, "box": box, "motion": motion, "historic": hist},
                                                    separators=(",", ":")) + ";\n")
    print("%.3f MB -> %s" % (os.path.getsize(a.out) / 1e6, a.out))
    for k, v in busiest.items():
        print("busiest box %-8s %s: %.1f h/season at %dN, lon %d  (all boxes %.0f h/season)" % (k[0], k[1], v[0], v[1], v[2], v[3]))
    print("motion boxes:", {s: {b: len(motion[s][b]) for b in BASINS} for s in motion}, " historic rows:", len(hist))


if __name__ == "__main__":
    main()
