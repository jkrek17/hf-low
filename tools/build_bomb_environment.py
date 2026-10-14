#!/usr/bin/env python3
"""Write docs/bomb-environment.html: the upper-level environment check and 24 h gale / storm / HF outlook.

    python3 tools/build_bomb_environment.py --fixes fixes_2004.csv.gz --env env_2004.csv.gz \
        --oof oof_outlook_full.csv.gz --cells diagram_cells.csv --family family.csv \
        --reliability reliability.csv --regions regions.csv --matches matches.csv --residual residual_unlisted_both_rules_2004.csv \
        --source "research/era5/bomb_env at <commit>"

Inputs are research/era5/intensity/results (fixes, env), research/era5/bomb_env/results (oof probabilities,
diagram cells, family, reliability), research/era5/gust_depth/results/matches.csv (OPC storm -> pipeline A track) and
the residual table of hf-probability-tracks (proxy storms OPC never listed). ERA5 PROXY, pipeline A, 2004-05 to
2025-26; the outlook probabilities are leave-one-season-out model output, not observations.
The page is self-contained (the flat/ copy has no data/ folder for pages other than index.html), so the data are
inlined from tools/bomb_environment.template.html at the marker below.

Wire format of window.HF_OUTLOOK (append-only):
  tracks[id] = [basin(0 atl,1 pac), season, fixes]; fixes = flat list, 12 numbers per fix:
    time YYYYMMDDHH, lat*10, lon*10, msl hPa*10, gust index kt*10, jet250 kt*10, vadv500 (1e-9 s-2)*100,
    eady (day-1)*1000, P(gale)*1000, P(storm)*1000, P(HF)*1000, 0 (spare)
    P(below gale) = 1 - the three. P is the held-out-season full-model probability that the peak gust index over
    the next 24 h (next two 00/12 UTC fixes) is in that class.
  diagrams[name] = {x, y, xe, ye, cells}: axis names, 7 quantile edges each, and 36 cells x-major
    [n, rapid-deepening rate*10000 (-1 when none scored), P(below gale), P(gale), P(storm), P(HF) of the NEXT-24 h class, all *10000].
    The third diagram, trough_baroclinic, has x = vadv500 and y = eady (added after the first two; same cell layout).
  regions = {cuts: [vadv500 cut, eady cut], rows: [[name, n, n_scored, rd rate, rd lo90, rd hi90, share of rapid-deepening fixes,
    share of fixes, P(below gale), P(gale), P(storm), P(HF)] x 4 in the order quiet, baroclinic only, forcing only, coupled]}.
    Region index of a fix = 2 * (vadv500 >= cuts[0]) + (eady >= cuts[1]); the cuts are the upper thirds over all fixes.
"""
import argparse
import json
import os

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MARK = "/*@@DATA@@*/"
RULE_MARK = "<!--@@GFSRULE@@-->"  # replaced by tools/gfs_gust_rule.html (made by tools/gfs_gust_rule.py)


def main():
    ap = argparse.ArgumentParser()
    for k in ("fixes", "env", "oof", "cells", "family", "reliability", "regions", "matches", "residual", "source"):
        ap.add_argument("--" + k, required=True)
    ap.add_argument("--skill", required=True, help="research/era5/bomb_env/results/skill.txt")
    ap.add_argument("--template", default=os.path.join(ROOT, "tools", "bomb_environment.template.html"))
    ap.add_argument("--out", default=os.path.join(ROOT, "docs", "bomb-environment.html"))
    a = ap.parse_args()

    m = pd.read_csv(a.matches)
    m = m[(m.peak == "HF") & m.track.notna()]
    tr = set(m.track.astype(int))
    res = pd.read_csv(a.residual)
    tr |= set(res.track.astype(int))
    F = pd.read_csv(a.fixes, dtype={"time": str})
    F = F[F.track.isin(tr)]
    E = pd.read_csv(a.env, dtype={"time": str})
    O = pd.read_csv(a.oof, dtype={"time": str})
    D = F.merge(E, on=["track", "time"]).merge(O, on=["track", "time"])
    D = D[D.season <= 2025].sort_values(["track", "time"])
    tracks = {}
    for tid, g in D.groupby("track"):
        flat = []
        for r in g.itertuples():
            flat += [int(r.time), round(r.lat * 10), round(r.lon * 10), round(r.msl * 10), round(r.g800 * 10),
                     round(r.jet250 * 10) if r.jet250 == r.jet250 else -1,
                     round(r.vadv500 * 100) if r.vadv500 == r.vadv500 else 0,
                     round(r.eady * 1000) if r.eady == r.eady else -1,
                     round(r.p1 * 1000), round(r.p2 * 1000), round(r.p3 * 1000), 0]
        tracks[str(int(tid))] = [0 if g.basin.iloc[0] == "atl" else 1, int(g.season.iloc[0]), flat]

    C = pd.read_csv(a.cells)
    dg = {}
    for name, g in C.groupby("diagram"):
        g = g.sort_values(["bx", "by"])
        xe = [round(float(v), 3) for v in list(g.drop_duplicates("bx").x_lo) + [g.x_hi.max()]]
        ye = [round(float(v), 3) for v in list(g.drop_duplicates("by").y_lo) + [g.y_hi.max()]]
        cells = [[int(r.n), int(round(r.rd_rate * 10000)) if r.rd_rate == r.rd_rate else -1,
                  int(round(r.p_below_gale * 10000)), int(round(r.p_gale * 10000)),
                  int(round(r.p_storm * 10000)), int(round(r.p_hf * 10000))] for r in g.itertuples()]
        dg[name] = {"x": g.xvar.iloc[0], "y": g.yvar.iloc[0], "xe": xe, "ye": ye, "cells": cells}
    R = pd.read_csv(a.regions)
    c = dg["trough_baroclinic"]
    regions = {"cuts": [c["xe"][4], c["ye"][4]],
               "rows": [[r.region, int(r.n), int(r.n_scored), round(r.rd_rate, 5), round(r.rd_lo90, 5), round(r.rd_hi90, 5),
                         round(r.share_of_rd, 4), round(r.share_of_fixes, 4), round(r.p_below_gale, 4), round(r.p_gale, 4),
                         round(r.p_storm, 4), round(r.p_hf, 4)] for r in R.itertuples()]}
    fam = pd.read_csv(a.family).round(4).to_dict("records")
    rel = pd.read_csv(a.reliability).round(4)
    skill = {}
    for line in open(a.skill):
        s = line.strip()
        for k in ("clim", "state", "A", "B", "full", "full_nogust"):
            if s.startswith(k + " "):
                skill.setdefault(k, []).append(s.split()[1])
    data = {"meta": {"what": "ERA5 PROXY, pipeline A (research/era5/hf_history). Held-out-season model output, not observations.",
                     "source": a.source, "tracks": len(tracks), "fixes": int(len(D)), "base_rd": 0.0343,
                     "cutoffs_kt": [38, 54, 71.7]},
            "family": fam, "skill": skill, "diagrams": dg, "regions": regions, "tracks": tracks}
    blob = "window.HF_OUTLOOK=" + json.dumps(data, separators=(",", ":")) + ";"
    html = open(a.template).read()
    assert MARK in html
    html = html.replace(RULE_MARK, open(os.path.join(ROOT, "tools", "gfs_gust_rule.html")).read().rstrip("\n"))
    open(a.out, "w").write(html.replace(MARK, blob))
    print(a.out, len(tracks), "tracks", len(D), "fixes", len(blob) // 1024, "KB data")


if __name__ == "__main__":
    main()
