"""POST HOC, secondary: does the hemispheric pattern index explain the Knox excess?

Uses the out-of-sample weekly pattern index committed by the hemispheric thread
(research/era5/hemispheric/results/oos_index_{atl,pac}.csv, hf-low PR 41; one value per season and
Oct-1-based 7-day week, each season's index from a fit on the other 21 seasons, built from ERA5 proxy
fields over days -7..-1 before the week). Not part of PREREGISTRATION.md; archive tier only, because the
index was built against the archive's weekly counts.

Times are permuted only among events in the same season-month AND the same tercile of the index of the
event's week. Same recipe as posthoc.py group "lag_tercile".
usage: posthoc_pattern.py OUT_DIR [--nperm N]
"""
import argparse
import json
import os

import numpy as np
import pandas as pd

import cluster as C
import posthoc as P

HEM = os.path.join(C.REPO, "research/era5/hemispheric/results")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--nperm", type=int, default=4000)
    a = ap.parse_args()
    rows, cfg = C.tier_data("arch")
    s0 = cfg["s0"]
    rng = np.random.default_rng(C.SEED + 77)
    res = {}
    for basin in ("atl", "pac"):
        ev, _ = C.window(rows, basin, cfg["s0"], cfg["s1"])
        tab = pd.read_csv(os.path.join(HEM, f"oos_index_{basin}.csv"))
        idx = tab.pivot(index="season", columns="week", values="idx")
        week = np.array([min(int(e["t"] // 168), 29) for e in ev])
        v = np.array([idx.loc[e["season"], w] for e, w in zip(ev, week)])
        q = np.quantile(v, [1 / 3, 2 / 3])
        terc = np.digitize(v, q)
        season = np.array([e["season"] for e in ev])
        t = np.array([e["t"] for e in ev])
        mon = C.MONIDX[np.clip((t // 24).astype(int), 0, C.NDAY - 1)]
        base = (season - season.min()) * 7 + mon
        out = {}
        for name, grp in (("base", base), ("pattern_tercile", base * 3 + terc)):
            o, m, tot = P.knox_bins(ev, grp, rng, a.nperm)
            out[name] = P.summarise(o, m, tot, rng)
        b, p = out["base"]["ratio_72h"], out["pattern_tercile"]["ratio_72h"]
        out["share_of_excess_explained"] = float(1 - (p - 1) / (b - 1))
        res[basin] = out
        print(basin, b, p, out["pattern_tercile"]["ratio_72h_ci"], out["share_of_excess_explained"], flush=True)
    json.dump(dict(nperm=a.nperm, basins=res), open(os.path.join(a.out, "arch_posthoc_pattern.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
