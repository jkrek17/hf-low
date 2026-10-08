"""Weekly table that lets the held-out numbers be recomputed without the ERA5 fields (about 0.5 MB).
One row per season-week 1979..2025 (Oct 1 + 30 weeks); PC scores are the primary split's (EOFs from 2004-05..2014-15,
unit variance in that period), window days -7..-1. Indices are raw (not standardised). usage: export_tables.py WORKDIR [OUT]"""
import os, sys, pickle
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hemlib as H
import run as R

work = sys.argv[1]
out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(H.HERE, "results")
ctx = R.Ctx(work)
fz = pickle.load(open(os.path.join(work, "frozen_primary.pkl"), "rb"))
seasons = list(range(1979, 2026))
pcs = fz["eof"].scores(ctx.maps_for("lag1", seasons))
cols = [f"{v[0]}{i+1}" for v in H.VARS for i in range(H.K)]
df = pd.DataFrame(pcs, columns=cols)
df.insert(0, "week", [k for s in seasons for k in range(H.NWEEK)])
df.insert(0, "season", [s for s in seasons for k in range(H.NWEEK)])
ix = ctx.index_table(seasons, "lag1") if min(seasons) >= 2001 else None
ixl = []
import datetime as dt
for s in seasons:
    for k in range(H.NWEEK):
        S = H.week_start(s, k)
        try:
            ixl.append(ctx.idx.window(S - dt.timedelta(days=7), S - dt.timedelta(days=1)))
        except Exception:
            ixl.append({c: np.nan for c in ("nao", "pna", "ao", "oni", "mjo1", "mjo2")})
df = pd.concat([df, pd.DataFrame(ixl)], axis=1)
arch = H.load_archive_counts(None)
for b in H.BASINS:
    t = H.weekly_table(arch, seasons, b)
    df[f"y_{b}"] = t.y.values; df[f"prev_{b}"] = t.prev.values
df["month"] = t.month.values
for kind in ("hf", "depth"):
    d = H.load_proxy_days(kind)
    for b in H.BASINS:
        t = H.weekly_table(d, seasons, b)
        df[f"pA{kind}_{b}"] = t.y.values; df[f"pA{kind}prev_{b}"] = t.prev.values
df.round(5).to_csv(os.path.join(out, "weekly_table.csv.gz"), index=False)
print(df.shape, os.path.getsize(os.path.join(out, "weekly_table.csv.gz")) / 1e6, "MB")
