"""In-situ check at Prins Christian Sund (ISD 043900, 60.03N 43.12W, 19 m), the
Cape Farewell station Moore and Renfrew used for Greenland tip jets.

Question: when pipeline A's gust maximum (ERA5 proxy) sits in the Cape
Farewell region, does the station see a strong wind at that hour, more than at
other Atlantic HF-strength times?

Groups (fixes from gustloc_fixes.csv, 2004-05 on, one row per time):
  CF   the gust maximum lies in 59-61.5N, 35-50W (any centre position)
  HI   a fix at or north of 60N, maximum elsewhere
  LO   all fixes south of 60N
Station value: the 10 m mean wind (ISD WND) and gust (OC1) in the report
nearest the time, within +-1 h. Quality codes 2, 3, 6, 7 are dropped.

Data: https://noaa-global-hourly-pds.s3.amazonaws.com/<year>/04390099999.csv,
cached in work/isd/ (ignored). Writes pcs-result.txt.
"""
import os, glob
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
KT = 1 / 0.514444
BAD = set("2367")
out = []


def say(*a):
    s = " ".join(str(x) for x in a); print(s); out.append(s)


def parse(path):
    d = pd.read_csv(path, usecols=lambda c: c in ("DATE", "WND", "OC1"), dtype=str)
    w = d.WND.str.split(",", expand=True)
    spd = pd.to_numeric(w[3], errors="coerce")
    ok = (spd != 9999) & ~w[4].isin(BAD)
    d["wind"] = np.where(ok, spd / 10 * KT, np.nan)
    d["wdir"] = np.where(ok & (w[0] != "999"), pd.to_numeric(w[0], errors="coerce"), np.nan)
    if "OC1" in d and d.OC1.notna().any():
        g = d.OC1.fillna("9999,9").str.split(",", expand=True)
        gs = pd.to_numeric(g[0], errors="coerce")
        d["gust"] = np.where((gs != 9999) & ~g[1].isin(BAD), gs / 10 * KT, np.nan)
    else:
        d["gust"] = np.nan
    d["t"] = pd.to_datetime(d.DATE)
    return d[["t", "wind", "wdir", "gust"]]


S = pd.concat([parse(p) for p in sorted(glob.glob(os.path.join(HERE, "work", "isd", "*.csv")))
               if os.path.getsize(p) > 1000]).dropna(subset=["wind"]).sort_values("t")
G = pd.read_csv(os.path.join(HERE, "gustloc_fixes.csv"))
G["t"] = pd.to_datetime(G.time.astype(str), format="%Y%m%d%H")
G["cf"] = (G.max_lat >= 59) & (G.max_lat <= 61.5) & (G.max_lon >= 310) & (G.max_lon <= 325)
T = G.groupby("t").agg(cf=("cf", "any"), hi=("fix_lat", lambda x: (x >= 60).any()),
                       gmax=("g800", "max")).reset_index()
T["grp"] = np.where(T.cf, "CF", np.where(T.hi, "HI", "LO"))
m = pd.merge_asof(T.sort_values("t"), S, on="t", direction="nearest", tolerance=pd.Timedelta("1h"))
say("PRINS CHRISTIAN SUND (ISD 043900) at Atlantic HF-strength times, pipeline A (ERA5 proxy), 2004-05 on")
say(f"station reports 2004-2025: {len(S)}; station climatology of mean wind: median {S.wind.median():.1f} kt, "
    f"95th pct {S.wind.quantile(.95):.1f}, 99th {S.wind.quantile(.99):.1f}")
for gname in ("CF", "HI", "LO"):
    x = m[m.grp == gname]
    v = x.wind.dropna()
    gg = x.gust.dropna()
    say(f"{gname}: times {len(x):4d}, with a report {len(v):4d}; mean wind median {v.median():5.1f} kt, "
        f">=34 kt {100 * (v >= 34).mean():5.1f}%, >=48 kt {100 * (v >= 48).mean():5.1f}%, "
        f">=64 kt {100 * (v >= 64).mean():5.1f}%; gust reported {len(gg)} (median {gg.median() if len(gg) else float('nan'):.1f})")
x = m[(m.grp == "CF")].dropna(subset=["wind"])
if len(x) > 10:
    say(f"CF: Spearman rank correlation, ERA5 index vs station mean wind: {x.gmax.rank().corr(x.wind.rank()):.2f} (n={len(x)})")
    wd = x.wdir.dropna()
    say(f"CF: station wind direction from W-NW (240-330) {100 * ((wd >= 240) & (wd <= 330)).mean():.1f}%, "
        f"from N-E (330-120) {100 * ((wd > 330) | (wd < 120)).mean():.1f}% (n={len(wd)})")
open(os.path.join(HERE, "pcs-result.txt"), "w").write("\n".join(out) + "\n")
