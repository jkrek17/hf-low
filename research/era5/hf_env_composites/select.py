"""Storm times for the HF-low composites (labels only; no fields are read).

usage: select.py FIXES_CSV_GZ OUT_CSV
FIXES_CSV_GZ is ../intensity/results/fixes_2004.csv.gz. One row per (track, anchor): anchor in
onset, peak, m12, m24, m48 (the same track's 00/12 UTC fix 12/24/48 h before onset, if it exists).
Population and definitions are in PREREGISTRATION.md.
"""
import sys
import numpy as np, pandas as pd

HF = 71.7
MONTHS = (10, 11, 12, 1, 2, 3, 4)
LAST = pd.Timestamp("2023-01-10")        # first time not in WeatherBench2


def build(F):
    F = F.copy()
    F["t"] = pd.to_datetime(F.time, format="%Y%m%d%H")
    hf = F[F.g800 >= HF].sort_values("t")
    on = hf.groupby("track").first().reset_index()
    pk = F.loc[hf.groupby("track").g800.idxmax()].set_index("track")
    ok = on.t.dt.month.isin(MONTHS) & (on.t < LAST)
    ok &= on.track.map(pk.t) < LAST
    on = on[ok]
    idx = F.set_index(["track", "t"])
    rows = []
    for r in on.itertuples():
        for anchor, t in [("onset", r.t), ("peak", pk.t[r.track]), ("m12", r.t - pd.Timedelta(hours=12)),
                          ("m24", r.t - pd.Timedelta(hours=24)), ("m48", r.t - pd.Timedelta(hours=48))]:
            if (r.track, t) not in idx.index:
                continue
            x = idx.loc[(r.track, t)]
            rows.append(dict(key=f"{r.track}_{anchor}", track=r.track, basin=r.basin, season=int(r.season),
                             anchor=anchor, time=t.strftime("%Y%m%d%H"), lat=x.lat, lon=x.lon,
                             heading=x.heading, g800=x.g800, onset_time=r.t.strftime("%Y%m%d%H")))
    return pd.DataFrame(rows)


if __name__ == "__main__":
    F = pd.read_csv(sys.argv[1], dtype={"time": str})
    D = build(F)
    D.to_csv(sys.argv[2], index=False)
    print(len(D), "rows;", D[D.anchor == "onset"].groupby("basin").size().to_dict(), "HF lows")
    print(D.groupby(["basin", "anchor"]).size().unstack())
    print("heading missing:", D.groupby("anchor").heading.apply(lambda s: int(s.isna().sum())).to_dict())
