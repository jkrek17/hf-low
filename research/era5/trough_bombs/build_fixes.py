"""Link the pressure-only lows into pipeline A tracks and build the 00/12 UTC forecast fixes (1979-2003 seasons).

ERA5 proxy, pipeline A. Outcome-free apart from the deepening rate that defines BOMB, which is computed here
(`ndr24`) but not read by anything before the pre-registered tracker check.

usage: build_fixes.py LOWS_DIR OUT_FIXES.csv.gz [--check]
  --check   only link and compare with the committed 1979-2003 catalog track points
            (`hf_history/results/era5_hf_catalog_tracks.csv`); prints the recovery numbers, writes nothing.
Needs the same linker as the committed catalog (`hf_history/track.py`: <= 900 km per 6 h, tracks >= 24 h).
"""
import os, sys, glob
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "hf_history"))
sys.path.insert(0, os.path.join(HERE, "..", "intensity"))
import track as TR
import fixes as FX


def tracks(lows_dir):
    files = sorted(glob.glob(f"{lows_dir}/*.csv"))
    df = pd.concat([pd.read_csv(f, dtype={"time": str}) for f in files], ignore_index=True)
    df = df[(df.lat >= 20) & (df.lat <= 75)]
    df = TR.link(df)
    n = df.groupby("track").size()
    df = df[df.track.isin(n[n >= TR.MIN_LEN].index)].copy()
    df["basin"] = [TR.basin(a, b) for a, b in zip(df.lat, df.lon)]
    df["season"] = [TR.season(t) for t in df.time]
    return df


def check(df):
    C = pd.read_csv(os.path.join(HERE, "..", "hf_history", "results", "era5_hf_catalog_tracks.csv"), dtype={"time": str})
    C = C[C.time < "2004060100"]
    m = C.merge(df[["time", "lat", "lon", "msl", "track"]], on="time", how="left", suffixes=("", "_r"))
    m = m[(m.lat_r.isna()) | ((abs(m.lat - m.lat_r) < 0.01) & (abs(m.lon - m.lon_r) < 0.01))]
    m = m.drop_duplicates(["track", "time"])
    n = len(C)
    found = m.lat_r.notna()
    same = found & ((m.msl - m.msl_r).abs() <= 0.2)
    print(f"catalog points before 2004-06: {n}; position found {found.sum()} ({found.mean():.4%}); "
          f"same pressure within 0.2 hPa {same.sum()} ({same.mean():.4%})")
    # does the re-linked track keep the catalog points together?
    g = m[found].groupby("track").track_r.agg(lambda s: s.nunique())
    print(f"catalog tracks {C.track.nunique()}; catalog tracks mapped to exactly one re-linked track: {(g == 1).sum()}")
    miss = C.merge(m[["track", "time"]].assign(ok=1), on=["track", "time"], how="left")
    print("missing points by season:", miss[miss.ok.isna()].time.str[:4].value_counts().sort_index().to_dict())


if __name__ == "__main__":
    lows, out = sys.argv[1], sys.argv[2]
    df = tracks(lows)
    print(len(df), "fixes", df.track.nunique(), "tracks")
    if "--check" in sys.argv:
        check(df)
        sys.exit()
    if len(sys.argv) > 3 and sys.argv[3].endswith(".csv.gz"):
        df.to_csv(sys.argv[3], index=False)              # linked points, kept in work/ for re-use
    df["g800"] = 0.0                      # gust is not used before 2004
    F = FX.build(df[["track", "time", "basin", "season", "lat", "lon", "msl", "g800"]])
    F = F.drop(columns=["g800", "hf_now", "hf24", "hf48", "left24"])
    F.to_csv(out, index=False)
    print(len(F), "forecast fixes", F.track.nunique(), "tracks; seasons", F.season.min(), "-", F.season.max())
