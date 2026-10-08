"""Select the HF and storm-force-only storms and their anchor times (labels only, no field is read).

ERA5 proxy, pipeline A (research/era5/hf_history). Seasons 2004-05 to 2025-26.
  HF   : track peak gust index g800 >= 71.7 kt (catalog events).
  SF   : track peak gust index in [SF_LO, 71.7) kt.
Anchors: HF onset = first 6-hourly in-domain point with g800 >= 71.7 (catalog tracks);
         peak    = pipeline A peak time and position (all_tracks.csv.gz), HF and SF.
Heading (bearing of motion) = bearing between the latest 00/12 UTC fix at or before t-6 h and the earliest
at or after t+6 h of the same track (intensity/results/fixes_2004.csv.gz); one-sided at a track end.
Writes results/storms.csv (every HF and SF storm, strata, anchors) and results/times_storm.csv (pull list).
usage: select_times.py REPO_ROOT [HF_CAP_PER_BASIN]
"""
import sys, os, numpy as np, pandas as pd
root = sys.argv[1]
cap = int(sys.argv[2]) if len(sys.argv) > 2 else 10**9
RES = os.path.join(root, "research/era5/hf_history/results")
OUT = os.path.join(root, "research/era5/hf_vs_storm/results")
HF, SF_LO, SF_HI2 = 71.7, 54.0, 60.0
SEED = 20261012
A = pd.read_csv(os.path.join(RES, "all_tracks.csv.gz"))
A = A[(A.season >= 2004) & A.basin.isin(["atl", "pac"])].copy()
A["grp"] = np.where(A.gust800_kt >= HF, "HF", np.where(A.gust800_kt >= SF_LO, "SF", "none"))
A = A[A.grp != "none"].copy()
A["mon"] = (A.peak_time // 10000) % 100
C = pd.read_csv(os.path.join(RES, "era5_hf_catalog_tracks.csv"))
C = C[C.g800 >= HF].sort_values(["track", "time"]).groupby("track").first().reset_index()
on = C.set_index("track")[["time", "lat", "lon"]].rename(columns=dict(time="onset_time", lat="onset_lat", lon="onset_lon"))
A = A.merge(on, left_on="track", right_index=True, how="left")
A.loc[A.grp == "SF", ["onset_time", "onset_lat", "onset_lon"]] = np.nan
# heading from 00/12 fixes
F = pd.read_csv(os.path.join(root, "research/era5/intensity/results/fixes_2004.csv.gz"), usecols=["track", "time", "lat", "lon"])
F["tm"] = pd.to_datetime(F.time.astype(str), format="%Y%m%d%H")
F = F.sort_values(["track", "tm"])
byt = {k: g for k, g in F.groupby("track")}
def bearing(la1, lo1, la2, lo2):
    p1, p2, dl = np.radians(la1), np.radians(la2), np.radians(lo2 - lo1)
    return np.degrees(np.arctan2(np.sin(dl) * np.cos(p2), np.cos(p1) * np.sin(p2) - np.sin(p1) * np.cos(p2) * np.cos(dl))) % 360
def heading(tr, t):
    g = byt.get(tr)
    if g is None: return np.nan
    t = pd.to_datetime(str(int(t)), format="%Y%m%d%H")
    a = g[g.tm <= t - pd.Timedelta(hours=6)]; b = g[g.tm >= t + pd.Timedelta(hours=6)]
    if len(a) == 0 and len(b) == 0: return np.nan
    p = a.iloc[-1] if len(a) else g[g.tm >= t].iloc[0]
    q = b.iloc[0] if len(b) else g[g.tm <= t].iloc[-1]
    if p.tm == q.tm: return np.nan
    return bearing(p.lat, p.lon, q.lat, q.lon)
A["hd_peak"] = [heading(tr, t) for tr, t in zip(A.track, A.peak_time)]
A["hd_onset"] = [heading(tr, t) if pd.notna(t) else np.nan for tr, t in zip(A.track, A.onset_time)]
A["age_peak_h"] = (pd.to_datetime(A.peak_time.astype(str), format="%Y%m%d%H") - pd.to_datetime(A.start.astype(str), format="%Y%m%d%H")).dt.total_seconds() / 3600
A["sf_hi"] = (A.grp == "SF") & (A.gust800_kt >= SF_HI2)
# storm-scale sample: HF cap per basin (random within basin), SF drawn 1:1 within basin x month, without replacement
rng = np.random.default_rng(SEED)
A["pull"] = False
hf_ix = []
for b, g in A[A.grp == "HF"].groupby("basin"):
    hf_ix += list(rng.permutation(g.index.values)[:min(cap, len(g))])
A.loc[hf_ix, "pull"] = True
for (b, m), g in A[A.pull & (A.grp == "HF")].groupby(["basin", "mon"]):
    pool = A[(A.grp == "SF") & (A.basin == b) & (A.mon == m)].index.values
    A.loc[rng.permutation(pool)[:len(g)], "pull"] = True
A.to_csv(os.path.join(OUT, "storms.csv"), index=False)
rows = []
for r in A[A.pull].itertuples():
    if r.grp == "HF":
        rows.append((r.track, "HF_onset", int(r.onset_time), r.onset_lat, r.onset_lon, r.hd_onset, r.basin))
        rows.append((r.track, "HF_peak", int(r.peak_time), r.peak_lat, r.peak_lon, r.hd_peak, r.basin))
    else:
        rows.append((r.track, "SF_peak", int(r.peak_time), r.peak_lat, r.peak_lon, r.hd_peak, r.basin))
T = pd.DataFrame(rows, columns=["track", "anchor", "time", "lat", "lon", "heading", "basin"])
T.to_csv(os.path.join(OUT, "times_storm.csv"), index=False)
print(A.groupby(["basin", "grp"]).size().unstack())
print("pull storms:", A[A.pull].groupby(["basin", "grp"]).size().to_dict())
print("anchors:", T.anchor.value_counts().to_dict(), "distinct times:", T.time.nunique(), "NaN heading:", T.heading.isna().sum())
