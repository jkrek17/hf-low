"""Test 3 stage B selection: 250 HF + 250 SF per basin from the 400+400 already pulled; 6 new lags (-48..-12, +12, +24 h) around the peak time.
Position = linear interpolation of the 00/12 UTC track positions (intensity/results/fixes_2004.csv.gz) to the lag time; a lag outside the track span is missing (counted).
usage: select_lc.py REPO_ROOT. Writes results/lc_sample.csv, results/times_lc.csv, results/lc_selection_summary.txt. Labels and positions only; no field read."""
import sys, os, numpy as np, pandas as pd
root = sys.argv[1]; E5 = os.path.join(root, "research/era5"); OUT = os.path.join(E5, "hf_lifecycle/results")
A = pd.read_csv(os.path.join(E5, "hf_vs_storm/results/storm/anchors.csv"))
A = A[A.anchor.isin(["HF_peak", "SF_peak"])]
rng = np.random.default_rng(20261012)
parts = []
for (b, a), g in A.groupby(["basin", "anchor"]):
    parts.append(g.sample(250, random_state=20261012))
S = pd.concat(parts)[["track", "basin", "anchor", "time", "lat", "lon", "heading", "peak_time", "season", "mon", "g800"]].copy()
S["grp"] = np.where(S.anchor == "HF_peak", 0, 1)
S.to_csv(os.path.join(OUT, "lc_sample.csv"), index=False)
F = pd.read_csv(os.path.join(E5, "intensity/results/fixes_2004.csv.gz"), usecols=["track", "time", "lat", "lon"])
F["tm"] = pd.to_datetime(F.time.astype(str), format="%Y%m%d%H")
F = F[F.track.isin(S.track)].sort_values(["track", "tm"])
byt = {t: g for t, g in F.groupby("track")}
LAGS = [-48, -36, -24, -12, 12, 24]
rows, miss = [], {l: 0 for l in LAGS}
for r in S.itertuples():
    g = byt.get(r.track)
    pk = pd.to_datetime(str(int(r.peak_time)), format="%Y%m%d%H") if not isinstance(r.peak_time, str) else pd.to_datetime(r.peak_time, format="%Y%m%d%H")
    for L in LAGS:
        t = pk + pd.Timedelta(hours=L)
        if g is None or len(g) < 2 or t < g.tm.iloc[0] or t > g.tm.iloc[-1]:
            miss[L] += 1; continue
        x = (g.tm - g.tm.iloc[0]).dt.total_seconds().values; xt = (t - g.tm.iloc[0]).total_seconds()
        lon = np.degrees(np.unwrap(np.radians(g.lon.values)))
        rows.append(dict(track=r.track, anchor=f"lag{L:+d}", time=t.strftime("%Y%m%d%H"), lat=float(np.interp(xt, x, g.lat.values)), lon=float(np.interp(xt, x, lon)) % 360,
                         heading=r.heading, basin=r.basin))
T = pd.DataFrame(rows); T.to_csv(os.path.join(OUT, "times_lc.csv"), index=False)
txt = [f"sample {len(S)} storms ({S.groupby(['basin','anchor']).size().to_dict()}); lags with a 00/12 track bracket: {len(T)} of {len(S)*len(LAGS)}; missing by lag {miss}; distinct times {T.time.nunique()}"]
open(os.path.join(OUT, "lc_selection_summary.txt"), "w").write("\n".join(txt) + "\n"); print(txt[0])
