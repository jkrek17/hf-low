"""POST HOC (not in the plan): among pipeline A Atlantic HF events (tele_intensity event table), does the lagged Greenland
high raise minimum pressure through the background (ring - clim) or through depth against the surroundings? ERA5 proxy.
usage: posthoc.py ALL_TRACKS CPC_DIR REPO_ROOT GH_DAILY FIXES ENV EVENT_TABLE OUT_TXT"""
import sys
import numpy as np, pandas as pd
import sg

tracks, cpc, repo, ghfile, fixes, env, evfile, out = sys.argv[1:9]
P = sg.prep(tracks, cpc, repo, ghfile, fixes, env)
W = P["Tw"]
E = pd.read_csv(evfile)
E = E[E.basin == "atl"]
M = W.merge(E[["track", "minp", "clim", "ring", "depth_ring", "bg", "anom_clim"]].rename(columns={"minp": "minp_ev"}), on="track", how="inner")
M = M[M.hf == 1].reset_index(drop=True)
mons = sorted(M.mon.unique())[1:]
Z = np.column_stack([(M.mon.values == m).astype(float) for m in mons] + [(M.season.values - M.season.mean()) / 10.0, np.ones(len(M))])
X = np.column_stack([M.zGH, M.zNAO, Z])
seasons = np.arange(sg.S0, sg.S1 + 1)
idx = {s: np.where(M.season.values == s)[0] for s in seasons}
rng = np.random.default_rng(sg.SEED)
outs = ["minp_ev", "bg", "depth_ring", "anom_clim", "gust800_kt"]
beta = lambda Xm, y: np.linalg.lstsq(Xm, y, rcond=None)[0]
est = {o: beta(X, M[o].values.astype(float))[0] for o in outs}
estn = {o: beta(X, M[o].values.astype(float))[1] for o in outs}
bs = {o: [] for o in outs}
for _ in range(2000):
    r = np.concatenate([idx[s] for s in rng.choice(seasons, len(seasons))])
    for o in outs:
        bs[o].append(beta(X[r], M[o].values.astype(float)[r])[0])
with open(out, "w") as f:
    f.write(f"POST HOC. Atlantic HF events (pipeline A, ERA5 proxy) with background pressure, 2004-05..2025-26: n = {len(M)} (HF events in the Oct-Apr window with event-table rows).\n")
    f.write("Effect of +1 SD lagged Greenland high, NAO in the model, month and trend terms, season-block bootstrap (2,000). Conditional on HF, so selection applies.\n")
    f.write("minp = track minimum pressure (hPa); bg = ring - clim (hPa, background anomaly at the deepest fix); depth_ring = minp - ring (hPa, negative = deeper than surroundings); anom_clim = minp - clim.\n\n")
    for o in outs:
        lo, hi = np.percentile(bs[o], [2.5, 97.5])
        f.write(f"{o:12s} GH|NAO {est[o]:+.3f} [{lo:+.3f}, {hi:+.3f}]   (NAO|GH {estn[o]:+.3f})\n")
print(open(out).read())
