"""Step 0: count candidate (parent, daughter) pairs and nothing else. No orientation, no null, no outcome."""
import numpy as np, pandas as pd
from common import *
T, F = load()
print("HF tracks", len(T), T.groupby("basin").size().to_dict())
print("gen_obs", T.groupby("basin").gen_obs.mean().round(3).to_dict(), " first fix within 12 h of start:", (T.lag0 <= 12).mean().round(3))
rows = []
for b, G in T.groupby("basin"):
    G = G.reset_index(drop=True)
    for j in range(len(G)):
        d = G.iloc[j]
        cand = np.where((G.t0.values < d.t0) & (G.t1.values >= d.t0) & (np.abs(d.t0 - G.tpeak.values) <= 48))[0]
        for i in cand:
            if i == j: continue
            p = G.iloc[i]
            pp = pos_at(F[p.track], d.t0)
            if pp is None: continue
            dist = hav(pp[0], pp[1], d.lat0, d.lon0)
            if dist <= 1500:
                rows.append((b, p.track, d.track, dist, d.gen_obs, p.gen_obs, d.lag0 <= 12))
R = pd.DataFrame(rows, columns="basin parent daughter dist dobs pobs dlag".split())
print("candidate pairs (P alive at D first fix, |t_D - peak_P| <= 48 h, <= 1500 km):")
print(R.groupby("basin").size().to_dict())
for lo in (0, 400):
    r = R[R.dist >= lo]; print(f" dist >= {lo} km:", r.groupby("basin").size().to_dict(), " D observed genesis:", r[r.dobs == True].groupby("basin").size().to_dict(),
          " D obs & lag<=12:", r[(r.dobs == True) & r.dlag].groupby("basin").size().to_dict())
print("daughters with >=1 parent:", R.groupby("basin").daughter.nunique().to_dict())
R.to_csv("work_candidates.csv", index=False)
