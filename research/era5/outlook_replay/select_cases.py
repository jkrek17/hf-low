"""Replay cases, fixed in PREREGISTRATION.md: per season the 2 cycles with the most OPC HF fixes in the following 96 h, 1 'no HF called'
cycle with the highest day-2 count of points in the two top gust bands, and 1 seeded random cycle. usage: select_cases.py"""
import os, glob, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, "results")
d = pd.read_csv(os.path.join(RES, "point_days.csv"), dtype={"cycle": str})
cyc = sorted(d.cycle.unique()); season = {c: (int(c[:4]) if int(c[4:6]) >= 10 else int(c[:4]) - 1) for c in cyc}
parts = []
for tag in ("Atl", "Pac"):
    a = pd.read_csv(f"/home/user/hf-low/data/hf_lows/HF_Data_-_{tag}.csv", dtype={"date": str})
    a = a[a.date.str.fullmatch(r"\d{10}", na=False)].copy(); a["t"] = pd.to_datetime(a.date, format="%Y%m%d%H", errors="coerce"); parts.append(a.dropna(subset=["t"]))
A = pd.concat(parts); H = A[A.Category == "HF"]
nhf = {}
for c in cyc:
    t0 = pd.Timestamp(c); nhf[c] = int(((H.t > t0) & (H.t <= t0 + pd.Timedelta(hours=96))).sum())
top = d[(d.day == 2) & (d.band >= 5)].groupby("cycle").n.sum().reindex(cyc).fillna(0)
rng = np.random.default_rng(20261014); rows = []
for s in sorted(set(season.values())):
    cs = [c for c in cyc if season[c] == s]; chosen = []
    # most HF fixes, cycles at least 4 days apart so the two are different storms
    for c in sorted(cs, key=lambda c: -nhf[c]):
        if all(abs((pd.Timestamp(c) - pd.Timestamp(x)).days) >= 4 for x in chosen) and len(chosen) < 2:
            chosen.append(c); rows.append((s, c, "most HF called", nhf[c], int(top[c])))
    z = [c for c in cs if nhf[c] == 0 and c not in chosen]
    c = max(z, key=lambda c: top[c]); chosen.append(c); rows.append((s, c, "none called", 0, int(top[c])))
    rest = [c for c in cs if c not in chosen]; c = rng.choice(rest); rows.append((s, str(c), "random", nhf[str(c)], int(top[str(c)])))
R = pd.DataFrame(rows, columns=["season", "cycle", "why", "opc_hf_fixes_96h", "gfs_day2_points_ge72kt"]); R.to_csv(os.path.join(RES, "replay_cases.csv"), index=False); print(R.to_string())
