"""Independent recomputation of tools/storm_review_facts.csv from a raw fixes CSV.
Usage: python3 -I verify_numbers.py FIXES_CSV FACTS_CSV [STORMS_CSV]
STORMS_CSV is optional; if given, rank/minp ordering is checked against it.
sw_max is NOT checked (needs ERA5 fields). hf_lat/hf_lon are checked against
several candidate definitions because the task did not define them."""
import sys
import numpy as np
import pandas as pd

TOL = {"minp": 0.05, "p_first": 0.05, "lat": 0.01, "lon": 0.01, "gust_idx": 0.05,
       "dp24_max": 0.05, "bomb_thr": 1e-6, "lat_dp24": 1e-6}
HF = 71.7


def wlon(x):
    return ((x + 180.0) % 360.0) - 180.0


def compute(g):
    g = g.sort_values("t").reset_index(drop=True)
    o = {}
    i = int(g.msl.values.argmin())
    o["t_minp"] = g.t[i]; o["minp"] = g.msl[i]
    o["lat_minp"] = g.lat[i]; o["lon_minp"] = wlon(g.lon[i])
    o["n_minp_ties"] = int((g.msl == g.msl.min()).sum())
    gg = g[g.g800 > 0]
    if len(gg):
        j = gg.g800.idxmax()
        o["t_gust"] = g.t[j]; o["gust_idx"] = g.g800[j]
    hf = g[g.g800 >= HF]
    o["n_hf"] = len(hf)
    if len(hf):
        o["first_hf"] = hf.t.min(); o["last_hf"] = hf.t.max()
        # candidate hf positions
        k = hf.g800.idxmax()
        o["hf_cand"] = {
            "first_hf": (hf.lat.iloc[0], wlon(hf.lon.iloc[0])),
            "last_hf": (hf.lat.iloc[-1], wlon(hf.lon.iloc[-1])),
            "max_gust_hf": (g.lat[k], wlon(g.lon[k])),
            "mean_hf": (hf.lat.mean(), wlon(np.degrees(np.arctan2(np.sin(np.radians(hf.lon)).mean(), np.cos(np.radians(hf.lon)).mean())))),
            "mean_hf_arith360": (hf.lat.mean(), wlon(hf.lon.mean())),
            "median_hf": (hf.lat.median(), wlon(hf.lon.median())),
            "minp_hf": (hf.lat[hf.msl.idxmin()], wlon(hf.lon[hf.msl.idxmin()])),
        }
    o["first_fix"] = g.t.iloc[0]; o["last_fix"] = g.t.iloc[-1]; o["p_first"] = g.msl.iloc[0]
    m = g.set_index("t").msl
    best = None
    for idx in range(len(g)):
        t = g.t[idx]; t0 = t - pd.Timedelta(hours=24)
        if t0 in m.index:
            d = m[t0] - g.msl[idx]
            if best is None or d > best[0]:
                best = (d, idx)
    if best:
        d, idx = best
        lm = g.lat.iloc[max(0, idx - 4):idx + 1].mean()
        o["dp24_max"] = d; o["t_dp24"] = g.t[idx]; o["lat_dp24"] = lm
        o["bomb_thr"] = 24 * np.sin(np.radians(lm)) / np.sin(np.radians(60))
    o["n_fixes"] = len(g)
    gaps = g.t.diff().dropna().dt.total_seconds() / 3600
    o["irregular_steps"] = int((gaps != 6).sum())
    return o


def main():
    fixes = pd.read_csv(sys.argv[1]); facts = pd.read_csv(sys.argv[2])
    storms = pd.read_csv(sys.argv[3]) if len(sys.argv) > 3 else None
    fixes["t"] = pd.to_datetime(fixes.time.astype(str), format="%Y%m%d%H")
    for c in ("t_minp", "t_gust", "first_hf", "last_hf", "first_fix", "last_fix", "t_dp24"):
        facts[c] = pd.to_datetime(facts[c])
    cols = ["t_minp", "minp", "lat_minp", "lon_minp", "t_gust", "gust_idx", "first_hf", "last_hf",
            "n_hf", "first_fix", "last_fix", "p_first", "dp24_max", "t_dp24", "lat_dp24", "bomb_thr"]
    tolkey = {"minp": "minp", "p_first": "p_first", "lat_minp": "lat", "lon_minp": "lon",
              "gust_idx": "gust_idx", "dp24_max": "dp24_max", "lat_dp24": "lat_dp24", "bomb_thr": "bomb_thr"}
    ok = {c: 0 for c in cols}; tot = {c: 0 for c in cols}
    hfc = {}; mism = []; odd = []
    for _, r in facts.iterrows():
        g = fixes[fixes.track == r.track]
        if g.empty:
            mism.append((r.track, "ALL", "no fixes", "")); continue
        o = compute(g)
        if o["n_minp_ties"] > 1: odd.append(f"track {r.track}: {o['n_minp_ties']} fixes tie for lowest msl (first used)")
        if o["irregular_steps"]: odd.append(f"track {r.track}: {o['irregular_steps']} non-6h steps in fixes")
        if o["n_hf"] == 0: odd.append(f"track {r.track}: no HF fix")
        if o.get("dp24_max", 0) < 0: odd.append(f"track {r.track}: negative dp24_max")
        if o.get("dp24_max") is not None and o["dp24_max"] < o.get("bomb_thr", 0): odd.append(f"track {r.track}: dp24_max {o['dp24_max']:.2f} < bomb_thr {o['bomb_thr']:.2f}")
        for c in cols:
            tot[c] += 1
            a, b = r[c], o.get(c)
            if b is None or pd.isna(a):
                good = (b is None) and pd.isna(a)
            elif c in tolkey:
                good = abs(float(a) - float(b)) <= TOL[tolkey[c]]
            elif c == "n_hf":
                good = int(a) == int(b)
            else:
                good = pd.Timestamp(a) == pd.Timestamp(b)
            if good: ok[c] += 1
            else: mism.append((r.track, c, a, b))
        if "hf_cand" in o and not pd.isna(r.hf_lat):
            for name, (la, lo) in o["hf_cand"].items():
                h = hfc.setdefault(name, [0, 0])
                h[1] += 1
                if abs(la - r.hf_lat) <= 0.01 and abs(wlon(lo - r.hf_lon)) <= 0.01: h[0] += 1
    n = len(facts)
    print(f"storms in facts: {n}")
    for c in cols: print(f"{c:10s} matched {ok[c]}/{tot[c]}")
    print("hf_lat/hf_lon candidate definitions (both within 0.01 deg):")
    for k, (a, b) in hfc.items(): print(f"  {k:18s} {a}/{b}")
    print("MISMATCHES:" if mism else "MISMATCHES: none")
    for m in mism: print("  ", m)
    # ordering
    print("ORDER CHECK (facts rank vs minp per basin):")
    for b, d in facts.groupby("basin"):
        d = d.sort_values("rank")
        mp = [compute(fixes[fixes.track == t].copy())["minp"] for t in d.track]
        print(f"  {b}: ranks {list(d['rank'])} minp {[round(x,1) for x in mp]} ascending={all(np.diff(mp) >= 0)}")
        if not all(np.diff(mp) >= 0): odd.append(f"{b}: recomputed minp not ascending with rank")
    if storms is not None:
        s = set(storms.track); f = set(facts.track)
        print(f"storm list vs facts tracks equal: {s == f}; only-in-list {sorted(s-f)}; only-in-facts {sorted(f-s)}")
        mg = storms.merge(facts, on="track", suffixes=("_s", "_f"))
        print(f"basin match {(mg.basin_s==mg.basin_f).sum()}/{len(mg)}; rank match {(mg.rank_s==mg.rank_f).sum()}/{len(mg)}")
    print("ODD:" if odd else "ODD: none")
    for x in odd: print("  ", x)
    print("sw_max: NOT CHECKED")


main()
