"""Step 1 (PREREGISTRATION.md): GFS-minus-ERA5 low position error by lead. ERA5 PROXY, pipeline A.
Usage: python3 step1_position_error.py   (run from this folder; reads two committed tables, writes results/)"""
import os, sys, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
GFS = os.path.join(ROOT, "research/era5/gfs_raw_gust_96h/results/raw_gust_fixes_96h.csv.gz")
ERA = os.path.join(ROOT, "research/era5/intensity/results/fixes_2004.csv.gz")
OUT = os.path.join(HERE, "results"); os.makedirs(OUT, exist_ok=True)
GATE_KM = 1200.0; NM = 1.852
g = pd.read_csv(GFS)
e = pd.read_csv(ERA, usecols=["time", "lat", "lon", "msl"])
g["valid"] = g.valid_time.astype(np.int64); e["time"] = e.time.astype(np.int64)
eg = {t: d[["lat", "lon"]].to_numpy() for t, d in e.groupby("time")}
def hav(la1, lo1, la2, lo2):
    p = np.pi / 180; a = np.sin((la2 - la1) * p / 2) ** 2 + np.cos(la1 * p) * np.cos(la2 * p) * np.sin((lo2 - lo1) * p / 2) ** 2
    return 2 * 6371.0 * np.arcsin(np.sqrt(a))
rows = []
for t, d in g.groupby("valid"):
    E = eg.get(t)
    for r in d.itertuples(index=False):
        if E is None or len(E) == 0:
            rows.append((r.lead, r.basin, r.G, np.nan, np.nan, np.nan)); continue
        dist = hav(r.lat, r.lon, E[:, 0], E[:, 1]); k = int(np.argmin(dist))
        if dist[k] > GATE_KM:
            rows.append((r.lead, r.basin, r.G, np.nan, np.nan, np.nan)); continue
        dlon = (r.lon - E[k, 1] + 180) % 360 - 180
        rows.append((r.lead, r.basin, r.G, r.lat - E[k, 0], dlon, dist[k] / NM))
m = pd.DataFrame(rows, columns=["lead", "basin", "G", "dlat", "dlon", "dist_nm"])
m["matched"] = m.dist_nm.notna()
m.to_csv(os.path.join(OUT, "matches.csv.gz"), index=False)
out = []
for grp, sel in (("G>=64kt (primary)", m.G >= 64), ("all fixes", m.G >= -1)):
    for basin in ("atl", "pac", "both"):
        s = m[sel & ((m.basin == basin) if basin != "both" else True)]
        for L, d in s.groupby("lead"):
            x = d[d.matched]
            q = x.dist_nm.quantile([.5, .75, .9]).values if len(x) else [np.nan] * 3
            out.append(dict(group=grp, basin=basin, lead=L, n=len(d), matched=len(x), matched_share=round(len(x) / len(d), 3),
                            median_nm=round(q[0], 1), p75_nm=round(q[1], 1), p90_nm=round(q[2], 1), mean_nm=round(x.dist_nm.mean(), 1),
                            mean_dlat=round(x.dlat.mean(), 3), mean_dlon=round(x.dlon.mean(), 3)))
t = pd.DataFrame(out); t.to_csv(os.path.join(OUT, "position_error_by_lead.csv"), index=False)
vec = m[(m.G >= 64) & m.matched][["lead", "basin", "dlat", "dlon"]]
vec.to_csv(os.path.join(OUT, "displacement_vectors_primary.csv"), index=False)
with open(os.path.join(OUT, "position_error.txt"), "w") as f:
    f.write("GFS minus ERA5 low position, nautical miles, matched within 1200 km at the same valid time (ERA5 PROXY, pipeline A)\n")
    f.write(t.to_string(index=False))
print(t[t.group.str.startswith("G")][["basin", "lead", "n", "matched_share", "median_nm", "p75_nm", "p90_nm"]].to_string(index=False))
