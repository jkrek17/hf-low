"""Independent recomputation of late_highlat numbers (RA-26). Own code; event table follows RA-7 filters."""
import os, numpy as np, pandas as pd
ERA = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
L = pd.read_csv(f"{ERA}/hf_history/results/lifecycle_events.csv")
M = L[(L.era == "2004+") & (~L.tc.astype(bool))].copy()
f = pd.read_csv(f"{ERA}/hf_structure/results/fixes.csv").sort_values(["track", "time"], kind="stable")
fo = f.drop_duplicates("track", keep="first").copy()          # earliest ROW per track
fo_first = f.groupby("track").first().reset_index()           # RA-7 style: first non-null per column
T = pd.read_csv(f"{ERA}/hf_history/results/era5_hf_catalog_tracks.csv")
C = pd.read_csv(f"{ERA}/hf_history/results/era5_hf_catalog.csv")
print("catalog duplicate tracks:", C.track.duplicated().sum(), "| fixes onset rows differing (row vs groupby.first) in gmax_r/gmax_lat/gmax_coast_km/terrain:",
      int((fo.set_index("track")[["gmax_r","gmax_lat","gmax_coast_km","terrain"]].astype(float).fillna(-999)
           != fo_first.set_index("track")[["gmax_r","gmax_lat","gmax_coast_km","terrain"]].astype(float).fillna(-999)).any(axis=1).sum()))
T = T[T.track.isin(M.track) & T.basin.notna() & (T.g800 >= 71.7)].copy()
T["t"] = pd.to_datetime(T.time.astype(str), format="%Y%m%d%H")
ton = T.groupby("track").t.min().rename("t_on").reset_index()

def build(fo):
    D = fo[["track","time","gmax_r","gmax_lat","gmax_coast_km","terrain"]].copy()
    D["t_on_f"] = pd.to_datetime(D.time.astype(str), format="%Y%m%d%H")
    X = M.merge(D, on="track", how="left").merge(ton, on="track", how="left").reset_index(drop=True)
    bad = (X.t_on_f != X.t_on) | X.gmax_r.isna()
    print(f"events {len(X)}; dropped (onset mismatch or no gmax_r): {int(bad.sum())}")
    X = X[~bad].copy()
    X["late"] = (X.h_on_minp >= 0).astype(int); X["late6"] = (X.h_on_minp >= 6).astype(int)
    X["terr"] = X.terrain.astype(bool)
    X["coast"] = X.gmax_coast_km <= 300; X["far"] = X.gmax_r > 400; X["near"] = X.gmax_r <= 400
    X["barrier"] = X.coast & X.far
    return X
E = build(fo)
E2 = build(fo_first)
A = E[E.basin == "atl"]
print("Atlantic events:", len(A), "| gmax_coast_km NaN:", int(A.gmax_coast_km.isna().sum()), "| gmax_lat NaN:", int(A.gmax_lat.isna().sum()), "| h_on_minp NaN:", int(A.h_on_minp.isna().sum()))
H = A[A.gmax_lat >= 60].copy(); S = A[A.gmax_lat < 60].copy()
# gust direction join
G = pd.read_csv(f"{ERA}/highlat/gustloc_fixes.csv")[["track","time","max_wdir"]]
print("gustloc duplicate track+time:", int(G.duplicated(["track","time"]).sum()))
H = H.merge(G, on=["track","time"], how="left")
H["north"] = np.where(H.max_wdir.isna(), np.nan, ((H.max_wdir >= 315) | (H.max_wdir <= 45)).astype(float))

def sh(d, col, lab):
    l, n = d[d.late == 1], d[d.late == 0]
    return l, n
def row(name, d, col, latecol="late"):
    d = d.dropna(subset=[col]); l = d[d[latecol] == 1][col].astype(float); n = d[d[latecol] == 0][col].astype(float)
    print(f"{name:42s} late {l.mean():.4f} ({int(l.sum())}/{len(l)})  non-late {n.mean():.4f} ({int(n.sum())}/{len(n)})  diff {l.mean()-n.mean():+.4f}")
    return l.mean() - n.mean()
print("\n(a) group n=%d LATE=%d nonlate=%d LATE6=%d terrain=%d (late&terrain=%d)" % (len(H), H.late.sum(), (1-H.late).sum(), H.late6.sum(), H.terr.sum(), (H.late.astype(bool)&H.terr).sum()))
print("\n(b)"); p1 = row("P1 BARRIER", H, "barrier")
print("(c)"); row("coast<=300", H, "coast"); row("gmax_r>400", H, "far"); row("gmax_r<=400", H, "near")
print("(d)"); row("NORTHERLY (with direction)", H, "north"); print("   direction coverage:", int(H.max_wdir.notna().sum()), "of", len(H), "| late", int(H[H.late==1].max_wdir.notna().sum()), "non-late", int(H[H.late==0].max_wdir.notna().sum()))
print("(e)"); row("BARRIER, LATE6 vs rest", H, "barrier", "late6")
print("(f)"); row("BARRIER, terrain-flagged removed", H[~H.terr], "barrier")
print("(g)"); row("BARRIER, Atlantic gmax_lat<60", S, "barrier")
# bootstrap
def boot(d, col, B=2000, seed=7):
    rng = np.random.default_rng(seed)
    seas = {s: g for s, g in d.groupby("season")}; ks = list(seas); out = []
    for _ in range(B):
        pick = rng.choice(len(ks), len(ks), replace=True)
        b = pd.concat([seas[ks[i]] for i in pick])
        l, n = b[b.late == 1][col], b[b.late == 0][col]
        if len(l) and len(n): out.append(l.mean() - n.mean())
    return np.percentile(out, [2.5, 97.5]), len(seas), len(out)
H["barrier"] = H.barrier.astype(float)
ci, ns, nb = boot(H, "barrier")
print(f"\nP1 season-block bootstrap: seasons={ns} draws={nb} 95% CI [{ci[0]:+.4f}, {ci[1]:+.4f}] (seed 7, numpy default_rng)")
for sd in (1, 2, 123):
    c2, _, _ = boot(H, "barrier", seed=sd); print(f"   seed {sd}: [{c2[0]:+.4f}, {c2[1]:+.4f}]")
print("seasons in group:", sorted(H.season.unique()))
print("\nVariant check, RA-7 groupby.first onset table:")
A2 = E2[E2.basin == "atl"]; H2 = A2[A2.gmax_lat >= 60]
print(" group n", len(H2), "late", H2.late.sum(), "P1:", end=" "); row("BARRIER", H2, "barrier")
