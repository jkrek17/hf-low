"""POST HOC checks added after the claim review (not registered): per-basin AUC, Atlantic split at 60N, candidate-set size, sample enrichment.
usage: post_hoc.py [footprints.csv.gz]   (default: results/scat_footprints.csv.gz)"""
import os, sys, numpy as np, pandas as pd
import analyse as A

here = os.path.dirname(os.path.abspath(__file__))
F = pd.read_csv(sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "results", "scat_footprints.csv.gz"))
T = pd.read_csv(os.path.join(here, "..", "hf_history", "results", "all_tracks.csv.gz"), dtype={"peak_time": str})
T = T[(T.season >= 2004) & T.peak_time.str[-2:].isin(["00", "12"]) & (T.gust800_kt > 0)]
print("candidate storms (00/12 UTC peak, stored index > 0): from 2004-05", len(T), " from 2006-07", int((T.season >= 2006).sum()))
G = F[F.cover >= 0.40].sort_values("cover", ascending=False).drop_duplicates("track").merge(T[["track", "peak_lat"]], on="track")
print("analysed storms", len(G), " share of 2006-07-on candidates", round(len(G) / int((T.season >= 2006).sum()), 3))
print("index >= 71.7: analysed", round(float((G.idx >= 71.7).mean()), 3), " candidates (2006-07 on)", round(float((T[T.season >= 2006].gust800_kt >= 71.7).mean()), 3))
print("passes in P4 set (cover >= 0.40):", int((F.cover >= 0.40).sum()), "from storms", F[F.cover >= 0.40].track.nunique())
print("\nPer-basin AUC (>=45 kt), season-block 95% CI:")
for b, g in G.groupby("basin"):
    bs = A.boot_seasons(g, A.auc)
    print(f"  {b}: {A.auc(g):.3f} [{np.percentile(bs, 2.5):.3f}, {np.percentile(bs, 97.5):.3f}] n={len(g)}")
print("\nAtlantic Spearman rho (footprint maxima) split at 60N:")
atl = G[G.basin == "atl"]
for name, g in (("60N and south", atl[atl.peak_lat <= 60]), ("north of 60N", atl[atl.peak_lat > 60])):
    bs = A.boot_seasons(g, A.rho)
    print(f"  {name}: {A.rho(g):.3f} [{np.percentile(bs, 2.5):.3f}, {np.percentile(bs, 97.5):.3f}] n={len(g)}  index>=71.7: {int((g.idx >= 71.7).sum())}")
print("\nP3 (index>=71.7 minus depth-matched), Atlantic south of 60N only (north of 60N dropped from the Atlantic):")
H = pd.concat([G[G.basin == "pac"], atl[atl.peak_lat <= 60]])
H = A.add_pdiff(H)
bs = A.boot_seasons(H, A.pdiff_mean)
print(f"  {A.pdiff_mean(H):.2f} kt [{np.nanpercentile(bs, 2.5):.2f}, {np.nanpercentile(bs, 97.5):.2f}] n={int(H.pdiff.notna().sum())}")
