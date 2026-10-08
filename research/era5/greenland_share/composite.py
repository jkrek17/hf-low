"""Descriptive composites for the Greenland-high share study (no tests): storm quantities by GH tercile within NAO tercile,
and a fix-density map. ERA5 PROXY, pipeline A. usage: composite.py ALL_TRACKS CPC_DIR REPO_ROOT GH_DAILY FIXES ENV OUT_DIR"""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sg


def main():
    tracks, cpc, repo, ghfile, fixes, env, out = sys.argv[1:8]
    os.makedirs(out, exist_ok=True)
    P = sg.prep(tracks, cpc, repo, ghfile, fixes, env)
    W, F = P["Tw"], P["F"]
    W = W[W.lat0.notna()].copy()
    W["tg"] = pd.qcut(W.zGH, 3, labels=[1, 2, 3]).astype(int)
    W["tn"] = pd.qcut(W.zNAO, 3, labels=[1, 2, 3]).astype(int)
    cols = {"n": ("track", "size"), "HF share %": ("hf", lambda x: 100 * x.mean()), "lat0": ("lat0", "mean"), "lon0": ("lon0", "mean"),
            "latm": ("latm", "mean"), "lonm": ("lonm", "mean"), "minp": ("minp", "mean"), "gust": ("gust800_kt", "mean"),
            "ndr_max": ("ndr_max", "mean"), "box %": ("box", lambda x: 100 * x.mean()), "jet0": ("jet0", "mean"), "eady0": ("eady0", "mean")}
    G = W.groupby(["tg", "tn"]).agg(**cols).round(2)
    G.index.names = ["GH tercile (1 low)", "NAO tercile (1 low)"]
    M = W.groupby("tg").agg(**cols).round(2)
    N = W.groupby("tn").agg(**cols).round(2)
    with open(os.path.join(out, "composites.txt"), "w") as f:
        f.write("Atlantic storms (pipeline A, ERA5 proxy), Oct-Apr 2004-05..2025-26. Terciles over storm-level lagged indices (days -10..-4 before genesis).\n")
        f.write(f"tercile cut points: GH z {W.zGH.quantile([1/3,2/3]).round(2).tolist()}, NAO z {W.zNAO.quantile([1/3,2/3]).round(2).tolist()}; r(GH,NAO) over storms = {np.corrcoef(W.zGH, W.zNAO)[0,1]:.2f}\n\n")
        f.write("By GH tercile (NAO not held fixed):\n" + M.to_string() + "\n\nBy NAO tercile (GH not held fixed):\n" + N.to_string() + "\n\n")
        f.write("By GH tercile within NAO tercile (cells under 150 storms are sparse):\n" + G.to_string() + "\n")
    # fix density maps
    F = F.merge(W[["track", "tg", "tn"]], on="track")
    F["lonS"] = sg.sl(F.lon)
    xb, yb = np.arange(-90, 41, 4), np.arange(30, 82, 2)

    def dens(sub, nst):
        h, _, _ = np.histogram2d(sub.lonS, sub.lat, bins=[xb, yb])
        return h.T / max(nst, 1)  # fixes per storm per cell
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6), sharey=True)
    vmax = None
    sets = [("All NAO", W, F), ("Middle NAO tercile", W[W.tn == 2], F[F.tn == 2]), ("Low-NAO tercile", W[W.tn == 1], F[F.tn == 1])]
    for a, (ttl, w, f) in zip(ax, sets):
        hi, lo = dens(f[f.tg == 3], (w.tg == 3).sum()), dens(f[f.tg == 1], (w.tg == 1).sum())
        d = hi - lo
        v = np.nanmax(np.abs(d)) * 0.9 + 1e-9
        pc = a.pcolormesh(xb, yb, d, cmap="RdBu_r", vmin=-v, vmax=v)
        a.set_title(f"{ttl}: strong minus weak GH\n(n storms {int((w.tg==3).sum())} vs {int((w.tg==1).sum())})", fontsize=10)
        a.set_xlabel("lon (deg E)"); a.plot([-50, -50, -15, -15, -50], [55, 67, 67, 55, 55], "k--", lw=0.8)
        plt.colorbar(pc, ax=a, label="fixes per storm per cell")
    ax[0].set_ylabel("lat (deg N)")
    fig.suptitle("Where lows sit (00/12 UTC fixes below 1010 hPa), top minus bottom tercile of lagged Greenland high. Pipeline A, ERA5 proxy. Exploratory.", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "composite_density.png"), dpi=110)
    print(open(os.path.join(out, "composites.txt")).read())


if __name__ == "__main__":
    main()
