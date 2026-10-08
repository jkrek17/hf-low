"""Composite MSLP maps for the cases (ERA5 proxy; pipeline A's lows; outcome G_T).

Uses the 50-80N, 80W-0 MSLP window stored at every candidate time by stage1.py (all 4,559
times, so no sampling). Panels: all cases, barrier cases, tip cases. Shading = case-mean
MSLP minus the mean over non-case times of the same months (weighted to the cases' month
mix); contours = case-mean MSLP. The black dashes mark the 300 km band around Greenland
inside which G_T is measured; the box marks the Cape Farewell tip-jet region.
Writes results/composite.png and results/composite.txt.
"""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(os.environ.get("ERA5_WORK", os.path.join(HERE, "work")), "greenland_jets")
RES = os.path.join(HERE, "results")
LAT = 80 - 0.25 * np.arange(120)                 # rows 40..160 of the 90..-90 grid
LON = 280 + 0.25 * np.arange(320)                # cols 1120..1440


def mean_field(times):
    acc = np.zeros((120, 320))
    for t in times:
        acc += np.load(f"{WORK}/mslp/{t}.npy") / 100.0
    return acc / len(times)


def main():
    d = pd.read_csv(f"{RES}/stage1_times.csv")
    d["month"] = (d.time // 10000) % 100
    nc = d[~d.case]
    out = []
    groups = {"all cases": d[d.case], "barrier cases": d[d.barrier], "tip cases": d[d.tip]}
    # month-matched reference: mean of non-case fields by month, weighted to each group's month mix
    mm = {m: mean_field(nc[nc.month == m].time.astype(str)) for m in sorted(d.month.unique())}
    st = np.load(f"{WORK}/static.npz")
    dgl = st["dgl"][40:160, 1120:1440]
    gl = st["gl"][40:160, 1120:1440]
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6), sharey=True)
    for a, (name, g) in zip(ax, groups.items()):
        f = mean_field(g.time.astype(str))
        mix = g.month.value_counts(normalize=True)
        ref = sum(mix[m] * mm[m] for m in mix.index)
        anom = f - ref
        im = a.pcolormesh(LON - 360, LAT, anom, cmap="RdBu_r", vmin=-14, vmax=14, shading="auto")
        cs = a.contour(LON - 360, LAT, f, levels=np.arange(960, 1040, 4), colors="k", linewidths=0.6)
        a.contour(LON - 360, LAT, gl.astype(float), levels=[0.5], colors="0.3", linewidths=0.9)
        a.contour(LON - 360, LAT, dgl, levels=[300], colors="k", linewidths=0.8, linestyles="--")
        a.add_patch(plt.Rectangle((-50, 59), 15, 2.5, fill=False, ec="m", lw=1.2))
        a.set_title(f"{name} (n={len(g)})", fontsize=10)
        a.set_xlim(-65, -5); a.set_ylim(52, 76)
        a.set_xlabel("longitude")
        out.append(f"{name}: n {len(g)}; mean MSLP anomaly over the Greenland ice sheet box (64-78N, 55-30W) "
                   f"{anom[(LAT[:, None] >= 64) & (LAT[:, None] <= 78) & (LON[None] - 360 >= -55) & (LON[None] - 360 <= -30) & gl].mean():+.1f} hPa; "
                   f"min anomaly {anom.min():+.1f} at {LAT[np.unravel_index(anom.argmin(), anom.shape)[0]]:.1f}N "
                   f"{LON[np.unravel_index(anom.argmin(), anom.shape)[1]] - 360:.1f}E; "
                   f"max anomaly {anom.max():+.1f} at {LAT[np.unravel_index(anom.argmax(), anom.shape)[0]]:.1f}N "
                   f"{LON[np.unravel_index(anom.argmax(), anom.shape)[1]] - 360:.1f}E")
    ax[0].set_ylabel("latitude")
    fig.colorbar(im, ax=ax, shrink=0.85, label="MSLP, cases minus same-month non-case times (hPa)")
    fig.suptitle("Mean sea-level pressure at times with a Greenland-coast gust of 71.7 kt or more (ERA5 proxy, "
                 "2004-05 to 2025-26, Nov-Mar)", fontsize=10)
    fig.savefig(f"{RES}/composite.png", dpi=130, bbox_inches="tight")
    open(f"{RES}/composite.txt", "w").write("\n".join(out) + "\n")
    print("\n".join(out))


if __name__ == "__main__":
    main()
