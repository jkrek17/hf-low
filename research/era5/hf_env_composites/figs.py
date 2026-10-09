"""Figures from results/comp_{rot,north}.npz: one file per frame x basin x variable group.

usage: figs.py     -> results/figs/{frame}_{basin}_{group}.png
Rows: raw mean (colour, with contours or wind vectors) and mean anomaly (colour, stippled where the season-block
bootstrap passes BH q < 0.05; stippling is descriptive because pixels are correlated). Columns: 48, 24, 12 h before
HF onset, onset, peak.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import BX, BY
from extract import FIELDS, BASE

HERE = os.path.dirname(os.path.abspath(__file__))
ANCH = ["m48", "m24", "m12", "onset", "peak"]
TITLE = {"m48": "48 h before onset", "m24": "24 h before onset", "m12": "12 h before onset", "onset": "HF onset", "peak": "HF peak"}
BASIN = {"atl": "North Atlantic", "pac": "North Pacific"}
GROUPS = {
    "z500": [("z500", "Z500 (m)", "raw", "contour"), ("z500", "Z500 anomaly (m)", "anom", None)],
    "wind250": [("ws250", "250 hPa wind speed (m/s)", "raw", "vec250"), ("ws250", "250 hPa speed anomaly (m/s)", "anom", None)],
    "wind500": [("ws500", "500 hPa wind speed (m/s)", "raw", "vec500"), ("ws500", "500 hPa speed anomaly (m/s)", "anom", None)],
    "div": [("d250", "250 hPa divergence (1e-5/s)", "raw", None), ("d250", "250 hPa divergence anomaly (1e-5/s)", "anom", None),
            ("d500", "500 hPa divergence (1e-5/s)", "raw", None), ("d500", "500 hPa divergence anomaly (1e-5/s)", "anom", None)],
}
CMAP = {"raw": "viridis", "anom": "RdBu_r"}


def panel(ax, Z, kind, vmin, vmax, extra=None, rej=None, xlab=None):
    im = ax.pcolormesh(BX, BY, Z, cmap=CMAP[kind], vmin=vmin, vmax=vmax, shading="auto", rasterized=True)
    if extra is not None:
        extra(ax)
    if rej is not None:
        yy, xx = np.where(rej[::3, ::3])
        ax.scatter(BX[::3, ::3][yy, xx], BY[::3, ::3][yy, xx], s=1.0, c="k", alpha=.55, linewidths=0)
    ax.plot(0, 0, "k+", ms=8, mew=1.5)
    ax.set_aspect("equal"); ax.set_xticks([-3000, 0, 3000]); ax.set_yticks([-3000, 0, 3000])
    ax.tick_params(labelsize=7)
    return im


def main():
    os.makedirs(os.path.join(HERE, "results", "figs"), exist_ok=True)
    for fr in ("rot", "north"):
        z = np.load(os.path.join(HERE, "results", f"comp_{fr}.npz"))
        xl = ("ahead of motion (km)", "left of motion (km)") if fr == "rot" else ("east (km)", "north (km)")
        for basin in ("atl", "pac"):
            for gname, rowsdef in GROUPS.items():
                nr = len(rowsdef)
                fig, axs = plt.subplots(nr, 5, figsize=(14, 2.55 * nr + 0.9), squeeze=False)
                ims = {}
                for r, (base, lab, kind, ex) in enumerate(rowsdef):
                    fi = FIELDS.index(base if kind == "raw" else "a_" + base)
                    allv = np.stack([z[f"{basin}_{a}_mean"][fi] for a in ANCH])
                    if kind == "raw":
                        vmin, vmax = (np.nanpercentile(allv, 1), np.nanpercentile(allv, 99))
                        if base.startswith("d"):
                            m = max(abs(vmin), abs(vmax)); vmin, vmax = -m, m
                    else:
                        m = np.nanpercentile(np.abs(allv), 99); vmin, vmax = -m, m
                    for c, a in enumerate(ANCH):
                        mean = z[f"{basin}_{a}_mean"]
                        extra = None
                        if ex == "contour":
                            extra = lambda ax, M=mean[fi]: ax.contour(BX, BY, M, levels=np.arange(4800, 6000, 100), colors="w", linewidths=.6)
                        if ex and ex.startswith("vec"):
                            lv = ex[3:]
                            U, V = mean[FIELDS.index("u" + lv)], mean[FIELDS.index("v" + lv)]
                            extra = lambda ax, U=U, V=V: ax.quiver(BX[::6, ::6], BY[::6, ::6], U[::6, ::6], V[::6, ::6], color="w", scale=420, width=.004)
                        rej = z[f"{basin}_{a}_rej"][BASE.index(base) if kind == "anom" else 0] if kind == "anom" else None
                        if kind == "anom":
                            rej = z[f"{basin}_{a}_rej"][BASE.index(base)]
                        im = panel(axs[r, c], mean[fi], kind, vmin, vmax, extra, rej)
                        if r == 0:
                            axs[r, c].set_title(f"{TITLE[a]}\nn = {int(z[f'{basin}_{a}_n'])}", fontsize=9)
                        if c == 0:
                            axs[r, c].set_ylabel(lab + "\n" + xl[1], fontsize=7)
                        if r == nr - 1:
                            axs[r, c].set_xlabel(xl[0], fontsize=7)
                    fig.colorbar(im, ax=axs[r, :].tolist(), shrink=.85, pad=.01, aspect=25).ax.tick_params(labelsize=7)
                fig.suptitle(f"{BASIN[basin]}, storm-centred mean of HF lows (ERA5 proxy, Oct-Apr onsets 2004-05 to 2022-23), "
                             f"{'rotated to the direction of motion' if fr == 'rot' else 'north-up'}", fontsize=10)
                fig.savefig(os.path.join(HERE, "results", "figs", f"{fr}_{basin}_{gname}.png"), dpi=110, bbox_inches="tight")
                plt.close(fig)
    print(sorted(os.listdir(os.path.join(HERE, "results", "figs"))))


if __name__ == "__main__":
    main()
