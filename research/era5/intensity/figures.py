"""Hart-style phase diagrams shaded by what happened next.

Two panels per outcome, B against -V_T lower and -V_T lower against -V_T
upper (Hart 2003's two views), each cell shaded by the observed frequency of
the outcome over the next 24 h among ERA5 fixes in that cell, seasons 2004-05 to 2025-26 (model.load applies the
project's fit-and-test floor).

usage: figures.py FIXES_CSV ENV_CSV OUTDIR
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import model

MIN_N = 30


def panel(ax, x, y, v, xe, ye, vmax, label):
    n, _, _ = np.histogram2d(x, y, [xe, ye])
    s, _, _ = np.histogram2d(x, y, [xe, ye], weights=v)
    with np.errstate(invalid="ignore"):
        f = np.where(n >= MIN_N, s / n, np.nan)
    im = ax.pcolormesh(xe, ye, f.T, cmap="viridis", vmin=0, vmax=vmax, shading="flat")
    ax.axhline(0, color="0.5", lw=0.8)
    ax.axvline(0, color="0.5", lw=0.8)
    return im


def main(fx, ev, outdir):
    os.makedirs(outdir, exist_ok=True)
    D = model.load(fx, ev)
    C = D[D.cls >= 0]
    H = D[~D.hf_now.astype(bool)]
    outcomes = [("rapid deepening (>= 1 Bergeron in 24 h)", C, (C.cls == 4).astype(float), 0.6, "phase_rapid_deepening.png"),
                ("rapid decay (<= -1 Bergeron in 24 h)", C, (C.cls == 0).astype(float), 0.15, "phase_rapid_decay.png"),
                ("HF-equivalent gust within 24 h, not HF now", H, H.hf24.astype(float), 0.5,
                 "phase_hf24_onset.png")]
    be, le, ue = np.arange(-70, 151, 10), np.arange(-560, 301, 40), np.arange(-560, 301, 40)
    for title, S, v, vmax, name in outcomes:
        fig, axs = plt.subplots(1, 2, figsize=(11, 4.6), constrained_layout=True)
        panel(axs[0], S.VTL.values, S.B.values, v.values, le, be, vmax, title)
        axs[0].set_xlabel("-V_T lower, 900-600 hPa (m)")
        axs[0].set_ylabel("B, 900-600 hPa thickness asymmetry (m)")
        axs[0].axhline(10, color="0.5", lw=0.8, ls="--")
        im = panel(axs[1], S.VTL.values, S.VTU.values, v.values, le, ue, vmax, title)
        axs[1].set_xlabel("-V_T lower, 900-600 hPa (m)")
        axs[1].set_ylabel("-V_T upper, 600-300 hPa (m)")
        fig.colorbar(im, ax=axs, label="observed frequency")
        fig.suptitle(f"ERA5 proxy, lows at 00/12 UTC, 2004-26: P({title}); cells with >= {MIN_N} fixes")
        fig.savefig(os.path.join(outdir, name), dpi=110)
        plt.close(fig)


if __name__ == "__main__":
    main(*sys.argv[1:4])
