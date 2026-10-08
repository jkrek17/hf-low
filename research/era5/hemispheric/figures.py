"""Figures from results/composite_maps.npz and the held-out JSON. usage: figures.py [RESULTS_DIR]"""
import os, sys, re, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
RES = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "results")
LAT = (-87.1875 + 5.625 * np.arange(32))[20:32]
LON = 5.625 * np.arange(64)
txt = open(os.path.join(REPO, "docs/data/coastlines.js")).read()
POLY = json.loads(txt[txt.index("{"): txt.rindex("}") + 1])["polygons"]


def xy(lat, lon):
    r = (90 - np.asarray(lat))
    th = np.radians(np.asarray(lon) - 90)          # 0E to the right-bottom; Greenwich at the bottom
    return r * np.cos(th), r * np.sin(th)


def base(ax):
    for p in POLY:
        a = np.array(p, float)
        la, lo = a[:, 1].copy(), a[:, 0] % 360
        la[la < 20] = np.nan
        x, y = xy(la, lo)
        # break at the dateline wrap
        wrap = np.where(np.abs(np.diff(lo)) > 180)[0] + 1
        x, y = np.insert(x, wrap, np.nan), np.insert(y, wrap, np.nan)
        ax.plot(x, y, lw=0.5, color="0.3")
    t = np.linspace(0, 2 * np.pi, 200)
    for la in (30, 60):
        ax.plot((90 - la) * np.cos(t), (90 - la) * np.sin(t), lw=0.3, color="0.6", ls=":")
    for lo, nm in ((0, "0"), (90, "90E"), (180, "180"), (270, "90W")):
        x, y = xy(np.array([20.0]), np.array([float(lo)]))
        ax.plot([0, x[0]], [0, y[0]], lw=0.3, color="0.6", ls=":")
        ax.text(x[0] * 1.07, y[0] * 1.07, nm, ha="center", va="center", fontsize=6.5, color="0.2")
    ax.set_xlim(-76, 76); ax.set_ylim(-76, 76); ax.set_aspect("equal"); ax.axis("off")


def panel(ax, field, q, unit, title, lim):
    lon = np.append(LON, 360.0)
    f = np.hstack([field, field[:, :1]]); qq = np.hstack([q, q[:, :1]])
    LO, LA = np.meshgrid(lon, np.append(LAT, 90.0))
    f = np.vstack([f, np.repeat(f[-1:], 1, 0)]); qq = np.vstack([qq, qq[-1:]])
    x, y = xy(LA, LO)
    cs = ax.contourf(x, y, f, levels=np.linspace(-lim, lim, 13), cmap="RdBu_r", extend="both")
    ax.contourf(x, y, np.where(qq < 0.05, 1.0, np.nan), levels=[0.5, 1.5], colors="none", hatches=["..."])
    base(ax)
    ax.set_title(title, fontsize=9)
    return cs


def main():
    d = np.load(os.path.join(RES, "composite_maps.npz"), allow_pickle=True)
    summ = json.load(open(os.path.join(RES, "composite_summary.json")))
    fig, axes = plt.subplots(2, 3, figsize=(10, 6.6))
    spec = [("z500", "m", "Z500 (m)", 60), ("u250", "m/s", "250 hPa zonal wind (m/s)", 8), ("sst", "K", "SST (K)", 0.6)]
    for i, (b, nm) in enumerate((("atl", "Atlantic"), ("pac", "Pacific"))):
        for j, (v, u, t, lim) in enumerate(spec):
            f, q = d[f"{b}_{v}_diff"], d[f"{b}_{v}_q"]
            cs = panel(axes[i, j], f, q, u, f"{nm}: {t}", lim)
            if i == 1:
                cb = fig.colorbar(cs, ax=axes[:, j], orientation="horizontal", fraction=0.04, pad=0.04, aspect=30)
                cb.set_label(f"top minus bottom quintile, {u}", fontsize=8); cb.ax.tick_params(labelsize=7)
    fig.suptitle("7-day mean anomaly in the week before, high minus low quintile of the out-of-sample pattern index\n"
                 "dots: BH-FDR 5% (season-block permutation). ERA5 proxy, 12 UTC, 2004-05 to 2025-26", fontsize=9)
    fig.savefig(os.path.join(RES, "pattern_maps.png"), dpi=130, bbox_inches="tight")
    print({b: {v: round(summ[b][v]["pattern_corr_halves"], 2) for v in summ[b]} for b in summ})


def skill_figure():
    H = {n: json.load(open(os.path.join(RES, f"heldout_{n}.json"))) for n in ("primary", "lag2", "conc", "swap")}
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.4))
    names = [("P|B1", "pattern (field EOFs)"), ("N|B1", "named indices"), ("NI|B1", "named + products"), ("SOM|B1", "SOM regimes")]
    w = 0.2
    for k, (key, lab) in enumerate(names):
        for j, b in enumerate(("atl", "pac")):
            ax[0].bar(j + (k - 1.5) * w, 100 * H["primary"][b]["SS"][key], w, label=lab if j == 0 else None,
                      color=["#1f4e79", "#c0504d", "#e6a09e", "#8a8a8a"][k])
    ax[0].axhline(0, color="k", lw=0.6); ax[0].set_xticks([0, 1]); ax[0].set_xticklabels(["Atlantic", "Pacific"])
    ax[0].set_ylabel("held-out skill vs B1 (%)"); ax[0].legend(fontsize=7, frameon=False)
    ax[0].set_title("2015-16 to 2025-26, models fixed on 2004-05 to 2014-15", fontsize=8)
    for j, (b, nm, c) in enumerate((("atl", "Atlantic", "#1f4e79"), ("pac", "Pacific", "#2e8b57"))):
        y = [100 * H[n][b]["SS"]["P|B1"] for n in ("lag2", "primary", "conc")]
        ax[1].plot([0, 1, 2], y, "o-", color=c, label=nm)
    ax[1].set_xticks([0, 1, 2]); ax[1].set_xticklabels(["days -14 to -8", "days -7 to -1\n(primary)", "same week\n(circular)"], fontsize=8)
    ax[1].axhline(0, color="k", lw=0.6); ax[1].set_ylabel("pattern skill vs B1 (%)"); ax[1].legend(fontsize=8, frameon=False)
    ax[1].set_title("how the signal fades with lead", fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(RES, "skill.png"), dpi=130)


if __name__ == "__main__":
    main(); skill_figure()
