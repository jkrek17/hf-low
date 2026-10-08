"""Figures for the P(HF)-along-tracks result (ERA5 PROXY, pipeline A). Reads results/ only.

usage: python3 -I figs.py [OUTDIR]       default OUTDIR = results/figs
"""
import os, sys, json, re
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(RES, "figs")
os.makedirs(OUT, exist_ok=True)
COAST = os.path.join(HERE, "..", "..", "..", "docs", "data", "coastlines.js")
EXT = {"atl": (-105, 15, 24, 72), "pac": (125, 250, 24, 72)}
NAME = {"atl": "North Atlantic", "pac": "North Pacific"}
LAND = "#e6e2da"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": "#888", "axes.linewidth": 0.6})


def coast():
    s = open(COAST).read()
    d = json.loads(s[s.index("{"):s.rindex("}") + 1])
    return d["polygons"]


POLY = coast()


def lon_for(b, lon):
    lon = np.asarray(lon, float) % 360
    return np.where(lon > 180, lon - 360, lon) if b == "atl" else lon


def base(ax, b):
    x0, x1, y0, y1 = EXT[b]
    for p in POLY:
        p = np.array(p)
        for sh in ((0,) if b == "atl" else (0, 360)):
            ax.fill(p[:, 0] + sh, p[:, 1], color=LAND, lw=0.3, ec="#b5afa3", zorder=0)
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
    ax.set_aspect(1 / np.cos(np.radians(50)))
    ax.set_facecolor("#f4f8fb")
    ax.set_title(NAME[b], fontsize=10, loc="left")
    ax.grid(color="white", lw=0.5, zorder=0.5)
    ax.tick_params(length=2)


def segs(fx, col):
    """12 h segments between consecutive 00/12 fixes of a track, coloured by P at the segment start."""
    fx = fx.sort_values(["track", "time"])
    t = pd.to_datetime(fx.time, format="%Y%m%d%H")
    same = (fx.track.values[1:] == fx.track.values[:-1]) & ((t.values[1:] - t.values[:-1]) == np.timedelta64(12, "h"))
    i = np.where(same)[0]
    return fx.iloc[i], fx.iloc[i + 1], fx[col].values[i]


def track_map(fx, col, title, fname, note):
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.6), gridspec_kw={"wspace": 0.06})
    cmap = plt.get_cmap("magma_r")
    for ax, b in zip(axs, ("atl", "pac")):
        base(ax, b)
        s = fx[fx.basin == b]
        a, c, p = segs(s, col)
        o = np.argsort(p)
        lo = lon_for(b, a.lon.values); hi = lon_for(b, c.lon.values)
        ok = np.abs(hi - lo) < 60
        L = np.stack([np.stack([lo, a.lat.values], 1), np.stack([hi, c.lat.values], 1)], 1)
        L, p = L[ok], p[ok]
        o = np.argsort(p)
        lc = LineCollection(L[o], array=p[o], cmap=cmap, norm=plt.Normalize(0, 1),
                            linewidths=np.where(p[o] > 0.5, 1.6, np.where(p[o] > 0.1, 0.9, 0.35)),
                            alpha=1.0, zorder=2)
        ax.add_collection(lc)
        ax.text(0.01, 0.02, f"{s.track.nunique():,} tracks, {len(s):,} fixes", transform=ax.transAxes, fontsize=8, color="#333")
    cb = fig.colorbar(lc, ax=axs, orientation="horizontal", fraction=0.04, pad=0.09, aspect=60)
    cb.set_label("P(HF-equivalent gust within 24 h) at the start of each 12 h segment")
    fig.suptitle(title, x=0.125, ha="left", fontsize=11)
    fig.text(0.125, -0.07, note, fontsize=7.5, color="#444", ha="left", va="top", wrap=True)
    fig.savefig(os.path.join(OUT, fname), dpi=140, bbox_inches="tight")
    plt.close(fig)


def main():
    f04 = pd.read_csv(os.path.join(RES, "fix_probs_2004.csv.gz"), dtype={"time": str})
    pre = pd.read_csv(os.path.join(RES, "fix_probs_pre2004.csv.gz"), dtype={"time": str})
    pre = pre[pre.basin.isin(["atl", "pac"])]
    track_map(pre, "P_N",
              "Fig 1  1979-2003 catalog tracks coloured by P(HF within 24 h), environment model N (ERA5 proxy, pipeline A)",
              "fig1_pre2004_tracks_P.png",
              "Events plus matched null cases only (4,311 of ~39,600 tracks), seasons 1979-80 to 2003-04. N = PR 12 full model without the gust index (environment pulled for these fixes). "
              "S-model version: fig1b. Not the archive. Counts before 2001 are not comparable with later seasons.")
    track_map(pre, "P_S",
              "Fig 1b  1979-2003 catalog tracks coloured by P(HF within 24 h), storm-state model S (ERA5 proxy, pipeline A)",
              "fig1b_pre2004_tracks_P_stateonly.png",
              "S = pressure, tendency, age, latitude, speed, day of year; no gust, no environment. Events plus matched null cases only.")
    track_map(f04, "P_N",
              "Fig 2  2004-2025 tracks coloured by P(HF within 24 h), environment model N, leave-one-season-out (ERA5 proxy, pipeline A)",
              "fig2_2004_2025_tracks_P.png",
              "Every low below 1010 hPa at 00/12 UTC (36,598 tracks), seasons 2004-05 to 2025-26. N = PR 12 full model without the current gust index; "
              "each season scored by a model that did not see it.")

    # ---------------- candidates (T2) and calibration (T1/T3)
    T = pd.read_csv(os.path.join(RES, "tracks_2004_summary.csv.gz"))
    R = pd.read_csv(os.path.join(RES, "residual_unlisted_2004.csv"))
    B = pd.read_csv(os.path.join(RES, "residual_unlisted_both_rules_2004.csv"))
    H = T[T.event & (T.Ppre_N >= 0.5)]
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.6), gridspec_kw={"wspace": 0.06})
    for ax, b in zip(axs, ("atl", "pac")):
        base(ax, b)
        for sel, col, sz, z, lab in ((H[H.listed & (H.basin == b)], "#7f8fa0", 9, 1, "archive-listed"),
                                     (R[R.basin == b], "#d9534f", 16, 2, "unlisted residual (strict rule)")):
            ax.scatter(lon_for(b, sel.peak_lon), sel.peak_lat, s=sz, color=col, alpha=0.75, lw=0, zorder=z, label=f"{lab} ({len(sel)})")
        st = B[(B.basin == b) & (B.margin >= 10)]
        ax.scatter(lon_for(b, st.peak_lon), st.peak_lat, s=46, facecolor="none", edgecolor="#111", lw=1.2, zorder=4,
                   label=f"unlisted under both rules, gust >= 81.7 kt ({len(st)})")
        ax.legend(loc="lower left", fontsize=7.5, frameon=True, framealpha=0.9)
    fig.suptitle("Fig 3  Where the gust peaks sit, 2004-2025 proxy events with P_pre >= 0.5 (model N): archive-listed and not (ERA5 proxy, pipeline A)", x=0.125, ha="left", fontsize=11)
    fig.text(0.125, -0.02, "Archive match = PR 52 peak-position rule. Residual = unlisted after removing tropical-cyclone-linked, Atlantic north of 60N and tracks of fewer than 3 fixes. "
             "Tropical-cyclone flags cover only IBTrACS points in tc_candidates.csv (a Gulf hurricane, 2005, is not flagged).", fontsize=7.5, color="#444", va="top")
    fig.savefig(os.path.join(OUT, "fig3_unlisted_high_P.png"), dpi=140, bbox_inches="tight")
    plt.close(fig)

    fig, axs = plt.subplots(1, 3, figsize=(13, 3.8), gridspec_kw={"wspace": 0.3})
    ax = axs[0]
    bins = [0, .1, .3, .5, .7, .9, 1.01]
    labs = ["0-.1", ".1-.3", ".3-.5", ".5-.7", ".7-.9", ".9-1"]
    Z = T[T.Ppre_N.notna()]
    allr, subr, evr = [], [], []
    for lo, hi in zip(bins[:-1], bins[1:]):
        s = Z[(Z.Ppre_N >= lo) & (Z.Ppre_N < hi)]
        allr.append(s.listed.mean()); subr.append(s[~s.event].listed.mean()); evr.append(s[s.event].listed.mean() if s.event.any() else np.nan)
    x = np.arange(len(labs))
    ax.bar(x - 0.27, evr, 0.26, color="#c0392b", label="gust >= 71.7 kt (proxy events)")
    ax.bar(x, subr, 0.26, color="#2e6f9e", label="gust < 71.7 kt")
    ax.bar(x + 0.27, allr, 0.26, color="#888", label="all tracks")
    ax.set_xticks(x); ax.set_xticklabels(labs); ax.set_xlabel("P_pre (model N), largest P before the track is HF")
    ax.set_ylabel("share archive-listed"); ax.legend(fontsize=7); ax.set_title("Archive-listed rate by P_pre, 2004-2025", fontsize=9, loc="left")
    ax = axs[1]
    Lg = H[H.listed & ~H.tc.astype(bool) & ~((H.basin == "atl") & (H.peak_lat >= 60)) & (H.n_dom >= 3)]
    bb = np.arange(0, 42, 2)
    ax.hist(Lg.margin.clip(upper=41), bb, density=True, color="#9aa7b4", alpha=0.8, label=f"listed ({len(Lg)})")
    ax.hist(R.margin.clip(upper=41), bb, density=True, color="#d9534f", alpha=0.6, label=f"unlisted residual ({len(R)})")
    ax.set_xlabel("gust index minus 71.7 kt"); ax.set_ylabel("density"); ax.legend(fontsize=7)
    ax.set_title("Unlisted high-P events sit near the threshold", fontsize=9, loc="left")
    ax = axs[2]
    ax.hist(Lg.minp, np.arange(920, 1000, 4), density=True, color="#9aa7b4", alpha=0.8, label="listed")
    ax.hist(R.minp, np.arange(920, 1000, 4), density=True, color="#d9534f", alpha=0.6, label="unlisted residual")
    ax.set_xlabel("minimum pressure (hPa)"); ax.legend(fontsize=7); ax.set_title("and are shallower", fontsize=9, loc="left")
    fig.suptitle("Fig 4  What P(HF) says about the archive, 2004-05 to 2025-26 (ERA5 proxy, pipeline A)", x=0.125, ha="left", fontsize=11, y=1.03)
    fig.savefig(os.path.join(OUT, "fig4_listed_rate_and_residual.png"), dpi=140, bbox_inches="tight")
    plt.close(fig)

    # ---------------- strongest storms
    S = pd.read_csv(os.path.join(RES, "strongest_storms.csv"))
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.6), gridspec_kw={"wspace": 0.06})
    for ax, b in zip(axs, ("atl", "pac")):
        base(ax, b)
        s = S[(S.basin == b) & (S.rank_by == "minp")]
        lon = lon_for(b, s.lon_minp.fillna(s.peak_lon))
        lat = s.lat_min.fillna(s.peak_lat)
        late = s.season >= 2001
        sc = ax.scatter(lon, lat, c=s.minp, cmap="viridis_r", s=70, edgecolor=np.where(late, "#111", "#c0392b"), lw=1.3, zorder=3, vmin=910, vmax=940)
        for r, x_, y_ in zip(s["rank"], lon, lat):
            if r <= 25:
                ax.text(x_ + 1.2, y_ + 0.8, str(r), fontsize=7, zorder=4)
    cb = fig.colorbar(sc, ax=axs, orientation="horizontal", fraction=0.04, pad=0.09, aspect=60)
    cb.set_label("ERA5 minimum central pressure (hPa)")
    fig.suptitle("Fig 5  The 25 deepest lows per basin, ERA5 1979-80 to 2025-26 (numbers are ranks; red outline = before 2001-02)", x=0.125, ha="left", fontsize=11)
    fig.text(0.125, -0.07, "ERA5 central pressure at the point of minimum pressure, pipeline A tracks. Several Atlantic ones are near Iceland and Greenland. Proxy, not the archive.",
             fontsize=7.5, color="#444", va="top")
    fig.savefig(os.path.join(OUT, "fig5_strongest_storms.png"), dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(os.listdir(OUT))


if __name__ == "__main__":
    main()
