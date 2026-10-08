"""Figures for the climatology atlas. Reads results/ and recomputes the box fields from the committed inputs.

usage: python3 -I research/era5/climatology_atlas/figures.py [outdir]     (default results/)

ARCHIVE = OPC archive, HF-window metrics. PROXY = ERA5 pipeline A, a proxy. Seasons 2004-05 to 2025-26.
"""
import json
import os
import re
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import atlas as A_  # noqa: E402

RES = os.path.join(HERE, "results")
OUT = sys.argv[1] if len(sys.argv) > 1 else RES
os.makedirs(OUT, exist_ok=True)
COL = {"archive": "#1f4e79", "proxy": "#d17a22"}
BN = A_.BN
EXT = {"atl": (-100, 30, 20, 75), "pac": (120, 250, 20, 75)}
MNAME = {10: "Oct", 11: "Nov", 12: "Dec", 1: "Jan", 2: "Feb", 3: "Mar", 4: "Apr"}


def coast():
    s = open(os.path.join(A_.ROOT, "docs", "data", "coastlines.js")).read()
    return json.loads(re.search(r"HF_COAST = (\{.*\});?\s*$", s, re.S).group(1))["polygons"]


POLYS = coast()


def draw_land(ax, basin):
    for poly in POLYS:
        p = np.array(poly)
        lon = np.degrees(np.unwrap(np.radians(p[:, 0])))     # continuous across the dateline
        for off in (-360, 0, 360):
            ax.fill(lon + off, p[:, 1], color="#d9d9d9", lw=0.3, ec="#888888", zorder=3)
    x0, x1, y0, y1 = EXT[basin]
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect(1 / np.cos(np.radians(47)))


def pcm(ax, basin, H, lat_e, lon_e, vmax):
    m = ax.pcolormesh(lon_e, lat_e, np.ma.masked_less_equal(H, 0), cmap="YlOrRd", vmin=0, vmax=vmax, zorder=2)
    draw_land(ax, basin)
    return m


def fig_centre_hours(AF, PF):
    fig, axs = plt.subplots(2, 2, figsize=(11, 8.2), constrained_layout=True)
    for j, b in enumerate(A_.BASINS):
        vmax = max(A_.box_hours(AF, b)[2].max(), A_.box_hours(PF, b)[2].max())
        for i, (src, fx) in enumerate((("archive", AF), ("proxy", PF))):
            la, lo, H = A_.box_hours(fx, b)
            m = pcm(axs[i, j], b, H, la, lo, vmax)
            axs[i, j].set_title(f"{BN[b]}: {src.upper()}  (total {H.sum():.0f} h per season)", fontsize=10)
        fig.colorbar(m, ax=axs[:, j], shrink=0.7, label="HF-centre hours per season per 5° x 10° box")
    fig.suptitle("Where hurricane-force lows sit: HF-centre hours per season, 2004-05 to 2025-26. Proxy = ERA5 pipeline A", fontsize=11)
    fig.savefig(os.path.join(OUT, "fig1_hf_centre_hours.png"), dpi=120)
    plt.close(fig)


def fig_by_month(AF):
    fig, axs = plt.subplots(2, 7, figsize=(17, 6.2), constrained_layout=True)
    for i, b in enumerate(A_.BASINS):
        f = AF[AF.basin == b]
        vmax = 0
        Hs = {}
        for m in A_.MONTHS:
            la, lo, H = A_.box_hours(f[f.t.dt.month == m], b)
            Hs[m] = (la, lo, H)
            vmax = max(vmax, H.max())
        for j, m in enumerate(A_.MONTHS):
            la, lo, H = Hs[m]
            mm = pcm(axs[i, j], b, H, la, lo, vmax)
            axs[i, j].set_title(f"{BN[b]} {MNAME[m]}", fontsize=9)
            axs[i, j].set_xticks([])
            axs[i, j].set_yticks([])
        fig.colorbar(mm, ax=axs[i, :], shrink=0.8, pad=0.01, label="h / season / box")
    fig.suptitle("ARCHIVE: HF-centre hours per season by calendar month (by month of the fix)", fontsize=11)
    fig.savefig(os.path.join(OUT, "fig2_hf_centre_hours_by_month.png"), dpi=110)
    plt.close(fig)


def fig_points(A, P):
    fig, axs = plt.subplots(2, 3, figsize=(15, 8), constrained_layout=True)
    for i, b in enumerate(A_.BASINS):
        for j, (what, lo_, la_, name) in enumerate((("first HF fix", "lon_on", "lat_on", "Where storms first reach HF"),
                                                     ("minimum pressure", "lon_mp", "lat_mp", "Where they reach minimum pressure"),
                                                     ("last HF fix", "lon_off", "lat_off", "Where they leave HF"))):
            ax = axs[i, j]
            for src, ev in (("archive", A), ("proxy", P)):
                d = ev[ev.basin == b].dropna(subset=[lo_, la_])
                ax.scatter(d[lo_], d[la_], s=5, alpha=0.25, color=COL[src], zorder=2, label=f"{src} (n={len(d)})", lw=0)
                ax.plot(d[lo_].median(), d[la_].median(), marker="*", ms=14, color=COL[src], mec="k", zorder=4)
            draw_land(ax, b)
            ax.set_title(f"{BN[b]}: {name}", fontsize=10)
            ax.legend(loc="lower left", fontsize=8, markerscale=3)
    fig.suptitle("Position of each HF low at three points of its life (star = median). Archive first fix is where OPC began warning, not genesis", fontsize=11)
    fig.savefig(os.path.join(OUT, "fig3_positions.png"), dpi=110)
    plt.close(fig)


def fig_monthly():
    mt = pd.read_csv(os.path.join(RES, "monthly_exposure.csv"))
    fig, axs = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)
    for j, b in enumerate(A_.BASINS):
        for i, (c, lo, hi, lab) in enumerate((("events_per_season", "ev_lo", "ev_hi", "HF events per season"),
                                               ("hf_hours_per_season", "h_lo", "h_hi", "HF hours per season"))):
            ax = axs[i, j]
            for k, src in enumerate(("archive", "proxy")):
                s = mt[(mt.basin == b) & (mt.source == src)].set_index("month").loc[A_.MONTHS]
                x = np.arange(len(A_.MONTHS)) + (k - 0.5) * 0.38
                ax.bar(x, s[c], 0.38, color=COL[src], label=src, yerr=[s[c] - s[lo], s[hi] - s[c]], capsize=2, error_kw=dict(lw=0.8))
            ax.set_xticks(range(len(A_.MONTHS)))
            ax.set_xticklabels([MNAME[m] for m in A_.MONTHS])
            ax.set_ylabel(lab)
            ax.set_title(BN[b], fontsize=10)
            ax.legend(fontsize=8)
    fig.suptitle("By month of first HF fix; bars = mean over 22 seasons, whiskers = 95% season-block bootstrap", fontsize=11)
    fig.savefig(os.path.join(OUT, "fig4_monthly_exposure.png"), dpi=120)
    plt.close(fig)


def fig_seasons():
    st = pd.read_csv(os.path.join(RES, "season_counts.csv"))
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.2), constrained_layout=True, sharey=False)
    for j, b in enumerate(A_.BASINS):
        for src in ("archive", "proxy"):
            s = st[(st.basin == b) & (st.source == src)]
            axs[j].plot(s.season, s.events, "-o", ms=4, color=COL[src], label=src)
        axs[j].set_title(BN[b], fontsize=10)
        axs[j].set_xlabel("season (start year; 1 June to 31 May)")
        axs[j].set_ylabel("HF events")
        axs[j].legend(fontsize=8)
        axs[j].set_ylim(bottom=0)
    fig.suptitle("Events per season. No trend is estimated or implied (decision 1); this is the spread to expect", fontsize=11)
    fig.savefig(os.path.join(OUT, "fig5_seasons.png"), dpi=120)
    plt.close(fig)


def fig_intensity(A, P):
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.2), constrained_layout=True)
    bins = np.arange(920, 1011, 4)
    for src, ev in (("archive", A), ("proxy", P)):
        for b, ls in (("atl", "-"), ("pac", "--")):
            v = ev[ev.basin == b].minp.dropna()
            axs[0].hist(v, bins=bins, histtype="step", density=True, color=COL[src], ls=ls, lw=1.5, label=f"{src} {BN[b]}")
    axs[0].set_xlabel("minimum pressure (hPa; archive = analysed, proxy = ERA5)")
    axs[0].set_ylabel("density")
    axs[0].legend(fontsize=7)
    for b, ls in (("atl", "-"), ("pac", "--")):
        v = P[P.basin == b].gust
        axs[1].hist(v, bins=np.arange(70, 130, 2.5), histtype="step", density=True, color=COL["proxy"], ls=ls, lw=1.5, label=f"proxy {BN[b]}")
    axs[1].set_xlabel("pipeline A 800 km gust index (kt); every event is at least 71.7")
    axs[1].legend(fontsize=8)
    for b, ls in (("atl", "-"), ("pac", "--")):
        v = P[P.basin == b].maxdeep.dropna()
        axs[2].hist(v, bins=np.arange(0, 4.1, 0.2), histtype="step", density=True, color=COL["proxy"], ls=ls, lw=1.5, label=f"proxy {BN[b]}")
    axs[2].axvline(1, color="k", lw=0.8)
    axs[2].set_xlabel("maximum 24 h deepening (Bergerons); 1 = bomb")
    axs[2].legend(fontsize=8)
    fig.suptitle("How strong: depth and wind at peak (proxy deepening needs a track surviving 24 h)", fontsize=11)
    fig.savefig(os.path.join(OUT, "fig6_intensity.png"), dpi=120)
    plt.close(fig)


def fig_duration(A, P):
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.2), constrained_layout=True)
    for ax, (c, lab, bins) in zip(axs, (("hf_h", "hours at HF (6 h per HF fix)", np.arange(0, 85, 6)),
                                         ("hf_spd", "translation speed while HF (kt, mean per event)", np.arange(0, 80, 3)),
                                         ("hf_dist_km", "distance covered while HF (km, archive)", np.arange(0, 3500, 150)))):
        for src, ev in (("archive", A), ("proxy", P)):
            if c not in ev:
                continue
            for b, ls in (("atl", "-"), ("pac", "--")):
                v = ev[ev.basin == b][c].dropna()
                if len(v):
                    ax.hist(v, bins=bins, histtype="step", density=True, color=COL[src], ls=ls, lw=1.5, label=f"{src} {BN[b]}")
        ax.set_xlabel(lab)
        ax.legend(fontsize=7)
    fig.suptitle("How long and how fast", fontsize=11)
    fig.savefig(os.path.join(OUT, "fig7_duration_speed.png"), dpi=120)
    plt.close(fig)


def fig_motion():
    mo = pd.read_csv(os.path.join(RES, "motion_steps.csv.gz"))
    fig, axs = plt.subplots(2, 2, figsize=(12, 8.2), constrained_layout=True)
    for j, b in enumerate(A_.BASINS):
        for i, src in enumerate(("archive", "proxy")):
            ax = axs[i, j]
            m = mo[(mo.source == src) & (mo.basin == b)].copy()
            m["bi"] = (m.latc // 5) * 5
            m["bj"] = (m.lonc // 10) * 10
            g = m.groupby(["bi", "bj"]).agg(u=("u", "mean"), v=("v", "mean"), n=("u", "size")).reset_index()
            g = g[g.n >= 25]
            q = ax.quiver(g.bj + 5, g.bi + 2.5, g.u, g.v, np.hypot(g.u, g.v), cmap="viridis", angles="xy", scale_units="xy",
                          scale=6, width=0.004, zorder=4, clim=(10, 45))
            draw_land(ax, b)
            ax.set_title(f"{BN[b]}: {src.upper()}  ({len(m)} six-hour steps)", fontsize=10)
        fig.colorbar(q, ax=axs[:, j], shrink=0.7, label="mean speed (kt)")
    fig.suptitle("Mean motion of HF lows while HF, in 5° x 10° boxes with at least 25 steps (arrow length = speed; direction is true)", fontsize=11)
    fig.savefig(os.path.join(OUT, "fig8_motion.png"), dpi=110)
    plt.close(fig)


def fig_simult_tc():
    sim = pd.read_csv(os.path.join(RES, "simultaneous.csv"))
    tc = pd.read_csv(os.path.join(RES, "tc_share_by_month.csv"))
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.2), constrained_layout=True)
    ax = axs[0]
    w = 0.2
    for i, (b, src) in enumerate([(b, s) for b in A_.BASINS for s in ("archive", "proxy")]):
        r = sim[(sim.basin == b) & (sim.source == src)].iloc[0]
        ax.bar(np.arange(4) + (i - 1.5) * w, [r.share_ge1, r.share_ge2, r.share_ge3, r.share_ge4], w,
               color=COL[src], alpha=1.0 if b == "atl" else 0.55, label=f"{BN[b]} {src}")
    ax.set_xticks(range(4))
    ax.set_xticklabels(["≥1", "≥2", "≥3", "≥4"])
    ax.set_xlabel("HF lows active at the same time in the basin")
    ax.set_ylabel("share of 6-hourly times, Oct-Apr")
    ax.legend(fontsize=8)
    ax = axs[1]
    for b, c in (("atl", "#1f4e79"), ("pac", "#7b3294")):
        s = tc[tc.basin == b]
        ax.plot(s.month.map(lambda m: (m - 8) % 12), s.tc / s.n, "-o", color=c, label=BN[b])
    ax.set_xticks(range(10))
    ax.set_xticklabels(["Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May"])
    ax.set_ylabel("share of proxy events within 400 km of a tropical cyclone")
    ax.legend(fontsize=8)
    fig.suptitle("Left: how many HF lows at once. Right: proxy events near a tropical cyclone, by month of first HF fix (first look)", fontsize=11)
    fig.savefig(os.path.join(OUT, "fig9_simultaneous_tc.png"), dpi=120)
    plt.close(fig)


def fig_latitude(A, P):
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.2), constrained_layout=True)
    for j, b in enumerate(A_.BASINS):
        for k, (src, ev) in enumerate((("archive", A), ("proxy", P))):
            xs, med, q1, q3 = [], [], [], []
            for i, m in enumerate(A_.MONTHS):
                v = ev[(ev.basin == b) & (ev.mon_on == m)].lat_on
                xs.append(i + (k - 0.5) * 0.15)
                med.append(v.median())
                q1.append(v.quantile(.25))
                q3.append(v.quantile(.75))
            axs[j].errorbar(xs, med, yerr=[np.array(med) - q1, np.array(q3) - med], fmt="o-", color=COL[src], capsize=3, label=src)
        axs[j].set_xticks(range(len(A_.MONTHS)))
        axs[j].set_xticklabels([MNAME[m] for m in A_.MONTHS])
        axs[j].set_ylabel("latitude of first HF fix (°N): median and quartiles")
        axs[j].set_title(BN[b], fontsize=10)
        axs[j].legend(fontsize=8)
    fig.suptitle("Seasonal migration of the place where storms first reach hurricane force", fontsize=11)
    fig.savefig(os.path.join(OUT, "fig10_latitude_by_month.png"), dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    A, AF = A_.load_archive()
    P, PF, _, _ = A_.load_proxy()
    fig_centre_hours(AF, PF)
    fig_by_month(AF)
    fig_points(A, P)
    fig_monthly()
    fig_seasons()
    fig_intensity(A, P)
    fig_duration(A, P)
    fig_motion()
    fig_simult_tc()
    fig_latitude(A, P)
    print("figures written to", OUT)
