#!/usr/bin/env python3
"""Storm review template v2: one standalone page per storm, with the meteorology behind the numbers.

usage: storm_review_v2.py WORKDIR FIXES_CSV OUTDIR TRACK [TRACK ...]
  WORKDIR     output of storm_review_extract.py (one npz per track) and, written here, env_<track>.csv
  FIXES_CSV   research/era5/hf_history/results/era5_hf_catalog_tracks.csv (pipeline A, every 6 h point)
  OUTDIR      where storm-review-<track>-*.png and storm-review-<track>.html are written
Everything is an ERA5 proxy (pipeline A tracks). The text of each page is generated from the numbers in this script's own
tables, not typed by hand. Template: tools/storm_review_template.md.
"""
import os, sys, html
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import build_storm_reviews as old
import hart_cps as hart
import storm_review_env as senv
import storm_review_extract as ex

HF = 71.7
KM_NM = 1.852
MAP_H = old.MAP_H
TERRAIN_M = 1000.0                              # m: terrain above this is masked on the low-level maps
PHASES = {"fcc": ("frontal, cold core", "#3b82f6"), "fwc": ("frontal, warm core", "#f2a93b"),
          "wsec": ("warm seclusion", "#e5502f"), "scc": ("symmetric, cold core", "#9ca3af")}
B_EDGE = 10.0                                   # m, Hart's asymmetry threshold
QUAD = (("front right", 1, 1), ("rear right", -1, 1), ("rear left", -1, -1), ("front left", 1, -1))
VARS = [("jet250", "strongest 250 hPa wind (the jet) within 540 nm", "kt", 0), ("div300", "upper-level divergence (300 hPa) within 270 nm", "&times;10<sup>&minus;5</sup> /s", 1),
        ("vadv500", "vorticity advection (500 hPa) within 270 nm", "&times;10<sup>&minus;9</sup> /s&sup2;", 1), ("eady", "Eady growth rate, the support for storm growth (850&ndash;500 hPa)", "/day", 2),
        ("tcwv", "moisture in the air column within 270 nm", "kg/m&sup2;", 0), ("sstgrad", "sharpest sea-surface temperature gradient within 270 nm", "K/100 km", 1)]


def phase_of(B, VTL):
    if not (np.isfinite(B) and np.isfinite(VTL)):
        return None
    return ("fcc" if VTL <= 0 else "fwc") if B >= B_EDGE else ("scc" if VTL <= 0 else "wsec")


# ------------------------------------------------------------------ data

def track_env(work, fx, track):
    """fx: the storm's prepped fixes (build_storm_reviews.prep). Adds the 1.5 degree diagnostics and the Hart class per fix."""
    path = os.path.join(work, f"env_{track}.csv")
    if not os.path.exists(path):
        h, _ = hart.track_motion(fx.lat.values, fx.lon.values, fx.t.astype("int64").values / 1e9)
        fixes = pd.DataFrame(dict(track=track, time=fx.t.dt.strftime("%Y%m%d%H"), lat=fx.lat, lon=fx.lon, heading=h))
        senv.env_for_fixes(fixes).to_csv(path, index=False)
    e = pd.read_csv(path, dtype={"time": str})
    e["t"] = pd.to_datetime(e.time, format="%Y%m%d%H")
    e["phase"] = [phase_of(b, v) for b, v in zip(e.B, e.VTL)]
    return fx.merge(e.drop(columns=["track", "time"]), on="t", how="left")


def opc_match(fx, basin):
    """Match the storm's pipeline A fixes to OPC archive fixes by time and position (within 200 nm). Returns None when nothing matches."""
    d = pd.read_csv(os.path.join(HERE, "..", "data", "hf_lows", "HF_Data_-_Atl.csv" if basin == "atl" else "HF_Data_-_Pac.csv"))
    d["t"] = pd.to_datetime(d.date.astype(str), format="%Y%m%d%H", errors="coerce")      # one Pacific row has an 11-digit date
    d["Pressure"] = pd.to_numeric(d.Pressure, errors="coerce")
    d = d.dropna(subset=["t", "Pressure"])
    d["Longitude"] = d.Longitude % 360
    m = d.merge(fx[["t", "lat", "lon", "msl", "gust"]], on="t")
    dl = np.radians(((m.Longitude - m.lon + 180) % 360) - 180); p1, p2 = np.radians(m.Latitude), np.radians(m.lat)
    m["dist"] = 2 * np.arcsin(np.sqrt(np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2)) * 3440.065
    m = m[m.dist <= 200]
    if m.empty:
        return None
    best = m.groupby("ID").size().idxmax()
    mm = m[m.ID == best]
    o = d[(d.ID == best) & d.t.isin(mm.t)].sort_values("t")        # only the OPC fixes that sit on this low (within 200 nm of the track)
    hf = o[o.Category == "HF"]
    era_hf = fx[fx.gust >= HF]
    both = set(hf.t) & set(era_hf.t)
    i = o.Pressure.idxmin()
    return dict(id=str(best), n_fix=len(o), first=o.t.iloc[0], last=o.t.iloc[-1], n_hf=len(hf), hf_first=hf.t.iloc[0] if len(hf) else None,
                hf_last=hf.t.iloc[-1] if len(hf) else None, minp=float(o.Pressure[i]), t_minp=o.t[i], n_matched=len(mm), dist_med=float(mm.dist.median()),
                n_era_hf=len(era_hf), n_both=len(both), era_hf_in_opc=len(set(era_hf.t) & set(o.t)))


_REF = None


def ref_rank(basin, var, value):
    global _REF
    if _REF is None:
        _REF = pd.read_csv(os.path.join(HERE, "storm_review_ref.csv"))
    r = _REF[(_REF.basin == basin) & (_REF["var"] == var)]
    if r.empty or not np.isfinite(value):
        return None
    q = r.iloc[0][[f"q{i}" for i in range(101)]].astype(float).values
    return float(np.interp(value, q, np.arange(101)))


def sr_data(z):
    """Storm-relative gust cells and per-time stats. Where the storm is nearly stationary (< 10 kt) the heading is unreliable, so
    the last heading from a time it was moving 10 kt or more is carried forward and the cells are rotated into that frame."""
    st = pd.DataFrame(z["sr_stats"], columns=[str(c) for c in z["sr_stats_cols"]])
    hd = st.heading.copy(); last = np.nan
    for i in range(len(st)):
        if st.speed_kt[i] >= 10:
            last = st.heading[i]
        hd[i] = last
    if hd.isna().any():
        hd = hd.bfill()
    st["frame_hd"] = hd
    c = pd.DataFrame(z["sr_cells"].T, columns=["h", "x", "y", "g"])
    rot = np.radians(st.set_index("h").heading.reindex(c.h).values - st.set_index("h").frame_hd.reindex(c.h).values)   # frame change
    x, y = c.x.values, c.y.values
    c["x"], c["y"] = x * np.cos(rot) + y * np.sin(rot), -x * np.sin(rot) + y * np.cos(rot)
    return st, c


def orography(work):
    """Surface height (m) on the 0.25 degree grid, fetched once (static field, about 1 MB)."""
    path = os.path.join(work, "orography_m.npy")
    if not os.path.exists(path):
        np.save(path, (ex.get(f"{ex.WB2}/geopotential_at_surface/0.0", (721, 1440)) / 9.80665).astype(np.float32))
    return np.load(path)


# ------------------------------------------------------------------ figures

def ptitle(ax, text, sub=""):
    ax.text(0, 1.04, text, transform=ax.transAxes, fontsize=9.5, fontweight="semibold", color=old.INK, va="bottom")
    if sub:
        ax.text(1, 1.04, sub, transform=ax.transAxes, fontsize=8, color=old.MUTED, va="bottom", ha="right")


def fig_series2(fx, d24, t0, st, title, path):
    hrs = (fx.t - t0) / pd.Timedelta(hours=1)
    j = int(np.nanargmin(d24)); end = hrs[j]
    fig, axs = plt.subplots(8, 1, figsize=(9.2, 16.2), sharex=True, gridspec_kw=dict(hspace=0.62))
    for a in axs:
        a.axvspan(end - 24, end, color="#fbe7c2", alpha=0.6, lw=0, zorder=0)
        a.axvline(0, color="#9ca3af", lw=0.9, ls=(0, (2, 3)), zorder=1); a.grid(axis="y", color=old.GRID, lw=0.8)
    BLUE, ORANGE, VIOLET, TEAL, RED, BROWN = "#2563eb", "#e5502f", "#6d5bd0", "#0f9d8a", "#d1495b", "#b45309"
    LW = dict(lw=2.2, solid_capstyle="round")
    a = axs[0]; ptitle(a, "Central pressure (hPa)")
    a.plot(hrs, fx.msl, color=BLUE, **LW); a.set_ylim(fx.msl.max() + 10, fx.msl.min() - 12)
    a.fill_between(hrs, fx.msl, fx.msl.max() + 10, color=BLUE, alpha=0.08, lw=0)
    i0 = int(fx.msl.values.argmin()); a.scatter([hrs[i0]], [fx.msl[i0]], s=46, color=BLUE, edgecolor="white", linewidth=1.5, zorder=5)
    a.annotate(f"{fx.msl[i0]:.0f} hPa", (hrs[i0], fx.msl[i0]), xytext=(12, -16), textcoords="offset points", fontsize=9, color=BLUE, fontweight="semibold", va="center")
    a = axs[1]; ptitle(a, "Fall over the previous 24 h (hPa)", "dashed: bomb yardstick at the storm's latitude")
    dd = -d24
    a.bar(hrs, dd, width=3.6, color=np.where(dd > 0, ORANGE, "#c7ccd4"), zorder=3)
    a.plot(hrs, 24 * np.sin(np.radians(fx.lat)) / np.sin(np.radians(60)), color=old.INK, lw=1.1, ls=(0, (4, 3)), zorder=4)
    a = axs[2]; ptitle(a, "Gust index (kt)", "highest ERA5 gust within 800 km; dashed: 71.7 kt")
    a.plot(hrs, fx.gust, color=ORANGE, **LW); a.axhline(HF, color=old.INK, lw=1, ls=(0, (4, 3)))
    ig = int(np.nanargmax(fx.gust.values)); a.scatter([hrs[ig]], [fx.gust[ig]], s=40, color=ORANGE, edgecolor="white", linewidth=1.4, zorder=5)
    a.annotate(f"{fx.gust[ig]:.0f} kt", (hrs[ig], fx.gust[ig]), xytext=(8, 4), textcoords="offset points", fontsize=9, color=ORANGE, fontweight="semibold")
    a = axs[3]; ptitle(a, "Ocean area with 71.7 kt+ gusts (thousand nm\u00b2)", "within 600 nm of the centre, 3-hourly")
    a.bar(st.h, st.area72_nm2 / 1000, width=2.4, color="#f08a63", zorder=3)
    a = axs[4]; ptitle(a, "250 hPa jet (kt)", "strongest wind within 1,000 km")
    a.plot(hrs, fx.jet250, color=VIOLET, **LW); a.fill_between(hrs, fx.jet250, fx.jet250.min() - 15, color=VIOLET, alpha=0.08, lw=0); a.set_ylim(fx.jet250.min() - 15, None)
    a = axs[5]; ptitle(a, "Upper-level lift", "")
    a.plot(hrs, fx.div300, color=TEAL, **LW); a.plot(hrs, fx.vadv500, color=RED, **LW); a.axhline(0, color="#9ca3af", lw=0.8)
    a.text(1, 1.04, "500 hPa vorticity advection (1e-9 /s\u00b2)", transform=a.transAxes, color=RED, fontsize=8, va="bottom", ha="right", fontweight="semibold")
    a.text(0.6, 1.04, "300 hPa divergence (1e-5 /s)", transform=a.transAxes, color=TEAL, fontsize=8, va="bottom", ha="right", fontweight="semibold")
    a = axs[6]; ptitle(a, "Hart B: thickness asymmetry (m)", "dashed: 10 m; above it the low is frontal")
    a.plot(hrs, fx.B, color=BROWN, **LW); a.axhline(B_EDGE, color=BROWN, lw=0.9, ls=(0, (4, 3)))
    a = axs[7]; ptitle(a, "Lower-level core (m)", "above zero the core is warm; strip below: structure class")
    a.plot(hrs, fx.VTL, color=TEAL, **LW); a.axhline(0, color=TEAL, lw=0.9, ls=(0, (4, 3)))
    lo, hi = a.get_ylim(); a.set_ylim(lo - 0.12 * (hi - lo), hi)
    for h_, p in zip(hrs, fx.phase):
        if isinstance(p, str):
            a.axvspan(h_ - 3, h_ + 3, ymin=0, ymax=0.07, color=PHASES[p][1], lw=0, zorder=4)
    a.set_xlabel(f"hours from lowest pressure ({t0:%d %b %Y %H}Z)", labelpad=6)
    fig.legend(handles=[Line2D([], [], color=c, lw=7, label=n) for n, c in PHASES.values()] + [plt.Rectangle((0, 0), 1, 1, color="#fbe7c2", label="24 h of fastest fall")],
               fontsize=8, loc="lower center", ncol=5, bbox_to_anchor=(0.5, 0.045), handlelength=1.2, columnspacing=1.2)
    fig.suptitle(title, fontsize=12, fontweight="semibold", x=0.01, ha="left", y=0.915)
    fig.savefig(path, dpi=100, bbox_inches="tight"); plt.close(fig)


def quarter_shares(cells):
    hf = cells[cells.g >= HF]
    return [((np.sign(hf.x) == sx) & (np.sign(hf.y) == sy)).mean() * 100 if len(hf) else np.nan for _, sx, sy in QUAD]


def fig_structure(fx, t0, st, cells, title, path):
    fig = plt.figure(figsize=(14.2, 5.0))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1.15, 0.8], wspace=0.28)
    ax = fig.add_subplot(gs[0]); e = fx.dropna(subset=["B", "VTL"])
    ax.axhline(B_EDGE, color="#d6d3cb", lw=1); ax.axvline(0, color="#d6d3cb", lw=1); ax.spines["bottom"].set_visible(True)
    ax.plot(e.VTL, e.B, color="#c4c2ba", lw=1.4, zorder=2)
    hrs = (e.t - t0) / pd.Timedelta(hours=1)
    sc = ax.scatter(e.VTL, e.B, c=hrs, cmap="coolwarm", vmin=-60, vmax=60, s=52, edgecolor="white", linewidth=1, zorder=3)
    i0 = (hrs.abs()).idxmin(); ax.scatter([e.VTL[i0]], [e.B[i0]], marker="*", s=300, color="#fff7d6", edgecolor=old.INK, linewidth=1.2, zorder=5)
    for k in e.index[::3]:
        ax.annotate(e.t[k].strftime("%d/%H"), (e.VTL[k], e.B[k]), xytext=(5, 5), textcoords="offset points", fontsize=7, color=old.MUTED)
    ax.text(0.02, 0.97, "frontal, cold core", transform=ax.transAxes, fontsize=8.5, va="top", color=PHASES["fcc"][1], fontweight="semibold")
    ax.text(0.98, 0.97, "frontal, warm core", transform=ax.transAxes, fontsize=8.5, va="top", ha="right", color=PHASES["fwc"][1], fontweight="semibold")
    ax.text(0.98, 0.03, "warm seclusion", transform=ax.transAxes, fontsize=8.5, ha="right", color=PHASES["wsec"][1], fontweight="semibold")
    ax.text(0.02, 0.03, "symmetric, cold core", transform=ax.transAxes, fontsize=8.5, color=PHASES["scc"][1], fontweight="semibold")
    ax.set_xlabel("lower-level core (m): warm to the right"); ax.set_ylabel("B (m): frontal above 10 m")
    cb = fig.colorbar(sc, ax=ax, location="bottom", shrink=0.7, pad=0.16, aspect=30); cb.set_label("hours from lowest pressure (star)"); cb.outline.set_visible(False); cb.ax.tick_params(length=0)
    ax.set_title("Where the structure went")
    ax = fig.add_subplot(gs[1]); R = 600
    hf = cells[cells.g >= HF]
    for r in (200, 400, 600):
        ax.add_patch(plt.Circle((0, 0), r, fill=False, ec="#d6d3cb", lw=0.8)); ax.text(r * 0.707 + 6, r * 0.707 + 6, f"{r} nm", fontsize=7, color=old.MUTED)
    bins = np.arange(-R, R + 1, 40)
    both = cells[cells.g >= 64]
    if len(both):
        ax.hist2d(both.x, both.y, bins=bins, cmap=LinearSegmentedColormap.from_list("o", ["#fff1de", "#ffc58a"]), cmin=1, vmin=0, vmax=max(30, np.percentile(np.histogram2d(both.x, both.y, bins=bins)[0], 99)) * 1.5)
    if len(hf):
        ax.hist2d(hf.x, hf.y, bins=bins, cmap=old.CM_GUST, cmin=1, vmin=-8, vmax=max(20, np.percentile(np.histogram2d(hf.x, hf.y, bins=bins)[0], 99)) * 1.1, zorder=3)
    ax.axhline(0, color="#d6d3cb", lw=0.8, zorder=1); ax.axvline(0, color="#d6d3cb", lw=0.8, zorder=1)
    ax.annotate("", xy=(R * 0.92, 0), xytext=(R * 0.6, 0), arrowprops=dict(arrowstyle="-|>", color=old.INK, lw=1.4), zorder=6); ax.text(R * 0.92, 14, "motion", fontsize=8, ha="right", color=old.INK)
    ax.scatter([0], [0], marker="o", s=60, color=old.INK, zorder=6)
    ax.set_xlim(-R, R); ax.set_ylim(-R, R); ax.set_aspect(1)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_xlabel("nm ahead (+) or behind (-)"); ax.set_ylabel("nm right (+) or left (-) of motion")
    ax.set_title("Where the strong gusts were (ocean, every 3 h)")
    ax.text(0.0, -0.17, "peach: 64 kt+    red: 71.7 kt+, darker = more samples", transform=ax.transAxes, fontsize=8, color=old.MUTED)
    ax = fig.add_subplot(gs[2]); ax.axis("off")
    ax.set_title("Share of 71.7 kt+ cells by quarter")
    n = len(hf); shares = []
    for name, sx, sy in QUAD:
        shares.append(((np.sign(hf.x) == sx) & (np.sign(hf.y) == sy)).mean() * 100 if n else np.nan)
    top = max([s for s in shares if np.isfinite(s)] + [1])
    for k, ((name, sx, sy), s) in enumerate(zip(QUAD, shares)):
        y0 = 0.78 - 0.17 * k; w = 0 if not np.isfinite(s) else s / 100
        ax.text(0.0, y0 + 0.075, name, transform=ax.transAxes, fontsize=9, color=old.MUTED)
        ax.add_patch(plt.Rectangle((0.0, y0), 0.78, 0.045, transform=ax.transAxes, color="#f1efe9"))
        ax.add_patch(plt.Rectangle((0.0, y0), 0.78 * w, 0.045, transform=ax.transAxes, color="#e5502f" if s == top else "#f0a98f"))
        ax.text(0.82, y0 + 0.005, "n/a" if not np.isfinite(s) else f"{s:.0f}%", transform=ax.transAxes, fontsize=10, fontweight="semibold")
    if n:
        rr = np.hypot(hf.x, hf.y)
        ax.text(0.0, 0.06, f"Median distance from the centre {np.median(rr):.0f} nm\n(90% within {np.percentile(rr, 90):.0f} nm)", transform=ax.transAxes, fontsize=8.5, va="top", color=old.MUTED)
    fig.suptitle(title, fontsize=12, fontweight="semibold", x=0.01, ha="left", y=1.03)
    fig.savefig(path, dpi=100, bbox_inches="tight"); plt.close(fig)
    return shares


def fig_maps_low(z, polys, title, path, fx, t0, oro):
    fig, axs = plt.subplots(2, 4, figsize=(17, 8.6), gridspec_kw=dict(wspace=0.08, hspace=0.22))
    axs = axs.ravel(); cf = None
    for k, h in enumerate(MAP_H):
        ax = axs[k]; la, lo = z[f"centre_{h}"]
        lats, lons = old.crop_axes(la, lo, old.DLAT, old.DLON)
        high = ex.crop(oro, la, lo, old.DLAT, old.DLON)[0] > TERRAIN_M          # same cells as the stored maps
        t850 = np.where(high, np.nan, z[f"t850_{h}"]); ws = np.where(high, np.nan, z[f"ws925_{h}"])
        cf = ax.contourf(lons, lats, t850, levels=np.arange(-30, 21, 2), cmap="coolwarm", extend="both", zorder=1)
        ax.contour(lons, lats, t850, levels=[0], colors="#555", linewidths=0.8, linestyles="--", zorder=2)
        ax.contourf(lons, lats, high.astype(float), levels=[0.5, 1.5], colors=["#d9d7d0"], zorder=1.5)
        ax.contour(lons, lats, z[f"msl_{h}"], levels=np.arange(900, 1060, 4), colors=old.INK, linewidths=0.55, alpha=0.8, zorder=3)
        if np.nanmax(ws) >= 50:
            ax.contour(lons, lats, ws, levels=[50, 64], colors=["#7c3aed", "#7c3aed"], linewidths=[1.2, 2.6], zorder=4)
        th = t0 + pd.Timedelta(hours=h); tracked = fx.t.iloc[0] <= th <= fx.t.iloc[-1]
        ax.plot([lo], [la], marker="x" if tracked else "+", color=old.INK if tracked else "#777", ms=9, mew=2, zorder=5)
        ax.set_xlim(lons[0], lons[-1]); ax.set_ylim(lats[-1], lats[0]); ax.set_aspect(1 / np.cos(np.radians(la)))
        old.draw_coast(ax, polys, 0.7)
        ax.set_title(f"{'peak' if h == 0 else f'{h:+d} h'}: lowest {np.nanmin(z[f'msl_{h}']):.0f} hPa, 925 hPa wind up to {np.nanmax(ws):.0f} kt", fontsize=8.8)
        ax.tick_params(labelsize=7); old.lonfmt(ax)
    ax = axs[7]; ax.axis("off")
    cb = fig.colorbar(cf, ax=ax, location="top", shrink=0.9, pad=0.0); cb.set_label("850 hPa temperature (deg C)")
    ax.text(0, 0.55, "colour: 850 hPa temperature (the dashed\nline is 0 C): the tight packing of colour\nmarks the frontal zone\nblack: sea-level pressure every 4 hPa\npurple: 925 hPa wind, thin 50 kt, thick 64 kt\ngrey: ground above 1,000 m, where the 850 and\n925 hPa levels are not meaningful\n\nERA5 proxy at 0.25 deg, 6-hourly\nx: storm position on the pipeline A track", transform=ax.transAxes, fontsize=8, va="top")
    fig.suptitle(title, fontsize=11, x=0.01, ha="left", y=0.955)
    fig.savefig(path, dpi=85, bbox_inches="tight"); plt.close(fig)


# ------------------------------------------------------------------ text

def hhmm(t):
    return pd.Timestamp(t).strftime("%d %b %H00Z")


def compass(x, y):
    fb = "ahead of" if x >= 0 else "behind"; lr = "right" if y >= 0 else "left"
    return f"{abs(x):.0f} nm {fb} and {abs(y):.0f} nm to the {lr} of the centre"


def meteo(f, fx, d24, st, cells, shares, opc, basin):
    """Return the generated text blocks of one review (HTML strings)."""
    t0 = pd.Timestamp(f["t_minp"]); j = int(np.nanargmin(d24)); tend = fx.t[j]; tstart = tend - pd.Timedelta(hours=24)
    win = fx[(fx.t >= tstart) & (fx.t <= tend)]
    out = {}
    # phase history
    runs = []
    for t, p in zip(fx.t, fx.phase):
        if isinstance(p, str):
            if not runs or runs[-1][0] != p:
                runs.append([p, t, t])
            else:
                runs[-1][2] = t
    first = fx.dropna(subset=["phase"]).iloc[0] if fx.phase.notna().any() else None
    ws = [r for r in runs if r[0] == "wsec"]
    pk = fx.iloc[(fx.t - t0).abs().argmin()]
    s = ""
    if first is not None:
        s += f"At its first tracked fix ({hhmm(first.t)}, {first.lat:.1f}&deg;N) the low was <b>{PHASES[first.phase][0]}</b> (Hart B {first.B:.0f} m, lower-level core {first.VTL:+.0f} m). "
    if isinstance(pk.phase, str):
        s += f"At its lowest pressure it was <b>{PHASES[pk.phase][0]}</b> (B {pk.B:.0f} m, lower-level core {pk.VTL:+.0f} m). "
    if ws:
        s += f"The first fix with a symmetric warm core (the warm-seclusion signature) is {hhmm(ws[0][1])}, {(ws[0][1] - t0) / pd.Timedelta(hours=1):+.0f} h from the lowest pressure. "
    else:
        s += "It never reached the symmetric warm-core signature of a warm seclusion in these 1.5&deg; diagnostics. "
    out["phase"] = s
    # forcing
    s = (f"During its fastest 24 h fall ({dt24(win)}) the strongest 250 hPa wind within about 540 nm averaged <b>{win.jet250.mean():.0f} kt</b>. "
         f"Upper-level divergence at 300 hPa over the centre averaged {win.div300.mean():+.1f} &times;10<sup>&minus;5</sup>/s (positive means air spreading out aloft, which draws air up from below), "
         f"and 500 hPa vorticity advection {win.vadv500.mean():+.1f} &times;10<sup>&minus;9</sup>/s&sup2; (positive means cyclonic spin arriving aloft). "
         f"The Eady growth rate (850&ndash;500 hPa), a measure of how favourable the temperature contrast and wind shear are for the storm to grow, averaged {win.eady.mean():.2f} per day, "
         f"and the moisture in the air column within about 270 nm averaged {win.tcwv.mean():.0f} kg/m&sup2;.")
    out["forcing"] = s
    # wind
    stx = st.copy(); pkr = stx.loc[stx.area72_nm2.idxmax()] if (stx.area72_nm2 > 0).any() else None
    s = ""
    if pkr is not None:
        s += (f"The ocean area with gusts of 71.7 kt or more within 600 nm peaked at <b>{pkr.area72_nm2:,.0f} nm&sup2;</b> at {hhmm(t0 + pd.Timedelta(hours=pkr.h))} "
              f"({pkr.h:+.0f} h from the lowest pressure). ")
        n72 = int((stx.area72_nm2 > 0).sum())
        s += f"Some such wind was present in {n72} of the {len(stx)} 3-hourly fields (about {n72 * 3} h). "
    mx = stx.loc[stx.maxg_kt.idxmax()]
    s += (f"The single highest ocean gust in a 3-hourly field was {mx.maxg_kt:.0f} kt at {hhmm(t0 + pd.Timedelta(hours=mx.h))}, "
          f"{compass(mx.max_fwd_nm, mx.max_right_nm)}. ")
    if np.all(np.isfinite(shares)):
        k = int(np.argmax(shares))
        s += f"Of the 71.7 kt+ cells, {shares[k]:.0f}% sat in the <b>{QUAD[k][0]}</b> quarter of the storm (frame fixed to the direction of motion; where the low was stalled the last direction it moved was used)."
    out["wind"] = s
    # setup percentile table
    hfx = fx[fx.gust >= HF]
    rows = []
    for var, lab, unit, nd in VARS:
        v = hfx[var].median() if len(hfx) >= 1 else np.nan
        r = ref_rank(basin, var, v)
        rows.append((lab, unit, v, nd, r))
    out["rank_rows"] = rows
    out["rank_n"] = len(hfx)
    # opc
    if opc:
        s = (f"OPC's archive carries this low as event {opc['id']}: {opc['n_fix']} of its fixes sit on this low (within 200 nm of the track), {opc['n_hf']} of them at hurricane force"
             + (f" ({hhmm(opc['hf_first'])} to {hhmm(opc['hf_last'])})" if opc["n_hf"] else "")
             + f". Its lowest analysed pressure was {opc['minp']:.0f} hPa at {hhmm(opc['t_minp'])}, against {f['minp']:.1f} hPa from ERA5 at {hhmm(f['t_minp'])}. "
             f"At the {opc['era_hf_in_opc']} of the proxy's {opc['n_era_hf']} hurricane-force-index fixes that OPC also has a fix for, OPC carried hurricane force on {opc['n_both']}. "
             f"Matched fixes are a median {opc['dist_med']:.0f} nm apart (OPC position against the pipeline A track).")
        out["opc"] = s
    return out


def dt24(win):
    return f"{hhmm(win.t.iloc[0])} to {hhmm(win.t.iloc[-1])}"


# ------------------------------------------------------------------ page

def fnum(v, nd):
    return "n/a" if not np.isfinite(v) else f"{v:.{nd}f}"


def ordinal(p):
    return f"higher than {p:.0f}% of them"


def page(f, fx, tx, opc, meteo_, shares, style, base):
    esc = html.escape
    cls = "atl" if f["basin"] == "atl" else "pac"
    b, lab = old.label(f)
    pre = pd.Timestamp(f["t_minp"]).year < 2004
    dur = ((pd.Timestamp(f["last_hf"]) - pd.Timestamp(f["first_hf"])) / pd.Timedelta(hours=1) + 6) if f["first_hf"] else 0
    nm = f' &mdash; the {esc(f["name"])} <span class="muted">(name inferred from date and position)</span>' if f["name"] else ""
    rows = "".join(f'<tr><td>{l}</td><td class="num">{fnum(v, nd)} {u}</td><td class="num">{"n/a" if r is None else ordinal(r)}</td></tr>' for l, u, v, nd, r in meteo_["rank_rows"])
    tiles = (f'<div class="tiles"><div><b>{f["minp"]:.1f}</b><span>lowest pressure (hPa)</span></div><div><b>{old.fall(f)}</b><span>fastest 24 h fall (hPa)</span></div>'
             f'<div><b>{f["gust_idx"]:.0f}</b><span>peak gust index (kt)</span></div><div><b>{dur:.0f}</b><span>hours at HF strength</span></div>'
             f'<div><b>{tx["area"]:,.0f}</b><span>peak ocean area at 71.7 kt+ (nm&sup2;)</span></div></div>')
    pre_note = ""
    if pre:
        pre_note = ('<p class="limit"><strong>Not direct observation.</strong> This storm is before the 2004&ndash;05 season, so no archive record covers it; every number and map '
                    'here is ERA5 reanalysis. ' + ('The gust index also runs high before 2001, so read its level with caution; central pressure is the better yardstick.' if f["gust_flag"] else '') + '</p>')
    opc_block = (f'<div class="blk"><div class="tag">What OPC had</div><p>{meteo_["opc"]}</p><p class="muted">OPC&rsquo;s call is analyst judgement of sustained wind; the gust index is a proxy for it, so &ldquo;OPC also carried it&rdquo; is agreement, not verification.</p></div>'
                 if "opc" in meteo_ else
                 '<div class="blk"><div class="tag">What OPC had</div><p>No OPC archive fix matches this storm (the archive starts with the 2004&ndash;05 season), so there is nothing to compare with.</p></div>')
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Storm review: {esc(lab)} - HF Low Archive</title>
<link rel="icon" type="image/svg+xml" href="favicon.svg">
{style}
</head>
<body>
<header class="masthead"><div class="masthead-in"><h1>Storm review: {esc(lab)}</h1><a href="storm-reviews.html">&larr; All storm reviews</a></div></header>
<main>
<div class="note"><strong>ERA5 proxy, not direct observation.</strong> ERA5 is a weather model re-run over past years and blended with the observations available. Tracks are the project&rsquo;s pipeline A tracks, found in ERA5 itself. The gust index is the strongest ERA5 10 m gust over the ocean within 800 km (about 430 nm) of the low, a stand-in (proxy) for the sustained wind OPC warns on. Structure and upper-level numbers are ERA5 at 1.5&deg;; maps and the wind field are ERA5 at 0.25&deg;. Distances are nautical miles.</div>
<section id="storm" class="storm {cls}">
  <h3><span class="rank">{b.replace("North ", "")} #{f["rank"]}</span><span class="basin-{cls}">{esc(lab)}</span>{nm}</h3>
  <div class="meta">Rank {f["rank"]} by lowest ERA5 central pressure in the {b}. Pipeline A track {f["track"]}.</div>
  {tiles}
  {pre_note}
  <div class="sum"><div class="tag">1. The storm in plain words</div><p>{old.summary(f)}</p></div>
</section>
<section><div class="tag">2. Track and wind field</div>
<figure class="fig"><div class="scroll"><a href="{base}-track.png"><img src="{base}-track.png" alt="Track and gust swath" loading="lazy"></a></div>
<figcaption>Left: 6-hourly fixes coloured by central pressure, red rings where the gust index reached 71.7 kt. Right: highest ERA5 10 m gust at each point from 48 h before to 24 h after the lowest pressure (land included; black contour 71.7 kt).</figcaption></figure></section>
<section><div class="tag">3. How it developed</div>
<p>{meteo_["phase"]}</p>
<p class="muted">How to read the structure numbers. Hart&rsquo;s B says how lopsided the low is: above about 10 m it is frontal, near zero it is symmetric. The lower-level core number says whether the 900 to 600 hPa layer is cold (negative) or warm (positive). A frontal low with a cold core is the early stage; a warm seclusion has warm air wrapped around the centre. The classes come from a coarse 1.5&deg; grid, so treat them as a guide.</p>
<figure class="fig"><div class="scroll"><a href="{base}-series.png"><img src="{base}-series.png" alt="Six-panel development series" loading="lazy"></a></div>
<figcaption>From the top: central pressure; the fall in pressure over the previous 24 h against the bomb yardstick (scaled by latitude); the gust index; the ocean area with gusts of 71.7 kt or more; the 250 hPa jet; the two upper-level lift terms (divergence and vorticity advection); and two structure measures (frontal against symmetric, cold core against warm core) with the structure class along the bottom. The tan band is the 24 h of fastest fall.</figcaption></figure></section>
<section><div class="tag">4. What drove it: upper-level support</div>
<p>{meteo_["forcing"]}</p>
<figure class="fig"><div class="scroll"><a href="{base}-maps.png"><img src="{base}-maps.png" alt="Upper-level maps" loading="lazy"></a></div>
<figcaption>Sea-level pressure (black), 500 hPa height (blue) and 250 hPa wind speed (shaded jet), every 12 h from 48 h before to 24 h after the lowest pressure.</figcaption></figure></section>
<section><div class="tag">5. Low-level structure: fronts and the low-level jet</div>
<figure class="fig"><div class="scroll"><a href="{base}-lower.png"><img src="{base}-lower.png" alt="850 hPa temperature and 925 hPa wind maps" loading="lazy"></a></div>
<figcaption>850 hPa temperature shows the fronts and the warm sector; the purple lines enclose 925 hPa wind of 50 and 64 kt, where the strongest low-level wind sat in each frame.</figcaption></figure></section>
<section><div class="tag">6. Where the strongest wind was</div>
<p>{meteo_["wind"]}</p>
<figure class="fig"><div class="scroll"><a href="{base}-structure.png"><img src="{base}-structure.png" alt="Hart phase space and storm-relative wind" loading="lazy"></a></div>
<figcaption>Left: the storm&rsquo;s path through Hart&rsquo;s cyclone phase space, a diagram of how lopsided the low is against how warm its core is. Middle: ocean cells at 64 kt or more, in a frame that moves and turns with the low. Right: how the 71.7 kt+ cells split between the four quarters.</figcaption></figure></section>
<section><div class="tag">7. Was its set-up unusual? Compared with other hurricane-force-strength lows</div>
<p>For this storm&rsquo;s {meteo_["rank_n"]} fixes with a gust index of 71.7 kt or more, the table gives the median of each ingredient and how it ranks among all such fixes in the same basin since 2004&ndash;05 (pipeline A, ERA5 proxy). Fixes from one storm are not independent, so read the rank as a rough guide, not a test.</p>
<div class="tablewrap"><table class="idx"><thead><tr><th>Ingredient</th><th class="num">This storm</th><th class="num">Higher than this share of HF-strength fixes</th></tr></thead><tbody>{rows}</tbody></table></div></section>
<section>{opc_block}</section>
<section><div class="tag">Read with care</div>
<p class="limit">ERA5 proxy throughout; never validated for storms before the archive begins. The structure measures (Hart&rsquo;s cyclone phase space) use seven pressure levels on a 1.5&deg; grid, coarser than the original method, so read the structure class as a guide, not an exact classification. Over Greenland and high ground ERA5&rsquo;s sea-level pressure is extrapolated, and gusts near coasts and sea ice are less reliable. The storm-relative frame uses the last direction of motion while the low is stalled. A storm is one case: nothing here is a trend or a forecast rule.</p></section>
</main></body></html>"""


INDEX_METHOD = """<section id="how">
  <details>
  <summary>How the storms were chosen and drawn (method)</summary>
  <p>The storms are the ones with the lowest ERA5 central pressure in each basin among the 1979&ndash;2025 tracks of pipeline A (the project&rsquo;s ERA5 storm tracker, <code>research/era5/hf_history</code>), taken in order of pressure from the lowest. Tropical-cyclone-flagged tracks are left out (three Pacific tracks, among them 30 Oct 1989 and 3 Oct 1981, which would otherwise rank second and third), and one Pacific storm from the 2024&ndash;25 season is left out because the weather-map fields used here end in January 2023. The rank shown is the rank among all tracks, so some numbers are missing. Ranking is by pressure alone, so the list is not a list of the strongest winds.</p>
  <p>Central pressure, position and gust index come from the 6-hourly pipeline A fixes. The gust index is the highest ERA5 10 m gust over ocean within 800 km (about 430 nm) of the low. Each review also uses WeatherBench2 ERA5 at 1.5&deg; for the 250 hPa jet, 300 hPa divergence, 500 hPa vorticity advection and Hart&rsquo;s cyclone-phase-space structure measures (frontal against symmetric, cold core against warm core), and at 0.25&deg; for the maps (sea-level pressure, 500 hPa height, 250 hPa wind, 850 hPa temperature, 925 hPa wind) re-centred on the storm. The gust swath and the storm-relative wind field use ERA5 hourly instantaneous gusts every 3 h over the 72 h window. The 24 h fall is the change in central pressure over the previous 24 h at each 6-hourly fix; the bomb line is 24 hPa in 24 h at 60&deg;N, scaled for other latitudes by sin(latitude)/sin(60&deg;). OPC warns on sustained winds; the gust index is a proxy for that, not the same thing. Where OPC&rsquo;s archive has the storm, its fixes are matched by time and position (within 200 nm) and shown as agreement, not verification.</p>
  <p class="limit"><strong>Read with care.</strong> ERA5 proxy throughout, never validated for storms before the archive begins. Before 2001 the gust index drifts upward at fixed storm depth, so gust numbers for the 1980s and 1990s storms are less reliable than their pressures. Not independently re-checked: how closely these ERA5 pressures match the observed minima of the named storms.</p>
  </details>
</section>"""


def write_index(rows, outdir):
    style = old.page_style(); esc = html.escape
    n = len(rows); na = sum(f["basin"] == "atl" for f in rows); npc = n - na
    npre = sum(pd.Timestamp(f["t_minp"]).year < 2004 for f in rows)
    deep = min(rows, key=lambda f: f["minp"]); fast = max(rows, key=lambda f: f["dp24_max"]); nhf = sum(1 for f in rows if f["first_hf"])
    gmax = max(f["gust_idx"] for f in rows); tr = []
    for f in sorted(rows, key=lambda f: (f["basin"], f["minp"])):
        b, lab = old.label(f); cls = "atl" if f["basin"] == "atl" else "pac"
        pre = pd.Timestamp(f["t_minp"]).year < 2004
        dur = ((pd.Timestamp(f["last_hf"]) - pd.Timestamp(f["first_hf"])) / pd.Timedelta(hours=1) + 6) if f["first_hf"] else 0
        nm = f' <span class="muted">{esc(f["name"].split(",")[0])}</span>' if f["name"] else ""
        tr.append(f'<tr data-basin="{cls}"><td><a href="storm-review-{f["track"]}.html"><span class="basin-{cls}">{esc(lab)}</span></a>{nm}</td><td class="num">{f["minp"]:.1f}</td>'
                  f'<td class="num">{old.fall(f)}</td><td class="num"><meter class="bar {cls}" min="0" max="{gmax:.0f}" value="{f["gust_idx"]:.0f}" aria-hidden="true"></meter> {f["gust_idx"]:.0f}</td>'
                  f'<td class="num">{dur:.0f}</td><td>{"reanalysis only" if pre else "archive era"}</td></tr>')
    page = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Storm reviews - Hurricane Force Extratropical Low Archive</title>
<link rel="icon" type="image/svg+xml" href="favicon.svg">
<meta name="description" content="Storm reviews of the {n} deepest extratropical lows in the ERA5 proxy record, {na} Atlantic and {npc} Pacific: development, jet and structure diagnostics, fronts, where the strong wind was, and what OPC had.">
<!-- generated by tools/storm_review_v2.py from tools/storm_review_facts.csv: edit the script, not this file -->
{style}
</head>
<body>
<header class="masthead"><div class="masthead-in">
  <h1>Storm reviews: the strongest lows, one at a time</h1>
  <a href="index.html">&larr; Back to the archive</a>
</div></header>
<main>
<div class="note"><strong>ERA5 proxy, not direct observation.</strong> ERA5 is a weather model re-run over past years and blended with the observations available. Every storm here was found in ERA5 by the project&rsquo;s pipeline A tracker, and every chart is drawn from ERA5 fields. {npre} of the {n} are before the 2004&ndash;05 season, where no archive record exists; the other {n - npre} are since then, and are still read here from ERA5. ERA5 smooths the very strongest winds and pressure minima, so real peaks may be stronger. These are case studies, not a trend: nothing here says whether such storms are becoming more or less common.</div>
<div class="answers">
  <div class="tag">What this page shows</div>
  <p>The {n} deepest lows ERA5 found, {na} in the Atlantic and {npc} in the Pacific, ranked by lowest central pressure. Each has its own page that walks through how the storm developed: how fast it deepened, what the upper-level pattern looked like, where the fronts and the low-level jet were, which part of the storm had the strongest wind, and whether OPC carried it as hurricane force. Use it to compare the setup of a storm you are watching with how the strongest storms on record developed. The list is ranked by pressure, so it is not a list of the strongest winds.</p>
</div>
<div class="hero">
  <div class="stat"><b>{n}</b><span>storms: {na} North Atlantic, {npc} North Pacific</span></div>
  <div class="stat"><b>{deep["minp"]:.1f} hPa</b><span>lowest: {esc(old.label(deep)[1])}</span></div>
  <div class="stat"><b>{fast["dp24_max"]:.0f} hPa</b><span>fastest 24 h fall: {esc(old.label(fast)[1])}</span></div>
  <div class="stat"><b>{nhf} of {n}</b><span>reach a gust index of 71.7 kt or more</span></div>
</div>
<div class="answers">
  <div class="tag">Open a review</div>
  <p class="muted">Click a storm to open its review. &ldquo;Gust index&rdquo; is the strongest ERA5 10 m gust over the ocean within about 430 nm of the low. &ldquo;Archive era&rdquo; means the storm is since the 2004&ndash;05 season, so OPC&rsquo;s archive can be compared; &ldquo;reanalysis only&rdquo; means it is earlier.</p>
  <div class="toolbar" role="group" aria-label="Show storms">
    <span>Show:</span>
    <button type="button" data-f="all" aria-pressed="true">All {n}</button>
    <button type="button" data-f="atl" aria-pressed="false">North Atlantic ({na})</button>
    <button type="button" data-f="pac" aria-pressed="false">North Pacific ({npc})</button>
  </div>
  <div class="tablewrap"><table class="idx">
  <thead><tr><th>Storm</th><th class="num">Lowest pressure (hPa)</th><th class="num">Fastest 24 h fall (hPa)</th><th class="num">Peak gust index (kt)</th><th class="num">Hours at HF strength</th><th>Source</th></tr></thead>
  <tbody id="rows">{"".join(tr)}</tbody></table></div>
  <p class="muted">Hours at HF strength is the number of 6-hourly positions with a gust index of at least 71.7 kt, times six. Peak gust index is the highest gust-index value over the storm&rsquo;s life.</p>
</div>
{INDEX_METHOD}
<footer><p class="muted">Seasons run 1 June to 31 May. Fields: WeatherBench2 ERA5 and ARCO-ERA5 (Google Cloud public data). <a href="index.html">Back to the archive</a></p></footer>
</main>
<script>
  (function () {{
    try {{ var saved = localStorage.getItem('hf-theme'); if (saved) document.documentElement.setAttribute('data-theme', saved); }} catch (err) {{ /* private mode */ }}
  }})();
  (function () {{
    var btns = document.querySelectorAll('.toolbar button'), rows = document.querySelectorAll('#rows tr');
    btns.forEach(function (b) {{
      b.addEventListener('click', function () {{
        btns.forEach(function (o) {{ o.setAttribute('aria-pressed', o === b ? 'true' : 'false'); }});
        rows.forEach(function (r) {{ r.hidden = b.dataset.f !== 'all' && r.dataset.basin !== b.dataset.f; }});
      }});
    }});
  }})();
</script>
</body>
</html>
"""
    with open(os.path.join(outdir, "storm-reviews.html"), "w") as fh:
        fh.write("\n".join(l.rstrip() for l in page.splitlines()) + "\n")


def quantize(outdir, rows):
    """160-colour PNGs keep the figure set small (the repository holds the PNGs)."""
    from PIL import Image
    for r in rows:
        for k in ("track", "series", "maps", "lower", "structure"):
            fp = os.path.join(outdir, f"storm-review-{r['track']}-{k}.png")
            Image.open(fp).convert("RGB").quantize(colors=160, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE).save(fp, optimize=True)


def build(work, fixes_csv, outdir, tracks):
    polys = old.coast()
    st_all = pd.read_csv(os.path.join(HERE, "storm_review_storms.csv"))
    fx_all = pd.read_csv(fixes_csv); fx_all["t"] = pd.to_datetime(fx_all.time.astype(str), format="%Y%m%d%H")
    style = old.page_style()
    oro = orography(work)
    rows = []
    for tr in tracks:
        row = st_all[st_all.track == int(tr)].iloc[0]
        fx = old.prep(fx_all[fx_all.track == row.track].copy())
        z = np.load(os.path.join(work, f"{row.track}.npz"))
        t0 = fx.t[fx.msl.idxmin()]
        fx = track_env(work, fx, int(row.track))
        d24 = old.deepening(fx)
        f = old.facts(row, fx, d24, t0, z["sw_gust"])
        st, cells = sr_data(z)
        base = os.path.join(outdir, f"storm-review-{row.track}")
        b, lab = old.label(f); head = f"{lab} (ERA5 proxy, pipeline A)"
        slon = np.degrees(np.unwrap(np.radians(z["sw_lons"]))); slon = slon + 360 * np.round((fx.lonu.mean() - slon.mean()) / 360)
        if os.environ.get("HF_HTML_ONLY"):                  # text only: figures already drawn
            shares = quarter_shares(cells)
        else:
            old.fig_track(fx, z["sw_gust"], z["sw_lats"], slon, polys, head, base + "-track.png")
            fig_series2(fx, d24, t0, st, head, base + "-series.png")
            old.fig_maps(z, polys, head + ": every 12 hours, 48 h before to 24 h after the lowest pressure", base + "-maps.png", fx, t0)
            fig_maps_low(z, polys, head + ": 850 hPa temperature and 925 hPa wind", base + "-lower.png", fx, t0, oro)
            shares = fig_structure(fx, t0, st, cells, head, base + "-structure.png")
        opc = opc_match(fx, f["basin"]) if t0.year >= 2004 else None
        m = meteo(f, fx, d24, st, cells, shares, opc, f["basin"])
        tx = dict(area=float(st.area72_nm2.max()))
        with open(base + ".html", "w") as fh:
            fh.write("\n".join(l.rstrip() for l in page(f, fx, tx, opc, m, shares, style, os.path.basename(base)).splitlines()) + "\n")
        rows.append(f); print(row.track, lab, flush=True)
    return rows


if __name__ == "__main__":
    if sys.argv[4:] == ["--all"]:                      # every storm in storm_review_storms.csv, then facts, index and small PNGs
        tracks = [str(t) for t in pd.read_csv(os.path.join(HERE, "storm_review_storms.csv")).track]
        rows = build(sys.argv[1], sys.argv[2], sys.argv[3], tracks)
        pd.DataFrame(rows).to_csv(os.path.join(HERE, "storm_review_facts.csv"), index=False)
        if not os.environ.get("HF_HTML_ONLY"):
            quantize(sys.argv[3], rows)
        write_index(rows, sys.argv[3])
    else:
        build(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:])
