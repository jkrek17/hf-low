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
PHASES = {"fcc": ("frontal, cold core", "#2a78d6"), "fwc": ("frontal, warm core", "#e0a21b"),
          "wsec": ("warm seclusion (symmetric, warm core)", "#d63c2a"), "scc": ("symmetric, cold core", "#7a7a76")}
B_EDGE = 10.0                                   # m, Hart's asymmetry threshold
QUAD = (("front right", 1, 1), ("rear right", -1, 1), ("rear left", -1, -1), ("front left", 1, -1))
VARS = [("jet250", "highest 250 hPa wind within 1,000 km", "kt", 0), ("div300", "300 hPa divergence within 500 km", "1e-5 /s", 1),
        ("vadv500", "500 hPa vorticity advection within 500 km", "1e-9 /s2", 1), ("eady", "Eady growth rate (850-500 hPa)", "/day", 2),
        ("tcwv", "column water vapour within 500 km", "kg/m2", 0), ("sstgrad", "largest SST gradient within 500 km", "K/100 km", 1)]


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
    o = d[d.ID == best].sort_values("t")
    mm = m[m.ID == best]
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

def shade_fall(axs, thr, hrs_end):
    for a in axs:
        a.axvspan(hrs_end - 24, hrs_end, color="#f3d9a8", alpha=0.45, lw=0, zorder=0)


def fig_series2(fx, d24, t0, st, title, path):
    hrs = (fx.t - t0) / pd.Timedelta(hours=1)
    j = int(np.nanargmin(d24)); end = hrs[j]
    fig, axs = plt.subplots(6, 1, figsize=(9.2, 13.2), sharex=True, gridspec_kw=dict(hspace=0.14, height_ratios=[1, 1, 1, 1, 1, 1.15]))
    shade_fall(axs, None, end)
    for a in axs:
        a.axvline(0, color="k", lw=0.8, ls=":"); a.grid(color="#ccc", lw=0.5)
    axs[0].plot(hrs, fx.msl, "-o", ms=3, color="#0a3f73"); axs[0].set_ylabel("central pressure\n(hPa)"); axs[0].invert_yaxis()
    axs[0].text(0.01, 0.06, "shaded: the 24 h of fastest fall", transform=axs[0].transAxes, fontsize=7.5, color="#8a6a1c")
    dd = -d24
    axs[1].bar(hrs, dd, width=4.5, color=np.where(dd > 0, "#c0392b", "#8aa"))
    axs[1].plot(hrs, 24 * np.sin(np.radians(fx.lat)) / np.sin(np.radians(60)), color="k", lw=1, ls="--", label="'bomb' yardstick at this latitude")
    axs[1].set_ylabel("fall over previous\n24 h (hPa)"); axs[1].legend(fontsize=7.5, loc="upper left")
    a = axs[2]; a.plot(hrs, fx.gust, "-o", ms=3, color="#b9460f"); a.axhline(HF, color="k", lw=1, ls="--")
    a.set_ylabel("gust index (kt)\nmax within 800 km"); a.text(hrs.min(), HF + 1.2, "71.7 kt", fontsize=7.5)
    a2 = a.twinx(); a2.bar(st.h, st.area72_nm2 / 1000, width=2.6, color="#e2b6a0", alpha=0.8, zorder=0); a2.set_ylabel("ocean area at 71.7 kt+\nwithin 600 nm (000 nm2)", color="#8a5a3c", fontsize=8)
    a2.spines["right"].set_visible(True); a.set_zorder(a2.get_zorder() + 1); a.patch.set_visible(False)
    a = axs[3]; a.plot(hrs, fx.jet250, "-o", ms=3, color="#5b3a8c"); a.set_ylabel("250 hPa jet\n(kt, within 1,000 km)")
    a = axs[4]; a.plot(hrs, fx.div300, "-o", ms=3, color="#1f6fbf", label="300 hPa divergence (1e-5 /s)")
    a.plot(hrs, fx.vadv500, "-s", ms=3, color="#c0392b", label="500 hPa vorticity advection (1e-9 /s2)"); a.axhline(0, color="k", lw=0.6)
    a.set_ylabel("upper-level lift\nterms"); a.legend(fontsize=7.5, loc="upper left")
    a = axs[5]; ph = fx.phase
    a.plot(hrs, fx.B, "-o", ms=3, color="#7a3b00", label="B: thickness asymmetry (m)"); a.axhline(B_EDGE, color="#7a3b00", lw=0.8, ls="--")
    a.set_ylabel("Hart B (m)", color="#7a3b00")
    b = a.twinx(); b.plot(hrs, fx.VTL, "-s", ms=3, color="#1c6b5a", label="lower-level core (VTL)"); b.axhline(0, color="#1c6b5a", lw=0.8, ls="--")
    b.set_ylabel("-VT lower (m)\n>0 warm core", color="#1c6b5a"); b.spines["right"].set_visible(True)
    for h_, p in zip(hrs, ph):
        if isinstance(p, str):
            a.axvspan(h_ - 3, h_ + 3, ymin=0, ymax=0.07, color=PHASES[p][1], lw=0)
    a.set_xlabel(f"hours from lowest pressure ({t0:%d %b %Y %H}Z)", labelpad=6)
    fig.legend(handles=[Line2D([], [], color=c, lw=6, label=n) for n, c in PHASES.values()], fontsize=8, loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.065), title="structure class (strip along the bottom of the last panel)", title_fontsize=8)
    fig.suptitle(title, fontsize=11, x=0.01, ha="left", y=0.905)
    fig.savefig(path, dpi=100, bbox_inches="tight"); plt.close(fig)


def fig_structure(fx, t0, st, cells, title, path):
    fig = plt.figure(figsize=(14.2, 4.9))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1.15, 0.85], wspace=0.3)
    # (a) Hart phase space
    ax = fig.add_subplot(gs[0]); e = fx.dropna(subset=["B", "VTL"])
    ax.axhline(B_EDGE, color="#888", lw=0.8); ax.axvline(0, color="#888", lw=0.8)
    ax.plot(e.VTL, e.B, color="#999", lw=1, zorder=2)
    hrs = (e.t - t0) / pd.Timedelta(hours=1)
    sc = ax.scatter(e.VTL, e.B, c=hrs, cmap="coolwarm", vmin=-60, vmax=60, s=40, edgecolor="k", linewidth=0.4, zorder=3)
    i0 = (hrs.abs()).idxmin(); ax.scatter([e.VTL[i0]], [e.B[i0]], marker="*", s=240, color="#fff", edgecolor="k", zorder=5)
    for k in e.index[::3]:
        ax.annotate(e.t[k].strftime("%d/%H"), (e.VTL[k], e.B[k]), xytext=(4, 4), textcoords="offset points", fontsize=6.5)
    xl = ax.get_xlim(); yl = ax.get_ylim()
    ax.text(0.02, 0.97, "frontal, cold core", transform=ax.transAxes, fontsize=8, va="top", color=PHASES["fcc"][1])
    ax.text(0.98, 0.97, "frontal, warm core", transform=ax.transAxes, fontsize=8, va="top", ha="right", color=PHASES["fwc"][1])
    ax.text(0.98, 0.03, "warm seclusion", transform=ax.transAxes, fontsize=8, ha="right", color=PHASES["wsec"][1])
    ax.text(0.02, 0.03, "symmetric, cold core", transform=ax.transAxes, fontsize=8, color=PHASES["scc"][1])
    ax.set_xlabel("-VT lower (m): warm core to the right"); ax.set_ylabel("B (m): asymmetry, frontal above 10 m")
    cb = fig.colorbar(sc, ax=ax, location="bottom", shrink=0.85, pad=0.16); cb.set_label("hours from lowest pressure")
    ax.set_title("Where the structure went (Hart phase space)\nstar = lowest pressure", fontsize=10)
    # (b) storm-relative gust cells
    ax = fig.add_subplot(gs[1]); R = 600
    hf = cells[cells.g >= HF]; sf = cells[(cells.g >= 64) & (cells.g < HF)]
    for r in (200, 400, 600):
        ax.add_patch(plt.Circle((0, 0), r, fill=False, ec="#aaa", lw=0.6, ls=":")); ax.text(0, r + 6, f"{r} nm", fontsize=6.5, color="#777", ha="center")
    bins = np.arange(-R, R + 1, 40)
    both = cells[cells.g >= 64]
    if len(both):
        ax.hist2d(both.x, both.y, bins=bins, cmap="Oranges", cmin=1, vmin=0, vmax=max(30, np.percentile(np.histogram2d(both.x, both.y, bins=bins)[0], 99)) * 1.3)
    if len(hf):
        ax.hist2d(hf.x, hf.y, bins=bins, cmap="Reds", cmin=1, vmin=-5, vmax=max(20, np.percentile(np.histogram2d(hf.x, hf.y, bins=bins)[0], 99)) * 1.1, zorder=3)
    ax.axhline(0, color="#888", lw=0.6); ax.axvline(0, color="#888", lw=0.6)
    ax.annotate("", xy=(R * 0.9, 0), xytext=(R * 0.55, 0), arrowprops=dict(arrowstyle="->", color="k", lw=1.6)); ax.text(R * 0.9, 12, "motion", fontsize=7.5, ha="right")
    ax.scatter([0], [0], marker="x", color="k", s=70, zorder=6)
    ax.set_xlim(-R, R); ax.set_ylim(-R, R); ax.set_aspect(1)
    ax.set_xlabel("nm ahead (+) or behind (-) the centre"); ax.set_ylabel("nm to the right (+) or left (-) of motion")
    ax.set_title("Where the strong gusts were, relative to the low\n(orange: 64 kt+; red: 71.7 kt+, darker = more 3-hourly samples; ocean only)", fontsize=10)
    # (c) sector shares and radius
    ax = fig.add_subplot(gs[2]); ax.axis("off")
    ax.set_title("Share of the 71.7 kt+ ocean cells\nby quarter of the storm", fontsize=10)
    n = len(hf); shares = []
    for name, sx, sy in QUAD:
        shares.append(((np.sign(hf.x) == sx) & (np.sign(hf.y) == sy)).mean() * 100 if n else np.nan)
    for k, ((name, sx, sy), s) in enumerate(zip(QUAD, shares)):
        y0 = 0.80 - 0.17 * k
        ax.text(0.0, y0 + 0.04, name, transform=ax.transAxes, fontsize=9)
        ax.add_patch(plt.Rectangle((0.38, y0 + 0.01), 0.55 * (0 if not np.isfinite(s) else s / 100), 0.06, transform=ax.transAxes, color="#d0452f"))
        ax.text(0.38 + 0.55 * (0 if not np.isfinite(s) else s / 100) + 0.02, y0 + 0.02, "n/a" if not np.isfinite(s) else f"{s:.0f}%", transform=ax.transAxes, fontsize=9)
    if n:
        rr = np.hypot(hf.x, hf.y)
        ax.text(0.0, 0.08, f"median distance of those cells from\nthe centre: {np.median(rr):.0f} nm (90% within {np.percentile(rr, 90):.0f} nm)\n{n:,} cell-times at 71.7 kt+", transform=ax.transAxes, fontsize=8.5, va="top")
    fig.suptitle(title, fontsize=11, x=0.01, ha="left", y=1.02)
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
        cf = ax.contourf(lons, lats, t850, levels=np.arange(-30, 21, 2), cmap="RdYlBu_r", extend="both", zorder=1)
        ax.contour(lons, lats, t850, levels=[0], colors="#444", linewidths=0.8, linestyles="--", zorder=2)
        ax.contourf(lons, lats, high.astype(float), levels=[0.5, 1.5], colors=["#c9c8c2"], zorder=1.5)
        ax.contour(lons, lats, z[f"msl_{h}"], levels=np.arange(900, 1060, 4), colors="k", linewidths=0.8, zorder=3)
        if np.nanmax(ws) >= 50:
            ax.contour(lons, lats, ws, levels=[50, 64], colors=["#6a1b9a", "#6a1b9a"], linewidths=[1.2, 2.4], zorder=4)
        th = t0 + pd.Timedelta(hours=h); tracked = fx.t.iloc[0] <= th <= fx.t.iloc[-1]
        ax.plot([lo], [la], marker="x" if tracked else "+", color="#111" if tracked else "#777", ms=9, mew=2, zorder=5)
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
        s += f"At its lowest pressure it was <b>{PHASES[pk.phase][0]}</b> (B {pk.B:.0f} m, {pk.VTL:+.0f} m). "
    if ws:
        s += f"The first fix with a symmetric warm core (the warm-seclusion signature) is {hhmm(ws[0][1])}, {(ws[0][1] - t0) / pd.Timedelta(hours=1):+.0f} h from the lowest pressure. "
    else:
        s += "It never reached the symmetric warm-core signature of a warm seclusion in these 1.5&deg; diagnostics. "
    out["phase"] = s
    # forcing
    s = (f"During its fastest 24 h fall ({dt24(win)}) the highest 250 hPa wind within 1,000 km averaged <b>{win.jet250.mean():.0f} kt</b>, "
         f"300 hPa divergence over the centre averaged {win.div300.mean():+.1f} (units 1e-5 /s; positive = air being drawn up) and "
         f"500 hPa vorticity advection {win.vadv500.mean():+.1f} (1e-9 /s&sup2;; positive = cyclonic spin arriving aloft). "
         f"The 850-500 hPa Eady growth rate, a measure of baroclinic (temperature-contrast) support, averaged {win.eady.mean():.2f} per day, "
         f"and column water vapour within 500 km averaged {win.tcwv.mean():.0f} kg/m&sup2;.")
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
        s = (f"OPC's archive carries this low as event {opc['id']}: {opc['n_fix']} fixes, {opc['n_hf']} of them at hurricane force"
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
    return f"above {p:.0f}% of them"


def page(f, fx, tx, opc, meteo_, shares, style, base):
    esc = html.escape
    cls = "atl" if f["basin"] == "atl" else "pac"
    b, lab = old.label(f)
    pre = pd.Timestamp(f["t_minp"]).year < 2004
    dur = ((pd.Timestamp(f["last_hf"]) - pd.Timestamp(f["first_hf"])) / pd.Timedelta(hours=1) + 6) if f["first_hf"] else 0
    nm = f' &mdash; the {esc(f["name"])} <span class="muted">(name inferred from date and position)</span>' if f["name"] else ""
    rows = "".join(f'<tr><td>{esc(l)}</td><td class="num">{fnum(v, nd)} {esc(u)}</td><td class="num">{"n/a" if r is None else ordinal(r)}</td></tr>' for l, u, v, nd, r in meteo_["rank_rows"])
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
<div class="note"><strong>ERA5 proxy, not direct observation.</strong> Tracks are pipeline A (found in ERA5 itself). The gust index is the highest ERA5 10 m gust over ocean within 800 km, a proxy for the sustained wind OPC warns on. Structure and forcing numbers are ERA5 at 1.5&deg; (WeatherBench2); maps and the wind field are ERA5 at 0.25&deg;. Distances are nautical miles.</div>
<section id="storm" class="storm {cls}" style="margin-top:0">
  <h3><span class="rank">{b.replace("North ", "")} #{f["rank"]}</span><span class="basin-{cls}">{esc(lab)}</span>{nm}</h3>
  <div class="meta">Rank {f["rank"]} by lowest ERA5 central pressure in the {b}. Pipeline A track {f["track"]}.</div>
  {tiles}
  {pre_note}
  <div class="sum"><div class="tag">1. The storm in plain words</div><p>{old.summary(f)}</p></div>
</section>
<section><div class="tag">2. Track and wind field</div>
<figure class="fig"><div class="scroll"><a href="{base}-track.png"><img src="{base}-track.png" alt="Track and gust swath" loading="lazy"></a>
<figcaption>Left: 6-hourly fixes coloured by central pressure, red rings where the gust index reached 71.7 kt. Right: highest ERA5 10 m gust at each point from 48 h before to 24 h after the lowest pressure (land included; black contour 71.7 kt).</figcaption></div></figure></section>
<section><div class="tag">3. How it developed</div>
<p>{meteo_["phase"]}</p>
<figure class="fig"><div class="scroll"><a href="{base}-series.png"><img src="{base}-series.png" alt="Six-panel development series" loading="lazy"></a></div>
<figcaption>From the top: central pressure; the fall over the previous 24 h against the latitude-scaled bomb yardstick; the gust index with the ocean area at 71.7 kt or more; the 250 hPa jet; the two upper-level lift terms; and Hart&rsquo;s structure parameters with the structure class along the bottom. The tan band is the 24 h of fastest fall.</figcaption></figure></section>
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
<figcaption>Left: the storm&rsquo;s path through Hart&rsquo;s cyclone phase space. Middle: ocean cells at 64 kt or more, in a frame that moves and turns with the low. Right: how the 71.7 kt+ cells split between the four quarters.</figcaption></figure></section>
<section><div class="tag">7. How its set-up compares with other hurricane-force-strength lows</div>
<p>Median over this storm&rsquo;s {meteo_["rank_n"]} fixes with a gust index of 71.7 kt or more, placed among all such fixes in the same basin since 2004&ndash;05 (pipeline A, ERA5 proxy; fixes from one storm are not independent, so read it as a rough rank, not a test).</p>
<div class="tablewrap"><table class="idx"><thead><tr><th>Ingredient</th><th class="num">This storm</th><th class="num">Compared with HF-strength fixes</th></tr></thead><tbody>{rows}</tbody></table></div></section>
<section>{opc_block}</section>
<section><div class="tag">Read with care</div>
<p class="limit">ERA5 proxy throughout; never validated for storms before the archive begins. The Hart terms use 900, 850, 700, 600, 500, 400 and 300 hPa heights on a 1.5&deg; grid, not the original 50 hPa spacing, so read the structure class as a guide, not a classification to the metre. Over Greenland and high ground ERA5&rsquo;s sea-level pressure is extrapolated, and gusts near coasts and sea ice are less reliable. The storm-relative frame uses the last direction of motion while the low is stalled. A storm is one case: nothing here is a trend or a forecast rule.</p></section>
</main></body></html>"""


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
    build(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:])
