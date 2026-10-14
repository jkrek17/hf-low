#!/usr/bin/env python3
"""Draw the charts and write docs/storm-reviews.html (the storm-review page).

usage: build_storm_reviews.py WORKDIR FIXES_CSV OUTDIR
  WORKDIR     output of storm_review_extract.py (one npz per track)
  FIXES_CSV   research/era5/hf_history/results/era5_hf_catalog_tracks.csv (pipeline A, every 6 h point)
  OUTDIR      docs/ (PNG files storm-review-<track>-*.png and storm-reviews.html are written here)
The storm list is tools/storm_review_storms.csv. Everything is an ERA5 proxy (pipeline A tracks); the numbers in the
page text come from the facts table written to tools/storm_review_facts.csv, not typed by hand.
"""
import os, sys, json, html
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
HF = 71.7
MAP_H = (-48, -36, -24, -12, 0, 12, 24)
DLAT, DLON = 22.0, 40.0
plt.rcParams.update({"font.size": 9.5, "axes.titlesize": 10.5, "axes.titleweight": "semibold", "figure.facecolor": "white", "savefig.facecolor": "white",
                     "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": "#8a8984", "axes.labelcolor": "#33322f",
                     "xtick.color": "#55544f", "ytick.color": "#55544f", "axes.grid": False, "legend.frameon": False})


def coast():
    with open(os.path.join(HERE, "..", "docs", "data", "coastlines.js")) as f:
        s = f.read()
    return json.loads(s[s.index("{"):s.rindex("}") + 1])["polygons"]


def draw_coast(ax, polys, lw=0.6):
    for p in polys:
        a = np.array(p)
        if np.abs(np.diff(a[:, 0])).max() > 180:
            continue                                     # skips the few polygons that straddle the antimeridian
        for sh in (-360, 0, 360):
            ax.plot(a[:, 0] + sh, a[:, 1], color="#555", lw=lw, zorder=3)


def crop_axes(la, lo, dlat, dlon):
    r0 = max(int(round((90 - (la + dlat)) / 0.25)), 0); r1 = min(int(round((90 - (la - dlat)) / 0.25)) + 1, 721)
    c0 = int(round((lo - dlon) / 0.25)); c1 = int(round((lo + dlon) / 0.25)) + 1
    return 90 - 0.25 * np.arange(r0, r1), 0.25 * np.arange(c0, c1)


def lonfmt(ax):
    def f(x, _):
        v = ((x + 180) % 360) - 180
        return f"{abs(v):.0f}\u00b0{'E' if 0 <= v < 180 else 'W'}"
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(f))


def prep(fx):
    fx = fx.sort_values("t").reset_index(drop=True)
    fx["lonu"] = np.degrees(np.unwrap(np.radians(fx.lon.values)))
    fx.loc[fx.lonu > 180, "lonu"] = fx.lonu          # keep unwrapped, shifted below if needed
    fx["gust"] = fx.g800.where(fx.g800 > 0)
    return fx


def deepening(fx):
    """Change over the previous 24 h in hPa per 24 h (negative = deepening) at each fix that has a fix 24 h earlier."""
    s = fx.set_index("t").msl
    prev = s.reindex(s.index - pd.Timedelta(hours=24)).values
    return (s.values - prev)


def facts(row, fx, d24, t0, sw):
    i0 = fx.msl.idxmin()
    ng = fx.gust.idxmax()
    hf = fx[fx.gust >= HF]
    j = int(np.nanargmin(d24))
    lat_mean = float(fx.lat.iloc[max(j - 4, 0):j + 1].mean())
    bomb = 24 * np.sin(np.radians(lat_mean)) / np.sin(np.radians(60))
    out = dict(track=row.track, basin=row.basin, name=row.name_inferred if isinstance(row.name_inferred, str) else "",
               rank=row["rank"], gust_flag=row.gust_flag if isinstance(row.gust_flag, str) else "",
               t_minp=str(fx.t[i0]), minp=float(fx.msl[i0]), lat_minp=float(fx.lat[i0]), lon_minp=float(((fx.lon[i0] + 180) % 360) - 180),
               t_gust=str(fx.t[ng]), gust_idx=float(fx.gust[ng]),
               first_hf=str(hf.t.iloc[0]) if len(hf) else "", last_hf=str(hf.t.iloc[-1]) if len(hf) else "",
               n_hf=int(len(hf)), hf_lat=float(hf.lat.iloc[0]) if len(hf) else np.nan,
               hf_lon=float(((hf.lon.iloc[0] + 180) % 360) - 180) if len(hf) else np.nan,
               first_fix=str(fx.t.iloc[0]), last_fix=str(fx.t.iloc[-1]), p_first=float(fx.msl.iloc[0]),
               dp24_max=float(-d24[j]), t_dp24=str(fx.t[j]), lat_dp24=lat_mean, bomb_thr=float(bomb),
               sw_max=float(np.nanmax(sw)))
    return out


def fig_track(fx, sw, slat, slon, polys, title, path):
    lon0 = fx.lonu
    x0, x1 = lon0.min() - 12, lon0.max() + 12
    y0, y1 = max(fx.lat.min() - 8, 5), min(fx.lat.max() + 8, 85)
    mid = 0.5 * (y0 + y1)
    hgt = float(np.clip(5.2 * (y1 - y0) / ((x1 - x0) * np.cos(np.radians(mid))) + 1.5, 3.8, 6.6))   # panel height follows the map's shape
    fig, axs = plt.subplots(1, 2, figsize=(13.5, hgt), gridspec_kw=dict(wspace=0.16))
    for ax in axs:
        ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.set_aspect(1 / np.cos(np.radians(mid)))
        draw_coast(ax, polys); ax.grid(color="#bbb", lw=0.4, ls=":"); ax.set_xlabel("longitude"); lonfmt(ax); ax.set_ylabel("latitude (N)")
    ax = axs[0]
    ax.plot(lon0, fx.lat, color="#888", lw=1, zorder=4)
    sc = ax.scatter(lon0, fx.lat, c=fx.msl, cmap="plasma_r", s=34, vmin=910, vmax=1010, zorder=5, edgecolor="k", linewidth=0.3)
    hf = fx[fx.gust >= HF]
    ax.scatter(hf.lonu, hf.lat, s=130, facecolor="none", edgecolor="#d00", linewidth=1.6, zorder=6)
    i0 = fx.msl.idxmin()
    ax.scatter([fx.lonu[i0]], [fx.lat[i0]], marker="*", s=260, color="#fff", edgecolor="k", zorder=7)
    for k, r in fx.iterrows():
        if r.t.hour == 0:
            ax.annotate(r.t.strftime("%d %b"), (r.lonu, r.lat), xytext=(4, 5), textcoords="offset points", fontsize=7, zorder=8)
    cb = fig.colorbar(sc, ax=ax, shrink=0.8, pad=0.02); cb.set_label("central pressure (hPa)")
    ax.legend(handles=[Line2D([], [], marker="o", ls="", markerfacecolor="none", markeredgecolor="#d00", markersize=10, label="HF-strength fix (gust index ≥ 71.7 kt)"),
                       Line2D([], [], marker="*", ls="", markerfacecolor="w", markeredgecolor="k", markersize=14, label="lowest pressure")],
              loc="lower left", fontsize=7.5, framealpha=0.9)
    ax.set_title("Track, coloured by central pressure (6-hourly fixes)")
    ax = axs[1]
    lonc = slon.copy()
    lev = [40, 50, 60, 64, 72, 80, 90, 100, 110]
    cf = ax.contourf(lonc, slat, sw, levels=lev, cmap="YlOrRd", extend="max", zorder=1)
    ax.contour(lonc, slat, sw, levels=[HF], colors="k", linewidths=1.4, zorder=2)
    ax.plot(lon0, fx.lat, color="k", lw=1, zorder=4); ax.scatter([fx.lonu[i0]], [fx.lat[i0]], marker="*", s=200, color="#fff", edgecolor="k", zorder=7)
    cb = fig.colorbar(cf, ax=ax, shrink=0.8, pad=0.02); cb.set_label("highest ERA5 10 m gust (kt)")
    ax.set_title("Gust swath, 48 h before to 24 h after lowest pressure\n(black contour = 71.7 kt)")
    fig.suptitle(title, fontsize=12, fontweight="semibold", x=0.01, ha="left", y=1.0)
    fig.savefig(path, dpi=100, bbox_inches="tight"); plt.close(fig)


def fig_series(fx, d24, t0, title, path):
    hrs = (fx.t - t0) / pd.Timedelta(hours=1)
    fig, axs = plt.subplots(3, 1, figsize=(9, 7.4), sharex=True, gridspec_kw=dict(hspace=0.12))
    axs[0].plot(hrs, fx.msl, "-o", ms=3, color="#0a3f73"); axs[0].set_ylabel("central pressure (hPa)"); axs[0].invert_yaxis()
    axs[0].axvline(0, color="k", lw=0.8, ls=":")
    dd = -d24
    axs[1].bar(hrs, dd, width=4.5, color=np.where(dd > 0, "#c0392b", "#8aa"))
    thr = 24 * np.sin(np.radians(fx.lat)) / np.sin(np.radians(60))
    axs[1].plot(hrs, thr, color="k", lw=1, ls="--", label="'bomb' threshold at the storm's latitude\n(24 hPa/24 h at 60°N, scaled by sin latitude)")
    axs[1].set_ylabel("deepening over the previous 24 h\n(hPa per 24 h)"); axs[1].legend(fontsize=7.5, loc="upper left"); axs[1].axvline(0, color="k", lw=0.8, ls=":")
    axs[2].plot(hrs, fx.gust, "-o", ms=3, color="#b9460f"); axs[2].axhline(HF, color="k", lw=1, ls="--")
    axs[2].text(hrs.min(), HF + 1.5, "71.7 kt: hurricane-force-equivalent", fontsize=7.5)
    axs[2].set_ylabel("gust index (kt)\nmax ERA5 gust within 800 km"); axs[2].axvline(0, color="k", lw=0.8, ls=":")
    axs[2].set_xlabel(f"hours from lowest pressure ({t0:%d %b %Y %H}Z)")
    for a in axs: a.grid(color="#ccc", lw=0.5)
    fig.suptitle(title, fontsize=11, x=0.01, ha="left", y=0.93)
    fig.savefig(path, dpi=100, bbox_inches="tight"); plt.close(fig)


def fig_maps(z, polys, title, path, fx, t0):
    fig, axs = plt.subplots(2, 4, figsize=(17, 8.6), gridspec_kw=dict(wspace=0.08, hspace=0.22))
    axs = axs.ravel()
    cf = None
    for k, h in enumerate(MAP_H):
        ax = axs[k]; la, lo = z[f"centre_{h}"]
        lats, lons = crop_axes(la, lo, DLAT, DLON)
        msl, z5, ws = z[f"msl_{h}"], z[f"z500_{h}"], z[f"ws250_{h}"]
        cf = ax.contourf(lons, lats, ws, levels=np.arange(30, 91, 10), cmap="BuPu", extend="max", zorder=1)
        cs = ax.contour(lons, lats, z5, levels=np.arange(4800, 6200, 60), colors="#1f6fbf", linewidths=0.7, zorder=2)
        ax.contour(lons, lats, msl, levels=np.arange(900, 1060, 4), colors="k", linewidths=0.8, zorder=3)
        ax.contour(lons, lats, msl, levels=[992], colors="k", linewidths=1.6, zorder=3)
        th = t0 + pd.Timedelta(hours=h)
        tracked = fx.t.iloc[0] <= th <= fx.t.iloc[-1]
        ax.plot([lo], [la], marker="x" if tracked else "+", color="#d00" if tracked else "#777", ms=9, mew=2, zorder=5)
        if not tracked:
            ax.text(0.02, 0.03, "not yet tracked: marker at the\nnearest tracked position", transform=ax.transAxes, fontsize=7, zorder=6, bbox=dict(fc="w", ec="none", alpha=0.8))
        ax.set_xlim(lons[0], lons[-1]); ax.set_ylim(lats[-1], lats[0]); ax.set_aspect(1 / np.cos(np.radians(la)))
        draw_coast(ax, polys, 0.7)
        mn = np.nanmin(msl)
        ax.set_title(f"{'peak' if h == 0 else f'{h:+d} h'}: lowest {mn:.0f} hPa", fontsize=9.5)
        ax.tick_params(labelsize=7); lonfmt(ax)
    ax = axs[7]; ax.axis("off")
    cb = fig.colorbar(cf, ax=ax, location="top", shrink=0.9, pad=0.0); cb.set_label("250 hPa wind speed (m/s)")
    ax.text(0, 0.55, "black: sea-level pressure every 4 hPa\n(thick = 992 hPa)\nblue: 500 hPa height every 60 m\nred ×: storm position (pipeline A track;\ngrey + where the tracker has not yet\npicked it up or has lost it)\n\nERA5 proxy at 0.25°, 6-hourly\nanalysis fields; the map moves\nwith the storm.", fontsize=9, va="top", transform=ax.transAxes)
    fig.suptitle(title, fontsize=11, x=0.01, ha="left", y=0.955)
    fig.savefig(path, dpi=85, bbox_inches="tight"); plt.close(fig)


MONTH = "January February March April May June July August September October November December".split()


def dt(s, fmt="%d %b %Y %H00Z"):
    return pd.Timestamp(s).strftime(fmt)


def label(f):
    t = pd.Timestamp(f["t_minp"])
    b = "North Atlantic" if f["basin"] == "atl" else "North Pacific"
    return b, f"{b}, {t.day} {MONTH[t.month - 1]} {t.year}"


def summary(f):
    t0 = pd.Timestamp(f["t_minp"]); b, _ = label(f)
    fh = pd.Timestamp(f["first_hf"]) if f["first_hf"] else None
    s = f"ERA5 puts its lowest pressure at {f['minp']:.1f} hPa on {dt(f['t_minp'])}, near {f['lat_minp']:.1f}&deg;N {abs(f['lon_minp']):.1f}&deg;{'W' if f['lon_minp'] < 0 else 'E'}. "
    if f['dp24_max'] > 0:
        s += f"Its fastest 24 h fall was {f['dp24_max']:.0f} hPa, ending {dt(f['t_dp24'])} at about {f['lat_dp24']:.0f}&deg;N"
        s += (f" (the bomb threshold there is {f['bomb_thr']:.0f} hPa, so it qualifies as a bomb by that yardstick). " if f['dp24_max'] >= f['bomb_thr']
              else f" (the bomb threshold there is {f['bomb_thr']:.0f} hPa, so it falls short of that yardstick). ")
    else:
        s += "The tracker follows it only after most of its deepening, so no 24 h fall can be measured from the track and the bomb test cannot be applied. "
    if fh is not None:
        dur = (pd.Timestamp(f["last_hf"]) - fh) / pd.Timedelta(hours=1) + 6
        s += f"It first reached hurricane-force-equivalent gusts (index at least 71.7 kt) on {dt(f['first_hf'])} near {f['hf_lat']:.1f}&deg;N {abs(f['hf_lon']):.1f}&deg;{'W' if f['hf_lon'] < 0 else 'E'}, " \
             f"was at that strength in {f['n_hf']} of its 6-hourly fixes (about {dur:.0f} h from first to last), and peaked at a gust index of {f['gust_idx']:.0f} kt on {dt(f['t_gust'])}. "
    pre_h = (t0 - pd.Timestamp(f["first_fix"])) / pd.Timedelta(hours=1)
    if pre_h < 48:
        s += f"The tracker only picks it up {pre_h:.0f} h before the lowest pressure, so its early development, and the 24 h fall above, are only partly captured. "
    s += f"The gust swath maximum anywhere in the window, land included, is {f['sw_max']:.0f} kt."
    return s


def main(work, fixes_csv, outdir):
    polys = coast()
    st = pd.read_csv(os.path.join(HERE, "storm_review_storms.csv"))
    fx_all = pd.read_csv(fixes_csv)
    fx_all["t"] = pd.to_datetime(fx_all.time.astype(str), format="%Y%m%d%H")
    rows = []
    for r in st.itertuples(index=False):
        row = pd.Series(r._asdict()); row["rank"] = row["rank"]
        fx = prep(fx_all[fx_all.track == row.track].copy())
        z = np.load(os.path.join(work, f"{row.track}.npz"))
        t0 = fx.t[fx.msl.idxmin()]
        d24 = deepening(fx)
        f = facts(row, fx, d24, t0, z["sw_gust"])
        b, lab = label(f)
        head = f"{lab} ({'ERA5 proxy, pipeline A'})"
        base = os.path.join(outdir, f"storm-review-{row.track}")
        # lon axis for the track map: make the swath longitudes continuous with the track's unwrapped longitudes
        slon = z["sw_lons"].copy()
        slon = np.unwrap(np.radians(slon)); slon = np.degrees(slon)
        shift = 360 * np.round((fx.lonu.mean() - slon.mean()) / 360); slon = slon + shift
        fig_track(fx, z["sw_gust"], z["sw_lats"], slon, polys, head, base + "-track.png")
        fig_series(fx, d24, t0, head, base + "-series.png")
        fig_maps(z, polys, head + ": every 12 hours from 48 h before to 24 h after the lowest pressure", base + "-maps.png", fx, t0)
        rows.append(f)
        print(row.track, lab, flush=True)
    pd.DataFrame(rows).to_csv(os.path.join(HERE, "storm_review_facts.csv"), index=False)
    from PIL import Image                                  # 160-colour PNGs: 10.8 MB -> 3.1 MB for the 30 figures
    for r in rows:
        for k in ("track", "series", "maps"):
            fp = os.path.join(outdir, f"storm-review-{r['track']}-{k}.png")
            Image.open(fp).convert("RGB").quantize(colors=160, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE).save(fp, optimize=True)
    write_page(rows, outdir)


def fall(f):
    return f"{f['dp24_max']:.0f}" if f["dp24_max"] > 0 else "n/a"


def page_style():
    with open(os.path.join(HERE, "..", "docs", "composites.html")) as f:
        c = f.read()
    style = c[c.index("<style>"):c.index("</style>") + 8]
    style = style.replace("</style>", """  .idx td.num, .idx th.num { text-align: right; }
  .hero { display: grid; gap: 10px; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); }
  .hero .stat { background: var(--surface); border: 1px solid var(--border); border-radius: 8px; padding: 10px 12px; }
  .hero .stat b { display: block; font-size: 24px; line-height: 1.15; font-variant-numeric: tabular-nums; }
  .hero .stat span { font-size: 12px; color: var(--ink-2); }
  .stack { display: flex; flex-direction: column; gap: 22px; }
  meter.bar { appearance: none; -webkit-appearance: none; width: 64px; height: 8px; vertical-align: middle; margin-right: 6px; border: 0; border-radius: 4px; background: var(--surface-sunk); }
  meter.bar::-webkit-meter-bar { background: var(--surface-sunk); border: 0; border-radius: 4px; height: 8px; }
  meter.bar.atl::-webkit-meter-optimum-value, meter.bar.atl::-webkit-meter-suboptimum-value, meter.bar.atl::-webkit-meter-even-less-good-value { background: var(--atl); border-radius: 4px; }
  meter.bar.pac::-webkit-meter-optimum-value, meter.bar.pac::-webkit-meter-suboptimum-value, meter.bar.pac::-webkit-meter-even-less-good-value { background: var(--pac); border-radius: 4px; }
  meter.bar.atl::-moz-meter-bar { background: var(--atl); } meter.bar.pac::-moz-meter-bar { background: var(--pac); }
  @media (max-width: 700px) { .toolbar { position: static !important; } }
  .toolbar { position: sticky; top: 0; z-index: 5; background: var(--page); border-bottom: 1px solid var(--border); padding: 8px 0; display: flex; flex-wrap: wrap; gap: 8px 14px; align-items: center; font-size: 14px; }
  .toolbar button, .toolbar select { font: inherit; padding: 5px 12px; border: 1px solid var(--border-strong); border-radius: 999px; background: var(--surface); color: var(--ink); cursor: pointer; }
  .toolbar button[aria-pressed="true"] { background: var(--navy); color: var(--navy-ink); border-color: var(--navy); }
  .toolbar button:focus-visible, .toolbar select:focus-visible, summary:focus-visible { outline: 2px solid var(--accent); outline-offset: 2px; }
  .storm { background: var(--surface); border: 1px solid var(--border); border-radius: 10px; padding: 16px 18px; border-top: 4px solid var(--border-strong); scroll-margin-top: 56px; }
  .storm.atl { border-top-color: var(--atl); } .storm.pac { border-top-color: var(--pac); }
  .storm[hidden] { display: none; }
  .storm h3 { font-size: 19px; margin: 0 0 2px; }
  .storm .rank { display: inline-block; font-size: 12px; font-weight: 650; padding: 1px 8px; border-radius: 999px; background: var(--surface-sunk); border: 1px solid var(--border); color: var(--ink-2); margin-right: 8px; vertical-align: 2px; }
  .storm .meta { font-size: 13px; color: var(--ink-2); margin-bottom: 10px; }
  .tiles { display: grid; gap: 8px; grid-template-columns: repeat(auto-fit, minmax(120px, 1fr)); margin: 0 0 12px; }
  .tiles div { background: var(--surface-sunk); border-radius: 8px; padding: 8px 10px; }
  .tiles b { display: block; font-size: 20px; line-height: 1.2; font-variant-numeric: tabular-nums; }
  .tiles span { font-size: 12px; color: var(--ink-2); }
  .storm .sum { border-left: 3px solid var(--accent); padding: 2px 0 2px 12px; margin-bottom: 12px; }
  .storm .sum p { margin: 0 0 6px; }
  .storm details { border-top: 1px solid var(--border); padding: 8px 0 2px; }
  .storm summary { cursor: pointer; font-weight: 650; padding: 4px 0; }
  .storm .fig { margin: 8px 0 10px; }
  .storm .fig img { border-radius: 6px; }
  .count { color: var(--ink-2); font-size: 13px; }
  .blk { background: var(--surface); border: 1px solid var(--border); border-radius: 8px; padding: 12px 14px; }
  section .tag { margin-bottom: 4px; }
</style>""")
    return style


def write_page(rows, outdir):
    style = page_style()
    esc = html.escape
    secs = []
    tr = []
    srt = sorted(rows, key=lambda f: (f["basin"], f["minp"]))
    gmax = max(f["gust_idx"] for f in rows)
    for f in srt:
        b, lab = label(f)
        cls = "atl" if f["basin"] == "atl" else "pac"
        pre = pd.Timestamp(f["t_minp"]).year < 2004
        anchor = f"storm-{f['track']}"
        dur = ((pd.Timestamp(f["last_hf"]) - pd.Timestamp(f["first_hf"])) / pd.Timedelta(hours=1) + 6) if f["first_hf"] else 0
        tr.append(f'<tr><td><a href="#{anchor}"><span class="basin-{cls}">{esc(lab)}</span></a></td><td class="num">{f["minp"]:.1f}</td>'
                  f'<td class="num">{fall(f)}</td><td class="num"><meter class="bar {cls}" min="0" max="{gmax:.0f}" value="{f["gust_idx"]:.0f}" aria-hidden="true"></meter> {f["gust_idx"]:.0f}</td>'
                  f'<td class="num">{dur:.0f}</td><td>{"reanalysis only" if pre else "archive era"}</td></tr>')
        nm = f' &mdash; the {esc(f["name"])} <span class="muted">(name inferred from date and position)</span>' if f["name"] else ""
        note = ""
        if pre:
            note = ('<p class="limit"><strong>Not direct observation.</strong> This storm is before the 2004&ndash;05 season. '
                    'No archive record covers it, so every number and map here is ERA5 reanalysis. '
                    + ('The gust index also runs high before 2001 (ERA5 gust drifts upward at fixed storm depth), so read its level with caution; central pressure is the better yardstick. ' if f["gust_flag"] else '')
                    + 'No trend is drawn from it.</p>')
        tiles = (f'<div class="tiles"><div><b>{f["minp"]:.1f}</b><span>lowest pressure (hPa)</span></div>'
                 f'<div><b>{fall(f)}</b><span>fastest 24 h fall (hPa)</span></div>'
                 f'<div><b>{f["gust_idx"]:.0f}</b><span>peak gust index (kt)</span></div>'
                 f'<div><b>{dur:.0f}</b><span>hours at HF strength</span></div></div>')
        sec = f"""<section id="{anchor}" class="storm {cls}" data-basin="{cls}" data-minp="{f["minp"]:.1f}" data-fall="{max(f["dp24_max"], 0):.1f}" data-gust="{f["gust_idx"]:.1f}" data-year="{pd.Timestamp(f["t_minp"]).year}">
  <h3><span class="rank">{b.replace("North ", "")} #{f["rank"]}</span><span class="basin-{cls}">{esc(lab)}</span>{nm}</h3>
  <div class="meta">Rank {f["rank"]} by lowest ERA5 central pressure in the {b}. ERA5 proxy, pipeline A track {f["track"]}.</div>
  {tiles}
  <div class="sum"><div class="tag">In plain words</div><p>{summary(f)}</p></div>
  {note}
<figure class="fig">
  <div class="scroll"><a href="storm-review-{f["track"]}-track.png"><img src="storm-review-{f["track"]}-track.png" alt="Track coloured by central pressure with hurricane-force fixes ringed, and a map of the highest ERA5 gust for {esc(lab)}" loading="lazy"></a></div>
  <figcaption>Where it went, and where the wind was. Left: 6-hourly pipeline A fixes coloured by central pressure; red rings mark fixes with a gust index of at least 71.7 kt. Right: the highest ERA5 10 m gust at each point from 48 h before to 24 h after the lowest pressure. It includes land and the strongest gusts around hills and coasts, so it is not the same as the ocean-only gust index.</figcaption>
</figure>
<details>
  <summary>How it developed</summary>
  <figure class="fig">
    <div class="scroll"><a href="storm-review-{f["track"]}-series.png"><img src="storm-review-{f["track"]}-series.png" alt="Time series of central pressure, 24-hour deepening and gust index for {esc(lab)}" loading="lazy"></a></div>
    <figcaption>Central pressure, the fall over the previous 24 h against the latitude-scaled bomb threshold, and the gust index against the 71.7 kt line.</figcaption>
  </figure>
</details>
<details>
  <summary>Weather maps every 12 hours</summary>
  <figure class="fig">
    <div class="scroll"><a href="storm-review-{f["track"]}-maps.png"><img src="storm-review-{f["track"]}-maps.png" alt="Maps of sea-level pressure, 500 hPa height and 250 hPa wind speed every 12 hours around {esc(lab)}" loading="lazy"></a></div>
    <figcaption>The pattern around it from 48 h before to 24 h after the lowest pressure: sea-level pressure (black), 500 hPa height (blue) and 250 hPa wind speed (shaded jet). Over Greenland and other high ground ERA5's sea-level pressure is extrapolated and noisy.</figcaption>
  </figure>
</details>
</section>"""
        secs.append(sec)
    n = len(rows); na = sum(f["basin"] == "atl" for f in rows); npc = n - na
    npre = sum(pd.Timestamp(f["t_minp"]).year < 2004 for f in rows); npost = n - npre
    deep = min(rows, key=lambda f: f["minp"])
    fast = max(rows, key=lambda f: f["dp24_max"])
    nhf = sum(1 for f in rows if f["first_hf"])
    page = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Storm reviews - Hurricane Force Extratropical Low Archive</title>
<link rel="icon" type="image/svg+xml" href="favicon.svg">
<meta name="description" content="Storm reviews of the {n} deepest extratropical lows in the ERA5 proxy record, {na} Atlantic and {npc} Pacific: track, development, weather maps every 12 hours and gust swath.">
<!-- generated by tools/build_storm_reviews.py from tools/storm_review_facts.csv: edit the script, not this file -->
{style}
</head>
<body>
<header class="masthead"><div class="masthead-in">
  <h1>Storm reviews: the strongest lows, one at a time</h1>
  <a href="index.html">&larr; Back to the archive</a>
</div></header>
<main>
<div class="note"><strong>ERA5 proxy, not direct observation.</strong> Every storm here was found in ERA5 reanalysis by the project's pipeline A tracker, and every chart is drawn from ERA5 fields (0.25&deg;, 6-hourly). {npre} of the {n} are before the 2004&ndash;05 season, where no archive record exists; the other {npost} are since then, and are still read here from ERA5. ERA5 smooths the very strongest winds and pressure minima, so real peaks may be stronger. These are case studies, not a trend: nothing here says whether such storms are becoming more or less common.</div>

<div class="hero">
  <div class="stat"><b>{n}</b><span>storms: {na} North Atlantic, {npc} North Pacific</span></div>
  <div class="stat"><b>{deep["minp"]:.1f} hPa</b><span>lowest: {esc(label(deep)[1])}</span></div>
  <div class="stat"><b>{fast["dp24_max"]:.0f} hPa</b><span>fastest 24 h fall: {esc(label(fast)[1])}</span></div>
  <div class="stat"><b>{nhf} of {n}</b><span>reach a gust index of 71.7 kt or more</span></div>
</div>

<div class="answers">
  <div class="tag">At a glance</div>
  <div class="tablewrap"><table class="idx">
  <thead><tr><th>Storm</th><th class="num">Lowest pressure (hPa)</th><th class="num">Fastest 24 h fall (hPa)</th><th class="num">Peak gust index (kt)</th><th class="num">Hours at HF strength</th><th>Source</th></tr></thead>
  <tbody>{"".join(tr)}</tbody></table></div>
  <p class="muted">Hours at HF strength is the number of 6-hourly fixes with a gust index of at least 71.7 kt times six, counted from the first to the last such fix. Peak gust index is the highest ERA5 gust within 800 km over the ocean, one number per fix. Click a storm name to jump to its review.</p>
</div>

<div class="toolbar" role="group" aria-label="Show storms">
  <span>Show:</span>
  <button type="button" data-f="all" aria-pressed="true">All {n}</button>
  <button type="button" data-f="atl" aria-pressed="false">North Atlantic ({na})</button>
  <button type="button" data-f="pac" aria-pressed="false">North Pacific ({npc})</button>
  <label>Sort by <select id="sortby">
    <option value="minp">Lowest pressure</option><option value="fall">Fastest 24 h fall</option><option value="gust">Peak gust index</option><option value="year">Newest first</option>
  </select></label>
  <span class="count" id="count" aria-live="polite"></span>
</div>

<div id="storms" class="stack">
{"".join(secs)}
</div>

<section id="how">
  <div class="tag">Method</div>
  <h2>How the storms were chosen and drawn</h2>
  <p>The storms are those with the lowest ERA5 central pressure in each basin among the 1979&ndash;2025 tracks of pipeline A (the research tracker, <code>research/era5/hf_history</code>), taken in order of pressure from the lowest. Tropical-cyclone-flagged tracks are left out (three Pacific tracks, among them 30 Oct 1989 and 3 Oct 1981, which would otherwise rank second and third), and one Pacific storm from the 2024&ndash;25 season is left out because the weather-map fields used here end in January 2023. The rank shown is the rank among all tracks, so some numbers are missing. Ranking is by pressure alone, so the list is not a list of the strongest winds.</p>
  <p>Central pressure, position and gust index come from the 6-hourly pipeline A fixes. The gust index is the highest ERA5 10 m gust over ocean within 800 km of the low. The weather maps use WeatherBench2 ERA5 at 0.25&deg;, 6-hourly (sea-level pressure, 500 hPa geopotential height, 250 hPa wind speed), re-centred on the storm at each map time. The swath uses ERA5 hourly instantaneous gusts every 3 h over the 72 h window. The 24 h fall is the change in central pressure over the previous 24 h at each 6-hourly fix; the bomb line scales 24 hPa per 24 h at 60&deg;N by sin(latitude)/sin(60&deg;). OPC warns on sustained winds; the gust index is a proxy for that, not the same thing.</p>
  <p class="limit"><strong>Read with care.</strong> ERA5 proxy throughout, never validated for storms before the archive begins. Before 2001 the gust index drifts upward at fixed storm depth, so gust numbers for the 1980s and 1990s storms are less reliable than their pressures. Not independently re-checked: how closely these ERA5 pressures match the observed minima of the named storms.</p>
</section>

<footer>
  <p class="muted">Seasons run 1 June to 31 May. Fields: WeatherBench2 ERA5 and ARCO-ERA5 (Google Cloud public data). <a href="index.html">Back to the archive</a></p>
</footer>
</main>
<script>
  /* Follow the theme chosen on the archive page, if one was saved. */
  (function () {{
    try {{
      var saved = localStorage.getItem('hf-theme');
      if (saved) document.documentElement.setAttribute('data-theme', saved);
    }} catch (err) {{ /* private mode: follow the system theme */ }}
  }})();
  /* Basin filter and sort. Without script every review simply shows in rank order. */
  (function () {{
    var box = document.getElementById('storms');
    var cards = Array.prototype.slice.call(box.querySelectorAll('.storm'));
    var btns = document.querySelectorAll('.toolbar button');
    var sel = document.getElementById('sortby');
    var count = document.getElementById('count');
    var filt = 'all';
    function render() {{
      var key = sel.value;
      cards.sort(function (a, b) {{
        var x = parseFloat(a.dataset[key]), y = parseFloat(b.dataset[key]);
        return key === 'minp' ? x - y : y - x;
      }});
      var shown = 0;
      cards.forEach(function (c) {{
        var on = filt === 'all' || c.dataset.basin === filt;
        c.hidden = !on; if (on) shown++;
        box.appendChild(c);
      }});
      count.textContent = shown + ' shown';
    }}
    function openHash() {{
      var el = location.hash && document.getElementById(location.hash.slice(1));
      if (el && el.classList.contains('storm')) {{
        if (el.hidden) {{ filt = 'all'; btns.forEach(function (b) {{ b.setAttribute('aria-pressed', b.dataset.f === 'all'); }}); render(); }}
        el.scrollIntoView();
      }}
    }}
    btns.forEach(function (b) {{
      b.addEventListener('click', function () {{
        filt = b.dataset.f;
        btns.forEach(function (o) {{ o.setAttribute('aria-pressed', o === b ? 'true' : 'false'); }});
        render();
      }});
    }});
    sel.addEventListener('change', render);
    window.addEventListener('hashchange', openHash);
    render(); openHash();
  }})();
</script>
</body>
</html>
"""
    with open(os.path.join(outdir, "storm-reviews.html"), "w") as fh:
        fh.write("\n".join(l.rstrip() for l in page.splitlines()) + "\n")


if __name__ == "__main__":
    main(*sys.argv[1:4])
