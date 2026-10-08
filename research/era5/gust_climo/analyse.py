"""Gust climatology maps and within-era gust rankings. ERA5 PROXY, not observation.

usage: ERA5_WORK=<dir> python3 -I research/era5/gust_climo/analyse.py    (after extract.py)

Inputs: the partial aggregates written by extract.py (ERA5 instantaneous 10 m gust, 12-hourly 00/12 UTC, October-April,
seasons 2004-05 to 2025-26, ocean points only here), results/ of hf_structure (location of the maximum gust at HF-strength
fixes, 2004-05 on) and hf_history (pipeline A catalog, for the within-era rankings).

Decision 1: the gust climatology maps use 2004-05 on only. The rankings are WITHIN eras (1979-2000 and 2004-2025 separately),
because pre-2001 gust at fixed storm depth drifts upward; a 1979-2000 gust is not comparable with a 2004+ gust. Sampling is
12-hourly, so maxima and exceedance frequencies are for the sampled times, not the true maxima. Gust is ocean-only here and is
ERA5's instantaneous 10 m gust: it is not a sustained wind.
"""
import glob
import os
import sys
import urllib.request
import numpy as np
import pandas as pd
import numcodecs
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "climatology_atlas"))
import figures as F_  # noqa: E402  (draw_land, coast)
import atlas as A_  # noqa: E402

WORK = os.environ.get("ERA5_WORK", os.path.join(HERE, "work"))
RES = os.path.join(HERE, "results")
HFH = os.path.join(HERE, "..", "hf_history", "results")
THRESH = (34.0, 50.0, 64.0, 71.7, 80.0, 90.0)
ARCO = "https://storage.googleapis.com/gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3"
WARN = "ERA5 PROXY (reanalysis instantaneous gust), NOT DIRECT OBSERVATION. 2004-05 to 2025-26, October-April, 12-hourly samples."
LAT = 90 - 0.25 * np.arange(361)
LON = 0.25 * np.arange(1440)
ROW0, ROW1 = 40, 301       # 80N .. 15N


def banner(fig, extra="", warn=WARN):
    fig.get_layout_engine().set(rect=(0, 0.04, 1, 0.96))
    fig.text(0.5, 0.005, warn + (" " + extra if extra else ""), ha="center", va="bottom", fontsize=8, color="#a40000", weight="bold")


def lsm():
    c = numcodecs.Blosc()
    with urllib.request.urlopen(f"{ARCO}/land_sea_mask/876576.0.0", timeout=120) as r:
        raw = r.read()
    return np.frombuffer(c.decode(raw), "<f4").reshape(721, 1440)[:361]


def combine():
    files = sorted(glob.glob(os.path.join(WORK, "chunk_*.npz")))
    n = 0
    sm = 0
    mx = 0
    ex = 0
    smax = {}
    nbytes = 0
    for f in files:
        z = np.load(f)
        n += int(z["n"])
        sm = sm + z["sm"]
        mx = np.maximum(mx, z["mx"])
        ex = ex + z["ex"].astype(np.uint32)
        nbytes += int(z["nbytes"])
        for s, a in zip(z["seasons"], z["smax"]):
            smax[int(s)] = np.maximum(smax.get(int(s), 0), a.astype(np.float32))
    return files, n, sm, mx, ex, smax, nbytes


def sub(a):
    return a[..., ROW0:ROW1, :]


def region_plot(ax, basin, field, vmin, vmax, cmap, ocean):
    x0, x1, y0, y1 = F_.EXT[basin]
    lat = LAT[ROW0:ROW1]
    if basin == "atl":
        sel = (LON >= 260) | (LON <= 30)
        lon = np.where(LON[sel] > 180, LON[sel] - 360, LON[sel])
        order = np.argsort(lon)
        lon = lon[order]
        z = field[:, sel][:, order]
    else:
        sel = (LON >= 120) & (LON <= 250)
        lon = LON[sel]
        z = field[:, sel]
    m = ax.pcolormesh(lon, lat, np.ma.masked_invalid(z), cmap=cmap, vmin=vmin, vmax=vmax, shading="auto", zorder=2)
    F_.draw_land(ax, basin)
    return m


def main():
    os.makedirs(RES, exist_ok=True)
    files, n, sm, mx, ex, smax, nbytes = combine()
    ocean = sub((lsm() < 0.5)[None])[0]
    mean = sub(sm[None])[0] / n
    mxs = sub(mx[None])[0]
    exs = sub(ex)
    annmax = np.mean([sub(smax[s][None])[0] for s in sorted(smax)], axis=0)
    freq = exs / n
    nan = lambda a: np.where(ocean, a, np.nan)
    lat = LAT[ROW0:ROW1]
    out = [WARN, "", f"Chunks {len(files)}, sampled times {n}, seasons {len(smax)}, streamed {nbytes/1e9:.1f} GB (global fields, northern hemisphere kept)"]
    names = {"annmax": nan(annmax), "absmax": nan(mxs), "mean": nan(mean)}
    for k, th in enumerate(THRESH):
        names[f"freq_ge{th:g}"] = nan(freq[k])
    domain = {}
    for b in ("atl", "pac"):
        x0, x1, y0, y1 = F_.EXT[b]
        if b == "atl":
            lm = ((LON >= 260) | (LON <= 30))
        else:
            lm = (LON >= 120) & (LON <= 250)
        la = (lat >= 20) & (lat <= 75)
        domain[b] = (la[:, None] & lm[None, :])
    for b in ("atl", "pac"):
        d = domain[b]
        a = names["annmax"]
        v = np.where(d, a, np.nan)
        i = np.nanargmax(v)
        r, c = np.unravel_index(i, v.shape)
        lon = LON[c] if b == "pac" else (LON[c] - 360 if LON[c] > 180 else LON[c])
        out.append(f"{b}: highest mean annual maximum gust {np.nanmax(v):.1f} kt at {lat[r]:.2f}N, {lon:.2f}E; ocean median of mean annual max {np.nanmedian(v):.1f} kt")
        fr = np.where(d, names["freq_ge64"], np.nan)
        i = np.nanargmax(fr)
        r, c = np.unravel_index(i, fr.shape)
        lon = LON[c] if b == "pac" else (LON[c] - 360 if LON[c] > 180 else LON[c])
        out.append(f"{b}: highest share of sampled times with gust >= 64 kt {np.nanmax(fr)*100:.2f}% at {lat[r]:.2f}N, {lon:.2f}E")
        fr = np.where(d, names["freq_ge71.7"], np.nan)
        i = np.nanargmax(fr)
        r, c = np.unravel_index(i, fr.shape)
        lon = LON[c] if b == "pac" else (LON[c] - 360 if LON[c] > 180 else LON[c])
        out.append(f"{b}: highest share with gust >= 71.7 kt {np.nanmax(fr)*100:.2f}% at {lat[r]:.2f}N, {lon:.2f}E")
        am = np.where(d, names["absmax"], np.nan)
        i = np.nanargmax(am)
        r, c = np.unravel_index(i, am.shape)
        lon = LON[c] if b == "pac" else (LON[c] - 360 if LON[c] > 180 else LON[c])
        out.append(f"{b}: highest single sampled gust {np.nanmax(am):.1f} kt at {lat[r]:.2f}N, {lon:.2f}E")
        # share of ocean cells in the domain with mean annual max >= 64 / 71.7 kt
        v = np.where(d & ocean, annmax, np.nan)
        tot = np.isfinite(v).sum()
        out.append(f"{b}: share of ocean cells (20-75N) whose mean annual maximum gust is >= 64 kt {np.nansum(v>=64)/tot*100:.1f}%, >= 71.7 kt {np.nansum(v>=71.7)/tot*100:.1f}%, >= 80 kt {np.nansum(v>=80)/tot*100:.1f}%")


    # open ocean: cells with no land within about 200 km (ellipse of 8 rows x 12 columns, a rough 200-km radius at 50-60N)
    from scipy import ndimage
    yy, xx = np.ogrid[-8:9, -12:13]
    foot = (yy / 8.0) ** 2 + (xx / 12.0) ** 2 <= 1
    near_land = ndimage.binary_dilation(~ocean, structure=foot, iterations=1)
    openo = ocean & ~near_land
    out.append("\n== Open ocean (no land within about 200 km) ==")
    for b in ("atl", "pac"):
        d = domain[b] & openo
        for nm, f in (("mean annual maximum gust (kt)", annmax), ("share of times with gust >= 64 kt (%)", freq[2] * 100), ("share of times with gust >= 50 kt (%)", freq[1] * 100)):
            v = np.where(d, f, np.nan)
            i = np.nanargmax(v)
            r, c = np.unravel_index(i, v.shape)
            lon = LON[c] if b == "pac" else (LON[c] - 360 if LON[c] > 180 else LON[c])
            out.append(f"{b}: highest {nm} {np.nanmax(v):.2f} at {lat[r]:.2f}N, {lon:.2f}E")
        v = np.where(d, annmax, np.nan)
        out.append(f"{b}: open-ocean cells with mean annual maximum >= 64 kt: {np.nansum(v>=64)/np.isfinite(v).sum()*100:.1f}%; >= 70 kt: {np.nansum(v>=70)/np.isfinite(v).sum()*100:.1f}%")

    # store a compact copy (15-80N, ocean only) for map layers; float16 / uint16, well under 10 MB
    np.savez_compressed(os.path.join(RES, "gust_climo_grid.npz"), lat=lat.astype(np.float32), lon=LON.astype(np.float32),
                        n_times=n, annual_max_mean_kt=np.nan_to_num(names["annmax"], nan=-1).astype(np.float16),
                        abs_max_kt=np.nan_to_num(names["absmax"], nan=-1).astype(np.float16),
                        freq_thresholds_kt=np.array(THRESH, np.float32),
                        freq=np.stack([np.nan_to_num(names[f"freq_ge{t:g}"], nan=-1) for t in THRESH]).astype(np.float16),
                        note="ERA5 proxy, not observation; 2004-05..2025-26 Oct-Apr 12-hourly; ocean only, -1 = land; 0.25 deg, rows 80N..15N")

    # ------------ figures
    for tag, field, vmin, vmax, lab in (("annmax", names["annmax"], 40, 100, "mean annual maximum gust (kt)"),
                                        ("freq_ge64", names["freq_ge64"] * 100, 0, 3.0, "% of sampled times with gust >= 64 kt"),
                                        ("freq_ge50", names["freq_ge50"] * 100, 0, 15, "% of sampled times with gust >= 50 kt"),
                                        ("freq_ge71.7", names["freq_ge71.7"] * 100, 0, 1.0, "% of sampled times with gust >= 71.7 kt")):
        fig, axs = plt.subplots(1, 2, figsize=(14, 5.2), constrained_layout=True)
        for ax, b in zip(axs, ("atl", "pac")):
            m = region_plot(ax, b, field, vmin, vmax, "magma_r", ocean)
            ax.set_title(A_.BN[b], fontsize=10)
        fig.colorbar(m, ax=axs, shrink=0.8, label=lab)
        fig.suptitle(f"ERA5 gust climatology: {lab}", fontsize=11)
        banner(fig)
        fig.savefig(os.path.join(RES, f"fig15_{tag}.png"), dpi=110)
        plt.close(fig)

    # ------------ where the maximum gust sits in HF storms (hf_structure, 2004-05 on)
    S = pd.read_csv(os.path.join(HERE, "..", "hf_structure", "results", "fixes.csv"))
    S = S.dropna(subset=["basin"])
    out.append("\n== Location of the maximum gust at HF-strength fixes (hf_structure fixes.csv, 2004-05 on) ==")
    fig, axs = plt.subplots(1, 2, figsize=(14, 5.2), constrained_layout=True)
    for ax, b in zip(axs, ("atl", "pac")):
        d = S[S.basin == b]
        lon = A_.lonfix(d.gmax_lon, b)
        lat_e = np.arange(20, 76, 5)
        lon_e = np.arange(-100, 31, 10) if b == "atl" else np.arange(120, 251, 10)
        H, _, _ = np.histogram2d(d.gmax_lat, lon, bins=[lat_e, lon_e])
        m = ax.pcolormesh(lon_e, lat_e, np.ma.masked_less_equal(H / 22.0 * 6, 0), cmap="YlOrRd", zorder=2)
        F_.draw_land(ax, b)
        fig.colorbar(m, ax=ax, shrink=0.8, label="hours per season with the storm's maximum gust in the box")
        k = np.unravel_index(H.argmax(), H.shape)
        out.append(f"{b}: {len(d)} HF fixes; busiest box for the maximum-gust point {lat_e[k[0]]:.0f}-{lat_e[k[0]]+5:.0f}N, {lon_e[k[1]]:.0f} to {lon_e[k[1]]+10:.0f}E "
                   f"({H[k]/22*6:.1f} h per season); median distance of the maximum from the centre {d.gmax_r.median():.0f} km; median maximum gust {d.gmax.median():.1f} kt; "
                   f"share of maxima within 100 km of a coast {float((d.gmax_coast_km < 100).mean()):.3f}")
        ax.set_title(A_.BN[b], fontsize=10)
    fig.suptitle("Where the strongest gust of an HF-strength storm sits (not the centre), 2004-05 to 2025-26", fontsize=11)
    banner(fig, "Storm-owned ocean gust within 1,200 km of the centre.", warn="ERA5 PROXY (pipeline A gust), NOT DIRECT OBSERVATION. HF-strength fixes, all months, 2004-05 to 2025-26.")
    fig.savefig(os.path.join(RES, "fig16_max_gust_location_hf_storms.png"), dpi=110)
    plt.close(fig)

    # ------------ within-era rankings by pipeline A gust index
    C = pd.read_csv(os.path.join(HFH, "era5_hf_catalog.csv"))
    E = C[C.role == "event"].copy()
    LC = pd.read_csv(os.path.join(HFH, "lifecycle_events.csv"))[["track", "tc"]]
    E = E.merge(LC, on="track", how="left")
    E["lon"] = np.where(E.basin == "atl", ((E.peak_lon + 180) % 360) - 180, E.peak_lon % 360)
    E["date"] = E.peak_time.astype(str).str[:8]
    rows = []
    for era, (a, b_) in {"E1 1979-2000": (1979, 2000), "E3 2004-2025": (2004, 2025)}.items():
        e = E[E.season.between(a, b_)].sort_values("gust800_kt", ascending=False).head(25)
        out.append(f"\n== Top 25 storms by pipeline A gust index within {era} (ranks are within this era only) ==")
        for i, (_, r) in enumerate(e.iterrows(), 1):
            out.append(f"{i:2d}. {r.date} {r.basin} gust index {r.gust800_kt:.1f} kt  min MSLP {r.minp:.0f} hPa  centre at peak ({r.peak_lat:.1f}N, {r.lon:.1f}E)  tc-linked={bool(r.tc)}")
            rows.append(dict(era=era, rank=i, date=r.date, basin=r.basin, gust_index_kt=r.gust800_kt, min_mslp_hpa=r.minp,
                             lat=r.peak_lat, lon=r.lon, tc_linked=bool(r.tc), track=int(r.track),
                             caveat="ERA5 proxy, not observation; rank within era only (gust drifts upward before 2001)"))
    pd.DataFrame(rows).to_csv(os.path.join(RES, "top_gust_storms_by_era.csv"), index=False, float_format="%.1f")
    open(os.path.join(RES, "gust_climo.txt"), "w").write("\n".join(out) + "\n")
    print("\n".join(out))


if __name__ == "__main__":
    main()
