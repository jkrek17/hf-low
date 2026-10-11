"""Figures and summary numbers for the basin composites. usage: ERA5_WORK=DIR figs.py   (needs DIR/basin_sums.npz)
Writes results/basin_atl.png, results/basin_pac.png, results/summary.txt, results/means.npz (mean fields and anomalies per group)."""
import os, sys, re, json
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(HERE, "..", "hf_env_composites"))
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from common import LAT, LON, WORK
sys.path.insert(0, HERE)
import extract as X
txt = open(os.path.join(REPO, "docs/data/coastlines.js")).read()
POLY = json.loads(txt[txt.index("{"): txt.rindex("}") + 1])["polygons"]
WIN = {"atl": (-100, 40, 20, 85), "pac": (100, 250, 20, 85)}
NAME = {"atl": "North Atlantic", "pac": "North Pacific"}
COLS = [("msl", "a_msl", "Sea-level pressure", "hPa", 4, 8), ("z500", "a_z500", "500 hPa height", "m", 60, 100), ("ws250", "a_ws250", "250 hPa wind speed", "m/s", 10, 12)]
IDX = {n: i for i, n in enumerate(["msl", "z500", "ws250", "a_msl", "a_z500", "a_ws250"])}


def sub(a, lon0, lon1, lat0, lat1):
    lons = np.arange(lon0, lon1 + 0.01, 1.5)
    ci = np.rint((lons % 360) / 1.5).astype(int) % 240
    la = (LAT >= lat0) & (LAT <= lat1)
    return a[np.ix_(la, ci)], lons, LAT[la]


def coast(ax, lon0, lon1):
    for p in POLY:
        a = np.array(p, float); lo, la = a[:, 0].copy(), a[:, 1]
        if lon1 > 180:
            lo = lo % 360
        wrap = np.where(np.abs(np.diff(lo)) > 180)[0] + 1
        ax.plot(np.insert(lo, wrap, np.nan), np.insert(la, wrap, np.nan), lw=0.5, color="0.25")


def main():
    z = np.load(os.path.join(WORK, "basin_sums.npz"))
    S = pd.read_csv(os.path.join(HERE, "..", "hf_env_composites", "results", "storm_times.csv"), dtype={"time": str})
    on = S[S.anchor == "onset"]
    out, lines, means = {}, [], {}
    # wind-speed anomaly = mean |V| minus the climatological mean |V| at each event's own day of year and hour
    # (the first version used |V_clim| of the mean vector, which is biased high everywhere)
    cws = np.load(os.path.join(WORK, "clim_ws250.npz"))["ws250"]
    EV = X.events(); csum = {}
    for r in EV.itertuples():
        w = np.datetime64(pd.Timestamp(r.when), "h"); h = int((w - w.astype("datetime64[D]")).astype(int)) // 6
        doy = int((w.astype("datetime64[D]") - w.astype("datetime64[Y]").astype("datetime64[D]")).astype(int)) + 1
        k = (r.basin, r.anchor); csum[k] = csum.get(k, 0.0) + cws[h, doy - 1].astype(np.float64)
    for b in ("atl", "pac"):
        lon0, lon1, lat0, lat1 = WIN[b]
        fig, axs = plt.subplots(2, 3, figsize=(15.5, 7.4), constrained_layout=True)
        for r, (anc, lab) in enumerate((("m48", "48 hours before onset"), ("onset", "HF onset"))):
            n = int(z[f"n_{b}|{anc}"]); M = z[f"s_{b}|{anc}"] / n
            M[IDX["a_ws250"]] = M[IDX["ws250"]] - csum[(b, anc)] / n
            means[f"{b}_{anc}"] = M.astype(np.float32)
            for c, (fk, ak, title, unit, cint, alim) in enumerate(COLS):
                ax = axs[r, c]
                A, lons, lats = sub(M[IDX[ak]], lon0, lon1, lat0, lat1); F = sub(M[IDX[fk]], lon0, lon1, lat0, lat1)[0]
                lim = alim
                im = ax.contourf(lons, lats, A, levels=np.linspace(-lim, lim, 21), cmap="RdBu_r", extend="both")
                lv = np.arange(np.floor(F.min() / cint) * cint, F.max() + cint, cint)
                if fk == "ws250":
                    lv = lv[lv >= 20]
                cs = ax.contour(lons, lats, F, levels=lv, colors="k", linewidths=0.6)
                ax.clabel(cs, fontsize=6, fmt="%d")
                coast(ax, lon0, lon1)
                if anc == "onset":
                    o = on[on.basin == b]; lo = o.lon.values.copy()
                    lo = np.where(lo > 180, lo - 360, lo) if lon1 <= 180 else lo
                    ax.plot(lo, o.lat.values, ".", ms=2.5, color="m", alpha=0.35)
                ax.set_xlim(lon0, lon1); ax.set_ylim(lat0, lat1); ax.set_aspect('auto')
                ax.set_title(f"{title}: {lab}\nn = {n} HF lows", fontsize=9)
                ax.set_xticks(np.arange(lon0, lon1 + 1, 20 if b == "atl" else 30)); ax.tick_params(labelsize=7)
                ax.set_xticklabels([f"{abs(((x + 180) % 360) - 180):.0f}{'E' if 0 < ((x + 180) % 360) - 180 < 180 else 'W'}" for x in ax.get_xticks()])
                if r == 1:
                    cb = fig.colorbar(im, ax=axs[:, c], shrink=0.6, location="bottom", pad=0.02)
                    cb.set_label(f"anomaly from 1990-2019 climate ({unit}); contours: the mean field", fontsize=8); cb.ax.tick_params(labelsize=7)
                # numbers
                sel = (lats >= 25) & (lats <= 75)
                A = np.where(sel[:, None], A, np.nan) if False else A
                if fk in ("z500", "msl"):
                    k = np.unravel_index(np.nanargmin(np.where(sel[:, None], A, np.nan)), A.shape)
                    lines.append(f"{b} {anc} {fk}: deepest mean anomaly {A[k]:+.1f} {unit} at {lats[k[0]]:.1f}N, {((lons[k[1]] + 180) % 360) - 180:.1f}E; mean field min {F.min():.1f}")
                else:
                    k = np.unravel_index(np.nanargmax(np.where(sel[:, None], A, np.nan)), A.shape)
                    lines.append(f"{b} {anc} {fk}: largest mean anomaly {A[k]:+.1f} {unit} at {lats[k[0]]:.1f}N, {((lons[k[1]] + 180) % 360) - 180:.1f}E; mean field max {F.max():.1f}")
        fig.suptitle(f"{NAME[b]}: typical pattern for hurricane-force (HF) lows, Oct-Apr onsets 2004-05 to 2022-23 (ERA5 proxy, pipeline A, 1.5 degrees). Shading: anomaly; magenta dots: where each HF onset happened", fontsize=10)
        fig.savefig(os.path.join(HERE, "results", f"basin_{b}.png"), dpi=88)
        plt.close(fig)
        o = on[on.basin == b]
        lines.append(f"{b} onset latitude: mean {o.lat.mean():.1f}N, 10th to 90th percentile {o.lat.quantile(.1):.1f} to {o.lat.quantile(.9):.1f}N; n {len(o)}")
    np.savez_compressed(os.path.join(HERE, "results", "means.npz"), **means)
    open(os.path.join(HERE, "results", "summary.txt"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
