"""Figures for hf_vs_storm from results/large and results/storm. usage: make_figs.py REPO_ROOT
ERA5 proxy, pipeline A. Land outline in the large-scale maps is the 0.25 degree land-sea mask contour (ARCO-ERA5)."""
import sys, os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
root = sys.argv[1]
RES = os.path.join(root, "research/era5/hf_vs_storm/results")
FIG = os.path.join(RES, "figs"); os.makedirs(FIG, exist_ok=True)
WORK = os.path.join(os.environ.get("ERA5_WORK", "/tmp/era5work"), "hf_vs_storm")
lsm = np.load(os.path.join(WORK, "lsm.npy"))[::4, ::4]                # 1 degree: rows 90..-90, lon 0..359
llat = 90 - np.arange(0, 721, 4) * 0.25; llon = np.arange(0, 1440, 4) * 0.25
BASIN = {"atl": "North Atlantic", "pac": "North Pacific"}
VAR = {"z500": ("Z500 anomaly", "m", 30), "mslp": ("MSLP anomaly", "hPa", 6), "u250": ("250 hPa zonal wind anomaly", "m/s", 6),
       "sst": ("SST anomaly", "K", 0.6), "tcwv": ("column water vapour anomaly", "kg/m2", 1.2)}
CFG = {"L": "days -10 to -4", "k0": "onset day", "k7": "day -7", "k4": "day -4", "k2": "day -2", "k1": "day -1"}
EXT = {"atl": (-100, 40, 25, 87), "pac": (110, 250, 25, 87)}

def shift(lon, b):
    return np.where(lon > 180, lon - 360, lon) if b == "atl" else lon

def lonaxis(lon, b):
    l = shift(lon, b)
    o = np.argsort(l)
    return l[o], o

def coast(ax, b):
    l = shift(llon, b); o = np.argsort(l)
    ax.contour(l[o], llat, (lsm[:, o] < 0.5).astype(float), levels=[0.5], colors="k", linewidths=0.5)

def largescale(cmpname="C1"):
    for b in ("atl", "pac"):
        M = np.load(os.path.join(RES, f"large/maps_{b}_{cmpname}.npz"))
        cfgs, vars_ = list(M["cfgs"]), list(M["vars"])
        lon = M["lon"]; lat = M["lat"]; lx, o = lonaxis(lon, b)
        for cfg in ("L", "k0"):
            fig, axes = plt.subplots(1, 5, figsize=(22, 3.9), constrained_layout=True)
            for ax, v in zip(axes, ["z500", "mslp", "u250", "sst", "tcwv"]):
                iv, ic = vars_.index(v), cfgs.index(cfg)
                d = M["diff"][ic, iv][:, o]; rej = M["rej"][ic, iv][:, o]
                t, u, s = VAR[v]
                im = ax.pcolormesh(lx, lat, d, cmap="RdBu_r", norm=TwoSlopeNorm(0, -s, s), shading="nearest")
                yy, xx = np.where(rej)
                ax.plot(lx[xx], lat[yy], "k.", ms=3)
                coast(ax, b); x0, x1, y0, y1 = EXT[b]; ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
                bx = M["box"][:, o]; ax.contour(lx, lat, bx.astype(float), levels=[0.5], colors="g", linewidths=1.2)
                ax.set_title(f"{t} ({u})", fontsize=10); fig.colorbar(im, ax=ax, shrink=0.8, pad=0.01)
            fig.suptitle(f"{BASIN[b]}: HF onset minus storm-force peak, {CFG[cfg]} before the anchor (ERA5 proxy, pipeline A; n HF {int(M['n_hf'])}, SF {int(M['n_sf'])}; dots: BH q<0.05; green: test box)", fontsize=11)
            fig.savefig(os.path.join(FIG, f"large_diff_{b}_{cmpname}_{cfg}.png"), dpi=90); plt.close(fig)
        # typical pattern: raw Z500 contours + U250 shading, HF and SF, L and k0
        fig, axes = plt.subplots(2, 3, figsize=(17, 7.6), constrained_layout=True)
        iz, iu = vars_.index("z500_raw"), vars_.index("u250_raw")
        for r, (key, lab) in enumerate((("hf", "HF onset"), ("sf", "storm-force peak"))):
            for c, cfg in enumerate(("L", "k0")):
                ax = axes[r, c]; ic = cfgs.index(cfg)
                u = M[key][ic, iu][:, o]; z = M[key][ic, iz][:, o]
                im = ax.pcolormesh(lx, lat, u, cmap="viridis", vmin=10, vmax=50, shading="nearest")
                cs = ax.contour(lx, lat, z, levels=np.arange(4800, 5900, 60), colors="w", linewidths=0.9)
                ax.clabel(cs, fontsize=6, fmt="%d"); coast(ax, b); x0, x1, y0, y1 = EXT[b]; ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
                ax.set_title(f"{lab}: {CFG[cfg]}", fontsize=10)
            # difference of raw Z500 at L for context column
            ic = cfgs.index("L")
        for r in range(2):
            ax = axes[r, 2]; key = ("hf", "sf")[r]; ic = cfgs.index("L")
            im2 = ax.pcolormesh(lx, lat, M[key][ic, vars_.index("mslp_raw")][:, o], cmap="coolwarm", vmin=990, vmax=1030, shading="nearest")
            coast(ax, b); x0, x1, y0, y1 = EXT[b]; ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.set_title(f"{('HF onset', 'storm-force peak')[r]}: MSLP, days -10 to -4", fontsize=10)
        fig.colorbar(im, ax=axes[:, :2], shrink=0.8, label="250 hPa zonal wind (m/s); white: Z500 (m)"); fig.colorbar(im2, ax=axes[:, 2], shrink=0.8, label="MSLP (hPa)")
        fig.savefig(os.path.join(FIG, f"large_typical_{b}.png"), dpi=90); plt.close(fig)
    S = pd.read_csv(os.path.join(RES, "large/series.csv"))
    for cmpname in ("C1",):
        fig, axes = plt.subplots(2, 5, figsize=(21, 6.2), constrained_layout=True, sharex=True)
        for r, b in enumerate(("atl", "pac")):
            for c, v in enumerate(["z500", "mslp", "u250", "sst", "tcwv"]):
                ax = axes[r, c]; d = S[(S.basin == b) & (S.cmp == cmpname) & (S["var"] == v)].sort_values("k")
                ax.fill_between(-d.k, d.lo, d.hi, color="C0", alpha=0.25); ax.plot(-d.k, d["diff"], "C0-o", ms=3); ax.axhline(0, color="k", lw=0.6)
                ax.set_title(f"{BASIN[b]}: {VAR[v][0]} ({VAR[v][1]})", fontsize=9); ax.axvspan(-10, -4, color="y", alpha=0.12)
        fig.supxlabel("days before the anchor (shaded: lead window); box mean, HF onset minus storm-force peak, 95% season-block bootstrap interval")
        fig.savefig(os.path.join(FIG, f"large_series_{cmpname}.png"), dpi=90); plt.close(fig)

def stormscale():
    X = np.arange(121) * 25 - 1500
    cfg = {"gust": ("10 m gust (kt)", "magma_r", 30, 90, 10), "ws": ("10 m wind speed (kt)", "magma_r", 15, 60, 8), "msl": ("MSLP (hPa)", "viridis", 940, 1020, 20),
           "d2m": ("2 m dewpoint (C)", "YlGnBu", -15, 18, 10)}
    for b in ("atl", "pac"):
        A = np.load(os.path.join(RES, "storm/comp_%s_C1.npz" % b)); Bc = np.load(os.path.join(RES, "storm/comp_%s_C2.npz" % b))
        fig, axes = plt.subplots(4, 5, figsize=(22, 16), constrained_layout=True)
        for r, f in enumerate(["gust", "ws", "msl", "d2m"]):
            lab, cm, v0, v1, _ = cfg[f]
            pan = [(A[f"{f}_m_hf"], "HF onset"), (Bc[f"{f}_m_hf"], "HF peak"), (A[f"{f}_m_sf"], "storm-force peak")]
            for c, (arr, t) in enumerate(pan):
                ax = axes[r, c]; im = ax.pcolormesh(X, X, arr, cmap=cm, vmin=v0, vmax=v1, shading="nearest"); ax.set_title(f"{t}: {lab}", fontsize=10)
                ax.plot(0, 0, "w+", ms=10); ax.set_aspect("equal")
            fig.colorbar(im, ax=axes[r, :3], shrink=0.8)
            dm = max(np.nanpercentile(np.abs(A[f"{f}_m_diff"]), 99), 1e-6)
            for c, (S_, t) in enumerate(((A, "HF onset minus SF peak"), (Bc, "HF peak minus SF peak"))):
                ax = axes[r, 3 + c]; im2 = ax.pcolormesh(X, X, S_[f"{f}_m_diff"], cmap="RdBu_r", norm=TwoSlopeNorm(0, -dm, dm), shading="nearest")
                yy, xx = np.where(S_[f"{f}_m_rej"]); sel = slice(None, None, 11)
                ax.plot(X[xx][sel], X[yy][sel], "k.", ms=1.5); ax.set_aspect("equal"); ax.set_title(t + ": " + lab, fontsize=10); fig.colorbar(im2, ax=ax, shrink=0.8)
            for ax in axes[r]: ax.set_xlabel("km ahead of motion (+x)"); ax.set_ylabel("km left of motion (+y)")
        fig.suptitle(f"{BASIN[b]}: storm-motion-relative composites, ERA5 proxy (pipeline A). n HF onset/peak {int(A['n_hf_rot'])}, SF {int(A['n_sf_rot'])} with a heading. Right of motion is down. "
                     "ERA5 reads low in extreme storms. Dots: BH q<0.05", fontsize=12)
        fig.savefig(os.path.join(FIG, f"storm_comp_{b}.png"), dpi=70); plt.close(fig)

if __name__ == "__main__":
    which = sys.argv[2] if len(sys.argv) > 2 else "all"
    if which in ("all", "large"): largescale()
    if which in ("all", "storm"): stormscale()
