"""Figures for the three follow-up tests (hf_lifecycle, hf_pattern_phase, hf_conversion). ERA5 proxy, pipeline A. usage: figs_followup.py REPO_ROOT"""
import sys, os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
root = sys.argv[1]; E5 = os.path.join(root, "research/era5")
BASIN = {"atl": "North Atlantic", "pac": "North Pacific"}; COL = {"atl": "#1f77b4", "pac": "#d95f02"}
AX = (np.arange(121) - 60) * 25.0

def series(R, vars_, labels, x, fn, title, ylab="difference in storm-to-storm SD", xlab="hours from peak"):
    fig, axes = plt.subplots(1, len(vars_), figsize=(4.4 * len(vars_), 3.6), constrained_layout=True, sharey=False)
    for ax, v, lab in zip(np.atleast_1d(axes), vars_, labels):
        for b in ("atl", "pac"):
            d = R[(R.basin == b) & (R["var"] == v)].sort_values(x)
            s = d["sd"].values
            ax.plot(d[x], d["diff"] / s, "o-", c=COL[b], label=BASIN[b]); ax.fill_between(d[x], d.lo / s, d.hi / s, color=COL[b], alpha=0.18)
            sig = d.q.values < 0.05; ax.plot(d[x][sig], (d["diff"] / s)[sig], "o", c=COL[b], ms=9, mfc="none", mew=1.5)
        ax.axhline(0, c="k", lw=0.7); ax.set_title(lab, fontsize=10); ax.set_xlabel(xlab); ax.grid(alpha=0.3)
    np.atleast_1d(axes)[0].set_ylabel(ylab); np.atleast_1d(axes)[0].legend(fontsize=8)
    fig.suptitle(title + " (rings: BH q<0.05; band: 95% season-bootstrap interval)", fontsize=10); fig.savefig(fn, dpi=90); plt.close(fig)

def forest(R, labels, fn, title, col="d_std", lo="lo", hi="hi", sdcol="sd"):
    fig, axes = plt.subplots(1, 2, figsize=(10, 0.42 * len(labels) + 1.6), constrained_layout=True, sharey=True)
    for ax, b in zip(axes, ("atl", "pac")):
        d = R[R.basin == b].set_index("var").loc[list(labels)]
        y = np.arange(len(d))[::-1]; s = d[sdcol].values
        ax.errorbar(d[col], y, xerr=[d[col] - d[lo] / s, d[hi] / s - d[col]], fmt="none", ecolor="grey")
        ax.scatter(d[col], y, c=np.where(d.q < 0.05, "k", "w"), edgecolors="k", zorder=3)
        ax.axvline(0, c="k", lw=0.7); ax.set_yticks(y); ax.set_yticklabels([labels[k] for k in d.index]); ax.set_title(BASIN[b]); ax.set_xlabel("difference in SD (filled: BH q<0.05)"); ax.grid(alpha=0.3)
    fig.suptitle(title, fontsize=10); fig.savefig(fn, dpi=90); plt.close(fig)

def mapgrid(C, rows, fn, title, ngrp):
    """rows: list of (key, label, unit, vmin_raw, vmax_raw, dlim, cmap). Columns: group A, group B, difference."""
    fig, axes = plt.subplots(len(rows), 3, figsize=(11.5, 3.5 * len(rows)), constrained_layout=True)
    for r, (k, lab, unit, v0, v1, dl, cm) in enumerate(rows):
        for c, (nm, tt) in enumerate(((ngrp[0], ngrp[2]), (ngrp[1], ngrp[3]), ("diff", "difference"))):
            ax = axes[r, c]; A = C[f"{k}_{nm}"]
            if nm == "diff":
                im = ax.pcolormesh(AX, AX, A, cmap="RdBu_r", norm=TwoSlopeNorm(0, -dl, dl), shading="nearest")
                yy, xx = np.where(C[f"{k}_rej"]); ax.plot(AX[xx[::5]], AX[yy[::5]], "k.", ms=1.5)
            else:
                im = ax.pcolormesh(AX, AX, A, cmap=cm, vmin=v0, vmax=v1, shading="nearest")
            ax.plot(0, 0, "k+", ms=10); ax.set_aspect("equal"); ax.set_title(f"{lab} ({unit}): {tt}", fontsize=9); fig.colorbar(im, ax=ax, shrink=0.8)
            ax.invert_yaxis() if False else None
    fig.suptitle(title + "  +x is the direction of motion; dots: pixel BH q<0.05 (spatially correlated, descriptive)", fontsize=10); fig.savefig(fn, dpi=85); plt.close(fig)

# ---------------- test 3
LC = os.path.join(E5, "hf_lifecycle/results")
b = pd.read_csv(os.path.join(LC, "stage_b_tests.csv"))
series(b, ["pc", "msl_grad", "g48_rmax", "gmax_r"], ["central MSLP (hPa)", "MSLP gradient (hPa/100 km)", "48-kt radius (km)", "radius of max gust (km)"], "lag",
       os.path.join(LC, "figs/lifecycle_stageB.png"), "Test 3 stage B: HF minus storm-force-only, hours around the peak (ERA5 proxy, pipeline A)")
a = pd.read_csv(os.path.join(LC, "stage_a_tests.csv")); a = a[a["sel"] == "all"]
series(a[a["var"].isin(["msl", "dp12"])], ["msl", "dp12"], ["central MSLP at 00/12 UTC (hPa)", "12 h pressure change (hPa)"], "lag", os.path.join(LC, "figs/lifecycle_stageA.png"),
       "Test 3 stage A: track table, all storms with a value at each lag")
# ---------------- test 2
PP = os.path.join(E5, "hf_pattern_phase/results")
P = pd.read_csv(os.path.join(PP, "pp_primary.csv"))
LAB = {"pc": "central MSLP", "msl_grad": "MSLP gradient", "gmax": "max gust", "wsmax": "max 10 m wind", "a_g48": "48-kt gust area", "gmax_r": "radius of max gust", "g48_rmax": "48-kt outer radius",
       "g48_right": "48-kt share right of motion", "d2m_500": "dewpoint (500 km)", "gust_factor": "gust factor", "a_hf": "HF-equivalent area"}
forest(P, LAB, os.path.join(PP, "figs/pp_forest.png"), "Test 2: HF onsets under the top minus bottom tercile of the pattern index, structure in SD (primary family, ERA5 proxy)")
for bsn in ("atl", "pac"):
    C = np.load(os.path.join(PP, f"comp_{bsn}.npz"))
    mapgrid(C, [("gust", "gust", "kt", 30, 90, 8, "viridis"), ("msl", "MSLP", "hPa", 960, 1010, 4, "viridis_r"), ("d2m", "2 m dewpoint", "degC", -5, 15, 2, "viridis")], os.path.join(PP, f"figs/pp_maps_{bsn}.png"),
            f"Test 2 {BASIN[bsn]}: HF onset, top (n {int(C['n_top'])}) vs bottom (n {int(C['n_bottom'])}) tercile", ("top", "bottom", "top tercile", "bottom tercile"))
# ---------------- test 1
CV = os.path.join(E5, "hf_conversion/results")
sa = pd.read_csv(os.path.join(CV, "stage_a_tests.csv")); sa["lagh"] = -sa.lag.map({"t0": 0, "t0-12h": 12, "t0-24h": 24})
series(sa, ["eady", "flux", "vadv500", "div300", "jet250", "tcwv"], ["Eady growth rate", "surface heat flux", "500 hPa vorticity advection", "300 hPa divergence", "250 hPa jet", "column water vapour"], "lagh",
       os.path.join(CV, "figs/conv_stageA.png"), "Test 1 stage A: converters minus non-converters, hours before t0 (env table)", xlab="hours from t0")
sb = pd.read_csv(os.path.join(CV, "stage_b_tests.csv"))
series(sb.rename(columns={"conv": "hf"}), ["jet_dist", "jet_cos", "thetae850", "stab", "baroc"], ["jet-max distance (km)", "jet bearing cos (+1 ahead)", "theta-e 850 (K)", "SST - T2m (K)", "|grad T850| (K/100 km)"], "lag",
       os.path.join(CV, "figs/conv_stageB.png"), "Test 1 stage B: converters minus non-converters (0.25 deg, 250 pairs/basin)", xlab="hours from t0")
for bsn in ("atl", "pac"):
    for L in (0, 24):
        p = os.path.join(CV, f"comp_{bsn}_L{L}.npz")
        if os.path.exists(p):
            C = np.load(p)
            mapgrid(C, [("ws250", "250 hPa wind speed", "m/s", 10, 70, 8, "viridis"), ("te", "850 hPa theta-e", "K", 270, 320, 5, "magma"), ("gt", "|grad T850|", "K/100 km", 0, 4, 0.6, "viridis")],
                    os.path.join(CV, f"figs/conv_maps_{bsn}_L{L}.png"), f"Test 1 {BASIN[bsn]}: t0{'-24 h' if L else ''}, converters (n {int(C['n_conv'])}) vs non-converters (n {int(C['n_non'])})", ("conv", "non", "converters", "non-converters"))
for bsn in ("atl", "pac"):
    for fld, lab, unit, v0, v1, dl, cm in (("gust", "gust", "kt", 30, 90, 10, "viridis"), ("msl", "MSLP", "hPa", 950, 1010, 8, "viridis_r")):
        fig_rows = []
        for L in (-48, -24, 0):
            C = np.load(os.path.join(LC, f"comp_{bsn}_lag{L}.npz")); fig_rows.append((L, C))
        fig, axes = plt.subplots(3, 3, figsize=(11.5, 10.5), constrained_layout=True)
        for r, (L, C) in enumerate(fig_rows):
            for c, (nm, tt) in enumerate((("hf", f"HF (n {int(C['n_hf'])})"), ("sf", f"storm-force only (n {int(C['n_sf'])})"), ("diff", "difference"))):
                ax = axes[r, c]; A = C[f"{fld}_{nm}"]
                if nm == "diff":
                    im = ax.pcolormesh(AX, AX, A, cmap="RdBu_r", norm=TwoSlopeNorm(0, -dl, dl), shading="nearest"); yy, xx = np.where(C[f"{fld}_rej"]); ax.plot(AX[xx[::5]], AX[yy[::5]], "k.", ms=1.5)
                else:
                    im = ax.pcolormesh(AX, AX, A, cmap=cm, vmin=v0, vmax=v1, shading="nearest")
                ax.plot(0, 0, "k+", ms=10); ax.set_aspect("equal"); ax.set_title(f"{L:+d} h: {lab} ({unit}), {tt}", fontsize=9); fig.colorbar(im, ax=ax, shrink=0.8)
        fig.suptitle(f"Test 3 {BASIN[bsn]}: {lab}, HF vs storm-force-only, hours from peak; +x is the direction of motion; dots: pixel BH q<0.05 (descriptive)", fontsize=10)
        fig.savefig(os.path.join(LC, f"figs/lifecycle_maps_{bsn}_{fld}.png"), dpi=80); plt.close(fig)
print("figures done")
