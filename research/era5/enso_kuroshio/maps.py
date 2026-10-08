"""Descriptive genesis-density maps (no significance claim): CP minus EP El Nino DJF genesis per winter,
and the EMI slope per SD (OLS with N34 and year). 5 degree bins, ocean first fixes, DJF window.
usage: maps.py TRACKS WINTER_INDICES LSM_NPY OUT_DIR
"""
import sys, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

T = pd.read_csv(sys.argv[1], dtype={"t0": str}); I = pd.read_csv(sys.argv[2]).set_index("winter")
lsm = np.load(sys.argv[3]); out = sys.argv[4]
T = T[T.ocean0 & (T.t0 >= T.winter.astype(str) + "120300") & (T.t0 < (T.winter + 1).astype(str) + "030100")]
la_e = np.arange(20, 56, 5); lo_e = np.arange(110, 191, 5)
H = np.zeros((47, len(la_e) - 1, len(lo_e) - 1))
for s, g in T.groupby("winter"):
    H[s - 1979] = np.histogram2d(g.lat0, g.lon0, bins=[la_e, lo_e])[0]
w = np.arange(1979, 2026)
en = I.n34 >= 0.5; cp = (en & (I.n4 > I.n3)).values; ep = (en & (I.n4 <= I.n3)).values
neu = (I.n34.abs() < 0.5).values
d_cp_ep = H[cp].mean(0) - H[ep].mean(0)
z = lambda v: (v - v.mean()) / v.std(ddof=1)
X = np.column_stack([np.ones(47), z(I.emi.values), z(I.n34.values), z(w.astype(float))])
beta = np.linalg.lstsq(X, H.reshape(47, -1), rcond=None)[0][1].reshape(H.shape[1:])
fig, ax = plt.subplots(1, 3, figsize=(15, 3.8), sharey=True)
def show(a, A, title, vmax, cmap="RdBu_r", unit=""):
    im = a.pcolormesh(lo_e, la_e, A, cmap=cmap, vmin=-vmax, vmax=vmax, shading="flat")
    a.contour(np.arange(0, 360, .25), 90 - .25 * np.arange(721), lsm, [0.5], colors="k", linewidths=.6)
    a.set_xlim(110, 190); a.set_ylim(20, 55); a.set_title(title, fontsize=10); plt.colorbar(im, ax=a, label=unit, shrink=.85)
    a.set_xlabel("lon (E)")
show(ax[0], H[neu].mean(0), f"mean DJF genesis, neutral winters (n={neu.sum()})", H[neu].mean(0).max(), "viridis" and "Reds", "genesis per winter per 5x5 deg")
ax[0].collections[0].set_clim(0, H[neu].mean(0).max())
show(ax[1], d_cp_ep, f"CP (n={cp.sum()}) minus EP (n={ep.sum()}) El Nino", max(1.0, np.abs(d_cp_ep).max()), unit="per winter")
show(ax[2], beta, "EMI slope per SD (N34, year held)", max(0.5, np.abs(beta).max()), unit="per winter per SD")
ax[0].set_ylabel("lat (N)")
for a, (a0, a1, o0, o1) in zip([ax[1], ax[2]] * 1, [(25, 35, 120, 142)] * 2): pass
fig.suptitle("Descriptive only (ERA5 proxy, MSLP lows, ocean first fix, 3 Dec-28/29 Feb). Not a significance map.", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(out, "genesis_maps.png"), dpi=130)
pd.DataFrame({"cp_minus_ep": d_cp_ep.ravel(), "emi_slope": beta.ravel(),
              "lat0": np.repeat(la_e[:-1], len(lo_e) - 1), "lon0": np.tile(lo_e[:-1], len(la_e) - 1)}).to_csv(os.path.join(out, "map_bins.csv"), index=False)
print("classes: CP", cp.sum(), "EP", ep.sum(), "neutral", neu.sum(), "max |CP-EP|", np.abs(d_cp_ep).max().round(2), "max |slope|", np.abs(beta).max().round(2))
