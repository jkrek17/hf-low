"""Q3: draw 12 cases (Storm Dennis and ex-Typhoon Nuri if in sample, rest random seed 20261009) with streak, quadrant and trough.

usage: ERA5_WORK=DIR qc_cases.py     pulls 12 chunks (~0.35 GB), writes results/qc_cases/*.png and results/qc_cases.png
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from features import *
from common import decode, WB2, WB2_T0, LEVELS
from qc import load

L = {l: LEVELS.index(l) for l in (250, 300, 500)}
SEED = 20261009


def pick(D):
    M = pd.read_csv(os.path.join(HERE, "..", "intensity", "results", "fixes_2004.csv.gz"), usecols=["track", "time", "basin", "msl"], dtype={"time": str})
    D = D.merge(M[["track", "time", "msl"]], on=["track", "time"])
    t = D.time
    names = []
    w = D[(D.basin == "atl") & (t >= "2020021400") & (t <= "2020021700")]
    names.append(("Storm Dennis (inferred)", w.loc[w.msl.idxmin()]))
    w = D[(D.basin == "pac") & (t >= "2014110600") & (t <= "2014110900")]
    names.append(("ex-Typhoon Nuri (inferred)", w.loc[w.msl.idxmin()]))
    pool = D[(D.nohead == 0) & (D.nojet == 0)]
    r = pool.sample(10, random_state=SEED)
    return names + [(f"random {i + 1}", row) for i, (_, row) in enumerate(r.iterrows())]


def main():
    D = load()
    mon = pd.to_datetime(D.time, format="%Y%m%d%H").dt.month
    D = D[mon.isin([10, 11, 12, 1, 2, 3, 4])].reset_index(drop=True)
    cases = pick(D)
    os.makedirs(os.path.join(HERE, "results", "qc_cases"), exist_ok=True)
    fig, axs = plt.subplots(6, 4, figsize=(18, 26))
    for k, (nm, r) in enumerate(cases):
        when = np.datetime64(f"{r.time[:4]}-{r.time[4:6]}-{r.time[6:8]}T{r.time[8:10]}")
        step = int((when - WB2_T0) / np.timedelta64(6, "h"))
        A = {tag: decode(WB2, var, f"{step // 8}.0.0.0")[0] for var, tag in (("u_component_of_wind", "u"), ("v_component_of_wind", "v"), ("geopotential", "z"))}
        ti = step % 8
        Fd = time_fields(A["u"][ti, L[250]], A["v"][ti, L[250]], A["u"][ti, L[300]], A["v"][ti, L[300]], A["z"][ti, L[500]])
        dbg = {}
        o = fix_features(Fd, r.lat, r.lon, r.heading, dbg)
        ext = [-4050, 4050, -4050, 4050]
        a1, a2 = axs[k // 2, (k % 2) * 2], axs[k // 2, (k % 2) * 2 + 1]
        if "V" in dbg:
            im = a1.imshow(dbg["V"], origin="lower", extent=ext, cmap="viridis", vmin=0, vmax=150)
            a1.contour(BX, BY, dbg["comp"].astype(float), [0.5], colors="w", linewidths=1.5)
            xm, ym = dbg["jmax"]
            a1.plot([xm], [ym], "r*", ms=14)
            a = dbg["a"]
            a1.annotate("", xy=(xm + 1500 * a[0], ym + 1500 * a[1]), xytext=(xm - 1500 * a[0], ym - 1500 * a[1]), arrowprops=dict(arrowstyle="->", color="r", lw=2))
        a1.plot([0], [0], "kx", ms=12, mew=3)
        q = [c for c in ("RE", "LE", "RX", "LX") if o[c] == 1]
        a1.set_title(f"{nm} {r.time} {r.basin} {r.lat:.0f}N heading {r.heading:.0f}\n250 hPa vortex-removed (kt), streak, quadrant {q[0] if q else 'none'}, s {o['s']:+.1f} n {o['n']:+.1f}", fontsize=8)
        if "zs" in dbg:
            a2.imshow(dbg["zs"], origin="lower", extent=ext, cmap="RdBu_r", vmin=-400, vmax=400)
            xt, yt = dbg["trough"]
            a2.plot([xt], [yt], "k^", ms=12)
            hd = np.radians(r.heading)
            a2.annotate("", xy=(1500 * np.sin(hd), 1500 * np.cos(hd)), xytext=(0, 0), arrowprops=dict(arrowstyle="->", color="g", lw=2))
        a2.plot([0], [0], "kx", ms=12, mew=3)
        a2.set_title(f"Z500' (m), trough depth {o['tdepth']:.0f} dist {o['tdist']:.1f} tilt {o['ttilt']:.0f} tphase {o['tphase']:+.1f}\n(green arrow = motion)", fontsize=8)
        for a_ in (a1, a2):
            a_.set_xlim(-4000, 4000); a_.set_ylim(-4000, 4000)
    plt.tight_layout()
    plt.savefig(os.path.join(HERE, "results", "qc_cases.png"), dpi=55)
    print("saved")
    pd.DataFrame([dict(case=n, **r[["track", "time", "basin", "lat", "lon", "heading"]].to_dict()) for n, r in cases]).to_csv(
        os.path.join(HERE, "results", "qc_cases.csv"), index=False)


if __name__ == "__main__":
    main()
