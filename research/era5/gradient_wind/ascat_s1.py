"""Secondary check S1: ASCAT (equivalent-neutral) wind against geostrophic wind from ERA5 MSLP, western Atlantic, by gradient bin.
Usage: python3 ascat_s1.py OUTFILE. ASCAT cells lie within 800 km of a pipeline A fix, so they are storm-neighbourhood cells; ASCAT saturates near 45-50 kt."""
import glob, os, sys, datetime as dt
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_pairs as bp

SH = bp.SH


def main(out):
    lsm = np.load(f"{SH}/era5-cache/lsm_nh_0p25.npy")
    land = (lsm >= 0.5).astype(float); c = np.cumsum(np.cumsum(np.pad(land, ((1, 0), (1, 0))), 0), 1)
    msl = {}
    for fn in sorted(glob.glob(f"{SH}/gulf-stream/msl/msl_*.npz")):
        z = np.load(fn)
        for t, a in zip(z["t"], z["msl"]): msl[int(t)] = a
    rows = []
    for fn in sorted(glob.glob(f"{SH}/ascat-check/cells/cells_*.npz")):
        r = np.load(fn)["rows"]
        lat = -89.875 + 0.25 * r[:, 0]; lon = 0.125 + 0.25 * r[:, 1]
        m = (lat >= 30) & (lat <= 45) & (lon >= 282) & (lon <= 305) & (r[:, 4] == 0) & (r[:, 3] <= 250)
        r = r[m]
        if len(r): rows.append(r)
    r = np.concatenate(rows); lat = -89.875 + 0.25 * r[:, 0]; lon = 0.125 + 0.25 * r[:, 1]
    t = r[:, 2].astype("int64") * 60
    tm = t.astype("datetime64[s]"); h = (tm.astype("datetime64[h]")).astype("int64"); hod = h % 24
    near = np.where(hod % 12 >= 6, hod + (12 - hod % 12), hod - hod % 12)                       # nearest 00/12 (hours)
    base = (h - hod) + near
    dtm = np.abs(t / 3600.0 - base)
    sel = dtm <= 1.5
    ts = base[sel].astype("datetime64[h]").astype(dt.datetime)
    keys = np.array([x.year * 1000000 + x.month * 10000 + x.day * 100 + x.hour for x in ts])
    r = r[sel]; lat = lat[sel]; lon = lon[sel]
    # ocean, 1 degree from land
    ri = np.rint((90 - lat) / 0.25).astype(int); ci = np.rint(lon / 0.25).astype(int)
    land_n = c[ri + 5, ci + 5] - c[ri - 4, ci + 5] - c[ri + 5, ci - 4] + c[ri - 4, ci - 4]
    ok = land_n == 0
    wind = r[:, 3] * 0.388768
    G = np.full(len(r), np.nan); grad = np.full(len(r), np.nan); win = np.zeros(len(r), int)
    for k in np.unique(keys):
        if k not in msl: continue
        g, GG, _, _ = bp.fields(msl[k].astype(float) + 85000.0)
        i = np.nonzero(keys == k)[0]
        a = np.rint((bp.ROW0 - lat[i]) / 0.25).astype(int); b = np.rint((lon[i] - bp.COL0) / 0.25).astype(int)
        G[i] = GG[a, b]; grad[i] = g[a, b]
        y, mo = k // 1000000, k // 10000 % 100
        win[i] = y if mo >= 6 else y - 1
    ok &= np.isfinite(G) & np.isfinite(grad)
    L = [f"ASCAT cells used {int(ok.sum())}, winters {len(np.unique(win[ok]))} ({win[ok].min()}-{win[ok].max()}); wind byte <=250, no rain flag, ocean >=1 degree from land, within 1.5 h of 00/12 UTC", ""]
    L.append("gradient bin (hPa/100nm) | n | winters | median ASCAT kt | p10-p90 | median ASCAT/G | share of cells at or above 45 kt")
    for a, b in zip([0, 2, 4, 6, 8, 10, 12], [2, 4, 6, 8, 10, 12, 99]):
        m = ok & (grad >= a) & (grad < b)
        if m.sum() < 200: L.append(f"{a}-{b}: thin"); continue
        w_ = wind[m]
        L.append(f"{a}-{b if b < 99 else '+'} | {int(m.sum())} | {len(np.unique(win[m]))} | {np.median(w_):.1f} | {np.percentile(w_,10):.1f}-{np.percentile(w_,90):.1f} | {np.median(w_/G[m]):.2f} | {100*(w_>=45).mean():.1f}%")
    open(out, "w").write("\n".join(L) + "\n"); print("\n".join(L))


if __name__ == "__main__":
    main(sys.argv[1])
