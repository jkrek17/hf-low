"""Composites, season-block bootstrap, descriptive scalars (PREREGISTRATION.md).

usage: ERA5_WORK=DIR analyse.py [NBOOT]
Reads DIR/boxes/*.npz and results/storm_times.csv; writes results/comp_rot.npz, comp_north.npz,
results/scalars.csv, results/summary.txt.
"""
import os, sys, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from common import WORK, BX, BY, BR
from extract import FIELDS, BASE

HERE = os.path.dirname(os.path.abspath(__file__))
ANCH = ["m48", "m24", "m12", "onset", "peak"]
SEED = 20261009
NF = len(FIELDS)


def load(frame):
    D = pd.read_csv(os.path.join(HERE, "results", "storm_times.csv"), dtype={"time": str}).set_index("key")
    boxes = {}
    for f in sorted(glob.glob(os.path.join(WORK, "boxes", "*.npz"))):
        z = np.load(f)
        keys, arr = z["keys"], z[frame]
        for i, k in enumerate(keys):
            boxes[k] = arr[i]
    missing = set(D.index) - set(boxes)
    assert not missing, f"{len(missing)} rows without boxes"
    return D, boxes


def bh(p, q=0.05):
    ok = np.isfinite(p)
    out = np.zeros(p.shape, bool)
    pv = p[ok]
    if pv.size == 0:
        return out
    o = np.argsort(pv)
    thr = q * (np.arange(1, pv.size + 1) / pv.size)
    below = pv[o] <= thr
    k = np.max(np.where(below)[0]) + 1 if below.any() else 0
    sel = np.zeros(pv.size, bool)
    sel[o[:k]] = True
    out[ok] = sel
    return out


def season_sums(X, seasons):
    """Per-season sums and counts of valid values: X [n, F, ny, nx] -> sums [S, F*ny*nx], counts, season ids."""
    ids = np.unique(seasons)
    flat = X.reshape(len(X), -1)
    ok = np.isfinite(flat)
    S = np.stack([np.where(ok[seasons == s], flat[seasons == s], 0).sum(0) for s in ids])
    C = np.stack([ok[seasons == s].sum(0) for s in ids])
    return S, C, ids


def boot_p(S, C, nboot, rng):
    """Two-sided pixelwise p of mean = 0 by resampling whole seasons."""
    ns = len(S)
    W = rng.multinomial(ns, np.ones(ns) / ns, size=nboot).astype(np.float32)
    num, den = W @ S.astype(np.float32), W @ C.astype(np.float32)
    with np.errstate(invalid="ignore", divide="ignore"):
        m = num / den
    pl, pg = np.nanmean(m <= 0, 0), np.nanmean(m >= 0, 0)
    p = np.minimum(1, 2 * np.minimum(pl, pg))
    p[den.min(0) == 0] = np.nan
    return p


def main(nboot):
    rng = np.random.default_rng(SEED)
    ny = nx = BX.shape[0]
    lines, rows = [], []
    out = {"rot": {}, "north": {}}
    for fr in ("rot", "north"):
        D, boxes = load(fr)
        for basin in ("atl", "pac"):
            for a in ANCH:
                keys = D.index[(D.basin == basin) & (D.anchor == a)]
                X = np.stack([boxes[k] for k in keys])
                good = np.isfinite(X[:, FIELDS.index("ws250"), ny // 2, nx // 2])
                X, kk = X[good], keys[good]
                seas = D.loc[kk, "season"].values
                mean = np.nanmean(X, 0)
                nvalid = int(len(X))
                tag = f"{basin}_{a}"
                out[fr][f"{tag}_mean"] = mean.astype(np.float32)
                out[fr][f"{tag}_n"] = np.array(nvalid)
                S, C, ids = season_sums(X[:, NF // 2:], seas)
                p = boot_p(S, C, nboot, rng).reshape(NF // 2, ny, nx)
                rej = np.stack([bh(p[j]) for j in range(NF // 2)])
                out[fr][f"{tag}_rej"] = rej
                lines.append(f"{fr} {basin} {a}: n={nvalid} storms, {len(ids)} seasons; pixels passing BH q<0.05: " +
                             ", ".join(f"{BASE[j]} {int(rej[j].sum())}" for j in range(NF // 2)))
                if fr == "rot":
                    for f in ("ws250", "ws500"):
                        rows.append(dict(frame=fr, basin=basin, anchor=a, scalar=f"S1 centre {f} (m/s)", n=nvalid,
                                         value=float(mean[FIELDS.index(f), ny // 2, nx // 2])))
                    for f, lim in (("ws250", 2500), ("ws500", 2500)):
                        m = np.where(BR <= lim, mean[FIELDS.index(f)], np.nan)
                        i, j = np.unravel_index(np.nanargmax(m), m.shape)
                        rows.append(dict(frame=fr, basin=basin, anchor=a, scalar=f"S2 max mean {f} within {lim} km (m/s)",
                                         n=nvalid, value=float(m[i, j]), x_km=float(BX[i, j]), y_km=float(BY[i, j])))
                m = np.where(BR <= 3000, mean[FIELDS.index("a_z500")], np.nan)
                i, j = np.unravel_index(np.nanargmin(m), m.shape)
                rows.append(dict(frame=fr, basin=basin, anchor=a, scalar="S4 min mean Z500 anomaly within 3000 km (m)",
                                 n=nvalid, value=float(m[i, j]), x_km=float(BX[i, j]), y_km=float(BY[i, j])))
                m = np.where(BR <= 2000, mean[FIELDS.index("d250")], np.nan)
                i, j = np.unravel_index(np.nanargmax(m), m.shape)
                rows.append(dict(frame=fr, basin=basin, anchor=a, scalar="S4 max mean 250 hPa divergence within 2000 km (1e-5/s)",
                                 n=nvalid, value=float(m[i, j]), x_km=float(BX[i, j]), y_km=float(BY[i, j])))
                # S3: per-storm location of the maximum raw 250 hPa speed within 2000 km
                ws = X[:, FIELDS.index("ws250")].copy()
                ws[:, BR > 2000] = np.nan
                flat = np.where(np.isfinite(ws), ws, -1).reshape(len(X), -1)
                am = flat.argmax(1)
                sx, sy = BX.ravel()[am], BY.ravel()[am]
                ahead = sx > 0
                ids2 = np.unique(seas)
                bs = []
                for _ in range(nboot):
                    pick = rng.choice(ids2, len(ids2))
                    sel = np.concatenate([np.where(seas == s)[0] for s in pick])
                    bs.append(ahead[sel].mean())
                lo, hi = np.percentile(bs, [5, 95])
                rows.append(dict(frame=fr, basin=basin, anchor=a, scalar="S3 share of storms whose max 250 hPa speed within 2000 km is ahead (x>0)",
                                 n=nvalid, value=float(ahead.mean()), lo90=float(lo), hi90=float(hi)))
                rows.append(dict(frame=fr, basin=basin, anchor=a, scalar="S3 median distance of that maximum (km)",
                                 n=nvalid, value=float(np.median(np.hypot(sx, sy)))))
    for fr in out:
        np.savez_compressed(os.path.join(HERE, "results", f"comp_{fr}.npz"), fields=np.array(FIELDS), **out[fr])
    pd.DataFrame(rows).to_csv(os.path.join(HERE, "results", "scalars.csv"), index=False)
    open(os.path.join(HERE, "results", "summary.txt"), "w").write(
        f"nboot {nboot}, seed {SEED}\n" + "\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 2000)
