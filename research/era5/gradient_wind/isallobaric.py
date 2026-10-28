"""Isallobaric wind test. Plan: PREREGISTRATION_ISALLOBARIC.md. Usage: python3 isallobaric.py PAIRS.npz OUTFILE (reads shared caches; no download)"""
import glob, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_pairs as bp

SH, KT, RHO = bp.SH, bp.KT, bp.RHO
B = 2000
rng = np.random.default_rng(20261028)
BLAT, BLON = np.arange(29.5, 41, 1.0), np.arange(280.5, 305, 1.0)


def vecs(p, ptend):
    la = bp.ROW0 - 0.25 * np.arange(p.shape[0]); f = (2 * 7.292e-5 * np.sin(np.radians(la)))[:, None]
    px, py, *_ = bp.derivs(p, la)
    tx, ty, *_ = bp.derivs(bp.box(ptend, 5), la)
    return -py / (RHO * f), px / (RHO * f), -tx / (RHO * f ** 2), -ty / (RHO * f ** 2), bp.box(ptend, 5)


def build():
    cells = [(a, b) for a in bp.LATS for b in bp.LONS] + [(a, b) for a in BLAT for b in BLON]
    r = np.rint((bp.ROW0 - np.array([c[0] for c in cells])) / 0.25).astype(int); c = np.rint((np.array([c[1] for c in cells]) - bp.COL0) / 0.25).astype(int)
    rows = []
    for fn in sorted(glob.glob(f"{SH}/gulf-stream/msl/msl_*.npz")):
        z = np.load(fn); t = z["t"]; a = z["msl"].astype(np.float64) + 85000.0
        h = np.array([bp.key_to_h1900(int(x)) for x in t])
        for i in range(1, len(t) - 1):
            if h[i + 1] - h[i] != 12 or h[i] - h[i - 1] != 12: continue                    # needs both neighbours 12 h away
            ptend = (a[i + 1] - a[i - 1]) / (24 * 3600.0)
            vgu, vgv, vau, vav, tend = vecs(a[i], ptend)
            rows.append(pd.DataFrame(dict(key=int(t[i]), lat=[x[0] for x in cells], lon=[x[1] for x in cells],
                                          vgu=vgu[r, c], vgv=vgv[r, c], vau=vau[r, c], vav=vav[r, c], tend=tend[r, c] * 43200.0 / 100.0)))
        print(os.path.basename(fn), flush=True)
    return pd.concat(rows, ignore_index=True)


def ranks(a): return np.argsort(np.argsort(a)).astype(float)


def main(pairs, out):
    iso = build(); iso["Vg"] = np.hypot(iso.vgu, iso.vgv) * KT; iso["Va"] = np.hypot(iso.vau, iso.vav) * KT
    iso["S"] = np.hypot(iso.vgu + iso.vau, iso.vgv + iso.vav) * KT
    iso["along"] = (iso.vau * iso.vgu + iso.vav * iso.vgv) / np.maximum(np.hypot(iso.vgu, iso.vgv), 1e-9) * KT
    d = pd.DataFrame({k: v for k, v in dict(np.load(pairs)).items()}); d["key"] = d.key.astype(np.int64)
    m = d.merge(iso, on=["key", "lat", "lon"], how="inner"); m = m[np.isfinite(m.S) & np.isfinite(m.gust)]
    L = [f"I1/I3 sample: {len(m)} cell-times with tendency, {m.winter.nunique()} winters"]
    # descriptive
    for a_, b_ in zip([0, 2, 4, 6, 8, 10, 12], [2, 4, 6, 8, 10, 12, 99]):
        s = m[(m.grad >= a_) & (m.grad < b_)]
        if len(s) > 200:
            L.append(f"  grad {a_}-{b_ if b_ < 99 else '+'}: |Va| median {s.Va.median():.2f} kt, p90 {s.Va.quantile(.9):.2f}, >5 kt {100*(s.Va>5).mean():.1f}%, along-geostrophic component median {s.along.median():+.2f} kt (p10 {s.along.quantile(.1):+.2f}, p90 {s.along.quantile(.9):+.2f}), n {len(s)}")
    # I1
    s = m[m.grad >= 4].reset_index(drop=True); rg, rVg, rS = ranks(s.gust.values), ranks(s.Vg.values), ranks(s.S.values)
    cor = lambda a, b, i: np.corrcoef(a[i], b[i])[0, 1]
    base = cor(rg, rS, slice(None)) - cor(rg, rVg, slice(None)); w = s.winter.values; u = np.unique(w); idx = [np.nonzero(w == x)[0] for x in u]; bs = []
    for _ in range(B):
        i = np.concatenate([idx[j] for j in rng.integers(0, len(u), len(u))]); bs.append(cor(rg, rS, i) - cor(rg, rVg, i))
    bs = np.array(bs); ci1 = np.percentile(bs, [2.5, 97.5]); p1 = float(2 * min((bs <= 0).mean(), (bs >= 0).mean()))
    L.append(f"I1 Spearman(gust, S) minus Spearman(gust, |Vg|), grad>=4: {base:+.4f} (95% {ci1[0]:+.4f} to {ci1[1]:+.4f}), n {len(s)}, {len(u)} winters, p {p1:.4f}")
    # I3
    g = m[(m.grad >= 6) & (m.grad < 10)].copy(); g["r"] = g.gust / g.G
    wm = lambda x: {int(k): float(np.median(v.r)) for k, v in x.groupby("winter") if len(v) >= 20}
    fa, ri = wm(g[g.tend <= -1]), wm(g[g.tend >= 1]); ws = sorted(set(fa) & set(ri)); df = np.array([fa[k] - ri[k] for k in ws])
    ix = rng.integers(0, len(df), (B, len(df))); b3 = df[ix].mean(1); ci3 = np.percentile(b3, [2.5, 97.5]); p3 = float(2 * min((b3 <= 0).mean(), (b3 >= 0).mean()))
    L.append(f"I3 falling (<= -1 hPa/12 h) minus rising (>= +1), median gust/G at grad 6-10: {df.mean():+.3f} (95% {ci3[0]:+.3f} to {ci3[1]:+.3f}), {len(ws)} winters, p {p3:.4f}; n falling {int((g.tend<=-1).sum())}, rising {int((g.tend>=1).sum())}")
    # I2
    a = pd.read_csv(f"{SH}/gulf-stream/tables/boxes_B.csv.gz"); a = a[a.winter <= 2024].dropna(subset=["v_ascat", "v_gfs", "dT"]).copy()
    hh = np.rint(a.tmin.values / 60.0 / 12.0) * 12.0
    a["key"] = pd.to_datetime((hh * 3600).astype("int64"), unit="s").strftime("%Y%m%d%H").astype(np.int64)
    a = a.merge(iso[["key", "lat", "lon", "Vg", "S"]], on=["key", "lat", "lon"], how="inner").reset_index(drop=True)
    def Xm(dd, c): return np.column_stack([np.ones(len(dd)), dd.v_gfs.values, dd[c].values])
    def lowo(c):
        p = np.empty(len(a))
        for w_ in a.winter.unique():
            tr, te = (a.winter != w_).values, (a.winter == w_).values
            p[te] = Xm(a[te], c) @ np.linalg.lstsq(Xm(a[tr], c), a.v_ascat.values[tr], rcond=None)[0]
        return p
    pV, pS = lowo("Vg"), lowo("S"); y = a.v_ascat.values
    days = a.day.values; ud, inv = np.unique(days, return_inverse=True)
    def sk(msk):
        e0, e1 = (y - pV) ** 2, (y - pS) ** 2; s0 = np.bincount(inv, e0 * msk, len(ud)); s1 = np.bincount(inv, e1 * msk, len(ud))
        ixx = rng.integers(0, len(ud), (B, len(ud))); b = 1 - s1[ixx].sum(1) / s0[ixx].sum(1)
        per = {int(w_): 1 - e1[msk & (a.winter.values == w_)].mean() / e0[msk & (a.winter.values == w_)].mean() for w_ in a.winter.unique()}
        return 1 - e1[msk].mean() / e0[msk].mean(), np.percentile(b, [2.5, 97.5]), per, b, np.sqrt(e0[msk].mean()), np.sqrt(e1[msk].mean())
    for lab, msk in [("I2 all", np.ones(len(a), bool)), ("I2 GFS wind >= 25 kt (secondary)", (a.v_gfs >= 25).values)]:
        s_, ci_, per, b_, r0, r1 = sk(msk); L.append(f"{lab}: adding S instead of |Vg| to GFS wind: skill {100*s_:+.2f}% (95% {100*ci_[0]:+.2f} to {100*ci_[1]:+.2f}), RMSE {r0:.3f} -> {r1:.3f} kt, n {int(msk.sum())}, by winter " + ", ".join(f"{k}:{100*v:+.1f}%" for k, v in per.items()))
        if lab == "I2 all": b2 = b_
    p2 = float(2 * min((b2 <= 0).mean(), (b2 >= 0).mean()))
    ps = np.maximum([p1, p2, p3], 1 / B); o = np.argsort(ps); q = np.empty(3); prev = 1.0
    for rk, i in zip(range(3, 0, -1), o[::-1]): prev = min(prev, ps[i] * 3 / rk); q[i] = prev
    L.append(f"BH over I1, I2, I3: p {ps[0]:.4f}, {ps[1]:.4f}, {ps[2]:.4f}; q {q[0]:.4f}, {q[1]:.4f}, {q[2]:.4f}")
    open(out, "w").write("\n".join(L) + "\n"); print("\n".join(L))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
