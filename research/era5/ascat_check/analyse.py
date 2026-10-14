"""Pre-registered tests P1-P4 (see PREREGISTRATION.md) on a footprint table from collocate.py.
P1 and P2 p-values are one-sided bootstrap p against the agreement threshold (rho 0.6, AUC 0.85); P3 two-sided vs 0; P4 permutation.
Season-block bootstrap (resample seasons) for intervals; BH-FDR over the 5 primary tests."""
import sys, numpy as np, pandas as pd

NB = 1000
LEVEL = 45.0
rng = np.random.default_rng(20261013)


def boot_seasons(F, stat, nb=NB):
    seas = F.season.unique()
    groups = {s: g for s, g in F.groupby("season")}
    out = []
    for _ in range(nb):
        pick = rng.choice(seas, len(seas))
        try:
            out.append(stat(pd.concat([groups[s] for s in pick])))
        except Exception:
            pass
    return np.array(out)


def p_two_sided(b):
    b = b[np.isfinite(b)]
    return float(min(1.0, 2 * min((b <= 0).mean(), (b >= 0).mean()) + 1 / (len(b) + 1)))


def rho(F):
    return F.era5_max_fp.rank().corr(F.scat_max.rank())


def auc(F):
    y = (F.scat_max >= LEVEL).astype(int)
    r = F.era5_max_fp.rank()
    n1, n0 = int(y.sum()), int((1 - y).sum())
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)      # Mann-Whitney AUC


def add_pdiff(F):
    """per HF storm: scat_max minus the mean scat_max of depth-matched sub-HF storms (same basin, minp within 5 hPa).
    Matching is fixed once on the full sample; the bootstrap then resamples seasons of these pairs."""
    lo = F[F.idx < 71.7]
    pd_ = np.full(len(F), np.nan)
    for k, (_, r) in enumerate(F.iterrows()):
        if r.idx >= 71.7:
            c = lo[(lo.basin == r.basin) & ((lo.minp - r.minp).abs() <= 5)]
            if len(c):
                pd_[k] = r.scat_max - c.scat_max.mean()
    return F.assign(pdiff=pd_)


def pdiff_mean(F):
    return float(F.pdiff.mean())


def matched_diff(F):
    """HF (index>=71.7) vs depth-matched sub-HF in the same basin, minp within 5 hPa, with replacement."""
    hf, lo = F[F.idx >= 71.7], F[F.idx < 71.7]
    d = []
    for _, r in hf.iterrows():
        c = lo[(lo.basin == r.basin) & ((lo.minp - r.minp).abs() <= 5)]
        if len(c):
            d.append(r.scat_max - c.scat_max.mean())
    return float(np.mean(d)) if d else np.nan


def hetero(F):
    G = F[F.era5_max_fp >= 54].copy()
    G["lr"] = np.log(G.scat_max / G.era5_max_fp)
    def Fstat(x):
        m = x.groupby("sensor").lr.agg(["mean", "size", "var"]).dropna()
        gm = np.average(m["mean"], weights=m["size"])
        return float((m["size"] * (m["mean"] - gm) ** 2).sum() / max(np.average(m["var"], weights=m["size"]), 1e-9))
    obs = Fstat(G)
    season_ix = [np.where(G.season.values == q)[0] for q in G.season.unique()]
    perm = []
    for _ in range(NB):
        H = G.copy()
        sv = H.sensor.values.copy()
        for ix in season_ix:
            sv[ix] = rng.permutation(sv[ix])
        H["sensor"] = sv
        perm.append(Fstat(H))
    return obs, float((np.array(perm) >= obs).mean() + 1 / (NB + 1)), G.groupby("sensor").lr.agg(["mean", "size"])


def bh(p):
    p = np.asarray(p); o = np.argsort(p); n = len(p)
    q = np.empty(n); prev = 1.0
    for k in range(n - 1, -1, -1):
        prev = min(prev, p[o[k]] * n / (k + 1)); q[o[k]] = prev
    return q


def main(path, mincov=0.40, quiet=False):
    F = pd.read_csv(path)
    F = F[F.cover >= mincov]
    # one pass per storm: the sensor with the highest coverage
    F = F.sort_values("cover", ascending=False).drop_duplicates("track")
    F["cover_ok"] = True
    rows = []
    for b in ("atl", "pac"):
        G = F[F.basin == b]
        if len(G) < 30:
            rows.append((f"P1 rho {b}", np.nan, np.nan, np.nan, np.nan, len(G))); continue
        bs = boot_seasons(G, rho)
        rows.append((f"P1 rho {b}", rho(G), *np.percentile(bs, [2.5, 97.5]), float((bs <= 0.6).mean()), len(G)))
    bs = boot_seasons(F, auc)
    rows.append(("P2 AUC (>=45 kt)", auc(F), *np.percentile(bs, [2.5, 97.5]), float((bs <= 0.85).mean()), len(F)))
    F = add_pdiff(F)
    bs = boot_seasons(F, pdiff_mean)
    rows.append(("P3 HF minus matched (kt)", pdiff_mean(F), *np.nanpercentile(bs, [2.5, 97.5]), p_two_sided(bs), int(F.pdiff.notna().sum())))
    allf = pd.read_csv(path); allf = allf[allf.cover >= mincov]
    obs, ph, tab = hetero(allf)
    rows.append(("P4 sensor heterogeneity F", obs, np.nan, np.nan, ph, len(allf)))
    R = pd.DataFrame(rows, columns=["test", "estimate", "lo95", "hi95", "p", "n"])
    R["q_BH"] = bh(R.p.fillna(1.0).values)
    if not quiet:
        print(R.to_string(index=False)); print(tab)
    return R, F


if __name__ == "__main__":
    main(sys.argv[1])
