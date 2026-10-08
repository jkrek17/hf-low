"""RA-16: SST-front strength vs HF-centre hours in 40-45N 160-170E (pipeline A proxy). Plan: PREREGISTRATION.md.
usage: python3 -I analysis.py FRONT_MONTHLY_CSV OUTDIR      (reads committed files only besides FRONT_MONTHLY_CSV)
"""
import sys, os, json
import numpy as np, pandas as pd
from scipy import stats
import statsmodels.api as sm

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
HFH = os.path.join(ROOT, "research", "era5", "hf_history", "results")
THR = 71.7
BOX = (40, 45, 160, 170)
NPERM, NBOOT = 50000, 10000
rng = np.random.default_rng(16)


def inbox(lat, lon):
    lon = np.asarray(lon, float) % 360
    return (lat >= BOX[0]) & (lat < BOX[1]) & (lon >= BOX[2]) & (lon < BOX[3])


def proxy_outcomes():
    C = pd.read_csv(os.path.join(HFH, "era5_hf_catalog.csv"))
    E = C[(C.role == "event") & (C.basin == "pac") & (C.season >= 2004) & (C.season <= 2025)].set_index("track")
    T = pd.read_csv(os.path.join(HFH, "era5_hf_catalog_tracks.csv"))
    T = T[T.track.isin(E.index)].copy()
    T["t"] = pd.to_datetime(T.time.astype(str), format="%Y%m%d%H")
    T["indom"] = T.basin.fillna("") != ""
    rows = []
    for tr, g in T.groupby("track"):
        g = g.sort_values("t").reset_index(drop=True)
        tmin = g.t[g.msl.idxmin()]
        hf = g[g.indom & (g.g800 >= THR)]
        hb = hf[inbox(hf.lat.values, hf.lon.values)]
        if hb.empty:
            continue
        rows.append(dict(track=tr, season=int(E.season[tr]), n=len(hb), n_dep=int((hb.t <= tmin).sum()),
                         n_mat=int((hb.t > tmin).sum())))
    D = pd.DataFrame(rows)
    S = D.groupby("season").agg(fix=("n", "sum"), dep=("n_dep", "sum"), mat=("n_mat", "sum"), C=("track", "nunique"))
    S = S.reindex(range(2004, 2026), fill_value=0)
    big = D.groupby("season").n.max() / D.groupby("season").n.sum()
    return S, D, big


def archive_outcomes():
    d = json.load(open(os.path.join(ROOT, "docs", "data", "hf-lows.json")))
    cnt = {s: 0 for s in range(2004, 2026)}
    for raw in d["lows"]:
        l = dict(zip(d["lowFields"], raw))
        if l["basin"] != "pac" or l["season"] < 2004 or l["season"] > 2025:
            continue
        f = pd.DataFrame(l["fixes"], columns=["date", "lat", "lon", "cat", "pres"])
        f["t"] = pd.to_datetime(f.date.astype("int64").astype(str), format="%Y%m%d%H")
        f = f.drop_duplicates("t")
        h = f[f.cat == "HF"]
        if len(h):
            cnt[l["season"]] += int(inbox(h.lat.values.astype(float), h.lon.values.astype(float)).sum())
    return pd.Series(cnt)


def predictors(path):
    M = pd.read_csv(path)
    out = {}
    for s in range(2004, 2026):
        son = M[(M.year == s) & M.month.isin([9, 10, 11])]
        djf = M[((M.year == s) & (M.month == 12)) | ((M.year == s + 1) & M.month.isin([1, 2]))]
        for k in ("front", "front_ke", "front_oe"):
            out.setdefault(f"{k}_son", {})[s] = son[k].mean()
            out.setdefault(f"{k}_djf", {})[s] = djf[k].mean()
    return pd.DataFrame(out)


def resid(v, Z):
    if Z is None:
        return v - v.mean()
    Z = np.column_stack([np.ones(len(v)), Z])
    b, *_ = np.linalg.lstsq(Z, v, rcond=None)
    return v - Z @ b


def spearman_perm(x, y, Z=None, nperm=NPERM):
    """Spearman (rank) correlation after removing covariates Z from the ranks; permutation p (two-sided)."""
    rx = resid(stats.rankdata(x).astype(float), Z); ry = resid(stats.rankdata(y).astype(float), Z)
    r = np.corrcoef(rx, ry)[0, 1]
    idx = np.argsort(rng.random((nperm, len(rx))), axis=1)
    rp = np.array([np.corrcoef(rx[i], ry)[0, 1] for i in idx])
    p = (1 + (np.abs(rp) >= abs(r) - 1e-12).sum()) / (nperm + 1)
    return float(r), float(p)


def fisher_ci(r, n, k=0):
    z = np.arctanh(r); se = 1 / np.sqrt(n - 3 - k)
    return float(np.tanh(z - 1.96 * se)), float(np.tanh(z + 1.96 * se))


def mde(n, alpha=0.05, power=0.8):
    return float(np.tanh((stats.norm.ppf(1 - alpha / 2) + stats.norm.ppf(power)) / np.sqrt(n - 3)))


def qp_slope(x, y):
    X = sm.add_constant((x - x.mean()) / x.std(ddof=1))
    m = sm.GLM(y, X, family=sm.families.Poisson()).fit(cov_type="HC1")
    b, se = m.params[1], m.bse[1]
    return float(np.exp(b)), float(np.exp(b - 1.96 * se)), float(np.exp(b + 1.96 * se))


def bh(p):
    p = np.asarray(p); n = len(p); o = np.argsort(p)
    q = np.empty(n); prev = 1.0
    for rank, i in zip(range(n, 0, -1), o[::-1]):
        prev = min(prev, p[i] * n / rank); q[i] = prev
    return q


def lag1(v):
    v = np.asarray(v, float); v = v - v.mean()
    return float((v[1:] * v[:-1]).sum() / (v * v).sum())


def main():
    out = sys.argv[2]; os.makedirs(out, exist_ok=True)
    S, D, big = proxy_outcomes()
    A = archive_outcomes()
    P = predictors(sys.argv[1])
    S["hours"] = S.fix * 6.0; S["dep_h"] = S.dep * 6.0; S["mat_h"] = S.mat * 6.0
    T = pd.concat([S, P, A.rename("arch_fix")], axis=1)
    T.to_csv(os.path.join(out, "season_table.csv"), float_format="%.4f")
    s07 = T.loc[2007:2025]; sall = T
    L = []; w = L.append
    w("RA-16: SST-front strength vs HF-centre hours, box 40-45N 160-170E. Pipeline A (ERA5 proxy). Plan: PREREGISTRATION.md.")
    w(f"Seasons: primary 2007-2025 (n={len(s07)}); sensitivity 2004-2025 (n={len(sall)}).")
    w(f"Box HF-centre hours per season, 2007-2025: mean {s07.hours.mean():.1f}, SD {s07.hours.std():.1f}, min {s07.hours.min():.0f}, max {s07.hours.max():.0f}; seasons with zero {int((s07.hours==0).sum())}.")
    w(f"Storms per season in box: mean {s07.C.mean():.1f}. Largest single-storm share of a season's box fixes: median {big.loc[2007:2025].median():.2f}, max {big.loc[2007:2025].max():.2f}.")
    w(f"Lag-1 autocorrelation 2007-2025: FRONT_SON {lag1(s07.front_son):.2f}, hours {lag1(s07.hours):.2f}. Corr(FRONT_SON, FRONT_DJF) = {np.corrcoef(s07.front_son, s07.front_djf)[0,1]:.2f}.")
    w(f"Analytic MDE |rho| (80% power, two-sided 0.05): n=19 {mde(19):.2f}, n=22 {mde(22):.2f}; with alpha 0.05/11: n=19 {mde(19, 0.05/11):.2f}.")
    tests = []

    def run(tid, desc, x, y, Z=None, n=None, k=0, extra=None):
        x = np.asarray(x, float); y = np.asarray(y, float)
        r, p = spearman_perm(x, y, Z)
        lo, hi = fisher_ci(r, len(x), k)
        row = dict(id=tid, test=desc, n=len(x), rho=r, lo=lo, hi=hi, p=p, pearson=float(np.corrcoef(x, y)[0, 1]))
        if extra:
            row.update(extra)
        tests.append(row)

    f07 = s07.front_son.values
    sl = qp_slope(f07, s07.fix.values)
    run("P1", "FRONT(SON) vs H, 2007-25", f07, s07.hours.values, extra=dict(qp_ratio=sl[0], qp_lo=sl[1], qp_hi=sl[2]))
    era = (sall.index < 2007).astype(float)
    run("S1", "FRONT(SON) vs H, 2004-25, product dummy removed", sall.front_son.values, sall.hours.values, Z=era, k=1)
    run("S2", "FRONT_KE(SON) vs H", s07.front_ke_son.values, s07.hours.values)
    run("S3", "FRONT_OE(SON) vs H", s07.front_oe_son.values, s07.hours.values)
    run("S4", "FRONT(DJF) concurrent vs H", s07.front_djf.values, s07.hours.values)
    run("S5", "FRONT(SON) vs H deepening-phase", f07, s07.dep_h.values)
    run("S6", "FRONT(SON) vs H mature/decay-phase", f07, s07.mat_h.values)
    run("S7", "FRONT(SON) vs storm count", f07, s07.C.values)
    run("S8", "FRONT(SON) vs archive box fixes, 2013-14 dummy removed", f07, s07.arch_fix.values, Z=(s07.index >= 2013).astype(float), k=1)
    run("S9", "FRONT(SON) vs H, year removed", f07, s07.hours.values, Z=s07.index.values.astype(float), k=1)
    # S10 bootstrap of rho_dep - rho_mat
    n = len(s07); d = []
    for _ in range(NBOOT):
        i = rng.integers(0, n, n)
        if np.ptp(f07[i]) == 0:
            continue
        a = stats.spearmanr(f07[i], s07.dep_h.values[i])[0]; b = stats.spearmanr(f07[i], s07.mat_h.values[i])[0]
        if np.isfinite(a) and np.isfinite(b):
            d.append(a - b)
    d = np.array(d)
    r5 = [t for t in tests if t["id"] == "S5"][0]["rho"]; r6 = [t for t in tests if t["id"] == "S6"][0]["rho"]
    p10 = float(min(1.0, 2 * min((d <= 0).mean(), (d >= 0).mean()) + 1 / len(d)))
    tests.append(dict(id="S10", test="rho(S5) - rho(S6), season bootstrap", n=n, rho=r5 - r6, lo=float(np.percentile(d, 2.5)),
                      hi=float(np.percentile(d, 97.5)), p=p10, pearson=np.nan))
    R = pd.DataFrame(tests)
    R["q_bh"] = bh(R.p.values)
    R.to_csv(os.path.join(out, "results.csv"), index=False, float_format="%.4f")
    w("\nid   n   rho    95% (Fisher / bootstrap for S10)   p(perm)  q(BH,11)  test")
    for _, t in R.iterrows():
        w(f"{t.id:4s} {int(t.n):2d} {t.rho:+.3f} [{t.lo:+.2f}, {t.hi:+.2f}]  {t.p:.4f}  {t.q_bh:.3f}  {t.test}")
    w(f"\nTests with q < 0.10: {int((R.q_bh < 0.10).sum())} of {len(R)}; with raw p < 0.05: {int((R.p < 0.05).sum())}.")
    p1 = R[R.id == "P1"].iloc[0]
    w(f"P1 effect size: quasi-Poisson ratio of box fixes per SD of FRONT(SON) {p1.qp_ratio:.2f} (95% robust {p1.qp_lo:.2f}-{p1.qp_hi:.2f}).")
    open(os.path.join(out, "summary.txt"), "w").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
