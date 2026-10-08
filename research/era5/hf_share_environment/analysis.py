"""Which environmental ingredients track the midwinter rise in the HF share? Plan: PREREGISTRATION.md.

Pipeline A tracks (ERA5 proxy); near-storm environment from research/era5/intensity/results/env_2004.csv.gz.

    python3 -I analysis.py <repo root> <out dir> [nboot]
"""
import os, sys, json
import numpy as np, pandas as pd
from multiprocessing import Pool

ROOT, OUT = sys.argv[1], sys.argv[2]
B = int(sys.argv[3]) if len(sys.argv) > 3 else 2000
SEED, THR, DEEP = 20261010, 71.7, -3.6
PRIMARY = ["jet250", "eady", "sstgrad", "sst_t500", "flux", "tcwv", "div300", "vadv500"]
SECOND = ["lat", "B", "VTL", "VTU", "sst"]
PER = {10: 0, 11: 0, 12: 1, 1: 1, 2: 1, 3: 2, 4: 2}          # 0 ON, 1 DJF, 2 MA
os.makedirs(OUT, exist_ok=True)


def load():
    t = pd.read_csv(f"{ROOT}/research/era5/hf_history/results/all_tracks.csv.gz")
    t = t[(t.season >= 2004) & (t.season <= 2025)].copy()
    t["gm"] = (t.start // 10000) % 100
    t = t[t.gm.isin(PER)].copy()
    t["per"] = t.gm.map(PER)
    t["hf"] = (t.gust800_kt >= THR).astype(int)
    f = pd.read_csv(f"{ROOT}/research/era5/intensity/results/fixes_2004.csv.gz")
    e = pd.read_csv(f"{ROOT}/research/era5/intensity/results/env_2004.csv.gz")
    f = f.merge(e, on=["track", "time"], how="left", validate="1:1")
    assert len(f) == 159430 and f.jet250.notna().all()
    f["fm"] = (f.time // 10000) % 100
    return t, f


def ref_tables(t, f):
    """one row per cyclone at three reference fixes"""
    ff = f[f.track.isin(t.track)].sort_values(["track", "time"])
    cols = PRIMARY + SECOND
    first = ff.groupby("track").head(1).set_index("track")[cols]
    deep1 = ff[ff.dp12 <= DEEP].groupby("track").head(1).set_index("track")[cols]
    mx = ff[ff.ndr24.notna()]
    mx = mx.loc[mx.groupby("track").ndr24.idxmax()].set_index("track")[cols]
    base = t.set_index("track")[["basin", "season", "gm", "per", "hf"]]
    base["deep"] = base.index.isin(deep1.index).astype(int)
    out = {}
    for name, src in (("primary", deep1), ("entry", first), ("maxdeep", mx)):
        have = base.join(src, how="inner")
        rest = base.loc[~base.index.isin(src.index)].join(first, how="inner")   # no such fix: first fix
        out[name] = pd.concat([have, rest]).reset_index().reset_index(drop=True)
    out["n_cyclones"] = len(t)
    out["n_with_fix"] = {k: len(v) for k, v in out.items() if isinstance(v, pd.DataFrame)}
    return out


def fix_table(t, f):
    ff = f[f.track.isin(t.track) & f.fm.isin(PER) & ~f.hf_now].copy()
    ff["per"] = ff.fm.map(PER)
    ff["hf"] = ff.hf24.astype(int)
    return ff[["track", "basin", "season", "per", "hf"] + PRIMARY + SECOND].reset_index(drop=True)


def fit(X, y, it=30):
    b = np.zeros(X.shape[1]); b[0] = np.log(y.mean() / (1 - y.mean()))
    for _ in range(it):
        p = 1 / (1 + np.exp(-X @ b)); w = p * (1 - p)
        s = np.linalg.solve(X.T @ (X * w[:, None]) + 1e-9 * np.eye(len(b)), X.T @ (y - p))
        b += s
        if np.abs(s).max() < 1e-8:
            break
    return b


def design(per, extra=None):
    cols = [np.ones(len(per)), (per == 1).astype(float), (per == 2).astype(float)]
    if extra is not None:
        cols += [extra[:, j] for j in range(extra.shape[1])]
    return np.column_stack(cols)


def contrasts(b):                                    # rise = DJF vs ON, fall = DJF vs MA
    return b[1], b[1] - b[2]


def compute(per, y, Z, names, deep=None, inter=False):
    r = {}
    n = len(y)
    r["n"] = n; r["hf_total"] = float(y.sum())
    for p, lab in ((0, "ON"), (1, "DJF"), (2, "MA")):
        m = per == p
        r[f"share_{lab}"] = y[m].mean() if m.any() else np.nan
    allok = ~np.isnan(Z).any(1)
    for j, x in enumerate(names):
        ok = ~np.isnan(Z[:, j])
        pp, yy, xx = per[ok], y[ok], Z[ok, j:j + 1]
        b0 = fit(design(pp), yy); b1 = fit(design(pp, xx), yy)
        c0r, c0f = contrasts(b0); c1r, c1f = contrasts(b1)
        z = xx[:, 0]
        r[f"delta_rise|{x}"] = z[pp == 1].mean() - z[pp == 0].mean()
        r[f"delta_fall|{x}"] = z[pp == 1].mean() - z[pp == 2].mean()
        r[f"beta|{x}"] = b1[3]
        r[f"c0_rise|{x}"] = c0r; r[f"c1_rise|{x}"] = c1r
        r[f"F_rise|{x}"] = 1 - c1r / c0r
        r[f"F_fall|{x}"] = 1 - c1f / c0f
        if inter:
            m2 = pp < 2
            Xi = np.column_stack([np.ones(m2.sum()), (pp[m2] == 1).astype(float), z[m2], z[m2] * (pp[m2] == 1)])
            r[f"inter|{x}"] = fit(Xi, yy[m2])[3]
    if allok.sum() > 100 and len(names) > 1:
        pp, yy, ZZ = per[allok], y[allok], Z[allok]
        b0 = fit(design(pp), yy); bj = fit(design(pp, ZZ), yy)
        c0r, c0f = contrasts(b0); c1r, c1f = contrasts(bj)
        r["F_rise|JOINT"] = 1 - c1r / c0r; r["F_fall|JOINT"] = 1 - c1f / c0f
        r["c0_rise|JOINT"] = c0r; r["n|JOINT"] = len(yy)
        for j, x in enumerate(names):
            keep = [k for k in range(len(names)) if k != j]
            bd = fit(design(pp, ZZ[:, keep]), yy)
            r[f"Funique_rise|{x}"] = r["F_rise|JOINT"] - (1 - contrasts(bd)[0] / c0r)
    if deep is not None:
        b = fit(design(per), deep.astype(float))
        r["deep_rise"], r["deep_fall"] = contrasts(b)
        m = deep == 1
        b = fit(design(per[m]), y[m])
        r["hfgd_rise"], r["hfgd_fall"] = contrasts(b)
        b = fit(design(per), y.astype(float))
        r["share_rise"], r["share_fall"] = contrasts(b)
        for p, lab in ((0, "ON"), (1, "DJF"), (2, "MA")):
            r[f"pdeep_{lab}"] = deep[per == p].mean()
            r[f"hfgd_{lab}"] = y[(per == p) & m].mean()
    return r


def run_job(job):
    label, basin, df, names, deep_flag, inter, nb, seed = job
    d = df[df.basin == basin].reset_index(drop=True)
    Zdf = d[names].astype(float)
    Z = ((Zdf - Zdf.mean()) / Zdf.std(ddof=0)).values
    per, y, sea = d.per.values, d.hf.values.astype(float), d.season.values
    deep = d.deep.values if deep_flag and "deep" in d else None
    pt = compute(per, y, Z, names, deep, inter)
    idx = {s: np.where(sea == s)[0] for s in np.unique(sea)}
    ss = np.array(list(idx))
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(nb):
        pick = rng.choice(ss, len(ss))
        ii = np.concatenate([idx[s] for s in pick])
        boots.append(compute(per[ii], y[ii], Z[ii], names, deep[ii] if deep is not None else None, inter))
    bt = pd.DataFrame(boots)
    # monthly profile (descriptive): mean z by genesis month, season-block interval
    return label, basin, pt, bt, len(d)


def main():
    t, f = load()
    R = ref_tables(t, f)
    fx = fix_table(t, f)
    print("cyclones in population", R["n_cyclones"], "with a 00/12 fix", R["n_with_fix"], "fix-level rows", len(fx), flush=True)
    jobs = []
    allcols = PRIMARY + SECOND
    pri = R["primary"]
    ok8 = pri[PRIMARY].notna().all(1)
    for basin in ("atl", "pac"):
        jobs += [
            ("primary", basin, pri, PRIMARY, True, True, B, SEED),
            ("S5_secondary", basin, pri, SECOND, False, False, B, SEED),
            ("S1_entry", basin, R["entry"], PRIMARY, False, False, B, SEED),
            ("S2_maxdeep", basin, R["maxdeep"], PRIMARY, False, False, B, SEED),
            ("S6a_commoncc", basin, pri[ok8], PRIMARY, False, False, B, SEED),
            ("S6b_deepening_only", basin, pri[pri.deep == 1], PRIMARY, False, False, B, SEED),
            ("S3_fixlevel", basin, fx, PRIMARY, False, False, min(B, 500), SEED),
        ]
    with Pool(4) as P:
        res = P.map(run_job, jobs, chunksize=1)
    rows, pts = [], {}
    for label, basin, pt, bt, n in res:
        pts[f"{label}|{basin}"] = {k: float(v) for k, v in pt.items()}
        bt.to_csv(f"{OUT}/boot_{label}_{basin}.csv.gz", index=False, float_format="%.5g") if label == "primary" else None
        for k, v in pt.items():
            lo, hi = np.nanpercentile(bt[k], [2.5, 97.5]) if k in bt else (np.nan, np.nan)
            if k in bt and k.split("|")[0] in ("delta_rise", "delta_fall", "beta", "inter"):
                nn = bt[k].dropna()
                p = min(1.0, 2 * (1 + min((nn >= 0).sum(), (nn <= 0).sum())) / (len(nn) + 1))
            else:
                p = np.nan
            q, x = (k.split("|") + [""])[:2]
            rows.append(dict(variant=label, basin=basin, quantity=q, ingredient=x, estimate=v, lo=lo, hi=hi, p=p, nboot=len(bt)))
    T = pd.DataFrame(rows)
    # Benjamini-Hochberg within a family
    T["q"] = np.nan
    fam = lambda v, qs: (T.variant == v) & T.quantity.isin(qs) & (T.ingredient != "JOINT") & T.p.notna()
    families = {"primary": ["delta_rise", "beta"], "S1_entry": ["delta_rise", "beta"], "S2_maxdeep": ["delta_rise", "beta"],
                "S3_fixlevel": ["delta_rise", "beta"], "S5_secondary": ["delta_rise", "beta"], "S6a_commoncc": ["delta_rise", "beta"],
                "S6b_deepening_only": ["delta_rise", "beta"]}
    for v, qs in families.items():
        m = fam(v, qs)
        p = T.loc[m, "p"].values; n = len(p); o = np.argsort(p)
        adj = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
        qv = np.empty(n); qv[o] = np.minimum(adj, 1.0)
        T.loc[m, "q"] = qv
    m = (T.variant == "primary") & (T.quantity == "delta_fall")
    p = T.loc[m, "p"].values; n = len(p); o = np.argsort(p)
    adj = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    qv = np.empty(n); qv[o] = np.minimum(adj, 1.0); T.loc[m, "q"] = qv
    m = T.quantity == "inter"
    p = T.loc[m, "p"].values; n = len(p); o = np.argsort(p)
    adj = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    qv = np.empty(n); qv[o] = np.minimum(adj, 1.0); T.loc[m, "q"] = qv
    T.to_csv(f"{OUT}/tests.csv", index=False, float_format="%.5g")
    json.dump({"n_cyclones": R["n_cyclones"], "n_with_fix": R["n_with_fix"], "fix_rows": len(fx)}, open(f"{OUT}/populations.json", "w"))
    # monthly profile at the primary reference fix
    prof = []
    for basin in ("atl", "pac"):
        d = pri[pri.basin == basin]
        for x in PRIMARY + SECOND:
            z = (d[x] - d[x].mean()) / d[x].std(ddof=0)
            for m_ in (10, 11, 12, 1, 2, 3, 4):
                s = d[d.gm == m_]
                prof.append(dict(basin=basin, ingredient=x, month=m_, mean_z=z[s.index].mean(), mean_raw=s[x].mean(), n=len(s),
                                 missing=s[x].isna().mean(), hf_share=s.hf.mean(), p_deep=s.deep.mean()))
    pd.DataFrame(prof).to_csv(f"{OUT}/monthly_profile.csv", index=False, float_format="%.5g")
    print("done")


if __name__ == "__main__":
    main()
