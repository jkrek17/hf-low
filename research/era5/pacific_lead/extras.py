"""RA-5 secondary S2 (direction difference), S5 (lag profile), S6 (same-week contrast), S7 (halves). Descriptive; no p-values
except the pre-registered permutation envelope of S5. See PREREGISTRATION.md."""
import json
import os

import numpy as np
import pandas as pd

import leadlib as L

OUT = os.path.join(L.HERE, "results")
rng = np.random.default_rng(L.SEED + 99)
LAGS = list(range(-3, 4))
W3 = slice(3, 27)  # weeks 3..26: common sample for lags -3..+3


def profile_r(ds, qidx, k, seas=None):
    """Partial r of Q(w-k) with A(w), controls month + A(w-1)."""
    seas = np.arange(ds.ns) if seas is None else np.asarray(seas)
    A = ds.idx["atl"][seas]
    y = A[:, W3].ravel()
    M = ds.M[seas][:, W3].reshape(-1, ds.M.shape[2])
    X = np.column_stack([M, A[:, 2:26].reshape(-1, 1), qidx[seas][:, 3 - k:27 - k].reshape(-1, 1)])
    return L.fit(y, X, "ols")[1]


def profile(ds, nperm=2000):
    res = []
    perms = [np.random.default_rng(L.SEED + 31 + i).permutation(ds.ns) for i in range(nperm)]
    for k in LAGS:
        r = profile_r(ds, ds.idx["pac"], k)
        nul = np.array([profile_r(ds, ds.idx["pac"][p], k) for p in perms])
        bt = np.array([profile_r(ds, ds.idx["pac"], k, rng.integers(0, ds.ns, ds.ns)) for _ in range(500)])
        res.append(dict(data=ds.label, lag=k, partial_r=r, boot_lo=np.percentile(bt, 2.5), boot_hi=np.percentile(bt, 97.5),
                        null_lo=np.percentile(nul, 2.5), null_hi=np.percentile(nul, 97.5)))
    return res


def halves(ds):
    out = []
    for nm, seas in (("2004-2014", np.arange(0, 11)), ("2015-2025", np.arange(11, 22))):
        for sp in (L.Spec("P1", "pois", 1, "idx"), L.Spec("P3", "ols", 1, "idx")):
            est, r = L.stat(ds, sp, seas=seas)
            bt = L.boot(ds, sp, rng, 1000, seas_pool=seas)
            lo, hi = np.percentile(bt[:, 0], [2.5, 97.5])
            f = np.exp if sp.kind == "pois" else (lambda v: v)
            out.append(dict(half=nm, test=sp.name, est=float(f(est)), lo=float(f(lo)), hi=float(f(hi)), partial_r=r))
    return out


def direction(ds):
    """Slopes of P3 (Pacific to Atlantic) and its reverse in the same bootstrap draws; difference with 95% interval."""
    fw, rv = L.Spec("P3", "ols", 1, "idx"), L.Spec("rev", "ols", 1, "idx", False, "pac", "atl")
    d, a, b = [], [], []
    for _ in range(2000):
        s = rng.integers(0, ds.ns, ds.ns)
        x, y = L.stat(ds, fw, seas=s)[0], L.stat(ds, rv, seas=s)[0]
        a.append(x); b.append(y); d.append(x - y)
    return dict(pac_to_atl=L.stat(ds, fw)[0], atl_to_pac=L.stat(ds, rv)[0],
                diff=L.stat(ds, fw)[0] - L.stat(ds, rv)[0], diff_lo=np.percentile(d, 2.5), diff_hi=np.percentile(d, 97.5))


def same_week(ds):
    """S6 labelled contrast: same-week correlation of the two indices and of the Pacific index with the Atlantic count."""
    a, p = ds.idx["atl"].ravel(), ds.idx["pac"].ravel()
    return dict(corr_idx_same_week=float(np.corrcoef(a, p)[0, 1]),
                corr_atl_idx_lag1_vs_pac_idx_lag1=float(np.corrcoef(ds.idx["atl"][:, 1:29].ravel(), ds.idx["pac"][:, 1:29].ravel())[0, 1]),
                corr_atl_persistence_lag1=float(np.corrcoef(ds.idx["atl"][:, 1:].ravel(), ds.idx["atl"][:, :-1].ravel())[0, 1]),
                corr_pac_persistence_lag1=float(np.corrcoef(ds.idx["pac"][:, 1:].ravel(), ds.idx["pac"][:, :-1].ravel())[0, 1]))


def main():
    a, r = L.load_2004("archive"), L.load_1979_2000()
    prof = pd.DataFrame(profile(a) + profile(r))
    prof.to_csv(os.path.join(OUT, "lag_profile.csv"), index=False)
    print(prof.round(3).to_string())
    ex = dict(halves=halves(a), direction_2004=direction(a), direction_1979=direction(r), same_week_2004=same_week(a),
              same_week_1979=same_week(r))
    json.dump(ex, open(os.path.join(OUT, "extras.json"), "w"), indent=1, default=float)
    print(json.dumps(ex, indent=1, default=float))


if __name__ == "__main__":
    main()
