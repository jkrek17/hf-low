"""Matched sensitivity: conditional logistic on matched case-control sets (PREREGISTRATION.md).

For each case in a half, the nearest 2 sampled controls of the same month and half by standardised distance on
the reference low's MSLP, latitude and longitude (greedy, no reuse, cases in a seeded random order). The
ingredients (GH, GRAD, STAB, NAO, motion) enter a conditional logit; matching replaces the baseline terms.
Season-block bootstrap (resampling the seasons of the cases; sets follow their case). ERA5 proxy.
Needs analyse.py to have run. Writes results/matched.txt.
"""
import os, sys
import numpy as np, pandas as pd
from statsmodels.discrete.conditional_models import ConditionalLogit

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(os.environ.get("ERA5_WORK", os.path.join(HERE, "work")), "greenland_jets")
RES = os.path.join(HERE, "results")
sys.path.insert(0, HERE)
import analyse as an
out = []


def say(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)


def match(z, rng, k=2):
    sets = []
    used = set()
    cases = z[z.case == 1].sample(frac=1.0, random_state=int(rng.integers(1 << 31)))
    ctrl = z[z.case == 0]
    for ci, c in cases.iterrows():
        pool = ctrl[(ctrl.month == c.month) & (~ctrl.index.isin(used))]
        if len(pool) == 0:
            continue
        d = np.sqrt(((pool[["ref_msl", "ref_lat", "ref_lon"]].values - c[["ref_msl", "ref_lat", "ref_lon"]].values.astype(float)) /
                     z[["ref_msl", "ref_lat", "ref_lon"]].std().values) ** 2).sum(1)
        pick = pool.index[np.argsort(d)[:k]]
        used.update(pick)
        sets.append((ci, list(pick)))
    return sets


def fit(z, sets, cols):
    rows, g, y = [], [], []
    for n, (ci, ps) in enumerate(sets):
        for j, i in enumerate([ci] + ps):
            rows.append(z.loc[i, cols].values.astype(float)); g.append(n); y.append(1 if j == 0 else 0)
    m = ConditionalLogit(np.array(y), np.array(rows), groups=np.array(g)).fit(disp=0)
    return m.params


def main():
    zall = pd.read_pickle(f"{WORK}/zall.pkl")
    cols = an.ING + an.MOT
    zst = pd.read_pickle(f"{WORK}/zstats.pkl")
    for half in ("heldout", "discovery"):
        z = zall[zall.half == half].reset_index(drop=True)
        rng = np.random.default_rng(an.SEED + 60)
        sets = match(z, rng)
        say(f"{half}: {len(sets)} matched sets (cases {int(z.case.sum())}), controls per set "
            f"{np.mean([len(p) for _, p in sets]):.2f}")
        dm = np.mean([z.loc[ci, "ref_msl"] - z.loc[p, "ref_msl"].mean() for ci, p in sets]) * zst["ref_msl"][1]
        dla = np.mean([z.loc[ci, "ref_lat"] - z.loc[p, "ref_lat"].mean() for ci, p in sets])
        dlo = np.mean([z.loc[ci, "ref_lon"] - z.loc[p, "ref_lon"].mean() for ci, p in sets])
        dgh = np.mean([z.loc[ci, "GH"] - z.loc[p, "GH"].mean() for ci, p in sets])
        say(f"  match quality, case minus matched controls: reference-low MSLP {dm:+.1f} hPa, lat {dla:+.2f}, lon {dlo:+.2f}; GH {dgh:+.2f} SD")
        b0 = fit(z, sets, cols)
        sidx = {s: [n for n, (ci, _) in enumerate(sets) if z.loc[ci, "season"] == s] for s in sorted(z.season.unique())}
        bs = []
        for _ in range(500):
            ss = rng.choice(list(sidx), len(sidx))
            sub = [sets[n] for s in ss for n in sidx[s]]
            try:
                bs.append(fit(z, sub, cols))
            except Exception:
                pass
        bs = np.array(bs)
        for i, c in enumerate(cols):
            say(f"  {c:5s} OR per SD {np.exp(b0[i]):.2f} [{np.exp(np.percentile(bs[:, i], 2.5)):.2f}, "
                f"{np.exp(np.percentile(bs[:, i], 97.5)):.2f}]  p {min(1, 2 * min((bs[:, i] >= 0).mean(), (bs[:, i] <= 0).mean())):.3f}")
    open(f"{RES}/matched.txt", "w").write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
