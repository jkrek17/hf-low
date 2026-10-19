"""Planted-effect power for T2 (LOSO Poisson deviance skill of the index over month + previous-week count, season-block permutation).
Run before any index-outcome association is read. Draws ignore the season-to-season variation in the index effect, so power is an upper bound."""
import sys, json
import numpy as np
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
import bomblib as B
def skill(df, y, idx):
    seasons = sorted(df.season.unique()); d1 = d0 = 0.0
    Xb = B.design(df); Xf = B.design(df, [idx])
    for s in seasons:
        tr = (df.season != s).values; te = ~tr
        for X, tag in ((Xb, 0), (Xf, 1)):
            b = B.poisson_fit(X[tr], y[tr]); mu = np.exp(np.clip(X[te] @ b, -20, 8)); dd = B.dev(y[te], mu)
            if tag: d1 += dd
            else: d0 += dd
    return 1 - d1 / d0
def perm_p(df, y, idx, nperm, rng, s_obs):
    seasons = sorted(df.season.unique()); blocks = {s: idx[(df.season == s).values] for s in seasons}; ge = 0
    for _ in range(nperm):
        order = rng.permutation(seasons); new = np.concatenate([blocks[o] for o in order]); 
        # blocks equal length (30 weeks), keep season labels in place
        ge += skill(df, y, new) >= s_obs
    return (ge + 1) / (nperm + 1)
if __name__ == "__main__":
    rng = np.random.default_rng(3); out = {}
    for b in ("atl", "pac"):
        df = B.weekly(b); idx = df.idx.values; z = (idx - idx.mean()) / idx.std()
        mm = df.groupby("month").bomb.transform("mean").values
        res = {}
        for rr in (1.0, 1.15, 1.3):
            hit = 0; nsim = 120
            for _ in range(nsim):
                y = rng.poisson(mm * rr ** z).astype(float); d2 = df.copy(); d2["bomb"] = y
                d2["prev"] = d2.groupby("season").bomb.shift(1).fillna(y.mean())
                sk = skill(d2, y, idx)
                if sk > 0 and perm_p(d2, y, idx, 60, rng, sk) < 0.05: hit += 1
            res[rr] = hit / nsim
        out[b] = res; print(b, res, flush=True)
    json.dump(out, open(__import__("os").path.join(B.HERE, "results/power.json"), "w"))
