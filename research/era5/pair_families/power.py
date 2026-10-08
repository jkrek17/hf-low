"""Step 0b: null mean and SD of candidate-pair counts under whole-track time permutation within season-month,
and minimum detectable excess (planted-daughter simulation). No observed orientation is looked at: only the observed
candidate COUNTS (already disclosed in count_pairs.py)."""
import sys, json
import numpy as np
from common import load
from engine import *
NPERM = int(sys.argv[1]) if len(sys.argv) > 1 else 500
T, F = load()
res = {}
for b, G in T.groupby("basin"):
    G = G.reset_index(drop=True)
    P = prep(G, F)
    st = strata(G)
    rng = np.random.default_rng(20261008)
    obs = summarise(pairs(P), P)
    null = [summarise(pairs(P, perm_shift(P, st, rng)), P) for _ in range(NPERM)]
    res[b] = dict(n_tracks=len(G), observed_counts=obs,
                  null_mean={k: float(np.mean([x[k] for x in null])) for k in obs},
                  null_sd={k: float(np.std([x[k] for x in null], ddof=1)) for k in obs})
    print(b, json.dumps(res[b]))
json.dump(res, open("results_power_counts.json", "w"), indent=1)
