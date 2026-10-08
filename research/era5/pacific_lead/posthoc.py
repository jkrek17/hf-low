"""RA-5 post hoc (not pre-registered): the naive month-only lagged correlation of the two indices, to show what the
same-week anticorrelation does to a lead-lag reading. Labelled post hoc in the README."""
import json
import os

import numpy as np

import leadlib as L

out = {}
for ds in (L.load_2004("archive"), L.load_1979_2000()):
    M = ds.M[:, 3:27].reshape(-1, ds.M.shape[2])
    P = M @ np.linalg.pinv(M)
    A = ds.idx["atl"][:, 3:27].ravel(); A = A - P @ A
    rows = {}
    for k in range(-3, 4):
        Q = ds.idx["pac"][:, 3 - k:27 - k].ravel(); Q = Q - P @ Q
        rows[k] = float(np.corrcoef(A, Q)[0, 1])
    out[ds.label] = rows
json.dump(out, open(os.path.join(L.HERE, "results", "posthoc_naive_xcorr.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
