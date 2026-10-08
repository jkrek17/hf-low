"""Run every analysis in PREREGISTRATION.md. usage: run.py OUT_DIR"""
import os, sys
import numpy as np, pandas as pd
import common as C

OUT = sys.argv[1]; os.makedirs(OUT, exist_ok=True)
T, I = C.load_tracks(), C.indices()
SPEC = [  # id, label, basin, s0, s1, now, basin_by
    ("P1", "primary: Atlantic, lagged, 1979-2000, minimum-pressure basin", "atl", 1979, 2000, False, "min"),
    ("S3", "basin at first fix", "atl", 1979, 2000, False, "gen"),
    ("S4", "same-time predictors (contrast)", "atl", 1979, 2000, True, "min"),
    ("S5", "bridge: source era 2001-2013 (positive control, not independent)", "atl", 2001, 2013, False, "min"),
    ("S6", "Pacific control, 1979-2000", "pac", 1979, 2000, False, "min"),
]
rows = []
for sid, lab, basin, s0, s1, now, bb in SPEC:
    D, E = C.table(T, I, basin, "NAO", "PNA", s0, s1, now=now, basin_by=bb)
    r = C.core.location_test(E)
    for k, nm in enumerate(("lon", "lat")):
        rows.append(dict(
                         test=sid, label=lab, outcome=nm, basin=basin, winters="%d-%d" % (s0, s1), gamma=r["gamma"][k],
                         se_cr1=r["se_gamma"][k], ci_lo=r["ci_gamma"][k, 0], ci_hi=r["ci_gamma"][k, 1],
                         p_two=r["p_lon"] if k == 0 else r["p_lat"], p_joint=r["p_joint"], n_tracks=r["n_events"],
                         n_winters=r["n_seasons"], betaNAO=r["A"][k], betaPNA=r["B"][k], y_sd=r["y_sd"][k]))
R = pd.DataFrame(rows)
g = R.gamma.values; p2 = R.p_two.values
R["p_one"] = np.where(g > 0, p2 / 2, 1 - p2 / 2)
# BH over the 7 tests of the plan: P1 lon, S1 = P1 lat, S2 = P1 joint, S3, S4, S5, S6 (lon, one-sided for P1/S3/S4/S5 predicted +; two-sided for S1, S2, S6)
sel = []
def pick(test, outcome): return R[(R.test == test) & (R.outcome == outcome)].iloc[0]
plan = [("P1", "P1 lon (one-sided)", pick("P1", "lon").p_one),
        ("S1", "S1 lat, two-sided", pick("P1", "lat").p_two),
        ("S2", "S2 joint 2 df", pick("P1", "lon").p_joint),
        ("S3", "S3 first-fix basin lon (one-sided)", pick("S3", "lon").p_one),
        ("S4", "S4 same-time lon (one-sided)", pick("S4", "lon").p_one),
        ("S5", "S5 bridge lon (one-sided)", pick("S5", "lon").p_one),
        ("S6", "S6 Pacific lon (two-sided)", pick("S6", "lon").p_two)]
P = pd.DataFrame(plan, columns=["id", "test", "p"])
P["q_BH"] = C.core.bh(P.p.values)
R.to_csv(os.path.join(OUT, "results.csv"), index=False)
P.to_csv(os.path.join(OUT, "plan_tests.csv"), index=False)
with open(os.path.join(OUT, "summary.txt"), "w") as f:
    f.write(R.round(4).to_string(index=False) + "\n\n" + P.round(4).to_string(index=False) + "\n")
print(R.round(3).to_string(index=False)); print(P.round(4).to_string(index=False))
