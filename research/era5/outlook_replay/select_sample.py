"""Systematic replay sample, fixed before looking at results: every 10th 00 UTC cycle (sorted order, offset 5) of the 1,061 cycles, 1 Oct-30 Apr
2021-22 to 2025-26. Unlike the 20 picked cases (replay_cases.csv) nothing here is chosen for interest, so coverage over this sample is an unbiased
read of the shading. Cycles already in replay_cases.csv are left out of this file. usage: select_sample.py"""
import os, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, "results")
d = pd.read_csv(os.path.join(RES, "point_days.csv"), dtype={"cycle": str}); cyc = sorted(d.cycle.unique())
have = set(pd.read_csv(os.path.join(RES, "replay_cases.csv"), dtype={"cycle": str}).cycle)
rows = [(int(c[:4]) if int(c[4:6]) >= 10 else int(c[:4]) - 1, c, "systematic") for c in cyc[5::10] if c not in have]
pd.DataFrame(rows, columns=["season", "cycle", "why"]).to_csv(os.path.join(RES, "replay_sample.csv"), index=False); print(len(rows), "sample cycles")
