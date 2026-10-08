"""Readable summary of results/*.csv into results/summary.txt. usage: report.py RESULTS_DIR"""
import os, sys
import pandas as pd
R = sys.argv[1]
pd.set_option("display.width", 200); pd.set_option("display.max_rows", 500)
T = pd.read_csv(os.path.join(R, "tests.csv")); D = pd.read_csv(os.path.join(R, "decomposition.csv"))
S = pd.read_csv(os.path.join(R, "secondary.csv")); P = pd.read_csv(os.path.join(R, "power_mde.csv"))
PP = pd.read_csv(os.path.join(R, "power_position.csv"))
out = ["HEMISPHERIC PATTERN CHANNELS, ERA5 proxy, pipeline A, seasons 2004-05..2025-26, per SD of the PR 41 out-of-sample index",
       "p: season-block permutation, two-sided, 2000 draws (floor 0.0005). CI: 95% season-block bootstrap (positions: clustered SE, t on 21 df).", ""]
out.append("PRIMARY (24 tests)"); out.append(T.drop(columns=["label"]).round(4).to_string(index=False))
n = int((T.q_all24 < 0.05).sum())
out.append(f"\nPass q_all24 < 0.05: {n} of {len(T)}; core six (T1-T3 x 2 basins) pass q_core6 < 0.05: {int((T.q_core6 < 0.05).sum())} of 6\n")
out.append("DECOMPOSITION (+1 SD, summed over 660 weeks; fractions of the change in HF lows)"); out.append(D.round(4).to_string(index=False))
out.append("\nSECONDARY (own BH family)"); out.append(S.round(4).to_string(index=False))
out.append("\nPOWER: RR per SD with 80% power (NaN = not reached by RR 1.20); positions: degrees per SD"); out.append(P.round(3).to_string(index=False)); out.append(PP.round(2).to_string(index=False))
open(os.path.join(R, "summary.txt"), "w").write("\n".join(out) + "\n")
print("\n".join(out[-0:3]))
