"""Combined Benjamini-Hochberg over this directory's primary tests and the two sibling threads' primary tests.

Sibling p-values are copied from their committed READMEs on the research branch (read 2026-10-08):
  enso_pna (ONI x PNA, Pacific HF count, 2004-05..2025-26): P1 interaction permutation p 0.71; P2 indirect effect p 0.036
  enso_kuroshio (EMI flavour, DJF Kuroshio-box genesis, 47 winters): P1 permutation p 0.895
Their families were corrected separately by their authors; this puts all eleven on one scale.
usage: family.py RESULTS_CSV OUT_TXT
"""
import sys
import pandas as pd
import core

R = pd.read_csv(sys.argv[1])
P = R[R.tier == "P"]
rows = []
for t, g in P.groupby("test"):
    lj = g[g.outcome.str.startswith("L:lon")].p_joint.iloc[0]
    cp = g[g.outcome.str.startswith("C")].p.iloc[0]
    rows += [(t + "L", lj), (t + "C", cp)]
rows += [("ENSOxPNA P1 (sibling)", 0.71), ("ENSOxPNA P2 (sibling)", 0.036), ("Kuroshio P1 (sibling)", 0.895)]
df = pd.DataFrame(rows, columns=["test", "p"])
df["q_combined"] = core.bh(df.p.values)
txt = "Combined BH over %d primary tests (this directory %d + 3 sibling)\n%s\n" % (len(df), len(df) - 3, df.round(3).to_string(index=False))
print(txt)
open(sys.argv[2], "w").write(txt)
