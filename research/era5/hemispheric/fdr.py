"""Benjamini-Hochberg over the pre-registered families (PREREGISTRATION.md). usage: fdr.py [RESULTS_DIR]"""
import os, sys, json
import numpy as np, pandas as pd
RES = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
J = lambda n: json.load(open(os.path.join(RES, f"heldout_{n}.json")))
rows = []
def add(fam, test, split, b, p, stat):
    rows.append(dict(family=fam, test=test, split=split, basin=b, p=p, stat=stat))
P = J("primary")
for b in ("atl", "pac"):
    for t in ("P1", "P3", "P4"):
        add("1", t, "primary", b, P[b][t]["p"], P[b][t].get("ss", P[b][t].get("delta")))
for nm in ("swap", "lag2"):
    S = J(nm)
    for b in ("atl", "pac"):
        for t in ("P1", "P3", "P4"):
            add("2", t, nm, b, S[b][t]["p"], S[b][t].get("ss", S[b][t].get("delta")))
for b in ("atl", "pac"):
    add("2", "S3_SOM", "primary", b, P[b]["SOM"]["p"], P[b]["SOM"]["ss"])
    add("2", "S4_proxy_transfer", "primary", b, P["S4_proxy"][b]["p"], P["S4_proxy"][b]["ss"])
    add("2", "S5_within_era", "1979-2000", b, P["S5_within_era"][b]["p"], P["S5_within_era"][b]["rr_per_sd"])
df = pd.DataFrame(rows)
def bh(p):
    p = np.asarray(p); n = len(p); o = np.argsort(p)
    q = np.empty(n); a = p[o] * n / (np.arange(n) + 1)
    q[o] = np.minimum(np.minimum.accumulate(a[::-1])[::-1], 1)
    return q
f1 = df[df.family == "1"].copy()
f1["q_family1"] = bh(f1.p.values)
df["q_family1"] = np.nan; df.loc[f1.index, "q_family1"] = f1.q_family1
df["q_family2_all"] = bh(df.p.values)
df.to_csv(os.path.join(RES, "fdr_table.csv"), index=False)
pd.set_option("display.width", 200)
print(df.round(4).to_string())
