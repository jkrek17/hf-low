"""RA-28 summary: primary family (4 tests) with BH q, verdicts from the pre-registered rules, global BH over all pre-registered p. No ERA5 access."""
import numpy as np, pandas as pd
T = pd.read_csv("results/conv_tests.csv"); P = pd.read_csv("results/paired_delta.csv"); Q = pd.read_csv("results/part2_tests.csv")
def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = np.empty(n); prev = 1.0
    for rank, i in zip(range(n, 0, -1), o[::-1]): prev = min(prev, p[i] * n / rank); q[i] = prev
    return q
out = []
say = lambda s="": (print(s), out.append(s))
prim = []
for b in ("atl", "pac"):
    r = T[(T.label == "P99") & (T.basin == b) & (T.test == "A2")].iloc[0]; prim.append((f"A2(P99) vs 1, {b}", r.p, r.rr, r.rr_lo, r.rr_hi))
for b in ("atl", "pac"):
    r = P[(P.label == "P99") & (P.basin == b)].iloc[0]; prim.append((f"Delta(P99-M), {b}", r.p, r.delta_rr, r.lo, r.hi))
q = bh([x[1] for x in prim])
say("PRIMARY FAMILY (BH within 4)")
for x, qq in zip(prim, q): say(f"  {x[0]:24s} est {x[2]:.3f} [{x[3]:.3f}, {x[4]:.3f}]  p {x[1]:.4f}  q {qq:.4f}")
say()
say("VERDICTS (Delta_RR = RR(P99) - RR(M); strengthens lo>0; weakens hi<0; keeps interval inside +-0.06; else can't tell)")
for b in ("atl", "pac"):
    r = P[(P.label == "P99") & (P.basin == b)].iloc[0]
    v = "strengthens" if r.lo > 0 else "weakens" if r.hi < 0 else "keeps (well-powered)" if (r.lo >= -0.06 and r.hi <= 0.06) else "can't tell"
    ref = {"atl": (1.252, 1.30), "pac": (1.137, 1.20)}[b]
    say(f"  {b}: {v}; RR M {r.rr_M:.3f}, P99 {r.rr_label:.3f}, Delta {r.delta_rr:+.3f} [{r.lo:+.3f}, {r.hi:+.3f}], MDE {r.mde_delta:.3f}; "
        f"P99 RR >= {ref[1]}: {r.rr_label >= ref[1]}; |P99 - published headline {ref[0]}| = {abs(r.rr_label - ref[0]):.3f} (within 0.03: {abs(r.rr_label - ref[0]) <= 0.03})")
allp = list(T.p) + list(P.p) + list(Q[Q.quantity == "beta_tot"].p) + list(Q[Q.quantity == "a_path"].p) + list(Q[Q.quantity == "b_path_OR"].p)
allp = np.array([x for x in allp if x == x]); qa = bh(allp)
say(); say(f"GLOBAL BH over all {len(allp)} pre-registered p-values (Part 1 A1-A3 x 4 labels x 2 basins, paired Deltas, Part 2 beta_tot, a-paths, b-paths): {(qa < 0.05).sum()} pass q < 0.05")
open("results/summary.txt", "w").write("\n".join(out) + "\n")
