"""Part A summary: BH over the defined placebo tests (P2 Pacific is degenerate and excluded), the pre-registered gate, readings, shuffle comparison."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import placelib as P
R = sys.argv[1]
D = pd.read_csv(f"{R}/partA_tests.csv"); I = pd.read_csv(f"{R}/partA_indices.csv"); S = pd.read_csv(f"{R}/partA_shuffle.csv")
bad = (D["index"] == "P2") & (D.basin == "pac")
D.loc[bad, ["rr", "rr_lo", "rr_hi", "p", "f", "f_lo", "f_hi", "A4_rr", "A4_lo", "A4_hi"]] = np.nan
D["q_placebo"] = np.nan
m = D["index"].isin(["P1", "P2", "P3"]) & ~bad
D.loc[m, "q_placebo"] = P.C.bh(D.p[m].values)
D.to_csv(f"{R}/partA_tests.csv", index=False)
out = []
out.append("PART A. placebo pattern indices (ERA5 proxy, pipeline A; per +1 SD, 22 seasons)\n")
out.append("Reproduction gate (r with committed oos_index >= 0.90 in both basins, A2 inside PR 85 interval):")
for b in ("atl", "pac"):
    r = I[(I["index"] == "P0") & (I.basin == b)].r_committed.iloc[0]
    a2 = D[(D["index"] == "P0") & (D.basin == b) & (D.test == "A2")].iloc[0]
    out.append(f"  {b}: r = {r:.3f}; A2 rate ratio {a2.rr:.3f} ({a2.rr_lo:.3f}-{a2.rr_hi:.3f})")
out.append("  Gate NOT met (Pacific r < 0.90). Per the plan, P0 refit is the reference; the committed index is not used in Part A.\n")
out.append("Placebo indices (A4 = share of cyclones reaching HF; A2 = P(HF | deepening); f_conv = A2/A4 in logs):")
for (ix, b), g in D[D["index"].isin(["P0", "P1", "P2", "P3"])].groupby(["index", "basin"], sort=False):
    a2 = g[g.test == "A2"].iloc[0]
    if np.isnan(a2.rr):
        out.append(f"  {ix} {b}: index degenerate (cross-validation chose no PCs in every fold); no placebo possible"); continue
    ok = (a2.A4_lo > 1 or a2.A4_hi < 1) and a2.frac_boot_ok >= 0.95
    reading = "no resolved share effect (reading 2)" if not ok else ("conversion dominates (reading 1)" if a2.f_lo > 0.5 else "mixed (reading 3)")
    out.append(f"  {ix} {b}: A4 {a2.A4_rr:.3f} ({a2.A4_lo:.3f}-{a2.A4_hi:.3f}); A2 {a2.rr:.3f} ({a2.rr_lo:.3f}-{a2.rr_hi:.3f}), q {a2.q_placebo if ix!='P0' else float('nan'):.3f}; f_conv {a2.f:.2f} ({a2.f_lo:.2f} to {a2.f_hi:.2f}); frac_boot_ok {a2.frac_boot_ok:.3f} -> {reading}")
out.append("\nShuffle null (100 season-shuffled P0 indices per basin): S = |A2| / (|A1|+|A2|+|A3|)")
for b in ("atl", "pac"):
    g = D[(D["index"] == "P0") & (D.basin == b)].set_index("test").rr
    a = np.abs(np.log(g[["A1", "A2", "A3"]].astype(float)))
    So = a.A2 / a.sum()
    allb = S[S.basin == b]
    s = allb.S.dropna()                                  # a shuffled index can be constant (no PCs chosen in any fold): S undefined
    out.append(f"  {b}: S(P0) = {So:.3f}; {len(s)} of {len(allb)} shuffles have a non-constant index; among those, median {s.median():.3f}, 95th pct {s.quantile(.95):.3f}; share >= observed {(s >= So).mean():.3f} ({int((s >= So).sum())} of {len(s)}); median |A4| {allb.A4.abs().median():.3f}")
open(f"{R}/partA-result.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out))
