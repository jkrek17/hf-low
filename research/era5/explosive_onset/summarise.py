"""Combine the Arm A and Arm B test tables (results/analysis_{B,A}_tests.csv, subsets) into the FDR accounting of the plan.
Run: python3 -I research/era5/explosive_onset/summarise.py <repo_root>"""
import sys, os, numpy as np, pandas as pd
root = sys.argv[1]; RES = os.path.join(root, "research/era5/explosive_onset/results")
def bh(p):
    p = np.asarray(p, float); o = np.argsort(p); m = len(p); q = np.empty(m); r = p[o] * m / (np.arange(m) + 1)
    q[o] = np.minimum.accumulate(r[::-1])[::-1]; return np.minimum(q, 1)
R = pd.concat([pd.read_csv(os.path.join(RES, f"analysis_{a}_tests.csv")) for a in "BA"], ignore_index=True)
S = pd.concat([pd.read_csv(os.path.join(RES, f"analysis_{a}_subsets.csv")) for a in "BA"], ignore_index=True)
P = R[R.baseline == "M1"].copy(); P["q_primary_all9"] = bh(P.p_decision)
L = ["PRIMARY FAMILY (within-time test vs the M1 tree baseline; p = larger of LR and season-bootstrap Wald):"]
L.append(P[["arm", "group", "k", "p_lr", "p_boot", "p_decision", "q_within_arm", "q_primary_all9", "dBSS_S", "dBSS_lo5", "dBSS_hi95"]].to_string(index=False, float_format=lambda x: f"{x:.4g}"))
L.append(f"primary tests passing q<0.05: within arm {int((P.q_within_arm<0.05).sum())} of {len(P)}; across all 9 {int((P.q_primary_all9<0.05).sum())} of {len(P)}")
rows = []
for _, r in R.iterrows():
    rows.append((f"{r.arm} {r.group} vs {r.baseline}: within-time LR/bootstrap", r.p_decision))
    rows.append((f"{r.arm} {r.group} vs {r.baseline}: held-out log-lik sign-flip", r.heldout_signflip_p))
    rows.append((f"{r.arm} {r.group} vs {r.baseline}: BSS sign-flip", r.dBSS_signflip_p))
for _, r in S.iterrows(): rows.append((f"{r.arm} 'all' model, subset {r.subset}: within-time", r.p_decision))
T = pd.DataFrame(rows, columns=["test", "p"]); T["q_all_tests"] = bh(T.p)
L.append(f"\nALL TESTS RUN ({len(T)}): passing q<0.05 over the whole set: {int((T.q_all_tests<0.05).sum())}")
L.append(T[T.q_all_tests < 0.05].to_string(index=False, float_format=lambda x: f"{x:.4g}"))
sec = R[R.baseline == "M0"]; L.append(f"\nSECONDARY (vs M0 logistic), within-time: {int((sec.q_secondary<0.05).sum())} of {len(sec)} pass q<0.05 within the secondary family")
T.to_csv(os.path.join(RES, "all_tests_fdr.csv"), index=False); P.to_csv(os.path.join(RES, "primary_family.csv"), index=False)
open(os.path.join(RES, "summary.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
