"""Test 2 maps: HF-onset composites (rotated to motion) under the top and bottom tercile of the pattern index, with pixelwise season-block bootstrap and BH. Descriptive.
usage: maps_pp.py ERA5_WORK [B_PIXEL]. Writes results/comp_<basin>.npz and pp_maps_summary.txt"""
import sys, os, glob, numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "hf_vs_storm"))
from strat import cells, boot, bh
W = sys.argv[1]; BP = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results"); rng = np.random.default_rng(20261012); SEASONS = np.arange(2004, 2026)
H = pd.read_csv(os.path.join(OUT, "hf_index.csv")); H = H[H.in_scope & H.tercile.isin(["top", "bottom"])]
boxes = {}
for d in ("hf_vs_storm", "hf_pattern_phase"):
    for f in glob.glob(os.path.join(W, d, "boxes", "*.npz")):
        z = np.load(f)
        for k in z.files:
            if "_HF_onset_" in k and k.endswith("_m"): boxes[k] = z[k]
lines = []
for b in ("atl", "pac"):
    x = H[H.basin == b]; x = x[[f"{t}_HF_onset_gust_m" in boxes for t in x.track]].reset_index(drop=True)
    x = x[[np.isfinite(boxes[f"{t}_HF_onset_gust_m"].astype(np.float32)).any() for t in x.track]].reset_index(drop=True)
    grp = (x.tercile == "bottom").astype(int).values; out = {"n_top": int((grp == 0).sum()), "n_bottom": int((grp == 1).sum())}
    for fld in ("gust", "ws", "msl", "d2m"):
        arr = np.stack([boxes[f"{t}_HF_onset_{fld}_m"].astype(np.float32).ravel() for t in x.track])
        sm, cn = cells(arr, grp, x.season.values, x.mon.values - 1, SEASONS, 12)
        obs, _, p = boot(sm, cn, BP, rng, keep=False); qq, rej = bh(p)
        for nm, v in (("diff", obs[0]), ("top", obs[1]), ("bottom", obs[2])): out[f"{fld}_{nm}"] = v.reshape(121, 121).astype(np.float32)
        out[f"{fld}_rej"] = rej.reshape(121, 121)
        lines.append(f"{b} {fld}: top {out['n_top']}, bottom {out['n_bottom']} rotated storms; {int(rej.sum())} of {rej.size} pixels pass BH q<0.05")
    np.savez_compressed(os.path.join(OUT, f"comp_{b}.npz"), **out)
open(os.path.join(OUT, "pp_maps_summary.txt"), "w").write("\n".join(lines) + "\n"); print("\n".join(lines))
