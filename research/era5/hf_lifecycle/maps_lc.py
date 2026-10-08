"""Test 3 stage B composites at lags -48, -24, 0 h (rotated to motion at the peak): gust and MSLP, HF vs storm-force-only, pixelwise season-block bootstrap + BH. Descriptive.
usage: maps_lc.py REPO_ROOT ERA5_WORK [B_PIXEL]. Writes results/comp_<basin>_lag<L>.npz"""
import sys, os, glob, numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "hf_vs_storm"))
from strat import cells, boot, bh
root, W = sys.argv[1], sys.argv[2]; BP = int(sys.argv[3]) if len(sys.argv) > 3 else 1000
OUT = os.path.join(root, "research/era5/hf_lifecycle/results"); rng = np.random.default_rng(20261012); SEASONS = np.arange(2004, 2026)
V = pd.read_csv(os.path.join(OUT, "lc_values.csv.gz")); V = V[V.snapped]
boxes = {}
for f in glob.glob(os.path.join(W, "hf_lifecycle", "boxes", "*.npz")):
    z = np.load(f)
    for k in z.files:
        if k.endswith("_gust") or k.endswith("_msl"): boxes[k] = z[k]
S = pd.read_csv(os.path.join(OUT, "lc_sample.csv")).set_index("track")
for f in glob.glob(os.path.join(W, "hf_vs_storm", "boxes", "*.npz")):
    z = np.load(f)
    for k in z.files:
        t, rest = k.split("_", 1)
        if rest in ("HF_peak_gust_m", "HF_peak_msl_m", "SF_peak_gust_m", "SF_peak_msl_m") and int(t) in S.index:
            boxes[f"{t}_lag+0_{rest.split('_')[2]}"] = z[k]
lines = []
for b in ("atl", "pac"):
    for L in (-48, -24, 0):
        x = V[(V.basin == b) & (V.lag == L)]; nm = f"lag{L:+d}"
        x = x[[f"{t}_{nm}_gust" in boxes and np.isfinite(boxes[f"{t}_{nm}_gust"].astype(np.float32)).any() for t in x.track]].reset_index(drop=True)
        out = {"n_hf": int((x.grp == 0).sum()), "n_sf": int((x.grp == 1).sum())}
        for fld in ("gust", "msl"):
            arr = np.stack([boxes[f"{t}_{nm}_{fld}"].astype(np.float32).ravel() for t in x.track])
            sm, cn = cells(arr, x.grp.values, x.season.values, x.mon.values - 1, SEASONS, 12)
            obs, _, p = boot(sm, cn, BP, rng, keep=False); qq, rej = bh(p)
            for n_, v in (("diff", obs[0]), ("hf", obs[1]), ("sf", obs[2])): out[f"{fld}_{n_}"] = v.reshape(121, 121).astype(np.float32)
            out[f"{fld}_rej"] = rej.reshape(121, 121)
            lines.append(f"{b} lag {L} h {fld}: HF {out['n_hf']}, SF {out['n_sf']} rotated storms; {int(rej.sum())} of {rej.size} pixels pass BH q<0.05")
        np.savez_compressed(os.path.join(OUT, f"comp_{b}_lag{L}.npz"), **out)
open(os.path.join(OUT, "stage_b_maps_summary.txt"), "w").write("\n".join(lines) + "\n"); print("\n".join(lines))
