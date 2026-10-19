#!/usr/bin/env python3
"""Freeze the numbers the weekly pattern loading note needs into data/pattern_loading/ (small, committed). Run once from a checkout of the research branch.

    python3 tools/pattern_loading_freeze.py RESEARCH_ROOT        (RESEARCH_ROOT = the research/era5 folder of claude/exciting-fermat-8vcvgq)

Reads: gfs_test/results/pr41_eofs_zu.npz (PR 41 EOFs, discovery seasons 2004-2014, z500 and u250), gfs_test/results/gefs_climatology.npz (GEFS v12 control model
climatology, one 12 x 64 map per week k and variable), hemispheric/results/frozen_primary.json (frozen betas of model P), basin_bomb_score/results/pattern_note_{table.csv,cuts.json}.
Writes frozen.npz (EOF and climatology arrays) and params.json (betas for the z and u PC blocks, class cut points, class table)."""
import json, os, sys
import numpy as np, pandas as pd

root = sys.argv[1]; out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "pattern_loading")
e = np.load(os.path.join(root, "gfs_test/results/pr41_eofs_zu.npz")); c = np.load(os.path.join(root, "gfs_test/results/gefs_climatology.npz"))
np.savez_compressed(os.path.join(out, "frozen.npz"), w=e["w"], mean_z500=e["mean_z500"], mean_u250=e["mean_u250"], vec_z500=e["vec_z500"], vec_u250=e["vec_u250"],
                    sd_z500=e["sd_z500"], sd_u250=e["sd_u250"], clim_z500=c["z500"], clim_u250=c["u250"])
fp = json.load(open(os.path.join(root, "hemispheric/results/frozen_primary.json")))
cuts = json.load(open(os.path.join(root, "basin_bomb_score/results/pattern_note_cuts.json")))
tab = pd.read_csv(os.path.join(root, "basin_bomb_score/results/pattern_note_table.csv"))
P = {"source": "research/era5 on claude/exciting-fermat-8vcvgq; index = frozen PR 41 model P z500 and u250 PC blocks (beta[8:18], beta[18:28]) on GEFS v12 control anomalies",
     "basins": {}}
for b in ("atl", "pac"):
    beta = fp["basins"][b]["models"]["P"]["beta"]
    P["basins"][b] = {"beta_z": beta[8:18], "beta_u": beta[18:28], **cuts[b],
                      "table": {r.cls: {k: float(getattr(r, k)) for k in ("weeks", "hf_per_week", "bomb_per_week", "share_weeks_hf", "share_weeks_bomb")} for r in tab[tab.basin == b].itertuples()}}
json.dump(P, open(os.path.join(out, "params.json"), "w"), indent=1)
print("wrote", os.listdir(out))
