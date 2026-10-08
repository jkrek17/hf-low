"""PR 14 primary decomposition (Atlantic NAO, Pacific PNA, lagged, Oct-Apr 2004-05..2025-26) on a derived track table.
Calls freq_split/split.py's decompose() unchanged, primary analyses only, same seed and order as split.py so that k = 1 reproduces PR 14.

usage: run_pr14.py ALL_TRACKS_VARIANT.csv.gz OUT_CSV [NPERM]
"""
import os, sys
import numpy as np, pandas as pd

ERA = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.join(ERA, "freq_split"))
import split as S  # noqa: E402

tracks, out = sys.argv[1], sys.argv[2]
nperm = int(sys.argv[3]) if len(sys.argv) > 3 else 2000
cpc = os.environ.get("CPC_DIR", "/mnt/project-files/teleconnection-test/cpc_indices")
repo = os.path.abspath(os.path.join(ERA, "..", ".."))
rng = np.random.default_rng(S.SEED)
T = S.load_tracks(tracks)
I = S.indices(cpc, repo)
R = []
for basin, idx in (("atl", "NAO"), ("pac", "PNA")):
    r = S.decompose(T, I, basin, idx, "octapr", 2004, 2025, nperm, rng)
    R.append(dict(r, analysis="primary", tracks_file=os.path.basename(tracks)))
    print(basin, idx, f"RR(HF) {np.exp(r['b_hf']):.3f} RR(all) {np.exp(r['b_all']):.3f} RR(share) {np.exp(r['b_share']):.3f} f {r['f_share']:.2f} n_hf {r['n_hf']}", flush=True)
pd.DataFrame(R).to_csv(out, index=False, float_format="%.6g")
