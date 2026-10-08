"""Export the per-storm analysis tables of the three follow-up tests (small, committed) so the numbers can be recomputed without the uncommitted ERA5 caches.
usage: export_followup_tables.py REPO_ROOT ERA5_WORK. ERA5 proxy, pipeline A."""
import sys, os, glob, numpy as np, pandas as pd
root, W = sys.argv[1], sys.argv[2]; E5 = os.path.join(root, "research/era5")
rd = lambda d: pd.concat([pd.read_csv(f) for f in glob.glob(os.path.join(W, d, "stats", "*.csv"))])
# test 3 stage B
S = pd.read_csv(os.path.join(E5, "hf_lifecycle/results/lc_sample.csv")); new = rd("hf_lifecycle"); old = rd("hf_vs_storm")
old = old[old.anchor.isin(["HF_peak", "SF_peak"])].merge(S[["track", "anchor"]], on=["track", "anchor"]); old["anchor"] = "lag+0"; old["snapped"] = old.match_km <= 25
A = pd.concat([new, old]); A["lag"] = A.anchor.str.replace("lag", "").astype(int)
A = A.merge(S[["track", "grp", "season", "mon"]], on="track")
A[["track", "basin", "grp", "season", "mon", "lag", "time", "snapped", "match_km", "g800", "pc", "msl_grad", "g48_rmax", "gmax_r"]].to_csv(os.path.join(E5, "hf_lifecycle/results/lc_values.csv.gz"), index=False)
# test 2
H = pd.read_csv(os.path.join(E5, "hf_pattern_phase/results/hf_index.csv")); H = H[H.in_scope]
st = pd.concat([rd("hf_vs_storm"), rd("hf_pattern_phase")]); st = st[st.anchor.isin(["HF_onset", "HF_peak"])].drop_duplicates(["track", "anchor"])
SC = ["pc", "msl_grad", "gmax", "wsmax", "a_g48", "gmax_r", "g48_rmax", "g48_right", "d2m_500", "gust_factor", "a_hf"]
st[["track", "anchor", "match_km", "g800", "g800_cat", "gmax_near_grn_ice"] + SC].merge(H[["track", "basin", "season", "mon", "z", "tercile", "onset_lat"]], on="track").to_csv(os.path.join(E5, "hf_pattern_phase/results/pp_values.csv.gz"), index=False)
# test 1 stage B
C = rd("hf_conversion_pull"); a = C.anchor.str.extract(r"g(\d)_p(\d+)_L(\d+)").astype(int); C["grp"], C["pair"], C["lag"] = a[0], a[1], a[2]
CS = pd.read_csv(os.path.join(E5, "hf_conversion/results/conv_sample.csv")).rename(columns={"grp": "g"})[["track", "pair", "g", "season", "mon"]]
C.merge(CS, left_on=["track", "pair", "grp"], right_on=["track", "pair", "g"]).drop(columns=["g", "anchor"]).to_csv(os.path.join(E5, "hf_conversion/results/conv_values.csv.gz"), index=False)
for f in ("hf_lifecycle/results/lc_values", "hf_pattern_phase/results/pp_values", "hf_conversion/results/conv_values"):
    print(f, os.path.getsize(os.path.join(E5, f + ".csv.gz")) / 1e6, "MB")
