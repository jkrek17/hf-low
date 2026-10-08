"""Storm-scale composites and structure scalars: HF onset / HF peak against storm-force-only peak (hf_vs_storm, PREREGISTRATION.md Part 2).

ERA5 proxy, pipeline A, ARCO-ERA5 0.25 degree. Reads $ERA5_WORK/hf_vs_storm/{stats,boxes} written by extract_storm.py.
usage: analyse_storm.py REPO_ROOT [B_SCALAR] [B_PIXEL]
Outputs results/storm/: anchors.csv (one row per anchor), scalars.csv, scalars_secondary.csv, comp_<basin>_<cmp>.npz, summary.txt
"""
import sys, os, glob, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from strat import cells, boot, bh
root = sys.argv[1]
BS = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
BP = int(sys.argv[3]) if len(sys.argv) > 3 else 1000
WORK = os.path.join(os.environ.get("ERA5_WORK", "/tmp/era5work"), "hf_vs_storm")
OUT = os.path.join(root, "research/era5/hf_vs_storm/results/storm"); os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(20261012)
SEASONS = np.arange(2004, 2026)

R = pd.concat(pd.read_csv(f) for f in sorted(glob.glob(os.path.join(WORK, "stats", "*.csv"))))
S = pd.read_csv(os.path.join(root, "research/era5/hf_vs_storm/results/storms.csv"))
R = R.merge(S[["track", "season", "mon", "sf_hi", "start", "onset_time", "peak_time"]], on="track")
R["grp"] = np.where(R.anchor == "SF_peak", 1, 0)
tt = pd.to_datetime(R.time.astype(str), format="%Y%m%d%H"); st = pd.to_datetime(R.start.astype(str), format="%Y%m%d%H")
R["age_h"] = (tt - st).dt.total_seconds() / 3600
lines = []
chk = R[R.g800_cat > 0]
lines.append(f"anchors {len(R)} ({R.anchor.value_counts().to_dict()}); re-detection match <= 25 km: {(R.match_km <= 25).mean():.4f}; "
             f"g800 within 0.5 kt of catalog: {(np.abs(chk.g800 - chk.g800_cat) <= 0.5).mean():.4f} of {len(chk)}; within 1 kt {(np.abs(chk.g800 - chk.g800_cat) <= 1).mean():.4f}")
R.to_csv(os.path.join(OUT, "anchors.csv"), index=False)

# ------------------------------------------------------------- scalars
SC = [("pc", "central MSLP (hPa)"), ("msl_grad", "MSLP gradient 0-500 km (hPa/100 km)"), ("gmax_r", "radius of max gust (km)"),
      ("g48_rmax", "outer radius of 48-kt gust area (km)"), ("g48_right", "share of 48-kt area right of motion"),
      ("d2m_500", "2 m dewpoint within 500 km (C)"), ("gust_factor", "gust factor (max gust / max 10 m wind)")]
DESC = [("gmax", "max gust (kt)"), ("wsmax", "max 10 m wind (kt)"), ("a_g48", "area gust >= 48 kt (km2)"), ("a_hf", "area gust >= 71.7 kt (km2)"),
        ("age_h", "storm age at anchor (h)"), ("own_ocean_frac", "owned ocean fraction")]

def scal(df, hfkind, cols, strat="mon", nstr=12, bins=None):
    d = df[df.anchor.isin([hfkind, "SF_peak"])]
    rows = []
    for b in ("atl", "pac"):
        x = d[d.basin == b]
        if strat == "pc":
            x = x.copy(); x["stratum"] = np.clip(np.floor((x.pc - 920) / 5), 0, 15).astype(int)
            # overlap: bins where both groups hold >= 5 storms
            cnt = x.groupby(["stratum", "grp"]).size().unstack(fill_value=0)
            keep = cnt.index[(cnt.get(0, 0) >= 5) & (cnt.get(1, 0) >= 5)]
            x = x[x.stratum.isin(keep)]
        else:
            x = x.assign(stratum=x.mon - 1)
        for c, lab in cols:
            v = x[c].values.astype(float)
            sm, cn = cells(v, x.grp.values, x.season.values, x.stratum.values, SEASONS, nstr)
            obs, reps, p = boot(sm, cn, BS, rng)
            ok = ~np.isnan(v)
            sd = np.sqrt((np.nanvar(v[(x.grp == 0).values & ok]) + np.nanvar(v[(x.grp == 1).values & ok])) / 2)
            rows.append(dict(basin=b, cmp="C1" if hfkind == "HF_onset" else "C2", scalar=c, label=lab, hf=obs[1], sf=obs[2], diff=obs[0],
                             lo=np.percentile(reps, 2.5), hi=np.percentile(reps, 97.5), se=reps.std(), p=p, mde80=2.8 * reps.std(),
                             sd=sd, d_std=obs[0] / sd if sd > 0 else np.nan, mde_std=2.8 * reps.std() / sd if sd > 0 else np.nan,
                             n_hf=int(((x.grp == 0) & ok).sum()), n_sf=int(((x.grp == 1) & ok).sum())))
    return pd.DataFrame(rows)

def with_q(T, key="p"):
    T = T.copy(); T["q"] = bh(T[key].values)[0]; return T

prim = with_q(pd.concat([scal(R, "HF_onset", SC), scal(R, "HF_peak", SC)]))
prim.to_csv(os.path.join(OUT, "scalars.csv"), index=False)
lines.append(f"primary scalar family: {int((prim.q < 0.05).sum())} of {len(prim)} pass BH q<0.05")
desc = pd.concat([scal(R, "HF_onset", DESC), scal(R, "HF_peak", DESC)]); desc.to_csv(os.path.join(OUT, "scalars_descriptive.csv"), index=False)
sec = []
# S1 depth-matched (central MSLP excluded: it is the matching variable)
s1 = with_q(pd.concat([scal(R, "HF_onset", SC[1:], "pc", 16), scal(R, "HF_peak", SC[1:], "pc", 16)])); s1["family"] = "S1_depth_matched"
# S2 SF-high only
R2 = R[(R.anchor != "SF_peak") | (R.sf_hi)]
s2 = with_q(pd.concat([scal(R2, "HF_onset", SC), scal(R2, "HF_peak", SC)])); s2["family"] = "S2_SF_high"
# S3 drop storms with the maximum gust within 100 km of Greenland / Iceland (Atlantic only)
R3 = R[~R.gmax_near_grn_ice.astype(bool)]
s3 = pd.concat([scal(R3, "HF_onset", SC), scal(R3, "HF_peak", SC)]); s3 = s3[s3.basin == "atl"]; s3 = with_q(s3); s3["family"] = "S3_no_Greenland_Iceland"
pd.concat([s1, s2, s3]).to_csv(os.path.join(OUT, "scalars_secondary.csv"), index=False)
for f, T in (("S1", s1), ("S2", s2), ("S3", s3)):
    lines.append(f"{f}: {int((T.q < 0.05).sum())} of {len(T)} pass q<0.05")
lines.append("share of storms whose max gust is within 100 km of Greenland/Iceland: " + str(R.groupby("anchor").gmax_near_grn_ice.mean().round(3).to_dict()))

# ------------------------------------------------------------- composites
FIELDS = ["gust", "ws", "msl", "d2m"]
boxes = {}
for f in sorted(glob.glob(os.path.join(WORK, "boxes", "*.npz"))):
    z = np.load(f)
    for k in z.files:
        boxes[k] = z[k]
print("boxes", len(boxes), flush=True)
for b in ("atl", "pac"):
    for cmpname, hfkind in (("C1", "HF_onset"), ("C2", "HF_peak")):
        x = R[(R.basin == b) & R.anchor.isin([hfkind, "SF_peak"])].reset_index(drop=True)
        out = {}
        for frame in ("m", "n"):
            for fld in FIELDS:
                arr = np.stack([boxes[f"{r.track}_{r.anchor}_{fld}_{frame}"].astype(np.float32).ravel() for r in x.itertuples()])
                sm, cn = cells(arr, x.grp.values, x.season.values, (x.mon.values - 1), SEASONS, 12)
                nb = BP if frame == "m" else 100
                obs, _, p = boot(sm, cn, nb, rng, keep=False)
                if frame == "m":
                    qq, rej = bh(p)
                else:
                    rej = np.zeros(p.shape, bool)
                out[f"{fld}_{frame}_diff"] = obs[0].reshape(121, 121).astype(np.float32)
                out[f"{fld}_{frame}_hf"] = obs[1].reshape(121, 121).astype(np.float32)
                out[f"{fld}_{frame}_sf"] = obs[2].reshape(121, 121).astype(np.float32)
                out[f"{fld}_{frame}_rej"] = rej.reshape(121, 121)
                out[f"{fld}_{frame}_p"] = p.reshape(121, 121).astype(np.float32)
                lines.append(f"{b} {cmpname} {fld} {frame}: {int(rej.sum())} of {rej.size} pixels pass BH q<0.05")
        nrot = x[x.heading.notna()]
        out["n_hf"] = int((x.grp == 0).sum()); out["n_sf"] = int((x.grp == 1).sum())
        out["n_hf_rot"] = int((nrot.grp == 0).sum()); out["n_sf_rot"] = int((nrot.grp == 1).sum())
        np.savez_compressed(os.path.join(OUT, f"comp_{b}_{cmpname}.npz"), **out)
open(os.path.join(OUT, "summary.txt"), "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
