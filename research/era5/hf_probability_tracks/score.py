"""Per-fix P(HF within 24 h) along ERA5 pipeline A tracks, 1979-2025 (ERA5 PROXY).

usage: python3 -I score.py        (reads committed files only; writes results/)

Models F, N, S as defined in PREREGISTRATION.md (the PR 12 code, unchanged):
  F  full (state + Hart + environment, includes g800)      2004+ only (needs the environment table)
  N  full_nogust                                           2004+ only
  S  state_nogust (storm state without g800)               2004+ (leave-one-season-out) and 1979-2003 (final fit)
2004-05..2025-26: leave-one-season-out probabilities, so no fix is scored by a model that saw its season.
1979-80..2003-04: S fitted on all 2004+ fixes, applied to the 00/12 UTC fixes of the catalog tracks
(events and matched null cases). Not the full low population; see the pre-registration.
"""
import os, sys, json
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
INT = os.path.join(HERE, "..", "intensity")
HFH = os.path.join(HERE, "..", "hf_history", "results")
sys.path.insert(0, INT)
import model as M
import fixes as FX

RES = os.path.join(HERE, "results")
STATE_NOGUST = [c for c in M.STATE if c != "g800"]
SETS = {"F": M.SETS["full"], "N": M.SETS["full_nogust"], "S": STATE_NOGUST}


def prepare(D):
    """The derived columns model.load() adds (copied so pre-2004 fixes get the same ones)."""
    D = D.copy()
    D["young"] = D.dp12.isna().astype(float)
    D["dp12"] = D.dp12.fillna(0.0)
    D["logage"] = np.log1p(D.age)
    D["pac"] = (D.basin == "pac").astype(float)
    doy = pd.to_datetime(D.time, format="%Y%m%d%H").dt.dayofyear.values
    D["doy_c"], D["doy_s"] = np.cos(2 * np.pi * doy / 365.25), np.sin(2 * np.pi * doy / 365.25)
    D["month"] = D.time.str[4:6].astype(int)
    return D


def brier(p, y):
    return (p - y) ** 2


def table(fc, ob):
    return M.table(fc, ob)


def main():
    os.makedirs(RES, exist_ok=True)
    D = M.load(os.path.join(INT, "results", "fixes_2004.csv.gz"), os.path.join(INT, "results", "env_2004.csv.gz"))
    D["hf24"] = D.hf24.astype(int)
    y = D.hf24.values
    seas = D.season.values
    rep = []
    w = rep.append
    w(f"Fixes {len(D)} on {D.track.nunique()} tracks, seasons {seas.min()}-{seas.max()} (n = {len(np.unique(seas))}); "
      f"hf24 base rate {y.mean():.4f}")
    # ---- leave-one-season-out for F, N, S and basin-month climatology
    P = {k: np.zeros(len(D)) for k in ["clim"] + list(SETS)}
    CUT = {k: np.zeros(len(D)) for k in SETS}
    for s in np.unique(seas):
        te = seas == s
        tr = ~te
        Dtr, Dte = D[tr], D[te]
        P["clim"][te] = M.clim_probs(Dtr, Dte.reset_index(drop=True).set_index(np.arange(te.sum())), "hf24", 2)[:, 1]
        for name, cols in SETS.items():
            prep = M.Prep(Dtr[cols].values.astype(float))
            m = M.fit(prep(Dtr[cols].values.astype(float)), Dtr.hf24.values)
            P[name][te] = m.predict_proba(prep(Dte[cols].values.astype(float)))[:, 1]
            ptr = m.predict_proba(prep(Dtr[cols].values.astype(float)))[:, 1]
            CUT[name][te] = np.quantile(ptr, 1 - Dtr.hf24.mean())
    # ---- reproduction of PR 12 (F) and skill of N, S
    bc = brier(P["clim"], y).sum()
    w("\nLeave-one-season-out Brier skill vs basin-month climatology, hf24 (PR 12 `full` = 0.423, HSS 0.59):")
    for name in SETS:
        bss = 1 - brier(P[name], y).sum() / bc
        t = table(P[name] >= CUT[name], y == 1)
        w(f"  {name}: BSS {bss:.4f}  POD {t['pod']:.2f} FAR {t['far']:.2f} CSI {t['csi']:.2f} HSS {t['hss']:.2f} bias {t['bias']:.2f}")
    # ---- per-fix table, 2004+
    out = D[["track", "time", "basin", "season", "lat", "lon", "msl", "g800", "hf_now", "hf24"]].copy()
    for name in SETS:
        out["P_" + name] = P[name].round(5)
        out["cut_" + name] = CUT[name].round(5)
    out.to_csv(os.path.join(RES, "fix_probs_2004.csv.gz"), index=False)
    # ---- final fit of S and N on all 2004+ fixes, applied to 1979-2003 catalog fixes
    ct = pd.read_csv(os.path.join(HFH, "era5_hf_catalog_tracks.csv"), dtype={"time": str})
    cat = pd.read_csv(os.path.join(HFH, "era5_hf_catalog.csv"))
    y_, m_ = ct.time.str[:4].astype(int), ct.time.str[4:6].astype(int)
    ct["season"] = np.where(m_ >= 6, y_, y_ - 1)
    Fx = FX.build(ct)
    Fx = prepare(Fx[Fx.season < 2004])
    prep = M.Prep(D[STATE_NOGUST].values.astype(float))
    m = M.fit(prep(D[STATE_NOGUST].values.astype(float)), D.hf24.values)
    Fx["P_S"] = m.predict_proba(prep(Fx[STATE_NOGUST].values.astype(float)))[:, 1]
    # the same model with g800, shown only as a labelled within-era contrast
    cg = M.SETS["state"]
    prep2 = M.Prep(D[cg].values.astype(float))
    m2 = M.fit(prep2(D[cg].values.astype(float)), D.hf24.values)
    Fx["P_Sg"] = m2.predict_proba(prep2(Fx[cg].values.astype(float)))[:, 1]
    role = cat.set_index("track").role
    Fx["role"] = Fx.track.map(role)
    keep = ["track", "time", "basin", "season", "lat", "lon", "msl", "g800", "hf_now", "hf24", "role", "P_S", "P_Sg"]
    Fx[keep].round({"P_S": 5, "P_Sg": 5}).to_csv(os.path.join(RES, "fix_probs_pre2004.csv.gz"), index=False)
    w(f"\nPre-2004: {len(Fx)} fixes on {Fx.track.nunique()} catalog tracks, seasons {int(Fx.season.min())}-{int(Fx.season.max())} "
      f"({len(cat)} catalog tracks over 1979-2025)")
    w(f"  mean P_S over fixes: events {Fx[Fx.role=='event'].P_S.mean():.3f}, null cases {Fx[Fx.role=='null_case'].P_S.mean():.3f}")
    w(f"  mean P_S - P_Sg (what the gust predictor adds): {np.mean(Fx.P_Sg - Fx.P_S):+.4f}; "
      f"1979-2000 {np.mean((Fx.P_Sg - Fx.P_S)[Fx.season <= 2000]):+.4f}, 2001-2003 {np.mean((Fx.P_Sg - Fx.P_S)[Fx.season >= 2001]):+.4f}")
    # N and F on the 2004+ catalog fixes, to compare with S on the same fixes (stand-in check)
    cat_tr = set(cat.track)
    sub = out[out.track.isin(cat_tr)]
    w(f"\n2004+ catalog fixes {len(sub)}: correlation of LOSO P_S with P_N {np.corrcoef(sub.P_S, sub.P_N)[0,1]:.3f}, "
      f"with P_F {np.corrcoef(sub.P_S, sub.P_F)[0,1]:.3f}")
    open(os.path.join(RES, "score.txt"), "w").write("\n".join(rep) + "\n")
    print("\n".join(rep))


if __name__ == "__main__":
    main()
