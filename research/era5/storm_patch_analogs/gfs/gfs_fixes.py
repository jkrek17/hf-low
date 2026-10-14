"""Matched lead 24 h GFS fixes (PR 185 rule) -> results/gfs_fixes_lead24.csv. Same code path as gfs_outlook/score_outlook_full.py:
GFS tracker fixes at valid 00/12 UTC, matched to the nearest pipeline A fix at the same valid time within 500 km; keeps those with all GFS
and ERA5 environment terms (the PR 185 scored sample, 18,167 expected). No GFS grid is read."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); RE = os.path.join(HERE, "..", "..")
sys.path.insert(0, os.path.join(RE, "gfs_outlook")); sys.path.insert(0, os.path.join(RE, "gfs_test")); sys.path.insert(0, os.path.join(RE, "intensity"))
import numpy as np, pandas as pd
_argv = sys.argv; sys.argv = [_argv[0]]
import score_outlook_full as SF
sys.argv = _argv
SO, AS, BS = SF.SO, SF.AS, SF.BS
HART, ENVC, ALLENV, REQ = SF.HART, SF.ENVC, SF.ALLENV, SF.REQ


def main():
    F, n_all = AS.gfs_features()
    F = F[F.lead.isin([0, 12, 24, 36, 48])].copy().reset_index(drop=True)
    D = BS.load(); E = D[D.season.isin(SF.SEASONS)].copy().reset_index(drop=True); E["time_i"] = E.time.astype(int)
    mi, md = AS.match_era5(F, E)
    F["era_i"] = mi; F["matched"] = mi >= 0
    Fm = F[F.matched].copy().reset_index(drop=True); Fm["cycle"] = Fm.cycle.astype(str)
    T1 = pd.read_csv(SF.ENV1, dtype={"cycle": str}); T2 = pd.read_csv(SF.ENVF, dtype={"cycle": str})
    Fm = Fm.merge(T1[["cycle", "lead", "track", "jet250", "vadv500", "eady"]], on=["cycle", "lead", "track"], how="left")
    Fm = Fm.merge(T2[["cycle", "lead", "track"] + HART + ["div300", "sst", "sstgrad", "sst_t500", "flux", "tcwv"]], on=["cycle", "lead", "track"], how="left")
    Fm["nosst"] = Fm.sst.isna().astype(float)
    flux_ok = (Fm.lead == 0) | Fm.flux.notna()
    Fm["gfs_ok"] = Fm[REQ].notna().all(axis=1) & flux_ok
    ei = Fm.era_i.values
    Ee = E.reindex(ei)
    Fm["era_ok"] = Ee[REQ + ["flux"]].notna().all(axis=1).values
    Fm["era_track"] = Ee.track.values; Fm["era_time"] = Ee.time.values
    Fm["y"] = Ee.gcls.values
    Fm["ok"] = Fm.gfs_ok & Fm.era_ok
    S = Fm[(Fm.lead == 24) & Fm.ok].copy().reset_index(drop=True)
    print("matched lead 24 scored sample:", len(S), "(PR 185: 18167); cycles", S.cycle.nunique(), "seasons", S.season.value_counts().sort_index().to_dict())
    S.to_csv(os.path.join(HERE, "results", "gfs_fixes_lead24.csv.gz"), index=False)
    Fm.to_pickle("/tmp/claude-0/work/Fm_all.pkl")


if __name__ == "__main__":
    main()
