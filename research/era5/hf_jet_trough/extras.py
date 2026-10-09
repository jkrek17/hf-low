"""After skill.py: HSS and basin split for groups that passed, sensitivity runs, corrected test counts.

usage: SCRATCH=DIR extras.py OUTDIR
Sensitivity (pre-registered, primary test only): J2 with vortex removal 1000 km / none, T2 with trough exclusion 500 / 1000 km, on the seasons
where those columns exist (the first 304 chunks of the pull, the newest seasons, were processed without them), base refitted on the same rows.
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np, pandas as pd
import skill as K
import evaluate as E


def main(out):
    D = pd.read_pickle(os.path.join(os.environ.get("SCRATCH", "/tmp"), "jt_D.pkl"))
    T = pd.read_csv(f"{out}/gains.csv")
    main_ = T[~T.target.str.contains("descriptive")]
    lines = [f"Corrected counts (non-descriptive tests): {int((main_.verdict == 'ADDS SKILL').sum())} of {len(main_)} meet the 'adds skill' rule with the within-family q; "
             f"q(all) < 0.05 in {int((main_.q_all < 0.05).sum())} of {len(main_)}; q(all) < 0.10 in {int((main_.q_all < 0.10).sum())}"]
    # HSS and basin split for passing groups
    ms = K.models_for("0,24")
    for tg, tname, msk in (("rd", "rapid deepening", (D.rd >= 0).values), ("hf48i", "hf48", (D.hf48i >= 0).values)):
        S, P, C, y = K.run(D, {f"{k} [lags 0,24]": v for k, v in ms.items()}, tg, msk)
        for k in ms:
            r = T[(T.target == tname) & (T.model == f"{k} [lags 0,24]")].iloc[0]
            if r.verdict != "ADDS SKILL":
                continue
            d, lo, hi, hb, hk = E.hss_gain(S, P, C, y, f"{k} [lags 0,24]")
            lines.append(f"{tname} {k}: HSS base {hb:.3f} -> {hk:.3f}, change {d:+.3f} [{lo:+.3f}, {hi:+.3f}] (yes/no at the training-count-matched cut)")
            for b in ("atl", "pac"):
                tb, _ = E.gains(S, P, C, y, tname, [f"{k} [lags 0,24]"], mask=(S.basin == b).values)
                rr = tb.iloc[0]
                lines.append(f"   {b}: gain {rr.gain:+.4f} [{rr.lo90:+.4f}, {rr.hi90:+.4f}] seasons better {int(rr.seasons_better)}/{int(rr.n_seasons)} (descriptive)")
    # sensitivity
    has = D.groupby("season").vmax2500_r1000.apply(lambda s: s.notna().mean())
    seas = sorted(has[has > 0.99].index)
    lines.append(f"Sensitivity seasons with the variant columns: {seas[0]}..{seas[-1]} ({len(seas)} seasons; {has.round(2).to_dict()})")
    Dv = D[D.season.isin(seas)].reset_index(drop=True)
    mods = {}
    for k, kw in (("default", {}), ("J2 rm1000", dict(suf_j="r1000")), ("J2 rm0", dict(suf_j="r0")),
                  ("T2 excl500", dict(suf_t="x500")), ("T2 excl1000", dict(suf_t="x1000"))):
        m = K.models_for("0,24", **kw)
        if k == "default":
            mods["J2 default"], mods["T2 default"], mods["J2+T2 default"] = m["J2"], m["T2"], m["J2+T2"]
        elif "J2" in k:
            mods[k] = m["J2"]
        else:
            mods[k] = m["T2"]
    for tg, tname, msk in (("hf24", "hf24", None), ("rd", "rapid deepening", (Dv.rd >= 0).values)):
        S, P, C, y = K.run(Dv, mods, tg, msk)
        tb = K.family(S, P, C, y, tname + " sensitivity", list(mods))
        lines.append(f"-- {tname}, {len(seas)} seasons, base BSS {tb.bss_base.iloc[0]:+.4f}")
        for r in tb.itertuples():
            lines.append(f"   {r.model:16s} gain {r.gain:+.4f} [{r.lo90:+.4f}, {r.hi90:+.4f}] seasons better {r.seasons_better}/{r.n_seasons}")
    open(f"{out}/extras.txt", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main(sys.argv[1])
