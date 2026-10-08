"""Print the headline numbers of every tier from results/*.json (what README.md quotes). usage: summarize.py RESULTS_DIR"""
import json, os, sys
R = sys.argv[1]
f = lambda x, n=3: "nan" if x is None else f"{x:.{n}f}"
for t in ("arch", "r1", "r2", "r3"):
    r = json.load(open(os.path.join(R, f"{t}_results.json")))
    for b, X in r["basins"].items():
        m = X["M1"]; ci = X["M1_ci"]; k = X["M3"]; g = X["M2"]
        print(f"[{t} {b}] events {X['n_events']} (outside window {X['n_outside_window']}); per season {X['per_season_mean']:.1f}")
        print(f"  W7  N1 E {f(m['N1']['E'][3])} CI {f(ci['N1']['W7'][0])}..{f(ci['N1']['W7'][1])} p {f(m['N1']['p'][3])} | N2 E {f(m['N2']['E'][3])} CI {f(ci['N2']['W7'][0])}..{f(ci['N2']['W7'][1])} | IDX E {f(m['IDX']['E'][3])}")
        print(f"  W30 N1 E {f(m['N1']['E'][5])} CI {f(ci['N1']['W30'][0])}..{f(ci['N1']['W30'][1])} p {f(m['N1']['p'][5])} | N2 E {f(m['N2']['E'][5])} CI {f(ci['N2']['W30'][0])}..{f(ci['N2']['W30'][1])} | IDX E {f(m['IDX']['E'][5])} IDX2 {f(m['IDX2']['E'][5])} SPLIT {f(m['SPLIT']['E'][5])}")
        print(f"  curve E (1,2,3,7,14,30 d, N1): {[round(x,3) for x in m['N1']['E']]}")
        print(f"  M2 gaps<=48h {g['obs']['G48']} vs {g['null_mean']['G48']:.1f} (ratio {g['ratio_G48']:.3f}) p {f(g['p']['G48'])}; CV {f(g['obs']['CV'])} (null {f(g['null_mean']['CV'])}) p {f(g['p']['CV'])}; KS {f(g['obs']['KS'])} p {f(g['p']['KS'])}; share<=72h {f(g['obs']['S72'])} (null {f(g['null_mean']['S72'])})")
        print(f"  M3 Knox 72h/1000km obs {k['obs'][1][1]} null {k['null_mean'][1][1]:.1f} ratio {f(k['ratio'][1][1],2)} p {f(k['p'][1][1],4)}; pairs<=2000km same season {k['pairs']}")
        c = X["index_coef"]
        print(f"  index {c['index']} RR/SD {f(c['RR'])} se(log) {f(c['se_cluster'])} perm p {f(c['p_perm'],4)}")
        if "M4" in X: print("  M4", X["M4"])
    print("  tests (basin, test, p, q):", [(x["basin"], x["test"], round(x["p"], 4), round(x["q_bh"], 4)) for x in r["tests"]])
