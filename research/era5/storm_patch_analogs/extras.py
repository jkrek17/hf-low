"""Post hoc descriptive checks and the analog display (no claims beyond what the pre-registration allows).
usage: ERA5_WORK=DIR extras.py   (run after score_patch.py)
 1. per-season gain of BASE+PC15 over BASE (from the stored out-of-fold probabilities' RPS needs the class probs)
 2. position control: BASE + smooth functions of lat/lon (is the gain just where the storm is?)
 3. analog display: five nearest other-season analogs for every HF-outcome fix; six case storms (lowest minimum pressure)
 4. compact patch product for /mnt/project-files/storm-patches/
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np, pandas as pd
import score_patch as SP
from score_patch import *          # BS, M, Feat, lr, STATE, build, WORK, OUT
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

S, X = SP.build()
y, season = S.gcls.values, S.season.values
Z = np.load(os.path.join(WORK, "oof_probs.npz"), allow_pickle=True)
L = {k: BS.rps(Z[k], y, 4) for k in ("clim", "base", "pc15")}
lines = []
def out(t=""):
    print(t); lines.append(t)
out("== Per-season RPSS gain of BASE+PC15 over BASE (post hoc, descriptive)")
for s in sorted(np.unique(season)):
    m = season == s
    out(f"  {s}-{str(s + 1)[2:]}  n {m.sum():5d}  base {BS.skill(L['base'][m], L['clim'][m]):+.3f}  pc15 {BS.skill(L['pc15'][m], L['clim'][m]):+.3f}  gain {(L['base'][m].sum() - L['pc15'][m].sum()) / L['clim'][m].sum():+.4f}")

# position control: LOSO, BASE + lat/lon polynomials (cos/sin lon harmonics x lat)
lon = np.radians(S.lon.values); lat = (S.lat.values - 50) / 15
pos = np.c_[np.cos(lon), np.sin(lon), np.cos(2 * lon), np.sin(2 * lon), lat * np.cos(lon), lat * np.sin(lon), lat ** 2]
B = S[STATE].values.astype(float)
P = np.zeros((len(S), 4))
for s in sorted(np.unique(season)):
    te = season == s
    P[te] = lr(np.c_[B, pos][~te], y[~te], np.c_[B, pos][te])
Lp = BS.rps(P, y, 4)
rng = np.random.default_rng(SP.SEED + 1)
for lab, a in (("BASE + position terms", Lp), ("BASE + PC15", L["pc15"])):
    pt, lo, hi, _ = BS.boot_gain(season, a, L["base"], L["clim"], rng)
    out(f"  position control: {lab} minus BASE {pt:+.4f} [{lo:+.4f}, {hi:+.4f}]")

# does the pattern add beyond position? (BASE + position terms + PC15, PCA refitted per fold)
P2 = np.zeros((len(S), 4))
fullm = np.ones((NP, NP), bool)
for s_ in sorted(np.unique(season)):
    te = season == s_
    pcs = Feat(X[~te], fullm, None).pcs(X)
    P2[te] = lr(np.c_[B, pos, pcs][~te], y[~te], np.c_[B, pos, pcs][te])
L2_ = BS.rps(P2, y, 4)
pt, lo, hi, _ = BS.boot_gain(season, L2_, Lp, L["clim"], rng)
out(f"  pattern beyond position: BASE + position + PC15 minus BASE + position {pt:+.4f} [{lo:+.4f}, {hi:+.4f}]")
pt, lo, hi, _ = BS.boot_gain(season, L2_, L["base"], L["clim"], rng)
out(f"  BASE + position + PC15 minus BASE {pt:+.4f} [{lo:+.4f}, {hi:+.4f}]")

# analog display
full = np.ones((NP, NP), bool)
F = Feat(X, full, None); W = F.pcs(X, whiten=True)
nn = NearestNeighbors(n_neighbors=400).fit(W)
hf = np.flatnonzero(y == 3)
dist, idx = nn.kneighbors(W[hf])
rows = []
for r, (d, ix) in zip(hf, zip(dist, idx)):
    ok = [j for j in ix if season[j] != season[r]][:5]
    dd = [d[list(ix).index(j)] for j in ok]
    row = dict(track=S.track[r], time=S.time[r], basin=S.basin[r], lat=S.lat[r], lon=S.lon[r], msl=S.msl[r])
    for n, (j, dj) in enumerate(zip(ok, dd), 1):
        row.update({f"a{n}_track": S.track[j], f"a{n}_time": S.time[j], f"a{n}_lat": S.lat[j], f"a{n}_lon": S.lon[j],
                    f"a{n}_class": BS.NAMES[y[j]], f"a{n}_dist": round(float(dj), 2)})
    rows.append(row)
A = pd.DataFrame(rows); A.to_csv(os.path.join(OUT, "analogs_hf_fixes.csv.gz"), index=False)
out(f"\n== Analog display: {len(A)} hurricane-force-outcome fixes, five other-season analogs each")
cls = np.array([[row.get(f"a{n}_class") for n in range(1, 6)] for row in rows])
out("  outcome class of all analogs: " + ", ".join(f"{n} {np.mean(cls == n):.3f}" for n in BS.NAMES))

# six case storms: lowest minimum pressure among the 8 seasons (chosen before looking at analogs)
mn = S.groupby("track").msl.idxmin()
best = S.loc[mn.values].sort_values("msl").head(6)
out("\n== Case storms (lowest minimum pressure at a 00/12 UTC fix)")
for i in best.index:
    ixs = np.flatnonzero(season != season[i]); d = np.linalg.norm(W[ixs] - W[i], axis=1); o = ixs[np.argsort(d)[:5]]
    out(f"  track {S.track[i]} {S.time[i]} {S.basin[i]} {S.msl[i]:.1f} hPa, next-24 h class {BS.NAMES[y[i]]}; analog outcomes: " + ", ".join(BS.NAMES[y[j]] for j in o))
    fig, ax = plt.subplots(1, 6, figsize=(18, 3.4))
    ext = [-500, 500, -500, 500]
    for a_, j, ttl in zip(ax, [i] + list(o), ["storm"] + [f"analog {n}" for n in range(1, 6)]):
        zc = X[j, 0].astype(float)
        im = a_.imshow(zc, origin="lower", extent=ext, cmap="RdBu_r", vmin=-250, vmax=250)
        a_.contour(np.linspace(-500, 500, NP), np.linspace(-500, 500, NP), X[j, 5].astype(float), colors="k", linewidths=0.6, levels=np.arange(940, 1040, 4))
        a_.set_title(f"{ttl}: {S.time[j]} {S.lat[j]:.0f}N {S.lon[j]:.0f}E\n{S.msl[j]:.0f} hPa, next 24 h: {BS.NAMES[y[j]]}", fontsize=8)
        a_.set_xticks([]); a_.set_yticks([])
    fig.suptitle("500 nm patch, north up. Colour: 500 hPa height anomaly (m); lines: sea-level pressure (hPa). ERA5 proxy, pipeline A", fontsize=9)
    fig.savefig(os.path.join(OUT, f"case_{S.track[i]}_{S.time[i]}.png"), dpi=70, bbox_inches="tight"); plt.close(fig)

# compact reusable product
dst = "/mnt/project-files/storm-patches"; os.makedirs(dst, exist_ok=True)
for s in sorted(np.unique(season)):
    m = season == s
    np.savez_compressed(os.path.join(dst, f"patches_{s}_{str(s + 1)[2:]}.npz"), keys=S.key.values[m], x=X[m], channels=np.array(SP.__dict__.get("CH", ["a_z500", "a_msl", "a_u250", "a_v250", "z500", "msl"])),
                        track=S.track.values[m], time=S.time.values[m], lat=S.lat.values[m], lon=S.lon.values[m], gcls=y[m])
open(os.path.join(OUT, "extras.txt"), "w").write("\n".join(lines) + "\n")
