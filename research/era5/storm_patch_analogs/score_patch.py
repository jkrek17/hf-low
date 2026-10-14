"""Pre-registered tests of the storm-centred 500 nm pattern (see PREREGISTRATION.md). ERA5 proxy, pipeline A.
usage: ERA5_WORK=DIR score_patch.py        (patches from extract_patches.py, merged into DIR/patches_all.npz)
Leave-one-season-out over the 8 seasons 2014-15..2021-22 (Oct-Apr). Outlook = 4 ordered classes of the next-24 h peak gust.
"""
import os, sys, glob
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "bomb_env"))
sys.argv = sys.argv[:1]
import numpy as np, pandas as pd
from joblib import Parallel, delayed
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import NearestNeighbors
import score as BS            # bomb_env/score.py: load, rps, brier, boot_gain, bh, skill, clim, NAMES
import model as M
from common_patch import WORK, PR, NP

SEED = 20261014
NPC, K = 15, 30
OUT = os.path.join(HERE, "results")
STATE = M.SETS["full"]


def build():
    D = BS.load()
    D["key"] = D.track.astype(str) + "|" + D.time
    F = pd.read_csv(os.path.join(OUT, "fixes_patch.csv"), dtype={"time": str})
    F["key"] = F.track.astype(str) + "|" + F.time
    P = {}
    for f in glob.glob(os.path.join(WORK, "patches", "*.npz")):
        z = np.load(f)
        for k, x in zip(z["keys"], z["x"]):
            P[str(k)] = x
    assert len(P) == len(F), (len(P), len(F))
    S = D[D.key.isin(set(F.key))].reset_index(drop=True)
    X = np.stack([P[k] for k in S.key]).astype(np.float32)          # [n, 6, 13, 13]
    return S, X


class Feat:
    """Scaler (per-channel SD), PCA and composite, fitted on the training rows only."""
    def __init__(self, X, mask, y, npc=NPC):
        self.sd = X.std((0, 2, 3), keepdims=True) + 1e-9
        self.mask = mask
        Z = self._flat(X)
        self.mu = Z.mean(0)
        self.pca = PCA(npc, svd_solver="full", random_state=0).fit(Z - self.mu)
        self.sdpc = np.sqrt(self.pca.explained_variance_)

    def _flat(self, X):
        Z = np.nan_to_num(X[:, :4] / self.sd[:, :4], nan=0.0)
        return Z[:, :, self.mask].reshape(len(X), -1)

    def pcs(self, X, whiten=False):
        s = self.pca.transform(self._flat(X) - self.mu)
        return s / self.sdpc if whiten else s

    def centred(self, X):
        return self._flat(X) - self.mu


def composite_feature(Xc_tr, hf_tr, season_tr, Xc_te):
    """Projection on the composite difference. Training rows: leave-own-season-out composite; test rows: all training."""
    d_all = Xc_tr[hf_tr].mean(0) - Xc_tr.mean(0)
    ftr = np.zeros(len(Xc_tr))
    for s in np.unique(season_tr):
        m = season_tr == s
        d = Xc_tr[~m & hf_tr].mean(0) - Xc_tr[~m].mean(0)
        ftr[m] = Xc_tr[m] @ d
    return ftr, Xc_te @ d_all


def analog_feats(Wtr, ytr, season_tr, Wte, k=K):
    """Outcome-class fractions of the k nearest library rows. Train rows: library = other training seasons.
    Test rows: library = all training seasons. Returns [n,4] fractions for train and test."""
    ftr = np.zeros((len(Wtr), 4))
    for s in np.unique(season_tr):
        m = season_tr == s
        nn = NearestNeighbors(n_neighbors=k).fit(Wtr[~m])
        idx = nn.kneighbors(Wtr[m], return_distance=False)
        ftr[m] = np.stack([np.bincount(ytr[~m][i], minlength=4) for i in idx]) / k
    nn = NearestNeighbors(n_neighbors=k).fit(Wtr)
    idx = nn.kneighbors(Wte, return_distance=False)
    fte = np.stack([np.bincount(ytr[i], minlength=4) for i in idx]) / k
    return ftr, fte


def lr(Xtr, ytr, Xte, ncls=4):
    prep = M.Prep(Xtr)
    m = LogisticRegression(C=1.0, max_iter=3000).fit(prep(Xtr), ytr)
    p = np.zeros((len(Xte), ncls)); p[:, m.classes_] = m.predict_proba(prep(Xte))
    return p


def fold(s, S, X, Xplc, full_mask, core_mask):
    te = (S.season == s).values; tr = ~te
    y = S.gcls.values; season = S.season.values
    B = S[STATE].values.astype(float); St = S[M.STATE].values.astype(float)
    R = {"clim": BS.clim(S, te, "gcls", 4), "base": lr(B[tr], y[tr], B[te]), "state": lr(St[tr], y[tr], St[te])}
    f_full = Feat(X[tr], full_mask, None)
    pc = f_full.pcs(X)
    R["pc15"] = lr(np.c_[B, pc][tr], y[tr], np.c_[B, pc][te])
    R["state_pc15"] = lr(np.c_[St, pc][tr], y[tr], np.c_[St, pc][te])
    f_core = Feat(X[tr], core_mask, None)
    pcc = f_core.pcs(X)
    R["s1_env"] = lr(np.c_[B, pcc][tr], y[tr], np.c_[B, pcc][te])
    Xc = f_full.centred(X)
    ftr, fte = composite_feature(Xc[tr], (y[tr] == 3), season[tr], Xc[te])
    R["s2_comp"] = lr(np.c_[B[tr], ftr], y[tr], np.c_[B[te], fte])
    W = f_full.pcs(X, whiten=True)
    atr, ate = analog_feats(W[tr], y[tr], season[tr], W[te])
    R["s3_analog"] = lr(np.c_[B[tr], atr[:, 1:]], y[tr], np.c_[B[te], ate[:, 1:]])
    R["analog_alone"] = ate                      # raw k-NN class fractions for the calibration check
    pcp = Feat(Xplc[tr], full_mask, None).pcs(Xplc)       # placebo patches (permuted among like fixes)
    R["s4_placebo"] = lr(np.c_[B, pcp][tr], y[tr], np.c_[B, pcp][te])
    return s, te, R


def main():
    S, X = build()
    rng = np.random.default_rng(SEED)
    # placebo patches: permute whole patches among fixes of the same basin, month, hour
    perm = np.arange(len(S)); pr = np.random.default_rng(SEED)
    grp = S.groupby(["basin", "month", S.time.str[8:10]]).indices
    for _, ix in grp.items():
        perm[ix] = pr.permutation(ix)
    Xplc = X[perm]
    full_mask = np.ones((NP, NP), bool); core_mask = PR >= 300.0
    assert np.isfinite(X[:, :4]).mean() > 0.999, "NaN patches"
    res = Parallel(n_jobs=4)(delayed(fold)(s, S, X, Xplc, full_mask, core_mask) for s in sorted(S.season.unique()))
    keys = [k for k in res[0][2]]
    P = {k: (np.zeros((len(S), 4)) if k != "analog_alone" else np.zeros((len(S), 4))) for k in keys}
    for s, te, R in res:
        for k in keys:
            P[k][te] = R[k]
    y = S.gcls.values; season = S.season.values
    L = {k: BS.rps(P[k], y, 4) for k in keys if k != "analog_alone"}
    lines = []
    def out(t=""):
        print(t); lines.append(t)
    out(f"Fixes {len(S)}, tracks {S.track.nunique()}, seasons {sorted(S.season.unique())} (n = {S.season.nunique()}); HF outcomes {(y==3).sum()} on {S.track[y==3].nunique()} tracks")
    out("Outlook class frequencies: " + ", ".join(f"{n} {f:.3f}" for n, f in zip(BS.NAMES, np.bincount(y, minlength=4) / len(y))))
    out("\n== Outlook RPSS vs basin-month climatology (LOSO, 8 seasons), 90% season-bootstrap CI")
    for k in L:
        if k == "clim":
            continue
        pt, lo, hi, _ = BS.boot_gain(season, L[k], L["clim"], L["clim"], rng)
        out(f"  {k:12s} {BS.skill(L[k], L['clim']):+.4f}  [{lo:+.4f}, {hi:+.4f}]")
    # reference: 22-season fit on the same fixes
    ref = pd.read_csv(os.path.join(HERE, "..", "bomb_env", "results", "oof_outlook_full.csv.gz"), dtype={"time": str})
    ref["key"] = ref.track.astype(str) + "|" + ref.time
    r = ref.set_index("key").loc[S.key]
    Lref = BS.rps(r[["p0", "p1", "p2", "p3"]].values, y, 4)
    out(f"  (reference) full model, 22-season fit, same fixes: {BS.skill(Lref, L['clim']):+.4f}  (not like-for-like: more training seasons)")
    rows = []
    def gain(tag, a, b):
        pt, lo, hi, p = BS.boot_gain(season, L[a], L[b], L["clim"], rng)
        rows.append(dict(test=tag, model=a, over=b, gain=pt, lo90=lo, hi90=hi, p=p)); return pt, lo, hi, p
    pt, lo, hi, p = gain("T1", "pc15", "base")
    verdict = "pattern adds skill" if (pt >= 0.010 and lo > 0) else ("no useful gain" if hi < 0.005 else "cannot tell")
    out(f"\n== PRIMARY T1: BASE + PC15 minus BASE: {pt:+.4f} [{lo:+.4f}, {hi:+.4f}]  verdict (pre-registered rule): {verdict}")
    sec = [("S1 environment only (r >= 300 km)", "s1_env"), ("S2 composite projection", "s2_comp"), ("S3 analog features", "s3_analog")]
    ps = []
    out("\n== Secondary (gain over BASE; BH-FDR over S1-S3)")
    tmp = [gain(n, k, "base") for n, k in sec]
    q = BS.bh([t[3] for t in tmp])
    for (n, k), t, qq in zip(sec, tmp, q):
        out(f"  {n:36s} {t[0]:+.4f} [{t[1]:+.4f}, {t[2]:+.4f}]  p {t[3]:.3f}  q {qq:.3f}")
        rows[[x['test'] for x in rows].index(n)]["q"] = qq
    t = gain("S4 placebo (not in FDR family)", "s4_placebo", "base")
    out(f"  S4 placebo (permuted patches)       {t[0]:+.4f} [{t[1]:+.4f}, {t[2]:+.4f}]  (over-fit cost of 15 extra columns)")
    t = gain("D state + PC15 over state", "state_pc15", "state")
    out(f"\n== Descriptive: state-only + PC15 over state-only {t[0]:+.4f} [{t[1]:+.4f}, {t[2]:+.4f}]")
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "tests.csv"), index=False)
    out("\n== Gain of BASE+PC15 over BASE by class (Brier skill gain) and by subset")
    for c in range(4):
        yc = (y == c).astype(float)
        a, b, cl = BS.brier(P["pc15"][:, c], yc), BS.brier(P["base"][:, c], yc), BS.brier(P["clim"][:, c], yc)
        pt, lo, hi, _ = BS.boot_gain(season, a, b, cl, rng)
        out(f"  {BS.NAMES[c]:16s} base BSS {BS.skill(b, cl):+.3f}  pc15 BSS {BS.skill(a, cl):+.3f}  gain {pt:+.4f} [{lo:+.4f}, {hi:+.4f}]")
    for lab, m in (("Atlantic", (S.basin == "atl").values), ("Pacific", (S.basin == "pac").values), ("track continues 24 h", (~S.left24).values)):
        pt, lo, hi, _ = BS.boot_gain(season[m], L["pc15"][m], L["base"][m], L["clim"][m], rng)
        out(f"  RPSS gain {lab:22s} n {m.sum():6d}  {pt:+.4f} [{lo:+.4f}, {hi:+.4f}]")
    # analog calibration
    pa = P["analog_alone"][:, 3]; yh = (y == 3).astype(float)
    pcl = P["clim"][:, 3]
    bss = 1 - ((pa - yh) ** 2).sum() / ((pcl - yh) ** 2).sum()
    pt, lo, hi, _ = BS.boot_gain(season, (pa - yh) ** 2, (pcl - yh) ** 2, (pcl - yh) ** 2, rng)
    pbase = P["base"][:, 3]
    bss_base = 1 - ((pbase - yh) ** 2).sum() / ((pcl - yh) ** 2).sum()
    out(f"\n== Analog calibration (k = {K}, stand-alone P(hurricane force next 24 h) = share of the {K} nearest past fixes that reached it)")
    out(f"  Brier skill vs basin-month climatology: analog alone {bss:+.3f} [{lo:+.3f}, {hi:+.3f}]; BASE logistic {bss_base:+.3f}")
    edges = [0, 1e-9, 0.034, 0.067, 0.1, 0.167, 0.267, 0.4, 1.01]
    b = np.digitize(pa, edges[1:-1])
    rel = []
    for j in range(len(edges) - 1):
        m = b == j
        if m.sum():
            rel.append(dict(bin=f"{edges[j]:.2f}-{edges[j+1]:.2f}", n=int(m.sum()), n_hf=int(yh[m].sum()), forecast=pa[m].mean(), observed=yh[m].mean()))
    rel = pd.DataFrame(rel); rel.to_csv(os.path.join(OUT, "analog_reliability.csv"), index=False)
    out(rel.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    slope = np.polyfit(rel.forecast, rel.observed, 1, w=np.sqrt(rel.n))[0]
    out(f"  reliability slope {slope:.2f} (1 = perfect)")
    c20 = ((pa >= 0.25) & (pa <= 0.35)); out(f"  'about 6 in 20' check: fixes where 25-35% of neighbours reached HF: n {c20.sum()}, observed HF rate {yh[c20].mean():.3f}")
    open(os.path.join(OUT, "skill.txt"), "w").write("\n".join(lines) + "\n")
    np.savez_compressed(os.path.join(WORK, "oof_probs.npz"), **{k: P[k] for k in P}, keys=S.key.values)
    S[["key", "track", "time", "basin", "season", "lat", "lon", "msl", "gcls"]].assign(**{f"p_{k}": P[k][:, 3] for k in ("base", "pc15", "s1_env", "s2_comp", "s3_analog", "s4_placebo")}).to_csv(os.path.join(OUT, "oof_hf_prob.csv.gz"), index=False)


if __name__ == "__main__":
    main()
